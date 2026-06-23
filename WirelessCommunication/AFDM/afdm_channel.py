"""
afdm_channel.py
================
Helper module for the AFDM modulation/demodulation notes notebook.

This module deliberately keeps everything that is *not* the AFDM
modulation/demodulation core, so that the companion notebook can focus on the
DAFT/IDAFT transmit-receive chain. It contains:

  * NR-like numerology  ............ nr_numerology()
  * Doubly-dispersive channel ...... gen_doubly_dispersive_channel(), apply_channel()
  * Pulse-shaping / interpolation .. frac_delay_kernel(), build_channel_matrix()
  * Constellation mappers .......... qam16_map/demap, qpsk_map/demap
  * DMRS (QPSK pilots) ............. gen_dmrs()
  * Noise + equalization ........... awgn(), lmmse_equalize()
  * DMRS channel estimation ........ estimate_channel_from_dmrs()

References (equation numbers) follow:
    H. S. Rou et al., "AFDM: Evolving OFDM Towards 6G+", IEEE OJ-COMS.

All comments/labels are kept in English to match the repository convention.
"""

import numpy as np

# ----------------------------------------------------------------------------
# 1. NR-like numerology (100 MHz @ 30 kHz SCS, FR1)
# ----------------------------------------------------------------------------
def nr_numerology():
    """Return a dict of 5G NR-like parameters for 100 MHz / 30 kHz SCS.

    273 resource blocks (RB) x 12 subcarriers = 3276 active subcarriers,
    placed inside a 4096-point FFT. Sampling rate fs = N_fft * SCS.
    """
    scs = 30e3                  # subcarrier spacing [Hz]
    N_fft = 4096                # FFT size
    N_rb = 273                  # resource blocks for 100 MHz @ 30 kHz
    N_sc = N_rb * 12            # active subcarriers = 3276
    fs = N_fft * scs           # sampling rate = 122.88 MHz
    return {
        "scs": scs,
        "bw": 100e6,
        "N_fft": N_fft,
        "N_rb": N_rb,
        "N_sc": N_sc,
        "fs": fs,
        "Ts": 1.0 / fs,
        "fc": 4.0e9,            # assumed carrier frequency [Hz]
    }


# ----------------------------------------------------------------------------
# 2. Doubly-dispersive channel (delay + Doppler)
# ----------------------------------------------------------------------------
def gen_doubly_dispersive_channel(N, fs, fc=4.0e9, velocity_kmh=500.0,
                                  scenario="highmobility", seed=0):
    """Generate a P-path doubly-dispersive channel (eq. (1)).

    Parameters
    ----------
    N  : int    block length in samples (used to normalize Doppler, eq. (4)).
    fs : float  sampling rate [Hz].
    fc : float  carrier frequency [Hz] (sets the max Doppler).
    velocity_kmh : float  mobile speed; large value -> significant Doppler so
                   that the AFDM advantage over OFDM is visible.
    scenario : "highmobility" uses a short TDL-like power-delay profile.
    seed : RNG seed for the per-path random Doppler angles.

    Returns
    -------
    paths : dict with physical and normalized per-path parameters
            h   : complex path gains          (P,)
            tau : path delays  [s]            (P,)
            nu  : Doppler shifts [Hz]         (P,)
            ell : normalized delays  l_p = tau/Ts        (P,)   eq.(4)
            fdop: normalized digital Doppler f_p = N*nu/fs (P,)  eq.(4)
    """
    rng = np.random.default_rng(seed)
    c = 3e8
    nu_max = velocity_kmh / 3.6 / c * fc      # maximum Doppler [Hz]

    # Short power-delay profile (relative power dB, delay in ns): TDL-like.
    if scenario == "highmobility":
        delays_ns = np.array([0.0, 30.0, 70.0, 90.0, 110.0, 190.0])
        powers_db = np.array([0.0, -1.5, -1.4, -3.6, -0.6, -9.1])
    else:
        delays_ns = np.array([0.0, 50.0, 120.0])
        powers_db = np.array([0.0, -3.0, -6.0])

    P = len(delays_ns)
    tau = delays_ns * 1e-9
    lin = 10 ** (powers_db / 20.0)
    # Rayleigh-distributed complex gains scaled by the PDP.
    phase = rng.uniform(0, 2 * np.pi, P)
    h = lin * np.exp(1j * phase)
    h = h / np.sqrt(np.sum(np.abs(h) ** 2))   # normalize total power to 1

    # Each path gets a Doppler nu_p = nu_max * cos(theta_p) (Jakes-like).
    theta = rng.uniform(0, 2 * np.pi, P)
    nu = nu_max * np.cos(theta)

    ell = tau * fs                              # normalized delay (eq. 4)
    fdop = N * nu / fs                          # normalized digital Doppler (eq. 4)

    return {
        "h": h, "tau": tau, "nu": nu,
        "ell": ell, "fdop": fdop,
        "nu_max": nu_max, "P": P,
    }


def make_paths(h, ell, fdop):
    """Build a paths dict directly from normalized arrays.

    Convenient for the illustrative/reproduction cases (e.g. the paper's Fig. 3
    setup) and for the BER study, where normalized delay l_p and digital
    Doppler f_p are specified directly rather than derived from a physical PDP.
    """
    h = np.asarray(h, dtype=complex)
    return {
        "h": h,
        "ell": np.asarray(ell, dtype=float),
        "fdop": np.asarray(fdop, dtype=float),
        "tau": np.zeros_like(h, dtype=float),
        "nu": np.zeros_like(h, dtype=float),
        "P": len(h),
    }


def apply_channel(s, paths, N, kernel="sinc", L=16):
    """Pass a time-domain block s through the doubly-dispersive channel.

    Implements the discrete input-output relation (eq. (7)):

        r[n] = sum_p h_p * exp(j 2 pi f_p n / N) * sum_m s[m] g((n-m)-l_p)

    using an efficient time-domain interpolation instead of an N x N matrix,
    so it scales to the NR FFT size N = 4096. A cyclic-prefix style circular
    indexing is assumed (the caller adds/removes the prefix).

    Parameters
    ----------
    s      : (N,) complex transmit samples (already prefix-removed reference).
    paths  : dict from gen_doubly_dispersive_channel().
    N      : block length.
    kernel : "sinc" (band-limited), "rc" (raised cosine), or "rect".
    L      : interpolation half-window (taps on each side) for fractional delay.

    Returns
    -------
    r : (N,) complex received samples (noise added separately via awgn()).
    """
    n = np.arange(N)
    r = np.zeros(N, dtype=complex)
    for hp, ell, fdop in zip(paths["h"], paths["ell"], paths["fdop"]):
        ell_int = int(np.floor(ell))
        frac = ell - ell_int
        # Fractional-delay FIR via the chosen interpolation kernel.
        taps = np.arange(-L, L + 1)
        g = frac_delay_kernel(taps - frac, kernel)
        g = g / np.sum(np.abs(g)) if kernel == "rect" else g
        # Delayed signal s[n - ell] with circular indexing.
        s_delayed = np.zeros(N, dtype=complex)
        for k, gk in zip(taps, g):
            if gk == 0:
                continue
            idx = (n - ell_int - k) % N
            s_delayed += gk * s[idx]
        doppler = np.exp(1j * 2 * np.pi * fdop * n / N)   # eq. (5) phase term
        r += hp * doppler * s_delayed
    return r


# ----------------------------------------------------------------------------
# 3. Pulse-shaping / fractional-delay interpolation kernel g(.) (eq. 6)
# ----------------------------------------------------------------------------
def frac_delay_kernel(x, kernel="sinc", alpha=0.5):
    """Effective discrete-time pulse kernel g(.) evaluated at offsets x.

    kernel : "sinc" -> ideal band-limited interpolation
             "rc"   -> raised cosine with roll-off alpha
             "rect" -> nearest-sample (rectangular, no pulse shaping)
    """
    x = np.asarray(x, dtype=float)
    if kernel == "sinc":
        return np.sinc(x)
    if kernel == "rc":
        # Raised-cosine interpolation kernel.
        sinc = np.sinc(x)
        denom = 1 - (2 * alpha * x) ** 2
        cos = np.cos(np.pi * alpha * x)
        out = sinc * np.where(np.abs(denom) < 1e-8, np.pi / 4, cos / denom)
        return out
    if kernel == "rect":
        return (np.abs(x) < 0.5).astype(float)
    raise ValueError(f"unknown kernel '{kernel}'")


def build_channel_matrix(paths, N, kernel="sinc"):
    """Build the explicit N x N channel matrix H = sum_p h_p V^{f_p} G(l_p).

    This is the dense formulation of eq. (14)/(16) and is intended for SMALL N
    visualization only (reproducing the heatmaps of Figs. 1-3); for large N use
    apply_channel() instead.
    """
    n = np.arange(N)[:, None]
    m = np.arange(N)[None, :]
    H = np.zeros((N, N), dtype=complex)
    for hp, ell, fdop in zip(paths["h"], paths["ell"], paths["fdop"]):
        # Interpolated (fractional) delay matrix G(l_p): [G]_{n,m}=g(n-m-l_p),
        # with circular (mod N) support so it matches the CP/CPP convolution.
        diff = ((n - m - ell + N / 2) % N) - N / 2
        G = frac_delay_kernel(diff, kernel)
        # Diagonal Doppler phase V^{f_p} = diag(exp(j 2 pi f_p n / N)).
        V = np.exp(1j * 2 * np.pi * fdop * n.ravel() / N)
        H += hp * (V[:, None] * G)
    return H


def per_path_matrix(hp, ell, fdop, N, kernel="sinc"):
    """Single-path channel matrix h_p V^{f_p} G(l_p) for Fig. 1/2 style plots."""
    return build_channel_matrix(
        {"h": np.array([hp]), "ell": np.array([ell]), "fdop": np.array([fdop])},
        N, kernel)


# ----------------------------------------------------------------------------
# 4. Constellation mappers (Gray-coded, unit average power)
# ----------------------------------------------------------------------------
_QPSK = np.array([1 + 1j, 1 - 1j, -1 + 1j, -1 - 1j]) / np.sqrt(2)


def qpsk_map(bits):
    """Map a bit array (length multiple of 2) to QPSK symbols (unit power)."""
    b = np.asarray(bits).reshape(-1, 2)
    idx = b[:, 0] * 2 + b[:, 1]
    return _QPSK[idx]


def qpsk_demap(sym):
    """Hard-decision QPSK demapping -> bit array."""
    sym = np.asarray(sym)
    b0 = (sym.real < 0).astype(int)
    b1 = (sym.imag < 0).astype(int)
    return np.column_stack([b0, b1]).reshape(-1)


# 16-QAM Gray map over levels {-3,-1,1,3}/sqrt(10).
_LEVELS = np.array([-3, -1, 1, 3]) / np.sqrt(10)
_GRAY = {0: 0, 1: 1, 3: 2, 2: 3}            # bits(2) -> level index
_GRAY_INV = {v: k for k, v in _GRAY.items()}


def qam16_map(bits):
    """Map a bit array (length multiple of 4) to 16-QAM symbols (unit power)."""
    b = np.asarray(bits).reshape(-1, 4)
    i_idx = _vec_gray(b[:, 0] * 2 + b[:, 1])
    q_idx = _vec_gray(b[:, 2] * 2 + b[:, 3])
    return _LEVELS[i_idx] + 1j * _LEVELS[q_idx]


def qam16_demap(sym):
    """Hard-decision 16-QAM demapping -> bit array."""
    sym = np.asarray(sym)
    i_idx = np.argmin(np.abs(sym.real[:, None] - _LEVELS[None, :]), axis=1)
    q_idx = np.argmin(np.abs(sym.imag[:, None] - _LEVELS[None, :]), axis=1)
    bi = _vec_gray_inv(i_idx)
    bq = _vec_gray_inv(q_idx)
    out = np.column_stack([bi // 2, bi % 2, bq // 2, bq % 2])
    return out.reshape(-1)


def _vec_gray(v):
    return np.array([_GRAY[int(x)] for x in v])


def _vec_gray_inv(v):
    return np.array([_GRAY_INV[int(x)] for x in v])


# ----------------------------------------------------------------------------
# 5. DMRS (QPSK reference symbols) + noise + equalization
# ----------------------------------------------------------------------------
def gen_dmrs(indices, seed=12345):
    """Generate QPSK DMRS symbols for the given resource indices.

    Uses a deterministic PRBS-like sequence so TX and RX agree.
    """
    rng = np.random.default_rng(seed)
    bits = rng.integers(0, 2, size=2 * len(indices))
    return qpsk_map(bits)


def awgn(x, snr_db, sig_power=1.0):
    """Add complex AWGN at the given SNR (dB) relative to sig_power."""
    snr = 10 ** (snr_db / 10.0)
    n0 = sig_power / snr
    noise = np.sqrt(n0 / 2) * (np.random.randn(*x.shape) + 1j * np.random.randn(*x.shape))
    return x + noise, n0


def lmmse_equalize(y, Xi, sigma2):
    """LMMSE equalization x_hat = (Xi^H Xi + sigma2 I)^{-1} Xi^H y."""
    N = Xi.shape[0]
    A = Xi.conj().T @ Xi + sigma2 * np.eye(N)
    return np.linalg.solve(A, Xi.conj().T @ y)


def estimate_channel_from_dmrs(y, dmrs_sym, dmrs_idx, N, smooth=11):
    """Diagonal (per affine-frequency) channel estimate from QPSK DMRS.

    LS estimate H_hat[k] = y[k]/dmrs[k] at pilot indices, then linearly
    interpolated and lightly smoothed over the full grid. Returns the diagonal
    estimate g_hat (length N) usable as a one-tap equalizer when the effective
    channel is approximately diagonal (the AFDM regime of interest).
    """
    h_ls = y[dmrs_idx] / dmrs_sym
    full = np.arange(N)
    g_re = np.interp(full, dmrs_idx, h_ls.real)
    g_im = np.interp(full, dmrs_idx, h_ls.imag)
    g_hat = g_re + 1j * g_im
    if smooth > 1:
        w = np.ones(smooth) / smooth
        g_hat = np.convolve(np.r_[g_hat, g_hat, g_hat], w, mode="same")[N:2 * N]
    return g_hat
