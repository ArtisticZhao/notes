"""Common helpers for the NR initial-access PRACH simulation notebook.

All functions are pure numpy / matplotlib so the notebook stays
"from first principles". The module covers four areas referenced
by the paper *On the Design Details of SS/PBCH, Signal Generation
and PRACH in 5G-NR* (Chakrapani, IEEE Access 2020):

    - ZC sequence generation, cyclic shift, cyclic xcorr
    - Sample-rate OFDM modulation / demodulation with PRACH tone mapping
    - Full PRACH TX / RX chain with optional coherent or non-coherent
      combining, common-phase de-rotation (V-C, eq. 14) and data-CP
      phase-ramp compensation (V-D, eq. 15)
    - PRACH cell-dimensioning (IV-C-2 four-step procedure, Table 2)
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    "PRACH_FORMATS",
    "zc_td", "zc_fd", "cyclic_shift", "cyclic_xcorr",
    "ofdm_modulate", "ofdm_demodulate", "tone_start_index",
    "prach_tx", "prach_channel", "prach_rx_pdp",
    "undo_common_phase", "compensate_cp_phase_ramp",
    "awgn", "detect_peaks_per_zone", "estimate_noise_var",
    "cell_dimensioning",
    "plot_pdp_with_zones", "plot_prach_tf_schematic",
]


# -----------------------------------------------------------------------------
# Parameter table -- short-sequence PRACH formats at fs = 122.88 Msps, kappa=64
# Values match 38.211 Table 6.3.3.1-2 scaled by mu (see paper Section IV).
# -----------------------------------------------------------------------------
PRACH_FORMATS = {
    # name : (L_RA, n_seq repetitions, N_u (sequence samples), N_cp (CP samples)
    #         at SCS=120 kHz (mu=3); has_guard_time)
    "A1": dict(L_RA=139, n_seq=2, N_u=2048, N_cp=288,  has_gt=False),
    "A2": dict(L_RA=139, n_seq=4, N_u=2048, N_cp=576,  has_gt=False),
    "A3": dict(L_RA=139, n_seq=6, N_u=2048, N_cp=864,  has_gt=False),
    "B1": dict(L_RA=139, n_seq=2, N_u=2048, N_cp=216,  has_gt=True),
    "B4": dict(L_RA=139, n_seq=12, N_u=2048, N_cp=936, has_gt=True),
    "C0": dict(L_RA=139, n_seq=1, N_u=2048, N_cp=1240, has_gt=True),
    "C2": dict(L_RA=139, n_seq=4, N_u=2048, N_cp=2048, has_gt=True),
}


# -----------------------------------------------------------------------------
# Zadoff-Chu sequences
# -----------------------------------------------------------------------------
def zc_td(u: int, L: int = 139) -> np.ndarray:
    """Time-domain ZC sequence, NR convention x_u(n) = exp(-j*pi*u*n*(n+1)/L)."""
    n = np.arange(L)
    return np.exp(-1j * np.pi * u * n * (n + 1) / L)


def zc_fd(u: int, L: int = 139) -> np.ndarray:
    """Pre-computed FFT of zc_td, useful as a matched-filter reference."""
    return np.fft.fft(zc_td(u, L))


def cyclic_shift(x: np.ndarray, n: int) -> np.ndarray:
    return np.roll(x, n)


def cyclic_xcorr(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Normalized cyclic cross-correlation; peak at the cyclic offset of a w.r.t. b."""
    L = len(a)
    return np.fft.ifft(np.fft.fft(a) * np.conj(np.fft.fft(b))) / L


# -----------------------------------------------------------------------------
# Sample-rate OFDM with PRACH tone mapping
# -----------------------------------------------------------------------------
def tone_start_index(n_fft: int, k1: int = 0, k_bar: int = 2,
                     K: int = 1) -> int:
    """Eq. (13) of the paper, simplified for k_mu_0 = 0 (single PRACH BWP).

    For PRACH A1 with delta_f_RA == delta_f_data, K = 1, k1 places the lowest
    PRB and k_bar = 2 is the sub-carrier offset inside the PRB.
    """
    return n_fft // 2 + K * k1 + k_bar


def ofdm_modulate(symbols_fd: np.ndarray, n_fft: int,
                  k_start: int, cp_len: int = 0) -> np.ndarray:
    """Map ``symbols_fd`` onto consecutive sub-carriers, IFFT, prepend CP.

    The grid is zero everywhere except indices ``k_start .. k_start+len-1``.
    Output is a single OFDM symbol (length ``cp_len + n_fft``).
    """
    grid = np.zeros(n_fft, dtype=complex)
    n = len(symbols_fd)
    grid[k_start:k_start + n] = symbols_fd
    # DC at bin 0 convention -> ifftshift after building a centred grid.
    td = np.fft.ifft(np.fft.ifftshift(grid)) * np.sqrt(n_fft)
    if cp_len:
        td = np.concatenate([td[-cp_len:], td])
    return td


def ofdm_demodulate(rx_td: np.ndarray, n_fft: int, cp_len: int,
                    k_start: int, n_tones: int) -> np.ndarray:
    """Remove CP, FFT, extract ``n_tones`` consecutive sub-carriers."""
    sym = rx_td[cp_len:cp_len + n_fft]
    grid = np.fft.fftshift(np.fft.fft(sym)) / np.sqrt(n_fft)
    return grid[k_start:k_start + n_tones]


# -----------------------------------------------------------------------------
# Full PRACH chain
# -----------------------------------------------------------------------------
def prach_tx(u: int, v: int, *, fmt: str = "A1",
             n_fft: int = 1024, k1: int = 0, k_bar: int = 2) -> np.ndarray:
    """Generate a single UE's baseband PRACH preamble for one TD occasion.

    Structure (no guard time formats are concatenated as CP || N_rep x SEQ):

        [ aggregate CP ][ SEQ_0 ][ SEQ_1 ] ... [ SEQ_{N_rep-1} ]

    Each SEQ is the IFFT of the frequency-domain ZC sequence placed on
    ``L_RA`` consecutive sub-carriers starting at ``tone_start_index``.
    """
    p = PRACH_FORMATS[fmt]
    L = p["L_RA"]
    n_rep = p["n_seq"]
    cp_len = p["N_cp"]
    Cv = v * _ncs_for_fmt(fmt)
    # FD reference (no shift) then apply phase ramp == cyclic shift in time domain.
    x_td = cyclic_shift(zc_td(u, L), Cv)
    x_fd = np.fft.fft(x_td)
    k_start = tone_start_index(n_fft, k1=k1, k_bar=k_bar)
    one_seq = ofdm_modulate(x_fd, n_fft, k_start, cp_len=0)  # CP appended below
    body = np.tile(one_seq, n_rep)
    cp = body[-cp_len:]
    return np.concatenate([cp, body])


def prach_channel(tx_list, taus, amps, *, snr_db: float,
                  rng: np.random.Generator | None = None) -> np.ndarray:
    """Sum UE signals with sample-level integer delays, then add AWGN.

    All ``tx_list[i]`` must have the same length. ``taus`` is in samples
    (>=0). Noise power is referenced to the *summed* signal mean power.
    """
    rng = rng or np.random.default_rng()
    n = len(tx_list[0])
    rx = np.zeros(n, dtype=complex)
    for s, tau, amp in zip(tx_list, taus, amps):
        rx += amp * np.roll(s, int(tau))
    return awgn(rx, snr_db, rng=rng)


def prach_rx_pdp(rx_td: np.ndarray, u_root: int, *, fmt: str = "A1",
                 n_fft: int = 1024, k1: int = 0, k_bar: int = 2,
                 n_ifft: int | None = None, combine: str = "noncoh"
                 ) -> np.ndarray:
    """Full PRACH receiver: CP strip -> FFT -> tone extract -> matched filter
    -> per-rep N_IFFT IFFT -> combine.

    Parameters
    ----------
    combine : 'noncoh' or 'coh'
        Non-coherent power sum, or coherent FD-sum then |IFFT|^2.
    n_ifft : optional zero-padded IFFT size; controls timing resolution
        (>= L_RA, default L_RA).
    """
    p = PRACH_FORMATS[fmt]
    L = p["L_RA"]
    n_rep = p["n_seq"]
    cp_len = p["N_cp"]
    n_ifft = n_ifft or L
    k_start = tone_start_index(n_fft, k1=k1, k_bar=k_bar)
    root_fd = zc_fd(u_root, L)

    # Strip aggregate CP, then for each SEQ do one N_FFT FFT.
    body = rx_td[cp_len:cp_len + n_rep * n_fft]
    seq_tones = np.empty((n_rep, L), dtype=complex)
    for r in range(n_rep):
        seg = body[r * n_fft:(r + 1) * n_fft]
        seq_tones[r] = ofdm_demodulate(np.concatenate([np.zeros(0), seg]),
                                       n_fft, cp_len=0,
                                       k_start=k_start, n_tones=L)

    corr_fd = seq_tones * np.conj(root_fd)  # FD matched filter, per rep
    # Zero-pad to n_ifft on each rep, then IFFT
    pad = np.zeros((n_rep, n_ifft), dtype=complex)
    pad[:, :L] = corr_fd
    td = np.fft.ifft(pad, axis=1)
    if combine == "noncoh":
        pdp = np.sum(np.abs(td) ** 2, axis=0)
    elif combine == "coh":
        pdp = np.abs(np.sum(td, axis=0)) ** 2
    else:
        raise ValueError(f"unknown combine={combine!r}")
    # Match the original notebook's amplitude scaling (per L^2) so peaks
    # land in a comparable range across L and N_IFFT.
    return pdp / (L ** 2)


# -----------------------------------------------------------------------------
# Receiver corrections (Section V-C and V-D of the paper)
# -----------------------------------------------------------------------------
def undo_common_phase(seq_tones: np.ndarray, f0_hz: float, mu: int,
                      symbol_indices) -> np.ndarray:
    """Remove the per-symbol common phase correction applied in the data DFE.

    For PRACH symbol ``l`` in the slot, the data-chain phase is
    ``exp(-i*2*pi*f0 * (t_start_l + T_CP_l))`` (eq. 14). Coherent combining
    requires undoing this rotation. The closed-form per-symbol value
    depends on the symbol layout in the slot; we accept a list of accumulated
    phases ``symbol_indices`` (already in seconds * f0 * 2*pi) for clarity.
    Returns ``seq_tones`` rotated symbol-by-symbol so they share a phase
    reference.
    """
    out = np.empty_like(seq_tones)
    for r, theta in enumerate(symbol_indices):
        out[r] = seq_tones[r] * np.exp(1j * 2 * np.pi * f0_hz * theta)
    return out


def compensate_cp_phase_ramp(seq_tones: np.ndarray, *, n_fft: int, k2: int,
                             N_RA_CP: int, N_CP_data: int) -> np.ndarray:
    """Apply eq. (15) phase-ramp correction per OFDM symbol of PRACH.

    When a common wideband FFT is used and the FFT window is aligned to the
    data CP (not the PRACH aggregate CP), each PRACH symbol's tones see a
    linear phase ramp rho_j across the L_RA sub-carriers. This function
    multiplies it out.
    """
    n_rep, L = seq_tones.shape
    out = np.empty_like(seq_tones)
    for j in range(n_rep):
        rho = (N_RA_CP - N_CP_data - j * N_CP_data) % n_fft
        n = np.arange(L)
        out[j] = seq_tones[j] * np.exp(1j * 2 * np.pi * (n + k2) * rho / n_fft)
    return out


# -----------------------------------------------------------------------------
# Noise / detection
# -----------------------------------------------------------------------------
def awgn(signal: np.ndarray, snr_db: float,
         rng: np.random.Generator | None = None) -> np.ndarray:
    rng = rng or np.random.default_rng()
    p_sig = np.mean(np.abs(signal) ** 2)
    p_n = p_sig / (10 ** (snr_db / 10))
    n = np.sqrt(p_n / 2) * (rng.standard_normal(signal.shape)
                            + 1j * rng.standard_normal(signal.shape))
    return signal + n


def estimate_noise_var(pdp: np.ndarray, n_cs: int) -> float:
    """Estimate noise floor from PDP samples that lie at the **end** of each
    ZCZ zone (least likely to contain a peak). Returns a scalar variance
    that the detector can use as a normaliser."""
    L = len(pdp)
    n_zones = L // n_cs
    # Take the last sample of each zone (and the leftover tail) as background.
    tail_idx = [z * n_cs + n_cs - 1 for z in range(n_zones)]
    background = pdp[tail_idx]
    return float(np.median(background))


def detect_peaks_per_zone(pdp: np.ndarray, n_cs: int,
                          *, threshold_ratio: float = 20.0,
                          noise_var: float | None = None) -> list[dict]:
    """For every ZCZ zone of length ``n_cs`` return a record dict.

    'detected' is True when ``peak / noise_var > threshold_ratio``. If
    ``noise_var`` is None, falls back to ``estimate_noise_var``.
    """
    if noise_var is None:
        noise_var = estimate_noise_var(pdp, n_cs)
    L = len(pdp)
    n_zones = L // n_cs
    out = []
    for z in range(n_zones):
        seg = pdp[z * n_cs:(z + 1) * n_cs]
        idx_in_zone = int(np.argmax(seg))
        peak = float(seg[idx_in_zone])
        out.append(dict(
            v=z,
            peak_idx=z * n_cs + idx_in_zone,
            tau_est=idx_in_zone,
            peak=peak,
            snr_db=10 * np.log10(peak / noise_var) if noise_var > 0 else np.inf,
            detected=(peak / noise_var) > threshold_ratio,
        ))
    return out


# -----------------------------------------------------------------------------
# Cell dimensioning (paper IV-C-2 four-step procedure -> Table 2)
# -----------------------------------------------------------------------------
# zeroCorrelationZoneConfig -> N'_CS for FR2 short formats (Table 6.3.3.1-7 of 38.211).
_N_CS_PRIME_FR2 = {
    "A1": 288, "A2": 576, "A3": 864,
    "B1": 216, "B4": 936,
    "C0": 1240, "C2": 2048,
}
# Quantization grid (subset of Table 6.3.3.1-7); finer table omitted for clarity.
_NCS_QUANTIZED_VALUES_FR2 = [0, 2, 4, 6, 8, 10, 12, 13, 15, 17, 19, 23, 27, 34, 46, 69]


def _quantize_ncs(n_cs_prime: int) -> int:
    """Pick the largest tabled NCS <= floor(N'_CS)."""
    candidates = [v for v in _NCS_QUANTIZED_VALUES_FR2 if v <= int(np.floor(n_cs_prime))]
    return max(candidates) if candidates else 0


def _ncs_for_fmt(fmt: str) -> int:
    p = PRACH_FORMATS[fmt]
    raw = _N_CS_PRIME_FR2[fmt] / p["N_u"] * p["L_RA"]
    return _quantize_ncs(raw)


def cell_dimensioning(fmt: str, mu: int = 3, *, n_preambles: int = 64,
                      tau_d_us: float = 4.69) -> dict:
    """4-step procedure from paper Section IV-C-2 (FR2, unrestricted set).

    Returns N_CS, cyclic-shifts/root, root sequences needed, max cell radius (m).
    """
    p = PRACH_FORMATS[fmt]
    L = p["L_RA"]
    n_cs_prime = _N_CS_PRIME_FR2[fmt] / p["N_u"] * L
    n_cs = _quantize_ncs(n_cs_prime)
    cv = L // n_cs if n_cs else 0
    roots = int(np.ceil(n_preambles / cv)) if cv else np.inf
    delta_f_RA = 120e3 * 2 ** (mu - 3)
    # Paper eq.: radius = (N_CS / (delta_f_RA * L_RA) - tau_d / 2^mu) * c / 2
    radius_m = (n_cs / (delta_f_RA * L) - tau_d_us * 1e-6 / 2 ** mu) * 3e8 / 2
    return dict(format=fmt, n_cs_prime=n_cs_prime, n_cs=n_cs,
                cv_per_root=cv, roots_needed=roots,
                cell_radius_m=max(0.0, radius_m))


# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------
def plot_pdp_with_zones(pdp: np.ndarray, n_cs: int, ues: list[dict] | None = None,
                        ax=None, title: str | None = None, normalize: bool = True):
    """Plot a PDP with zero-correlation-zone shading and optional UE markers."""
    if ax is None:
        _, ax = plt.subplots(figsize=(13, 4.5))
    y = pdp / pdp.max() if normalize else pdp
    ax.plot(y, "b-", linewidth=1.2)
    L = len(pdp)
    n_zones = L // n_cs
    used_v = {u["v"] for u in (ues or [])}
    for z in range(n_zones):
        color = "lightgreen" if z in used_v else "lightgray"
        ax.axvspan(z * n_cs, (z + 1) * n_cs, color=color, alpha=0.25)
        ax.text(z * n_cs + n_cs / 2, 1.07, f"v={z}",
                ha="center", fontsize=8)
    for ue in ues or []:
        expected = (ue["v"] * n_cs + ue["tau"]) % L
        ax.axvline(expected, color="red", linestyle=":", alpha=0.6)
        ax.annotate(f"{ue['name']}\nv={ue['v']}, tau={ue['tau']}",
                    xy=(expected, 0.92), xytext=(expected, 0.58),
                    ha="center", fontsize=9,
                    bbox=dict(boxstyle="round", facecolor="yellow", alpha=0.85),
                    arrowprops=dict(arrowstyle="->", color="red"))
    ax.set_xlim(0, L); ax.set_ylim(0, 1.2)
    ax.set_xlabel("PDP index"); ax.set_ylabel("Normalized power")
    if title:
        ax.set_title(title)
    return ax


def plot_prach_tf_schematic(fmt: str = "A1", *, m_fd: int = 2, n_td: int = 6,
                            slot_us: float = 125.0, ax=None):
    """Schematic of multiple TD x FD occasions inside one PRACH slot."""
    if ax is None:
        _, ax = plt.subplots(figsize=(13, 4.2))
    p = PRACH_FORMATS[fmt]
    L = p["L_RA"]; n_rep = p["n_seq"]
    delta_f_RA = 120e3  # FR2 mu=3
    seq_us = 1e6 / delta_f_RA
    cp_us = p["N_cp"] / (delta_f_RA * p["N_u"]) * 1e6
    total_per_occ = cp_us + n_rep * seq_us
    prach_bw_mhz = L * delta_f_RA / 1e6
    fd_y = [i * (prach_bw_mhz + 1.0) for i in range(m_fd)]
    for fd_idx, y0 in enumerate(fd_y):
        t = 1.0
        for i in range(n_td):
            ax.add_patch(plt.Rectangle((t, y0), cp_us, prach_bw_mhz,
                                       facecolor="orange",
                                       label="CP" if (i == 0 and fd_idx == 0) else None,
                                       edgecolor="black", linewidth=0.4))
            for r in range(n_rep):
                ax.add_patch(plt.Rectangle((t + cp_us + r * seq_us, y0),
                                           seq_us, prach_bw_mhz,
                                           facecolor="steelblue",
                                           label="SEQ" if (i == 0 and r == 0 and fd_idx == 0) else None,
                                           edgecolor="black", linewidth=0.4))
            ax.text(t + total_per_occ / 2, y0 + prach_bw_mhz / 2,
                    f"TD#{i}", ha="center", va="center",
                    fontsize=9, color="white")
            t += total_per_occ
        ax.text(-2, y0 + prach_bw_mhz / 2, f"FD#{fd_idx}",
                ha="right", va="center", fontsize=10, fontweight="bold")
    ax.set_xlim(-3, slot_us + 2)
    ax.set_ylim(-1, fd_y[-1] + prach_bw_mhz + 1)
    ax.set_xlabel("time within slot (us)")
    ax.set_ylabel("frequency (MHz, schematic)")
    ax.set_title(f"PRACH {fmt} in {slot_us:.0f}us slot: "
                 f"{n_td} TD x {m_fd} FD; "
                 f"each TD = CP + {n_rep}xSEQ; "
                 f"PRACH BW/FD = {prach_bw_mhz:.2f} MHz")
    ax.legend(loc="upper right")
    return ax


# -----------------------------------------------------------------------------
# Self test
# -----------------------------------------------------------------------------
def _selftest():
    L = 139
    x = zc_td(1, L)
    assert np.allclose(np.abs(x), 1.0)

    # OFDM roundtrip
    n_fft = 1024
    fd = np.fft.fft(x)
    k0 = tone_start_index(n_fft)
    sym = ofdm_modulate(fd, n_fft, k0, cp_len=128)
    rec = ofdm_demodulate(sym, n_fft, cp_len=128, k_start=k0, n_tones=L)
    assert np.max(np.abs(rec - fd)) < 1e-9, f"roundtrip err={np.max(np.abs(rec - fd))}"

    # PRACH chain: a single UE, no noise, no delay -> peak at index 0
    tx = prach_tx(u=1, v=0, fmt="A1", n_fft=n_fft)
    rx = prach_channel([tx], taus=[0], amps=[1.0], snr_db=60)
    pdp = prach_rx_pdp(rx, u_root=1, fmt="A1", n_fft=n_fft)
    assert int(np.argmax(pdp)) == 0

    # cell dimensioning A1 numbers
    d = cell_dimensioning("A1")
    assert d["n_cs"] == 19
    assert d["cv_per_root"] == 7
    assert d["roots_needed"] == 10
    print("nr_init_access self-test ok:", d)


if __name__ == "__main__":
    _selftest()
