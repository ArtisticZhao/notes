from matplotlib import pyplot as plt
from p2ch10.vis import findPositiveSamples, showCandidate
positiveSample_list = findPositiveSamples()

series_uid = positiveSample_list[11][2]
showCandidate(series_uid)

plt.show()