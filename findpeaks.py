import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_prominences
from scipy.stats import sem

from singlefig import Chrom


def find_chrom_peaks(chrom: Chrom, output='indices'):
    """ Finds peaks from a chromatogram using the mAU values, returns dictionary of peaks"""
    var_dict = {}
    peaks_dict = {}

    for x in chrom.df.columns:

        if 'mAU' in x and len(x) < 7:
            ''' Should put try/except here '''

            label = x + 'peaks'
            avg_x = np.convolve(chrom.df[x], np.ones(3), 'valid') / 3
            peaks, _ = find_peaks(avg_x, height=(0, 3000), prominence=1.5)
            chrom.df[label] = pd.Series(peaks)


def plot_peaks(chrom: Chrom):
    for a in chrom.df.columns:
        if 'peaks' in a:
            ''' Should make these aliases for speed and coherence '''
            mL_col = 'ml' + a.rstrip('peaks').lstrip('mAU')
            mAU_col = a.rstrip('peaks')

            x_vals = []
            y_vals = []

            for b in chrom.df[a]:
                if np.isnan(b):
                    break
                else:
                    x_vals.append(chrom.df[mL_col][b])
                    y_vals.append(chrom.df[mAU_col][b])

            chrom.Axes.plot(x_vals, y_vals, marker='D', linestyle='', markersize='8', color='k')

            for c, txt in enumerate(x_vals):
                an = chrom.Axes.annotate(np.around(txt, decimals=2), (x_vals[c], y_vals[c]),
                                         xycoords='data', xytext=(x_vals[c], y_vals[c]+np.log(np.max(y_vals))),
                                         fontweight='bold')
                an.draggable()

    # for var in df.columns:
    #     for key in keys:
    #         if key in var:
    #             var_dict[var] = np.array(df[var])
    #     if 'mAU' in var:
    #         var_dict[var + '_movmean'] = np.convolve(var_dict[var], np.ones((50,)) / 50)
    #
    # for signal in var_dict.keys():
    #     if '_movmean' in signal:
    #         peaks, _ = find_peaks(var_dict[signal], prominence=1)
    #         peaks_dict[signal.split('_')[0]+'_peaks'] = peaks
    #
    # if output == 'values':
    #     for var in peaks_dict.keys():
    #         peaks_dict[var] = var_dict['ml'+var.split('U')[1].split('_')[0]][peaks_dict[var]]
    #     return peaks_dict
    # else:
    #     return peaks_dict


# def compare_peaks(value, array2: np.array):
#
#     len_arr1 = len(array2)
#
#     if sem([value, array2[0]])/np.mean([value, array2[0]]) < 0.05:
#         return np.mean([value, array2[0]])
#     elif value < array2[0] - (0.05*array2[0]):
#         print('Value << Array')
#         return None
#     else:
#         try:
#             return compare_peaks(value, array2[1:])
#         except IndexError:
#             print('No values found')
#             return None


# chrom = Chrom(os.getcwd()+'/TestDataFiles/07222020 HIV Assembly #21 (WTSso7d+p75+U5-25_9Cy5) 001 copy.xls')
#
# peaks = find_chrom_peaks(chrom, 'values')
# print(peaks)
# plt.plot(chrom.df['ml280'], chrom.df['mAU280'], color='blue')
# plt.plot(chrom.df['ml254'], chrom.df['mAU254'], color='red')
# plt.plot(chrom.df['ml646'], chrom.df['mAU646'], color='m')
# plt.plot(chrom.df['ml280'][peaks['mAU280_peaks']], chrom.df['mAU280'][peaks['mAU280_peaks']], 'x', color='blue')
# plt.plot(chrom.df['ml254'][peaks['mAU254_peaks']], chrom.df['mAU254'][peaks['mAU254_peaks']], 'x', color='red')
# plt.plot(chrom.df['ml646'][peaks['mAU646_peaks']], chrom.df['mAU646'][peaks['mAU646_peaks']], 'x', color='m')
plt.show()

# figure2 = plt.figure()
#
# peaks, _ = find_peaks(x, prominence=1)
#
# prominences = peak_prominences(x, peaks)[0]
#
# contour_heights = x[peaks] - prominences
#
# plt.plot(x)
# plt.plot(peaks, x[peaks], "x")
# plt.plot(np.zeros_like(x), "--", color="gray")
# plt.vlines(x=peaks, ymin=contour_heights, ymax=x[peaks])
#
# plt.show()
