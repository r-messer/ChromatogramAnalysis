import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import time
import os

cwd = os.getcwd()

''' Assemlby #16 '''
path1 = '/Users/ryan/Desktop/PFV Grant Renewal/MMTV WT#5 unlabeled 061819 001.xlsm'


''' Assemlby #19 '''
path2 = '/Users/ryan/Desktop/PFV Grant Renewal/01232020 PFV #93 WT 616-Cy3_675-ddA 001.xlsm'


t0 = time.perf_counter()

dataframe1 = pd.read_excel(path1, header=None)
dataframe2 = pd.read_excel(path2, header=None)

x_reg = list(dataframe1[10][3:97])
y_reg = list(dataframe1[11][3:97])

frac_m, frac_b = np.polyfit(x_reg, y_reg, 1)

masterlist = [dataframe1]


# Initialize Parellel Processing


def reindex_chrom(df):
    """ Reindexes top 3 rows of data frame to the format *signal* + *axis label*,
     deletes rows that do not contain numbers (top 3 rows) """
    new_index1 = []
    new_index2 = []

    for i in df.iloc[1]:
        if isinstance(i, str):
            try:
                new_index1.append(i.split('_', 1)[1])
            except IndexError:
                new_index1.append(i.split()[0])

    for i in df.iloc[2]:
        if isinstance(i, str):
            new_index2.append(i.split()[0])

    final_index = []

    x = 0
    for i in new_index1:
        """ This is an interesting iteration, what would it look like recursively? """
        y = 0

        while y < 2:
            final_index.append(new_index2[x] + i)
            x += 1
            y += 1

    df.columns = final_index

    df.drop([0, 1, 2], inplace=True)
    df.reset_index(drop=True, inplace=True)
    print(df.columns)
    return df


def frac_reg(x):
    """ Linear equation for mLFraction, FractionFraction regression """
    return (frac_m * x) + frac_b


def inv_frac_reg(x):
    """ Invers of frac_reg """
    return (x / frac_m) - frac_b


def plot_chrom(df, path):
    """ Plots data based on indices in dataframe. Will only plot if index starts with 'ml' and
     is < 6. """

    fig, ax = plt.subplots(constrained_layout=True, figsize=size)

    # Default plt settings
    ax.set_prop_cycle(mpl.cycler(color=['r', 'b', 'm', 'c']))
    mpl.rc('lines', linewidth=4)

    plt.xlim([-0.007904, 25.998215])

    # Title
    ax.set_title(path.split('/')[-1], size=18, weight='bold')

    ax.set_xlabel('Retention Volume', size=12)
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.set_ylabel('mAU', size=12)

    ax.tick_params('both', labelsize=12, width=3, length=5)

    ax.spines["top"].set_linewidth(3)
    ax.spines["left"].set_linewidth(3)
    ax.spines["right"].set_linewidth(3)

    secax = ax.secondary_xaxis('bottom', functions=(frac_reg, inv_frac_reg))
    secax.set_xlabel('Fraction', size=12)
    secax.tick_params('both', labelsize=12, width=3, length=4)

    secax.spines["bottom"].set_linewidth(3)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, fontsize=12, bbox_to_anchor=(1.01, 1.0), loc='upper left', borderaxespad=0.)

    for x in df.columns:
        if x[1] == 'l' and len(x) < 6:
            try:
                label = x.split('l')[1]
                ax.plot(df[x], df['mAU'+label], label=label)
            except KeyError as c:
                print('{} is unable to be plotted'.format(c))


# Plotting Functions
fact = 1.8
size = [6.4*fact, 4.8*fact]

# Call reindex_chrom
for df in masterlist:
    reindex_chrom(df)

# Plot reindexed chroms


plot_chrom(dataframe1, path1)
plot_chrom(dataframe2, path2)


t1 = time.perf_counter()

print(f'\nCompleted in: {t1-t0} seconds')


plt.show()



