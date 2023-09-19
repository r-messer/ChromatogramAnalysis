from abc import ABC
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import os
import time
from tkinter import Tk
from tkinter.filedialog import askopenfilename


# Initialize Parellel Processing


class Chrom(pd.DataFrame, ABC):
    """ Class Chrom inherits from pd.Dataframe, stores chromatogram data which can be
     reindexed, plotted. Takes variable 'path' as str where an excel doc should be located"""

    def __init__(self, path: str):
        super().__init__()
        self.path = path
        self.df = pd.read_excel(path, header=None)
        self.reindex_chrom()

    def reindex_chrom(self):
        """ Reindexes top 3 rows of data frame to the format *signal* + *axis label*,
         deletes rows that do not contain numbers (top 3 rows) """
        new_index1 = []
        new_index2 = []

        for i in self.df.iloc[1]:
            if isinstance(i, str):
                try:
                    new_index1.append(i.split('_', 1)[1])
                except IndexError:
                    new_index1.append(i.split()[0])

        for i in self.df.iloc[2]:
            if isinstance(i, str):
                new_index2.append(i.split()[0])

        final_index = []

        x = 0
        for i in new_index1:
            y = 0
            while y < 2:
                final_index.append(new_index2[x] + i)
                x += 1
                y += 1

        self.df.columns = final_index

        self.df.drop([0, 1, 2], inplace=True)
        self.df.reset_index(drop=True, inplace=True)

        return self.df

    def plot_chrom(self, size: int = 1.8, graphstyle='full', concB: bool = False):
        """ Plots data based on indices in dataframe. Will only plot if index starts with 'ml' and
         is < 6. """

        x_reg = list([x for x in self.df['mlFraction'] if not np.isnan(x)])
        """ Generate a list of x values from mLs of the Fraction """

        y_reg = list(range(1, len(x_reg) + 1))
        """Evenly spaced integers from 1 to length of mLFraction; Represents evenly spaced fractions of total 
        retention """

        frac_m, frac_b = np.polyfit(x_reg, y_reg, 1)

        def frac_reg(x):
            """ Linear equation for mlFraction, FractionFraction regression """
            return (frac_m * x) + frac_b

        def inv_frac_reg(x):
            """ Invers of frac_reg """
            return (x / frac_m) - frac_b

        fig, ax = plt.subplots(constrained_layout=True, figsize=[6.4 * size, 4.8 * size])

        # Default plt settings
        ax.set_prop_cycle(mpl.cycler(color=['r', 'b', 'm', 'c']))
        mpl.rc('lines', linewidth=4)

        # Title
        title = self.path.split('/')[-1]
        ax.set_title(title, size=18, weight='bold')

        # Axis Labels
        ax.set_xlabel('Retention Volume (mL)', size=12)
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()
        ax.set_ylabel('mAU', size=12)

        ax.tick_params('both', labelsize=12, width=3, length=5)
        # ax.set_xticks(np.arange(0,x_reg[-1], step=int(len(x_reg)/12))) auto ticks are better than mine,
        # maybe revisit later. Would need to update renderer when zooming in on graph - annoying

        ax.spines["top"].set_linewidth(3)
        ax.spines["left"].set_linewidth(3)
        ax.spines["right"].set_linewidth(3)

        xlim = plt.xlim([-0.1, x_reg[-1]])
        """ Crops out equilibration after run """
        # ylim = plt.ylim([-10, np.array(self.df['mAU276']).max()])

        # Secondary Axis
        secax = ax.secondary_xaxis('bottom', functions=(frac_reg, inv_frac_reg))
        secax.set_xlabel('Fraction', size=12)
        secax.tick_params('both', labelsize=10, width=3, length=4, labelrotation=70)
        secax.set_xticks(np.arange(0, len(x_reg)))
        secax.set_xticklabels(self.df['FractionFraction'][0:len(x_reg)])
        """ Replaces numbers with Fraction strings from FractionFraction """

        secax.spines["bottom"].set_linewidth(3)

        # Baseline
        hline = plt.hlines(0, -1, 30, linestyles=':')
        hline.set_linewidth(2)
        # Plot Data
        for x in self.df.columns:
            if x[1] == 'l' and len(x) < 6:
                """ Excludes series whose titles do not have an 'l' AND
                have less than 6 characters
                i.e. mAU280 will not be plotted as x values """

                try:
                    label = x.split('l')[1]
                    ax.plot(self.df[x], self.df['mAU' + label], label=label)
                except KeyError as c:
                    print(f'{c} is unable to be plotted')
            # elif x[0] == '%' and concB is True:
            #     """ Plot Concentration %B Curve if concB is True """
            #     ax.plot(self.df[x], self.df)




        # Legend
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, fontsize=12)
