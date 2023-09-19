from abc import ABC
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import os
import time



# Initialize Parellel Processing


class Chrom(pd.DataFrame, ABC):
    """ Class Chrom inherits from pd.Dataframe, stores chromatogram data which can be
     reindexed, plotted. Takes variable 'path' as str where an excel doc should be located"""

    def __init__(self, path: str):
        super().__init__()
        self.path = path

        self.df = pd.read_excel(path, header=None)

        self.reindex_chrom()
        self.Figure = plt.Figure
        self.Axes = plt.Axes


    def reindex_chrom(self):
        """ Reindexes top 3 rows of data frame to the format *signal* + *axis label*,
         deletes rows that do not contain numbers (top 3 rows) """
        new_index1 = []
        new_index2 = []

        for i in self.df.iloc[1]:
            if isinstance(i, str):
                try:
                    new_index1.append(i.split('_')[-1])
                except IndexError:
                    new_index1.append(i.split()[0])

        for i in self.df.iloc[2]:
            if isinstance(i, str):
                n = i.split()[0]
                if 'Frac' in n:
                    new_index2.append('Fraction')
                else:
                    new_index2.append(n)
        print("New Index1:  ", new_index1, '\n', "New Index2:  ", new_index2)

        final_index = []

        x = 0
        for i in new_index1:
            y = 0
            if 'Frac' in i:
                i = 'Fraction'
            elif 'nm' in i:
                i = i.rstrip('nm')


            while y < 2:
                final_index.append(new_index2[x] + i)
                x += 1
                y += 1
        print("Final Index:  ", final_index)
        self.df.columns = final_index

        self.df.drop([0, 1, 2], inplace=True)
        self.df.reset_index(drop=True, inplace=True)

        return self

    def plot_chrom(self, size: int = 1.8, graphstyle='full'):
        """ Plots data based on indices in dataframe.Returns objects (figure, axes) """
        ''' Will only plot if index starts with 'ml' and
         is < 6. '''

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
            """ Inverse of frac_reg """
            return (x / frac_m) - frac_b

        self.Figure, self.Axes = plt.subplots(constrained_layout=True, figsize=[6.4 * size, 4.8 * size])


        # Default plt settings
        self.Axes.set_prop_cycle(mpl.cycler(color=['r', 'b', 'm', 'c']))
        mpl.rc('lines', linewidth=4)

        # Title
        # title = self.path.split('/')[-1].split('.')[0]
        # self.Axes.set_title(title, size=18, weight='bold')

        # Axis Labels
        self.Axes.set_xlabel('Retention Volume (mL)', size=12)
        self.Axes.xaxis.set_label_position('top')
        self.Axes.xaxis.tick_top()
        self.Axes.set_ylabel('mAU', size=12)

        self.Axes.tick_params('both', labelsize=12, width=3, length=5)
        # ax.set_xticks(np.arange(0,x_reg[-1], step=int(len(x_reg)/12))) auto ticks are better than mine,
        # maybe revisit later. Would need to update renderer when zooming in on graph - annoying

        self.Axes.spines["top"].set_linewidth(3)
        self.Axes.spines["left"].set_linewidth(3)
        self.Axes.spines["right"].set_linewidth(3)

        plt.xlim([-0.1, 24])
        """ Crops out equilibration after run """
        # ylim = plt.ylim([-10, np.array(self.df['mAU276']).max()])

        # Secondary Axis
        secax = self.Axes.secondary_xaxis('bottom', functions=(frac_reg, inv_frac_reg))
        secax.set_xlabel('Fractions', size=12)
        secax.tick_params('both', labelsize=10, width=3, length=4, labelrotation=70)
        secax.set_xticks(np.arange(0, len(x_reg)))
        secax.set_xticklabels(self.df['FractionFraction'][0:len(x_reg)])

        #Set Bottom Labels to every ther
        for label in secax.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        """ Replaces numbers with Fraction strings from FractionFraction """

        secax.spines["bottom"].set_linewidth(3)

        # Baseline
        hline = plt.hlines(0, -1, np.max(self.df['mlFraction']), linestyles=':')
        hline.set_linewidth(2)
        # Plot Data
        for x in self.df.columns:
            if x[1] == 'l' and len(x) < 6:
                """ Excludes series whose titles do not have an 'l' AND
                have less than 6 characters
                i.e. mAU280 will not be plotted as x values """

                try:
                    label = x.split('l')[1]
                    self.Axes.plot(self.df[x], self.df['mAU' + label], label=label)
                except KeyError as c:
                    print(f'{c} is unable to be plotted')

        # Legend
        handles, labels = self.Axes.get_legend_handles_labels()
        self.Axes.legend(handles, labels, fontsize=12)
        plt.grid(1,axis='y')
        return self.Figure, self.Axes
