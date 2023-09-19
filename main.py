from singlefig import *
from findpeaks import *
from tkinter import Tk
from tkinter.filedialog import askopenfilename

Tk().withdraw()
path1 = askopenfilename()
#path1 = '/Volumes/RYANSUSB/Chromatograms/2022/05252022 No 1 PARP+PAR Only Control Superdex 200 inc 10-300 PAR Histone Separation 001.csv'

t0 = time.perf_counter()

C = Chrom(path1)
print(C.path.split('/')[-1])

fig, ax = C.plot_chrom()

find_chrom_peaks(C)

plot_peaks(C)

# C.df.to_csv(path_or_buf='~/Desktop/chromcheck.csv')

t1 = time.perf_counter()

plt.show()

print(f'Finished in:  {t1-t0} second(s)')