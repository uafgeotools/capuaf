#-----------------------------------------------------------
# 
# script to plot scatter and histograms for comparisons bet vel models and bet
# inversions
# 
#-----------------------------------------------------------
# 

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, MultipleLocator
import numpy as np
import sys

# verify and read from inputfile
if len(sys.argv)!=2:
    sys.stderr.write('usage: %s datafile.txt\n' % sys.argv[0])
    sys.exit('%s: datafile = AMP summary from CAP out file' % sys.argv[0])

inputfile = str(sys.argv[1])
try:
    data = np.genfromtxt(inputfile, dtype=None, names = ['theta', 'dvr', 'dummy'])
except:
    sys.exit('\n%s: abort. error reading file: %s\n' % (sys.argv[0], inputfile))

sys.stderr.write("%s: processing file: %s\n" % (sys.argv[0], inputfile))

# put into arrays
data_x = data['theta']
data_y = data['dvr']

# axes dimensions
#left, width = 0.1, 0.65
left, width = 0.1, 0.35
bottom, height = 0.1, 0.35
bottom_h = left_h = left+width+0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.12]
rect_histy = [left_h, bottom, 0.12, height]

# rectangular Figure
plt.figure(1, figsize=(6, 6))

# handles 
axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)

# histograms, set bar ranges
xlim0, xlimf, dbinx = -5, 70, 5      # all figures
ylim0, ylimf, dbiny = -14, 1, 1     # NP vs P (0d, 1d)  # Figure 9, S2a
#ylim0, ylimf, dbiny = -14, 13, 1    # NP vs P (0d, 1d)  # Figure S2b
#ylim0, ylimf, dbiny = -12.5, 1, 1   # NP vs P (0d, 1d)
#ylim0, ylimf, dbiny = -25, 25, 5    # 0d vs 1d (NP, P)
#ylim0, ylimf, dbiny = -10.5, 2.5, 1 # 0d vs 1d (NP, P)
#ylim0, ylimf = -30, 30, 5

binset_x = np.arange(xlim0, xlimf, dbinx)   # top hist
binset_y = np.arange(ylim0, ylimf, dbiny)   # side hist -- NOTE -1 OFFSET

#binset_x = np.linspace(xlim0, xlimf, 15)     # top hist
#binset_y = np.linspace(ylim0, ylimf, 11)     # side hist

hist_x, bin_edges_x = np.histogram(data_x, bins=binset_x)   # top hist
hist_y, bin_edges_y = np.histogram(data_y, bins=binset_y)   # side hist

#print binset_x
#print bin_edges_x
#print binset_y
print bin_edges_y
print hist_y

# plot
color_bars=[0.3, 1.0, 0.3]

# scatter
axScatter.scatter(data_x, data_y, facecolor=color_bars, s=30)

# histograms
axHistx.bar(bin_edges_x[:-1], hist_x, width=dbinx, color=color_bars)
axHisty.barh(bin_edges_y[:-1], hist_y, height=dbiny, color=color_bars)

# no labels for histograms
nullfmt   = NullFormatter() 
axScatter.set_xlabel(r'$\theta$')

# Figure 9, S2a
# NP vs P, utuhalf theta_dvr_eid_utuhalf
# NP vs P. utu1D theta_dvr_eid_utu1d
axScatter.set_ylabel(r'$VR_{\mathrm{np}}\/-\/VR_{\mathrm{p}}$')

# utuhalf vs utu1d, both with polarities, file theta_dvr_eid_utuhalf_utu1d_pol
# Figure S2b
#axScatter.set_ylabel(r'$VR_{\mathrm{half}}\/-\/VR_{\mathrm{1d}}$')

#axScatter.axhline(y=0)

# ticks
minorLocator_scatter_x = MultipleLocator(5)
majorLocator_scatter_x = MultipleLocator(20)
minorLocator_scatter_y = MultipleLocator(1) # y-ticks scatter + side hist 1 or 5
majorLocator_scatter_y = MultipleLocator(5) # y-labels scatter. 10 or 5

minorLocator_histo_x = MultipleLocator(10)  # horiz ticks side histo
majorLocator_histo_x = MultipleLocator(50)  # 50 vs 30
minorLocator_histo_y = MultipleLocator(10)  # y-ticks top histo
majorLocator_histo_y = MultipleLocator(50)  # 50 vs 30

axHistx.xaxis.set_major_formatter(nullfmt)
axHistx.yaxis.set_major_locator(majorLocator_histo_y)   # y-labels top hist
axHistx.yaxis.set_minor_locator(minorLocator_histo_y)   # y-ticks top hist

axHisty.xaxis.set_major_locator(majorLocator_histo_x)   # x-labels side hist
axHisty.xaxis.set_minor_locator(minorLocator_histo_x)   # x-ticks side hist
axHisty.yaxis.set_major_formatter(nullfmt)              # y-labels side hist
axHisty.yaxis.set_minor_locator(minorLocator_scatter_y) # y-ticks side hist
axHisty.yaxis.set_major_locator(majorLocator_scatter_y) # y-ticks side hist

axScatter.xaxis.set_minor_locator(minorLocator_scatter_x) # x-ticks scatter
axScatter.xaxis.set_major_locator(majorLocator_scatter_x) # y-labels scatter
axScatter.yaxis.set_minor_locator(minorLocator_scatter_y) # x-ticks scatter
axScatter.yaxis.set_major_locator(majorLocator_scatter_y) # y-labels scatter

# histogram ranges
# 0d vs 1d   == (0,30)
# NP vs P 0d == (0,50)
axHistx.set_ylim( (0, 50) )     # top histo
axHisty.set_xlim( (0, 50) )     # side histo
#axHistx.set_ylim( (0, 30) )     # top histo
#axHisty.set_xlim( (0, 30) )     # side histo

# axes limits
#xlim0, xlimf = -10, 60
axScatter.set_xlim((xlim0, xlimf))    # x-range scatter
axScatter.set_ylim((ylim0, ylimf))  # y-range scatter
axHistx.set_xlim((xlim0, xlimf))      # x-range top hist
axHisty.set_ylim((ylim0, ylimf))    # y-range side hist

plt.savefig("output.eps", orientation='portrait', format='eps',
        bbox_inches='tight')

