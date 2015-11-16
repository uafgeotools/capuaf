#!/home/alvizuri/unixutils/python/local/bin/python2.7
#####!/opt/antelope/python2.7.2/bin/python2.7
#----------------------------------------------------------
# script to compute stats mean/median/stdev and histogram
#
# input: 2 files 
#
#   NAME                    X               Y           Z      MEDIAN    MEAN     STDEV   NPTS      
#   20100516063454464 -6.718560e+01 -2.226000e+01 6.399000e-01 PLLO_XP -67.0793 -22.3337 4.5790 2.830000e-06 1.900000e-06
#           col0            col1          col2         col3    col4       col5    col6    col7    amp1          amp2
#
# output
#   a "stats" file that aggregates data 
#   histogram and output values: evid, npts, mean, ave, stdev
#
# 20150111 cralvizuri - 
#----------------------------------------------------------

import sys
import os   # to verify empty files
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# verify that inputfile exists, then extract data from it
if len(sys.argv)!=3:
    sys.stderr.write('usage: %s datafile1.txt datafile2.txt\n' % sys.argv[0])
    sys.stderr.write('NOTE datafiles have amplitudes for each station and component' % sys.argv[0])
    sys.exit()

input_ratios_V = str(sys.argv[1])
input_ratios_R = str(sys.argv[2])

def get_histograms(inputfile):

    if os.stat(inputfile).st_size==0:
        sys.exit('%s: abort. file is empty: %s\n' % (sys.argv[0], inputfile))

    try:
        data = np.genfromtxt(inputfile, dtype=None, names = ['col0', 'col1',
            'col2', 'col3', 'col4','col5', 'col6', 'col7', 'amp1', 'amp2'])
    except:
        sys.exit('\n%s: abort. error reading file: %s\n' % (sys.argv[0], inputfile))

    sys.stderr.write("\n%s: processing file: %s\n" % (sys.argv[0], inputfile))

    # extract data into arrays
    amp1    = data['amp1']
    amp2    = data['amp2']

    # if single or empty fields, then abort
    try:
        # will throw exception if empty field
        npts = len(data)
    except TypeError:
        npts = 1
        sys.exit("\n%s: abort. %s has single or empty rows\n" % (sys.argv[0], inputfile))

    # if colums are constant then use them as control label
    # NOTE POSSIBLE PROBLEM: when there are multiple files (eg depth test)
    if data[0][0] == data[1][0]:
        sys.stderr.write("%s: found constant column: EVID\n" % sys.argv[0])
        fixedname = data['col0'][0]
        fixedlo   = data['col1'][0]
        fixedla   = data['col2'][0]
        fixeddp   = data['col3'][0]

    elif data[0][4] == data[1][4]:
        sys.stderr.write("%s: found constant column: STA\n" % sys.argv[0])
        colname   = data['col4'][0]
        fixedlo   = data['col5'][0]
        fixedla   = data['col6'][0]
        fixeddp   = data['col7'][0]

        fixedname = colname.replace("_XP", "")

    else:
        # if both columns constant then something may be wrong...
        sys.exit('%s: abort. input data constant in two columns...' % (sys.argv[0]))

    # compute log ratio for non nans
    logratio = np.log(amp1 / amp2)
    logratio = logratio[~np.isnan(logratio)]

    # exit if NPTS (non nan) = 0 
    npts = np.count_nonzero(~np.isnan(logratio))
    if npts < 1:
        sys.exit('\n%s: abort. npts = %d. all values are NAN: %s\n' % (sys.argv[0],
            npts, inputfile))

    # compute mean, std. output result
    logratio_mean = np.nanmean(logratio, dtype=np.float64)
    logratio_med = np.median(logratio)
    logratio_std = np.nanstd(logratio, dtype=np.float64, ddof=1)

    #if np.isnan(logratio_std):
    #    sys.exit('%s: abort. STDEV is nan' % (sys.argv[0]))

    # output format
    sys.stderr.write('output format: name lon lat dep median mean stdev npts\n')

    # output data to screen
    sys.stderr.write('%s %f %f %f %f %f %f %03d \n' % (
        fixedname, fixedlo, fixedla, fixeddp,
        logratio_med, logratio_mean, logratio_std, npts))

    # output data to file
#   fid1.write('%s %f %f %f %f %f %f %03d \n' % (
#       fixedname, fixedlo, fixedla, fixeddp,
#       logratio_med, logratio_mean, logratio_std, npts))
#   fid1.close()

    # histogram parameters
    # define number of bins and bin edges (important)
    xmin = -5
    xmax = 5
    dbin = 0.5
    bin0 = xmin + (dbin/2)
    binf = xmax + dbin
    binsize = np.arange(bin0, binf, dbin)

    # compute histograms
    hist_logratio, bin_edges = np.histogram(logratio, bins=binsize)

# logratio_mean, logratio_med,  logratio_std 

    return (hist_logratio, bin_edges, fixedname, npts, dbin, xmin, xmax, logratio_mean, logratio_med,  logratio_std)


histoV, bin_edgesV, fixedname, npts, dbin, xmin, xmax, logratio_V_mean, logratio_V_med,  logratio_V_std = get_histograms(input_ratios_V)
histoR, bin_edgesR, fixedname, npts, dbin, xmin, xmax, logratio_R_mean, logratio_R_med,  logratio_R_std = get_histograms(input_ratios_R)

ymax_logratio_V = max(histoV)
ymax_logratio_R = max(histoR)

# nticks, starting tick
main_title = "%s npts=%02d" % (fixedname, npts)
#sys.stdout.write("NOTE data appended to table.amp_ratio.stations_utu60_PV.txt / table.amp_ratio.stations_utu60_PR.txt\n")

# figure out ytick spacing (depends on amount of data)
yticks=2
y0=1

if ymax_logratio_V >= ymax_logratio_R:
    ymax = ymax_logratio_V
else:
    ymax = ymax_logratio_R

if (ymax > 10):
    yticks = 10
    y0 = 10
elif (ymax >= 5) & (ymax <= 10):
    yticks = 2
    y0 = 2
else:
    yticks = 1
    y0 = 1

#---------------------------------------------------------- 
# begin plotting
#---------------------------------------------------------- 
# set output files
outplot="hist_log_ratios_%s.eps" % fixedname

# begin plotting
fig = plt.figure(figsize=(3,2))
plt.suptitle(main_title, fontsize=12)

# histogram colors
plt.rcParams.update({'font.size': 8})
mygray    = [0.7, 0.7, 0.7]
col_bar_1  = [0.5, 0.7, 1.0]
col_bar_2  = [0.5, 1.0, 0.7]
col_vline = [1.0, 0.0, 0.0]

alpha_var = 0.98

#---------------------------------------------------------- 
# histogram log ratio 1
#---------------------------------------------------------- 
ax1 = plt.subplot(2,1,1)
plt.bar(bin_edgesV[:-1], histoV, width=dbin, color=col_bar_1)

# axes range and ticks
plt.axvline(x=0, color=col_vline, linestyle='-', linewidth=1.1)
plt.xlim([xmin, xmax])
plt.ylim([0, ymax])
plt.yticks(np.arange(y0, ymax, yticks))

# annotations
annot_logratio_V = "m=%5.3f\na=%5.3f\ns=%5.3f" % (logratio_V_med, logratio_V_mean, logratio_V_std)
ax1.text(0.05, 0.9, annot_logratio_V, horizontalalignment='left',
        verticalalignment='top', transform = ax1.transAxes)

# remove middle x-axis, reduce spacing bet plots
plt.setp(ax1.get_xticklabels(), visible=False)
plt.subplots_adjust(hspace=0.0, bottom=0.1)

#---------------------------------------------------------- 
# histogram log ratio 2
#---------------------------------------------------------- 
ax2 = plt.subplot(2,1,2, sharex=ax1)
plt.bar(bin_edgesR[:-1], histoR, width=dbin, color=col_bar_2)

# axes range and ticks
plt.axvline(x=0, color=col_vline, linestyle='-', linewidth=1.1)
plt.xlim([xmin, xmax])
plt.ylim([0, ymax])
plt.yticks(np.arange(y0, ymax, yticks))

# annotations
annot_logratio_R = "m=%5.3f\na=%5.3f\ns=%5.3f" % (logratio_R_med, logratio_R_mean, logratio_R_std)
ax2.text(0.05, 0.9, annot_logratio_R, horizontalalignment='left',
        verticalalignment='top', transform = ax2.transAxes)


sys.stderr.write("%s: done. output file: %s\n" % (sys.argv[0] , outplot))
#plt.savefig(outplot, orientation='portrait', papertype=None, format='pdf')
plt.savefig(outplot, orientation='portrait', format='eps')

