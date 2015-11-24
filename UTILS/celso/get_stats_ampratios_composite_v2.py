#!/home/alvizuri/unixutils/python/local/bin/python2.7
#!/opt/antelope/python2.7.2/bin/python2.7
#----------------------------------------------------------
# this script plots histograms for each PLUTONS station.
# data comes from compilation table with
#
# input: table with station amp (max) and stdevs
# fixedname         lon                 lat      varname           lo             la        ln(Vobs/Vsyn) Vobs     Vsyn       ln(Robs/Rsyn) Robs       Rsyn
# 0                  1                  2           3              4               5           6          7           8        9           10          11
# PL03              -66.9451        -22.0156 20100601073604690 -6.714330e+01 -2.240580e+01 1.580000 1.810000e-06 3.720000e-07 -0.460000 6.970000e-07 1.100000e-06
# PL03              -66.9451        -22.0156 20100601221352679 -6.715120e+01 -2.241310e+01 1.630000 1.390000e-06 2.740000e-07 -0.390000 5.280000e-07 7.810000e-07
#
#            OR
#
# 20120202125816847 -6.702280e+01 -2.238860e+01 PLWB_XP        -66.8347 -22.3537 4.390000 7.890000e-06 9.820000e-08 1.370000 6.990000e-07 1.770000e-07
# 20120202125816847 -6.702280e+01 -2.238860e+01 PLRR_XP        -66.8820 -22.2611 4.310000 6.450000e-06 8.660000e-08 3.670000 6.150000e-06 1.560000e-07
#
#
# output: hist_station_amplification_utu60.pdf
#
# 20140728 celsoa - 
#----------------------------------------------------------

#-----------------------------------------------------------
# switches

#compare_components=1    # 0 => V2V R2R
                        # 1 => V2R

#-----------------------------------------------------------

# compute histogram for each station
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# check data from input
if len(sys.argv)!=4:
    sys.exit('usage: %s data_component_V data_component_R [a2a] [v2r] \n' % sys.argv[0])

# switch for histogram type
# a2a => V2V + R2R
# v2r => V2R_obs + V2R_syn + scatter
compare_components=sys.argv[3]
if compare_components!="a2a" and compare_components!="v2r":
    sys.exit('%s: abort. histogram type \"%s\" not recognized.' %
            (sys.argv[0],compare_components))

# function to read input files
def read_data(inputfile):
 
    # check that file exists
    sys.stderr.write("%s: working on inputfile: %s\n" % (sys.argv[0],  inputfile))
    # 20100516063454464 -6.718560e+01 -2.226000e+01 6.399000e-01 PL03_XP -66.9451 -22.0156 4.6285 8.540000e-06 3.130000e-06
    # 20100516085611725 -6.718490e+01 -2.226290e+01 5.061200e+00 PL03_XP -66.9451 -22.0156 4.6285 3.220000e-06 1.050000e-06
    # 20100601073604690 -6.714330e+01 -2.240580e+01 1.089700e+00 PL03_XP -66.9451 -22.0156 4.6285 4.960000e-07 1.910000e-07
    data = np.genfromtxt(inputfile, dtype=None, names = ['evid', 'evlo', 'evla',
        'evdp','stname', 'stlo', 'stla', 'stdp', 'amp1', 'amp2'])

    ramp1 = data['amp1']
    ramp2 = data['amp2']

    # check that both files have same number of points
    npts1 = len(ramp1)
    npts2 = len(ramp2)
    if (npts1 == npts2):
        sys.stderr.write("%s: npts (vertical, radial)= (%d, %d)\n\n" %
                (sys.argv[0], npts1, npts2))
    else:
        sys.error("%s: abort. %d != %d !\n" % (npts1, npts2))

    return(ramp1, ramp2, npts1)

[Vobs, Vsyn, npts] = read_data(sys.argv[1])
[Robs, Rsyn, npts] = read_data(sys.argv[2])

# file_data_vertical = str(sys.argv[1])
# file_data_radial = str(sys.argv[2])
# 
# sys.stderr.write("%s: working on inputfile: %s\n\n" % (sys.argv[0],  file_data_vertical))
# sys.stderr.write("%s: working on inputfile: %s\n\n" % (sys.argv[0],  file_data_radial))
# 
# data_ver = np.genfromtxt(file_data_vertical, dtype=None, names = ['evid', 'evlo', 'evla', 'evdp','stname', 'stlo', 'stla', 'stdp', 'amp1', 'amp2'])
# data_rad = np.genfromtxt(file_data_radial,   dtype=None, names = ['evid', 'evlo', 'evla', 'evdp','stname', 'stlo', 'stla', 'stdp', 'amp1', 'amp2'])
# 
# # extract data into arrays
# Vobs = data['vobs']
# Robs = data['robs']
# 
# Vsyn = data['vsyn']
# Rsyn = data['rsyn']

#sys.exit()

#-----------------------------------------------------------
# compute ratios
# note: simply swap commands below for V/R or V/V
#-----------------------------------------------------------
if compare_components == "a2a":
    # get amplitudes between observed and synthetic
    sys.stderr.write("%s: computing V/V, R/R\n" % sys.argv[0])
    ratio_Vo_Vs = np.log(abs(Vobs)/abs(Vsyn))
    ratio_Ro_Rs = np.log(abs(Robs)/abs(Rsyn))
else:
    # get amplitudes between Vertical and Radial
    sys.stderr.write("%s: computing V/R (obs), V/R (syn)\n" % sys.argv[0])
    ratio_Vo_Vs = np.log(abs(Vobs)/abs(Robs))
    ratio_Ro_Rs = np.log(abs(Vsyn)/abs(Rsyn))

# histogram parameters
plt.rcParams.update({'font.size': 8})
mygray    = [0.7, 0.7, 0.7]
col_bar1  = [0.5, 0.7, 1.0]
col_bar1  = [0.5, 1.0, 0.7]
col_vline = [1.0, 0.0, 0.0]
#alpha_var = 0.98

# define number of bins and bin edges (important)
# plot bars BOUNDED by limits defined in arange
xmax = 6   # 10 or 6
xmin = -xmax
dbin = 0.5
binsize = np.arange(xmin, xmax, dbin)

# plot bars CENTERED at limits defined in arange
#bin0 = xmin + (dbin/2)
#binf = xmax + dbin
#dbin = 0.5
#binsize = np.arange(bin0, binf, dbin)
#binsize = 20

sys.stderr.write("min V= %f max V=%f \n" % (min(ratio_Vo_Vs), max(ratio_Vo_Vs)))
sys.stderr.write("min R= %f max R=%f \n\n" % (min(ratio_Ro_Rs), max(ratio_Ro_Rs)))

# check to ensure that histogram calculation includes full range of data!!
if (min(ratio_Vo_Vs) < xmin) or  (max(ratio_Vo_Vs) > xmax) \
        or (min(ratio_Ro_Rs) < xmin) or (max(ratio_Ro_Rs) > xmax):
        sys.stderr.write("\n*** WARNING ** data ranges larger than histogram ranges!:\n")
        sys.stderr.write("histogram ranges:\nmin %f max %f\nData ranges:\n" % (xmin, xmax))
        sys.stderr.write("min V %f max V %f \n" % (min(ratio_Vo_Vs), max(ratio_Vo_Vs)))
        sys.stderr.write("min R %f max R %f \n\n" % (min(ratio_Ro_Rs), max(ratio_Ro_Rs)))

# max amp for histogram plots
if compare_components == "a2a":
    ymax_hist = 250     # V2V, R2R
else:
    ymax_hist = 400     # V/R

xmax_hist = 6
xmax_hist_scatter = xmax_hist

# output, title...
main_title ="npts = %d" % (npts)

# begin canvas (?)
# key: make adjust figsize so all plots fit on a 
# squares of equal sizes (regardless individual amplitudes)
#fig = plt.figure(figsize=(8.0, 2.25))
#plt.subplots_adjust(wspace = 0.2, bottom = 0.15)
fig = plt.figure(figsize=(5.5, 2.6))

# calc histograms
hist_Vo_Vs, hist_edges_Vo_Vs = np.histogram(ratio_Vo_Vs, bins=binsize)
hist_Ro_Rs, hist_edges_Ro_Rs = np.histogram(ratio_Ro_Rs, bins=binsize)

#---------------------------------------------------------- 
# plot histogram obs PR/PV
#---------------------------------------------------------- 
#ax1 = plt.subplot(1,3,1, adjustable='box', aspect=0.015)
ax1 = plt.subplot(1,2,1, adjustable='box')
#ax1 = plt.subplot(1,2,1)
plt.axvline(color=[0,0,0], linewidth=0.1)
plt.bar(hist_edges_Vo_Vs[:-1], hist_Vo_Vs, width=dbin, color=col_bar1, linewidth=0.5)

if compare_components == "a2a":
    plt.xlabel('ln(Vobs/Vsyn)')     # Vo / Vs
else:
    plt.xlabel('Observed ln(V/R)')  # V / R

plt.ylabel('Count')
plt.ylim([0, ymax_hist])
plt.xlim([-xmax_hist, xmax_hist])

ax1.xaxis.set_minor_locator(MultipleLocator(1))
ax1.xaxis.set_major_locator(MultipleLocator(5))

#---------------------------------------------------------- 
# plot histogram syn PR/PV
#---------------------------------------------------------- 
#ax2 = plt.subplot(1,3,2, adjustable='box', aspect=0.02)
ax2 = plt.subplot(1,2,2)
plt.axvline(color=[0,0,0], linewidth=0.1)
plt.bar(hist_edges_Ro_Rs[:-1], hist_Ro_Rs, width=dbin, color=col_bar1, linewidth=0.5)

if compare_components == "a2a":
    plt.xlabel('ln(Robs/Rsyn)')     # Vo / Vs
else:
    plt.xlabel('Synthetic ln(V/R)') # V / R

plt.ylim([0, ymax_hist])
plt.xlim([-xmax_hist, xmax_hist])

ax2.xaxis.set_minor_locator(MultipleLocator(1))
ax2.xaxis.set_major_locator(MultipleLocator(5))

plt.tight_layout(pad=0.4, w_pad=1.5)    # this is needed!!
plt.subplots_adjust(top=0.89)
plt.suptitle(main_title)

# NOTE ps2raster trims too much if saving in ps
outfile="amp_ratios_plot_histogram.eps"
plt.savefig(outfile, orientation='portrait', format='eps')
sys.stderr.write("outfile: %s\n" % outfile)

#---------------------------------------------------------- 
# plot scatter x=syn, y=obs
#---------------------------------------------------------- 
fig = plt.figure(figsize=(2.3, 2.2))
plt.subplots_adjust(wspace = 0.3, bottom = 0.2)
ax3 = plt.subplot(1,1,1, aspect='equal')
plt.axvline(color=[0,0,0], linewidth=0.1)
plt.axhline(color=[0,0,0], linewidth=0.1)
plt.plot([-xmax_hist_scatter, xmax_hist_scatter],[-xmax_hist_scatter, xmax_hist_scatter], '-', color=[1, 0, 0], linewidth=1)
plt.plot(ratio_Ro_Rs, ratio_Vo_Vs, 'o', markersize=3.0, color=col_bar1)
plt.xlim([-xmax_hist, xmax_hist])
plt.ylim([-xmax_hist, xmax_hist])

ax3.xaxis.set_minor_locator(MultipleLocator(1))
ax3.xaxis.set_major_locator(MultipleLocator(5))
ax3.yaxis.set_minor_locator(MultipleLocator(1))
ax3.yaxis.set_major_locator(MultipleLocator(5))

# labels
if compare_components == "a2a":
    plt.xlabel('ln(Vobs/Vsyn)')
    plt.ylabel('ln(Robs/Rsyn)')
else:
    plt.xlabel('Synthetic ln(V/R)') # V / R
    plt.ylabel('Observed ln(V/R)')  # V / R

# count and print points on each quadrant
npts_ru = 0
npts_rd = 0
npts_lu = 0
npts_ld = 0

for i in np.arange(0, len(ratio_Vo_Vs)):
    if   ((ratio_Ro_Rs[i]>0) & (ratio_Vo_Vs[i]>0)):
        npts_ru +=1
    elif ((ratio_Ro_Rs[i]>0) & (ratio_Vo_Vs[i]<0)):
        npts_rd +=1
    elif ((ratio_Ro_Rs[i]<0) & (ratio_Vo_Vs[i]>0)):
        npts_lu +=1
    else:
        npts_ld +=1

text_quad_ru = "%03d" % npts_ru
text_quad_rd = "%03d" % npts_rd
text_quad_lu = "%03d" % npts_lu
text_quad_ld = "%03d" % npts_ld
props = dict(facecolor='white', linewidth=0 )
#plt.text(0.80, 0.8, text_quad_ru, transform=ax3.transAxes, bbox=props)  #box around label
plt.text(0.80, 0.8, text_quad_ru, transform=ax3.transAxes)  #box around label
plt.text(0.80, 0.1, text_quad_rd, transform=ax3.transAxes)
plt.text(0.10, 0.8, text_quad_lu, transform=ax3.transAxes)
plt.text(0.10, 0.1, text_quad_ld, transform=ax3.transAxes)

plt.suptitle(main_title)
#debug
#plt.show()

# save pdf
outfile="amp_ratios_plot_scatter.eps"
plt.savefig(outfile, orientation='portrait', format='eps')
sys.stderr.write("outfile: %s\n" % outfile)

