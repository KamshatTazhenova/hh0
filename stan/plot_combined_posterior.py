import numpy as np
import numpy.random as npr
import scipy.stats as ss
import scipy.optimize as so
import os
notk = False
if 'DISPLAY' not in os.environ.keys():
		notk = True
elif os.environ['DISPLAY'] == ':0':
		notk = True
if notk:
	import matplotlib
	matplotlib.use('Agg')
import matplotlib.pyplot as mp
import matplotlib.cm as mpcm
import matplotlib.colors as mpc
import pickle
import pystan as ps
import fnmatch
import getdist as gd
import getdist.plots as gdp
import h5py

def pretty_hist(data, bins, axis, color, density=False, fill=True, \
				ls='-', zorder=None, label=None):

	hist, bin_edges = np.histogram(data, bins=bins, density=density)
	bins_to_plot = np.append(bins, bins[-1])
	hist_to_plot = np.append(np.insert(hist, 0, 0.0), 0.0)
	if zorder is not None:
		if label is not None:
			axis.step(bins_to_plot, hist_to_plot, where='pre', \
					  color=color, linestyle=ls, zorder=zorder, \
					  label=label)
		else:
			axis.step(bins_to_plot, hist_to_plot, where='pre', \
					  color=color, linestyle=ls, zorder=zorder)
		if fill:
			axis.fill_between(bins_to_plot, hist_to_plot, \
							  color=color, alpha=0.7, step='pre', \
							  zorder=zorder)
	else:
		if label is not None:
			axis.step(bins_to_plot, hist_to_plot, where='pre', \
					  color=color, linestyle=ls, label=label)
		else:
			axis.step(bins_to_plot, hist_to_plot, where='pre', \
					  color=color, linestyle=ls)
		if fill:
			axis.fill_between(bins_to_plot, hist_to_plot, \
							  color=color, alpha=0.7, step='pre')

def norm_dist(dist, delta_x):
	norm = np.sum(dist)
	return dist / norm
	

# plotting settings
lw = 1.5
mp.rc('font', family = 'serif')
mp.rcParams['text.latex.preamble'] = [r'\boldmath']
mp.rcParams['axes.linewidth'] = lw
mp.rcParams['lines.linewidth'] = lw
cm = mpcm.get_cmap('plasma')
cols = [cm(x) for x in np.linspace(0.1, 0.9, 10)]

n_event = 1
n_rpt = 1
h_0_true = 70.0
q_0_true = -0.5

fixed_n_bns = True
z = False
vary_m_c = False
inc_rate_redshift = False
cut_bad_runs = True
n_overlay = 1
base = 'bias_test'
pars = ['h_0']
par_names = ['H_0']
par_ranges = {}
if z:
	
	pars += ['z']
	par_names += ['z']
else:
	q_0_true = None
n_pars = len(pars)
if vary_m_c:
	base += '_vary_m_c'
if inc_rate_redshift:
	base += '_rr'

n = 20   #number of individual posteriors
n_samples = 5000
n_pars = 1
n_rpt = 1
n_chains = 4
samples = np.zeros((n_chains * n_samples, n_pars, n_rpt, n))

for k in range(n):
	filename = 'runs/' + str(k+1) + '/bias_test_samples_1_rpts.h5'
	f = h5py.File(filename, 'r')
	raw_samples = f['samples'][:]
	for i in range(0, n_chains):
		for j in range(0, n_pars):
			samples[i * n_samples: (i + 1) * n_samples, j, :, k] = \
				raw_samples[:, i, j, :]
mp.figure(1)
n_grid = 1500
h_0_min_plot = 10
h_0_max_plot = 150
h_0_grid = np.linspace(h_0_min_plot, h_0_max_plot, n_grid)
d_h_0 = h_0_grid[1] - h_0_grid[0]
h_0_post_product = 1.0

print 'mean = ', np.mean(samples)

if (n > 100):
	
	m = int(n/100)
	print 'm=', m
	for p in range(m):
		for i in range(100*p, 100*(p+1)):	
			gd_samps = gd.MCSamples(samples=samples[:,:,0, i], \
								names=pars, labels=par_names, \
								ranges=par_ranges)
			h_0_post = norm_dist(gd_samps.get1DDensity('h_0').Prob(h_0_grid), d_h_0)
			h_0_post_product = h_0_post_product*h_0_post
			mp.plot(h_0_grid, h_0_post, alpha=0.5)
		h_0_post_product = h_0_post_product/sum(h_0_post_product)
	for j in range(100*m, n):
		gd_samps = gd.MCSamples(samples=samples[:,:,0, j], \
							names=pars, labels=par_names, \
							ranges=par_ranges)
		h_0_post = norm_dist(gd_samps.get1DDensity('h_0').Prob(h_0_grid), d_h_0)
		h_0_post_product = h_0_post_product*h_0_post
		mp.plot(h_0_grid, h_0_post, alpha=0.5)
		
	h_0_post_final = h_0_post_product/sum(h_0_post_product)
	mp.plot(h_0_grid, h_0_post_final, alpha=1.0, color='black')
	mp.ylim(0.0)
	print sum(h_0_post_final)
	
else:	
	for i in range(n):
		gd_samps = gd.MCSamples(samples=samples[:,:,0, i], \
								names=pars, labels=par_names, \
								ranges=par_ranges)
		h_0_post = norm_dist(gd_samps.get1DDensity('h_0').Prob(h_0_grid), d_h_0)
		h_0_post_product = h_0_post_product*h_0_post
		#print 'product', h_0_post_product 
		mp.plot(h_0_grid, h_0_post, alpha=0.5)
		
	h_0_post_final = h_0_post_product/sum(h_0_post_product)
	mp.plot(h_0_grid, h_0_post_final, alpha=1.0, color='blue')
	mp.ylim(0.0)
	print sum(h_0_post_final)

	
max_h_0_index = np.argmax(h_0_post_final)
print 'max prob = ', h_0_post_final[max_h_0_index]
max_h_0 = h_0_grid[max_h_0_index]
print 'max_post = ', max_h_0

hdr = 0.0 #highest density region
i = len(h_0_post_final) - 1
hist1_sorted = sorted(h_0_post_final)
print 'sorted = ', hist1_sorted[i]
while(hdr < 0.68):
	hdr = hdr + hist1_sorted[i]
	i=i-1
print 'hdr = ', hdr	
print i+1
print sum(hist1_sorted[i+1:len(h_0_post_final)])

print 'prob line on = ', hist1_sorted[i+1]

#searching for index
almost_prob_line = abs(h_0_post_final - hist1_sorted[i+1])
index = np.argmin(almost_prob_line)
#print hist1[index]

k_left = 0
while(h_0_post_final[k_left] < hist1_sorted[i+1]):
	k_left = k_left+1

print 'left_h_0 = ', h_0_grid[k_left] - max_h_0

k_right = len(h_0_post_final)-1
while(h_0_post_final[k_right] < hist1_sorted[i+1]):
	k_right = k_right - 1

print 'right_h_0 = ', h_0_grid[k_right] - max_h_0

mp.fill_between(h_0_grid[k_left:k_right], h_0_post_final[k_left:k_right], color='blue', alpha = 0.3)
mp.xlabel(r'$H_0\,{\rm (km/s/Mpc)}$', fontsize=16)
mp.ylabel(r'${\rm P}(H_0|{\rm sample})$', fontsize=16)
mp.legend(['$H_0 =$' + str(round(max_h_0, 2)) + '(+' +  str(round(h_0_grid[k_right] - max_h_0, 2)) +')' + '(' +  str(round(h_0_grid[k_left] - max_h_0, 2)) + ')'])
mp.show()

