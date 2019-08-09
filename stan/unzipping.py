import bz2file
import pandas as pd
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
import pylab as plb
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.stats import norm
import h5py

#catalog_filename='4724.csv.bz2'
catalog_filename='5563.csv.bz2'
with bz2file.BZ2File(catalog_filename) as galaxy_fd:
    galaxy_sample = pd.read_csv(galaxy_fd, sep=',', comment='#', na_values = '\N')


print np.min(galaxy_sample, axis=0)
arr = galaxy_sample.to_numpy()

#cosmological redshift
max_z_cos = np.max(arr[:, 4])
min_z_cos = np.min(arr[:, 4])

plt.figure(1)
bins_obs = np.linspace(min_z_cos, max_z_cos, 100)
plt.hist(arr[:,4], bins_obs, density='true')
plt.xlabel('z', fontsize=12)
plt.ylabel('Cosmological redshift', fontsize=12)

#redshift with peculiar velocity
max_z_true = np.max(arr[:, 5])
min_z_true = np.min(arr[:, 5])

plt.figure(2)
bins_obs = np.linspace(min_z_true, max_z_true, 100)
plt.hist(arr[:,5], bins_obs, density='true')
plt.xlabel('z', fontsize=12)
plt.ylabel('True redshift', fontsize=12)

#observed photmetric redshift
max_z_obs = np.max(arr[:, 6])
min_z_obs = np.min(arr[:, 6])

plt.figure(3)
bins_obs = np.linspace(min_z_obs, max_z_obs, 100)
plt.hist(arr[:,6], bins_obs, density='True')
plt.xlabel('z', fontsize=12)
plt.ylabel('Photometric redshift', fontsize=12)

#observed magnitude
plt.figure(4)
plt.scatter(arr[:,4], arr[:,1], alpha=0.4, s=0.2)
#observed mag cut
mag = arr[(arr[:,1] <= 23.9), :]
plt.scatter(mag[:,4], mag[:,1], alpha=0.4, s=0.2)
plt.xlabel('Cosmo z', fontsize=12)
plt.ylabel('Observed mag', fontsize=12)


#absolute magnitude
plt.figure(5)
plt.scatter(arr[:,4], arr[:,0], alpha=0.4, s=0.2)
plt.gca().invert_yaxis()
plt.xlabel('Cosmo z', fontsize=12)
plt.ylabel('Absolute mag', fontsize=12)

#absolute mag cut
plt.figure(6)
plt.scatter(mag[:,4], mag[:,0], alpha=0.2, s=0.2, color = 'grey')
abs_mag = mag[(mag[:,0] <= -15.5), :]
#plt.scatter(abs_mag[:,4], abs_mag[:,0], alpha=0.4, s=0.2)
plt.gca().invert_yaxis()
plt.xlabel('z', fontsize=12)
plt.xlim(0.06, 0.3)
plt.ylabel('Absolute mag, $i$', fontsize=12)


#observed redshift after cut
max_z_obs = np.max(abs_mag[:, 6])
min_z_obs = np.min(abs_mag[:, 6])

plt.figure(3)
bins_obs = np.linspace(min_z_obs, max_z_obs, 100)
plt.hist(abs_mag[:,6], bins_obs, alpha = 0.6, density='True')
plt.xlabel('z', fontsize=12)
plt.ylabel('Density', fontsize=12)

#fitting a stright line
hist, bin_edges = np.histogram(abs_mag[:,6], bins_obs, density = True)
x = bins_obs[0:len(bins_obs)-1] + (bins_obs[-1] - bins_obs[0])/200
fit = np.polyfit(x, hist, 1)
y = fit[1] + fit[0]*x
print 'fit[0] = ', fit[0]
print 'fit[1] = ', fit[1]
#plt.plot(x, y, color='g', label = 'straight line fit')



#new catalog of observed redshift
c = 2.998e5
arr = arr[(arr[:,4] <= 0.3),:]
mag = arr[(arr[:,1] <= 23.9), :]
abs_mag = mag[(mag[:,0] <= -15.5), :]

sig_v_pec = 500.0
sig_z_obs = 0.0
z_with_v_pec = np.zeros(len(abs_mag))
z_ph_obs = np.zeros(len(abs_mag))
for i in range(len(abs_mag)):
	v_pec_true = npr.randn(1) * (sig_v_pec + sig_z_obs*c)
	z_ph_obs[i] = abs_mag[i,4] + v_pec_true / c
	#z_ph_obs[i] = z_with_v_pec[i] + npr.randn(1) * sig_z_obs

max_z_obs = np.max(z_ph_obs)
min_z_obs = np.min(z_ph_obs)
bins_obs = np.linspace(min_z_obs, max_z_obs, 100)

base = 'bias_test_'
with h5py.File(base + 'z_obs1.h5', 'w') as f:
		f.create_dataset('z_ph_obs', data=z_ph_obs)

plt.figure(7)
plt.hist(z_ph_obs, bins_obs, alpha = 0.6, density='True')
plt.xlabel('z', fontsize=12)
plt.ylabel('Density', fontsize=12)

print 'z_min_obs = ', min_z_obs

filename = 'bias_test_z_obs1.h5'
f = h5py.File(filename, 'r')
raw_samples = f['z_ph_obs']
obs_z_from_file = np.zeros(len(abs_mag))
for i in range(0, len(abs_mag)):
	obs_z_from_file[i] = raw_samples[i]

plt.figure(8)
max_z_obs = np.max(obs_z_from_file)
min_z_obs = np.min(obs_z_from_file)
bins_obs = np.linspace(min_z_obs, max_z_obs, 100)
plt.hist(obs_z_from_file, bins_obs, alpha = 0.6, density='True')
plt.xlabel('z', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.show()


