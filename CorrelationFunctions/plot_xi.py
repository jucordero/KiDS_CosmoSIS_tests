import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# This routine plots the two point shear correlations funcionts xi_plus and xi_minus as a funcion of angle, for each thomographic combination of bins.
# The plot follows the same style used on the KiDS-450 cosmological parameter paper. Two png files are created, one for each correlation function. 

xi1 = np.genfromtxt('xi_1-1')
xi2 = np.genfromtxt('xi_1-2')
xi3 = np.genfromtxt('xi_1-3')
xi4 = np.genfromtxt('xi_1-4')
xi5 = np.genfromtxt('xi_2-2')
xi6 = np.genfromtxt('xi_2-3')
xi7 = np.genfromtxt('xi_2-4')
xi8 = np.genfromtxt('xi_3-3')
xi9 = np.genfromtxt('xi_3-4')
xi10 = np.genfromtxt('xi_4-4')

xi = [xi1,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9,xi10]

f, ((axoff6,axoff5,axoff4,ax44),(axoff3,axoff2,ax33,ax34),(axoff1,ax22,ax23,ax24),(ax11,ax12,ax13,ax14)) = plt.subplots(4,4,sharex='col',sharey='row')
ax = (ax11,ax12,ax13,ax14,ax22,ax23,ax24,ax33,ax34,ax44,axoff1,axoff2,axoff3,axoff4,axoff5,axoff6)

for i in range(10):
	x = xi[i][3:9,0]
	y = xi[i][3:9,2]
	ax[i].plot(x,x*y*10000,'.',color='k', marker='o')
	ax[i].yaxis.tick_right()
	ax[i].set_ylim(-1,4)
	ax[i].set_xlim(4,300)
	ax[i].set_yticks([0,1,2,3])
	ax[i].set_xticks([10,100])
	ax[i].set_xscale("log", nonposx='clip')

for i in range(10,16):
	ax[i].axis('off')

f.text(0.5, 0.04, r'$\theta$', ha='center')
f.text(0.96, 0.5, r'$\xi\theta_{-}$[$10^{-4}$ arcmin]', va='center', rotation='vertical')

f.subplots_adjust(hspace=0, wspace=0)

plt.savefig('correlation_xim_ATHENA.pdf', dpi=200)

f, ((ax14,ax24,ax34,ax44),(ax13,ax23,ax33,axoff6),(ax12,ax22,axoff4,axoff5),(ax11,axoff1,axoff2,axoff3)) = plt.subplots(4,4,sharex='col',sharey='row')
ax = (ax11,ax12,ax13,ax14,ax22,ax23,ax24,ax33,ax34,ax44,axoff1,axoff2,axoff3,axoff4,axoff5,axoff6)

for i in range(10):
	x = xi[i][:7,0]
	y = xi[i][:7,1]
	ax[i].plot(x,x*y*10000,'.',color='k', marker='o')
	if i == 0 or i == 4 or i == 7 or i == 9:
		ax[i].xaxis.tick_bottom()
		ax[i].set_xlabel(r'$\theta$')
	ax[i].set_ylim(-1,4)
	ax[i].set_xlim(0.5,80)
	ax[i].set_yticks([0,1,2,3])
	ax[i].set_xticks([1,10])
	ax[i].set_xscale("log", nonposx='clip')

for i in range(10,16):
	ax[i].axis('off')

f.text(0.04, 0.5, r'$\xi\theta_{+}$[$10^{-4}$ arcmin]', va='center', rotation='vertical')

f.subplots_adjust(hspace=0, wspace=0)

plt.savefig('correlation_xip_ATHENA.pdf', dpi=200)
