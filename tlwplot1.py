import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.gridspec as gridspec
import numpy as np
from cmath import sqrt
import cmath as cm

#%-------------------------------------------------------------
#% Matlab Subroutine tlwplot.M
#% Called from tlwmenu.M--interactive menu system.
#%
#% By Robert Hart for Meteo 574 / Fall 1995
#% Penn State University Meteorology
#% Updated 1 March 2018 to increase npts given faster cpus
#%
#% Program which analyzes and simulates the flow over
#% an isolated mountain in the presence of a two-
#% layer atmosphere.
#%
#% Parameters received from menu:
#%-------------------------------
#% Lupper - Scorer parameter of upper layer
#% Llower - Scorer parameter of lower layer
#% U      - Surface wind speed (in m/s)
#% H      - Height above ground of 2-layer interface (in m)
#% a      - Half-width of mountain (in m)
#% ho     - maximum height of mountain
#% xdom   - Horizontal extent of domain (m)
#% zdom   - Vertical extent of domain (m)
#% mink   - Minimum wavenumber in Fourier analysis of flow
#% maxk   - Maximum wavenumber in Fourier analysis of flow
#%
#%--------------------------------------------------------------

Lupper=0.0004;
Llower=0.1;
U=20;
H=3.5;
a=0.5;
h0=3;
xdom=40;
zdom=10;
mink=0;
maxk=30;



def tlwplot(Lupper,Llower,U,H,a,h0,xdom,zdom,mink,maxk):
	npts = 100;  # cells in each direction
	dk = 0.367 / a; # % wavenumber interval size (smaller the better, but  need to watch calc.time!)
	nk = (maxk - mink) / dk;   # loops for wave# integration

	minx = -.25 * xdom; # % leftmost limit of domain
	maxx = .75 * xdom; #% rightmost limit of domain
	minz = 0; #% lower limit of domain
	maxz = zdom; #% upper limit of domain


	matrix1 = np.zeros([npts + 1, npts + 1]); #% temp matrix used in integration
	matrix2 = np.zeros([npts + 1, npts + 1]); #% temp matrix used in integration
	matrix3 = np.zeros([npts , npts ]); #% sum matrix used in integration
	umesh=U*np.ones(matrix3.shape)
	dx = (maxx - minx) / npts; #% grid cell size in horizontal
	dz = (maxz - minz) / npts; #% grid cell size in vertical

	x = np.arange(minx, maxx,dx); #% array of x gridpoints
	z = np.arange(minz, maxz,dz); #% array of z gridpoints
	k = np.arange(mink, maxk,dk); #% array of wavenumbers

	ht=0.
	x,z=np.meshgrid(x,z)
	nk=np.int(np.round(nk,decimals=0))

	for kloop in range(1,nk):

		kk = k[kloop]; #% horiz         wavenumber
		print(kk)
		m = sqrt((Llower * Llower - kk * kk)); #% vert         wave  # in lower layer
		
		n = sqrt((kk * kk - Lupper * Lupper)); #% vert     wave  # in upper layer
		print(m,n)

		if ((m+1j*n) == 0): # % check for divide by zero
			r = 9e99;
		else:
			r = (m-1j*n) / (m+1j*n); # % reflection coefficient end;

		R = r * cm.exp(1j*2* m * H);# % reflection     calculation

		A = ((1 + r) * (cm.exp(H*n+1j*H*m))) / (1 + R);

		C = 1 / (1 + R);

		D = R*C;
		print(dk,'dk')
		hs =  np.pi*a*h0 *cm.exp(-a * np.abs(kk));
		ht = ht + np.pi * dk * a * cm.exp(-a * np.abs(kk)); 
		print(hs,ht,'hsht')
		aboveH = (A * np.exp(-z * n)) ;# % calculate     w in each     layer
		belowH = (C * np.exp(1j*z* m)+D*np.exp(-1j*z*m))

		aboveidx=np.where((z>H))
		belowidx=np.where((z<=H))

		zmat=np.ones(aboveH.shape)
		zmat[aboveidx]=aboveH[aboveidx]*z[aboveidx]
		zmat[belowidx]=belowH[belowidx]*z[belowidx]
		matrix2 = zmat*(-1j*kk*hs*U)*np.exp(-1j*x*kk)

		if kloop > 1: #% trapezoidal integration / summation
			matrix3 = matrix3 + .5 * (matrix1+matrix2) * dk;
		else:
			matrix1 = matrix2;

	w = np.real(matrix3/ ht);
	plt.figure()
	plt.clf()
	plt.close('all')
   # fig = plt.figure(figsize=(1, 1))
   # gs = gridspec.GridSpec(nrows=1, ncols=1, height_ratios=[1, 1, 2])

    #  Varying density along a streamline
   # ax0 = fig.add_subplot(gs[0, 0])

    #plt.streamplot(Z[0,100:200],X[:,100:200],U,w)
    #plt.contourf(Z[:,0:100],X[:,0:100],w)
	Z = np.arange(minz,maxz,dz)
	X = np.arange(minx,maxx,dx)
	
	plt.contourf(x,z,w)
	plt.streamplot(x,z,umesh,w, density=[0.5, 1],color='k')
	plt.ylim(0,10)
	Hline = H * np.ones([1, npts + 1]); #% create    interface     line
	plt.plot( Hline, 'm--'); #% draw     interface   line in magenta!
	plt.savefig('mw_fig.png')
	return (w)
w=tlwplot(Lupper,Llower,U,H,a,h0,xdom,zdom,mink,maxk)

