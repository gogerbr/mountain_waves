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
Llower=0.001;
U=20.;
H=3.5;
a=0.71;
h0=0.5;
xdom=40;
zdom=10;
mink=0.;
maxk=30.;





def tlwplot(Lupper,Llower,U,H,a,h0,xdom,zdom,mink,maxk):
	npts = 1000;  # cells in each direction
	dk = 0.367 / a; # % wavenumber interval size (smaller the better, but  need to watch calc.time!)
	nk = (maxk - mink) / dk;   # loops for wave# integration

	minx = -0.25 * xdom; # % leftmost limit of domain
	maxx = 0.75 * xdom; #% rightmost limit of domain
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
		m = cm.sqrt((Llower * Llower - kk * kk)); #% vert         wave  # in lower layer
		
		n = cm.sqrt((kk * kk - Lupper * Lupper)); #% vert     wave  # in upper layer

		if ((m+1j*n) == 0): # % check for divide by zero
			r = 9e99;
		else:
			r = (m-1j*n) / (m+1j*n); # % reflection coefficient end;

		R = r*np.exp(1j*2.*m*H);# % reflection     calculation

		A = ((1 + r) * (cm.exp(H*n+1j*H*m))) / (1 + R);
		C = 1 / (1 + R)
		D = R*C;

		hs =  np.pi*a*h0 *np.exp(-a * kk);
		ht = ht + np.pi * dk * a * np.exp(-a * kk); 
		print(hs,'hs')

	
		aboveH = (A * np.exp(-z * n))*(z>H) ;# % calculate     w in each     layer

		belowH = (C * np.exp(1j*z*m)+D*np.exp(-1j*m*z))*(z<=H)


	  	

		matrix2 = (-1j*U*kk*hs*(aboveH+belowH))*np.exp(kk*-1j*x)

		if kloop > 1: #% trapezoidal integration / summation
			matrix3 = matrix3 + .5 * (matrix2) * dk;
			print(np.min(matrix1),np.min(matrix2),np.min(matrix3))
			print(np.max(matrix1),np.max(matrix2),np.max(matrix3))
		else:
			matrix1 = matrix2;

	w = np.real(matrix3/ ht);
	plt.figure()


	Z = np.arange(minz,maxz,dz)
	X = np.arange(minx,maxx,dx)
	
	plt.contourf(x,z,w,cmap='coolwarm')
	plt.colorbar()


	#Hline = H * np.ones([1, npts + 1]); 

	#plt.plot( Hline, 'm--'); #% draw     interface   line in magenta!
	plt.savefig('mw_fig.png')
	plt.close('all')

	return (w)

w=tlwplot(Lupper,Llower,U,H,a,h0,xdom,zdom,mink,maxk)
print(np.max(w),np.min(w))
