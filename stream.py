#!/usr/bin/env python
import sys
import matplotlib
import numpy as np
#mountain wave stream function; python version of stream.m
#%-------------------------------------------------------%
###############				%
#% Called by:  tlwplot.m					%
#%							%
#% by Robert Hart for Meteo 574 / Fall 95		%
#% Penn State University Meteorology                     %
#%							%
#% This subroutine performs a streamline analysis	%
#% of the wind field.  					%
#%							%
#% Parameters passed to this routine:			%
#% ----------------------------------			%
#% x	- array of x-gridpoints				%
#% y	- array of y-gridpoints				%
#% U 	- speed of x-direction wind (constant)		%
#% v	- array of y-direction wind			%
#% num	- # of evenly spaced streamlines to draw	%
#%							%
#% NOTE:  Subroutine is written for a constant   	%
#% 	 horizontal wind velocity.			%
#%-------------------------------------------------------%
#hold;				% hold figure%-------------------------------------------------------%
#% Matlab subroutine stream.m				%
#% Called by:  tlwplot.m					%
#%							%
#% by Robert Hart for Meteo 574 / Fall 95		%
#% Penn State University Meteorology                     %
#%							%
#% This subroutine performs a streamline analysis	%
#% of the wind field.  					%
#%							%
#% Parameters passed to this routine:			%
#% ----------------------------------			%
#% x	- array of x-gridpoints				%
#% y	- array of y-gridpoints				%
#% U 	- speed of x-direction wind (constant)		%
#% v	- array of y-direction wind			%
#% num	- # of evenly spaced streamlines to draw	%
#%							%
#% NOTE:  Subroutine is written for a constant   	%
#% 	 horizontal wind velocity.			%
#%-------------------------------------------------------%
#hold;				% hold figure

#Hold(True)


x=10000
y=10000
U=10.
v=2.
num=50


def streamfunction(x,y,U,v,num):
    print("hello!")
  #  global x,y,U,v,num
  #  xsize = len(x);
    xsize = np.arange(1,x);#% size  of x - grid
    ysize = np.arange(1,y);
   # ysize = ysize[2];# % size     of     y - grid
    miny = np.min(ysize);# % min.y - value
    maxy = np.max(ysize);# % max.y - value
    minx = np.min(xsize);# % min.x - value
    maxx = np.max(xsize);# % max.x - value
    print(miny,maxy,minx,maxx)
    dx = (maxx - minx) / xsize; #% x - dir gridspacing
    dy = (maxy - miny) / ysize; #% y - dir grid    spacing
    dh = ysize / num; #% streamline spacing
    print(dx,dy,dh)

    tstep = dx / U; #% time to cross cell

   # mtncolor = [.02 .77 .02]; % color of mountain(lt  green here)

    ycell = 1; #% initial location

    for jj in range(0,num):
        ycell=1+dh*(jj-1)
        if np.any(ycell)<1:
            ycell=1
        elif np.any(ycell)>y:
            ycell=ysize
        print(ycell)

vars=streamfunction(x,y,U,v,num)
