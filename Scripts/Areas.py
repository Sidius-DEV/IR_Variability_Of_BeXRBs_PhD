from scipy.integrate import quad
import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#This script calculates the blue and red surface areas and these are compared with the measured equivalent widths
#Equations A.17 and A.18

file = np.loadtxt('HR_ave_velocities.dat')
file2 = np.loadtxt('HR_ave_EW.dat')


MJD = file[:,0]
#cycles = file[:,1]
v_blue = file[:,2]
v_red = file[:,3]

EW_blue = file2[:,2]
EW_red = file2[:,3]

MJD_eph = 43366.275


def integrand(x,eccen):
    return 1.0/((1.0 + eccen*np.cos(x))**2)

ep = 0.4

#print I[0]

a_in = 1.0/(1.0 - ep)
inc = 0.52
v_crit = 525.0
degtorad = np.pi/180.0
radtodeg = 180.0/np.pi

g = open("HR_Areas_ave_out.txt","w")
for i in range(len(MJD)):
    t1 = ((2.0*v_crit)/(v_red[i] - v_blue[i]))**2
    t2 = ((np.sin(inc))**2)/(1.0 - ep**2)
    a_p = t1*t2
    ts1 = 0.5*(a_p**2 - a_in**2)*(1.0 - ep**2)*np.cos(inc)
    MJD_conv = MJD[i] - 2400000.5
    ratio = ((v_red[i]+v_blue[i])/(v_red[i]-v_blue[i]))
    om = np.arccos(ratio*(1.0/ep))
    cos_om = ratio*(1.0/ep)
    f01 = np.arccos(-ep*cos_om) - om
    f01_deg = f01*radtodeg
    
    f02 = (2.0*np.pi - np.arccos(-ep*cos_om)) - om
    f02_deg = f02*radtodeg
    
    
    I_b = quad(integrand,f01,f02,args=(ep))
    ts2 = I_b[0]
    S_blue = ts1*ts2
    I_r = quad(integrand,f02,f01+2.0*np.pi,args=(ep))
    ts3 = I_r[0]
    S_red = ts1*ts3
    ratio_areas = S_blue/S_red
    ratio_EW = EW_blue[i]/EW_red[i]
    g.write("%0.3f  %0.3f   %0.3f   %0.3f   %0.3f   %0.3f   %0.3f   %0.3f   %0.3f   %0.3f   %0.3f\n" %(MJD[i],S_blue,S_red,EW_blue[i],EW_red[i],ratio_areas,ratio_EW,I_b[0],I_r[0],f01_deg,f02_deg))
g.close()


