from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

# Script to perform multiple Gaussian fits on the double peak H-alpha line (on the blue, center and red peak)

file = np.loadtxt("spec_eg6.txt")
x = file[:,0]
y = file[:,1]

#plt.plot(x,y)
#plt.show()


def func(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y



guess = [6557.5, 2.0, 8.0, 6561.0, -0.4, 3.0, 6565.0, 1.1, 2.0]
    #for i in range(1):
#guess += [60+80*i, 46000, 25]


popt, pcov = curve_fit(func, x, y, p0=guess)
print(popt)
fit = func(x, *popt)

plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Normalised flux')
plt.plot(x, y, 'rd', )
plt.plot(x, fit , 'b-', linewidth=2)
#plt.show()
plt.savefig('Gaussian_demo.png')