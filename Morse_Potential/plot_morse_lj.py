import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Argon parameters
eps = 0.2379
D = 0.2379
sigma = 3.405
re = pow(2,1/6)*sigma
a = 1.07

r = np.arange(3,10,0.00001)
lj = 4 * eps * ((np.power(sigma, 12) / np.power(r, 12)) - (np.power(sigma, 6) / np.power(r, 6)))
morse = D * ((np.exp(-2*a*(r-re))) - 2*(np.exp(-a*(r-re))))

plt.title("Plotting Morse and LJ Potential")
plt.ylim(-2*eps, +5*eps)
plt.ylabel('Energy')
plt.xlabel('Distance')
plt.text(8,0.7, 'a = 1.07', fontsize=11, bbox=dict(facecolor='red', alpha=0.3))
plt.plot(r,morse,'b', label='Morse Potential')
plt.plot(r,lj,'g', label='Lennard-Jones Potential')
plt.grid()
plt.legend()
plt.show()

