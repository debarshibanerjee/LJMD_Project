import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def morse(r,a):
    D = 0.2379
    re = 3.82198
    # return 2 * a * D * ((np.exp(-2*a*(r-re))) - 1*(np.exp(-a*(r-re))))
    return D * ((np.exp(-2*a*(r-re))) - 2*(np.exp(-a*(r-re))))

def lj(r):
    eps = 0.2379
    sigma = 3.405
    # return 48 * eps * (np.power(sigma, 12) / np.power(r, 13)) - 24 * eps * (np.power(sigma, 6) / np.power(r, 7))
    return 4 * eps * ((np.power(sigma, 12) / np.power(r, 12)) - (np.power(sigma, 6) / np.power(r, 6)))

# Argon parameters
D = eps = 0.2379
sigma = 3.405
re = pow(2,1/6)*sigma

r = np.arange(3,10,0.00001)
# lj_y = 48 * eps * (np.power(sigma, 12) / np.power(r, 13)) - 24 * eps * (np.power(sigma, 6) / np.power(r, 7))
lj_y = 4 * eps * ((np.power(sigma, 12) / np.power(r, 12)) - (np.power(sigma, 6) / np.power(r, 6)))

popt, pcov = curve_fit(morse,r,lj_y)

print("a (popt fitted) =", popt)
print("pcov (fitted) =", pcov)

morse_vals = morse(r,popt[0])

plt.title("Estimating appropriate Morse parameters from LJ results")
plt.ylim(-2*eps, +5*eps)
plt.plot(r,morse_vals,'b', label='Morse Potential')
plt.plot(r,lj_y,'g',label='Lennard-Jones Potential')
plt.grid()
plt.legend()
plt.show()
