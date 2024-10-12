import numpy as np
import matplotlib.pyplot as plt
import math as Math
import seaborn as sns
from scipy.stats import gaussian_kde
import scipy
import math



_sigmaS = 2.6
sigmaA = 0.0041
_sigmaT = _sigmaS + sigmaA

# Define the function f(r)
def Rdipole(r):
    eta = 1.5
    r2 = r*r
    g = 0

    # Compute isotropic phase function
    sigmaTr = np.sqrt(3 * sigmaA * _sigmaT)
    _alpha = _sigmaS / _sigmaT

    # Effective transport coefficient
    sigmaTr = np.sqrt(3 * sigmaA * _sigmaT)
    
    # Aproximation for the diffuse reflectance (fresnel)
    Fdr = (-1.440 / eta*eta) + (0.710 / eta) + 0.668 + 0.0636 * eta
    A = (1 + Fdr) / (1 - Fdr)    # Boundary condition for the change between refraction indexes

    lu = 1 / _sigmaT
    zr = lu
    zv = -lu * (1.0 + (4.0/3.0) * A) # Negative because it is inside the surface

    dr = np.sqrt(r2 + zr*zr) 
    dv = np.sqrt(r2 + zv*zv) 

    sigmaTrDr = sigmaTr * dr
    sigmaTrDv = sigmaTr * dv

    # Compute main formula
    C1 = zr * (1 + sigmaTrDr) * np.exp(-sigmaTrDr) / (dr*dr*dr)
    C2 = zv * (1 + sigmaTrDv) * np.exp(-sigmaTrDv) / (dv*dv*dv)

    result = _alpha * (C1 - C2) / (4*np.pi)

    return result


def expDecay(epsilon):
    # Define the decay constant
    decay_constant = _sigmaT

    # Compute the exponential decay
    decay = - np.log(1 - epsilon) / decay_constant

    return decay 

def F(r, xi):

    return 1.0 - 0.25 * math.exp(-r) - 0.75 * math.exp(-r/3) - xi

  

steps = 1024
x0 = 0

R = np.array([])
XI= np.array([])


def computeInverse(xi):
    global x0, R, XI

    r = scipy.optimize.newton(F,x0,args=(xi,))
    x0 = r

    XI = np.append(XI, xi)
    R = np.append(R, r)

    return r

for i in range(steps):
    computeInverse(i / steps)

full = 1 - (1 / (steps)) / 4
computeInverse(full)

#print(R * 1/sigmaTr)


def sampleF(rnd):

    # Compute index
    i_xi = np.floor(rnd * steps).astype(int)

    # Find previous R
    r1 = R[i_xi]
    r2 = R[i_xi + 1]

    x1 = XI[i_xi]
    x2 = XI[i_xi + 1]

    # Linearly interpolate
    r = r1 + (r2 - r1) * (rnd - x1) / (x2 - x1)

    D = (2 * sigmaA) / (3*_sigmaT*_sigmaT)
    sigmaTr = math.sqrt(sigmaA/D)
    d = 1 / sigmaTr

    return r * d






# Example function to generate random samples
def generate_random_samples(size=10000):
    # Generates random samples from a uniform [0-1) distribution
    samples = np.random.rand(size)

    return sampleF(samples)

# Generate random samples
samples = generate_random_samples()

# Normalize y so it integrates to 1
sns.kdeplot(samples, color='g', label='KDE')

# Create an array of r values (e.g., from 0 to 10 with 1000 points)
r = np.linspace(0, 10, 1000)

# Calculate the corresponding f(r) values
f_r = Rdipole(r)

# Plot the function over the previous 
plt.plot(r, f_r, label='f(r)')

# Add labels and title
plt.xlabel('r')
plt.ylabel('f(r)')
plt.title('Plot of f(r)')
plt.legend()

# plot y log
#plt.yscale('log')

# Display the plot
plt.grid(True)
plt.show()