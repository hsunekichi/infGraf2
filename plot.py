import numpy as np
import matplotlib.pyplot as plt
import math as Math

# Define the function f(r)
def f(r):
    eta = 1.5
    r2 = r*r
    sigmaS = 2.6
    sigmaA = 0.0041
    g = 0

    # Compute isotropic phase function
    _sigmaS = (1 - g) * sigmaS
    _sigmaT = _sigmaS + sigmaA
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


# Create an array of r values (e.g., from 0 to 10 with 1000 points)
r = np.linspace(0, 10, 1000)

# Calculate the corresponding f(r) values
f_r = f(r)

# Plot the function
plt.plot(r, f_r, label='f(r) = sin(r)')  # Replace label with your actual function

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