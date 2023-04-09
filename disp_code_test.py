import numpy as np
from matplotlib import pyplot as plt

#   real z_disp_giants(real plx, vector r, vector sunpos, real amplitude_z, real R_scale) {
  
#   real disp_mean;
#   vector[3] starpos;
#   real Rstar;

#   starpos = (1000.0/plx)*r + sunpos;
#   Rstar = sqrt(starpos[1]^2+starpos[2]^2);
#   disp_mean = amplitude_z * exp(-Rstar/R_scale);
#   return disp_mean;

def z_disp_giants(Rstar, amplitude_z, R_scale):
    return amplitude_z * np.exp( -Rstar/R_scale)

A = 808.69  # Uit de eerdere berekeningen
A = 80.869
R = np.arange(14000)
R_scale = 12e3  # Uit de slechte fit. Vervangen

disp_values = z_disp_giants(R, A, R_scale)

plt.figure(1)
plt.plot(R, disp_values)
plt.xlabel("$R_{star}$")
plt.ylabel("disp")
plt.ylim(0,90)
plt.show()