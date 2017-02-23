import numpy as np
import matplotlib.pyplot as plt

# BEAM
rho = 265.e-10              # resistivity of the beam (Ohm-meter)
mu_0 = 4.e-7 * np.pi        # vacuum permeability (N/A^2)
h_0 = 25.e-3                # thickness of beam (m); first and last 61m of track
h = 15.e-3                  # thickness of beam (m); rest of track
z_0 = 10.e-3                # distance from beam (m)
# MAGNETS
V = 0.0254 ** 3             # volume of neodymium magnets (m^3)
M = 7580.e3                 # magnetization of neodymium magnets (A/m)
m = V * M                   # magnetic vertical dipole moment (A-m^2)
w = 2 * rho / (mu_0 * h)    # velocity of magnetic propagation through beam (m/s)
# WHEEL


F_lift = lambda v: ((3 * mu_0 * m ** 2) / (32 * np.pi * z_0 ** 4)) * (1 - (w / np.sqrt(v ** 2 + w ** 2)))
F_drag = lambda v: (w / v) * F_lift(v)

# Compute the x and y coordinates for points on sine and cosine curves
x = np.arange(0.001, 30)
F_lift_plot = F_lift(x)
F_drag_plot = F_drag(x)

# Plot the points using matplotlib
plt.plot(x, F_lift_plot)
plt.plot(x, F_drag_plot)
plt.xlabel('Velocity (m/s)')
plt.ylabel('Force (N)')
plt.title('Force vs Velocity')
plt.legend(['Lift Force', 'Drag Force'])
plt.show()