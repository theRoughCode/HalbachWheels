import numpy as np
import matplotlib.pyplot as plt

# BEAM
rho = 265.e-10              # resistivity of the beam (Ohm-meter)
mu_0 = 4.e-7 * np.pi        # vacuum permeability (N/A^2)
h_0 = 25.e-3                # thickness of beam (m); first and last 61m of track
h = 15.e-3                  # thickness of beam (m); rest of track
z_0 = 10.e-3                # distance from beam (m)
# MAGNETS
sets = 5                    # number of halbach sets (4 magnets)
num_magnets = sets * 4      # total number of magnets
spacing = 0.0127            # space between each magnet (m)
side_len = 0.0254           # side length of each magnet (m)
V = side_len ** 3           # volume of neodymium magnets (m^3)
M = 7580.e3                 # magnetization of neodymium magnets (A/m)
m = V * M                   # magnetic vertical dipole moment (A-m^2)
w = 2 * rho / (mu_0 * h)    # velocity of magnetic propagation through beam (m/s)
# WHEEL
r = ((side_len + spacing) * num_magnets) / (2 * np.pi)  # radius of wheel


# CONVERSIONS
rpm_to_angular_vel = lambda rpm: 30 * rpm / np.pi
tangential_vel = lambda angular_vel: r * angular_vel    # tangential velocity of wheel given
rpm_to_tan_vel = lambda rpm: tangential_vel(rpm_to_angular_vel(rpm))

# FORCES (RPM -> Force)
F_lift = lambda rpm: ((3 * mu_0 * m ** 2) / (32 * np.pi * z_0 ** 4)) * (1 - (w / np.sqrt(rpm_to_tan_vel(rpm) ** 2 + w ** 2)))
F_drag = lambda rpm: (w / rpm_to_tan_vel(rpm)) * F_lift(rpm)

def getPOI(f, g, x0):
    tolerance = 0.0001
    h = lambda x: f(x) - g(x)
    dh = lambda x : (h(x+tolerance) - h(x)) / tolerance
    x1 = x0 - h(x0)/dh(x0)
    if h(x1) > tolerance:
        return getPOI(f, g, x1)
    else:
        return x1


# Compute the x and y coordinates for points on sine and cosine curves
x = np.arange(0.001, 30)
F_lift_plot = F_lift(x)
F_drag_plot = F_drag(x)
max_drag_x = np.argmax(F_drag_plot)
max_drag_y = F_drag(max_drag_x)
POI_x = getPOI(F_lift, F_drag, 2)
POI_y = F_drag(POI_x)

# Plot the points using matplotlib
plt.plot(x, F_lift_plot)
plt.plot(x, F_drag_plot)
plt.plot(POI_x, POI_y, 'ro')
plt.plot(max_drag_x, max_drag_y, 'bD')
plt.xlabel('RPM (cycles/min)')
plt.ylabel('Force (N)')
plt.title('Force vs RPM')
plt.legend(['Lift Force', 'Drag Force'])
#plt.annotate('POI = (' + str('{0:.2f}'.format(POI_x)) + ', ' + str('{0:.2f}'.format(POI_y)) + ')', xy=(2, 1), xytext=(0.5, POI_y))
#plt.annotate('Max F_drag = (' + str('{0:.2f}'.format(max_drag_x)) + ', ' + str('{0:.2f}'.format(max_drag_y)) + ')',
 #            xy=(2, 1), xytext=(max_drag_x + 0.1, max_drag_y))
plt.show()