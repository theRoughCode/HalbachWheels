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
M = 1138.e3                 # magnetization of neodymium magnets (A/m)
m = V * M                   # magnetic vertical dipole moment (A-m^2)
w = 2 * rho / (mu_0 * h)    # velocity of magnetic propagation through beam (m/s)
# WHEEL
r = ((side_len + spacing) * num_magnets) / (2 * np.pi)  # radius of wheel

# CONVERSIONS
rpm_to_angular_vel = lambda rpm: np.pi * rpm / 30
tangential_vel = lambda angular_vel: r * angular_vel    # tangential velocity of wheel given
rpm_to_tan_vel = lambda rpm: tangential_vel(rpm_to_angular_vel(rpm))

# FORCES (v_relative -> Force)
F_lift_v = lambda v: ((3 * mu_0 * m ** 2) / (32 * np.pi * z_0 ** 4)) * (1 - (w / np.sqrt(v ** 2 + w ** 2)))
F_drag_v = lambda v: (w / v) * F_lift_v(v)

# FORCES (RPM -> Force)
F_lift_rpm = lambda rpm: ((3 * mu_0 * m ** 2) / (32 * np.pi * z_0 ** 4)) * (1 - (w / np.sqrt(rpm_to_tan_vel(rpm) ** 2 + w ** 2)))
F_drag_rpm = lambda rpm: (w / rpm_to_tan_vel(rpm)) * F_lift_rpm(rpm)


def Newton(f, x0, tolerance, derivative=False):
    if derivative:
        h = lambda x: abs((f(x + tolerance) - f(x)) / tolerance)
    else:
        h = f
    dh = lambda x: (h(x + tolerance) - h(x)) / tolerance
    x1 = x0 - h(x0) / dh(x0)
    if h(x1) > tolerance:
        return Newton(f, x1, tolerance, derivative)
    else:
        return x1


def getPOI(f, g):
    tolerance = 0.0001
    x0 = 250
    return Newton(lambda x: abs(f(x) - g(x)), x0, tolerance)


def getMax(f):
    tolerance = 0.01
    x0 = 270
    return Newton(f, x0, tolerance, True)


# CONTROLS
display_points = True       # display POI and max/min
display_labels = False      # display labels for points
ZOOM = False                # zoom in to key points
F_var = "rpm"               # v or rpm
START = 0.001
END = 2000                  # RPM (default)
if ZOOM:
    END /= 16
if F_var == "v":
    END = rpm_to_tan_vel(END)
STEP = 0.1

# Compute the x and y coordinates
x = np.arange(START, END, STEP)
if F_var == "v":
    F_lift = F_lift_v
    F_drag = F_drag_v
else:
    F_lift = F_lift_rpm
    F_drag = F_drag_rpm
F_lift_plot = F_lift(x)
F_drag_plot = F_drag(x)
if display_points:
    max_drag_x = getMax(F_drag)
    max_drag_y = F_drag(max_drag_x)
    POI_x = getPOI(F_lift, F_drag)
    POI_y = F_drag(POI_x)

# Plot the points using matplotlib
plt.plot(x, F_lift_plot)
plt.plot(x, F_drag_plot)
plt.xlabel('RPM (cycles/min)' if F_var == "rpm" else 'V (m/s)')
plt.ylabel('Force (N)')
plt.title('Force vs RPM' if F_var == "rpm" else 'Force vs V')
plt.legend(['Lift Force', 'Drag Force'])
if display_points:
    plt.plot(POI_x, POI_y, 'ro')
    plt.plot(max_drag_x, max_drag_y, 'bD')
    if display_labels:
        plt.annotate('POI = (' + str('{0:.2f}'.format(POI_x)) + ', ' + str('{0:.2f}'.format(POI_y)) + ')', xy=(2, 1), xytext=(POI_x - 150, POI_y))
        plt.annotate('Max F_drag = (' + str('{0:.2f}'.format(max_drag_x)) + ', ' + str('{0:.2f}'.format(max_drag_y)) + ')',
                     xy=(2, 1), xytext=(max_drag_x - 30, max_drag_y + 20))
plt.show()