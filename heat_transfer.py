import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

# User parameters
flibe_velocity = 2.0  # m/s
flibe_inlet_temp = 900  # K
pipe_diameter = 0.02  # m (typical channel width)

# Domain setup
dx = 1e-4
lengths = {
    'w': 0.003,
    'v1': 0.01,
    'flibe1': 0.02,
    'v2': 0.03,
    'flibe2': 1.0,
    'v3': 0.03
}
regions = ['w', 'v1', 'flibe1', 'v2', 'flibe2', 'v3']
materials = {'w': 'tungsten', 'v1': 'vcrti', 'flibe1': 'flibe',
             'v2': 'vcrti', 'flibe2': 'flibe', 'v3': 'vcrti'}

# Material properties
props = {
    'tungsten': {'k': 170, 'cp': 134, 'rho': 19250, 'q_gen': 1e6},
    'vcrti': {'k': 25, 'cp': 500, 'rho': 6000, 'q_gen': 1e6},
    'flibe': {
    'k': 1.0,             # W/m-K
    'cp': 2410.0,         # J/kg-K
    'rho': 1940.0,        # kg/m³
    'mu': 6.0e-3,         # Pa·s
    'q_gen': 1e6
}
}

# Convection coefficient (Dittus-Boelter)
Re = props['flibe']['rho'] * flibe_velocity * pipe_diameter / props['flibe']['mu']
Pr = props['flibe']['cp'] * props['flibe']['mu'] / props['flibe']['k']
Nu = 0.023 * Re**0.8 * Pr**0.4
h = Nu * props['flibe']['k'] / pipe_diameter
print(f"Convective heat transfer coefficient h = {h:.2f} W/m²-K")

# Boundary temperatures
T_left = 1000  # K
T_right = 800  # K

# Grid setup
x = []
k = []
cp = []
rho = []
q_gen = []
region_labels = []

for r in regions:
    n = int(lengths[r] / dx)
    mat = materials[r]
    x.extend(np.linspace(len(x)*dx, len(x)*dx + lengths[r], n, endpoint=False))
    k.extend([props[mat]['k']] * n)
    cp.extend([props[mat]['cp']] * n)
    rho.extend([props[mat]['rho']] * n)
    q_gen.extend([props[mat]['q_gen']] * n)
    region_labels.extend([r] * n)

x = np.array(x)
N = len(x)
T = np.ones(N) * flibe_inlet_temp
T[0] = T_left
T[-1] = T_right

# Time integration (transient FDM)
dt = 0.01
time = 20
nt = int(time / dt)

for tstep in range(nt):
    T_new = T.copy()
    for i in range(1, N-1):
        alpha = k[i] / (rho[i] * cp[i])
        # Standard conduction update
        conduction = alpha * dt / dx**2 * (T[i+1] - 2*T[i] + T[i-1])
        source = dt * q_gen[i] / (rho[i] * cp[i])
        T_new[i] = T[i] + conduction + source

        # Convection at FLiBe-metal interfaces
        if region_labels[i] == 'flibe1' and region_labels[i+1] != 'flibe1':
            # Right edge of flibe1
            T_new[i] += dt * h * (flibe_inlet_temp - T[i]) / (rho[i] * cp[i] * dx)
        elif region_labels[i] == 'flibe2' and region_labels[i+1] != 'flibe2':
            # Right edge of flibe2
            T_new[i] += dt * h * (flibe_inlet_temp - T[i]) / (rho[i] * cp[i] * dx)
        elif region_labels[i] != 'flibe1' and region_labels[i+1] == 'flibe1':
            # Left edge of flibe1
            T_new[i+1] += dt * h * (flibe_inlet_temp - T[i+1]) / (rho[i+1] * cp[i+1] * dx)
        elif region_labels[i] != 'flibe2' and region_labels[i+1] == 'flibe2':
            # Left edge of flibe2
            T_new[i+1] += dt * h * (flibe_inlet_temp - T[i+1]) / (rho[i+1] * cp[i+1] * dx)

    # Apply boundary conditions
    T_new[0] = T_left
    T_new[-1] = T_right
    T = T_new.copy()

# Plot result
plt.plot(x, T)
plt.xlabel('x [m]')
plt.ylabel('Temperature [K]')
plt.title('1D Temperature Profile with Convection')
plt.grid()
plt.tight_layout()
plt.show()