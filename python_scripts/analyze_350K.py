import numpy as np
import matplotlib.pyplot as plt

# 1. Load the GROMACS data (skipping the @ and # header lines)
data = np.loadtxt('density_350K.xvg', comments=['@', '#'])
z_axis = data[:, 0]  # Distance along the Z-axis (nanometers)
density = data[:, 1] # Density (kg/m^3)

# 2. Extract the Gas Phase Density (The vacuum edges of our box)
# We look at the first 2 nm and the last 2 nm of the box where the steam is
vapor_region = (z_axis < 2.0) | (z_axis > 13.0)
gas_density_kg_m3 = np.mean(density[vapor_region])

# 3. Calculate Vapor Pressure using the Ideal Gas Law (P = (density / Molar Mass) * R * T)
R = 8.314 # Ideal gas constant (J / mol K)
T = 350.0 # Temperature (Kelvin)
M_w = 0.018015 # Molar mass of water (kg / mol)

# Calculate Pressure in Pascals, then convert to kiloPascals (kPa)
pressure_pa = (gas_density_kg_m3 / M_w) * R * T
pressure_kpa = pressure_pa / 1000.0

print("=======================================")
print(f"Simulation Temperature: {T} K")
print(f"Calculated Gas Density: {gas_density_kg_m3:.4f} kg/m^3")
print(f"Calculated Vapor Pressure: {pressure_kpa:.4f} kPa")
print("=======================================")

# 4. Plot the Liquid Slab Profile
plt.figure(figsize=(8, 5))
plt.plot(z_axis, density, color='blue', linewidth=2, label='Water Density')
plt.fill_between(z_axis, density, color='blue', alpha=0.3)

plt.title('Density Profile of the Water Slab (350K)', fontsize=14)
plt.xlabel('Z-axis position (nm)', fontsize=12)
plt.ylabel('Density (kg / m$^3$)', fontsize=12)
plt.axhline(y=gas_density_kg_m3, color='red', linestyle='--', label=f'Gas Density (avg)')
plt.legend()
plt.grid(True, alpha=0.5)

# Save the plot as an image
plt.savefig('slab_profile_350K.png', dpi=300)
print("Plot saved as 'slab_profile_350K.png'")

