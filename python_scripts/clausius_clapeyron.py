import numpy as np
import matplotlib.pyplot as plt

# Your extracted simulation data
temperatures = np.array([300.0, 350.0])
pressures_kpa = np.array([8.5921, 86.0148]) 

# Convert to standard units for math
R = 8.314 # J / mol K
pressures_pa = pressures_kpa * 1000

# X and Y axis for the Clausius-Clapeyron plot
x_inv_T = 1.0 / temperatures
y_ln_P = np.log(pressures_pa)

# Calculate slope (m)
slope = (y_ln_P[1] - y_ln_P[0]) / (x_inv_T[1] - x_inv_T[0])

# Calculate dH_vap (Slope = -dH_vap / R)
dH_vap_J_mol = -slope * R
dH_vap_kJ_mol = dH_vap_J_mol / 1000.0

print("=======================================")
print(f"Calculated Enthalpy of Vaporization (dH_vap): {dH_vap_kJ_mol:.2f} kJ/mol")
print(f"Textbook value for Water: 40.65 kJ/mol")
error = abs((dH_vap_kJ_mol - 40.65) / 40.65) * 100
print(f"Percentage Error: {error:.2f}%")
print("=======================================")

# Plotting the Graph
plt.figure(figsize=(8, 5))
plt.plot(x_inv_T, y_ln_P, marker='o', color='red', linestyle='-', linewidth=2, markersize=8)
plt.title('Clausius-Clapeyron Plot from MD Simulation', fontsize=14)
plt.xlabel('1 / Temperature (K$^{-1}$)', fontsize=12)
plt.ylabel('ln(Vapor Pressure)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)

# Add text box with the result
plt.text(0.0031, 9.5, f'$\Delta H_{{vap}}$ = {dH_vap_kJ_mol:.2f} kJ/mol', 
         fontsize=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='black'))

plt.savefig('clausius_clapeyron_plot.png', dpi=300)
print("Plot saved as 'clausius_clapeyron_plot.png'")

