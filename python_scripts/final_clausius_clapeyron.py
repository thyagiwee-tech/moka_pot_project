import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import os

# 1. The 7 temperatures we simulated
temperatures = np.array([300.0, 315.0, 330.0, 350.0, 365.0, 380.0, 400.0])
pressures_kpa =[]
M_w = 0.018015 # kg/mol
R = 8.314 # J/mol K

# 2. Automatically read the GROMACS data for all 7 simulations
for T in temperatures:
    filename = f'../01_water_slab/density_{int(T)}K.xvg'
    if os.path.exists(filename):
        data = np.loadtxt(filename, comments=['@', '#'])
        z_axis, density = data[:, 0], data[:, 1]
        
        # Extract gas phase density (edges of the box)
        gas_density = np.mean(density[(z_axis < 2.0) | (z_axis > 13.0)])
        
        # Convert to pressure (kPa)
        P_kpa = ((gas_density / M_w) * R * T) / 1000.0
        pressures_kpa.append(P_kpa)
    else:
        print(f"Warning: Data for {int(T)}K not found!")

if len(pressures_kpa) > 0:
    pressures_kpa = np.array(pressures_kpa)
    pressures_pa = pressures_kpa * 1000

    # 3. Clausius-Clapeyron Math (1/T vs ln(P))
    x_inv_T = 1.0 / temperatures[:len(pressures_kpa)]
    y_ln_P = np.log(pressures_pa)

    # Perform Linear Regression (Line of Best Fit)
    slope, intercept, r_value, p_value, std_err = linregress(x_inv_T, y_ln_P)

    # Calculate dH_vap
    dH_vap_kJ_mol = (-slope * R) / 1000.0
    r_squared = r_value**2

    print("=======================================")
    print(f"Data Points: {len(pressures_kpa)}")
    for T, P in zip(temperatures[:len(pressures_kpa)], pressures_kpa):
        print(f"T = {T}K  -->  P = {P:.2f} kPa")
    print("---------------------------------------")
    print(f"Calculated dH_vap: {dH_vap_kJ_mol:.2f} kJ/mol")
    print(f"R-squared value:   {r_squared:.4f}")
    print("=======================================")

    # 4. Plotting the 7-Point Graph
    plt.figure(figsize=(9, 6))

    # Plot the raw data points
    plt.scatter(x_inv_T, y_ln_P, color='red', s=80, label='Simulated Data (MD)', zorder=5)

    # Plot the Line of Best Fit
    line_of_best_fit = slope * x_inv_T + intercept
    plt.plot(x_inv_T, line_of_best_fit, color='blue', linestyle='--', linewidth=2, label='Linear Regression')

    plt.title('Clausius-Clapeyron Plot (7 Data Points)', fontsize=14, fontweight='bold')
    plt.xlabel(r'1 / Temperature (K$^{-1}$)', fontsize=12)
    plt.ylabel('ln(Vapor Pressure in Pa)', fontsize=12)
    plt.grid(True, linestyle=':', alpha=0.8)
    plt.legend(fontsize=11)

    # Add a text box with the final physics results
    text_box = (f"$\\Delta H_{{vap}}$ = {dH_vap_kJ_mol:.2f} kJ/mol\n"
                f"$R^2$ = {r_squared:.4f}")
    plt.text(0.05, 0.15, text_box, transform=plt.gca().transAxes, fontsize=12,
             bbox=dict(facecolor='white', alpha=0.9, edgecolor='black'))

    plt.tight_layout()
    plt.savefig('../python_scripts/final_clausius_plot.png', dpi=300)
    print("Plot saved as 'final_clausius_plot.png' in the python_scripts folder!")
