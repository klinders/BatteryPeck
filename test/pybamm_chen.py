# pybamm_chen.py
# Compare PyBaMM SPMe simulation to experimental Chen dataset across 4 segments

import pybamm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Load experimental data
df = pd.read_csv("data/Chen2020/LGM50_cell03.csv", skiprows=13)

# STRIP hidden spaces from all column names to prevent KeyErrors
df.columns = df.columns.str.strip()

# Process current
df["Md"] = df["Md"].astype(str).str.strip()
df["Processed Current [A]"] = 0.0
df.loc[df["Md"] == "C", "Processed Current [A]"] = -df["Current [A]"]
df.loc[df["Md"] == "D", "Processed Current [A]"] = df["Current [A]"]

# Load lookup table and create interpolator
lut = pd.read_csv("data/Chen2020/soc_ocv_lut.csv")
v_to_soc = interp1d(lut["Voltage"], lut["SoC"], fill_value="extrapolate")

def extract_discharge_segments(df_in):
    # Extract isolated R to D to R segments from dataframe
    segments = []
    n = len(df_in)
    i = 0
    md = df_in["Md"].values
    
    while i < n:
        if md[i] == "D":
            # Find start of preceding R block
            start_idx = i
            while start_idx > 0 and md[start_idx - 1] == "R":
                start_idx -= 1
                
            # Find end of D block
            end_idx = i
            while end_idx < n - 1 and md[end_idx + 1] == "D":
                end_idx += 1
                
            # Find end of succeeding R block
            while end_idx < n - 1 and md[end_idx + 1] == "R":
                end_idx += 1
                
            segments.append(df_in.iloc[start_idx:end_idx + 1].copy())
            i = end_idx
        i += 1
        
    # Ignore segments that do not start at full charge
    return [s for s in segments if s["Voltage [V]"].iloc[0] >= 4.1]

segments = extract_discharge_segments(df)
print(f"Found {len(segments)} discharge segments to evaluate.\n")

# Set up matplotlib figure
fig, (ax_volt, ax_curr) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
ax_volt.set_ylabel("Voltage (V)")
ax_curr.set_ylabel("Current (A)")
ax_curr.set_xlabel("Time (s)")

colours = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"]

for idx, seg in enumerate(segments):
    print(f"Simulating segment {idx}...")
    
    # Extract segment arrays
    t_exp = seg["Test Time [s]"].values
    t_sim = t_exp - t_exp[0]
    v_exp = seg["Voltage [V]"].values
    i_exp = seg["Processed Current [A]"].values
    
    # Read starting voltage to set accurate initial SOC
    v_init = v_exp[0]
    soc_init = float(v_to_soc(v_init))
    
    # Build standard isothermal SPMe
    model = pybamm.lithium_ion.SPMe()
    param = pybamm.ParameterValues("Chen2020")

    # Decrease voltage cut-off to prevent early termination
    param["Lower voltage cut-off [V]"] = 2.0
    
    # Feed experimental current profile into PyBaMM
    current_interpolant = pybamm.Interpolant(t_sim, i_exp, pybamm.t)
    param["Current function [A]"] = current_interpolant
    
    # Run simulation
    sim = pybamm.Simulation(model, parameter_values=param)
    
    try:
        # Solve using experimental time steps
        sol = sim.solve(t_eval=t_sim, initial_soc=soc_init)
        
        # Extract raw solver outputs
        t_sol = sol["Time [s]"].entries
        v_sol = sol["Terminal voltage [V]"].entries
        
        # Interpolate PyBaMM results back onto the exact experimental time grid
        # This handles both shape mismatches and early solver termination
        v_interp_func = interp1d(t_sol, v_sol, bounds_error=False, fill_value=np.nan)
        v_sim_mapped = v_interp_func(t_sim)
        
        # Filter out NaNs in case the PyBaMM simulation terminated early (e.g. hit 2.5V limit)
        valid = ~np.isnan(v_sim_mapped)
        v_sim_clean = v_sim_mapped[valid]
        v_exp_clean = v_exp[valid]
        
        # Calculate metrics on the overlapping region
        rmse_v = np.sqrt(np.mean((v_sim_clean - v_exp_clean)**2))
        ss_res_v = np.sum((v_exp_clean - v_sim_clean)**2)
        ss_tot_v = np.sum((v_exp_clean - np.mean(v_exp_clean))**2)
        r2_v = 1.0 - (ss_res_v / ss_tot_v) if ss_tot_v != 0 else 0.0
        
        # Print metrics
        print(f"Segment {idx} Results:")
        print(f"  Starting voltage | {v_init:.4f} V")
        print(f"  Voltage          | RMSE: {rmse_v * 1000:.2f} mV | R²: {r2_v:.4f}")
        if not valid.all():
            print(f"  Note: PyBaMM terminated early at {t_sol[-1]:.1f}s (Exp segment length: {t_sim[-1]:.1f}s)")
        print()
        
        # Plot experimental traces
        t_shift = t_exp[0]
        c = colours[idx % len(colours)]
        
        ax_volt.plot(t_exp, v_exp, color="black", linestyle="--", 
                     label="Exp. voltage" if idx == 0 else "")
        ax_curr.plot(t_exp, i_exp, color="black", linestyle="--", 
                     label="Exp. current" if idx == 0 else "")
        
        # Plot simulation traces (using the raw solver time for the plot)
        ax_volt.plot(t_sol + t_shift, v_sol, color=c, 
                     label="PyBaMM voltage" if idx == 0 else "")
        ax_volt.scatter(t_shift, v_init, color="black", marker="x", zorder=5)
        
    except pybamm.SolverError as e:
        print(f"Warning: Simulation {idx} failed to solve. {e}\n")

ax_volt.legend(loc="lower left")
ax_curr.legend(loc="lower left")
plt.suptitle("PyBaMM SPMe validation", fontsize=14)
plt.tight_layout()
plt.show()