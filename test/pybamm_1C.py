import pybamm
import pandas as pd

# Load standard isothermal SPMe and Chen2020 parameters
model = pybamm.lithium_ion.SPMe()
param = pybamm.ParameterValues("Chen2020")

# Define a 1C discharge (Chen2020 nominal capacity is exactly 5 Ah, so 1C = 5A)
experiment = pybamm.Experiment([
    "Discharge at 5 A until 2.5 V"
])

# Run simulation
sim = pybamm.Simulation(model, parameter_values=param, experiment=experiment)
sim.solve()

# Extract variables
t = sim.solution["Time [s]"].entries
v = sim.solution["Terminal voltage [V]"].entries

# Export to CSV
df = pd.DataFrame({"Time [s]": t, "Voltage [V]": v})
df.to_csv("pybamm_spme_1c.csv", index=False)
print("PyBaMM baseline saved to pybamm_spme_1c.csv")