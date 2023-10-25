- The file "Modified_VLM_Parameters.py" contains all the geometry and propellant properties for the thruster. 
- "Density_model.py" contains an estimation for the density of water vapor based ont he steam tables which is called within the file "Modified_VLM_model_ver3.py".

- "Modified_VLM_model_ver3.py" runs the 1D model for the heating chamber of the thruster. It is fed the input operating conditions along with the geometry and properties of the propellant. It is then called in the script "Modified_VLM_nozzle_model".

- "Modified_VLM_nozzle_model" runs the model for the nozzle section and returns the output parameters such as thrust, Isp, etc. 

- "Coupled_model.py" calls both the heating chmaber and nozzle model functions iteratively until the mass flow rate converges (It couples the two models).

- "op_conds_sensitivity_HC.py" calls the coupled model function for a range of input operating conditions and saves the output data in respective files.

# Run only the file named "op_conds_sensitivity_HC.py" 