# Glycolysis

Package for modeling mammalian glycolysis activity

You can use it to simulate glycolysis using the following code:
```julia
using Glycolysis, DifferentialEquations, CairoMakie

# Set ATPase rate to be equal to 10% of glycolysis Vmax
model_params.ATPase_Vmax = 0.1 * 2 * model_params.HK1_Vmax * model_params.HK1_Conc
prob = ODEProblem(glycolysis_ODEs, initial_concentrations, (0,100), model_params)

#Simulate glycolysis activity
sol = solve(
    prob,
    Rodas4(),
    abstol = 1e-12,
    reltol = 1e-5,
)

#Plot [ATP] over time
lines(sol.t, [sol.u[i].ATP for i in 1:length(sol)], axis=(limits=((nothing),(0,10e-3)),))

![image](https://user-images.githubusercontent.com/75404066/205253830-57abfa66-d48e-4025-98ad-d3f19980f197.png)
```
