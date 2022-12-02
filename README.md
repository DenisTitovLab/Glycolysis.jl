# Glycolysis

Package for modeling mammalian glycolysis activity

You can use it to simulate glycolysis using the following code:
```julia
using Glycolysis, DifferentialEquations, CairoMakie, Revise

# Set ATPase rate to be equal to 10% of glycolysis Vmax which is set by HK1 Vmax
model_params.ATPase_Vmax = 0.1 * 2 * model_params.HK1_Vmax * model_params.HK1_Conc
prob = ODEProblem(glycolysis_ODEs, initial_concentrations, (0, 100), model_params)

#Simulate glycolysis activity
sol = solve(prob, Rodas4(), abstol = 1e-12, reltol = 1e-5)

#Plot all [Metabolite] over time
fig = Figure()
ax = Axis(fig[1, 1], yscale = log10, xlabel = "Time, min", ylabel = "[Metabolite], M")
for metabolite in [:ATP, :ADP, :Phosphate]
    lines!(ax, sol.t, sol[metabolite, :], label = String(metabolite))
end
Legend(fig[1, 2], ax)
fig
```
<img src="https://user-images.githubusercontent.com/75404066/205350460-f73de619-14c3-42a4-95d7-7fa0d6e077ee.png" width="500">
