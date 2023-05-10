# Glycolysis

Package for modeling mammalian glycolysis activity.  
See our [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2022.12.28.522046) for examples of how to use this package to study glycolysis regulation.  
Code to reproduce all figures in the preprint are in _examples/biorxiv figures_ folder.  

You can use it to simulate glycolysis using the following code:
```julia
using Glycolysis, DifferentialEquations, CairoMakie

# Set ATPase rate to be equal to 10% of glycolysis Vmax which is set by HK1 Vmax
glycolysis_params.ATPase_Vmax = 0.1 * 2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc
prob = ODEProblem(glycolysis_ODEs, glycolysis_init_conc, (0, 100), glycolysis_params)

#Simulate glycolysis activity
sol = solve(prob, Rodas4(), abstol = 1e-12, reltol = 1e-5)

#List of all metabolites in the model can be found by running
propertynames(glycolysis_init_conc)

#Plot all [Metabolite] over time
fig = Figure()
ax = Axis(fig[1, 1], yscale = log10, xlabel = "Time, min", ylabel = "[Metabolite], M")
#Substitute :ATP, :ADP, :Phosphate for any metabolite(s) from propertynames(glycolysis_init_conc)
for metabolite in [:ATP, :ADP, :Phosphate]
    lines!(ax, sol.t, sol[metabolite, :], label = String(metabolite))
end
Legend(fig[1, 2], ax)
fig
```
<img src="https://user-images.githubusercontent.com/75404066/205350460-f73de619-14c3-42a4-95d7-7fa0d6e077ee.png" width="500">