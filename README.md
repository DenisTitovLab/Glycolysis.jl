# Glycolysis

Package for modeling mammalian glycolysis activity

'''julia
using Glycolysis, DifferentialEquations

prob = ODEProblem(glycolysis_ODEs, initial_concentrations, (0,10), model_params)
sol = solve(
    prob,
    Rodas4(),
    abstol = 1e-12,
    reltol = 1e-5,
)
'''
