using Glycolysis
using DifferentialEquations, ProgressMeter
using CairoMakie, Dates, Printf, Statistics, StatsBase
using DataFrames, CSV

##
# Precalculate output of complete model and model without regulation
function find_ATP_at_ATPase_range(glycolysis_params,
        glycolysis_init_conc;
        n_Vmax_ATPase_values = 1000)
    tspan = (0.0, 1e8)
    pathway_Vmax = 2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc
    ATPases = 10 .^ range(log10(0.003), log10(1.0), n_Vmax_ATPase_values) .* pathway_Vmax
    prob = ODEProblem(glycolysis_ODEs, glycolysis_init_conc, tspan, glycolysis_params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(ensemble_prob,
        Rodas5P(),
        EnsembleThreads(),
        trajectories = n_Vmax_ATPase_values,
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
        save_start = false)
    ATP_prod_rate = [Glycolysis.conc_to_rates(sol.u[end], sol.prob.p).ATPprod /
                     sol.prob.p.ATPase_Vmax for sol in sim if sol.retcode == ReturnCode.Success]
    ATP_energy = -log.([Glycolysis.conc_to_disequilibrium_ratios(sol.u[end],
        sol.prob.p).Q_Keq_ATPase for sol in sim if sol.retcode == ReturnCode.Success])
    ATP_conc = [sol.u[end].ATP for sol in sim if sol.retcode == ReturnCode.Success]
    return (ATP_conc = ATP_conc,
        ATP_prod_rate = ATP_prod_rate,
        ATP_energy = ATP_energy,
        ATPase_Vmax = ATPases /
                      (2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc))
end

glycolysis_params.ATPase_Km_ATP = 1e-9
glycolysis_params_copy = deepcopy(glycolysis_params)
Complete_Model_Simulation_Data = @time find_ATP_at_ATPase_range(glycolysis_params_copy,
    glycolysis_init_conc)

no_reg_params = deepcopy(glycolysis_params)
no_reg_params.HK1_K_a_G6P_cat = Inf
no_reg_params.HK1_K_i_G6P_reg = Inf
no_reg_params.HK1_K_a_Pi = Inf
no_reg_params.PFKP_L = 0.0
no_reg_params.GAPDH_L = 0.0
no_reg_params.PKM2_L = 0.0
no_reg_params.ATPase_Km_ATP = 1e-9

No_Reg_Model_Simulation_Data = find_ATP_at_ATPase_range(no_reg_params, glycolysis_init_conc)

tspan = (0.0, 240.0 / 4)
Initial_ATPase_Vmax_frac = 0.03
High_ATPase_Vmax_frac = Initial_ATPase_Vmax_frac * 2
Low_ATPase_Vmax_frac = Initial_ATPase_Vmax_frac / 2
ATPase_change_time1 = 60 / 4
ATPase_change_time2 = 120 / 4
ATPase_change_time3 = 180 / 4
glycolysis_params.ATPase_Vmax = Initial_ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc *
                                glycolysis_params.HK1_Vmax
glycolysis_params.ATPase_Km_ATP = 1e-9

function affect1!(integrator)
    integrator.p.ATPase_Vmax = High_ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc *
                               glycolysis_params.HK1_Vmax
end
function affect2!(integrator)
    integrator.p.ATPase_Vmax = Low_ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc *
                               glycolysis_params.HK1_Vmax
end
function affect3!(integrator)
    integrator.p.ATPase_Vmax = Initial_ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc *
                               glycolysis_params.HK1_Vmax
end

PresetTime_cb1 = PresetTimeCallback(ATPase_change_time1, affect1!)
PresetTime_cb2 = PresetTimeCallback(ATPase_change_time2, affect2!)
PresetTime_cb3 = PresetTimeCallback(ATPase_change_time3, affect3!)
cb_set = CallbackSet(PresetTime_cb1, PresetTime_cb2, PresetTime_cb3)

init_cond_prob = ODEProblem(glycolysis_ODEs,
    glycolysis_init_conc,
    (0, 1e8),
    glycolysis_params)
init_cond_sol = solve(init_cond_prob,
    Rodas5P(),
    abstol = 1e-15,
    reltol = 1e-8,
    save_everystep = false)
new_init_cond = init_cond_sol.u[end]
prob = ODEProblem(glycolysis_ODEs,
    new_init_cond,
    tspan,
    glycolysis_params,
    callback = cb_set)
sol = solve(prob,
    Rodas5P(),
    abstol = 1e-15,
    reltol = 1e-8,
    saveat = [k for k in tspan[1]:((tspan[2] - tspan[1]) / 10_000):tspan[2]])

timepoints = sol.t
ATPprod = [Glycolysis.conc_to_rates(conc, glycolysis_params).ATPprod for conc in sol.u] ./
          (2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax)
ATPenergy = [Glycolysis.conc_to_disequilibrium_ratios(conc, glycolysis_params).Q_Keq_ATPase
             for conc in sol.u]
ATPase = Initial_ATPase_Vmax_frac * ones(length(sol))
ATPase[sol.t .> ATPase_change_time1 .&& sol.t .< ATPase_change_time2] .= High_ATPase_Vmax_frac
ATPase[sol.t .> ATPase_change_time2 .&& sol.t .< ATPase_change_time3] .= Low_ATPase_Vmax_frac
ATP = [conc.ATP for conc in sol.u]

# Complete_Model_Simulation_Data = find_ATP_at_ATPase_range(glycolysis_params, glycolysis_init_conc)

no_reg_glycolysis_params = deepcopy(glycolysis_params)
no_reg_glycolysis_params.HK1_K_a_G6P_cat = Inf
no_reg_glycolysis_params.HK1_K_i_G6P_reg = Inf
no_reg_glycolysis_params.HK1_K_a_Pi = Inf
no_reg_glycolysis_params.PFKP_L = 0.0
no_reg_glycolysis_params.GAPDH_L = 0.0
no_reg_glycolysis_params.PKM2_L = 0.0

no_reg_glycolysis_params.ATPase_Vmax = Initial_ATPase_Vmax_frac * 2 *
                                       no_reg_glycolysis_params.HK1_Conc *
                                       no_reg_glycolysis_params.HK1_Vmax
no_reg_glycolysis_params.ATPase_Km_ATP = 1e-9

no_reg_init_cond_prob = ODEProblem(glycolysis_ODEs,
    glycolysis_init_conc,
    (0, 1e8),
    no_reg_glycolysis_params)
no_reg_init_cond_sol = solve(no_reg_init_cond_prob,
    Rodas5P(),
    abstol = 1e-15,
    reltol = 1e-8,
    save_everystep = false)
no_reg_new_init_cond = no_reg_init_cond_sol.u[end]
no_reg_prob = ODEProblem(glycolysis_ODEs,
    no_reg_new_init_cond,
    tspan,
    no_reg_glycolysis_params,
    callback = cb_set)
no_reg_sol = solve(no_reg_prob,
    Rodas5P(),
    abstol = 1e-15,
    reltol = 1e-8,
    saveat = [k for k in tspan[1]:((tspan[2] - tspan[1]) / 10_000):tspan[2]])
no_reg_timepoints = no_reg_sol.t
no_reg_ATPprod = [Glycolysis.conc_to_rates(conc, no_reg_glycolysis_params).ATPprod
                  for conc in no_reg_sol.u] /
                 (2 * no_reg_glycolysis_params.HK1_Conc * no_reg_glycolysis_params.HK1_Vmax)
no_reg_ATPenergy = [Glycolysis.conc_to_disequilibrium_ratios(conc,
    glycolysis_params).Q_Keq_ATPase for conc in no_reg_sol.u]
no_reg_ATPase = Initial_ATPase_Vmax_frac * ones(length(no_reg_sol))
no_reg_ATPase[no_reg_sol.t .> ATPase_change_time1 .&& no_reg_sol.t .< ATPase_change_time2] .= High_ATPase_Vmax_frac
no_reg_ATPase[no_reg_sol.t .> ATPase_change_time2 .&& no_reg_sol.t .< ATPase_change_time3] .= Low_ATPase_Vmax_frac
no_reg_ATP = [conc.ATP for conc in no_reg_sol.u]

##
# Plot the results
size_inches = (6.5, 5)
size_pt = 72 .* size_inches
set_theme!(Theme(fontsize = 6,
    Axis = (xticksize = 1,
        yticksize = 1,
        # xticklabelsize = 6,
        # yticklabelsize = 6,
        yticklabelpad = 1,
        ylabelpadding = 3)))
fig = Figure(size = size_pt)

# full_color = :black
# no_reg_color = :red

full_color = :Black
no_reg_color = :Red
line_ATPase_width = 2
line_ATPase_color = :grey
no_reg_linestyle = :dash
line_models_width = 4

# Plot ATP production
ax_ATP_prod = Axis(fig[1, 1],
    limits = (nothing, (0.0, 1.5 * maximum(ATPase))),
    xlabel = "Time, min",
    ylabel = rich("ATP production rate, relative to glycolysis V", subscript("max")),
    title = "Dynamically matching\nATP supply and demand")
ax_ATPase = Axis(fig[1, 1],
    limits = (nothing, (0.0, 1.5 * maximum(ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)
ATP_prod_line = lines!(ax_ATP_prod, timepoints, ATPprod, color = full_color)
no_reg_ATP_prod_line = lines!(ax_ATP_prod,
    no_reg_timepoints,
    no_reg_ATPprod,
    linestyle = no_reg_linestyle,
    color = no_reg_color)
ATPase_line = lines!(ax_ATPase,
    timepoints,
    ATPase,
    linestyle = :dot,
    color = line_ATPase_color,
    linewidth = line_ATPase_width)

axislegend(ax_ATP_prod,
    [ATP_prod_line, no_reg_ATP_prod_line, ATPase_line],
    ["Full Model", "Model w/o allost.", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (0, -4, 0, -4),
    patchsize = (10, 5))

# Plot [ATP]
ax_ATP_conc = Axis(fig[1, 2],
    limits = (nothing, (2e-5, 15e-3)),
    xlabel = "Time, min",
    ylabel = "[ATP], mM",
    title = "Dynamically maintaining\nATP concentration",
    yscale = log10,
    ytickformat = ys -> ["$(round(y*1000, sigdigits = 3))" for y in ys])
ax_ATPase = Axis(fig[1, 2],
    limits = (nothing, (0.0, 1.5 * maximum(ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)
ATP_line = lines!(ax_ATP_conc, timepoints, ATP, color = full_color)
no_reg_ATP_line = lines!(ax_ATP_conc,
    no_reg_timepoints,
    no_reg_ATP,
    linestyle = no_reg_linestyle,
    color = no_reg_color)
ATPase_line = lines!(ax_ATPase,
    timepoints,
    ATPase,
    linestyle = :dot,
    color = line_ATPase_color,
    linewidth = line_ATPase_width)

lines!(ax_ATP_conc,
    [tspan[1], tspan[2]],
    repeat([glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP],
        2),
    color = :grey,
    linestyle = :dash)
text!(ax_ATP_conc,
    0,
    1.05 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = :grey)
axislegend(ax_ATP_conc,
    [ATP_line, no_reg_ATP_line, ATPase_line],
    ["Full Model", "Model w/o allost.", "ATPase rate"],
    # position = :rt,
    position = (1, 0.88),
    rowgap = 1,
    framevisible = false,
    padding = (0, -4, 0, -4),
    patchsize = (10, 5))

# Plot energy of ATPase reaction
ax_ATP_energy = Axis(fig[1, 3],
    limits = (nothing, (0, 35)),
    xlabel = "Time, min",
    ylabel = rich("Energy of ATPase reaction, k", subscript("B"), "T"),
    title = "Dynamically maintaining\nATP energy",
    ytickformat = ys -> ["$(Int(round(y)))" for y in ys])
ax_ATPase = Axis(fig[1, 3],
    limits = (nothing, (0.0, 1.5 * maximum(ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)
ATPase_energy_line = lines!(ax_ATP_energy, timepoints, -log.(ATPenergy), color = full_color)
no_reg_ATPase_energy_line = lines!(ax_ATP_energy,
    no_reg_timepoints,
    -log.(no_reg_ATPenergy),
    linestyle = no_reg_linestyle,
    color = no_reg_color)
ATPase_line = lines!(ax_ATPase,
    timepoints,
    ATPase,
    linestyle = :dot,
    color = line_ATPase_color,
    linewidth = line_ATPase_width)
axislegend(ax_ATP_energy,
    [ATPase_energy_line, no_reg_ATPase_energy_line, ATPase_line],
    ["Full Model", "Model w/o allost.", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (0, -4, 0, -4),
    patchsize = (10, 5))

# ATPase tracking with ATP production
ax_ATPase_range = Axis(fig[2, 1],
    limits = ((0.003, 1.0), (0, 1.2)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "Ratio of ATP prod.\n to ATPase rates",
    title = "Matching steady state\nATP supply and demand",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9],
        ["0.3%", "1%", "3%", "10%", "30%", "90%"]),
    ytickformat = ys -> ["$(round(y, sigdigits=1))" for y in ys])
ATP_line = lines!(ax_ATPase_range,
    Complete_Model_Simulation_Data.ATPase_Vmax,
    Complete_Model_Simulation_Data.ATP_prod_rate,
    color = full_color)
ATP_line_no_reg = lines!(ax_ATPase_range,
    No_Reg_Model_Simulation_Data.ATPase_Vmax,
    No_Reg_Model_Simulation_Data.ATP_prod_rate,
    color = no_reg_color,
    linestyle = :dash)
axislegend(ax_ATPase_range,
    [ATP_line, ATP_line_no_reg],
    ["Full model", "Model w/o allost."],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (0, -4, 0, -4),
    patchsize = (10, 5))

# Plot [ATP] maintenance
ax_ATPase_range = Axis(fig[2, 2],
    limits = ((0.003, 1.0), (-0.15e-3, 11e-3)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "[ATP],mM",
    title = "Maintaining steady state\nATP concentration",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9],
        ["0.3%", "1%", "3%", "10%", "30%", "90%"]),
    ytickformat = ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys])
ATP_line = lines!(ax_ATPase_range,
    Complete_Model_Simulation_Data.ATPase_Vmax,
    Complete_Model_Simulation_Data.ATP_conc,
    color = full_color)
ATP_line_no_reg = lines!(ax_ATPase_range,
    No_Reg_Model_Simulation_Data.ATPase_Vmax,
    No_Reg_Model_Simulation_Data.ATP_conc,
    color = no_reg_color,
    linestyle = no_reg_linestyle)
lines!(ax_ATPase_range,
    [0.003, 1.0],
    repeat([glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP],
        2),
    color = :grey,
    linestyle = :dash)
text!(ax_ATPase_range,
    0.0033,
    1.02 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = :grey)
axislegend(ax_ATPase_range,
    [ATP_line, ATP_line_no_reg],
    ["Full model", "Model w/o allost."],
    position = :rt,
    # position = (1.075, 1.05),
    rowgap = 1,
    framevisible = false,
    padding = (0, -4, 0, -4),
    patchsize = (10, 5))

# Plot ATPase energy
ax_ATPase_range = Axis(fig[2, 3],
    limits = ((0.003, 1.0), (0, 35)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "ATP hydrolysis energy, kᵦT",
    title = "Maintaining steady state\nATP energy",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9],
        ["0.3%", "1%", "3%", "10%", "30%", "90%"]),
    ytickformat = ys -> ["$(Int(round(y, sigdigits=2)))" for y in ys]
    # yticklabelcolor = full_color,
    # ylabelcolor = full_color,
)
ATP_line = lines!(ax_ATPase_range,
    Complete_Model_Simulation_Data.ATPase_Vmax,
    Complete_Model_Simulation_Data.ATP_energy,
    color = full_color)
ATP_line_no_reg = lines!(ax_ATPase_range,
    No_Reg_Model_Simulation_Data.ATPase_Vmax,
    No_Reg_Model_Simulation_Data.ATP_energy,
    color = no_reg_color,
    linestyle = :dash)
axislegend(ax_ATPase_range,
    [ATP_line, ATP_line_no_reg],
    ["Full model", "Model w/o allost."],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (0, -4, 0, -4),
    patchsize = (10, 5))

colgap!(fig.layout, 7.5)
rowgap!(fig.layout, 7.5)

label_a = fig[1, 1, TopLeft()] = Label(fig,
    "A",
    fontsize = 12,
    halign = :right,
    padding = (0, 5, 5, 0))
label_b = fig[1, 2, TopLeft()] = Label(fig,
    "B",
    fontsize = 12,
    halign = :right,
    padding = (0, 5, 5, 0))
label_c = fig[1, 3, TopLeft()] = Label(fig,
    "C",
    fontsize = 12,
    halign = :right,
    padding = (0, 5, 5, 0))
label_d = fig[2, 1, TopLeft()] = Label(fig,
    "D",
    fontsize = 12,
    halign = :right,
    padding = (0, 5, 5, 0))
label_e = fig[2, 2, TopLeft()] = Label(fig,
    "E",
    fontsize = 12,
    halign = :right,
    padding = (0, 5, 5, 0))
label_f = fig[2, 3, TopLeft()] = Label(fig,
    "F",
    fontsize = 12,
    halign = :right,
    padding = (0, 5, 5, 0))

fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_Fig3_allostery_required_for_ATP_maintenence.png", fig, px_per_unit = 4)
