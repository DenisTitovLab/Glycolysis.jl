using Glycolysis
using DifferentialEquations, ProgressMeter
using CairoMakie, Dates, Printf, Statistics, StatsBase
using DataFrames, CSV

##
# Precalculate output of our model
function find_ATP_at_ATPase_range(glycolysis_params,
        glycolysis_init_conc;
        n_Vmax_ATPase_values = 1000)
    tspan = (0.0, 1e8)
    pathway_Vmax = 2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc
    ATPases = 10 .^ range(log10(0.003), log10(0.9), n_Vmax_ATPase_values) .* pathway_Vmax
    prob = ODEProblem(glycolysis_ODEs, glycolysis_init_conc, tspan, glycolysis_params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(ensemble_prob,
        Rodas5P(),
        EnsembleThreads(),
        trajectories = 1_000,
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
        save_start = false)
    ATP_prod_rate = [Glycolysis.conc_to_rates(sol.u[end], sol.prob.p).ATPprod /
                     sol.prob.p.ATPase_Vmax
                     for sol in sim if sol.retcode == ReturnCode.Success]
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
Our_Model_Simulation_Data = find_ATP_at_ATPase_range(glycolysis_params_copy,
    glycolysis_init_conc)
supported_ATPase_range = Our_Model_Simulation_Data.ATPase_Vmax[Our_Model_Simulation_Data.ATP_conc .> 0.5 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP)]
if !isempty(supported_ATPase_range)
    supported_ATPase_fold_range = round(supported_ATPase_range[end] /
                                        supported_ATPase_range[1],
        sigdigits = 2)
    println("Our model ATPase range where ATP>0.5*(ATP+ADP+AMP) is $(supported_ATPase_fold_range) fold")
else
    println("Our model ATPase range where ATP>0.5*(ATP+ADP+AMP) is 0 fold")
end

# Precalculate output of Mulquiney 1999 model
include("for_FigS2_Mulquiney1999Model.jl")
# mulquiney_glycolysis_init_conc = deepcopy(glycolysis_init_conc)
function mulquiney_find_ATP_at_ATPase_range(mulquiney_glycolysis_params,
        mulquiney_glycolysis_init_conc;
        n_Vmax_ATPase_values = 1000)
    tspan = (0.0, 1e8)
    pathway_Vmax = 2 * mulquiney_glycolysis_params.HK1_kcat_f *
                   mulquiney_glycolysis_params.HK1_e0
    ATPases = 10 .^ range(log10(0.003), log10(0.9), n_Vmax_ATPase_values) .* pathway_Vmax
    prob = ODEProblem(mulquiney_glycolysis_ODEs,
        mulquiney_glycolysis_init_conc,
        tspan,
        mulquiney_glycolysis_params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(ensemble_prob,
        Rodas5P(),
        EnsembleThreads(),
        trajectories = 1_000,
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
        save_start = false)
    ATP_prod_rate = [mulquiney_conc_to_rates(sol.u[end], sol.prob.p).ATPprod /
                     sol.prob.p.ATPase_Vmax
                     for sol in sim if sol.retcode == ReturnCode.Success]
    ATP_energy = -log.([mulquiney_conc_to_disequilibrium_ratios(sol.u[end],
        sol.prob.p).Q_Keq_ATPase for sol in sim if sol.retcode == ReturnCode.Success])
    ATP_conc = [sol.u[end].ATP for sol in sim if sol.retcode == ReturnCode.Success]
    return (ATP_conc = ATP_conc,
        ATP_prod_rate = ATP_prod_rate,
        ATP_energy = ATP_energy,
        ATPase_Vmax = ATPases / (2 * mulquiney_glycolysis_params.HK1_kcat_f *
                       mulquiney_glycolysis_params.HK1_e0))
end
mulquiney_glycolysis_params.ATPase_Km_ATP = 1e-9
mulquiney_glycolysis_params_copy = deepcopy(mulquiney_glycolysis_params)
mulquiney_Model_Simulation_Data = mulquiney_find_ATP_at_ATPase_range(mulquiney_glycolysis_params_copy,
    mulquiney_glycolysis_init_conc)
supported_ATPase_range = mulquiney_Model_Simulation_Data.ATPase_Vmax[mulquiney_Model_Simulation_Data.ATP_conc .> 0.5 * (mulquiney_glycolysis_init_conc.ATP + mulquiney_glycolysis_init_conc.ADP + mulquiney_glycolysis_init_conc.AMP)]
if !isempty(supported_ATPase_range)
    supported_ATPase_fold_range = round(supported_ATPase_range[end] /
                                        supported_ATPase_range[1],
        sigdigits = 2)
    println("Mulquiney model ATPase range where ATP>0.5*(ATP+ADP+AMP) is $(supported_ATPase_fold_range) fold")
else
    println("Mulquiney model ATPase range where ATP>0.5*(ATP+ADP+AMP) is 0 fold")
end
no_reg_mulquiney_glycolysis_params = deepcopy(mulquiney_glycolysis_params)
no_reg_mulquiney_glycolysis_params.HK1_Ki_G6P = Inf
no_reg_mulquiney_glycolysis_params.PFKP_L = 0.0
no_reg_mulquiney_glycolysis_params.PKM2_L = 0.0
no_reg_mulquiney_glycolysis_params.ATPase_Km_ATP = 1e-9
no_reg_mulquiney_Model_Simulation_Data = mulquiney_find_ATP_at_ATPase_range(no_reg_mulquiney_glycolysis_params,
    mulquiney_glycolysis_init_conc)

# Precalculate output of vanHeerden 2014 model
include("for_FigS2_vanHeerden2014Model.jl")
# van_heerden_glycolysis_init_conc = deepcopy(glycolysis_init_conc)
function van_heerden_find_ATP_at_ATPase_range(van_heerden_glycolysis_params,
        van_heerden_glycolysis_init_conc;
        n_Vmax_ATPase_values = 1000)
    tspan = (0.0, 1e8)
    pathway_Vmax = 2 * van_heerden_glycolysis_params.PFKP_Vmax
    ATPases = 10 .^ range(log10(0.003), log10(0.9), n_Vmax_ATPase_values) .* pathway_Vmax
    prob = ODEProblem(van_heerden_glycolysis_ODEs,
        van_heerden_glycolysis_init_conc,
        tspan,
        van_heerden_glycolysis_params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(ensemble_prob,
        Rodas4P2(),
        EnsembleThreads(),
        trajectories = 1_000,
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
        save_start = false)
    ATP_prod_rate = [van_heerden_conc_to_rates(sol.u[end], sol.prob.p).ATPprod /
                     sol.prob.p.ATPase_Vmax
                     for sol in sim if sol.retcode == ReturnCode.Success]
    ATP_energy = -log.([van_heerden_conc_to_disequilibrium_ratios(sol.u[end],
        sol.prob.p).Q_Keq_ATPase for sol in sim if sol.retcode == ReturnCode.Success])
    ATP_conc = [sol.u[end].ATP for sol in sim if sol.retcode == ReturnCode.Success]
    return (ATP_conc = ATP_conc,
        ATP_prod_rate = ATP_prod_rate,
        ATP_energy = ATP_energy,
        ATPase_Vmax = ATPases / pathway_Vmax)
end
van_heerden_glycolysis_params.ATPase_Km_ATP = 1e-9
van_heerden_glycolysis_params_copy = deepcopy(van_heerden_glycolysis_params)
van_heerden_Model_Simulation_Data = van_heerden_find_ATP_at_ATPase_range(van_heerden_glycolysis_params_copy,
    van_heerden_glycolysis_init_conc)
supported_ATPase_range = van_heerden_Model_Simulation_Data.ATPase_Vmax[van_heerden_Model_Simulation_Data.ATP_conc .> 0.5 * (van_heerden_glycolysis_init_conc.ATP + van_heerden_glycolysis_init_conc.ADP + van_heerden_glycolysis_init_conc.AMP)]
if !isempty(supported_ATPase_range)
    supported_ATPase_fold_range = round(supported_ATPase_range[end] /
                                        supported_ATPase_range[1],
        sigdigits = 2)
    println("van Heerden model ATPase range where ATP>0.5*(ATP+ADP+AMP) is $(supported_ATPase_fold_range) fold")
else
    println("van Heerden model ATPase range where ATP>0.5*(ATP+ADP+AMP) is 0 fold")
end
no_reg_van_heerden_glycolysis_params = deepcopy(van_heerden_glycolysis_params)
no_reg_van_heerden_glycolysis_params.HK1_KiG6P = Inf
no_reg_van_heerden_glycolysis_params.PFKP_L0 = 0.0
no_reg_van_heerden_glycolysis_params.PKM2_L0 = 0.0
no_reg_van_heerden_glycolysis_params.ATPase_Km_ATP = 1e-9
no_reg_van_heerden_Model_Simulation_Data = van_heerden_find_ATP_at_ATPase_range(no_reg_van_heerden_glycolysis_params,
    van_heerden_glycolysis_init_conc)

# Precalculate output of Shestov 2014 model
include("for_FigS2_Shestov2014Model.jl")
# shestov_glycolysis_init_conc = deepcopy(glycolysis_init_conc)
shestov_ATPase_Vmax_shift = 10000
function shestov_find_ATP_at_ATPase_range(shestov_glycolysis_params,
        shestov_glycolysis_init_conc;
        n_Vmax_ATPase_values = 1000,
        shestov_ATPase_Vmax_shift = shestov_ATPase_Vmax_shift)
    tspan = (0.0, 1e8)
    # pathway_Vmax = 2 * shestov_glycolysis_params.HK1_Vmax * shestov_glycolysis_params.HK1_Conc
    pathway_Vmax = shestov_glycolysis_params.MCT_Vmax_ltr

    ATPases = 10 .^ range(log10(0.003), log10(0.9), n_Vmax_ATPase_values) .* pathway_Vmax ./
              shestov_ATPase_Vmax_shift
    prob = ODEProblem(shestov_glycolysis_ODEs,
        shestov_glycolysis_init_conc,
        tspan,
        shestov_glycolysis_params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax_f = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(ensemble_prob,
        Rodas5P(),
        EnsembleThreads(),
        trajectories = 1_000,
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
        save_start = false)
    ATP_prod_rate = [shestov_conc_to_rates(sol.u[end], sol.prob.p).ATPprod /
                     sol.prob.p.ATPase_Vmax_f
                     for sol in sim if sol.retcode == ReturnCode.Success]
    ATP_energy = -log.([shestov_conc_to_disequilibrium_ratios(sol.u[end],
        sol.prob.p).Q_Keq_ATPase for sol in sim if sol.retcode == ReturnCode.Success])
    ATP_conc = [sol.u[end].ATP for sol in sim if sol.retcode == ReturnCode.Success]
    return (ATP_conc = ATP_conc,
        ATP_prod_rate = ATP_prod_rate,
        ATP_energy = ATP_energy,
        ATPase_Vmax = ATPases / pathway_Vmax)
end
shestov_glycolysis_params_copy = deepcopy(shestov_glycolysis_params)
shestov_glycolysis_params_copy.ATPase_Km_f = 1e-9
shestov_Model_Simulation_Data = shestov_find_ATP_at_ATPase_range(shestov_glycolysis_params_copy,
    shestov_glycolysis_init_conc)
supported_ATPase_range = shestov_Model_Simulation_Data.ATPase_Vmax[shestov_Model_Simulation_Data.ATP_conc .> 0.5 * (shestov_glycolysis_init_conc.ATP + shestov_glycolysis_init_conc.ADP + shestov_glycolysis_init_conc.AMP)]
if !isempty(supported_ATPase_range)
    supported_ATPase_fold_range = round(supported_ATPase_range[end] /
                                        supported_ATPase_range[1],
        sigdigits = 2)
    println("Shestov model ATPase range where ATP>0.5*(ATP+ADP+AMP) is $(supported_ATPase_fold_range) fold")
else
    println("Shestov model ATPase range where ATP>0.5*(ATP+ADP+AMP) is 0 fold")
end
no_reg_shestov_glycolysis_params = deepcopy(shestov_glycolysis_params)
no_reg_shestov_glycolysis_params.HK1_L = 0.0
no_reg_shestov_glycolysis_params.PFKP_Ka = 0.0
no_reg_shestov_glycolysis_params.PKM2_Ka = 0.0
no_reg_shestov_glycolysis_params.ATPase_Km_f = 1e-9
no_reg_shestov_Model_Simulation_Data = shestov_find_ATP_at_ATPase_range(no_reg_shestov_glycolysis_params,
    shestov_glycolysis_init_conc)

##
# Plot the results
size_inches = (6.5, 5)
size_pt = 72 .* size_inches

CairoMakie.set_theme!(CairoMakie.Theme(fontsize = 6,
    Axis = (xticksize = 1,
        yticksize = 1,
        # xticklabelsize = 6,
        # yticklabelsize = 6,
        yticklabelpad = 1,
        ylabelpadding = 3)))
fig = Figure(size = size_pt)

# our_model_color = :black
# mulquiney_color = :red

our_model_color = :Grey
our_model_linestyle = :solid

van_heerden_color = :Black
van_heerden_linestyle = :dash
no_reg_van_heerden_color = :Red
no_reg_van_heerden_linestyle = :dot

shestov_color = :Black
shestov_linestyle = :dash
no_reg_shestov_color = :Red
no_reg_shestov_linestyle = :dot

mulquiney_color = :Black
mulquiney_linestyle = :dash
no_reg_mulquiney_color = :Red
no_reg_mulquiney_linestyle = :dot

# ATPase tracking with ATP production
ax_ATPase_range = Axis(fig[1, 1],
    limits = ((0.003, 0.9), (0, 1.9)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "Ratio of ATP prod.\n to ATPase rates",
    title = "Matching steady state\nATP supply and demand",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9]),
    xtickformat = xs -> ["$(round(x*100, sigdigits=1))" for x in xs],
    ytickformat = ys -> ["$(round(y, sigdigits=1))" for y in ys])

Our_model_line = lines!(ax_ATPase_range,
    Our_Model_Simulation_Data.ATPase_Vmax,
    Our_Model_Simulation_Data.ATP_prod_rate,
    color = our_model_color,
    linestyle = our_model_linestyle)
mulquiney_model_line = lines!(ax_ATPase_range,
    mulquiney_Model_Simulation_Data.ATPase_Vmax,
    mulquiney_Model_Simulation_Data.ATP_prod_rate,
    color = mulquiney_color,
    linestyle = mulquiney_linestyle)
no_reg_mulquiney_model_line = lines!(ax_ATPase_range,
    no_reg_mulquiney_Model_Simulation_Data.ATPase_Vmax,
    no_reg_mulquiney_Model_Simulation_Data.ATP_prod_rate,
    color = no_reg_mulquiney_color,
    linestyle = no_reg_mulquiney_linestyle)

axislegend(ax_ATPase_range,
    [mulquiney_model_line, no_reg_mulquiney_model_line],
    ["Full Model", "Model w/o allost."],
    "Mulquiney Model",
    position = :ct,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (-5, -5, 0, -4),
    patchsize = (10, 5))

axislegend(ax_ATPase_range,
    [Our_model_line],
    ["Our Model"],
    position = :cc,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (0, -5, 0, -24),
    patchsize = (10, 5))

# Plot [ATP] maintenance
Our_ATP_pool_size = glycolysis_init_conc.ATP + glycolysis_init_conc.ADP +
                    glycolysis_init_conc.AMP
mulquiney_ATP_pool_size = mulquiney_glycolysis_init_conc.ATP +
                          mulquiney_glycolysis_init_conc.ADP +
                          mulquiney_glycolysis_init_conc.AMP
van_heerden_ATP_pool_size = van_heerden_glycolysis_init_conc.ATP +
                            van_heerden_glycolysis_init_conc.ADP +
                            van_heerden_glycolysis_init_conc.AMP
shestov_ATP_pool_size = shestov_glycolysis_init_conc.ATP +
                        shestov_glycolysis_init_conc.ADP + shestov_glycolysis_init_conc.AMP
ax_ATPase_range = Axis(fig[1, 2],
    limits = ((0.003, 0.9), (-0.01, 2.0)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "[ATP]/([ATP]+[ADP]+[AMP])",
    title = "Maintaining steady state\nATP concentration",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9]),
    xtickformat = xs -> ["$(round(x*100, sigdigits=1))" for x in xs]
    # ytickformat = ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys],
)

Our_model_line = lines!(ax_ATPase_range,
    Our_Model_Simulation_Data.ATPase_Vmax,
    Our_Model_Simulation_Data.ATP_conc ./ Our_ATP_pool_size,
    color = our_model_color,
    linestyle = our_model_linestyle)
mulquiney_model_line = lines!(ax_ATPase_range,
    mulquiney_Model_Simulation_Data.ATPase_Vmax,
    mulquiney_Model_Simulation_Data.ATP_conc ./ mulquiney_ATP_pool_size,
    color = mulquiney_color,
    linestyle = mulquiney_linestyle)
no_reg_mulquiney_model_line = lines!(ax_ATPase_range,
    no_reg_mulquiney_Model_Simulation_Data.ATPase_Vmax,
    no_reg_mulquiney_Model_Simulation_Data.ATP_conc ./ mulquiney_ATP_pool_size,
    color = no_reg_mulquiney_color,
    linestyle = no_reg_mulquiney_linestyle)

lines!(ax_ATPase_range,
    [0.003, 1.0],
    repeat([1.0], 2),
    color = (:grey, 0.5),
    linestyle = :dot)
text!(ax_ATPase_range,
    0.0033,
    1.02,
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = (:grey, 0.5))
axislegend(ax_ATPase_range,
    [mulquiney_model_line, no_reg_mulquiney_model_line],
    ["Full Model", "Model w/o allost."],
    "Mulquiney Model",
    position = :ct,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (-5, -5, 0, -4),
    patchsize = (10, 5))

axislegend(ax_ATPase_range,
    [Our_model_line],
    ["Our Model"],
    position = :cc,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (0, -5, 0, -24),
    patchsize = (10, 5))

# Plot ATPase energy
ax_ATPase_range = Axis(fig[1, 3],
    limits = ((0.003, 0.9), (0, 45)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "ATP hydrolysis energy, kᵦT",
    title = "Maintaining steady state\nATP energy",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9]),
    xtickformat = xs -> ["$(round(x*100, sigdigits=1))" for x in xs],
    ytickformat = ys -> ["$(Int(round(y, sigdigits=2)))" for y in ys]
    # yticklabelcolor = our_model_color,
    # ylabelcolor = our_model_color,
)

Our_model_line = lines!(ax_ATPase_range,
    Our_Model_Simulation_Data.ATPase_Vmax,
    Our_Model_Simulation_Data.ATP_energy,
    color = our_model_color,
    linestyle = our_model_linestyle)
mulquiney_model_line = lines!(ax_ATPase_range,
    mulquiney_Model_Simulation_Data.ATPase_Vmax,
    mulquiney_Model_Simulation_Data.ATP_energy,
    color = mulquiney_color,
    linestyle = mulquiney_linestyle)
no_reg_mulquiney_model_line = lines!(ax_ATPase_range,
    no_reg_mulquiney_Model_Simulation_Data.ATPase_Vmax,
    no_reg_mulquiney_Model_Simulation_Data.ATP_energy,
    color = no_reg_mulquiney_color,
    linestyle = no_reg_mulquiney_linestyle)

axislegend(ax_ATPase_range,
    [mulquiney_model_line, no_reg_mulquiney_model_line],
    ["Full Model", "Model w/o allost."],
    "Mulquiney Model",
    position = :ct,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (-5, -5, 0, -4),
    patchsize = (10, 5))

axislegend(ax_ATPase_range,
    [Our_model_line],
    ["Our Model"],
    position = :cc,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (0, -5, 0, -24),
    patchsize = (10, 5))

axislegend(ax_ATPase_range,
    [Our_model_line],
    ["Our Model"],
    position = :cc,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (0, -5, 0, -24),
    patchsize = (10, 5))

# ATPase tracking with ATP production
ax_ATPase_range = Axis(fig[2, 1],
    limits = ((0.003, 0.9) ./ shestov_ATPase_Vmax_shift, (0, 1.9)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "Ratio of ATP prod.\n to ATPase rates",
    title = "Matching steady state\nATP supply and demand",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9] ./ shestov_ATPase_Vmax_shift),
    xtickformat = xs -> ["$(round(x*100, sigdigits=1))" for x in xs],
    ytickformat = ys -> ["$(round(y, sigdigits=1))" for y in ys])

# Our_model_line = lines!(
#     ax_ATPase_range,
#     Our_Model_Simulation_Data.ATPase_Vmax,
#     Our_Model_Simulation_Data.ATP_prod_rate,
#     color = our_model_color,
#     linestyle = our_model_linestyle
# )
shestov_model_line = lines!(ax_ATPase_range,
    shestov_Model_Simulation_Data.ATPase_Vmax,
    shestov_Model_Simulation_Data.ATP_prod_rate,
    color = shestov_color,
    linestyle = shestov_linestyle)
no_reg_shestov_model_line = lines!(ax_ATPase_range,
    no_reg_shestov_Model_Simulation_Data.ATPase_Vmax,
    no_reg_shestov_Model_Simulation_Data.ATP_prod_rate,
    color = no_reg_shestov_color,
    linestyle = no_reg_shestov_linestyle)

axislegend(ax_ATPase_range,
    [shestov_model_line, no_reg_shestov_model_line],
    ["Full Model", "Model w/o allost."],
    "Shestov Model",
    position = :ct,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (-5, -5, 0, -4),
    patchsize = (10, 5))

# Plot [ATP] maintenance
Our_ATP_pool_size = glycolysis_init_conc.ATP + glycolysis_init_conc.ADP +
                    glycolysis_init_conc.AMP
mulquiney_ATP_pool_size = mulquiney_glycolysis_init_conc.ATP +
                          mulquiney_glycolysis_init_conc.ADP +
                          mulquiney_glycolysis_init_conc.AMP
van_heerden_ATP_pool_size = van_heerden_glycolysis_init_conc.ATP +
                            van_heerden_glycolysis_init_conc.ADP +
                            van_heerden_glycolysis_init_conc.AMP
shestov_ATP_pool_size = shestov_glycolysis_init_conc.ATP +
                        shestov_glycolysis_init_conc.ADP + shestov_glycolysis_init_conc.AMP
ax_ATPase_range = Axis(fig[2, 2],
    limits = ((0.003, 0.9) ./ shestov_ATPase_Vmax_shift, (-0.01, 2.0)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "[ATP]/([ATP]+[ADP]+[AMP])",
    title = "Maintaining steady state\nATP concentration",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9] ./ shestov_ATPase_Vmax_shift),
    xtickformat = xs -> ["$(round(x*100, sigdigits=1))" for x in xs]
    # ytickformat = ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys],
)

# Our_model_line = lines!(
#     ax_ATPase_range,
#     Our_Model_Simulation_Data.ATPase_Vmax,
#     Our_Model_Simulation_Data.ATP_conc ./ Our_ATP_pool_size,
#     color = our_model_color,
#     linestyle = our_model_linestyle
# )
shestov_model_line = lines!(ax_ATPase_range,
    shestov_Model_Simulation_Data.ATPase_Vmax,
    shestov_Model_Simulation_Data.ATP_conc ./ shestov_ATP_pool_size,
    color = shestov_color,
    linestyle = shestov_linestyle)
no_reg_shestov_model_line = lines!(ax_ATPase_range,
    no_reg_shestov_Model_Simulation_Data.ATPase_Vmax,
    no_reg_shestov_Model_Simulation_Data.ATP_conc ./ shestov_ATP_pool_size,
    color = no_reg_shestov_color,
    linestyle = no_reg_shestov_linestyle)

lines!(ax_ATPase_range,
    [0.003, 1.0] ./ shestov_ATPase_Vmax_shift,
    repeat([1.0], 2),
    color = (:grey, 0.5),
    linestyle = :dot)
text!(ax_ATPase_range,
    0.0033 / shestov_ATPase_Vmax_shift,
    1.02,
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = (:grey, 0.5))
axislegend(ax_ATPase_range,
    [shestov_model_line, no_reg_shestov_model_line],
    ["Full Model", "Model w/o allost."],
    "Shestov Model",
    position = :ct,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (-5, -5, 0, -4),
    patchsize = (10, 5))

# Plot ATPase energy
ax_ATPase_range = Axis(fig[2, 3],
    limits = ((0.003, 0.9) ./ shestov_ATPase_Vmax_shift, (0, 45)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "ATP hydrolysis energy, kᵦT",
    title = "Maintaining steady state\nATP energy",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9] ./ shestov_ATPase_Vmax_shift),
    xtickformat = xs -> ["$(round(x*100, sigdigits=1))" for x in xs],
    ytickformat = ys -> ["$(Int(round(y, sigdigits=2)))" for y in ys]
    # yticklabelcolor = our_model_color,
    # ylabelcolor = our_model_color,
)

# Our_model_line = lines!(
#     ax_ATPase_range,
#     Our_Model_Simulation_Data.ATPase_Vmax,
#     Our_Model_Simulation_Data.ATP_energy,
#     color = our_model_color,
#     linestyle = our_model_linestyle
# )
shestov_model_line = lines!(ax_ATPase_range,
    shestov_Model_Simulation_Data.ATPase_Vmax,
    shestov_Model_Simulation_Data.ATP_energy,
    color = shestov_color,
    linestyle = shestov_linestyle)
no_reg_shestov_model_line = lines!(ax_ATPase_range,
    no_reg_shestov_Model_Simulation_Data.ATPase_Vmax,
    no_reg_shestov_Model_Simulation_Data.ATP_energy,
    color = no_reg_shestov_color,
    linestyle = no_reg_shestov_linestyle)

axislegend(ax_ATPase_range,
    [shestov_model_line, no_reg_shestov_model_line],
    ["Full Model", "Model w/o allost."],
    "Shestov Model",
    position = :ct,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (-5, -5, 0, -4),
    patchsize = (10, 5))

# ATPase tracking with ATP production
ax_ATPase_range = Axis(fig[3, 1],
    limits = ((0.003, 0.9), (0, 1.9)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "Ratio of ATP prod.\n to ATPase rates",
    title = "Matching steady state\nATP supply and demand",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9]),
    xtickformat = xs -> ["$(round(x*100, sigdigits=1))" for x in xs],
    ytickformat = ys -> ["$(round(y, sigdigits=1))" for y in ys])

Our_model_line = lines!(ax_ATPase_range,
    Our_Model_Simulation_Data.ATPase_Vmax,
    Our_Model_Simulation_Data.ATP_prod_rate,
    color = our_model_color,
    linestyle = our_model_linestyle)
van_heerden_model_line = lines!(ax_ATPase_range,
    van_heerden_Model_Simulation_Data.ATPase_Vmax,
    van_heerden_Model_Simulation_Data.ATP_prod_rate,
    color = van_heerden_color,
    linestyle = van_heerden_linestyle)
no_reg_van_heerden_model_line = lines!(ax_ATPase_range,
    no_reg_van_heerden_Model_Simulation_Data.ATPase_Vmax,
    no_reg_van_heerden_Model_Simulation_Data.ATP_prod_rate,
    color = no_reg_van_heerden_color,
    linestyle = no_reg_van_heerden_linestyle)

axislegend(ax_ATPase_range,
    [van_heerden_model_line, no_reg_van_heerden_model_line],
    ["Full Model", "Model w/o allost."],
    "van Heerden Model",
    position = :ct,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (-5, -5, 0, -4),
    patchsize = (10, 5))

axislegend(ax_ATPase_range,
    [Our_model_line],
    ["Our Model"],
    position = :cc,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (0, -5, 0, -24),
    patchsize = (10, 5))

# Plot [ATP] maintenance
Our_ATP_pool_size = glycolysis_init_conc.ATP + glycolysis_init_conc.ADP +
                    glycolysis_init_conc.AMP
mulquiney_ATP_pool_size = mulquiney_glycolysis_init_conc.ATP +
                          mulquiney_glycolysis_init_conc.ADP +
                          mulquiney_glycolysis_init_conc.AMP
van_heerden_ATP_pool_size = van_heerden_glycolysis_init_conc.ATP +
                            van_heerden_glycolysis_init_conc.ADP +
                            van_heerden_glycolysis_init_conc.AMP
shestov_ATP_pool_size = shestov_glycolysis_init_conc.ATP +
                        shestov_glycolysis_init_conc.ADP + shestov_glycolysis_init_conc.AMP
ax_ATPase_range = Axis(fig[3, 2],
    limits = ((0.003, 0.9), (-0.01, 2.0)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "[ATP]/([ATP]+[ADP]+[AMP])",
    title = "Maintaining steady state\nATP concentration",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9]),
    xtickformat = xs -> ["$(round(x*100, sigdigits=1))" for x in xs]
    # ytickformat = ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys],
)

Our_model_line = lines!(ax_ATPase_range,
    Our_Model_Simulation_Data.ATPase_Vmax,
    Our_Model_Simulation_Data.ATP_conc ./ Our_ATP_pool_size,
    color = our_model_color,
    linestyle = our_model_linestyle)
van_heerden_model_line = lines!(ax_ATPase_range,
    van_heerden_Model_Simulation_Data.ATPase_Vmax,
    van_heerden_Model_Simulation_Data.ATP_conc ./ van_heerden_ATP_pool_size,
    color = van_heerden_color,
    linestyle = van_heerden_linestyle)
no_reg_van_heerden_model_line = lines!(ax_ATPase_range,
    no_reg_van_heerden_Model_Simulation_Data.ATPase_Vmax,
    no_reg_van_heerden_Model_Simulation_Data.ATP_conc ./ van_heerden_ATP_pool_size,
    color = no_reg_van_heerden_color,
    linestyle = no_reg_van_heerden_linestyle)

lines!(ax_ATPase_range,
    [0.0003, 1.0] .* 2,
    repeat([1.0], 2),
    color = (:grey, 0.5),
    linestyle = :dot)
text!(ax_ATPase_range,
    0.0033 * 2,
    1.02,
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = (:grey, 0.5))
axislegend(ax_ATPase_range,
    [van_heerden_model_line, no_reg_van_heerden_model_line],
    ["Full Model", "Model w/o allost."],
    "van Heerden Model",
    position = :ct,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (-5, -5, 0, -4),
    patchsize = (10, 5))

axislegend(ax_ATPase_range,
    [Our_model_line],
    ["Our Model"],
    position = :cc,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (0, -5, 0, -24),
    patchsize = (10, 5))

# Plot ATPase energy
ax_ATPase_range = Axis(fig[3, 3],
    limits = ((0.003, 0.9), (0, 45)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "ATP hydrolysis energy, kᵦT",
    title = "Maintaining steady state\nATP energy",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9]),
    xtickformat = xs -> ["$(round(x*100, sigdigits=1))" for x in xs],
    ytickformat = ys -> ["$(Int(round(y, sigdigits=2)))" for y in ys]
    # yticklabelcolor = our_model_color,
    # ylabelcolor = our_model_color,
)

Our_model_line = lines!(ax_ATPase_range,
    Our_Model_Simulation_Data.ATPase_Vmax,
    Our_Model_Simulation_Data.ATP_energy,
    color = our_model_color,
    linestyle = our_model_linestyle)
van_heerden_model_line = lines!(ax_ATPase_range,
    van_heerden_Model_Simulation_Data.ATPase_Vmax,
    van_heerden_Model_Simulation_Data.ATP_energy,
    color = van_heerden_color,
    linestyle = van_heerden_linestyle)
no_reg_van_heerden_model_line = lines!(ax_ATPase_range,
    no_reg_van_heerden_Model_Simulation_Data.ATPase_Vmax,
    no_reg_van_heerden_Model_Simulation_Data.ATP_energy,
    color = no_reg_van_heerden_color,
    linestyle = no_reg_van_heerden_linestyle)

axislegend(ax_ATPase_range,
    [van_heerden_model_line, no_reg_van_heerden_model_line],
    ["Full Model", "Model w/o allost."],
    "van Heerden Model",
    position = :ct,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (-5, -5, 0, -4),
    patchsize = (10, 5))

axislegend(ax_ATPase_range,
    [Our_model_line],
    ["Our Model"],
    position = :cc,
    titlegap = 1,
    orientation = :horizontal,
    rowgap = 1,
    colgap = 7,
    framevisible = false,
    padding = (0, -5, 0, -24),
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
label_g = fig[3, 1, TopLeft()] = Label(fig,
    "G",
    fontsize = 12,
    halign = :right,
    padding = (0, 5, 5, 0))
label_h = fig[3, 2, TopLeft()] = Label(fig,
    "H",
    fontsize = 12,
    halign = :right,
    padding = (0, 5, 5, 0))
label_i = fig[3, 3, TopLeft()] = Label(fig,
    "I",
    fontsize = 12,
    halign = :right,
    padding = (0, 5, 5, 0))
fig

# uncomment the line below to save the plot
save("Results/$(Dates.format(now(),"mmddyy"))_Fig3-fig_suppl1_Mulquiney_Shestov_van_Heerden_models.png", fig, px_per_unit = 4)
