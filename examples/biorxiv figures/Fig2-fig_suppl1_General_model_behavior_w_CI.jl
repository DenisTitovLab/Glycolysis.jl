# it take about 3 min to run this code on 8-core computer
# the code is parallelized so the time will scale with the number of CPU cores

using Glycolysis
using DifferentialEquations
using CairoMakie, DataFrames, Dates, Printf, CSV, Statistics


##
# Precalculate output of complete model
using Distributed
addprocs(8; exeflags = "--project")
@everywhere using Glycolysis
@everywhere using DifferentialEquations, ProgressMeter
using Distributions, LabelledArrays, Statistics, StatsBase, DataFrames, CSV
using DataFrames, CSV, Dates

@everywhere Initial_ATPase_Vmax_frac = 0.06
@everywhere High_ATPase_Vmax_frac = Initial_ATPase_Vmax_frac * 2
@everywhere Low_ATPase_Vmax_frac = Initial_ATPase_Vmax_frac / 2
@everywhere ATPase_change_time1 = 60 / 4
@everywhere ATPase_change_time2 = 120 / 4
@everywhere ATPase_change_time3 = 180 / 4
@everywhere tspan = (0.0, 240.0 / 4)

@everywhere function glyc_step_response(
    params,
    init_conc;
    default_params = glycolysis_params,
    Initial_ATPase_Vmax_frac = Initial_ATPase_Vmax_frac,
    High_ATPase_Vmax_frac = High_ATPase_Vmax_frac,
    Low_ATPase_Vmax_frac = Low_ATPase_Vmax_frac,
    ATPase_change_time1 = ATPase_change_time1,
    ATPase_change_time2 = ATPase_change_time2,
    ATPase_change_time3 = ATPase_change_time3,
    tspan = tspan,
)
    # below is due to bug in rand where it's hard to sample 0.0 distribution
    init_conc.Citrate = 0.0
    init_conc.F26BP = 0.0
    init_conc.Phenylalanine = 0.0
    default_params.CK_Vmax = 0.0
    default_params.NDPK_Vmax = 0.0

    # pathway_Vmax = 2 * default_params.HK1_Conc * default_params.HK1_Vmax
    # pathway_Vmax = min(2 * params.HK1_Vmax * params.HK1_Conc, 2 * params.PFKP_Vmax * params.PFKP_Conc, 2 * params.ALDO_Vmax * params.ALDO_Conc)
    # pathway_Vmax = 2 * params.HK1_Vmax * params.HK1_Conc
    pathway_Vmax = min(2 * params.HK1_Vmax * params.HK1_Conc, 2 * params.PFKP_Vmax * params.PFKP_Conc)
    params.ATPase_Vmax = Initial_ATPase_Vmax_frac * pathway_Vmax
    params.ATPase_Km_ATP = 1e-6
    function affect1!(integrator)
        integrator.p.ATPase_Vmax = High_ATPase_Vmax_frac * pathway_Vmax
    end
    function affect2!(integrator)
        integrator.p.ATPase_Vmax = Low_ATPase_Vmax_frac * pathway_Vmax
    end
    function affect3!(integrator)
        integrator.p.ATPase_Vmax = Initial_ATPase_Vmax_frac * pathway_Vmax
    end

    PresetTime_cb1 = PresetTimeCallback(ATPase_change_time1, affect1!)
    PresetTime_cb2 = PresetTimeCallback(ATPase_change_time2, affect2!)
    PresetTime_cb3 = PresetTimeCallback(ATPase_change_time3, affect3!)
    cb_set = CallbackSet(PresetTime_cb1, PresetTime_cb2, PresetTime_cb3)

    init_cond_prob = ODEProblem(glycolysis_ODEs, init_conc, (0, 1e8), params)
    init_cond_sol = solve(init_cond_prob, Rodas5P(), abstol = 1e-15, reltol = 1e-8, save_everystep = false)
    new_init_cond = init_cond_sol.u[end]
    prob = ODEProblem(glycolysis_ODEs, new_init_cond, tspan, params, callback = cb_set)
    sol = solve(
        prob,
        Rodas5P(),
        abstol = 1e-15,
        reltol = 1e-8,
        saveat = [k for k = tspan[1]:((tspan[2]-tspan[1])/10_000):tspan[2]],
    )

    timepoints = sol.t
    ATPprod = [Glycolysis.conc_to_rates(conc, params).ATPprod for conc in sol.u] / pathway_Vmax
    ATPenergy = [Glycolysis.conc_to_disequilibrium_ratios(conc, params).Q_Keq_ATPase for conc in sol.u]
    ATPase = Initial_ATPase_Vmax_frac * ones(length(sol))
    ATPase[sol.t.>ATPase_change_time1.&&sol.t.<ATPase_change_time2] .= High_ATPase_Vmax_frac
    ATPase[sol.t.>ATPase_change_time2.&&sol.t.<ATPase_change_time3] .= Low_ATPase_Vmax_frac
    ATP = [conc.ATP for conc in sol.u]
    ADP = [conc.ADP for conc in sol.u]
    AMP = [conc.AMP for conc in sol.u]
    if sol.retcode == ReturnCode.Success && init_cond_sol.retcode == ReturnCode.Success
        return (
            timepoints = timepoints,
            ATPprod = ATPprod,
            ATPenergy = ATPenergy,
            ATPase = ATPase,
            ATP = ATP,
            normATP = ATP ./ (ATP + ADP + AMP),
            # params = params,
            # init_conc = init_conc,
        )
    else
        return 0.0
    end
end

default_params_ics_res = glyc_step_response(glycolysis_params, glycolysis_init_conc)

number_boostraps = 10_000
glycolysis_params_copy = deepcopy(glycolysis_params)
glycolysis_params_copy.CK_Vmax = 1.0
glycolysis_params_copy.NDPK_Vmax = 1.0
bootstrap_glycolysis_params = [
    rand.(truncated.(Normal.(glycolysis_params_copy, glycolysis_params_uncertainty); lower = 0.0)) for
    i = 1:number_boostraps
]
# next line take about 1 min to execute on an 8-core computer
bootstrap_param_step_response_list = @showprogress "Computing Parameters Sensitivity" pmap(
    x -> glyc_step_response(x, glycolysis_init_conc),
    bootstrap_glycolysis_params,
)

bootstrap_enzyme_conc = []
fold_range = 3
for i = 1:number_boostraps
    var_enzyme_glycolysis_params = deepcopy(glycolysis_params)
    var_enzyme_glycolysis_params.GLUT_Conc =
        sqrt(fold_range) .^ (-2 .+ 4 * rand()) * glycolysis_params.GLUT_Conc
    var_enzyme_glycolysis_params.HK1_Conc =
        sqrt(fold_range) .^ (-2 .+ 4 * rand()) * glycolysis_params.HK1_Conc
    var_enzyme_glycolysis_params.GPI_Conc =
        sqrt(fold_range) .^ (-2 .+ 4 * rand()) * glycolysis_params.GPI_Conc
    var_enzyme_glycolysis_params.PFKP_Conc =
        sqrt(fold_range) .^ (-2 .+ 4 * rand()) * glycolysis_params.PFKP_Conc
    var_enzyme_glycolysis_params.ALDO_Conc =
        sqrt(fold_range) .^ (-2 .+ 4 * rand()) * glycolysis_params.ALDO_Conc
    var_enzyme_glycolysis_params.TPI_Conc =
        sqrt(fold_range) .^ (-2 .+ 4 * rand()) * glycolysis_params.TPI_Conc
    var_enzyme_glycolysis_params.GAPDH_Conc =
        sqrt(fold_range) .^ (-2 .+ 4 * rand()) * glycolysis_params.GAPDH_Conc
    var_enzyme_glycolysis_params.PGK_Conc =
        sqrt(fold_range) .^ (-2 .+ 4 * rand()) * glycolysis_params.PGK_Conc
    var_enzyme_glycolysis_params.PGM_Conc =
        sqrt(fold_range) .^ (-2 .+ 4 * rand()) * glycolysis_params.PGM_Conc
    var_enzyme_glycolysis_params.PKM2_Conc =
        sqrt(fold_range) .^ (-2 .+ 4 * rand()) * glycolysis_params.PKM2_Conc
    var_enzyme_glycolysis_params.ENO_Conc =
        sqrt(fold_range) .^ (-2 .+ 4 * rand()) * glycolysis_params.ENO_Conc
    var_enzyme_glycolysis_params.LDH_Conc =
        sqrt(fold_range) .^ (-2 .+ 4 * rand()) * glycolysis_params.LDH_Conc
    var_enzyme_glycolysis_params.MCT_Conc =
        sqrt(fold_range) .^ (-2 .+ 4 * rand()) * glycolysis_params.MCT_Conc
    push!(bootstrap_enzyme_conc, var_enzyme_glycolysis_params)
end

# next line take about 0.5 min to execute on an 8-core computer
bootstrap_enzyme_conc_step_response_list = @showprogress "Computing Enzyme Level Sensitivity" pmap(
    x -> glyc_step_response(x, glycolysis_init_conc),
    bootstrap_enzyme_conc,
)

glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
glycolysis_init_conc_copy.Citrate = 1.0
glycolysis_init_conc_copy.F26BP = 1.0
glycolysis_init_conc_copy.Phenylalanine = 1.0
bootstrap_init_cond = [
    rand.(
        truncated.(
            Normal.(glycolysis_init_conc_copy, Glycolysis.glycolysis_init_conc_uncertainty);
            lower = 0.0,
        )
    ) for i = 1:number_boostraps
]

# next line take about 0.5 min to execute on an 8-core computer
bootstrap_init_cond_step_response_list = @showprogress "Computing Initial Conditions Sensitivity" pmap(
    x -> glyc_step_response(glycolysis_params, x),
    bootstrap_init_cond,
)

# Remove all the workers
for n in workers()
    rmprocs(n)
end


##
# Precalculate mean and CI

param_bootstrap = DataFrame()
for (i, res) in enumerate(bootstrap_param_step_response_list)
    if res != 0.0 #&& res.ATPprod[end] / Initial_ATPase_Vmax_frac .> 0.99
        append!(param_bootstrap, DataFrame(merge((trajectory = repeat([i], length(res.ATP)),), res)))
    end
end
nrow(param_bootstrap) / (10000 * 10004)
init_cond_bootstrap = DataFrame()
for (i, res) in enumerate(bootstrap_init_cond_step_response_list)
    if res != 0.0 #&& res.ATPprod[end] / Initial_ATPase_Vmax_frac .> 0.99
        append!(init_cond_bootstrap, DataFrame(merge((trajectory = repeat([i], length(res.ATP)),), res)))
    end
end
nrow(init_cond_bootstrap) / (10000 * 10004)

enzyme_conc_bootstrap = DataFrame()
for (i, res) in enumerate(bootstrap_enzyme_conc_step_response_list)
    if res != 0.0 #&& res.ATPprod[end] / Initial_ATPase_Vmax_frac .> 0.99
        append!(enzyme_conc_bootstrap, DataFrame(merge((trajectory = repeat([i], length(res.ATP)),), res)))
    end
end
nrow(enzyme_conc_bootstrap) / (10000 * 10004)

init_cond_bootstrap = DataFrame()
for (i, res) in enumerate(bootstrap_init_cond_step_response_list)
    if res != 0.0 #&& res.ATPprod[end] / Initial_ATPase_Vmax_frac .> 0.99
        append!(init_cond_bootstrap, DataFrame(merge((trajectory = repeat([i], length(res.ATP)),), res)))
    end
end
nrow(init_cond_bootstrap) / (10000 * 10004)

qlow(x) = percentile(x, 2.5)
qhigh(x) = percentile(x, 97.5)
processed_param_bootstrap = combine(
    groupby(param_bootstrap, :timepoints),
    Not(:timepoints) .=> median .=> names(param_bootstrap[:, Not(:timepoints)]),
    Not(:timepoints) .=>
        (x -> qlow(x)) .=> [name * "_qlow" for name in names(param_bootstrap[:, Not(:timepoints)])],
    Not(:timepoints) .=>
        (x -> qhigh(x)) .=> [name * "_qhigh" for name in names(param_bootstrap[:, Not(:timepoints)])],
)

processed_enzyme_conc_bootstrap = combine(
    groupby(enzyme_conc_bootstrap, :timepoints),
    Not(:timepoints) .=> median .=> names(enzyme_conc_bootstrap[:, Not(:timepoints)]),
    Not(:timepoints) .=>
        (x -> qlow(x)) .=> [name * "_qlow" for name in names(enzyme_conc_bootstrap[:, Not(:timepoints)])],
    Not(:timepoints) .=>
        (x -> qhigh(x)) .=>
            [name * "_qhigh" for name in names(enzyme_conc_bootstrap[:, Not(:timepoints)])],
)

processed_init_cond_bootstrap = combine(
    groupby(init_cond_bootstrap, :timepoints),
    Not(:timepoints) .=> median .=> names(init_cond_bootstrap[:, Not(:timepoints)]),
    Not(:timepoints) .=>
        (x -> qlow(x)) .=> [name * "_qlow" for name in names(init_cond_bootstrap[:, Not(:timepoints)])],
    Not(:timepoints) .=>
        (x -> qhigh(x)) .=> [name * "_qhigh" for name in names(init_cond_bootstrap[:, Not(:timepoints)])],
)


##
# Plot the results
size_inches = (6.5, 6)
size_pt = 72 .* size_inches
set_theme!(Theme(fontsize = 6, Axis = (
    xticksize = 1,
    yticksize = 1,
    # xticklabelsize = 6,
    # yticklabelsize = 6,
    yticklabelpad = 1,
    ylabelpadding = 3,
)))
fig = Figure(size = size_pt)




# Plot enzyme level bootstrapping results
# Plot ATP production
ax_ATP_prod = Axis(
    fig[1, 1:2],
    # limits = (nothing, (0.5, 2)),
    limits = (nothing, (0.0, 1.5 * maximum(default_params_ics_res.ATPase))),
    # limits = (nothing, (1/1.2, 1.2)),
    # yscale = log10,
    xlabel = "Time, min",
    ylabel = "ATP production rate\nrelative to glycolysis Vmax",
    yticklabelcolor = Makie.wong_colors()[3],
    ylabelcolor = Makie.wong_colors()[3],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    title = "Matching ATP supply and demand",
)
ax_ATPase = Axis(
    fig[1, 1:2],
    limits = (nothing, (0.0, 1.5 * maximum(default_params_ics_res.ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)

ATP_prod_line = lines!(
    ax_ATP_prod,
    processed_enzyme_conc_bootstrap.timepoints,
    (processed_enzyme_conc_bootstrap.ATPprod),
    color = Makie.wong_colors()[3],
)
band!(
    ax_ATP_prod,
    processed_enzyme_conc_bootstrap.timepoints,
    processed_enzyme_conc_bootstrap.ATPprod_qlow,
    processed_enzyme_conc_bootstrap.ATPprod_qhigh,
    color = Makie.wong_colors(0.5)[3],
)
ATPase_line = lines!(
    ax_ATPase,
    processed_enzyme_conc_bootstrap.timepoints,
    processed_enzyme_conc_bootstrap.ATPase,
    linestyle = :dot,
    color = :Black,
    linewidth = 1,
)

axislegend(
    ax_ATP_prod,
    [ATP_prod_line, ATPase_line],
    ["ATP production", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -3, -8, -6),
    patchsize = (7.5, 7.5),
)

# Plot [ATP]
ax_ATP_conc = Axis(
    fig[1, 3:4],
    limits = (nothing, (0, 14e-3)),
    xlabel = "Time, min",
    ylabel = "[ATP], mM",
    ytickformat = ys -> ["$(round(y*1000, sigdigits = 3))" for y in ys],
    yticklabelcolor = Makie.wong_colors()[1],
    ylabelcolor = Makie.wong_colors()[1],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    title = "Maintaining ATP concentration",
    width = 85,
)
ax_ATPase = Axis(
    fig[1, 3:4],
    limits = (nothing, (0.0, 1.5 * maximum(default_params_ics_res.ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
    width = 85,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)

ATP_line = lines!(
    ax_ATP_conc,
    processed_enzyme_conc_bootstrap.timepoints,
    processed_enzyme_conc_bootstrap.ATP,
    color = Makie.wong_colors()[1],
)
band!(
    ax_ATP_conc,
    processed_enzyme_conc_bootstrap.timepoints,
    processed_enzyme_conc_bootstrap.ATP_qlow,
    processed_enzyme_conc_bootstrap.ATP_qhigh,
    color = Makie.wong_colors(0.5)[1],
)
ATPase_line = lines!(
    ax_ATPase,
    processed_enzyme_conc_bootstrap.timepoints,
    processed_enzyme_conc_bootstrap.ATPase,
    linestyle = :dot,
    color = :Black,
    linewidth = 1,
)
lines!(
    ax_ATP_conc,
    [tspan[1], tspan[2]],
    repeat([glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP], 2),
    color = :grey,
    linestyle = :dash,
)
text!(
    ax_ATP_conc,
    tspan[1],
    1.01 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = :grey,
)
axislegend(
    ax_ATP_conc,
    [ATP_line, ATPase_line],
    ["[ATP]", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -3, -8, -6),
    patchsize = (7.5, 7.5),
)

# Plot energy of ATPase reaction
ax_ATP_energy = Axis(
    fig[1, 5:6],
    limits = (nothing, (0, 33)),
    xlabel = "Time, min",
    ylabel = rich("Energy of ATPase reaction, k", subscript("B"), "T"),
    # yscale = log10,
    ytickformat = ys -> ["$(Int(round(y)))" for y in ys],
    yticklabelcolor = Makie.wong_colors()[4],
    ylabelcolor = Makie.wong_colors()[4],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    title = "Maintaining ATP energy",
)
ax_ATPase = Axis(
    fig[1, 5:6],
    limits = (nothing, (0.0, 1.5 * maximum(default_params_ics_res.ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)

ATPase_energy_line = lines!(
    ax_ATP_energy,
    processed_enzyme_conc_bootstrap.timepoints,
    -log.(processed_enzyme_conc_bootstrap.ATPenergy),
    color = Makie.wong_colors()[4],
)
band!(
    ax_ATP_energy,
    processed_enzyme_conc_bootstrap.timepoints,
    -log.(processed_enzyme_conc_bootstrap.ATPenergy_qlow),
    -log.(processed_enzyme_conc_bootstrap.ATPenergy_qhigh),
    color = Makie.wong_colors(0.5)[4],
)
ATPase_line = lines!(
    ax_ATPase,
    processed_enzyme_conc_bootstrap.timepoints,
    processed_enzyme_conc_bootstrap.ATPase,
    linestyle = :dot,
    color = :Black,
    linewidth = 1,
)
axislegend(
    ax_ATP_energy,
    [ATPase_energy_line, ATPase_line],
    ["ATPase energy", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -3, -8, -6),
    patchsize = (7.5, 7.5),
)



# Plot parameter bootstrapping results
# Plot ATP production
ax_ATP_prod = Axis(
    fig[2, 1:2],
    # limits = (nothing, (0.5, 2)),
    limits = (nothing, (0.0, 1.5 * maximum(default_params_ics_res.ATPase))),
    # limits = (nothing, (1/1.2, 1.2)),
    # yscale = log10,
    xlabel = "Time, min",
    ylabel = "ATP production rate\nrelative to glycolysis Vmax",
    yticklabelcolor = Makie.wong_colors()[3],
    ylabelcolor = Makie.wong_colors()[3],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    title = "Matching ATP supply and demand",
)
ax_ATPase = Axis(
    fig[2, 1:2],
    limits = (nothing, (0.0, 1.5 * maximum(default_params_ics_res.ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)

ATP_prod_line = lines!(
    ax_ATP_prod,
    processed_param_bootstrap.timepoints,
    (processed_param_bootstrap.ATPprod),
    color = Makie.wong_colors()[3],
)
band!(
    ax_ATP_prod,
    processed_param_bootstrap.timepoints,
    processed_param_bootstrap.ATPprod_qlow,
    processed_param_bootstrap.ATPprod_qhigh,
    color = Makie.wong_colors(0.5)[3],
)
ATPase_line = lines!(
    ax_ATPase,
    processed_param_bootstrap.timepoints,
    processed_param_bootstrap.ATPase,
    linestyle = :dot,
    color = :Black,
    linewidth = 1,
)

axislegend(
    ax_ATP_prod,
    [ATP_prod_line, ATPase_line],
    ["ATP production", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -3, -8, -6),
    patchsize = (7.5, 7.5),
)

# Plot [ATP]
ax_ATP_conc = Axis(
    fig[2, 3:4],
    limits = (nothing, (0, 14e-3)),
    xlabel = "Time, min",
    ylabel = "[ATP], mM",
    ytickformat = ys -> ["$(round(y*1000, sigdigits = 3))" for y in ys],
    yticklabelcolor = Makie.wong_colors()[1],
    ylabelcolor = Makie.wong_colors()[1],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    title = "Maintaining ATP concentration",
    width = 85,
)
ax_ATPase = Axis(
    fig[2, 3:4],
    limits = (nothing, (0.0, 1.5 * maximum(default_params_ics_res.ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
    width = 85,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)

ATP_line = lines!(
    ax_ATP_conc,
    processed_param_bootstrap.timepoints,
    processed_param_bootstrap.ATP,
    color = Makie.wong_colors()[1],
)
band!(
    ax_ATP_conc,
    processed_param_bootstrap.timepoints,
    processed_param_bootstrap.ATP_qlow,
    processed_param_bootstrap.ATP_qhigh,
    color = Makie.wong_colors(0.5)[1],
)
ATPase_line = lines!(
    ax_ATPase,
    processed_param_bootstrap.timepoints,
    processed_param_bootstrap.ATPase,
    linestyle = :dot,
    color = :Black,
    linewidth = 1,
)
lines!(
    ax_ATP_conc,
    [tspan[1], tspan[2]],
    repeat([glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP], 2),
    color = :grey,
    linestyle = :dash,
)
text!(
    ax_ATP_conc,
    tspan[1],
    1.01 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = :grey,
)
axislegend(
    ax_ATP_conc,
    [ATP_line, ATPase_line],
    ["[ATP]", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -3, -8, -6),
    patchsize = (7.5, 7.5),
)

# Plot energy of ATPase reaction
ax_ATP_energy = Axis(
    fig[2, 5:6],
    limits = (nothing, (0, 33)),
    xlabel = "Time, min",
    ylabel = rich("Energy of ATPase reaction, k", subscript("B"), "T"),
    # yscale = log10,
    ytickformat = ys -> ["$(Int(round(y)))" for y in ys],
    yticklabelcolor = Makie.wong_colors()[4],
    ylabelcolor = Makie.wong_colors()[4],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    title = "Maintaining ATP energy",
)
ax_ATPase = Axis(
    fig[2, 5:6],
    limits = (nothing, (0.0, 1.5 * maximum(default_params_ics_res.ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)

ATPase_energy_line = lines!(
    ax_ATP_energy,
    processed_param_bootstrap.timepoints,
    -log.(processed_param_bootstrap.ATPenergy),
    color = Makie.wong_colors()[4],
)
band!(
    ax_ATP_energy,
    processed_param_bootstrap.timepoints,
    -log.(processed_param_bootstrap.ATPenergy_qlow),
    -log.(processed_param_bootstrap.ATPenergy_qhigh),
    color = Makie.wong_colors(0.5)[4],
)
ATPase_line = lines!(
    ax_ATPase,
    processed_param_bootstrap.timepoints,
    processed_param_bootstrap.ATPase,
    linestyle = :dot,
    color = :Black,
    linewidth = 1,
)
axislegend(
    ax_ATP_energy,
    [ATPase_energy_line, ATPase_line],
    ["ATPase energy", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -3, -8, -6),
    patchsize = (7.5, 7.5),
)




# Plot initial condition bootstrapping results
# Plot ATP production
ax_ATP_prod = Axis(
    fig[3, 1:2],
    # limits = (nothing, (0.5, 2)),
    limits = (nothing, (0.0, 1.5 * maximum(default_params_ics_res.ATPase))),
    # limits = (nothing, (1/1.2, 1.2)),
    # yscale = log10,
    xlabel = "Time, min",
    ylabel = "ATP production rate\nrelative to glycolysis Vmax",
    yticklabelcolor = Makie.wong_colors()[3],
    ylabelcolor = Makie.wong_colors()[3],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    title = "Matching ATP supply and demand",
)
ax_ATPase = Axis(
    fig[3, 1:2],
    limits = (nothing, (0.0, 1.5 * maximum(default_params_ics_res.ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)

ATP_prod_line = lines!(
    ax_ATP_prod,
    processed_init_cond_bootstrap.timepoints,
    (processed_init_cond_bootstrap.ATPprod),
    color = Makie.wong_colors()[3],
)
band!(
    ax_ATP_prod,
    processed_init_cond_bootstrap.timepoints,
    processed_init_cond_bootstrap.ATPprod_qlow,
    processed_init_cond_bootstrap.ATPprod_qhigh,
    color = Makie.wong_colors(0.5)[3],
)
ATPase_line = lines!(
    ax_ATPase,
    processed_init_cond_bootstrap.timepoints,
    processed_init_cond_bootstrap.ATPase,
    linestyle = :dot,
    color = :Black,
    linewidth = 1,
)

axislegend(
    ax_ATP_prod,
    [ATP_prod_line, ATPase_line],
    ["ATP production", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -3, -8, -6),
    patchsize = (7.5, 7.5),
)

# Plot [ATP]
ax_ATP_conc = Axis(
    fig[3, 3:4],
    limits = (nothing, (0, 1.3)),
    xlabel = "Time, min",
    ylabel = "[ATP], normalized to adenine pool",
    ytickformat = ys -> ["$(round(y, sigdigits = 3))" for y in ys],
    yticklabelcolor = Makie.wong_colors()[1],
    ylabelcolor = Makie.wong_colors()[1],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    title = "Maintaining ATP concentration",
    width = 85,
)
ax_ATPase = Axis(
    fig[3, 3:4],
    limits = (nothing, (0.0, 1.5 * maximum(default_params_ics_res.ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
    width = 85,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)

ATP_line = lines!(
    ax_ATP_conc,
    processed_init_cond_bootstrap.timepoints,
    processed_init_cond_bootstrap.normATP,
    color = Makie.wong_colors()[1],
)
band!(
    ax_ATP_conc,
    processed_init_cond_bootstrap.timepoints,
    processed_init_cond_bootstrap.normATP_qlow,
    processed_init_cond_bootstrap.normATP_qhigh,
    color = Makie.wong_colors(0.5)[1],
)
ATPase_line = lines!(
    ax_ATPase,
    processed_init_cond_bootstrap.timepoints,
    processed_init_cond_bootstrap.ATPase,
    linestyle = :dot,
    color = :Black,
    linewidth = 1,
)
# lines!(
#     ax_ATP_conc,
#     [tspan[1], tspan[2]],
#     repeat([glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP], 2),
#     color = :grey,
#     linestyle = :dash,
# )
# text!(
#     ax_ATP_conc,
#     tspan[1],
#     1.01 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
#     text = "Adenine pool size",
#     align = (:left, :bottom),
#     color = :grey,
# )
axislegend(
    ax_ATP_conc,
    [ATP_line, ATPase_line],
    ["[ATP]", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -3, -8, -6),
    patchsize = (7.5, 7.5),
)

# Plot energy of ATPase reaction
ax_ATP_energy = Axis(
    fig[3, 5:6],
    limits = (nothing, (0, 33)),
    xlabel = "Time, min",
    ylabel = rich("Energy of ATPase reaction, k", subscript("B"), "T"),
    # yscale = log10,
    ytickformat = ys -> ["$(Int(round(y)))" for y in ys],
    yticklabelcolor = Makie.wong_colors()[4],
    ylabelcolor = Makie.wong_colors()[4],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    title = "Maintaining ATP energy",
)
ax_ATPase = Axis(
    fig[3, 5:6],
    limits = (nothing, (0.0, 1.5 * maximum(default_params_ics_res.ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)

ATPase_energy_line = lines!(
    ax_ATP_energy,
    processed_init_cond_bootstrap.timepoints,
    -log.(processed_init_cond_bootstrap.ATPenergy),
    color = Makie.wong_colors()[4],
)
band!(
    ax_ATP_energy,
    processed_init_cond_bootstrap.timepoints,
    -log.(processed_init_cond_bootstrap.ATPenergy_qlow),
    -log.(processed_init_cond_bootstrap.ATPenergy_qhigh),
    color = Makie.wong_colors(0.5)[4],
)
ATPase_line = lines!(
    ax_ATPase,
    processed_init_cond_bootstrap.timepoints,
    processed_init_cond_bootstrap.ATPase,
    linestyle = :dot,
    color = :Black,
    linewidth = 1,
)
axislegend(
    ax_ATP_energy,
    [ATPase_energy_line, ATPase_line],
    ["ATPase energy", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -3, -8, -6),
    patchsize = (7.5, 7.5),
)

# Final Figure edits
colgap!(fig.layout, 10)
rowgap!(fig.layout, 5)
resize_to_layout!(fig)

label_a = fig[1, 1, TopLeft()] = Label(fig, "A", fontsize = 12, halign = :right, padding = (0, 15, 10, 0))
label_b = fig[1, 3, TopLeft()] = Label(fig, "B", fontsize = 12, halign = :right, padding = (0, 5, 10, 0))
label_c = fig[1, 5, TopLeft()] = Label(fig, "C", fontsize = 12, halign = :right, padding = (0, 5, 10, 0))
label_d = fig[2, 1, TopLeft()] = Label(fig, "D", fontsize = 12, halign = :right, padding = (0, 15, 10, 0))
label_e = fig[2, 3, TopLeft()] = Label(fig, "E", fontsize = 12, halign = :right, padding = (0, 5, 10, 0))
label_f = fig[2, 5, TopLeft()] = Label(fig, "F", fontsize = 12, halign = :right, padding = (0, 5, 10, 0))
label_g = fig[3, 1, TopLeft()] = Label(fig, "G", fontsize = 12, halign = :right, padding = (0, 15, 10, 0))
label_h = fig[3, 3, TopLeft()] = Label(fig, "H", fontsize = 12, halign = :right, padding = (0, 5, 10, 0))
label_i = fig[3, 5, TopLeft()] = Label(fig, "I", fontsize = 12, halign = :right, padding = (0, 5, 10, 0))


fig
# uncomment the line below to save the plot
save("Results/$(Dates.format(now(),"mmddyy"))_Fig2-fig_suppl1_model_behavior_and_validation_w_CI.png", fig, px_per_unit = 4)
