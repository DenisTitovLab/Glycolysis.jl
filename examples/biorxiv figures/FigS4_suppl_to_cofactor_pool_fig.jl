using Glycolysis
using DifferentialEquations, ProgressMeter
using CairoMakie, Dates, Printf, Statistics, StatsBase
using DataFrames, CSV

#= 
Global sensitivity analysis (GSA) is a computationally intensive task that
needs to performed on a computing cluster with many cores. Here we load the results
of such a computation for plotting. The code to run this analysis on the cluster is 
available in gsa_cluster_code/ folder.
=#
#Download sens indexes for init cond and pool sizes
GSA_ATP_AUC_cofactor_pools_sens_ind =
    CSV.read("gsa_cluster_code/050723_ATP_AUC_gsa_sobol_cofactor_pool_20000_3x_range.csv", DataFrame)
hist_ATP_AUC_cofactor_pools =
    CSV.read("gsa_cluster_code/050623_hist_ATP_AUC_10000_runs_3x_pool_sizes.csv", DataFrame)

##
# Precalculate and save outputs for model with different Pyruvate, NAD+NADH, Lactate_media, Glucose_media

function find_ATP_at_ATPase_range(glycolysis_params, glycolysis_init_conc; n_Vmax_ATPase_values = 1000)
    tspan = (0.0, 1e8)
    pathway_Vmax = 2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc
    ATPases = 10 .^ range(log10(0.003), log10(0.5), n_Vmax_ATPase_values) .* pathway_Vmax
    prob = ODEProblem(glycolysis_ODEs, glycolysis_init_conc, tspan, glycolysis_params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(
        ensemble_prob,
        Rodas4(),
        EnsembleThreads(),
        trajectories = n_Vmax_ATPase_values,
        abstol = 1e-12,
        reltol = 1e-5,
        save_everystep = false,
        save_start = false,
    )
    ATP_conc = [sol.u[end].ATP for sol in sim]
    return (
        ATP_conc = ATP_conc,
        ATPase_Vmax = ATPases / (2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc),
    )
end

function find_ATP_at_ATPase_range_const_Pi(
    glycolysis_params,
    glycolysis_init_conc;
    n_Vmax_ATPase_values = 1000,
)
    function glycolysis_ODEs_const_Pi(ds, s, glycolysis_params, t)
        ds.Glucose_media = 0
        ds.Glucose =
            Glycolysis.rate_GLUT(s.Glucose_media, s.Glucose, glycolysis_params) -
            Glycolysis.rate_HK1(s.Glucose, s.G6P, s.ATP, s.ADP, s.Phosphate, glycolysis_params)
        ds.G6P =
            Glycolysis.rate_HK1(s.Glucose, s.G6P, s.ATP, s.ADP, s.Phosphate, glycolysis_params) -
            Glycolysis.rate_GPI(s.G6P, s.F6P, glycolysis_params)
        ds.F6P = (
            Glycolysis.rate_GPI(s.G6P, s.F6P, glycolysis_params) - Glycolysis.rate_PFKP(
                s.F6P,
                s.ATP,
                s.F16BP,
                s.ADP,
                s.Phosphate,
                s.Citrate,
                s.F26BP,
                glycolysis_params,
            )
        )
        ds.F16BP = (
            Glycolysis.rate_PFKP(
                s.F6P,
                s.ATP,
                s.F16BP,
                s.ADP,
                s.Phosphate,
                s.Citrate,
                s.F26BP,
                glycolysis_params,
            ) - Glycolysis.rate_ALDO(s.F16BP, s.GAP, s.DHAP, glycolysis_params)
        )
        ds.GAP = (
            Glycolysis.rate_ALDO(s.F16BP, s.GAP, s.DHAP, glycolysis_params) +
            Glycolysis.rate_TPI(s.GAP, s.DHAP, glycolysis_params) -
            Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, glycolysis_params)
        )
        ds.DHAP =
            Glycolysis.rate_ALDO(s.F16BP, s.GAP, s.DHAP, glycolysis_params) -
            Glycolysis.rate_TPI(s.GAP, s.DHAP, glycolysis_params)
        ds.BPG =
            Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, glycolysis_params) -
            Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, glycolysis_params)
        ds.ThreePG =
            Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, glycolysis_params) -
            Glycolysis.rate_PGM(s.ThreePG, s.TwoPG, glycolysis_params)
        ds.TwoPG =
            Glycolysis.rate_PGM(s.ThreePG, s.TwoPG, glycolysis_params) -
            Glycolysis.rate_ENO(s.TwoPG, s.PEP, glycolysis_params)
        ds.PEP =
            Glycolysis.rate_ENO(s.TwoPG, s.PEP, glycolysis_params) -
            Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, glycolysis_params)
        ds.Pyruvate =
            Glycolysis.rate_PKM2(
                s.PEP,
                s.ADP,
                s.Pyruvate,
                s.ATP,
                s.F16BP,
                s.Phenylalanine,
                glycolysis_params,
            ) - Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, glycolysis_params)
        ds.Lactate =
            Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, glycolysis_params) -
            Glycolysis.rate_MCT(s.Lactate, s.Lactate_media, glycolysis_params)
        ds.Lactate_media = 0
        ds.ATP = (
            -Glycolysis.rate_HK1(s.Glucose, s.G6P, s.ATP, s.ADP, s.Phosphate, glycolysis_params) -
            Glycolysis.rate_PFKP(
                s.F6P,
                s.ATP,
                s.F16BP,
                s.ADP,
                s.Phosphate,
                s.Citrate,
                s.F26BP,
                glycolysis_params,
            ) +
            Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, glycolysis_params) +
            Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, glycolysis_params) - Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, glycolysis_params) +
            Glycolysis.rate_AK(s.ATP, s.ADP, s.AMP, glycolysis_params)
        )
        ds.ADP = (
            Glycolysis.rate_HK1(s.Glucose, s.G6P, s.ATP, s.ADP, s.Phosphate, glycolysis_params) +
            Glycolysis.rate_PFKP(
                s.F6P,
                s.ATP,
                s.F16BP,
                s.ADP,
                s.Phosphate,
                s.Citrate,
                s.F26BP,
                glycolysis_params,
            ) - Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, glycolysis_params) -
            Glycolysis.rate_PKM2(
                s.PEP,
                s.ADP,
                s.Pyruvate,
                s.ATP,
                s.F16BP,
                s.Phenylalanine,
                glycolysis_params,
            ) + Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, glycolysis_params) -
            2 * Glycolysis.rate_AK(s.ATP, s.ADP, s.AMP, glycolysis_params)
        )
        ds.AMP = Glycolysis.rate_AK(s.ATP, s.ADP, s.AMP, glycolysis_params)
        # ds.Phosphate =
        #     Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, glycolysis_params) -
        #     Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, glycolysis_params)
        ds.Phosphate = 0

        ds.NAD =
            Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, glycolysis_params) -
            Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, glycolysis_params)
        ds.NADH =
            Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, glycolysis_params) -
            Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, glycolysis_params)
        ds.F26BP = 0
        ds.Citrate = 0
        ds.Phenylalanine = 0
    end
    tspan = (0.0, 1e8)
    pathway_Vmax = 2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc
    ATPases = 10 .^ range(log10(0.003), log10(0.5), n_Vmax_ATPase_values) .* pathway_Vmax
    prob = ODEProblem(glycolysis_ODEs_const_Pi, glycolysis_init_conc, tspan, glycolysis_params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(
        ensemble_prob,
        Rodas4(),
        # RadauIIA5(),
        EnsembleThreads(),
        trajectories = n_Vmax_ATPase_values,
        abstol = 1e-12,
        reltol = 1e-5,
        save_everystep = false,
        save_start = false,
    )
    ATP_prod_rate = [
        Glycolysis.conc_to_rates(sol.u[end], sol.prob.p).ATPprod / sol.prob.p.ATPase_Vmax for
        sol in sim if Int(sol.retcode) == 1 && all(sol.u[end] .>= 0.0)
    ]
    ATP_energy = [
        Glycolysis.conc_to_disequilibrium_ratios(sol.u[end], sol.prob.p).Q_Keq_ATPase for
        sol in sim if Int(sol.retcode) == 1 && all(sol.u[end] .>= 0.0)
    ]
    ATP_conc = [sol.u[end].ATP for sol in sim if Int(sol.retcode) == 1 && all(sol.u[end] .>= 0.0)]
    ATPase_Vmax = [
        ATPase for (i, ATPase) in enumerate(ATPases) if Int(sim[i].retcode) == 1 && all(sim[i].u[end] .>= 0.0)
    ]
    return (
        ATP_conc = ATP_conc,
        ATP_prod_rate = ATP_prod_rate,
        ATP_energy = ATP_energy,
        ATPase_Vmax = ATPase_Vmax / (2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc),
    )
end

Complete_Model_Simulation_Data = find_ATP_at_ATPase_range(glycolysis_params, glycolysis_init_conc)

# init_cond_list_NADpool = []
# model_names_list_NADpool = []
# for multiplier in [10.0, 3.0, 1.0, 0.33, 0.1]
#     glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
#     glycolysis_init_conc_copy.NAD = multiplier * glycolysis_init_conc.NAD
#     glycolysis_init_conc_copy.NADH = multiplier * glycolysis_init_conc.NADH
#     push!(init_cond_list_NADpool, glycolysis_init_conc_copy)
#     push!(model_names_list_NADpool, "$(multiplier)x")
# end
# Simulation_Data_NADpool = []
# for model_glycolysis_init_conc in init_cond_list_NADpool
#     res = find_ATP_at_ATPase_range(glycolysis_params, model_glycolysis_init_conc)
#     push!(Simulation_Data_NADpool, res)
# end

init_cond_list_Adenine = []
model_names_list_Adenine = []
# ATP_multiplier_normalizer = [4.0, 2.0, 1.0, 0.5, 0.25]
ATP_multiplier_normalizer = [2.0, 1.0, 0.5]
Adenine_pool_size = glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP
for multiplier in ATP_multiplier_normalizer
    glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
    # glycolysis_init_conc_copy.Pyruvate = multiplier * glycolysis_init_conc.Pyruvate
    glycolysis_init_conc_copy.ATP = multiplier * glycolysis_init_conc.ATP
    glycolysis_init_conc_copy.ADP = multiplier * glycolysis_init_conc.ADP
    glycolysis_init_conc_copy.AMP = multiplier * glycolysis_init_conc.AMP

    push!(init_cond_list_Adenine, glycolysis_init_conc_copy)
    push!(model_names_list_Adenine, "$(Int(round(multiplier * Adenine_pool_size * 1000)))mM")
end
Simulation_Data_Adenine = []
for model_glycolysis_init_conc in init_cond_list_Adenine
    res = find_ATP_at_ATPase_range(glycolysis_params, model_glycolysis_init_conc)
    push!(Simulation_Data_Adenine, res)
end

init_cond_list_Pi = []
model_names_list_Pi = []
for multiplier in [2.0, 1.0, 0.5]
    glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
    glycolysis_init_conc_copy.G6P = multiplier * glycolysis_init_conc.G6P
    glycolysis_init_conc_copy.F6P = multiplier * glycolysis_init_conc.F6P
    glycolysis_init_conc_copy.F16BP = multiplier * glycolysis_init_conc.F16BP
    glycolysis_init_conc_copy.GAP = multiplier * glycolysis_init_conc.GAP
    glycolysis_init_conc_copy.DHAP = multiplier * glycolysis_init_conc.DHAP
    glycolysis_init_conc_copy.BPG = multiplier * glycolysis_init_conc.BPG
    glycolysis_init_conc_copy.ThreePG = multiplier * glycolysis_init_conc.ThreePG
    glycolysis_init_conc_copy.TwoPG = multiplier * glycolysis_init_conc.TwoPG
    glycolysis_init_conc_copy.PEP = multiplier * glycolysis_init_conc.PEP
    # glycolysis_init_conc_copy.Phosphate = 10 * glycolysis_init_conc.Phosphate
    glycolysis_init_conc_copy.Phosphate = 5e-3
    push!(init_cond_list_Pi, glycolysis_init_conc_copy)
    push!(model_names_list_Pi, "$(multiplier)x")
end
Simulation_Data_Pi = []
for model_glycolysis_init_conc in init_cond_list_Pi
    res = find_ATP_at_ATPase_range_const_Pi(
        glycolysis_params,
        model_glycolysis_init_conc;
        n_Vmax_ATPase_values = 100,
    )
    push!(Simulation_Data_Pi, res)
end

##

# Plot results

size_inches = (6.5, 5)
size_pt = 72 .* size_inches
set_theme!(Theme(fontsize = 6, Axis = (
    xticksize = 1,
    yticksize = 1,
    # xticklabelsize = 5,
    # yticklabelsize = 5,
    yticklabelpad = 1,
    ylabelpadding = 3,
)))
fig = Figure(resolution = size_pt)


#Plot hist ans sens indexes for cofactor pools
#Plot [ATP] AUC histogram data
ax_ATP_AUC_hist = Axis(
    fig[1, 1],
    limits = ((nothing), (0, nothing)),
    title = "Variance of ∑[ATP] at a 9x\nrange of cofactor pool sizes",
    xlabel = "∑[ATP], normalized to adenine pool size",
    ylabel = "Fraction of counts",
    xtickformat = xs -> ["$(round(x,sigdigits=1))" for x in xs],
)

hist_all_glycolysis_params = hist!(
    ax_ATP_AUC_hist,
    hist_ATP_AUC_cofactor_pools.all_params,
    color = (Makie.wong_colors()[3], 0.5),
    normalization = :probability,
    bins = range(0, 1, 30),
)
text!(
    ax_ATP_AUC_hist,
    0.5,
    0.08;
    text = "CV=$(round(std(hist_ATP_AUC_cofactor_pools.all_params) / mean(hist_ATP_AUC_cofactor_pools.all_params), sigdigits=2))",
    fontsize = 8,
)
# Plot sensitivity indexes
barplot_df = GSA_ATP_AUC_cofactor_pools_sens_ind
barplot_df = rename(barplot_df, "Row" => "sens_ind")
barplot_df = stack(barplot_df; variable_name = :Parameters, value_name = :sens_ind_value)
barplot_df = combine(
    groupby(barplot_df, :Parameters),
    :Parameters,
    :sens_ind,
    :sens_ind_value,
    groupindices => :group_index,
)
barplot_df = combine(
    groupby(barplot_df, :sens_ind),
    :Parameters,
    :sens_ind,
    :sens_ind_value,
    :group_index,
    groupindices => :dodge_index,
)

# replace.(String.(glycolysis_params_group_names), " " => "\n", count=1)
ax_Sens_Ind = Axis(
    fig[1, 2],
    limits = ((nothing), (0, 1.0)),
    title = "Sens. indxs for ∑[ATP] at a 9x\nrange of cofactor pool sizes",
    ylabel = "Sensitivity Index Values",
    xticks = (
        1:length(unique(barplot_df.Parameters)),
        replace.(
            unique(barplot_df.Parameters),
            "Adenine_pool_scaler" => "Adenine\npool",
            "NAD_pool_scaler" => "NAD(H)\npool",
            "Pi_scaler" => "Phosphate\npool",
        ),
    ),
    xticklabelrotation = π / 2,
)
colors = Makie.wong_colors()
barplot!(
    ax_Sens_Ind,
    barplot_df.group_index,
    barplot_df.sens_ind_value,
    dodge = barplot_df.dodge_index,
    color = colors[barplot_df.dodge_index],
)
labels = ["S1", "ST"]
elements = [PolyElement(polycolor = colors[i]) for i = 1:length(labels)]
axislegend(
    ax_Sens_Ind,
    elements,
    labels,
    # [hist_all_glycolysis_params, hist_fixed_HK_Km_Vmax],
    # ["All glycolysis_params varied", "All glycolysis_params varied\nexcept HK1 Km, Vmax"],
    position = :lt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
    patchsize = (10, 5),
)

#Plot Pi pool variation
const_Pi_Panel = fig[2, 1] = GridLayout()
ax_ATP_conc = Axis(
    const_Pi_Panel[1, 1],
    limits = ((0.003, 0.3), (0.0, 0.010)),
    xscale = log10,
    xlabel = rich("ATPase, % of pathway V",subscript("max")),
    ylabel = "[ATP],mM",
    title = "Effect of Phosphate pool changes\nat constant [Phosphate] = 5mM",
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3]),
    xtickformat = xs -> [x*100 >=1 ? "$(Int(round(x*100, sigdigits=1)))%" : "$(round(x*100, sigdigits=1))%" for x in xs],
    ytickformat = ys -> ["$(round(y*1000, sigdigits=3))" for y in ys],
)

color = Makie.wong_colors()
linewidth = 5
for (i, model_res) in enumerate(Simulation_Data_Pi)
    lines!(
        ax_ATP_conc,
        model_res.ATPase_Vmax,
        model_res.ATP_conc,
        color = color[i],
        label = model_names_list_Pi[i],
        linewidth = linewidth
        # label = replace(
        #     model_names_list_NADpool[i],
        #     "$(Int(round(Pi_pool*1000)))mM" => "$(Int(round(Pi_pool*1000)))mM\nDefault Model",
        # ),
    )
    linewidth *= 0.5
end
lines!(
    ax_ATP_conc,
    [0.001, 1.0],
    repeat([glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP], 2),
    color = :grey,
    linestyle = :dash,
)
text!(
    ax_ATP_conc,
    0.004,
    1.03 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = :grey,
)
Legend(
    const_Pi_Panel[1, 1, Right()],
    ax_ATP_conc,
    "Phosphate Pool\nfold-change\nrelative to default\nmodel",
    rowgap = 2,
    padding = (5, 2, 4, 2),
    patchsize = (10, 5),
    patchlabelgap = 5,
    framevisible = false,
    halign = :center,
    valign = :center,
)
rowgap!(const_Pi_Panel, 5)
colgap!(const_Pi_Panel, 5)

#Plot model with more or less ATP pool
Adenine_Panel = fig[2, 2] = GridLayout()
ax_ATP_conc = Axis(
    Adenine_Panel[1, 1],
    limits = ((0.003, 0.3), (0.0, 0.010)),
    xscale = log10,
    xlabel = rich("ATPase, % of pathway V", subscript("max")),
    ylabel = "[ATP],mM (Normalized by Pool Size Change)",
    title = "Effect of adenine pool size\non maintenance of ATP",
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3]),
    xtickformat = xs -> [x*100 >=1 ? "$(Int(round(x*100, sigdigits=1)))%" : "$(round(x*100, sigdigits=1))%" for x in xs],
    ytickformat = ys -> ["$(round(y*1000, sigdigits=3))" for y in ys],
)

color = Makie.wong_colors()

for (i, model_res) in enumerate(Simulation_Data_Adenine)
    supported_ATPase_range = model_res.ATPase_Vmax[model_res.ATP_conc ./ ATP_multiplier_normalizer[i] .> 0.5*(glycolysis_init_conc.ATP+glycolysis_init_conc.ADP+glycolysis_init_conc.AMP)]
    supported_ATPase_fold_range = round(supported_ATPase_range[end]/supported_ATPase_range[1], sigdigits = 2)
    lines!(
        ax_ATP_conc,
        model_res.ATPase_Vmax,
        model_res.ATP_conc ./ ATP_multiplier_normalizer[i],
        color = model_names_list_Adenine[i] != "9mM" ? color[i] : :Black,
        label = replace(
            model_names_list_Adenine[i],
            "$(Int(round(1.0 * Adenine_pool_size * 1000)))mM" => "$(Int(round(1.0 * Adenine_pool_size * 1000)))mM\n($(supported_ATPase_fold_range)-fold)\nDefault Model",
            "mM" => "mM\n($(supported_ATPase_fold_range)-fold)"
        ),
    )
end
lines!(
    ax_ATP_conc,
    [0.001, 1.0],
    repeat([glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP], 2),
    color = :grey,
    linestyle = :dash,
)
text!(
    ax_ATP_conc,
    0.004,
    1.01 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = :grey,
)

Legend(
    Adenine_Panel[1, 1, Right()],
    ax_ATP_conc,
    "Adenine Pool, mM\n(fold ATPase range where\n[ATP] > 1/2 Adenine Pool)",
    rowgap = 2,
    padding = (5, 2, 4, 2),
    patchsize = (10, 5),
    patchlabelgap = 5,
    titlegap = 3,
    framevisible = false,
    # tellheight = false,
    # tellwidth = false,
    halign = :center,
    valign = :center,
)
rowgap!(Adenine_Panel, 5)
colgap!(Adenine_Panel, 5)


colgap!(fig.layout, 5)
rowgap!(fig.layout, 5)


label_a = fig[1, 1, TopLeft()] = Label(fig, "A", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_b = fig[1, 2, TopLeft()] = Label(fig, "B", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_c = fig[2, 1, TopLeft()] = Label(fig, "C", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_d = fig[2, 2, TopLeft()] = Label(fig, "D", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))

fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_FigS4_gsa_and_sims_for_cofactor_pool_sizes.png", fig, px_per_unit = 4)
