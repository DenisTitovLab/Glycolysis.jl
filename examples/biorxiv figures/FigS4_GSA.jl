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
ATP_AUC_hist_data =
    CSV.read("gsa_cluster_code/041124_hist_ATP_AUC_10000runs.csv", DataFrame)
ATP_energy_AUC_hist_data =
    CSV.read("gsa_cluster_code/041124_hist_ATP_energy_AUC_10000runs.csv", DataFrame)
ATP_prod_AUC_hist_data =
    CSV.read("gsa_cluster_code/041124_hist_ATP_prod_AUC_10000runs.csv", DataFrame)

#Precalculate default model Values
function find_ATP_at_ATPase_range(
    glycolysis_params,
    glycolysis_init_conc;
    n_Vmax_ATPase_values = 100,
)
    tspan = (0.0, 1e8)
    pathway_Vmax = 2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc
    ATPases = 10 .^ range(log10(0.001), log10(1.0), n_Vmax_ATPase_values) .* pathway_Vmax
    prob = ODEProblem(glycolysis_ODEs, glycolysis_init_conc, tspan, glycolysis_params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(
        ensemble_prob,
        Rodas5P(),
        EnsembleThreads(),
        trajectories = n_Vmax_ATPase_values,
        abstol = 1e-13,
        reltol = 1e-6,
        save_everystep = false,
        save_start = false,
    )
    ATP_prod_rate = [
        Glycolysis.conc_to_rates(sol.u[end], sol.prob.p).ATPprod / sol.prob.p.ATPase_Vmax for sol in sim if sol.retcode == ReturnCode.Success
    ]
    ATP_energy =
        -log.([
            Glycolysis.conc_to_disequilibrium_ratios(sol.u[end], sol.prob.p).Q_Keq_ATPase
            for sol in sim if sol.retcode == ReturnCode.Success
        ])
    ATP_conc = [sol.u[end].ATP for sol in sim if sol.retcode == ReturnCode.Success]
    return (
        ATP_conc = ATP_conc,
        ATP_prod_rate = ATP_prod_rate,
        ATP_energy = ATP_energy,
        ATPase_Vmax = ATPases /
                      (2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc),
    )
end

glycolysis_params.ATPase_Km_ATP = 1e-9
glycolysis_params_copy = deepcopy(glycolysis_params)
Complete_Model_Simulation_Data =
    @time find_ATP_at_ATPase_range(glycolysis_params_copy, glycolysis_init_conc)

Complete_model_mean_ATP =
    mean(Complete_Model_Simulation_Data.ATP_conc) /
    (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP)
Complete_model_mean_ATPprod = mean(Complete_Model_Simulation_Data.ATP_prod_rate)
Complete_model_mean_ATP_energy = mean(Complete_Model_Simulation_Data.ATP_energy)

##
#Processes sensitivity index data
GSA_ATP_AUC =
    CSV.read("gsa_cluster_code/041124_ATP_AUC_gsa_sobol_10000_3x_range.csv", DataFrame)
GSA_ATP_energy_AUC = CSV.read(
    "gsa_cluster_code/041124_ATPase_energy_AUC_gsa_sobol_10000_3x_range.csv",
    DataFrame,
)
GSA_ATP_prod_AUC =
    CSV.read("gsa_cluster_code/041124_ATPprod_AUC_gsa_sobol_10000_3x_range.csv", DataFrame)

HK1_Km_Vmax = [:HK1_K_Glucose, :HK1_K_G6P, :HK1_Conc, :HK1_Vmax]
HK1_reg = [
    :HK1_K_a_ATP,
    :HK1_β_Glucose_ATP,
    :HK1_K_a_ADP,
    :HK1_K_i_G6P_reg,
    :HK1_K_a_G6P_cat,
    :HK1_K_a_Pi,
]
PFKP_Km_Vmax = [:PFKP_K_ATP, :PFKP_K_F16BP, :PFKP_K_ADP, :PFKP_Conc, :PFKP_Vmax]
PFKP_reg = [:PFKP_L, :PFKP_K_a_F6P, :PFKP_K_i_ATP_reg, :PFKP_K_a_ADP_reg, :PFKP_K_Phosphate]
HK1_PFKP_reg = [HK1_reg; PFKP_reg]
HK1_PFKP_Km_Vmax = [HK1_Km_Vmax; PFKP_Km_Vmax]
PKM2_reg = [:PKM2_a_KdF16BP, :PKM2_i_KdF16BP]
GAPDH = [
    param for param in propertynames(glycolysis_params) if occursin("GAPDH", string(param))
]

Keqs =
    [param for param in propertynames(glycolysis_params) if occursin("Keq", string(param))]
Other = [
    param for param in propertynames(glycolysis_params) if
    param ∉ [HK1_PFKP_Km_Vmax; HK1_PFKP_reg; GAPDH]
]

params_group_names = [:All, :HK1_PFKP_Km_Vmax, :HK1_PFKP_reg, :GAPDH, :All_Other]
params_names = [
    [name for name in propertynames(glycolysis_params)],
    HK1_PFKP_Km_Vmax,
    HK1_PFKP_reg,
    GAPDH,
    Other,
]


ATP_AUC_Sens_Ind = DataFrame(
    Parameters = replace.(String.(params_group_names), "_" => " "),
    S1 = [
        sum(Iterators.filter(!isnan, GSA_ATP_AUC[1, params_name])) for
        params_name in params_names
    ],
    ST = [
        sum(Iterators.filter(!isnan, GSA_ATP_AUC[2, params_name])) for
        params_name in params_names
    ],
)
ATP_energy_AUC_Sens_Ind = DataFrame(
    Parameters = replace.(String.(params_group_names), "_" => " "),
    S1 = [
        sum(Iterators.filter(!isnan, GSA_ATP_energy_AUC[1, params_name])) for
        params_name in params_names
    ],
    ST = [
        sum(Iterators.filter(!isnan, GSA_ATP_energy_AUC[2, params_name])) for
        params_name in params_names
    ],
)
ATP_prod_AUC_Sens_Ind = DataFrame(
    Parameters = replace.(String.(params_group_names), "_" => " "),
    S1 = [
        sum(Iterators.filter(!isnan, GSA_ATP_prod_AUC[1, params_name])) for
        params_name in params_names
    ],
    ST = [
        sum(Iterators.filter(!isnan, GSA_ATP_prod_AUC[2, params_name])) for
        params_name in params_names
    ],
)

##
# Plot the results
size_inches = (6.5, 5)
size_pt = 72 .* size_inches
set_theme!(
    Theme(
        fontsize = 6,
        Axis = (
            xticksize = 1,
            yticksize = 1,
            # xticklabelsize = 5,
            # yticklabelsize = 5,
            yticklabelpad = 1,
            ylabelpadding = 3,
        ),
    ),
)
fig = Figure(size = size_pt)

"function to calculate error bar position for stacked bar plots"
function compute_x(x, dodge; width = 1.0, gap = 0.2, dodge_gap = 0.03)
    scale_width(dodge_gap, n_dodge) = (1 - (n_dodge - 1) * dodge_gap) / n_dodge
    function shift_dodge(i, dodge_width, dodge_gap)
        (dodge_width - 1) / 2 + (i - 1) * (dodge_width + dodge_gap)
    end
    width *= 1 - gap
    n_dodge = maximum(dodge)
    dodge_width = scale_width(dodge_gap, n_dodge)
    shifts = shift_dodge.(dodge, dodge_width, dodge_gap)
    return x .+ width .* shifts
end

#Plot [ATP] AUC histogram data
ax_ATP_AUC_hist = Axis(
    fig[1, 2],
    limits = ((nothing), (0, 0.10)),
    title = "Variance of mean [ATP] in Fig. 3E\n at a 9x range of model parameters",
    xlabel = "mean [ATP], normalized to adenine pool size",
    ylabel = "Fraction of counts",
    xtickformat = xs -> ["$(round(x,sigdigits=1))" for x in xs],
)
hist_all_params = hist!(
    ax_ATP_AUC_hist,
    ATP_AUC_hist_data.all_params,
    color = (Makie.wong_colors()[3], 0.5),
    normalization = :probability,
    bins = range(0, 1, 50),
)
text!(
    ax_ATP_AUC_hist,
    0.0,
    0.08;
    text = "CV=$(round(std(ATP_AUC_hist_data.all_params) / mean(ATP_AUC_hist_data.all_params), sigdigits=2))",
    fontsize = 8,
)
vlines!(ax_ATP_AUC_hist, Complete_model_mean_ATP, color = :red, linestyle = :dash)
text!(
    ax_ATP_AUC_hist,
    Complete_model_mean_ATP,
    0.02;
    text = "Complete Model",
    color = :red,
    fontsize = 8,
    rotation = π / 2,
)

# Plot sensitivity indexes
barplot_df = stack(
    ATP_AUC_Sens_Ind[2:end, [:Parameters, :S1, :ST]],
    [:S1, :ST];
    variable_name = :sens_ind,
    value_name = :sens_ind_value,
)
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

ax_Sens_Ind = Axis(
    fig[2, 2],
    limits = ((nothing), (0, 0.9)),
    title = "Sens. indxs for mean [ATP] in Fig. 3E\n at a 9x range of model parameters",
    ylabel = "Sensitivity Index Values",
    xticks = (
        1:length(unique(barplot_df.Parameters)),
        replace.(
            unique(barplot_df.Parameters),
            "HK1 PFKP Km Vmax" => "HK1 & PFKP\nKm & Vmax",
            "HK1 PFKP" => "HK1 & PFKP",
            "Km Vmax" => "Km & Vmax",
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
vspan!(ax_Sens_Ind, 1.5, 2.5; ymin = 0.0, ymax = 0.35, color = (:red, 0.1))


labels = ["S1", "ST"]
elements = [PolyElement(polycolor = colors[i]) for i = 1:length(labels)]
axislegend(
    ax_Sens_Ind,
    elements,
    labels,
    # [hist_all_params, hist_fixed_HK_Km_Vmax],
    # ["All glycolysis_params varied", "All glycolysis_params varied\nexcept HK1 Km, Vmax"],
    position = (1.1, 1.1),
    rowgap = 1,
    framevisible = false,
    patchsize = (10, 5),
)


#Plot AUC histogram data
ax_ATP_AUC_hist = Axis(
    fig[1, 3],
    limits = ((nothing), (0, 0.10)),
    title = "Variance of mean Energy in Fig. 3F\n at a 9x range of model parameters",
    xlabel = "mean Energy",
    ylabel = "Fraction of counts",
    xtickformat = xs -> ["$(round(x,sigdigits=1))" for x in xs],
)
hist_all_params = hist!(
    ax_ATP_AUC_hist,
    ATP_energy_AUC_hist_data.all_params,
    color = (Makie.wong_colors()[3], 0.5),
    normalization = :probability,
    bins = range(0, 25, 50),
)
text!(
    ax_ATP_AUC_hist,
    0.0,
    0.08;
    text = "CV=$(round(std(ATP_energy_AUC_hist_data.all_params) / mean(ATP_energy_AUC_hist_data.all_params), sigdigits=2))",
    fontsize = 8,
)
vlines!(ax_ATP_AUC_hist, Complete_model_mean_ATP_energy, color = :red, linestyle = :dash)
text!(
    ax_ATP_AUC_hist,
    Complete_model_mean_ATP_energy,
    0.02;
    text = "Complete Model",
    color = :red,
    fontsize = 8,
    rotation = π / 2,
)

# Plot sensitivity indexes
barplot_df = stack(
    ATP_energy_AUC_Sens_Ind[2:end, [:Parameters, :S1, :ST]],
    [:S1, :ST];
    variable_name = :sens_ind,
    value_name = :sens_ind_value,
)
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

ax_Sens_Ind = Axis(
    fig[2, 3],
    limits = ((nothing), (0, 0.9)),
    title = "Sens. indxs for mean Energy in Fig. 3F\n at a 9x range of model parameters",
    ylabel = "Sensitivity Index Values",
    xticks = (
        1:length(unique(barplot_df.Parameters)),
        replace.(
            unique(barplot_df.Parameters),
            "HK1 PFKP Km Vmax" => "HK1 & PFKP\nKm & Vmax",
            "HK1 PFKP" => "HK1 & PFKP",
            "Km Vmax" => "Km & Vmax",
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
vspan!(ax_Sens_Ind, 1.5, 2.5; ymin = 0.0, ymax = 0.35, color = (:red, 0.1))


labels = ["S1", "ST"]
elements = [PolyElement(polycolor = colors[i]) for i = 1:length(labels)]
axislegend(
    ax_Sens_Ind,
    elements,
    labels,
    # [hist_all_params, hist_fixed_HK_Km_Vmax],
    # ["All glycolysis_params varied", "All glycolysis_params varied\nexcept HK1 Km, Vmax"],
    position = (1.1, 1.1),
    rowgap = 1,
    framevisible = false,
    patchsize = (10, 5),
)


#Plot ATP prod AUC histogram data
ax_ATP_AUC_hist = Axis(
    fig[1, 1],
    limits = ((nothing), (0, 0.10)),
    title = "Variance of mean ATP prod/ATPase in Fig. 3D\n at a 9x range of model parameters",
    xlabel = "mean ATP prod/ATPase",
    ylabel = "Fraction of counts",
    xtickformat = xs -> ["$(round(x,sigdigits=1))" for x in xs],
)
hist_all_params = hist!(
    ax_ATP_AUC_hist,
    ATP_prod_AUC_hist_data.all_params,
    color = (Makie.wong_colors()[3], 0.5),
    normalization = :probability,
    bins = range(0, 1, 50),
)
text!(
    ax_ATP_AUC_hist,
    0.0,
    0.08;
    text = "CV=$(round(std(ATP_prod_AUC_hist_data.all_params) / mean(ATP_prod_AUC_hist_data.all_params), sigdigits=2))",
    fontsize = 8,
)
vlines!(ax_ATP_AUC_hist, Complete_model_mean_ATPprod, color = :red, linestyle = :dash)
text!(
    ax_ATP_AUC_hist,
    Complete_model_mean_ATPprod,
    0.02;
    text = "Complete Model",
    color = :red,
    fontsize = 8,
    rotation = π / 2,
)

# Plot sensitivity indexes
barplot_df = stack(
    ATP_prod_AUC_Sens_Ind[2:end, [:Parameters, :S1, :ST]],
    [:S1, :ST];
    variable_name = :sens_ind,
    value_name = :sens_ind_value,
)
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

ax_Sens_Ind = Axis(
    fig[2, 1],
    limits = ((nothing), (0, 0.9)),
    title = "Sens. indxs for mean ATP prod/ATPase in Fig. 3D\n at a 9x range of model parameters",
    ylabel = "Sensitivity Index Values",
    xticks = (
        1:length(unique(barplot_df.Parameters)),
        replace.(
            unique(barplot_df.Parameters),
            "HK1 PFKP Km Vmax" => "HK1 & PFKP\nKm & Vmax",
            "HK1 PFKP" => "HK1 & PFKP",
            "Km Vmax" => "Km & Vmax",
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
vspan!(ax_Sens_Ind, 1.5, 2.5; ymin = 0.0, ymax = 0.35, color = (:red, 0.1))


labels = ["S1", "ST"]
elements = [PolyElement(polycolor = colors[i]) for i = 1:length(labels)]
axislegend(
    ax_Sens_Ind,
    elements,
    labels,
    # [hist_all_params, hist_fixed_HK_Km_Vmax],
    # ["All glycolysis_params varied", "All glycolysis_params varied\nexcept HK1 Km, Vmax"],
    position = (1.1, 1.1),
    rowgap = 1,
    framevisible = false,
    patchsize = (10, 5),
)

colgap!(fig.layout, 7.5)
rowgap!(fig.layout, 5)

label_a =
    fig[1, 1, TopLeft()] =
        Label(fig, "A", fontsize = 12, halign = :right, padding = (0, 15, 5, 0))
label_b =
    fig[1, 2, TopLeft()] =
        Label(fig, "B", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_c =
    fig[1, 3, TopLeft()] =
        Label(fig, "C", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_d =
    fig[2, 1, TopLeft()] =
        Label(fig, "D", fontsize = 12, halign = :right, padding = (0, 15, 5, 0))
label_e =
    fig[2, 2, TopLeft()] =
        Label(fig, "E", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_f =
    fig[2, 3, TopLeft()] =
        Label(fig, "F", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))

fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_FigS4_global_sensitivity_analysis.png", fig, px_per_unit = 4)
