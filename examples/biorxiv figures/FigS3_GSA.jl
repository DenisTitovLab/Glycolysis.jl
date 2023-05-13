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
ATP_AUC_data = CSV.read("gsa_cluster_code/050623_hist_ATP_AUC_10000runs_only_HK_PFK_Vmax.csv", DataFrame)
ATP_energy_AUC_data =
    CSV.read("gsa_cluster_code/050623_hist_ATP_energy_AUC_10000runs_only_HK_PFK_Vmax.csv", DataFrame)
ATP_prod_AUC_data = CSV.read("gsa_cluster_code/050623_hist_ATP_prod_AUC_10K_only_HK_PFK_Vmax.csv", DataFrame)

std(ATP_AUC_data.all_params) / mean(ATP_AUC_data.all_params)
std(ATP_energy_AUC_data.all_params) / mean(ATP_energy_AUC_data.all_params)
std(ATP_prod_AUC_data.all_params) / mean(ATP_prod_AUC_data.all_params)

ATP_AUC_Sens_Ind = CSV.read("gsa_cluster_code/050723_ATP_AUC_gsa_sobol_sens_ind_HK_PFK_Vmax.csv", DataFrame)
ATP_energy_AUC_Sens_Ind =
    CSV.read("gsa_cluster_code/050723_ATP_energy_AUC_gsa_sobol_sens_ind_HK_PFK_Vmax.csv", DataFrame)
ATP_prod_AUC_Sens_Ind =
    CSV.read("gsa_cluster_code/050723_ATP_prod_AUC_gsa_sobol_sens_ind_HK_PFK_Vmax.csv", DataFrame)


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
fig = Figure(resolution = size_pt)


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
    ATP_AUC_data.all_params,
    color = (Makie.wong_colors()[3], 0.5),
    normalization = :probability,
    bins = range(0, 1, 50),
)
text!(
    ax_ATP_AUC_hist,
    0.0,
    0.08;
    text = "CV=$(round(std(ATP_AUC_data.all_params) / mean(ATP_AUC_data.all_params), sigdigits=2))",
    fontsize = 8,
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
# replace.(String.(params_group_names), " " => "\n", count=1)

ax_Sens_Ind = Axis(
    fig[2, 2],
    limits = ((nothing), (0, 0.9)),
    title = "Sens. indxs for mean [ATP] in Fig. 3E\n at a 9x range of model parameters",
    ylabel = "Sensitivity Index Values",
    xticks = (1:length(unique(barplot_df.Parameters)), replace.(unique(barplot_df.Parameters), "HK1 PFKP" => "HK1 & PFKP", "Km Vmax" => "Km & Vmax")),
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
vspan!(ax_Sens_Ind, 1.5, 2.5; ymin = 0.0, ymax = 0.55, color = (:red, 0.1))


labels = ["S1", "ST"]
elements = [PolyElement(polycolor = colors[i]) for i = 1:length(labels)]
axislegend(
    ax_Sens_Ind,
    elements,
    labels,
    # [hist_all_params, hist_fixed_HK_Km_Vmax],
    # ["All glycolysis_params varied", "All glycolysis_params varied\nexcept HK1 Km, Vmax"],
    position = :lt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
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
    ATP_energy_AUC_data.all_params,
    color = (Makie.wong_colors()[3], 0.5),
    normalization = :probability,
    bins = range(0, 25, 50),
)
text!(
    ax_ATP_AUC_hist,
    0.0,
    0.08;
    text = "CV=$(round(std(ATP_energy_AUC_data.all_params) / mean(ATP_energy_AUC_data.all_params), sigdigits=2))",
    fontsize = 8,
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
    xticks = (1:length(unique(barplot_df.Parameters)), replace.(unique(barplot_df.Parameters), "HK1 PFKP" => "HK1 & PFKP", "Km Vmax" => "Km & Vmax")),
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
vspan!(ax_Sens_Ind, 1.5, 2.5; ymin = 0.0, ymax = 0.55, color = (:red, 0.1))


labels = ["S1", "ST"]
elements = [PolyElement(polycolor = colors[i]) for i = 1:length(labels)]
axislegend(
    ax_Sens_Ind,
    elements,
    labels,
    # [hist_all_params, hist_fixed_HK_Km_Vmax],
    # ["All glycolysis_params varied", "All glycolysis_params varied\nexcept HK1 Km, Vmax"],
    position = :lt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
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
    ATP_prod_AUC_data.all_params,
    color = (Makie.wong_colors()[3], 0.5),
    normalization = :probability,
    bins = range(0, 1, 50),
)
text!(
    ax_ATP_AUC_hist,
    0.0,
    0.08;
    text = "CV=$(round(std(ATP_prod_AUC_data.all_params) / mean(ATP_prod_AUC_data.all_params), sigdigits=2))",
    fontsize = 8,
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
    xticks = (1:length(unique(barplot_df.Parameters)), replace.(unique(barplot_df.Parameters), "HK1 PFKP" => "HK1 & PFKP", "Km Vmax" => "Km & Vmax")),
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
vspan!(ax_Sens_Ind, 1.5, 2.5; ymin = 0.0, ymax = 0.55, color = (:red, 0.1))


labels = ["S1", "ST"]
elements = [PolyElement(polycolor = colors[i]) for i = 1:length(labels)]
axislegend(
    ax_Sens_Ind,
    elements,
    labels,
    # [hist_all_params, hist_fixed_HK_Km_Vmax],
    # ["All glycolysis_params varied", "All glycolysis_params varied\nexcept HK1 Km, Vmax"],
    position = :lt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
    patchsize = (10, 5),
)

colgap!(fig.layout, 7.5)
rowgap!(fig.layout, 5)

label_a = fig[1, 1, TopLeft()] = Label(fig, "A", fontsize = 12, halign = :right, padding = (0, 15, 5, 0))
label_b = fig[1, 2, TopLeft()] = Label(fig, "B", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_c = fig[1, 3, TopLeft()] = Label(fig, "C", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_d = fig[2, 1, TopLeft()] = Label(fig, "D", fontsize = 12, halign = :right, padding = (0, 15, 5, 0))
label_e = fig[2, 2, TopLeft()] = Label(fig, "E", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_f = fig[2, 3, TopLeft()] = Label(fig, "F", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))

fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_FigS3_global_sensitivity_analysis.png", fig, px_per_unit = 4)
