using Glycolysis
using DifferentialEquations, BenchmarkTools
using CairoMakie, DataFrames, Dates, Printf, CSV, XLSX, Statistics, Measurements
using SwarmMakie
# CairoMakie.activate!(type = "svg")

# M/min * (µl/mg) = µmol/min/mg
0.301 *
Glycolysis.glycolysis_params_w_uncertainty.HK1_Conc *
Glycolysis.glycolysis_params_w_uncertainty.HK1_Vmax *
(Glycolysis.cell_volume_correction / Glycolysis.cell_protein_density)

##
# Plot the results
size_inches = (4, 3)
size_pt = 72 .* size_inches
set_theme!(
    Theme(
        fontsize = 6,
        Axis = (
            xticksize = 1,
            yticksize = 1,
            # xticklabelsize = 6,
            # yticklabelsize = 6,
            yticklabelpad = 1,
            ylabelpadding = 3,
        ),
    ),
)
fig = Figure(size = size_pt)

#Plot [Metabolite] of model vs data
Model_Result_bootstrap = CSV.read(
    # "Results/071324_Glycolysis_Processed_Total_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25.csv",
    "Results/092424_Glycolysis_Processed_Total_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25.csv",
    DataFrame,
)
# Model_Result_bootstrap = CSV.read(
#     "Results/071324_Glycolysis_Processed_Free_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25.csv",
#     DataFrame)

Model_Result_bootstrap_no_allostery = CSV.read(
    # "Results/081524_Glycolysis_Processed_Total_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25_no_allo_reg.csv",
    "Results/092424_Glycolysis_Processed_Total_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25_no_allostery.csv",
    DataFrame,
)
# Model_Result_bootstrap = CSV.read(
#     "Results/081524_Glycolysis_Processed_Free_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25_no_allo_reg.csv",
#     DataFrame)


Model_Result_bootstrap =
    Model_Result_bootstrap[:, Not(r"AMP", r"NTP", r"NDP", r"Phosphocreatine", r"Creatine")]
Model_Result_bootstrap_no_allostery = Model_Result_bootstrap_no_allostery[
    :,
    Not(r"AMP", r"NTP", r"NDP", r"Phosphocreatine", r"Creatine"),
]

Experimental_Data = DataFrame(
    XLSX.readtable(
        "Data/Supplementary File 1. Levels of enzymes, metabolites and isotope tracing.xlsx",
        "Metabolite concentrations";
        infer_eltypes = true,
    ),
)

ax1 = Axis(
    fig[1, 1],
    limits = (nothing, (2e-9, 2e-1)),
    yscale = log10,
    xticklabelrotation = pi / 2,
    ylabel = "Cytosolic [Metabolite], M",
    # include log ticks
    yticks = LogTicks(LinearTicks(8)),
)
tick_labels = []
column_names =
    replace.(
        names(
            Model_Result_bootstrap[
                :,
                Cols(
                    Between(:Glucose_median, :Lactate_median),
                    Between(:ATP_median, :NADH_median),
                ),
            ],
        ),
        "_median" => "",
    )
column_names = filter(x -> !occursin("BPG", x), column_names)
color_index_nt = NamedTuple()

for (i, name) in enumerate(column_names)
    points = Model_Result_bootstrap[:, name*"_median"]
    points_qlow = Model_Result_bootstrap[:, name*"_qlow"]
    points_qhigh = Model_Result_bootstrap[:, name*"_qhigh"]
    n_points = nrow(Model_Result_bootstrap)
    positions = i .+ (collect(1:n_points) .- (n_points ÷ 2)) ./ (1.5 * n_points)
    experimental_data =
        collect(skipmissing(Experimental_Data[:, name])) ./
        Glycolysis.cell_volume_correction
    global color_index_nt = merge(color_index_nt, (; Symbol(name) => i % 7))
    band!(
        positions,
        points_qlow,
        points_qhigh,
        color = Makie.wong_colors(0.3)[(i) % 7 == 0 ? 7 : (i) % 7],
    )
    # CairoMakie.scatter!(ax1, positions, points, markersize = 5)
    CairoMakie.lines!(
        ax1,
        positions,
        points,
        color = Makie.wong_colors(1)[(i) % 7 == 0 ? 7 : (i) % 7],
    )
    beeswarm!(
        ax1,
        i * ones(length(experimental_data)),
        experimental_data,
        algorithm = QuasirandomJitter(; jitter_width = 1.0),
        gutter = 1.0,
        markersize = 3,
        strokewidth = 0.1,
        color = (:white, 0.1),
    )
    push!(tick_labels, name)
end

ax1.xticks = (1:length(column_names), tick_labels)

for (i, name) in enumerate(column_names)
    points = Model_Result_bootstrap_no_allostery[:, name*"_median"]
    points_qlow = Model_Result_bootstrap_no_allostery[:, name*"_qlow"]
    points_qhigh = Model_Result_bootstrap_no_allostery[:, name*"_qhigh"]
    n_points = nrow(Model_Result_bootstrap_no_allostery)
    positions = i .+ (collect(1:n_points) .- (n_points ÷ 2)) ./ (1.5 * n_points)
    experimental_data =
        collect(skipmissing(Experimental_Data[:, name])) ./
        Glycolysis.cell_volume_correction
    global color_index_nt = merge(color_index_nt, (; Symbol(name) => i % 7))
    band!(positions, points_qlow, points_qhigh, color = (:Grey, 0.1))
    # CairoMakie.scatter!(ax1, positions, points, markersize = 5)
    CairoMakie.lines!(
        ax1,
        positions,
        points,
        # color = Makie.wong_colors(1)[(i) % 7 == 0 ? 7 : (i) % 7],
        color = (:Grey, 0.5),
        # linestyle = :dot,
    )
    # push!(tick_labels, name)
end

# hideydecorations!(inset_box)
axislegend(
    ax1,
    [
        [
            # PolyElement(color = Makie.wong_colors(0.3)[1]),
            LineElement(
                points = Point2f.(range(0, 1, 7), 0.5 .* ones(7)),
                color = Makie.wong_colors(0.3)[1:7],
                linewidth = 5,
            ),
            LineElement(
                points = Point2f.(range(0, 1, 7), 0.5 .* ones(7)),
                color = Makie.wong_colors()[1:7],
            ),
        ],
        [PolyElement(color = (:Grey, 0.1)), LineElement(color = (:Grey, 0.5))],
        [
            MarkerElement(
                marker = :circle,
                markersize = 4,
                strokewidth = 0.15,
                color = (:white, 0.1),
            ),
        ],
    ],
    ["Model ± 95%CI", "Model wo allostery ± 95%CI", "Data"],
    position = :lb,
    patchsize = (10.0f0, 5.0f0),
    groupgap = 8,
    padding = (2.0f0, 2.0f0, 2.0f0, 2.0f0),
    patchlabelgap = 3,
    titlegap = 2,
    # framevisible = false,
    framewidth = 0.1,
    backgroundcolor = (:white, 0.75),
)

# Final Figure edits
colgap!(fig.layout, 5)
rowgap!(fig.layout, 5)
resize_to_layout!(fig)

fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_FigS2_effect_of_allostery_on_metab_levels.png", fig, px_per_unit = 4)
