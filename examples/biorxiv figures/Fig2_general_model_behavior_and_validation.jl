using Glycolysis
using OrdinaryDiffEq, DiffEqCallbacks, BenchmarkTools
using CairoMakie, DataFrames, DataFramesMeta, Dates, Printf, CSV, XLSX, Statistics, Measurements, StatsBase
using SwarmMakie

##
# Precalculate output of complete model

tspan = (0.0, 240.0 / 4)
Initial_ATPase_Vmax_frac = 0.05
High_ATPase_Vmax_frac = Initial_ATPase_Vmax_frac * 2
Low_ATPase_Vmax_frac = Initial_ATPase_Vmax_frac / 2
ATPase_change_time1 = 60 / 4
ATPase_change_time2 = 120 / 4
ATPase_change_time3 = 180 / 4
glycolysis_params.ATPase_Vmax =
    Initial_ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax

function affect1!(integrator)
    integrator.p.ATPase_Vmax =
        High_ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax
end
function affect2!(integrator)
    integrator.p.ATPase_Vmax =
        Low_ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax
end
function affect3!(integrator)
    integrator.p.ATPase_Vmax =
        Initial_ATPase_Vmax_frac *
        2 *
        glycolysis_params.HK1_Conc *
        glycolysis_params.HK1_Vmax
end

PresetTime_cb1 = PresetTimeCallback(ATPase_change_time1, affect1!)
PresetTime_cb2 = PresetTimeCallback(ATPase_change_time2, affect2!)
PresetTime_cb3 = PresetTimeCallback(ATPase_change_time3, affect3!)
cb_set = CallbackSet(PresetTime_cb1, PresetTime_cb2, PresetTime_cb3)

init_cond_prob =
    ODEProblem(glycolysis_ODEs, glycolysis_init_conc, (0, 1e8), glycolysis_params)
init_cond_sol =
    solve(init_cond_prob, Rodas5P(), abstol = 1e-15, reltol = 1e-8, save_everystep = false)
new_init_cond = init_cond_sol.u[end]
prob =
    ODEProblem(glycolysis_ODEs, new_init_cond, tspan, glycolysis_params, callback = cb_set)
sol = solve(
    prob,
    Rodas5P(),
    abstol = 1e-15,
    reltol = 1e-8,
    saveat = [k for k = tspan[1]:((tspan[2]-tspan[1])/10_000):tspan[2]],
)

timepoints = sol.t
ATPprod =
    [Glycolysis.conc_to_rates(conc, glycolysis_params).ATPprod for conc in sol.u] / (2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax)
ATPenergy = [
    Glycolysis.conc_to_disequilibrium_ratios(conc, glycolysis_params).Q_Keq_ATPase for
    conc in sol.u
]
ATPase = Initial_ATPase_Vmax_frac * ones(length(sol))
ATPase[sol.t.>ATPase_change_time1.&&sol.t.<ATPase_change_time2] .= High_ATPase_Vmax_frac
ATPase[sol.t.>ATPase_change_time2.&&sol.t.<ATPase_change_time3] .= Low_ATPase_Vmax_frac
ATP = [conc.ATP for conc in sol.u]

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
            # xticklabelsize = 6,
            # yticklabelsize = 6,
            yticklabelpad = 1,
            ylabelpadding = 3,
        ),
    ),
)
fig = Figure(size = size_pt)

# Plot ATP production
ax_ATP_prod = Axis(
    fig[1, 1:3],
    # limits = (nothing, (0.9, 1.1)),
    limits = (nothing, (0.0, 1.5 * maximum(ATPase))),
    xlabel = "Time, min",
    ylabel = "ATP prod. rate, relative to glycolysis Vmax",
    yticklabelcolor = Makie.wong_colors()[3],
    ylabelcolor = Makie.wong_colors()[3],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    title = "Matching ATP supply\nand demand",
    width = 75,
    yticklabelrotation = pi / 2,
    yticks = LinearTicks(3),
)
ax_ATPase = Axis(
    fig[1, 1:3],
    limits = (nothing, (0.0, 1.5 * maximum(ATPase))),
    ylabel = "ATPase rate, relative to glycolysis Vmax",
    yaxisposition = :right,
    ygridvisible = false,
    width = 75,
    yticklabelrotation = pi / 2,
    yticks = LinearTicks(3),
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)
ATP_prod_line = lines!(ax_ATP_prod, timepoints, ATPprod, color = Makie.wong_colors()[3])
ATPase_line =
    lines!(ax_ATPase, timepoints, ATPase, linestyle = :dot, color = :Black, linewidth = 1)
axislegend(
    ax_ATP_prod,
    [ATP_prod_line, ATPase_line],
    ["ATP prod.", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (0, -4, 0, -4),
    patchsize = (7.5, 7.5),
)

# Plot dynamic [ATP] top inset
adenine_pool_size =
    glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP
Pane_B_ATP = fig[1, 4:6] = GridLayout(2, 1)
ax_ATP_conc = Axis(
    Pane_B_ATP[1, 1],
    # fig[1, 2][1, 1],
    limits = (nothing, (0.998 * adenine_pool_size, 1.00055 * adenine_pool_size)),
    xlabel = "Time, min",
    ylabel = "[ATP], mM",
    ytickformat = ys -> ["$(round(y*1000, sigdigits = 4))" for y in ys],
    yticklabelcolor = Makie.wong_colors()[1],
    ylabelcolor = Makie.wong_colors()[1],
    yticks = [10.29e-3, 10.31e-3],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    title = "Maintaining ATP concentration",
    width = 75,
    yticklabelrotation = pi / 2,
)
hidexdecorations!(ax_ATP_conc)
# hidespines!(ax_ATPase)
ATP_line = lines!(ax_ATP_conc, timepoints, ATP, color = Makie.wong_colors()[1])
lines!(
    ax_ATP_conc,
    [tspan[1], tspan[2]],
    repeat(
        [glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP],
        2,
    ),
    color = :grey,
    linestyle = :dash,
)
text!(
    ax_ATP_conc,
    0,
    (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = :grey,
)
axislegend(
    ax_ATP_conc,
    [ATP_line],
    ["[ATP]"],
    position = (1, -0.1),
    colgap = 5,
    framevisible = false,
    padding = (0, -4, 0, -4),
    patchsize = (7.5, 7.5),
    orientation = :horizontal,
)
# Plot dynamic [ATP] bottom
ax_ATP_conc = Axis(
    Pane_B_ATP[2:3, 1],
    limits = (nothing, (0, 1.3 * adenine_pool_size)),
    xlabel = "Time, min",
    ylabel = "[ATP], mM",
    ytickformat = ys -> ["$(round(y*1000, sigdigits = 3))" for y in ys],
    yticklabelcolor = Makie.wong_colors()[1],
    ylabelcolor = Makie.wong_colors()[1],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    width = 75,
    yticklabelrotation = pi / 2,
    tellwidth = true,
)
ax_ATPase = Axis(
    Pane_B_ATP[2:3, 1],
    limits = (nothing, (0.0, 1.5 * maximum(ATPase))),
    ylabel = "ATPase rate, relative to glycolysis Vmax",
    yaxisposition = :right,
    ygridvisible = false,
    width = 75,
    yticklabelrotation = pi / 2,
    yticks = LinearTicks(3),
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)
ATP_line = lines!(ax_ATP_conc, timepoints, ATP, color = Makie.wong_colors()[1])
ATPase_line =
    lines!(ax_ATPase, timepoints, ATPase, linestyle = :dot, color = :Black, linewidth = 1)
lines!(
    ax_ATP_conc,
    [tspan[1], tspan[2]],
    repeat(
        [glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP],
        2,
    ),
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
    position = :ct,
    colgap = 5,
    patchlabelgap = 2,
    framevisible = false,
    padding = (0, -4, 0, -4),
    patchsize = (5, 5),
    orientation = :horizontal,
)
rowgap!(Pane_B_ATP, 7.5)


# Plot energy of ATPase reaction
ax_ATP_energy = Axis(
    fig[1, 7:9],
    limits = (nothing, (0, 32)),
    xlabel = "Time, min",
    ylabel = rich("Energy of ATPase reaction, k", subscript("B"), "T"),
    # yscale = log10,
    ytickformat = ys -> ["$(Int(round(y)))" for y in ys],
    yticklabelcolor = Makie.wong_colors()[4],
    ylabelcolor = Makie.wong_colors()[4],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    title = "Maintaining ATP energy",
    width = 75,
    yticklabelrotation = pi / 2,
)
ax_ATP_energy_ATPase = Axis(
    fig[1, 7:9],
    limits = (nothing, (0.0, 1.5 * maximum(ATPase))),
    ylabel = "ATPase rate, relative to glycolysis Vmax",
    yaxisposition = :right,
    ygridvisible = false,
    width = 75,
    yticklabelrotation = pi / 2,
    yticks = LinearTicks(3),
)
hidespines!(ax_ATP_energy_ATPase)
hidexdecorations!(ax_ATP_energy_ATPase)
ATPase_energy_line =
    lines!(ax_ATP_energy, timepoints, -log.(ATPenergy), color = Makie.wong_colors()[4])
ATPase_line = lines!(
    ax_ATP_energy_ATPase,
    timepoints,
    ATPase,
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
    padding = (0, -4, 0, -4),
    patchsize = (7.5, 7.5),
)

#Plot [Metabolite] of model vs data
Model_Result_bootstrap = CSV.read(
    "Results/092424_Glycolysis_Processed_Total_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25.csv",
    # "Results/092824_Glycolysis_Processed_Total_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25_10uMF26BP.csv",
    DataFrame,
)
# Model_Result_bootstrap = CSV.read(
#     "Results/092424_Glycolysis_Processed_Free_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25.csv",
#     DataFrame)

# Model_Result_bootstrap = CSV.read(
#     "Results/081524_Glycolysis_Processed_Total_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25_no_allo_reg.csv",
#     DataFrame,
# )
# Model_Result_bootstrap = CSV.read(
#     "Results/081524_Glycolysis_Processed_Free_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25_no_allo_reg.csv",
#     DataFrame)


Model_Result_bootstrap =
    Model_Result_bootstrap[:, Not(r"AMP", r"NTP", r"NDP", r"Phosphocreatine", r"Creatine")]

Experimental_Data = DataFrame(
    XLSX.readtable(
        "Data/Supplementary File 1. Levels of enzymes, metabolites and isotope tracing.xlsx",
        "Metabolite concentrations";
        infer_eltypes = true,
    ),
)

ax_glyc_rates = Axis(
    fig[1, 10],
    limits = (nothing, (0.2e-4, 1e-1)),
    ylabel = "Glycolysis rate, µmol/min per mg cell protein",
    # yticks = LinearTicks(3),
    yscale = log10,
    # yticks = LogTicks(LinearTicks(3)),
    yticks = [1e-4, 1e-3, 1e-2],
    ygridvisible = false,
    width = 10,
    yticklabelrotation = pi / 2,
)
max_model_rate_w_uncertainty =
    Glycolysis.glycolysis_params_w_uncertainty.HK1_Conc *
    Glycolysis.glycolysis_params_w_uncertainty.HK1_Vmax *
    (Glycolysis.cell_volume_correction / Glycolysis.cell_protein_density)

band!(
    ax_glyc_rates,
    [0.6, 1.4],
    0.301 .* [
        Measurements.value(max_model_rate_w_uncertainty) +
        Measurements.uncertainty(max_model_rate_w_uncertainty),
    ],
    0.009 .* [
        Measurements.value(max_model_rate_w_uncertainty) -
        Measurements.uncertainty(max_model_rate_w_uncertainty),
    ],
    color = Makie.wong_colors(0.3)[1],
)
hidexdecorations!(ax_glyc_rates)

glyc_rate_data = [
    0.058583333,
    0.031083333,
    0.028833333,
    0.038083333,
    0.027583333,
    0.03275,
    0.029333333,
    0.003872222,
    0.004333333,
    0.003916111,
    0.0103,
    0.006188889,
    0.006105556,
    0.011277778,
    0.012311111,
    0.006851667,
    0.000256111,
    0.004001852,
    0.010777778,
    0.0023,
    0.003705556,
    0.003716667,
    0.000222083,
    0.000238681,
    0.000236319,
    0.012416667,
    0.002062745,
    0.002413333,
]
sort(glyc_rate_data)
median(glyc_rate_data) / Measurements.value(max_model_rate_w_uncertainty)
beeswarm!(
    ax_glyc_rates,
    1.0 .* ones(length(glyc_rate_data)),
    glyc_rate_data,
    algorithm = QuasirandomJitter(; jitter_width = 1.0),
    gutter = 1.0,
    markersize = 3,
    strokewidth = 0.1,
    color = (:white, 0.1),
)
Legend(
    fig[1, 10, Right()],
    [
        [PolyElement(color = (Makie.wong_colors()[1], 0.3))],
        [
            MarkerElement(
                marker = :circle,
                markersize = 4,
                strokewidth = 0.15,
                color = (:white, 0.1),
            ),
        ],
    ],
    ["Model", "Data"],
    # position = :cb,
    # patchsize = (5.0f0, 5.0f0),
    # groupgap = 1,
    # padding = (1.0f0, 1.0f0, 1.0f0, 1.0f0),
    # patchlabelgap = 1,
    # titlegap = 2,
    # framevisible = false,
    # framewidth = 0.1,
    # backgroundcolor = (:white, 0.75),
    patchsize = (5.0f0, 5.0f0),
    groupgap = 8,
    padding = (2.0f0, 2.0f0, 2.0f0, 2.0f0),
    patchlabelgap = 3,
    titlegap = 2,
    framevisible = false,
    alignmode = Mixed(right = 10),
)








#Plot data vs model metabolite concentrations
ax_concs = Axis(
    fig[2, 1:5],
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
    # CairoMakie.scatter!(ax_concs, positions, points, markersize = 5)
    CairoMakie.lines!(
        ax_concs,
        positions,
        points,
        color = Makie.wong_colors(1)[(i) % 7 == 0 ? 7 : (i) % 7],
    )
    beeswarm!(
        ax_concs,
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

ax_concs.xticks = (1:length(column_names), tick_labels)

#Draw an inset
high_bound = 15e-3
low_bound = 1e-3
left_bound = 11.6
right_bound = 12.4
inset_box = Axis(
    fig[2, 1:5];
    halign = :center,
    valign = :bottom,
    width = Relative(0.2),
    height = Relative(0.4),
    alignmode = Mixed(left = 4, right = -16, bottom = 4, top = 2),
    xlabel = rich("ATPase,% of pathway V", subscript("max")),
    ylabel = "[ATP], M",
    xscale = log10,
    yscale = log10,
    limits = ((0.02 / 1.5, 0.2 * 1.5), (low_bound, high_bound)),
    xticks = ([0.02, 0.06, 0.2], string.(Int.(100 .* [0.02, 0.06, 0.2]))),
    # yticks = ([0.01, 0.1, 0.2], string.(Int.(100 .* [0.01, 0.1, 0.2]))),
    yticks = LogTicks(LinearTicks(1)),
    ygridvisible = false,
    xgridvisible = false,
    xlabelpadding = 0,
    xticklabelpad = 0,
    ylabelpadding = 0,
    yticklabelpad = 0,
    xlabelsize = 5,
    ylabelsize = 5,
    xticklabelsize = 5,
    yticklabelsize = 5,
    spinewidth = 0.5,
)
translate!(inset_box.blockscene, 0, 0, 1000)
lines!(
    ax_concs,
    [left_bound, left_bound, right_bound, right_bound, left_bound],
    [low_bound, high_bound, high_bound, low_bound, low_bound],
    color = :grey,
    linestyle = :dot,
)
lines!(ax_concs, [left_bound, 9.15], [low_bound, 1e-6], color = :grey, linestyle = :dot)
lines!(ax_concs, [right_bound, 12.25], [low_bound, 1e-6], color = :grey, linestyle = :dot)
CairoMakie.lines!(
    inset_box,
    Model_Result_bootstrap.ATPase_Vmax_frac,
    Model_Result_bootstrap.ATP_median,
    color = Makie.wong_colors(1)[color_index_nt.ATP],
)
band!(
    Model_Result_bootstrap.ATPase_Vmax_frac,
    Model_Result_bootstrap[:, "ATP_qlow"],
    Model_Result_bootstrap[:, "ATP_qhigh"],
    color = Makie.wong_colors(0.3)[color_index_nt.ATP],
)
# hideydecorations!(inset_box)
axislegend(
    ax_concs,
    [
        [
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
        [
            MarkerElement(
                marker = :circle,
                markersize = 4,
                strokewidth = 0.15,
                color = (:white, 0.1),
            ),
        ],
    ],
    ["Model±95%CI", "Data"],
    position = :lb,
    patchsize = (7.5f0, 5.0f0),
    groupgap = 8,
    padding = (1.0f0, 1.0f0, 1.0f0, 1.0f0),
    patchlabelgap = 1,
    titlegap = 2,
    # framevisible = false,
    framewidth = 0.1,
    backgroundcolor = (:white, 0.75),
)










#Plot disequilibrium ratios of model vs data
#Load data
Model_Result_bootstrap_free = CSV.read(
    "Results/092424_Glycolysis_Free_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25.csv",
    DataFrame,
)
Model_Result_bootstrap_no_allostery_free = CSV.read(
    "Results/092424_Glycolysis_Free_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25_no_allostery.csv",
    DataFrame,
)
Experimental_Data = DataFrame(
    XLSX.readtable(
        "Data/Supplementary File 1. Levels of enzymes, metabolites and isotope tracing.xlsx",
        "Metabolite concentrations";
        infer_eltypes = true,
    ),
)

#Process data to extract Q/Keq ratios for each reaction
Disequilibrium_Ratios = DataFrame()
for row in eachrow(Model_Result_bootstrap_free)
    push!(
        Disequilibrium_Ratios,
        merge(
            (ATPase_Vmax_frac = row.ATPase_Vmax_frac,),
            (ATPprod_ATPase_ratios = row.ATPprod_ATPase_ratios,),
            convert(
                NamedTuple,
                Glycolysis.conc_to_disequilibrium_ratios(row, Glycolysis.glycolysis_params),
            ),
        ),
    )
end
qlow(x) = percentile(x, 2.5)
qhigh(x) = percentile(x, 97.5)
Processed_Disequilibrium_Ratios = @chain Disequilibrium_Ratios begin
    @rsubset!(1.01 >= :ATPprod_ATPase_ratios >= 0.99)
    @by(
        :ATPase_Vmax_frac,
        :count = length(:Q_Keq_GLUT),
        $(propertynames(Disequilibrium_Ratios) .=> median),
        $(propertynames(Disequilibrium_Ratios) .=> qlow),
        $(propertynames(Disequilibrium_Ratios) .=> qhigh),
    )
end

#Process no allostery data to extract Q/Keq ratios for each reaction
Disequilibrium_Ratios_no_allostery = DataFrame()
for row in eachrow(Model_Result_bootstrap_no_allostery_free)
    push!(
        Disequilibrium_Ratios_no_allostery,
        merge(
            (ATPase_Vmax_frac = row.ATPase_Vmax_frac,),
            (ATPprod_ATPase_ratios = row.ATPprod_ATPase_ratios,),
            convert(
                NamedTuple,
                Glycolysis.conc_to_disequilibrium_ratios(row, Glycolysis.glycolysis_params),
            ),
        ),
    )
end
qlow(x) = percentile(x, 2.5)
qhigh(x) = percentile(x, 97.5)
Processed_Disequilibrium_Ratios_no_allostery =
    @chain Disequilibrium_Ratios_no_allostery begin
        @rsubset!(1.01 >= :ATPprod_ATPase_ratios >= 0.99)
        @by(
            :ATPase_Vmax_frac,
            :count = length(:Q_Keq_GLUT),
            $(propertynames(Disequilibrium_Ratios_no_allostery) .=> median),
            $(propertynames(Disequilibrium_Ratios_no_allostery) .=> qlow),
            $(propertynames(Disequilibrium_Ratios_no_allostery) .=> qhigh),
        )
    end

#Process experimental data to extract Q/Keq ratios for each reaction

function sample_non_missing(df::DataFrame)
    sampled_row = NamedTuple()
    for col in names(df)
        non_missing_values = skipmissing(df[!, col])
        sampled_value = sample(collect(non_missing_values), 1)[1]
        sampled_row = merge(sampled_row, (Symbol(col) => sampled_value,))
    end
    return sampled_row
end

Disequilibrium_Ratios_data = DataFrame()
glycolysis_params_modified_keq = deepcopy(Glycolysis.glycolysis_params)
glycolysis_params_modified_keq.ALDO_Keq = 1.3e-4
glycolysis_params_modified_keq.TPI_Keq = 0.11
glycolysis_params_modified_keq.GAPDH_Keq = 2.0
n_bootstrap = 1000
for i = 1:n_bootstrap
    #sample one non missing value from each column of the experimental data
    bootstrap = sample_non_missing(Experimental_Data[:, Between(:Glucose_media, :NADH)])
    push!(
        Disequilibrium_Ratios_data,
        convert(
            NamedTuple,
            Glycolysis.conc_to_disequilibrium_ratios(
                bootstrap,
                glycolysis_params_modified_keq,
            ),
        ),
    )
end
# Filter NaN from DataFrame
filter!(row -> !any(isnan(x) for x in row), Disequilibrium_Ratios_data)

ax_disequil_ratios = Axis(
    fig[2, 6:7],
    # limits = (nothing, (2e-9, 2e-1)),
    yscale = log10,
    xticklabelrotation = pi / 2,
    ylabel = "Disequilibrium Ratio",
    # include log ticks
    yticks = LogTicks(LinearTicks(8)),
)
tick_labels = []
column_names = names(Disequilibrium_Ratios[:, Between(:Q_Keq_HK1, :Q_Keq_ATPase)])
filter!(x -> x != "Q_Keq_MCT", column_names)
for (i, name) in enumerate(column_names)
    points = Processed_Disequilibrium_Ratios[:, name*"_median"]
    points_qlow = Processed_Disequilibrium_Ratios[:, name*"_qlow"]
    points_qhigh = Processed_Disequilibrium_Ratios[:, name*"_qhigh"]
    n_points = nrow(Processed_Disequilibrium_Ratios)
    positions = i .+ (collect(1:n_points) .- (n_points ÷ 2)) ./ (1.5 * n_points)
    CairoMakie.boxplot!(
        ax_disequil_ratios,
        i * ones(nrow(Disequilibrium_Ratios_data)),
        Disequilibrium_Ratios_data[!, name],
        whiskerwidth = 0.5,
        show_outliers = false,
        width = 0.75,
        strokewidth = 1,
        whiskerlinewidth = 1,
        medianlinewidth = 1,
        strokecolor = (:grey, 0.5),
        mediancolor = (:grey, 0.5),
        whiskercolor = (:grey, 0.5),
        color = (:grey, 0.1),
        label = "Experimental Data",
    )
    CairoMakie.lines!(
        ax_disequil_ratios,
        positions,
        points,
        color = Makie.wong_colors(1)[(i) % 7 == 0 ? 7 : (i) % 7],
        label = "Model Predictions ± 95%CI",
    )
    band!(
        positions,
        points_qlow,
        points_qhigh,
        color = Makie.wong_colors(0.3)[(i) % 7 == 0 ? 7 : (i) % 7],
        label = "Model Predictions ± 95%CI",
    )
    push!(tick_labels, replace(name, "Q_Keq_" => ""))
end
ax_disequil_ratios.xticks = (1:length(column_names), tick_labels)
axislegend(
    ax_disequil_ratios,
    [
        [
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
        [
            PolyElement(
                # marker = :circle,
                # markersize = 4,
                strokewidth = 0.15,
                color = (:grey, 0.1),
            ),
            LineElement(
                # points = Point2f.(range(0, 1, 7), 0.5 .* ones(7)),
                color = (:grey, 0.5),
            ),
        ],
    ],
    ["Model±95%CI", "Data"],
    # position = :rb,
    position = (0, 0),
    patchsize = (7.5f0, 5.0f0),
    groupgap = 8,
    padding = (1.0f0, 1.0f0, 1.0f0, 1.0f0),
    patchlabelgap = 1,
    titlegap = 2,
    # framevisible = false,
    framewidth = 0.1,
    backgroundcolor = (:white, 0.75),
)








# Plot 13C-tracing of model vs data
Tracing_Panel = fig[2, 8:10] = GridLayout()

#Load Lactate data
Lactate_Model_results =
    CSV.read("Results/071324_13C_lactate_labeling_w_CI_0.15.csv", DataFrame)

Lactate_Experimental_data = DataFrame(
    XLSX.readtable(
        "Data/Supplementary File 1. Levels of enzymes, metabolites and isotope tracing.xlsx",
        "13C Lactate tracing";
        infer_eltypes = true,
    ),
)

#Load Glucose data
Glucose_Model_results =
    CSV.read("Results/071324_13C_glucose_labeling_w_CI_0.15.csv", DataFrame)

Glucose_Experimental_data = DataFrame(
    XLSX.readtable(
        "Data/Supplementary File 1. Levels of enzymes, metabolites and isotope tracing.xlsx",
        "13C Glucose tracing";
        infer_eltypes = true,
    ),
)

Model_results = Glucose_Model_results
Experimental_data = Glucose_Experimental_data
markersize = 5

ax_13C_Glucose = Axis(
    Tracing_Panel[1, 1],
    limits = ((-2, 32), (-0.05, 1.15)),
    xlabel = "Time, min",
    ylabel = "Fraction ¹³C labeling",
    title = "[U-¹³C]Glucose",
)
column_names =
    replace.(names(Model_results[:, Between(:Glucose_mean, :Lactate_mean)]), "_mean" => "")
column_names = ["Glucose", "G6P", "F16BP", "PEP", "Pyruvate", "Lactate"]
# column_names = ["Lactate"]
for (i, name) in enumerate(column_names)
    time = Model_results.time
    points = Model_results[:, name*"_mean"]
    points_qlow = Model_results[:, name*"_qlow"]
    points_qhigh = Model_results[:, name*"_qhigh"]
    n_points = length(points)
    lines!(time, points, label = "$name")
    band!(
        time,
        points_qlow,
        points_qhigh,
        color = Makie.wong_colors(0.1)[i % 7 == 0 ? 7 : i % 7],
    )
    experimental_time = [0; unique(Experimental_data[:, "time (min)"])]
    experimental_points = [
        0
        [
            mean(
                skipmissing(
                    vec(
                        Matrix(
                            Experimental_data[
                                Experimental_data[:, "time (min)"].==timepoint,
                                Regex("$name"),
                            ],
                        ),
                    ),
                ),
            ) for timepoint in experimental_time[2:end]
        ]
    ]
    experimental_error = [
        0
        [
            std(
                skipmissing(
                    vec(
                        Matrix(
                            Experimental_data[
                                Experimental_data[:, "time (min)"].==timepoint,
                                Regex("$name"),
                            ],
                        ),
                    ),
                ),
            ) for timepoint in experimental_time[2:end]
        ]
    ]
    CairoMakie.scatter!(
        ax_13C_Glucose,
        experimental_time,
        experimental_points,
        markersize = markersize,
        color = Makie.wong_colors(1)[i % 7 == 0 ? 7 : i % 7],
        label = "$name",
    )
    CairoMakie.errorbars!(
        ax_13C_Glucose,
        experimental_time,
        experimental_points,
        experimental_error,
        color = Makie.wong_colors(1)[i % 7 == 0 ? 7 : i % 7],
    )
end

Model_results = Lactate_Model_results
Experimental_data = Lactate_Experimental_data

ax_13C_Lactate = Axis(
    Tracing_Panel[1, 2],
    limits = ((-2, 32), (-0.1, 1.1)),
    xlabel = "Time, min",
    ylabel = "Fraction ¹³C labeling",
    title = "[U-¹³C]Lactate",
)
column_names =
    replace.(names(Model_results[:, Between(:Glucose_mean, :Lactate_mean)]), "_mean" => "")
column_names = ["Glucose", "G6P", "F16BP", "PEP", "Pyruvate", "Lactate"]
# column_names = ["Lactate"]
for (i, name) in enumerate(column_names)
    time = Model_results.time
    points = Model_results[:, name*"_mean"]
    points_qlow = Model_results[:, name*"_qlow"]
    points_qhigh = Model_results[:, name*"_qhigh"]
    n_points = length(points)
    lines!(time, points, label = "$name")
    band!(
        time,
        points_qlow,
        points_qhigh,
        color = Makie.wong_colors(0.1)[i % 7 == 0 ? 7 : i % 7],
    )
    experimental_time = [0; unique(Experimental_data[:, "time (min)"])]
    experimental_points = [
        0
        [
            mean(
                skipmissing(
                    vec(
                        Matrix(
                            Experimental_data[
                                Experimental_data[:, "time (min)"].==timepoint,
                                Regex("$name"),
                            ],
                        ),
                    ),
                ),
            ) for timepoint in experimental_time[2:end]
        ]
    ]
    experimental_error = [
        0
        [
            std(
                skipmissing(
                    vec(
                        Matrix(
                            Experimental_data[
                                Experimental_data[:, "time (min)"].==timepoint,
                                Regex("$name"),
                            ],
                        ),
                    ),
                ),
            ) for timepoint in experimental_time[2:end]
        ]
    ]
    CairoMakie.scatter!(
        ax_13C_Lactate,
        experimental_time,
        experimental_points,
        markersize = markersize,
        color = Makie.wong_colors(1)[i % 7 == 0 ? 7 : i % 7],
        label = "$name",
    )
    CairoMakie.errorbars!(
        ax_13C_Lactate,
        experimental_time,
        experimental_points,
        experimental_error,
        color = Makie.wong_colors(1)[i % 7 == 0 ? 7 : i % 7],
    )
end
Model_elem = [[PolyElement(color = (:black, 0.3)), LineElement()]]
Data_elem = [
    MarkerElement(
        color = :black,
        marker = '●',
        markersize = markersize,
        strokecolor = :black,
    ),
]
Metabolite_elem = [
    [
        LineElement(color = Makie.wong_colors()[i % 7 == 0 ? 7 : i % 7]),
        MarkerElement(
            color = Makie.wong_colors()[i % 7 == 0 ? 7 : i % 7],
            marker = '●',
            markersize = markersize,
            strokecolor = :black,
        ),
    ] for i = 1:length(column_names)
]
Legend(
    Tracing_Panel[1, 2, Right()],
    [[Model_elem, Data_elem], Metabolite_elem],
    [["Modeln95%CI", "Data"], column_names],
    ["Labels", "Metabolites"],
    patchsize = (10.0f0, 5.0f0),
    groupgap = 8,
    padding = (2.0f0, 2.0f0, 2.0f0, 2.0f0),
    patchlabelgap = 3,
    titlegap = 2,
    framevisible = false,
)
linkyaxes!(ax_13C_Glucose, ax_13C_Lactate)
hideydecorations!(ax_13C_Lactate, grid = false)
colgap!(Tracing_Panel, 5)

# Final Figure edits
Pane_B_ATP.alignmode = Mixed(left = -8)
# ax_ATP_energy.alignmode = Mixed(left = -1)
# ax_ATP_energy_ATPase.alignmode = Mixed(left = -1)


colgap!(fig.layout, 7.5)
rowgap!(fig.layout, 5)
resize_to_layout!(fig)




label_a =
    fig[1, 1, TopLeft()] =
        Label(fig, "A", fontsize = 12, halign = :right, padding = (0, 10, 0, 0))
label_b =
    fig[1, 4, TopLeft()] =
        Label(fig, "B", fontsize = 12, halign = :right, padding = (0, 5, 0, 0))
label_c =
    fig[1, 7, TopLeft()] =
        Label(fig, "C", fontsize = 12, halign = :right, padding = (0, 5, 0, 0))
label_d =
    fig[1, 10, TopLeft()] =
        Label(fig, "D", fontsize = 12, halign = :right, padding = (0, 10, 0, 0))
label_e =
    fig[2, 1, TopLeft()] =
        Label(fig, "E", fontsize = 12, halign = :right, padding = (0, 15, 0, 0))
label_f =
    fig[2, 6, TopLeft()] =
        Label(fig, "F", fontsize = 12, halign = :right, padding = (0, 15, 0, 0))
label_g =
    fig[2, 8, TopLeft()] =
        Label(fig, "G", fontsize = 12, halign = :right, padding = (0, 10, 0, 0))

fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_Fig2_model_behavior_and_validation_Lact_media_0_w_disequilibrium.pdf", fig, pt_per_unit = 1)
