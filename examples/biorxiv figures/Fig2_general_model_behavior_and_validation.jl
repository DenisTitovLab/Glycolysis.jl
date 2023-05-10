using Glycolysis
using DifferentialEquations
using CairoMakie, DataFrames, Dates, Printf, CSV, XLSX, Statistics

##
# Precalculate output of complete model


tspan = (0.0, 240.0 / 4)
Initial_ATPase_Vmax_frac = 0.03
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
        Initial_ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax
end

PresetTime_cb1 = PresetTimeCallback(ATPase_change_time1, affect1!)
PresetTime_cb2 = PresetTimeCallback(ATPase_change_time2, affect2!)
PresetTime_cb3 = PresetTimeCallback(ATPase_change_time3, affect3!)
cb_set = CallbackSet(PresetTime_cb1, PresetTime_cb2, PresetTime_cb3)

init_cond_prob = ODEProblem(glycolysis_ODEs, glycolysis_init_conc, (0, 1e8), glycolysis_params)
init_cond_sol = solve(init_cond_prob, Rodas5P(), abstol = 1e-12, reltol = 1e-5, save_everystep = false)
new_init_cond = init_cond_sol.u[end]
prob = ODEProblem(glycolysis_ODEs, new_init_cond, tspan, glycolysis_params, callback = cb_set)
sol = solve(
    prob,
    Rodas5P(),
    abstol = 1e-15,
    reltol = 1e-5,
    saveat = [k for k = tspan[1]:((tspan[2]-tspan[1])/10_000):tspan[2]],
)

timepoints = sol.t
ATPprod =
    [Glycolysis.conc_to_rates(conc, glycolysis_params).ATPprod for conc in sol.u] / (2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax)
ATPenergy = [Glycolysis.conc_to_disequilibrium_ratios(conc, glycolysis_params).Q_Keq_ATPase for conc in sol.u]
ATPase = Initial_ATPase_Vmax_frac * ones(length(sol))
ATPase[sol.t.>ATPase_change_time1.&&sol.t.<ATPase_change_time2] .= High_ATPase_Vmax_frac
ATPase[sol.t.>ATPase_change_time2.&&sol.t.<ATPase_change_time3] .= Low_ATPase_Vmax_frac
ATP = [conc.ATP for conc in sol.u]



##
# Plot the results
size_inches = (6.5, 5)
size_pt = 72 .* size_inches
set_theme!(Theme(fontsize = 6, Axis = (
    xticksize = 1,
    yticksize = 1,
    # xticklabelsize = 6,
    # yticklabelsize = 6,
    yticklabelpad = 1,
    ylabelpadding = 3,
)))
fig = Figure(resolution = size_pt)

# Plot ATP production
ax_ATP_prod = Axis(
    fig[1, 1:2],
    # limits = (nothing, (0.9, 1.1)),
    limits = (nothing, (0.0, 1.5 * maximum(ATPase))),
    xlabel = "Time, min",
    ylabel = rich("ATP production rate, relative to glycolysis V", subscript("max")),
    yticklabelcolor = Makie.wong_colors()[3],
    ylabelcolor = Makie.wong_colors()[3],
    yticklabelfont = :bold,
    ylabelfont = :bold,
    title = "Matching ATP supply and demand",
)
ax_ATPase = Axis(
    fig[1, 1:2],
    limits = (nothing, (0.0, 1.5 * maximum(ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)
ATP_prod_line = lines!(ax_ATP_prod, timepoints, ATPprod, color = Makie.wong_colors()[3])
ATPase_line = lines!(ax_ATPase, timepoints, ATPase, linestyle = :dot, color = :Black, linewidth = 1)
axislegend(
    ax_ATP_prod,
    [ATP_prod_line, ATPase_line],
    ["ATP production", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, -8, -8),
    patchsize = (7.5, 7.5),
)

# Plot [ATP]
ax_ATP_conc = Axis(
    fig[1, 3:4],
    limits = (nothing, (0, 11.5e-3)),
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
    limits = (nothing, (0.0, 1.5 * maximum(ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
    width = 85,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)
ATP_line = lines!(ax_ATP_conc, timepoints, ATP, color = Makie.wong_colors()[1])
ATPase_line = lines!(ax_ATPase, timepoints, ATPase, linestyle = :dot, color = :Black, linewidth = 1)
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
    padding = (-5, -5, -8, -8),
    patchsize = (7.5, 7.5),
)

# Plot energy of ATPase reaction
ax_ATP_energy = Axis(
    fig[1, 5:6],
    limits = (nothing, (0, 30)),
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
    limits = (nothing, (0.0, 1.5 * maximum(ATPase))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)
ATPase_energy_line = lines!(ax_ATP_energy, timepoints, -log.(ATPenergy), color = Makie.wong_colors()[4])
ATPase_line = lines!(ax_ATPase, timepoints, ATPase, linestyle = :dot, color = :Black, linewidth = 1)
axislegend(
    ax_ATP_energy,
    [ATPase_energy_line, ATPase_line],
    ["ATPase energy", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, -8, -8),
    patchsize = (7.5, 7.5),
)

#Plot [Metabolite] of model vs data
Model_Result_bootstrap = CSV.read(
    "Results/042023_SteadyState_Metabolites_at_ATPase_range_10000bootstraps_w_CI.csv",
    DataFrame,
)
Model_Result_bootstrap = Model_Result_bootstrap[:, Not(r"AMP")]

Experimental_Data = DataFrame(XLSX.readtable("Data/Data S1. Levels of enzymes, metabolites and isotope tracing.xlsx", "Metabolite concentrations"; infer_eltypes=true))

ax1 = Axis(
    fig[2, 1:3],
    limits = (nothing, (1e-9, 1e-1)),
    yscale = log10,
    xticklabelrotation = pi / 2,
    ylabel = "[Metabolite], M",
)
tick_labels = []
column_names =
    replace.(
        names(
            Model_Result_bootstrap[
                :,
                Cols(Between(:Glucose_median, :Lactate_median), Between(:ATP_median, :NADH_median)),
            ],
        ),
        "_median" => "",
    )
for (i, name) in enumerate(column_names)
    points = Model_Result_bootstrap[:, name*"_median"]
    points_qlow = Model_Result_bootstrap[:, name*"_qlow"]
    points_qhigh = Model_Result_bootstrap[:, name*"_qhigh"]
    n_points = nrow(Model_Result_bootstrap)
    positions = i .+ (collect(1:n_points) .- (n_points ÷ 2)) ./ (1.5 * n_points)
    experimental_data = collect(skipmissing(Experimental_Data[:, name]))
    # CairoMakie.boxplot!(
    #     i * ones(length(experimental_data)),
    #     experimental_data,
    #     color = Makie.wong_colors(0.6)[i % 7 == 0 ? 7 : i % 7],
    #     show_outliers = false,
    # )
    band!(positions, points_qlow, points_qhigh, color = Makie.wong_colors(0.3)[(i+3) % 7 == 0 ? 7 : (i+3) % 7])
    # CairoMakie.scatter!(ax1, positions, points, markersize = 5)
    CairoMakie.lines!(ax1, positions, points, color = Makie.wong_colors(1)[(i+3) % 7 == 0 ? 7 : (i+3) % 7])
    CairoMakie.scatter!(
        i * ones(length(experimental_data)) .- 0.25 .+ 0.5 .* rand(length(experimental_data)),
        experimental_data,
        markersize = 3,
        strokewidth = 0.1,
        color = (:white, 0.1),
        # color = Makie.wong_colors(0.6)[i % 7 == 0 ? 7 : i % 7],
        # show_outliers = false,
    )
    push!(tick_labels, name)
end
ax1.xticks = (1:length(column_names), tick_labels)

#Draw an inset
high_bound = 10e-3
low_bound = 1e-3
left_bound = 12.6
right_bound = 13.4
inset_box = Axis(
    fig[2, 1:3];
    halign = :center,
    valign = :bottom,
    width = Relative(0.2),
    height = Relative(0.4),
    alignmode = Mixed(left = 8, right = -20, bottom = 4, top = 2),
    xlabel = rich("ATPase,% of pathway V", subscript("max")),
    ylabel = "[ATP], M",
    xscale = log10,
    yscale = log10,
    limits = ((0.01, 0.225), (low_bound, high_bound)),
    xticks = ([0.01, 0.03, 0.1], string.(Int.(100 .* [0.01, 0.03, 0.1]))),
    # yticks = ([0.01, 0.1, 0.2], string.(Int.(100 .* [0.01, 0.1, 0.2]))),
    ygridvisible = false,
    xgridvisible = false,
    xlabelpadding = 0,
    xticklabelpad = 0,
    ylabelpadding = 0,
    yticklabelpad = 0,
    xlabelsize=5,
    ylabelsize=5,
    xticklabelsize=5,
    yticklabelsize=5,
    spinewidth=0.5
)
translate!(inset_box.blockscene, 0, 0, 1000)
lines!(ax1, [left_bound, left_bound, right_bound, right_bound, left_bound], [low_bound, high_bound, high_bound, low_bound, low_bound], color = :grey, linestyle = :dot)
lines!(ax1, [left_bound, 9.9], [low_bound, 1e-6], color = :grey, linestyle = :dot)
lines!(ax1, [right_bound, 13], [low_bound, 1e-6], color = :grey, linestyle = :dot)
CairoMakie.lines!(inset_box, Model_Result_bootstrap.ATPase_Vmax_frac, Model_Result_bootstrap.ATP_median, color = Makie.wong_colors(1)[3])
band!(Model_Result_bootstrap.ATPase_Vmax_frac, Model_Result_bootstrap[:, "ATP_qlow"], Model_Result_bootstrap[:, "ATP_qhigh"], color = Makie.wong_colors(0.3)[3])
# hideydecorations!(inset_box)


# Plot 13C-tracing of model vs data
Tracing_Panel = fig[2, 4:6] = GridLayout()

#Load Lactate data
Lactate_Model_results = CSV.read(
    "Results/042023_13C_lactate_labeling_w_CI_0.15.csv",
    DataFrame,
);
Lactate_Experimental_data = DataFrame(XLSX.readtable("Data/Data S1. Levels of enzymes, metabolites and isotope tracing.xlsx", "13C Lactate tracing"; infer_eltypes=true))

#Load Glucose data
Glucose_Model_results = CSV.read(
    "Results/042023_13C_glucose_labeling_w_CI_0.15.csv",
    DataFrame,
);
Glucose_Experimental_data = DataFrame(XLSX.readtable("Data/Data S1. Levels of enzymes, metabolites and isotope tracing.xlsx", "13C Glucose tracing"; infer_eltypes=true))


Model_results = Glucose_Model_results
Experimental_data = Glucose_Experimental_data
Experimental_data = Experimental_data[Experimental_data.Cell.=="C2C12", :]
markersize = 5

ax_13C_Glucose = Axis(
    Tracing_Panel[1, 1],
    limits = ((-2, 32), (-0.05, 1.15)),
    xlabel = "Time, min",
    ylabel = "Fraction ¹³C labeling",
    title = "[U-¹³C]Glucose",
)
column_names = replace.(names(Model_results[:, Between(:Glucose_mean, :Lactate_mean)]), "_mean" => "")
column_names = ["Glucose", "G6P", "F16BP", "PEP", "Pyruvate", "Lactate"]
for (i, name) in enumerate(column_names)
    time = Model_results.time
    points = Model_results[:, name*"_mean"]
    points_qlow = Model_results[:, name*"_qlow"]
    points_qhigh = Model_results[:, name*"_qhigh"]
    n_points = length(points)
    lines!(time, points, label = "$name")
    band!(time, points_qlow, points_qhigh, color = Makie.wong_colors(0.1)[i % 7 == 0 ? 7 : i % 7])
    experimental_points =
        [0; [mean(skipmissing(Experimental_data[timepoint, Regex("$name")])) for timepoint = 1:4]]
    experimental_time = [0; unique(Experimental_data[:, "time (min)"])]
    experimental_error =
        [0; [2 * std(skipmissing(Experimental_data[timepoint, Regex("$name")])) for timepoint = 1:4]]
    CairoMakie.scatterlines!(
        ax_13C_Glucose,
        experimental_time,
        experimental_points,
        markersize = markersize,
        markercolor = Makie.wong_colors(1)[i % 7 == 0 ? 7 : i % 7],
        linestyle = :dash,
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
Experimental_data = Experimental_data[Experimental_data.Cell.=="HeLa", :]

ax_13C_Lactate = Axis(
    Tracing_Panel[1, 2],
    limits = ((-2, 32), (-0.1, 1.1)),
    xlabel = "Time, min",
    ylabel = "Fraction ¹³C labeling",
    title = "[U-¹³C]Lactate",
)
column_names = replace.(names(Model_results[:, Between(:Glucose_mean, :Lactate_mean)]), "_mean" => "")
column_names = ["Glucose", "G6P", "F16BP", "PEP", "Pyruvate", "Lactate"]
for (i, name) in enumerate(column_names)
    time = Model_results.time
    points = Model_results[:, name*"_mean"]
    points_qlow = Model_results[:, name*"_qlow"]
    points_qhigh = Model_results[:, name*"_qhigh"]
    n_points = length(points)
    lines!(time, points, label = "$name")
    band!(time, points_qlow, points_qhigh, color = Makie.wong_colors(0.1)[i % 7 == 0 ? 7 : i % 7])
    experimental_points =
        [0; [mean(skipmissing(Experimental_data[timepoint, Regex("$name")])) for timepoint = 1:4]]
    experimental_time = [0; unique(Experimental_data[:, "time (min)"])]
    experimental_error =
        [0; [2 * std(skipmissing(Experimental_data[timepoint, Regex("$name")])) for timepoint = 1:4]]
    CairoMakie.scatterlines!(
        ax_13C_Lactate,
        experimental_time,
        experimental_points,
        markersize = markersize,
        markercolor = Makie.wong_colors(1)[i % 7 == 0 ? 7 : i % 7],
        linestyle = :dash,
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

Model_elem = [LineElement(color = :black, linestyle = nothing)]
Data_elem = [
    LineElement(color = :black, linestyle = :dash),
    MarkerElement(color = :black, marker = '●', markersize = markersize, strokecolor = :black),
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
    [["Model", "Data"], column_names],
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
colgap!(fig.layout, 10)
rowgap!(fig.layout, 5)
resize_to_layout!(fig)

label_a = fig[1, 1, TopLeft()] = Label(fig, "A", fontsize = 12, halign = :right, padding = (0, 15, 10, 0))
label_b = fig[1, 3, TopLeft()] = Label(fig, "B", fontsize = 12, halign = :right, padding = (0, 5, 10, 0))
label_c = fig[1, 5, TopLeft()] = Label(fig, "C", fontsize = 12, halign = :right, padding = (0, 5, 10, 0))
label_d = fig[2, 1, TopLeft()] = Label(fig, "D", fontsize = 12, halign = :right, padding = (0, 15, 0, 0))
label_e = fig[2, 4, TopLeft()] = Label(fig, "E", fontsize = 12, halign = :right, padding = (0, 5, 0, 0))

fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_Fig2_model_behavior_and_validation.png", fig, px_per_unit = 4)
