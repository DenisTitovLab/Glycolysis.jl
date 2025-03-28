using Glycolysis
using OrdinaryDiffEq, ProgressMeter
using CairoMakie, Dates, Printf, Statistics, StatsBase
using DataFrames, XLSX

# glycolysis_params.AK_Vmax = 0.03
# glycolysis_params.NDPK_Vmax = 0.03
# glycolysis_params.CK_Vmax = 0.0

##
# Precalculate and save outputs for each model
using Distributed
addprocs()

@everywhere using OrdinaryDiffEq, Glycolysis

start = now()

function calculate_ATP_at_range_ATPases(glycolysis_params, glycolysis_init_conc)
    tspan = (0.0, 1e8)
    prob = ODEProblem(glycolysis_ODEs, glycolysis_init_conc, tspan, glycolysis_params)
    n_Vmax_ATPase_values = 1000
    Pathway_Vmax = 2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax
    ATPases = 10 .^ range(log10(0.003), log10(0.45), n_Vmax_ATPase_values) .* Pathway_Vmax
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(
        ensemble_prob,
        Rodas5P(),
        EnsembleDistributed(),
        trajectories = n_Vmax_ATPase_values,
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
    )

    ATP_conc = [sol.u[end].ATP for sol in sim if sol.retcode == ReturnCode.Success]
    ATP_energy =
        -log.([
            Glycolysis.conc_to_disequilibrium_ratios(sol.u[end], sol.prob.p).Q_Keq_ATPase
            for sol in sim if sol.retcode == ReturnCode.Success
        ])
    Converged_sol_ATPases =
        [ATPases[i] for (i, sol) in enumerate(sim) if sol.retcode == ReturnCode.Success]

    return (
        range_ATPases = Converged_sol_ATPases ./ Pathway_Vmax,
        ATP_conc = ATP_conc,
        ATP_energy = ATP_energy,
    )
end

Complete_Model_Simulation_Data =
    calculate_ATP_at_range_ATPases(glycolysis_params, glycolysis_init_conc)

#Precalculate results at different enzyme levels
glycolysis_params_list = []
model_names_list = []
for multiplier in [10, 1, 0.1, 0.02]
    local glycolysis_params_copy = deepcopy(glycolysis_params)
    #glycolysis_params_copy.ALDO_Conc = multiplier * glycolysis_params.ALDO_Conc
    glycolysis_params_copy.GAPDH_Conc = multiplier * glycolysis_params.GAPDH_Conc
    glycolysis_params_copy.TPI_Conc = multiplier * glycolysis_params.TPI_Conc
    glycolysis_params_copy.PGK_Conc = multiplier * glycolysis_params.PGK_Conc
    glycolysis_params_copy.ENO_Conc = multiplier * glycolysis_params.ENO_Conc
    glycolysis_params_copy.PGM_Conc = multiplier * glycolysis_params.PGM_Conc
    glycolysis_params_copy.PKM2_Conc = multiplier * glycolysis_params.PKM2_Conc
    glycolysis_params_copy.LDH_Conc = multiplier * glycolysis_params.LDH_Conc
    #glycolysis_params_copy.HK1_Conc = multiplier * glycolysis_params.HK1_Conc
    #glycolysis_params_copy.PFKP_Conc = multiplier * glycolysis_params.PFKP_Conc

    push!(glycolysis_params_list, glycolysis_params_copy)
    push!(model_names_list, "$(multiplier)x Lower Glycolysis")
end
model_names_list
Simulation_Data = []
for model_glycolysis_params in glycolysis_params_list
    res = calculate_ATP_at_range_ATPases(model_glycolysis_params, glycolysis_init_conc)
    push!(Simulation_Data, res)
end
Simulation_Data

elapsed_minutes = round(now() - start, Minute)
println("Duration: ", elapsed_minutes)
#free workers
rmprocs(workers())

##
# Plot results

#Precalculate enzyme rates
Geiger_data_summary = DataFrame(
    XLSX.readtable(
        "Data/Supplementary File 1. Levels of enzymes, metabolites and isotope tracing.xlsx",
        "Enzyme concentrations";
        infer_eltypes = true,
    ),
)

glycolysis_enzymes_Vmax = [
    glycolysis_params.GLUT_Vmax,
    glycolysis_params.HK1_Vmax,
    glycolysis_params.GPI_Vmax,
    glycolysis_params.PFKP_Vmax,
    glycolysis_params.ALDO_Vmax,
    glycolysis_params.TPI_Vmax,
    glycolysis_params.GAPDH_Vmax,
    glycolysis_params.PGK_Vmax,
    glycolysis_params.PGM_Vmax,
    glycolysis_params.ENO_Vmax,
    glycolysis_params.PKM2_Vmax_a,
    glycolysis_params.LDH_Vmax,
    glycolysis_params.MCT_Vmax,
]
glycolysis_enzymes_MW =
    1000 .* [
        glycolysis_params.GLUT_MW,
        glycolysis_params.HK1_MW,
        glycolysis_params.GPI_MW,
        glycolysis_params.PFKP_MW,
        glycolysis_params.ALDO_MW,
        glycolysis_params.TPI_MW,
        glycolysis_params.GAPDH_MW,
        glycolysis_params.PGK_MW,
        glycolysis_params.PGM_MW,
        glycolysis_params.ENO_MW,
        glycolysis_params.PKM2_MW,
        glycolysis_params.LDH_MW,
        glycolysis_params.MCT_MW,
    ]
glycolysis_gene_names = [
    "SLC2A1",
    "SLC2A2",
    "SLC2A3",
    "SLC2A4",
    "HK1",
    "HK2",
    "HK3",
    "HK4",
    "GPI",
    "PFKP",
    "PFKM",
    "PFKL",
    "ALDOA",
    "ALDOB",
    "ALDOC",
    "TPI1",
    "GAPDH",
    "PGK1",
    "PGK2",
    "PGAM1",
    "PGAM2",
    "PGAM3P",
    "PGAM4",
    #"BPGM",
    #"TIGAR",
    "ENO1",
    "ENO2",
    "ENO3",
    "PKL",
    "PKR",
    "PKM1",
    "PKM2",
    "LDHA",
    "LDHB",
    "LDHC",
    "SLC16A1",
    "SLC16A3",
    "SLC16A7",
    "SLC16A8",
]
cell_protein_density = 0.2 # mg protein/ul intracellular volume

# Plot
size_inches = (6.5, 2.5)
size_pt = 72 .* size_inches
set_theme!(
    Theme(
        fontsize = 6,
        Axis = (xticksize = 1, yticksize = 1, yticklabelpad = 1, ylabelpadding = 3),
    ),
)
fig = Figure(size = size_pt)
CairoMakie.activate!(type = "svg")

#Plot enzyme concentrations
ax_enz_Vmax =
    fig[1, 1] = Axis(
        fig,
        yscale = log10,
        xticklabelrotation = pi / 2,
        ylabel = "Cellular Enzyme Concentration, M",
        # limits = ((0, 14), (1e0, 1e4)),
        # alignmode=Mixed(right=-70),
    )
tick_labels = []
global position_counter = 1
for j = 1:nrow(Geiger_data_summary)
    enzyme_conc = [Geiger_data_summary[j, k] for k = 5:ncol(Geiger_data_summary)]
    #Skip gene if it has no data
    if all(ismissing.(enzyme_conc))
        continue
        #Skip gene if it has more than the fraction below of missing data
    elseif sum(ismissing.(enzyme_conc)) / length(enzyme_conc) > 0.7
        continue
    end
    enzyme_conc = 1_000 .* collect(skipmissing(enzyme_conc)) ./ glycolysis_enzymes_MW[j]
    positions = (position_counter) * ones(length(enzyme_conc))
    boxplot!(ax_enz_Vmax, positions, enzyme_conc, show_outliers = false)
    scatter!(
        ax_enz_Vmax,
        positions .- 0.2 .+ 0.4 .* rand(length(enzyme_conc)),
        enzyme_conc,
        strokewidth = 0.1,
        markersize = 3,
    )
    global position_counter += 1
    push!(tick_labels, Geiger_data_summary[j, :Enzyme])
end

ax_enz_Vmax.xticks = (1:position_counter-1, tick_labels)



#Plot maximal enzyme activity data
ax_enz_Vmax =
    fig[1, 2] = Axis(
        fig,
        yscale = log10,
        xticklabelrotation = pi / 2,
        ylabel = "Maximal Cellular Enzyme Activity, mM/min",
        limits = ((0, 14), (1e0, 1e4)),
        # alignmode=Mixed(right=-70),
    )
global position_counter = 1
tick_labels = []
for j = 1:nrow(Geiger_data_summary)
    enzyme_conc = [Geiger_data_summary[j, k] for k = 5:ncol(Geiger_data_summary)]
    #Skip gene if it has no data
    if all(ismissing.(enzyme_conc))
        continue
        #Skip gene if it has more than the fraction below of missing data
    elseif sum(ismissing.(enzyme_conc)) / length(enzyme_conc) > 0.7
        continue
    end
    enzyme_conc =
        1000 .* collect(skipmissing(enzyme_conc)) .* cell_protein_density .*
        glycolysis_enzymes_Vmax[j]
    positions = (position_counter) * ones(length(enzyme_conc))
    boxplot!(ax_enz_Vmax, positions, enzyme_conc, show_outliers = false)
    scatter!(
        ax_enz_Vmax,
        positions .- 0.2 .+ 0.4 .* rand(length(enzyme_conc)),
        enzyme_conc,
        strokewidth = 0.1,
        markersize = 3,
    )
    global position_counter += 1
    push!(tick_labels, Geiger_data_summary[j, :Enzyme])
end
HK_mean =
    1000 *
    median(
        Array(
            Geiger_data_summary[Geiger_data_summary.Enzyme.=="HK", Between(:A549, :U2OS)],
        ),
    ) *
    cell_protein_density *
    glycolysis_params.HK1_Vmax
lines!(
    0:14,
    [HK_mean .* ones(7); 2 .* HK_mean .* ones(8)],
    color = "Grey",
    linestyle = :dash,
)
ax_enz_Vmax.xticks = (1:position_counter-1, tick_labels)


#Plot model with more or less enzymes in lower glycolysis
ax_ATP_conc = Axis(
    fig[1, 3],
    limits = (
        (0.003, 0.45),
        (
            0.0,
            1.6 * (
                glycolysis_init_conc.ATP +
                glycolysis_init_conc.ADP +
                glycolysis_init_conc.AMP
            ),
        ),
    ),
    xscale = log10,
    title = "Effect of Lower Glycolysis Enzyme Levels\non Maintenance of ATP",
    xlabel = rich("ATPase, % of pathway V", subscript("max")),
    ylabel = "[ATP],mM",
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3]),
    xtickformat = xs -> [
        x * 100 >= 1 ? "$(Int(round(x*100, sigdigits=1)))%" :
        "$(round(x*100, sigdigits=1))%" for x in xs
    ],
    ytickformat = ys -> ["$(round(y*1000, sigdigits=3))" for y in ys],
)

global color = Makie.wong_colors()
global linewidth = 2
for (i, model_res) in enumerate(Simulation_Data)
    supported_ATPase_range =
        model_res.range_ATPases[model_res.ATP_conc.>0.5*(glycolysis_init_conc.ATP+glycolysis_init_conc.ADP+glycolysis_init_conc.AMP)]
    supported_ATPase_fold_range =
        !isempty(supported_ATPase_range) ?
        round(supported_ATPase_range[end] / supported_ATPase_range[1], sigdigits = 2) : 0.0
    lines!(
        ax_ATP_conc,
        model_res.range_ATPases,
        model_res.ATP_conc,
        linewidth = linewidth,
        color = model_names_list[i] != "1.0x Lower Glycolysis" ? color[i] : :Black,
        label = replace(model_names_list[i], " Lower Glycolysis" => ""),
        # label = replace(
        #     model_names_list[i],
        #     "1.0x Lower Glycolysis" => "Default Model\n($(supported_ATPase_fold_range)-fold)",
        #     " Lower Glycolysis" => "\n($(supported_ATPase_fold_range)-fold)",
        # ),
    )
    global linewidth *= 1
end
lines!(
    ax_ATP_conc,
    [0.003, 1.0],
    repeat(
        [glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP],
        2,
    ),
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
axislegend(
    ax_ATP_conc,
    ax_ATP_conc,
    # "Fold change of\nenzyme concentrations\nfrom TPI to LDH\nrelative to default model\n(fold ATPase range where\n[ATP] > 1/2 Adenine Pool)",
    "Fold change of [TPI], [GAPDH],\n[PGK], [PGM], [ENO], [PKM2], [LDH]",
    titlegap = 2,
    colgap = 3,
    patchlabelgap = 1,
    position = (0.5, 1.125),
    patchsize = (7.5, 2.5),
    orientation = :horizontal,
    framevisible = false,
)

colgap!(fig.layout, 10)
rowgap!(fig.layout, 5)

label_a =
    fig[1, 1, TopLeft()] =
        Label(fig, "A", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_b =
    fig[1, 2, TopLeft()] =
        Label(fig, "B", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_c =
    fig[1, 3, TopLeft()] =
        Label(fig, "C", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))

fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_Fig6_high_lower_glyc_enz_level_required_for_ATP_maintenance.pdf", fig, pt_per_unit = 1)
