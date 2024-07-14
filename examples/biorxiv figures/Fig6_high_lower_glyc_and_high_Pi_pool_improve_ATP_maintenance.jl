using Glycolysis
using DifferentialEquations, ProgressMeter
using CairoMakie, Dates, Printf, Statistics, StatsBase
using DataFrames, XLSX

# glycolysis_params.AK_Vmax = 0.03
# glycolysis_params.NDPK_Vmax = 0.03
# glycolysis_params.CK_Vmax = 0.0

##
# Precalculate and save outputs for each model
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
        EnsembleThreads(),
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

#Precalculate results at different Pi pool sizes
Pi_pool = (
    glycolysis_init_conc.G6P +
    glycolysis_init_conc.F6P +
    2 * glycolysis_init_conc.F16BP +
    glycolysis_init_conc.GAP +
    glycolysis_init_conc.DHAP +
    2 * glycolysis_init_conc.BPG +
    glycolysis_init_conc.ThreePG +
    glycolysis_init_conc.TwoPG +
    glycolysis_init_conc.PEP +
    glycolysis_init_conc.Phosphate
)

init_cond_list_Pi = []
model_names_list_Pi = []
# Pi_conc_list = [0.100, 0.050, Pi_pool, 0.010, 0.005, 0.00]
# Pi_conc_list = [0.050, Pi_pool, 0.010, 0.00]
Pi_conc_list = [2 * Pi_pool, Pi_pool, Pi_pool / 2, 0.00]
for conc in Pi_conc_list
    local multiplier = conc / Pi_pool
    local glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
    glycolysis_init_conc_copy.G6P = multiplier * glycolysis_init_conc.G6P
    glycolysis_init_conc_copy.F6P = multiplier * glycolysis_init_conc.F6P
    glycolysis_init_conc_copy.F16BP = multiplier * glycolysis_init_conc.F16BP
    glycolysis_init_conc_copy.GAP = multiplier * glycolysis_init_conc.GAP
    glycolysis_init_conc_copy.DHAP = multiplier * glycolysis_init_conc.DHAP
    glycolysis_init_conc_copy.BPG = multiplier * glycolysis_init_conc.BPG
    glycolysis_init_conc_copy.ThreePG = multiplier * glycolysis_init_conc.ThreePG
    glycolysis_init_conc_copy.TwoPG = multiplier * glycolysis_init_conc.TwoPG
    glycolysis_init_conc_copy.PEP = multiplier * glycolysis_init_conc.PEP
    glycolysis_init_conc_copy.Phosphate = multiplier * glycolysis_init_conc.Phosphate

    push!(init_cond_list_Pi, glycolysis_init_conc_copy)
    push!(model_names_list_Pi, "$(Int(round(conc*1000)))mM")
end

Simulation_Data_Pi = []
for model_glycolysis_init_conc in init_cond_list_Pi
    res = calculate_ATP_at_range_ATPases(glycolysis_params, model_glycolysis_init_conc)
    push!(Simulation_Data_Pi, res)
end


#Precalculate results for model with +/- AK, NPDK, CK
model_names_list_AK_NDPK_CK = ["None", "AK", "NDPK", "CK"]
glycolysis_params_no_AK_NDPK_CK = deepcopy(glycolysis_params)
glycolysis_params_no_AK_NDPK_CK.AK_Vmax = 0.0
glycolysis_params_no_AK_NDPK_CK.NDPK_Vmax = 0.0
glycolysis_params_no_AK_NDPK_CK.CK_Vmax = 0.0
param_list_AK_NDPK_CK = []
for enzyme in model_names_list_AK_NDPK_CK
    local glycolysis_params_copy = deepcopy(glycolysis_params)
    if enzyme == "None"
        glycolysis_params_copy.AK_Vmax = 0.0
        glycolysis_params_copy.NDPK_Vmax = 0.0
        glycolysis_params_copy.CK_Vmax = 0.0
    elseif enzyme == "AK"
        glycolysis_params_copy.AK_Vmax = 0.003
        glycolysis_params_copy.NDPK_Vmax = 0.0
        glycolysis_params_copy.CK_Vmax = 0.0
    elseif enzyme == "NDPK"
        glycolysis_params_copy.AK_Vmax = 0.0
        glycolysis_params_copy.NDPK_Vmax = 0.03
        glycolysis_params_copy.CK_Vmax = 0.0
    elseif enzyme == "CK"
        glycolysis_params_copy.AK_Vmax = 0.0
        glycolysis_params_copy.NDPK_Vmax = 0.0
        glycolysis_params_copy.CK_Vmax = 0.03
    end
    push!(param_list_AK_NDPK_CK, glycolysis_params_copy)
end

multiplier = 1.0
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
glycolysis_init_conc_copy.Phosphate = multiplier * glycolysis_init_conc.Phosphate

Simulation_Data_AK_NDPK_CK = []
for model_glycolysis_params in param_list_AK_NDPK_CK
    res = calculate_ATP_at_range_ATPases(model_glycolysis_params, glycolysis_init_conc_copy)
    push!(Simulation_Data_AK_NDPK_CK, res)
end






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
size_inches = (6.5, 5)
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
            1.4 * (
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

#Plot a pie chart of phosphate containing metabolite levels
ax_Pi_pie = Axis(
    fig[2, 1],
    title = "Cellular Phosphate Pools",
    autolimitaspect = 1,
    limits = (nothing, (-5.5, 5.5)),
)

#Download experimental data and process it
Experimental_Data = DataFrame(
    XLSX.readtable(
        "Data/Supplementary File 1. Levels of enzymes, metabolites and isotope tracing.xlsx",
        "Metabolite concentrations";
        infer_eltypes = true,
    ),
)
# List metabolites that contain Pi
Pi_metabolites_names = [
    "G6P",
    "F6P",
    "F16BP",
    "GAP",
    "DHAP",
    "BPG",
    "ThreePG",
    "TwoPG",
    "PEP",
    "Phosphate",
    "ATP",
    "ADP",
]
Experimental_Data[:, :ATP] *= 3
Experimental_Data[:, :ADP] *= 2
Experimental_Data[:, :AMP] *= 1
Experimental_Data[:, :F16BP] *= 2
Experimental_Data[:, :BPG] *= 2
# Create a new DataFrame with the mean of each column
Experimental_Data_Mean_DF = combine(
    Experimental_Data[:, Pi_metabolites_names],
    names(Experimental_Data[:, Pi_metabolites_names]) .=>
        mean ∘ skipmissing .=> names(Experimental_Data[:, Pi_metabolites_names]),
)
#Stack DataFrame
Experimental_Data_Mean_DF_sorted = sort(
    stack(
        Experimental_Data_Mean_DF,
        Pi_metabolites_names,
        variable_name = :metab_name,
        value_name = :conc,
    ),
    :conc,
    rev = true,
)
Experimental_Data_Mean_DF_sorted.conc /= Glycolysis.cell_volume_correction
#Make a new DF with the sum of the three Pi pools
Pi_pool_summary_df = DataFrame(
    pool_names = [
        "Phosphate in\nATP+ADP+AMP",
        "Inorganic\nPhosphate",
        "Phosphate in\nGlycolytic\nIntermediates",
    ],
    conc = [
        sum(
            Experimental_Data_Mean_DF_sorted[
                Experimental_Data_Mean_DF_sorted.metab_name.∈Ref(["ATP", "ADP", "AMP"]),
                :conc,
            ],
        ),
        sum(
            Experimental_Data_Mean_DF_sorted[
                Experimental_Data_Mean_DF_sorted.metab_name.∈Ref(["Phosphate"]),
                :conc,
            ],
        ),
        sum(
            Experimental_Data_Mean_DF_sorted[
                Experimental_Data_Mean_DF_sorted.metab_name.∈Ref([
                    "G6P",
                    "F6P",
                    "F16BP",
                    "GAP",
                    "DHAP",
                    "BPG",
                    "ThreePG",
                    "TwoPG",
                    "PEP",
                ]),
                :conc,
            ],
        ),
    ],
)
#Actually plot the data
r1 = 6
r0 = 2
pie!(
    ax_Pi_pie,
    Pi_pool_summary_df.conc,
    color = [Makie.wong_colors()[1+i%7] for i = 1:length(Pi_pool_summary_df.conc)],
    radius = r1,
    inner_radius = r0,
    strokecolor = :white,
    strokewidth = 1,
)

rot(θ) = 0 < θ < 120 ? θ * π / 180 - π / 2 : θ * π / 180 + π / 2
θ =
    (cumsum(Pi_pool_summary_df.conc) - Pi_pool_summary_df.conc / 2) .* 360 /
    sum(Pi_pool_summary_df.conc)
scθ = sincosd.(θ)
for (li, θi, sci) in zip(Pi_pool_summary_df.pool_names, θ, scθ)
    # ri = sqrt(r0 * r1)
    ri = median([r0, r1])
    text!(
        ax_Pi_pie,
        li *
        "\n$(round(1000*Pi_pool_summary_df.conc[findfirst(Pi_pool_summary_df.pool_names .== li)], sigdigits = 2))mM",
        position = ri .* (sci[2], sci[1]),
        rotation = rot(θi),
        justification = :center,
        align = (:center, :center),
    )
end
# print sum of pool in the middle of plot
text!(
    ax_Pi_pie,
    "Total\nPhosphate\n$(Int(round(1000*sum(Pi_pool_summary_df.conc), sigdigits = 2)))mM",
    align = (:center, :center),
    justification = :center,
)
hidedecorations!(ax_Pi_pie)
hidespines!(ax_Pi_pie)


#Plot model with more or less phosphate pool
ax_ATP_conc = Axis(
    fig[2, 2],
    limits = (
        (0.003, 0.45),
        (
            0.0,
            1.4 * (
                glycolysis_init_conc.ATP +
                glycolysis_init_conc.ADP +
                glycolysis_init_conc.AMP
            ),
        ),
    ),
    xscale = log10,
    title = "Effect of Non-Adenine Phosphate Pool Size\non Maintenance of ATP",
    xlabel = rich("ATPase, % of pathway V", subscript("max")),
    ylabel = "[ATP],mM",
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3]),
    xtickformat = xs -> [
        x * 100 >= 1 ? "$(Int(round(x*100, sigdigits=1)))%" :
        "$(round(x*100, sigdigits=1))%" for x in xs
    ],
    ytickformat = ys -> ["$(round(y*1000, sigdigits=3))" for y in ys],
)
color = Makie.wong_colors()
for (i, model_res) in enumerate(Simulation_Data_Pi)
    supported_ATPase_range =
        model_res.range_ATPases[model_res.ATP_conc.>0.5*(glycolysis_init_conc.ATP+glycolysis_init_conc.ADP+glycolysis_init_conc.AMP)]
    supported_ATPase_fold_range =
        round(supported_ATPase_range[end] / supported_ATPase_range[1], sigdigits = 2)
    lines!(
        ax_ATP_conc,
        model_res.range_ATPases,
        model_res.ATP_conc,
        color = model_names_list_Pi[i] != "20mM" ? color[i] : :Black,
        label = model_names_list_Pi[i],
        # label = replace(
        #     model_names_list_Pi[i],
        #     "$(Int(round(Pi_pool*1000)))mM" => "Default Model\n$(Int(round(Pi_pool*1000)))mM\n($(supported_ATPase_fold_range)-fold)",
        #     "mM" => "mM\n($(supported_ATPase_fold_range)-fold)",
        # ),
    )
end
lines!(
    ax_ATP_conc,
    [0.001, 1.0],
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
    # "Phosphate Pool, mM\n(fold ATPase range where\n[ATP] > 1/2 Adenine Pool)",
    "Non-Adenine Phosphate Pool, mM",
    titlegap = 2,
    colgap = 3,
    patchlabelgap = 1,
    position = (0.5, 1.125),
    patchsize = (7.5, 2.5),
    orientation = :horizontal,
    framevisible = false,
)

#Plot model with more or less phosphate pool and +/- AK, NPDK, CK
ax_ATP_conc = Axis(
    fig[2, 3],
    limits = (
        (0.003, 0.45),
        (
            0.0,
            1.4 * (
                glycolysis_init_conc.ATP +
                glycolysis_init_conc.ADP +
                glycolysis_init_conc.AMP
            ),
        ),
    ),
    xscale = log10,
    title = "Effect of AK, NDPK, CK\non Maintenance of ATP",
    xlabel = rich("ATPase, % of pathway V", subscript("max")),
    ylabel = "[ATP],mM",
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3]),
    xtickformat = xs -> [
        x * 100 >= 1 ? "$(Int(round(x*100, sigdigits=1)))%" :
        "$(round(x*100, sigdigits=1))%" for x in xs
    ],
    ytickformat = ys -> ["$(round(y*1000, sigdigits=3))" for y in ys],
)
color = Makie.wong_colors()
for (i, model_res) in enumerate(Simulation_Data_AK_NDPK_CK)
    # supported_ATPase_range =
    #     model_res.range_ATPases[model_res.ATP_conc.>0.5*(glycolysis_init_conc.ATP+glycolysis_init_conc.ADP+glycolysis_init_conc.AMP)]
    # supported_ATPase_fold_range =
    #     round(supported_ATPase_range[end] / supported_ATPase_range[1], sigdigits = 2)
    lines!(
        ax_ATP_conc,
        model_res.range_ATPases,
        model_res.ATP_conc,
        color = model_names_list_AK_NDPK_CK[i] != "AK" ? color[i] : :Black,
        label = model_names_list_AK_NDPK_CK[i],
        # label = replace(
        #     model_names_list_AK_NDPK_CK[i],
        #     "$(Int(round(Pi_pool*1000)))mM" => "Default Model\n$(Int(round(Pi_pool*1000)))mM\n($(supported_ATPase_fold_range)-fold)",
        #     "mM" => "mM\n($(supported_ATPase_fold_range)-fold)",
        # ),
    )
end
lines!(
    ax_ATP_conc,
    [0.001, 1.0],
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
    "Included Enzymes",
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
label_d =
    fig[2, 1, TopLeft()] =
        Label(fig, "D", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_e =
    fig[2, 2, TopLeft()] =
        Label(fig, "E", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_f =
    fig[2, 3, TopLeft()] =
        Label(fig, "F", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))

fig

# uncomment the line below to save the plot
save("Results/$(Dates.format(now(),"mmddyy"))_Fig6_high_lower_glyc_and_Pi_pool_required_for_ATP_maintenance.png", fig, px_per_unit = 4)
