using Glycolysis
using DifferentialEquations, ProgressMeter
using CairoMakie, Dates, Printf, Statistics, StatsBase
using DataFrames, CSV


##
# Make vector of glycolysis_params of models to be plotted
# params_list = [glycolysis_params]
# model_names_list = ["Complete_Model"]
params_list = []
model_names_list = []
params_copy = deepcopy(glycolysis_params)
params_copy.HK1_K_a_G6P_cat = Inf
params_copy.HK1_K_i_G6P_reg = Inf
params_copy.HK1_K_a_Pi = Inf
push!(params_list, params_copy)
push!(model_names_list, "No_HK1_allost._reg.")

params_copy = deepcopy(glycolysis_params)
params_copy.PFKP_L = 0.0
push!(params_list, params_copy)
push!(model_names_list, "No_PFKP_allost._reg.")

params_copy = deepcopy(glycolysis_params)
params_copy.GAPDH_L = 0.0
params_copy.PKM2_L = 0.0
push!(params_list, params_copy)
push!(model_names_list, "No_GAPDH_allost._reg._and_No_PKM2_allost._reg.")

params_copy = deepcopy(glycolysis_params)
params_copy.HK1_K_a_G6P_cat = Inf
params_copy.HK1_K_i_G6P_reg = Inf
params_copy.HK1_K_a_Pi = Inf
params_copy.PFKP_L = 0.0
push!(params_list, params_copy)
push!(model_names_list, "No_HK1_allost._reg._and_No_PFKP_allost._reg.")

params_copy = deepcopy(glycolysis_params)
params_copy.HK1_K_a_G6P_cat = Inf
push!(params_list, params_copy)
push!(model_names_list, "No_HK1_cat_G6P_inh.")

params_copy = deepcopy(glycolysis_params)
params_copy.HK1_K_i_G6P_reg = Inf
push!(params_list, params_copy)
push!(model_names_list, "No_HK1_allo_G6P_inh.")

params_copy = deepcopy(glycolysis_params)
params_copy.PFKP_K_i_ATP_reg = Inf
push!(params_list, params_copy)
push!(model_names_list, "No_PFKP_ATP_inh.")

params_copy = deepcopy(glycolysis_params)
params_copy.HK1_K_i_G6P_reg = Inf
params_copy.PFKP_K_i_ATP_reg = Inf
params_copy.HK1_K_a_G6P_cat = Inf
push!(params_list, params_copy)
push!(model_names_list, "No_HK1_G6P_inh._and_No_PFKP_ATP_inh.")

params_copy = deepcopy(glycolysis_params)
params_copy.HK1_K_a_Pi = Inf
push!(params_list, params_copy)
push!(model_names_list, "No_HK1_Pi_act.")
params_copy = deepcopy(glycolysis_params)
params_copy.PFKP_K_Phosphate = Inf
push!(params_list, params_copy)
push!(model_names_list, "No_PFKP_Pi_act.")
params_copy = deepcopy(glycolysis_params)
params_copy.PFKP_K_a_ADP_reg = Inf
push!(params_list, params_copy)
push!(model_names_list, "No_PFKP_ADP_act.")

params_copy = deepcopy(glycolysis_params)
params_copy.HK1_K_a_Pi = Inf
params_copy.PFKP_K_Phosphate = Inf
params_copy.PFKP_K_a_ADP_reg = Inf
push!(params_list, params_copy)
push!(model_names_list, "No_HK1_Pi_act._and_No_PFKP_ADP,_Pi_act.")
@assert length(model_names_list) == length(params_list)

##
# Precalculate and save outputs for each model
# This calculation can takes ~10 minutes on an 8 core machine

function find_ATP_at_ATPase_range(params, init_conc; n_Vmax_ATPase_values = 1000)
    tspan = (0.0, 1e8)
    pathway_Vmax = 2 * params.HK1_Vmax * params.HK1_Conc
    ATPases = 10 .^ range(log10(0.001), log10(1.0), n_Vmax_ATPase_values) .* pathway_Vmax
    prob = ODEProblem(glycolysis_ODEs, init_conc, tspan, params)
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
        save_start = false,
    )
    ATP_conc = [sol.u[end].ATP for sol in sim if sol.retcode == ReturnCode.Success]
    ATPase_Vmax = [ATPases[i] / (2 * params.HK1_Vmax * params.HK1_Conc) for (i, sol) in enumerate(sim) if sol.retcode == ReturnCode.Success]
    return (ATP_conc = ATP_conc, ATPase_Vmax = ATPase_Vmax)
end

Complete_Model_Simulation_Data = find_ATP_at_ATPase_range(glycolysis_params, glycolysis_init_conc)
Complete_Model_Simulation_Data
no_reg_params = deepcopy(glycolysis_params)
no_reg_params.HK1_K_a_G6P_cat = Inf
no_reg_params.HK1_K_i_G6P_reg = Inf
no_reg_params.HK1_K_a_Pi = Inf
no_reg_params.PFKP_L = 0.0
no_reg_params.GAPDH_L = 0.0
no_reg_params.PKM2_L = 0.0
No_Reg_Model_Simulation_Data = find_ATP_at_ATPase_range(no_reg_params, glycolysis_init_conc)
No_Reg_Model_Simulation_Data
Simulation_Data = []
for (i, model_params) in enumerate(params_list)
    println(model_names_list[i])
    res = find_ATP_at_ATPase_range(model_params, glycolysis_init_conc)
    push!(Simulation_Data, res)
end
Simulation_Data


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
fig = Figure(size = size_pt)

#Loop through models and create all the axis for plotting

traces_panels = Vector{typeof(GridLayout())}(undef, length(model_names_list))
# default_color = Makie.wong_colors()[1]
default_color = :Black
default_linestyle = :dot
variant_color = Makie.wong_colors()[3]
variant_linestyle = nothing
# no_reg_color = Makie.wong_colors()[6]
no_reg_color = :Red
no_reg_linestyle = :dot

#Loop through models and plot everything
adenine_pool_size = glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP
for (i, model_name) in enumerate(model_names_list)
    ax_ATP_conc = Axis(
        if i <= 4
            fig[1, i]
        elseif i > 4 && i <= 8
            fig[2, i-4]
        elseif i > 8 && i <= 12
            fig[3, i-8]
        else
            @error "More than 6 plots won't fit"
        end,
        limits = ((0.001, 1.0), (-0.3e-3, adenine_pool_size * 1.5)),
        xscale = log10,
        xlabel = "ATPase, % of pathway Vmax",
        ylabel = "[ATP],mM",
        xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9], ["0.3%", "1%", "3%", "10%", "30%", "90%"]),
        ytickformat = ys -> ["$(Int(round(y*1000, sigdigits=3)))" for y in ys],
        xlabelpadding = 0.0,
        xticklabelpad = 0.0,
    )

    lines!(
        ax_ATP_conc,
        Simulation_Data[i].ATPase_Vmax,
        Simulation_Data[i].ATP_conc,
        linestyle = variant_linestyle,
        color = variant_color,
    )
    lines!(
        ax_ATP_conc,
        Complete_Model_Simulation_Data.ATPase_Vmax,
        Complete_Model_Simulation_Data.ATP_conc,
        linestyle = default_linestyle,
        color = default_color,
    )
    lines!(
        ax_ATP_conc,
        No_Reg_Model_Simulation_Data.ATPase_Vmax,
        No_Reg_Model_Simulation_Data.ATP_conc,
        linestyle = no_reg_linestyle,
        color = no_reg_color,
    )
    axislegend(
        ax_ATP_conc,
        [[LineElement(linestyle = variant_linestyle, color = variant_color)]],
        [replace(model_name, "_and_" => "\n", "_" => " ", "wo" => "w/o")],
        position = :lt,
        rowgap = 2,
        padding = (0, 0, 0, -4),
        patchsize = (10, 5),
        patchlabelgap = 2,
        framevisible = false,
    )
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
end

Legend(
    fig[0, 1:4],
    [
        [LineElement(linestyle = default_linestyle, color = default_color)],
        [LineElement(linestyle = no_reg_linestyle, color = no_reg_color)],
    ],
    ["Full Model", "Model w/o allost."],
    rowgap = 2,
    padding = (2, 2, 2, 2),
    patchsize = (10, 5),
    patchlabelgap = 2,
    # framevisible = false,
    halign = :center,
    valign = :bottom,
    orientation = :horizontal,
)

rowgap!(fig.layout, 5)
colgap!(fig.layout, 10)


label_a = fig[1, 1, TopLeft()] = Label(fig, "A", fontsize = 12, halign = :right, padding = (0, 7, -7, 0))
label_b = fig[1, 2, TopLeft()] = Label(fig, "B", fontsize = 12, halign = :right, padding = (0, 7, -7, 0))
label_c = fig[1, 3, TopLeft()] = Label(fig, "C", fontsize = 12, halign = :right, padding = (0, 7, -7, 0))
label_d = fig[1, 4, TopLeft()] = Label(fig, "D", fontsize = 12, halign = :right, padding = (0, 7, -7, 0))
label_e = fig[2, 1, TopLeft()] = Label(fig, "E", fontsize = 12, halign = :right, padding = (0, 7, -7, 0))
label_f = fig[2, 2, TopLeft()] = Label(fig, "F", fontsize = 12, halign = :right, padding = (0, 7, -7, 0))
label_g = fig[2, 3, TopLeft()] = Label(fig, "G", fontsize = 12, halign = :right, padding = (0, 7, -7, 0))
label_h = fig[2, 4, TopLeft()] = Label(fig, "H", fontsize = 12, halign = :right, padding = (0, 7, -7, 0))
label_k = fig[3, 1, TopLeft()] = Label(fig, "K", fontsize = 12, halign = :right, padding = (0, 7, -7, 0))
label_l = fig[3, 2, TopLeft()] = Label(fig, "L", fontsize = 12, halign = :right, padding = (0, 7, -7, 0))
label_m = fig[3, 3, TopLeft()] = Label(fig, "M", fontsize = 12, halign = :right, padding = (0, 7, -7, 0))
label_n = fig[3, 4, TopLeft()] = Label(fig, "N", fontsize = 12, halign = :right, padding = (0, 7, -7, 0))

fig

# uncomment the line below to save the plot
save("Results/$(Dates.format(now(),"mmddyy"))_Fig4_ATP_stability_wo_specific_reg.png", fig, px_per_unit = 4)
