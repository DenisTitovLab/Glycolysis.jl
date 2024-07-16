using Glycolysis
using DifferentialEquations, ProgressMeter, LabelledArrays
using CairoMakie, Dates, Printf, Statistics, StatsBase
using DataFrames, CSV
using FileIO

##
# Precalculate output of complete model and model without regulation
function glycolysis_ODEs_fixed_Pi(ds, s, params, t)
    ds.Glucose_media = 0.0
    ds.Glucose = Glycolysis.rate_GLUT(s, params) -
                 Glycolysis.rate_HK1(s, params)
    ds.G6P = Glycolysis.rate_HK1(s, params) -
             Glycolysis.rate_GPI(s, params)
    ds.F6P = (
        Glycolysis.rate_GPI(s, params) -
        Glycolysis.rate_PFKP(s, params)
    )
    ds.F16BP = (
        Glycolysis.rate_PFKP(s, params) -
        Glycolysis.rate_ALDO(s, params)
    )
    ds.GAP = (
        Glycolysis.rate_ALDO(s, params) + Glycolysis.rate_TPI(s, params) -
        Glycolysis.rate_GAPDH(s, params)
    )
    ds.DHAP = Glycolysis.rate_ALDO(s, params) - Glycolysis.rate_TPI(s, params)
    ds.BPG = Glycolysis.rate_GAPDH(s, params) -
             Glycolysis.rate_PGK(s, params)
    ds.ThreePG = Glycolysis.rate_PGK(s, params) -
                 Glycolysis.rate_PGM(s, params)
    ds.TwoPG = Glycolysis.rate_PGM(s, params) - Glycolysis.rate_ENO(s, params)
    ds.PEP = Glycolysis.rate_ENO(s, params) -
             Glycolysis.rate_PKM2(s, params)
    ds.Pyruvate = Glycolysis.rate_PKM2(s, params) -
                  Glycolysis.rate_LDH(s, params)
    ds.Lactate = Glycolysis.rate_LDH(s, params) -
                 Glycolysis.rate_MCT(s, params)
    ds.Lactate_media = 0.0
    ds.ATP = (
        -Glycolysis.rate_HK1(s, params) -
        Glycolysis.rate_PFKP(s, params) +
        Glycolysis.rate_PGK(s, params) +
        Glycolysis.rate_PKM2(s, params) -
        Glycolysis.rate_ATPase(s, params) +
        Glycolysis.rate_AK(s, params)
    )
    ds.ADP = (
        Glycolysis.rate_HK1(s, params) +
        Glycolysis.rate_PFKP(s, params) -
        Glycolysis.rate_PGK(s, params) -
        Glycolysis.rate_PKM2(s, params) +
        Glycolysis.rate_ATPase(s, params) -
        2 * Glycolysis.rate_AK(s, params)
    )
    ds.AMP = Glycolysis.rate_AK(s, params)
    ds.Phosphate = 0.0
    ds.NAD = Glycolysis.rate_LDH(s, params) -
             Glycolysis.rate_GAPDH(s, params)
    ds.NADH = Glycolysis.rate_GAPDH(s, params) -
              Glycolysis.rate_LDH(s, params)
    ds.F26BP = 0.0
    ds.Citrate = 0.0
    ds.Phenylalanine = 0.0
end

function find_ATP_at_ATPase_range(
        ODE_func,
        params,
        init_conc;
        n_Vmax_ATPase_values = 1000,
        min_ATPase = 0.001,
        max_ATPase = 1.0
)
    tspan = (0.0, 1e8)
    pathway_Vmax = 2 * params.HK1_Vmax * params.HK1_Conc
    ATPases = 10 .^ range(log10(min_ATPase), log10(max_ATPase), n_Vmax_ATPase_values) .*
              pathway_Vmax
    prob = ODEProblem(ODE_func, init_conc, tspan, params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(
        ensemble_prob,
        RadauIIA5(),
        EnsembleSerial(),
        trajectories = n_Vmax_ATPase_values,
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
        save_start = false,
        maxiters = 1e3
    )
    ATP_conc = [sol.u[end].ATP for sol in sim if sol.retcode == ReturnCode.Success]
    Glucose_conc = [sol.u[end].Glucose for sol in sim if sol.retcode == ReturnCode.Success]
    G6P_conc = [sol.u[end].G6P for sol in sim if sol.retcode == ReturnCode.Success]
    F6P_conc = [sol.u[end].F6P for sol in sim if sol.retcode == ReturnCode.Success]
    F16BP_conc = [sol.u[end].F16BP for sol in sim if sol.retcode == ReturnCode.Success]
    GAP_conc = [sol.u[end].GAP for sol in sim if sol.retcode == ReturnCode.Success]
    DHAP_conc = [sol.u[end].DHAP for sol in sim if sol.retcode == ReturnCode.Success]
    BPG_conc = [sol.u[end].BPG for sol in sim if sol.retcode == ReturnCode.Success]
    ThreePG_conc = [sol.u[end].ThreePG for sol in sim if sol.retcode == ReturnCode.Success]
    TwoPG_conc = [sol.u[end].TwoPG for sol in sim if sol.retcode == ReturnCode.Success]
    PEP_conc = [sol.u[end].PEP for sol in sim if sol.retcode == ReturnCode.Success]
    Pyruvate_conc = [sol.u[end].Pyruvate
                     for sol in sim if sol.retcode == ReturnCode.Success]
    Lactate_conc = [sol.u[end].Lactate for sol in sim if sol.retcode == ReturnCode.Success]

    ATPase_Vmax = [ATPases[i]
                   for (i, sol) in enumerate(sim) if sol.retcode == ReturnCode.Success] ./
                  (2 * params.HK1_Vmax * params.HK1_Conc)
    Q_Keq_PFKP = [Glycolysis.conc_to_disequilibrium_ratios(sol.u[end], params).Q_Keq_PFKP
                  for
                  sol in sim if sol.retcode == ReturnCode.Success]
    Q_Keq_HK1 = [Glycolysis.conc_to_disequilibrium_ratios(sol.u[end], params).Q_Keq_HK1
                 for
                 sol in sim if sol.retcode == ReturnCode.Success]
    Q_Keq_PGK = [Glycolysis.conc_to_disequilibrium_ratios(sol.u[end], params).Q_Keq_PGK
                 for
                 sol in sim if sol.retcode == ReturnCode.Success]
    Q_Keq_PKM2 = [Glycolysis.conc_to_disequilibrium_ratios(sol.u[end], params).Q_Keq_PKM2
                  for
                  sol in sim if sol.retcode == ReturnCode.Success]
    return (
        ATP_conc = ATP_conc,
        ATPase_Vmax = ATPase_Vmax,
        Glucose_conc = Glucose_conc,
        G6P_conc = G6P_conc,
        F6P_conc = F6P_conc,
        F16BP_conc = F16BP_conc,
        GAP_conc = GAP_conc,
        DHAP_conc = DHAP_conc,
        BPG_conc = BPG_conc,
        ThreePG_conc = ThreePG_conc,
        TwoPG_conc = TwoPG_conc,
        PEP_conc = PEP_conc,
        Pyruvate_conc = Pyruvate_conc,
        Lactate_conc = Lactate_conc,
        Q_Keq_PFKP = Q_Keq_PFKP,
        Q_Keq_HK1 = Q_Keq_HK1,
        Q_Keq_PGK = Q_Keq_PGK,
        Q_Keq_PKM2 = Q_Keq_PKM2
    )
end

glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
glycolysis_params.ATPase_Km_ATP = 1e-9
Complete_Model_Simulation_Data = find_ATP_at_ATPase_range(
    glycolysis_ODEs, glycolysis_params, glycolysis_init_conc_copy)

no_reg_params = deepcopy(glycolysis_params)
no_reg_params.HK1_K_a_G6P_cat = Inf
no_reg_params.HK1_K_i_G6P_reg = Inf
no_reg_params.HK1_K_a_Pi = Inf
no_reg_params.PFKP_L = 0.0
no_reg_params.GAPDH_L = 0.0
no_reg_params.PKM2_L = 0.0
no_reg_params.ATPase_Km_ATP = 1e-9

glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
No_Reg_Model_Simulation_Data = find_ATP_at_ATPase_range(
    glycolysis_ODEs, no_reg_params, glycolysis_init_conc_copy)

glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
glycolysis_init_conc_copy.Phosphate = 1e-3
Const_Pi_No_Reg_Model_Simulation_Data = find_ATP_at_ATPase_range(
    glycolysis_ODEs_fixed_Pi, no_reg_params, glycolysis_init_conc_copy)

##
# Plot the results
size_inches = (5, 3)
size_pt = 72 .* size_inches
set_theme!(Theme(fontsize = 6,
    Axis = (
        xticksize = 1,
        yticksize = 1,
        xticklabelsize = 8,
        yticklabelsize = 8,
        yticklabelpad = 1,
        ylabelpadding = 3
    )))
fig = Figure(size = size_pt)

full_color = Makie.wong_colors()[1]
no_reg_color = Makie.wong_colors()[6]
full_color = :Black
no_reg_color = :Red
condition_no_reg_color = Makie.wong_colors()[3]
full_linestyle = :dot
no_reg_linestyle = :dot
condition_no_reg_linestyle = nothing
adenine_pool_color = :Grey
adenine_pool_linestyle = :dash

line_ATPase_width = 2
line_ATPase_color = :grey
# no_allo_line_models_style = :dot
# no_reg_linestyle = [0.5, 1.5, 2.5, 3.5]
no_allo_line_models_style = [0.5, 1, 1.5, 2] .* 2

# Plot [ATP] maintenance with constant Phosphate
ax_ATPase_range = Axis(
    fig[1, 1],
    limits = ((0.001, 1.0), (1e-3, nothing)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "[Metabolite],M",
    title = "Metabolite concentration changes with model\nwithout allostery at constant [Phosphate]=1mM",
    titlesize = 10,
    xlabelsize = 12,
    ylabelsize = 12,
    xscale = log10,
    yscale = log10,
    xticks = ([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.9],
        ["0.1%", "0.3%", "1%", "3%", "10%", "30%", "90%"])
)
unphosphorylated_color = :Grey
unphosphorylated_style = :dot
after_gapdh_color = :Grey
after_gapdh_style = :dash

lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    Const_Pi_No_Reg_Model_Simulation_Data.Glucose_conc,
    label = "Glucose",
    color = unphosphorylated_color,
    linestyle = unphosphorylated_style
)
lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    Const_Pi_No_Reg_Model_Simulation_Data.G6P_conc,
    label = "G6P"
)
lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    Const_Pi_No_Reg_Model_Simulation_Data.F6P_conc,
    label = "F6P"
)
lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    Const_Pi_No_Reg_Model_Simulation_Data.F16BP_conc,
    label = "F16BP"
)
lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    abs.(Const_Pi_No_Reg_Model_Simulation_Data.GAP_conc),
    label = "GAP"
)
lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    abs.(Const_Pi_No_Reg_Model_Simulation_Data.DHAP_conc),
    label = "DHAP"
)
lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    abs.(Const_Pi_No_Reg_Model_Simulation_Data.BPG_conc),
    label = "BPG",
    color = after_gapdh_color,
    linestyle = after_gapdh_style
)
lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    abs.(Const_Pi_No_Reg_Model_Simulation_Data.ThreePG_conc),
    label = "3PG",
    color = after_gapdh_color,
    linestyle = after_gapdh_style
)
lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    abs.(Const_Pi_No_Reg_Model_Simulation_Data.TwoPG_conc),
    label = "2PG",
    color = after_gapdh_color,
    linestyle = after_gapdh_style
)
lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    abs.(Const_Pi_No_Reg_Model_Simulation_Data.PEP_conc),
    label = "PEP",
    color = after_gapdh_color,
    linestyle = after_gapdh_style
)
lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    Const_Pi_No_Reg_Model_Simulation_Data.Pyruvate_conc,
    label = "Pyruvate",
    color = unphosphorylated_color,
    linestyle = unphosphorylated_style
)
lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    Const_Pi_No_Reg_Model_Simulation_Data.Lactate_conc,
    label = "Lactate",
    color = unphosphorylated_color,
    linestyle = unphosphorylated_style
)
Legend(
    fig[1, 2],
    ax_ATPase_range,
    # position = :lt,
    # position = (1.075, 1.05),
    rowgap = 1,
    # framevisible = false,
    # padding = (-5, -5, 0, -8),
    patchsize = (20, 5),
    labelsize = 8
)

colgap!(fig.layout, 10)
rowgap!(fig.layout, 10)

# label_a = fig[1, 1, TopLeft()] = Label(fig, "A", fontsize = 12, halign = :right, padding = (0, 5, 5, 0))

fig

# uncomment the line below to save the plot
save("Results/$(Dates.format(now(),"mmddyy"))_Fig5_fig_suppl1_glyc_intermediate_at_constant_Pi.png", fig, px_per_unit = 4)
