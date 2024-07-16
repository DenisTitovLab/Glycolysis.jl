using Glycolysis
using DifferentialEquations, ProgressMeter, LabelledArrays
using CairoMakie, Dates, Printf, Statistics, StatsBase
using DataFrames, CSV
using FileIO


##
# Precalculate output of complete model and model without regulation

function glycolysis_ODEs_fixed_Pi(ds, s, params, t)
    ds.Glucose_media = 0.0
    ds.Glucose =
        Glycolysis.rate_GLUT(s.Glucose_media, s.Glucose, params) -
        Glycolysis.rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params)
    ds.G6P =
        Glycolysis.rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) -
        Glycolysis.rate_GPI(s.G6P, s.F6P, params)
    ds.F6P = (
        Glycolysis.rate_GPI(s.G6P, s.F6P, params) -
        Glycolysis.rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params)
    )
    ds.F16BP = (
        Glycolysis.rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) -
        Glycolysis.rate_ALDO(s.F16BP, s.GAP, s.DHAP, params)
    )
    ds.GAP = (
        Glycolysis.rate_ALDO(s.F16BP, s.GAP, s.DHAP, params) + Glycolysis.rate_TPI(s.GAP, s.DHAP, params) -
        Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    )
    ds.DHAP =
        Glycolysis.rate_ALDO(s.F16BP, s.GAP, s.DHAP, params) - Glycolysis.rate_TPI(s.GAP, s.DHAP, params)
    ds.BPG =
        Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params) -
        Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params)
    ds.ThreePG =
        Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) -
        Glycolysis.rate_PGM(s.ThreePG, s.TwoPG, params)
    ds.TwoPG = Glycolysis.rate_PGM(s.ThreePG, s.TwoPG, params) - Glycolysis.rate_ENO(s.TwoPG, s.PEP, params)
    ds.PEP =
        Glycolysis.rate_ENO(s.TwoPG, s.PEP, params) -
        Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params)
    ds.Pyruvate =
        Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) -
        Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params)
    ds.Lactate =
        Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params) -
        Glycolysis.rate_MCT(s.Lactate, s.Lactate_media, params)
    ds.Lactate_media = 0.0
    ds.ATP = (
        -Glycolysis.rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) -
        Glycolysis.rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) +
        Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) +
        Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) -
        Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) +
        Glycolysis.rate_AK(s.ATP, s.ADP, s.AMP, params)
    )
    ds.ADP = (
        Glycolysis.rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) +
        Glycolysis.rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) -
        Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) -
        Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) +
        Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) -
        2 * Glycolysis.rate_AK(s.ATP, s.ADP, s.AMP, params)
    )
    ds.AMP = Glycolysis.rate_AK(s.ATP, s.ADP, s.AMP, params)
    # ds.Phosphate =
    #     Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) -
    #     Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    ds.Phosphate = 0.0
    ds.NAD =
        Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params) -
        Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    ds.NADH =
        Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params) -
        Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params)
    ds.F26BP = 0.0
    ds.Citrate = 0.0
    ds.Phenylalanine = 0.0
end

function glycolysis_ODEs_Vmax_HK_PFK_eq_ATPase(ds, s, params, t)
    # 0.5 * params.ATPase_Vmax * (1 - s.G6P * s.ADP / (params.HK1_Keq * s.Glucose * s.ATP))
    # 0.5 * params.ATPase_Vmax * (1 - s.F16BP * s.ADP / (params.PFKP_Keq * s.F6P * s.ATP))
    # params.ATPase_Vmax * (1 - s.Phosphate * s.ADP / (params.ATPase_Keq * s.ATP))
    # Glycolysis.rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params)
    # Glycolysis.rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params)
    # Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params)
    ds.Glucose_media = 0.0
    ds.Glucose =
        Glycolysis.rate_GLUT(s.Glucose_media, s.Glucose, params) -
        0.5 * Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params)
    ds.G6P =
        0.5 * Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) -
        Glycolysis.rate_GPI(s.G6P, s.F6P, params)
    ds.F6P = (
        Glycolysis.rate_GPI(s.G6P, s.F6P, params) -
        0.5 * Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params)
    )
    ds.F16BP = (
        0.5 * Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) -
        Glycolysis.rate_ALDO(s.F16BP, s.GAP, s.DHAP, params)
    )
    ds.GAP = (
        Glycolysis.rate_ALDO(s.F16BP, s.GAP, s.DHAP, params) + Glycolysis.rate_TPI(s.GAP, s.DHAP, params) -
        Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    )
    ds.DHAP =
        Glycolysis.rate_ALDO(s.F16BP, s.GAP, s.DHAP, params) - Glycolysis.rate_TPI(s.GAP, s.DHAP, params)
    ds.BPG =
        Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params) -
        Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params)
    ds.ThreePG =
        Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) -
        Glycolysis.rate_PGM(s.ThreePG, s.TwoPG, params)
    ds.TwoPG = Glycolysis.rate_PGM(s.ThreePG, s.TwoPG, params) - Glycolysis.rate_ENO(s.TwoPG, s.PEP, params)
    ds.PEP =
        Glycolysis.rate_ENO(s.TwoPG, s.PEP, params) -
        Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params)
    ds.Pyruvate =
        Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) -
        Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params)
    ds.Lactate =
        Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params) -
        Glycolysis.rate_MCT(s.Lactate, s.Lactate_media, params)
    ds.Lactate_media = 0.0
    ds.ATP = (
        -0.5 * Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) -
        0.5 * Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) +
        Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) +
        Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) -
        Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) +
        Glycolysis.rate_AK(s.ATP, s.ADP, s.AMP, params)
    )
    ds.ADP = (
        0.5 * Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) +
        0.5 * Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) -
        Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) -
        Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) +
        Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) -
        2 * Glycolysis.rate_AK(s.ATP, s.ADP, s.AMP, params)
    )
    ds.AMP = Glycolysis.rate_AK(s.ATP, s.ADP, s.AMP, params)
    ds.Phosphate =
        Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) -
        Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    ds.NAD =
        Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params) -
        Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    ds.NADH =
        Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params) -
        Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params)
    ds.F26BP = 0.0
    ds.Citrate = 0.0
    ds.Phenylalanine = 0.0
end

function glycolysis_ODEs_fixed_BPG_NAD_NADH(ds, s, params, t)
    ds.Glucose_media = 0.0
    ds.Glucose =
        Glycolysis.rate_GLUT(s.Glucose_media, s.Glucose, params) -
        Glycolysis.rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params)
    ds.G6P =
        Glycolysis.rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) -
        Glycolysis.rate_GPI(s.G6P, s.F6P, params)
    ds.F6P = (
        Glycolysis.rate_GPI(s.G6P, s.F6P, params) -
        Glycolysis.rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params)
    )
    ds.F16BP = (
        Glycolysis.rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) -
        Glycolysis.rate_ALDO(s.F16BP, s.GAP, s.DHAP, params)
    )
    ds.GAP = (
        Glycolysis.rate_ALDO(s.F16BP, s.GAP, s.DHAP, params) + Glycolysis.rate_TPI(s.GAP, s.DHAP, params) -
        Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    )
    ds.DHAP =
        Glycolysis.rate_ALDO(s.F16BP, s.GAP, s.DHAP, params) - Glycolysis.rate_TPI(s.GAP, s.DHAP, params)
    ds.BPG = 0.0
    # ds.BPG =
    #     Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params) -
    #     Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params)
    ds.ThreePG =
        Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) -
        Glycolysis.rate_PGM(s.ThreePG, s.TwoPG, params)
    ds.TwoPG = Glycolysis.rate_PGM(s.ThreePG, s.TwoPG, params) - Glycolysis.rate_ENO(s.TwoPG, s.PEP, params)
    ds.PEP =
        Glycolysis.rate_ENO(s.TwoPG, s.PEP, params) -
        Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params)
    ds.Pyruvate =
        Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) -
        Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params)
    ds.Lactate =
        Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params) -
        Glycolysis.rate_MCT(s.Lactate, s.Lactate_media, params)
    ds.Lactate_media = 0.0
    ds.ATP = (
        -Glycolysis.rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) -
        Glycolysis.rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) +
        Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) +
        Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) -
        Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) +
        Glycolysis.rate_AK(s.ATP, s.ADP, s.AMP, params)
    )
    ds.ADP = (
        Glycolysis.rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) +
        Glycolysis.rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) -
        Glycolysis.rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) -
        Glycolysis.rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) +
        Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) -
        2 * Glycolysis.rate_AK(s.ATP, s.ADP, s.AMP, params)
    )
    ds.AMP = Glycolysis.rate_AK(s.ATP, s.ADP, s.AMP, params)
    ds.Phosphate =
        Glycolysis.rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) -
        Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    ds.NAD = 0.0
    # ds.NAD =
    #     Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params) -
    #     Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    ds.NADH = 0.0
    # ds.NADH =
    #     Glycolysis.rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params) -
    #     Glycolysis.rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params)
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
    max_ATPase = 1.0,
)
    tspan = (0.0, 1e8)
    pathway_Vmax = 2 * params.HK1_Vmax * params.HK1_Conc
    ATPases = 10 .^ range(log10(min_ATPase), log10(max_ATPase), n_Vmax_ATPase_values) .* pathway_Vmax
    prob = ODEProblem(ODE_func, init_conc, tspan, params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(
        ensemble_prob,
        Rodas5P(),
        EnsembleThreads(),
        trajectories = 1_000,
        abstol = 1e-15,
        reltol = 1e-5,
        save_everystep = false,
        save_start = false,
    )
    ATP_conc = [sol.u[end].ATP for sol in sim if sol.retcode == ReturnCode.Success]
    ATPase_Vmax = [ATPases[i] for (i, sol) in enumerate(sim) if sol.retcode == ReturnCode.Success] ./ (2 * params.HK1_Vmax * params.HK1_Conc)
    return (ATP_conc = ATP_conc, ATPase_Vmax = ATPase_Vmax)
end

glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
glycolysis_params.ATPase_Km_ATP = 1e-9
Complete_Model_Simulation_Data =
    find_ATP_at_ATPase_range(glycolysis_ODEs, glycolysis_params, glycolysis_init_conc_copy)

no_reg_params = deepcopy(glycolysis_params)
no_reg_params.HK1_K_a_G6P_cat = Inf
no_reg_params.HK1_K_i_G6P_reg = Inf
no_reg_params.PFKP_L = 0.0
no_reg_params.GAPDH_L = 0.0
no_reg_params.PKM2_L = 0.0
no_reg_params.ATPase_Km_ATP = 1e-9

glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
No_Reg_Model_Simulation_Data =
    find_ATP_at_ATPase_range(glycolysis_ODEs, no_reg_params, glycolysis_init_conc_copy)

glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
glycolysis_init_conc_copy.Phosphate = 1e-3
Const_Pi_No_Reg_Model_Simulation_Data =
    find_ATP_at_ATPase_range(glycolysis_ODEs_fixed_Pi, no_reg_params, glycolysis_init_conc_copy)

# glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
# glycolysis_init_conc_copy.BPG = 1e-5
# glycolysis_init_conc_copy.NADH = 1e-5
# glycolysis_init_conc_copy.NAD = 1e-3
# Const_BPG_NAD_NADH_No_Reg_Model_Simulation_Data =
#     find_ATP_at_ATPase_range(glycolysis_ODEs_fixed_BPG_NAD_NADH, no_reg_params, glycolysis_init_conc_copy)

no_reg_Keq_params = deepcopy(glycolysis_params)
no_reg_Keq_params.HK1_K_a_G6P_cat = Inf
no_reg_Keq_params.HK1_K_i_G6P_reg = Inf
no_reg_Keq_params.PFKP_L = 0.0
no_reg_Keq_params.GAPDH_L = 0.0
no_reg_Keq_params.PKM2_L = 0.0
no_reg_Keq_params.ATPase_Km_ATP = 1e-9
no_reg_Keq_params.HK1_Keq = 0.05
no_reg_Keq_params.PFKP_Keq = 0.05
no_reg_Keq_params.PGK_Keq = 5.0
no_reg_Keq_params.GAPDH_Keq = 0.5

Keq_HK_PFK_GAPDH_PGK_No_Reg_Model_Simulation_Data = find_ATP_at_ATPase_range(
    glycolysis_ODEs,
    no_reg_Keq_params,
    glycolysis_init_conc;
    min_ATPase = 0.001,
    max_ATPase = 1.0,
)


glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
Vmax_HK_PFK_eq_ATPase_No_Reg_Model_Simulation_Data =
    find_ATP_at_ATPase_range(glycolysis_ODEs_Vmax_HK_PFK_eq_ATPase, no_reg_params, glycolysis_init_conc_copy)




# Precalculate dynamic change in [Metabolite] after allostery removal
glycolysis_params_copy = deepcopy(glycolysis_params)
glycolysis_params_copy_affect2 = deepcopy(glycolysis_params)

tspan = (0.0, 60.0)
callback_time1 = 30
callback_time2 = 120
Initial_ATPase_Vmax_frac = 0.03
glycolysis_params_copy.ATPase_Vmax =
    Initial_ATPase_Vmax_frac * 2 * glycolysis_params_copy.HK1_Conc * glycolysis_params_copy.HK1_Vmax
glycolysis_params_copy_affect2.ATPase_Vmax =
    Initial_ATPase_Vmax_frac *
    2 *
    glycolysis_params_copy_affect2.HK1_Conc *
    glycolysis_params_copy_affect2.HK1_Vmax
no_reg_params.ATPase_Vmax = Initial_ATPase_Vmax_frac * 2 * no_reg_params.HK1_Conc * no_reg_params.HK1_Vmax
function affect1!(integrator)
    integrator.p .= no_reg_params
    # integrator.p.HK1_K_a_G6P_cat = Inf
    # integrator.p.HK1_K_i_G6P_reg = Inf
    # integrator.p.PFKP_L = 0.0
    # integrator.p.GAPDH_L = 0.0
    # integrator.p.PKM2_L = 0.0
    # integrator.p.ATPase_Km_ATP = 1e-9
end
function affect2!(integrator)
    integrator.p .= glycolysis_params_copy_affect2
end


PresetTime_cb1 = PresetTimeCallback(callback_time1, affect1!)
PresetTime_cb2 = PresetTimeCallback(callback_time2, affect2!)
cb_set = CallbackSet(PresetTime_cb1, PresetTime_cb2)

init_cond_prob = ODEProblem(glycolysis_ODEs, glycolysis_init_conc, (0, 1e8), glycolysis_params_copy)
init_cond_sol = solve(init_cond_prob, Rodas5P(), abstol = 1e-15, reltol = 1e-5, save_everystep = false)
new_init_cond = init_cond_sol.u[end]
prob = ODEProblem(glycolysis_ODEs, new_init_cond, tspan, glycolysis_params_copy, callback = PresetTime_cb1)
sol = solve(
    prob,
    Rodas5P(),
    abstol = 1e-15,
    reltol = 1e-5,
    saveat = [k for k = tspan[1]:((tspan[2]-tspan[1])/10_000):tspan[2]],
)

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


futile_cycle_schematic = load(
    "Glycolysis_schematic_futile_cycle.png",
)

ax_futile_cycle_schematic, im =
    image(fig[1, 1:3], rotr90(futile_cycle_schematic), axis = (aspect = DataAspect(),))
ax_futile_cycle_schematic.alignmode = Mixed(top = -20, bottom = -20, left = -25, right = -15)
# ax_futile_cycle_schematic.width = 80
ax_futile_cycle_schematic.tellwidth = false
hidedecorations!(ax_futile_cycle_schematic)
hidespines!(ax_futile_cycle_schematic)

# Plot changing in Pi, intermediate before GAPDH and after GAPDH

# Plot metabolite with switch off/on of allostery
ax_allost_on_off = Axis(
    fig[1, 4:6],
    # limits = ((0.001, 1.0), (-0.15e-3, 13e-3)),
    xlabel = "Time, min",
    ylabel = "[Metabolite],M",
    title = "Metabolite concentration changes\nin response to removal of HK and PFK allostery",
    yscale = log10,
    yticks = ([1e-8, 1e-6, 1e-4, 1e-2], ["0.01µM", "1µM", "0.1mM", "10mM"]),
    # ytickformat = ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys],
)
allost_on_off_lines = []
for metabolite in [:ATP, :ADP, :Phosphate, :G6P, :F6P, :F16BP, :GAP, :DHAP, :BPG, :TwoPG, :ThreePG, :PEP]
    if metabolite ∈ [:G6P, :F6P, :F16BP, :GAP, :DHAP]
        color = Makie.wong_colors()[1]
        linestyle = :dash
    elseif metabolite ∈ [:BPG, :TwoPG, :ThreePG, :PEP]
        color = Makie.wong_colors()[3]
        linestyle = :dash
    elseif metabolite ∈ [:ATP]
        # color = Makie.wong_colors()[1]
        color = :Black
        linestyle = nothing
    elseif metabolite ∈ [:ADP]
        # color = Makie.wong_colors()[1]
        color = :Grey
        linestyle = nothing
    elseif metabolite ∈ [:Phosphate]
        color = Makie.wong_colors()[5]
        color = :Red
        linestyle = nothing
    end
    allost_on_off_line = lines!(
        ax_allost_on_off,
        sol.t,
        sol[metabolite, :],
        label = String(metabolite),
        color = color,
        linestyle = linestyle,
    )
    push!(allost_on_off_lines, allost_on_off_line)
end
axislegend(
    ax_allost_on_off,
    [
        allost_on_off_lines[1],
        allost_on_off_lines[2],
        allost_on_off_lines[3],
        allost_on_off_lines[4],
        allost_on_off_lines[end],
    ],
    ["ATP", "ADP", "Phosphate", "Metabolites before GAPDH", "Metabolites after GAPDH"],
    position = :lb,
    # position = (1.075, 1.05),
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
    patchsize = (10, 5),
)


# Plot [ATP] maintenance with constant Phosphate
ax_ATPase_range = Axis(
    fig[2, 1:2],
    limits = ((0.001, 1.0), (-0.2e-3, 15e-3)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "[ATP],mM",
    title = "Metabolite concentration changes\nwith constant [Phosphate]",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.9], ["0.1%", "0.3%", "1%", "3%", "10%", "30%", "90%"]),
    ytickformat = ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys],
)
lines!(
    ax_ATPase_range,
    Complete_Model_Simulation_Data.ATPase_Vmax,
    Complete_Model_Simulation_Data.ATP_conc,
    color = full_color,
    linestyle = full_linestyle,
    label = "Full model",
)
lines!(
    ax_ATPase_range,
    No_Reg_Model_Simulation_Data.ATPase_Vmax,
    No_Reg_Model_Simulation_Data.ATP_conc,
    color = no_reg_color,
    linestyle = no_reg_linestyle,
    label = "Model w/o allost.",
)
lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    Const_Pi_No_Reg_Model_Simulation_Data.ATP_conc,
    color = condition_no_reg_color,
    linestyle = condition_no_reg_linestyle,
    label = "Model w/o allost.\n[Phosphate]=1mM",
)
lines!(
    ax_ATPase_range,
    [0.001, 1.0],
    repeat([glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP], 2),
    color = adenine_pool_color,
    linestyle = adenine_pool_linestyle,
)
text!(
    ax_ATPase_range,
    0.0033,
    1.02 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = adenine_pool_color,
)
axislegend(
    ax_ATPase_range,
    position = :lt,
    # position = (1.075, 1.05),
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
    patchsize = (10, 5),
)


# Plot [ATP] maintenance with different Keq HK PFK GAPDH PGK
ax_ATPase_range = Axis(
    fig[2, 3:4],
    limits = ((0.001, 1.0), (-0.2e-3, 15e-3)),
    # limits = ((0.0001, 1.0), (-0.2e-3, 14e-3)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "[ATP],mM",
    title = "Maintaining ATP concentration with lower\nKeq for HK1, PFKP, GAPDH, and PGK reactions",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.9], ["0.1%", "0.3%", "1%", "3%", "10%", "30%", "90%"]),
    ytickformat = ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys],
)
lines!(
    ax_ATPase_range,
    Complete_Model_Simulation_Data.ATPase_Vmax,
    Complete_Model_Simulation_Data.ATP_conc,
    color = full_color,
    linestyle = full_linestyle,
    label = "Full model",
)
lines!(
    ax_ATPase_range,
    No_Reg_Model_Simulation_Data.ATPase_Vmax,
    No_Reg_Model_Simulation_Data.ATP_conc,
    color = no_reg_color,
    linestyle = no_reg_linestyle,
    label = "Model w/o allost.",
)
lines!(
    ax_ATPase_range,
    Keq_HK_PFK_GAPDH_PGK_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    Keq_HK_PFK_GAPDH_PGK_No_Reg_Model_Simulation_Data.ATP_conc,
    color = condition_no_reg_color,
    linestyle = condition_no_reg_linestyle,
    label = "Model w/o allost.\nKeq HK1=Keq PFKP=0.05\nKeq GAPDH=0.5, Keq PGK=5",
)
lines!(
    ax_ATPase_range,
    [0.001, 1.0],
    repeat([glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP], 2),
    color = adenine_pool_color,
    linestyle = adenine_pool_linestyle,
)
text!(
    ax_ATPase_range,
    0.0033,
    1.02 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = adenine_pool_color,
)
axislegend(
    ax_ATPase_range,
    position = :lt,
    # position = (1.075, 1.05),
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
    patchsize = (10, 5),
)
ax_ATPase_range.width = 120


# Plot [ATP] maintenance with HK and PFK rate sub'd to ATPase rate
ax_ATPase_range = Axis(
    fig[2, 5:6],
    limits = ((0.001, 1.0), (-0.2e-3, 15e-3)),
    xlabel = "ATPase, % of pathway Vmax",
    ylabel = "[ATP],mM",
    title = "Maintaining ATP concentration\nwith 0.5V(HK1)=0.5V(PFKP)=V(ATPase)",
    xscale = log10,
    # yscale = log10,
    xticks = ([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.9], ["0.1%", "0.3%", "1%", "3%", "10%", "30%", "90%"]),
    ytickformat = ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys],
)
lines!(
    ax_ATPase_range,
    Complete_Model_Simulation_Data.ATPase_Vmax,
    Complete_Model_Simulation_Data.ATP_conc,
    color = full_color,
    linestyle = full_linestyle,
    label = "Full model",
)
lines!(
    ax_ATPase_range,
    No_Reg_Model_Simulation_Data.ATPase_Vmax,
    No_Reg_Model_Simulation_Data.ATP_conc,
    color = no_reg_color,
    linestyle = no_reg_linestyle,
    label = "Model w/o allost.",
)
lines!(
    ax_ATPase_range,
    Vmax_HK_PFK_eq_ATPase_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    Vmax_HK_PFK_eq_ATPase_No_Reg_Model_Simulation_Data.ATP_conc,
    color = condition_no_reg_color,
    linestyle = condition_no_reg_linestyle,
    label = "Model w/o allost.\n0.5V(HK1)=0.5V(PFKP)=V(ATPase)",
)
lines!(
    ax_ATPase_range,
    [0.001, 1.0],
    repeat([glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP], 2),
    color = adenine_pool_color,
    linestyle = adenine_pool_linestyle,
)
text!(
    ax_ATPase_range,
    0.0033,
    1.02 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = adenine_pool_color,
)
axislegend(
    ax_ATPase_range,
    # position = :lt,
    position = (0.5, 1.0),
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
    patchsize = (10, 5),
)

colgap!(fig.layout, 0)
rowgap!(fig.layout, 10)


label_a = fig[1, 1, TopLeft()] = Label(fig, "A", fontsize = 12, halign = :right, padding = (0, 5, 5, 0))
label_b = fig[1, 4, TopLeft()] = Label(fig, "B", fontsize = 12, halign = :right, padding = (0, 5, 5, 0))
label_c = fig[2, 1, TopLeft()] = Label(fig, "C", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_d = fig[2, 3, TopLeft()] = Label(fig, "D", fontsize = 12, halign = :right, padding = (0, 0, 5, 0))
label_e = fig[2, 5, TopLeft()] = Label(fig, "E", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))

fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_Fig5_mechanism_of_ATP_maintenence_by_allost.png", fig, px_per_unit = 4)

