using Glycolysis
using DifferentialEquations, ProgressMeter
using CairoMakie, Dates, Printf, Statistics, StatsBase
using DataFrames, CSV

#=
Running this code might take up to 20 min on 8 core computer.
Most of the time spent running simulations at constant phosphate.
You can decrease the running time by setting n_Vmax_ATPase_values keyword argument of
find_ATP_at_ATPase_range() to a smaller number thatn default 1000 value.
The downside of this is that the resulting plots will be less smooth as some of the simulations
with constant phosphate fail.
=#

#=
Global sensitivity analysis (GSA) is a computationally intensive task that
needs to performed on a computing cluster with many cores. Here we load the results
of such a computation for plotting. The code to run this analysis on the cluster is
available in gsa_cluster_code/ folder.
=#
#Download sens indexes for init cond and pool sizes
GSA_ATP_AUC_cofactor_pools_sens_ind = CSV.read(
    "gsa_cluster_code/041124_ATP_AUC_gsa_sobol_cofactor_pool_10000_3x_range.csv",
    DataFrame,
)
hist_ATP_AUC_cofactor_pools =
    CSV.read("gsa_cluster_code/041124_hist_ATP_AUC_10000_runs_3x_pool_sizes.csv", DataFrame)
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
# Precalculate and save outputs for model with different Phosphate pool and Adenine pool sizes

function glycolysis_ODEs_fixed_Pi(ds, s, params, t)
    ds.Glucose_media = 0.0
    ds.Glucose = Glycolysis.rate_GLUT(s, params) - Glycolysis.rate_HK1(s, params)
    ds.G6P = Glycolysis.rate_HK1(s, params) - Glycolysis.rate_GPI(s, params)
    ds.F6P = Glycolysis.rate_GPI(s, params) - Glycolysis.rate_PFKP(s, params)
    ds.F16BP = Glycolysis.rate_PFKP(s, params) - Glycolysis.rate_ALDO(s, params)
    ds.GAP = (
        Glycolysis.rate_ALDO(s, params) + Glycolysis.rate_TPI(s, params) -
        Glycolysis.rate_GAPDH(s, params)
    )
    ds.DHAP = Glycolysis.rate_ALDO(s, params) - Glycolysis.rate_TPI(s, params)
    ds.BPG = Glycolysis.rate_GAPDH(s, params) - Glycolysis.rate_PGK(s, params)
    ds.ThreePG = Glycolysis.rate_PGK(s, params) - Glycolysis.rate_PGM(s, params)
    ds.TwoPG = Glycolysis.rate_PGM(s, params) - Glycolysis.rate_ENO(s, params)
    ds.PEP = Glycolysis.rate_ENO(s, params) - Glycolysis.rate_PKM2(s, params)
    ds.Pyruvate = Glycolysis.rate_PKM2(s, params) - Glycolysis.rate_LDH(s, params)
    ds.Lactate = Glycolysis.rate_LDH(s, params) - Glycolysis.rate_MCT(s, params)
    ds.Lactate_media = 0.0
    ds.ATP = (
        -Glycolysis.rate_HK1(s, params) - Glycolysis.rate_PFKP(s, params) +
        Glycolysis.rate_PGK(s, params) +
        Glycolysis.rate_PKM2(s, params) - Glycolysis.rate_ATPase(s, params) +
        Glycolysis.rate_AK(s, params)
    )
    ds.ADP = (
        Glycolysis.rate_HK1(s, params) + Glycolysis.rate_PFKP(s, params) -
        Glycolysis.rate_PGK(s, params) - Glycolysis.rate_PKM2(s, params) +
        Glycolysis.rate_ATPase(s, params) - 2 * Glycolysis.rate_AK(s, params)
    )
    ds.AMP = Glycolysis.rate_AK(s, params)
    ds.Phosphate = 0.0
    ds.NAD = Glycolysis.rate_LDH(s, params) - Glycolysis.rate_GAPDH(s, params)
    ds.NADH = Glycolysis.rate_GAPDH(s, params) - Glycolysis.rate_LDH(s, params)
    ds.F26BP = 0.0
    ds.Citrate = 0.0
    ds.Phenylalanine = 0.0
end

function find_ATP_at_ATPase_range(
    ODE_func,
    params,
    init_conc;
    n_Vmax_ATPase_values = 1000,
    min_ATPase = 0.003,
    max_ATPase = 0.3,
)
    tspan = (0.0, 1e8)
    pathway_Vmax = 2 * params.HK1_Vmax * params.HK1_Conc
    ATPases =
        10 .^ range(log10(min_ATPase), log10(max_ATPase), n_Vmax_ATPase_values) .*
        pathway_Vmax
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
        trajectories = n_Vmax_ATPase_values,
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
        save_start = false,
    )
    ATP_conc = [sol.u[end].ATP for sol in sim if sol.retcode == ReturnCode.Success]

    ATPase_Vmax =
        [ATPases[i] for (i, sol) in enumerate(sim) if sol.retcode == ReturnCode.Success] ./ (2 * params.HK1_Vmax * params.HK1_Conc)

    return (ATP_conc = ATP_conc, ATPase_Vmax = ATPase_Vmax)
end

no_nonadenine_pi_pool_glycolysis_init_conc = deepcopy(glycolysis_init_conc)
no_nonadenine_pi_pool_glycolysis_init_conc.G6P = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.F6P = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.F16BP = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.GAP = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.DHAP = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.BPG = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.ThreePG = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.TwoPG = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.PEP = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.Phosphate = 0.0
Complete_Model_Simulation_Data = find_ATP_at_ATPase_range(
    glycolysis_ODEs,
    glycolysis_params,
    no_nonadenine_pi_pool_glycolysis_init_conc,
)



#MARK: Precalculate results for model with different adenine pool sizes
init_cond_list_Adenine = []
model_names_list_Adenine = []
# ATP_multiplier_normalizer = [4.0, 2.0, 1.0, 0.5, 0.25]
ATP_multiplier_normalizer = [2.0, 1.0, 0.5]
Adenine_pool_size =
    glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP
for multiplier in ATP_multiplier_normalizer
    glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
    # glycolysis_init_conc_copy.Pyruvate = multiplier * glycolysis_init_conc.Pyruvate
    glycolysis_init_conc_copy.ATP = multiplier * glycolysis_init_conc.ATP
    glycolysis_init_conc_copy.ADP = multiplier * glycolysis_init_conc.ADP
    glycolysis_init_conc_copy.AMP = multiplier * glycolysis_init_conc.AMP

    push!(init_cond_list_Adenine, glycolysis_init_conc_copy)
    push!(
        model_names_list_Adenine,
        "$(Int(round(multiplier * Adenine_pool_size * 1000)))mM",
    )
end
Simulation_Data_Adenine = []
for model_glycolysis_init_conc in init_cond_list_Adenine
    res = find_ATP_at_ATPase_range(
        glycolysis_ODEs,
        glycolysis_params,
        model_glycolysis_init_conc,
        min_ATPase = 0.003,
        max_ATPase = 0.6,
    )
    push!(Simulation_Data_Adenine, res)
end




#Precalculate results for model with fixed Pi
init_cond_list_Pi = []
model_names_list_Pi = []
for multiplier in [2.0, 1.0, 0.5]
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
    # glycolysis_init_conc_copy.Phosphate = 10 * glycolysis_init_conc.Phosphate
    glycolysis_init_conc_copy.Phosphate = 10e-3
    push!(init_cond_list_Pi, glycolysis_init_conc_copy)
    push!(model_names_list_Pi, "$(multiplier)x")
end
Simulation_Data_Pi = []
for model_glycolysis_init_conc in init_cond_list_Pi
    res = find_ATP_at_ATPase_range(
        glycolysis_ODEs_fixed_Pi,
        glycolysis_params,
        model_glycolysis_init_conc,
        # Uncomment n_Vmax_ATPase_values and set it to a smaller number (e.g. 100) to decrease running time
        # n_Vmax_ATPase_values = 300,
    )
    push!(Simulation_Data_Pi, res)
end

#Precalculate results for model without CK, AK, NDPK
glycolysis_params_no_AK_NDPK_CK = deepcopy(glycolysis_params)
glycolysis_params_no_AK_NDPK_CK.AK_Vmax = 0.0
glycolysis_params_no_AK_NDPK_CK.NDPK_Vmax = 0.0
glycolysis_params_no_AK_NDPK_CK.CK_Vmax = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc = deepcopy(glycolysis_init_conc)
no_nonadenine_pi_pool_glycolysis_init_conc.G6P = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.F6P = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.F16BP = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.GAP = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.DHAP = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.BPG = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.ThreePG = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.TwoPG = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.PEP = 0.0
no_nonadenine_pi_pool_glycolysis_init_conc.Phosphate = 0.0
no_AK_NDPK_CK_Model_Simulation_Data = find_ATP_at_ATPase_range(
    glycolysis_ODEs,
    glycolysis_params_no_AK_NDPK_CK,
    no_nonadenine_pi_pool_glycolysis_init_conc,
)

#Precalculate results for model with CK with different Creatine Phosphate pool sizes and no non-adenine phosphate pool
Init_cond_creatine_pool_concs_list = []
Creatine_pool_concs_list = [0, 0.2e-3, 2e-3, 20e-3]
multiplier = 0.0
for creatine_pool in Creatine_pool_concs_list
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
    glycolysis_init_conc_copy.Creatine = creatine_pool * 0.1
    glycolysis_init_conc_copy.Phosphocreatine = creatine_pool * 0.9
    push!(Init_cond_creatine_pool_concs_list, glycolysis_init_conc_copy)
end
Simulation_Data_CK = []
Init_cond_creatine_pool_concs_list
glycolysis_params_copy = deepcopy(glycolysis_params)
glycolysis_params_copy.AK_Vmax = 0.0
glycolysis_params_copy.NDPK_Vmax = 0.0
glycolysis_params_copy.CK_Vmax = 0.03
for model_init_cond in Init_cond_creatine_pool_concs_list
    res = find_ATP_at_ATPase_range(glycolysis_ODEs, glycolysis_params_copy, model_init_cond)
    push!(Simulation_Data_CK, res)
end


#Precalculate results for model with NDPK with different Creatine Phosphate pool sizes and no non-adenine phosphate pool
Init_cond_NTP_pool_concs_list = []
NTP_pool_concs_list = [0, 0.2e-3, 2e-3, 20e-3]
multiplier = 0.0
for NTP_pool in NTP_pool_concs_list
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
    glycolysis_init_conc_copy.NTP = NTP_pool * 0.9
    glycolysis_init_conc_copy.NDP = NTP_pool * 0.1
    push!(Init_cond_NTP_pool_concs_list, glycolysis_init_conc_copy)
end
Simulation_Data_NDPK = []
glycolysis_params_copy = deepcopy(glycolysis_params)
glycolysis_params_copy.AK_Vmax = 0.0
glycolysis_params_copy.NDPK_Vmax = 0.03
glycolysis_params_copy.CK_Vmax = 0.0
for model_init_cond in Init_cond_NTP_pool_concs_list
    res = find_ATP_at_ATPase_range(glycolysis_ODEs, glycolysis_params_copy, model_init_cond)
    push!(Simulation_Data_NDPK, res)
end



#Precalculate dynamic changes of ATP +/- AK, NDPK, CK
glycolysis_params_copy = deepcopy(glycolysis_params)
glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
glycolysis_init_conc_copy.ATP +
glycolysis_init_conc_copy.ADP +
glycolysis_init_conc_copy.AMP
glycolysis_init_conc_copy.ATP += glycolysis_init_conc_copy.AMP
glycolysis_init_conc_copy.AMP = 0.0

tspan = (0.0, 240.0 / 4)
Initial_ATPase_Vmax_frac = 0.03
High_ATPase_Vmax_frac = Initial_ATPase_Vmax_frac * 2
Low_ATPase_Vmax_frac = Initial_ATPase_Vmax_frac / 2
ATPase_change_time1 = 60 / 4
ATPase_change_time2 = 120 / 4
ATPase_change_time3 = 180 / 4
glycolysis_params_copy.ATPase_Vmax =
    Initial_ATPase_Vmax_frac *
    2 *
    glycolysis_params_copy.HK1_Conc *
    glycolysis_params_copy.HK1_Vmax

function affect1!(integrator)
    integrator.p.ATPase_Vmax =
        High_ATPase_Vmax_frac *
        2 *
        glycolysis_params_copy.HK1_Conc *
        glycolysis_params_copy.HK1_Vmax
end
function affect2!(integrator)
    integrator.p.ATPase_Vmax =
        Low_ATPase_Vmax_frac *
        2 *
        glycolysis_params_copy.HK1_Conc *
        glycolysis_params_copy.HK1_Vmax
end
function affect3!(integrator)
    integrator.p.ATPase_Vmax =
        Initial_ATPase_Vmax_frac *
        2 *
        glycolysis_params_copy.HK1_Conc *
        glycolysis_params_copy.HK1_Vmax
end

PresetTime_cb1 = PresetTimeCallback(ATPase_change_time1, affect1!)
PresetTime_cb2 = PresetTimeCallback(ATPase_change_time2, affect2!)
PresetTime_cb3 = PresetTimeCallback(ATPase_change_time3, affect3!)
cb_set = CallbackSet(PresetTime_cb1, PresetTime_cb2, PresetTime_cb3)

Dynamic_ATP_w_wo_AK_NDPK_CK = []
Dynamic_ATPprod_w_wo_AK_NDPK_CK = []
Dynamic_ATPenergy_w_wo_AK_NDPK_CK = []
Dynamic_ATPase_w_wo_AK_NDPK_CK = []
Dynamic_timepoints_w_wo_AK_NDPK_CK = []
enzyme_combos_AK_CK_NDPK = ["+ AK, CK, NDPK", "- AK, CK, NDPK"]#, "no AK", "no CK", "no NDPK"]
for (i, enzyme) in enumerate(enzyme_combos_AK_CK_NDPK)
    if enzyme == "- AK, CK, NDPK"
        glycolysis_params_copy.AK_Vmax = 0.0
        glycolysis_params_copy.CK_Vmax = 0.0
        glycolysis_params_copy.NDPK_Vmax = 0.0
    elseif enzyme == "+ AK, CK, NDPK"
        glycolysis_params_copy.AK_Vmax = 0.03
        glycolysis_params_copy.CK_Vmax = 0.03
        glycolysis_params_copy.NDPK_Vmax = 0.03
    elseif enzyme == "no AK"
        glycolysis_params_copy.AK_Vmax = 0.0
        glycolysis_params_copy.CK_Vmax = 0.03
        glycolysis_params_copy.NDPK_Vmax = 0.03
    elseif enzyme == "no CK"
        glycolysis_params_copy.AK_Vmax = 0.03
        glycolysis_params_copy.CK_Vmax = 0.0
        glycolysis_params_copy.NDPK_Vmax = 0.03
    elseif enzyme == "no NDPK"
        glycolysis_params_copy.AK_Vmax = 0.03
        glycolysis_params_copy.CK_Vmax = 0.03
        glycolysis_params_copy.NDPK_Vmax = 0.0
    end
    println(
        "AK = $(glycolysis_params_copy.AK_Vmax), CK = $(glycolysis_params_copy.CK_Vmax), NDPK = $(glycolysis_params_copy.NDPK_Vmax)",
    )

    init_cond_prob = ODEProblem(
        glycolysis_ODEs,
        glycolysis_init_conc_copy,
        (0, 1e8),
        glycolysis_params_copy,
    )
    init_cond_sol = solve(
        init_cond_prob,
        Rodas5P(),
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
    )
    new_init_cond = init_cond_sol.u[end]
    prob = ODEProblem(
        glycolysis_ODEs,
        new_init_cond,
        tspan,
        glycolysis_params_copy,
        callback = cb_set,
    )
    sol = solve(
        prob,
        Rodas5P(),
        abstol = 1e-15,
        reltol = 1e-8,
        saveat = [k for k = tspan[1]:((tspan[2]-tspan[1])/10_000):tspan[2]],
    )

    push!(
        Dynamic_ATPprod_w_wo_AK_NDPK_CK,
        [Glycolysis.conc_to_rates(conc, glycolysis_params_copy).ATPprod for conc in sol.u] / (2 * glycolysis_params_copy.HK1_Conc * glycolysis_params_copy.HK1_Vmax),
    )
    push!(
        Dynamic_ATPenergy_w_wo_AK_NDPK_CK,
        [
            Glycolysis.conc_to_disequilibrium_ratios(
                conc,
                glycolysis_params_copy,
            ).Q_Keq_ATPase for conc in sol.u
        ],
    )
    ATPase = Initial_ATPase_Vmax_frac * ones(length(sol))
    ATPase[sol.t.>ATPase_change_time1.&&sol.t.<ATPase_change_time2] .= High_ATPase_Vmax_frac
    ATPase[sol.t.>ATPase_change_time2.&&sol.t.<ATPase_change_time3] .= Low_ATPase_Vmax_frac
    push!(Dynamic_ATPase_w_wo_AK_NDPK_CK, ATPase)
    push!(Dynamic_ATP_w_wo_AK_NDPK_CK, [conc.ATP for conc in sol.u])
    push!(Dynamic_timepoints_w_wo_AK_NDPK_CK, sol.t)
end


##
# Plot results

size_inches = (6.5, 7)
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


#Plot hist ans sens indexes for cofactor pools
#Plot [ATP] AUC histogram data
stacked_A_B = fig[1, 1] = GridLayout()
ax_ATP_AUC_hist = Axis(
    # fig[1, 1][1,1],
    stacked_A_B[1, 1],
    limits = ((nothing), (0, nothing)),
    title = "Variance of ∑[ATP] at a 9x\nrange of cofactor pool sizes",
    xlabel = "∑[ATP], normalized to adenine pool size",
    ylabel = "Fraction of counts",
    xtickformat = xs -> ["$(round(x,sigdigits=1))" for x in xs],
)

hist_all_glycolysis_params = hist!(
    ax_ATP_AUC_hist,
    hist_ATP_AUC_cofactor_pools.all_params,
    color = (Makie.wong_colors()[3], 0.5),
    normalization = :probability,
    bins = range(0, 1, 50),
)
text!(
    ax_ATP_AUC_hist,
    0.0,
    0.1;
    text = "CV=$(round(std(hist_ATP_AUC_cofactor_pools.all_params) / mean(hist_ATP_AUC_cofactor_pools.all_params), sigdigits=2))",
    fontsize = 8,
)
vlines!(ax_ATP_AUC_hist, Complete_model_mean_ATP, color = :red, linestyle = :dash)
text!(
    ax_ATP_AUC_hist,
    Complete_model_mean_ATP * 1.03,
    0.02;
    text = "Complete\nModel",
    color = :red,
    fontsize = 8,
    # rotation = π / 2,
)
# Plot sensitivity indexes
barplot_df = GSA_ATP_AUC_cofactor_pools_sens_ind
barplot_df = stack(barplot_df; variable_name = :Parameters, value_name = :sens_ind_value)
barplot_df = rename(barplot_df, "Row" => "sens_ind")
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

# replace.(String.(glycolysis_params_group_names), " " => "\n", count=1)
ax_Sens_Ind = Axis(
    # fig[1, 1][2,1],
    stacked_A_B[2, 1],
    limits = ((nothing), (0, 1.0)),
    title = "Sens. indxs for ∑[ATP] at a 9x\nrange of cofactor pool sizes",
    ylabel = "Sensitivity Index",
    xticks = (
        1:length(unique(barplot_df.Parameters)),
        replace.(
            unique(barplot_df.Parameters),
            "Adenine_pool_scaler" => "Adenine\npool",
            "NAD_pool_scaler" => "NAD(H)\npool",
            "Pi_scaler" => "Phosphate\npool",
        ),
    ),
    # xticklabelrotation = π / 2,
)
colors = Makie.wong_colors()
barplot!(
    ax_Sens_Ind,
    barplot_df.group_index,
    barplot_df.sens_ind_value,
    dodge = barplot_df.dodge_index,
    color = colors[barplot_df.dodge_index],
)
labels = ["S1", "ST"]
elements = [PolyElement(polycolor = colors[i]) for i = 1:length(labels)]
axislegend(
    ax_Sens_Ind,
    elements,
    labels,
    # [hist_all_glycolysis_params, hist_fixed_HK_Km_Vmax],
    # ["All glycolysis_params varied", "All glycolysis_params varied\nexcept HK1 Km, Vmax"],
    position = :lt,
    rowgap = 1,
    colgap = 7,
    patchlabelgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -5),
    patchsize = (10, 5),
    orientation = :horizontal,
)

#Plot Pi pool variation
ax_ATP_conc = Axis(
    fig[1, 2],
    limits = (
        (0.003, 0.9),
        (
            0.0,
            1.5 * (
                glycolysis_init_conc.ATP +
                glycolysis_init_conc.ADP +
                glycolysis_init_conc.AMP
            ),
        ),
    ),
    xscale = log10,
    xlabel = rich("ATPase, % of pathway V", subscript("max")),
    ylabel = "[ATP],mM",
    title = "Effect of Non-Adenine Phosphate Pool\nchanges at constant [Phosphate] = 10mM",
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9]),
    xtickformat = xs -> [
        x * 100 >= 1 ? "$(Int(round(x*100, sigdigits=1)))%" :
        "$(round(x*100, sigdigits=1))%" for x in xs
    ],
    ytickformat = ys -> ["$(Int(round(y*1000, sigdigits=3)))" for y in ys],
)

color = Makie.wong_colors()
linewidth = 5
for (i, model_res) in enumerate(Simulation_Data_Pi)
    lines!(
        ax_ATP_conc,
        model_res.ATPase_Vmax,
        model_res.ATP_conc,
        color = color[i],
        label = model_names_list_Pi[i],
        linewidth = linewidth,
        # label = replace(
        #     model_names_list_NADpool[i],
        #     "$(Int(round(Pi_pool*1000)))mM" => "$(Int(round(Pi_pool*1000)))mM\nDefault Model",
        # ),
    )
    linewidth *= 0.5
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
    1.03 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = :grey,
)
axislegend(
    ax_ATP_conc,
    ax_ATP_conc,
    "Non-adenine phosphate pool\nrelative to 20mM in full model",
    titlegap = 2,
    colgap = 3,
    patchlabelgap = 1,
    position = (0.5, 1.1),
    patchsize = (7.5, 2.5),
    orientation = :horizontal,
    framevisible = false,
)

#Plot model with more or less ATP pool
ax_ATP_conc = Axis(
    fig[1, 3],
    limits = ((0.003, 0.9), (0.0, 1.5)),
    xscale = log10,
    xlabel = rich("ATPase, % of pathway V", subscript("max")),
    ylabel = "[ATP], normalized by adenine pool",
    title = "Effect of adenine pool size\non maintenance of ATP",
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9]),
    xtickformat = xs -> [
        x * 100 >= 1 ? "$(Int(round(x*100, sigdigits=1)))%" :
        "$(round(x*100, sigdigits=1))%" for x in xs
    ],
    # ytickformat = ys -> ["$(round(y*1000, sigdigits=3))" for y in ys],
)

color = Makie.wong_colors()

for (i, model_res) in enumerate(Simulation_Data_Adenine)
    supported_ATPase_range =
        model_res.ATPase_Vmax[model_res.ATP_conc./ATP_multiplier_normalizer[i].>0.5*(glycolysis_init_conc.ATP+glycolysis_init_conc.ADP+glycolysis_init_conc.AMP)]
    supported_ATPase_fold_range =
        round(supported_ATPase_range[end] / supported_ATPase_range[1], sigdigits = 2)
    lines!(
        ax_ATP_conc,
        model_res.ATPase_Vmax,
        model_res.ATP_conc ./ (
            ATP_multiplier_normalizer[i] * (
                glycolysis_init_conc.ATP +
                glycolysis_init_conc.ADP +
                glycolysis_init_conc.AMP
            )
        ),
        color = model_names_list_Adenine[i] != "9mM" ? color[i] : :Black,
        # label = replace(
        #     model_names_list_Adenine[i],
        #     "$(Int(round(1.0 * Adenine_pool_size * 1000)))mM" => "$(Int(round(1.0 * Adenine_pool_size * 1000)))mM\n($(supported_ATPase_fold_range)-fold)\nDefault Model",
        #     "mM" => "mM\n($(supported_ATPase_fold_range)-fold)",
        # ),
        label = "$(ATP_multiplier_normalizer[i])x",
    )
end
lines!(ax_ATP_conc, [0.001, 1.0], repeat([1.0], 2), color = :grey, linestyle = :dash)
text!(
    ax_ATP_conc,
    0.004,
    1.01,
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = :grey,
)
axislegend(
    ax_ATP_conc,
    ax_ATP_conc,
    "Adenine pool size\nrelative to 9mM in default model",
    titlegap = 2,
    colgap = 3,
    patchlabelgap = 1,
    position = (0.5, 1.1),
    patchsize = (7.5, 2.5),
    orientation = :horizontal,
    framevisible = false,
)




#Plot dynamic response to ATPase change +/- AK, NPDK, CK
line_ATPase_width = 2
line_ATPase_color = :grey

# Plot ATP production
ax_ATP_prod = Axis(
    fig[2, 1],
    limits = (nothing, (0.0, 1.5 * maximum(Dynamic_ATPase_w_wo_AK_NDPK_CK[1]))),
    xlabel = "Time, min",
    ylabel = rich("ATP production rate, relative to glycolysis V", subscript("max")),
    title = "Dynamically matching\nATP supply and demand",
)
ax_ATPase = Axis(
    fig[2, 1],
    limits = (nothing, (0.0, 1.5 * maximum(Dynamic_ATPase_w_wo_AK_NDPK_CK[1]))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)

for (i, enzyme_name) in enumerate(enzyme_combos_AK_CK_NDPK)
    lines!(
        ax_ATP_prod,
        Dynamic_timepoints_w_wo_AK_NDPK_CK[i],
        Dynamic_ATPprod_w_wo_AK_NDPK_CK[i],
        label = enzyme_name,
        color = enzyme_name != "+ AK, CK, NDPK" ? :Red : :Black,
        linestyle = enzyme_name != "+ AK, CK, NDPK" ? :dash : nothing,
    )
end

ATPase_line = lines!(
    ax_ATPase,
    Dynamic_timepoints_w_wo_AK_NDPK_CK[1],
    Dynamic_ATPase_w_wo_AK_NDPK_CK[1],
    linestyle = :dot,
    color = line_ATPase_color,
    linewidth = line_ATPase_width,
)

axislegend(
    ax_ATP_prod,
    # [ATP_prod_line, no_reg_ATP_prod_line, ATPase_line],
    # ["Full Model", "Model w/o allost.", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (0, -4, 0, -4),
    patchsize = (10, 5),
)

# Plot [ATP]
ax_ATP_conc = Axis(
    fig[2, 2],
    limits = (nothing, (8.665e-3, 8.685e-3)),
    xlabel = "Time, min",
    ylabel = "[ATP], mM",
    title = "Dynamically maintaining\nATP concentration",
    yscale = log10,
    yticks = ([8.67e-3, 8.68e-3]),
    ytickformat = ys -> ["$(round(y*1000, sigdigits = 3))" for y in ys],
)
ax_ATPase = Axis(
    fig[2, 2],
    limits = (nothing, (0.0, 1.5 * maximum(Dynamic_ATPase_w_wo_AK_NDPK_CK[1]))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)

for (i, enzyme_name) in enumerate(enzyme_combos_AK_CK_NDPK)
    lines!(
        ax_ATP_conc,
        Dynamic_timepoints_w_wo_AK_NDPK_CK[i],
        Dynamic_ATP_w_wo_AK_NDPK_CK[i],
        color = enzyme_name != "+ AK, CK, NDPK" ? :Red : :Black,
        linestyle = enzyme_name != "+ AK, CK, NDPK" ? :dash : nothing,
        label = enzyme_name,
    )
end

ATPase_line = lines!(
    ax_ATPase,
    Dynamic_timepoints_w_wo_AK_NDPK_CK[1],
    Dynamic_ATPase_w_wo_AK_NDPK_CK[1],
    linestyle = :dot,
    color = line_ATPase_color,
    linewidth = line_ATPase_width,
)

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
    1.05 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = :grey,
)
axislegend(
    ax_ATP_conc,
    # [ATP_line, no_reg_ATP_line, ATPase_line],
    # ["Full Model", "Model w/o allost.", "ATPase rate"],
    # position = :rt,
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (0, -4, 0, -4),
    patchsize = (10, 5),
)

# Plot energy of ATPase reaction
ax_ATP_energy = Axis(
    fig[2, 3],
    limits = (nothing, (0, 31)),
    xlabel = "Time, min",
    ylabel = rich("Energy of ATPase reaction, k", subscript("B"), "T"),
    title = "Dynamically maintaining\nATP energy",
    ytickformat = ys -> ["$(Int(round(y)))" for y in ys],
)
ax_ATPase = Axis(
    fig[2, 3],
    limits = (nothing, (0.0, 1.5 * maximum(Dynamic_ATPase_w_wo_AK_NDPK_CK[1]))),
    ylabel = rich("ATPase rate, relative to glycolysis V", subscript("max")),
    yaxisposition = :right,
    ygridvisible = false,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)

for (i, enzyme_name) in enumerate(enzyme_combos_AK_CK_NDPK)
    lines!(
        ax_ATP_energy,
        Dynamic_timepoints_w_wo_AK_NDPK_CK[i],
        -log.(Dynamic_ATPenergy_w_wo_AK_NDPK_CK[i]),
        color = enzyme_name != "+ AK, CK, NDPK" ? :Red : :Black,
        linestyle = enzyme_name != "+ AK, CK, NDPK" ? :dash : nothing,
        label = enzyme_name,
    )
end

ATPase_energy_line = lines!(
    ax_ATPase,
    Dynamic_timepoints_w_wo_AK_NDPK_CK[1],
    Dynamic_ATPase_w_wo_AK_NDPK_CK[1],
    linestyle = :dot,
    color = line_ATPase_color,
    linewidth = line_ATPase_width,
)

axislegend(
    ax_ATP_energy,
    # [ATPase_energy_line, no_reg_ATPase_energy_line, ATPase_line],
    # ["Full Model", "Model w/o allost.", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (0, -4, 0, -4),
    patchsize = (10, 5),
)

#Plot model with more or less phosphate pool and +/- AK, NPDK, CK
ax_ATP_conc = Axis(
    fig[3, 1],
    limits = (
        (0.003, 0.9),
        (
            0.0,
            1.5 * (
                glycolysis_init_conc.ATP +
                glycolysis_init_conc.ADP +
                glycolysis_init_conc.AMP
            ),
        ),
    ),
    xscale = log10,
    title = "Effect of AK\non Maintenance of ATP\nwithout non-adenine phosphate pool",
    xlabel = rich("ATPase, % of pathway V", subscript("max")),
    ylabel = "[ATP],mM",
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9]),
    xtickformat = xs -> [
        x * 100 >= 1 ? "$(Int(round(x*100, sigdigits=1)))%" :
        "$(round(x*100, sigdigits=1))%" for x in xs
    ],
    ytickformat = ys -> ["$(round(y*1000, sigdigits=3))" for y in ys],
)
color = Makie.wong_colors()

lines!(
    ax_ATP_conc,
    no_AK_NDPK_CK_Model_Simulation_Data.ATPase_Vmax,
    no_AK_NDPK_CK_Model_Simulation_Data.ATP_conc,
    # color=enzyme_name != "AK" ? color[i] : :Black,
    label = "without AK",
)
lines!(
    ax_ATP_conc,
    Complete_Model_Simulation_Data.ATPase_Vmax,
    Complete_Model_Simulation_Data.ATP_conc,
    # color=enzyme_name != "AK" ? color[i] : :Black,
    label = "with AK",
)

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
    # "Included Enzymes",
    titlegap = 2,
    colgap = 3,
    patchlabelgap = 1,
    position = (0.5, 1.1),
    patchsize = (7.5, 2.5),
    orientation = :horizontal,
    framevisible = false,
)


#Plot model with different Creatine pool with CK
ax_ATP_conc = Axis(
    fig[3, 2],
    limits = (
        (0.003, 0.9),
        (
            0.0,
            1.5 * (
                glycolysis_init_conc.ATP +
                glycolysis_init_conc.ADP +
                glycolysis_init_conc.AMP
            ),
        ),
    ),
    xscale = log10,
    title = "Effect of Creatine pool size and CK\non Maintenance of ATP\nwithout non-adenine phosphate pool",
    xlabel = rich("ATPase, % of pathway V", subscript("max")),
    ylabel = "[ATP],mM",
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9]),
    xtickformat = xs -> [
        x * 100 >= 1 ? "$(Int(round(x*100, sigdigits=1)))%" :
        "$(round(x*100, sigdigits=1))%" for x in xs
    ],
    ytickformat = ys -> ["$(round(y*1000, sigdigits=3))" for y in ys],
)
color = Makie.wong_colors()
for (i, model_res) in enumerate(Simulation_Data_CK)
    supported_ATPase_range =
        model_res.ATPase_Vmax[model_res.ATP_conc.>0.5*(glycolysis_init_conc.ATP+glycolysis_init_conc.ADP+glycolysis_init_conc.AMP)]
    supported_ATPase_fold_range =
        round(supported_ATPase_range[end] / supported_ATPase_range[1], sigdigits = 2)
    lines!(
        ax_ATP_conc,
        model_res.ATPase_Vmax,
        model_res.ATP_conc,
        label = string(1000 * Creatine_pool_concs_list[i]),
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
    "Creatine pool size, mM",
    titlegap = 2,
    colgap = 3,
    patchlabelgap = 1,
    position = (0.5, 1.1),
    patchsize = (7.5, 2.5),
    orientation = :horizontal,
    framevisible = false,
)

#Plot model with more or less phosphate pool and +/- AK, NPDK, CK
ax_ATP_conc = Axis(
    fig[3, 3],
    limits = (
        (0.003, 0.9),
        (
            0.0,
            1.5 * (
                glycolysis_init_conc.ATP +
                glycolysis_init_conc.ADP +
                glycolysis_init_conc.AMP
            ),
        ),
    ),
    xscale = log10,
    title = "Effect of NTP pool size and NDPK\non Maintenance of ATP\nwithout non-adenine phosphate pool",
    xlabel = rich("ATPase, % of pathway V", subscript("max")),
    ylabel = "[ATP],mM",
    xticks = ([0.003, 0.01, 0.03, 0.1, 0.3, 0.9]),
    xtickformat = xs -> [
        x * 100 >= 1 ? "$(Int(round(x*100, sigdigits=1)))%" :
        "$(round(x*100, sigdigits=1))%" for x in xs
    ],
    ytickformat = ys -> ["$(round(y*1000, sigdigits=3))" for y in ys],
)
color = Makie.wong_colors()
for (i, model_res) in enumerate(Simulation_Data_NDPK)
    supported_ATPase_range =
        model_res.ATPase_Vmax[model_res.ATP_conc.>0.5*(glycolysis_init_conc.ATP+glycolysis_init_conc.ADP+glycolysis_init_conc.AMP)]
    supported_ATPase_fold_range =
        round(supported_ATPase_range[end] / supported_ATPase_range[1], sigdigits = 2)
    lines!(
        ax_ATP_conc,
        model_res.ATPase_Vmax,
        model_res.ATP_conc,
        label = string(1000 * NTP_pool_concs_list[i]),
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
    "NTP pool size, mM",
    titlegap = 2,
    colgap = 3,
    patchlabelgap = 1,
    position = (0.5, 1.1),
    patchsize = (7.5, 2.5),
    orientation = :horizontal,
    framevisible = false,
)
rowgap!(stacked_A_B, 5)
colgap!(fig.layout, 5)
rowgap!(fig.layout, 5)

label_a =
    fig[1, 1][1, 1, TopLeft()] =
        Label(fig, "A", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_b =
    fig[1, 1][2, 1, TopLeft()] =
        Label(fig, "B", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_c =
    fig[1, 2, TopLeft()] =
        Label(fig, "C", fontsize = 12, halign = :right, padding = (0, 15, 5, 0))
label_d =
    fig[1, 3, TopLeft()] =
        Label(fig, "D", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_e =
    fig[2, 1, TopLeft()] =
        Label(fig, "E", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_f =
    fig[2, 2, TopLeft()] =
        Label(fig, "F", fontsize = 12, halign = :right, padding = (0, 15, 5, 0))
label_g =
    fig[2, 3, TopLeft()] =
        Label(fig, "G", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_h =
    fig[3, 1, TopLeft()] =
        Label(fig, "H", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))
label_i =
    fig[3, 2, TopLeft()] =
        Label(fig, "I", fontsize = 12, halign = :right, padding = (0, 15, 5, 0))
label_j =
    fig[3, 3, TopLeft()] =
        Label(fig, "J", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))

fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_FigS5_gsa_and_sims_for_cofactor_pool_sizes_and_AK_CK_NDPK.png", fig, px_per_unit = 4)
