cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.instantiate()

#=
Comments about duration and memory use:
atol and btol affect the memory usage and time of the calculation.
atol = 1e-15 and btol = 1e-8 are 10s per run and 200MB per run
atol = 1e-14 and btol = 1e-7 are 5s per run and 100MB per run
atol = 1e-13 and btol = 1e-6 are 2s per run and 60MB per run
It might be beneficial to decrease tolerance from default atol = 1e-15 and btol = 1e-8 to ensure that
the calculation on a cluster is not killed due to memory usage and doesn;t take too many hours.
=#

using Distributed, ClusterManagers

#Increase time before worker terminates without hearing from master process. Useful for large number of cores.
ENV["JULIA_WORKER_TIMEOUT"] = 600.0

# ENV["SLURM_JOB_CPUS_PER_NODE"] give number of cores in "40(x4),32,20" format
# "40(x2),32,20" can be parsed to 132 using code below to get total number of allocated cores

subs = Dict("x" => "*", "(" => "", ")" => "");
np = sum(eval(Meta.parse(replace(ENV["SLURM_JOB_CPUS_PER_NODE"], r"x|\(|\)" => s -> subs[s]))))
addprocs(SlurmManager(np); exeflags = "--project")

# # comment all lines above and uncomment two line below if running on local computer
# using Distributed
# addprocs(8; exeflags = "--project")

@everywhere using Glycolysis
@everywhere using DifferentialEquations
using GlobalSensitivity, QuasiMonteCarlo, LabelledArrays, DataFrames, CSV, Dates

2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax
# write a function outputting AUC of [ATP] estimated as the sum of [ATP] over log scale
@everywhere function glycolysis_output(
    params_scaler,
    glycolysis_params,
    glycolysis_init_conc,
    n_Vmax_ATPase_values,
)
    scaled_params = params_scaler .* glycolysis_params
    tspan = (0.0, 1e8)

    # Glycolysis_Vmax = min(
    #     # 2 * scaled_params.GLUT_Vmax * scaled_params.GLUT_Conc,
    #     2 * scaled_params.HK1_Vmax * scaled_params.HK1_Conc,
    #     # 2 * scaled_params.GPI_Vmax * scaled_params.GPI_Conc,
    #     2 * scaled_params.PFKP_Vmax * scaled_params.PFKP_Conc,
    #     # 2 * scaled_params.ALDO_Vmax * scaled_params.ALDO_Conc,
    #     # 2 * scaled_params.TPI_Vmax * scaled_params.TPI_Conc,
    #     # scaled_params.GAPDH_Vmax * scaled_params.GAPDH_Conc,
    #     # scaled_params.PGK_Vmax * scaled_params.PGK_Conc,
    #     # scaled_params.PGM_Vmax * scaled_params.PGM_Conc,
    #     # scaled_params.ENO_Vmax * scaled_params.ENO_Conc,
    #     # scaled_params.PKM2_Vmax * scaled_params.PKM2_Conc,
    #     # scaled_params.LDH_Vmax * scaled_params.LDH_Conc,
    #     # scaled_params.MCT_Vmax * scaled_params.MCT_Conc,
    # )

    #Value below has to be calculated from 2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc to prevent the latter from influencing analysis through affect of ATPases rate
    Glycolysis_Vmax = 2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc
    ATPases = 10 .^ range(log10(0.001), log10(1.0), n_Vmax_ATPase_values) .* Glycolysis_Vmax

    scaled_params.ATPase_Keq = glycolysis_params.ATPase_Keq
    scaled_params.AK_Km_ADP = glycolysis_params.AK_Km_ADP
    scaled_params.AK_Km_ATP = glycolysis_params.AK_Km_ATP
    scaled_params.AK_Km_AMP = glycolysis_params.AK_Km_AMP
    scaled_params.AK_Vmax = glycolysis_params.AK_Vmax
    scaled_params.AK_Keq = glycolysis_params.AK_Keq
    scaled_params.AK_MW = glycolysis_params.AK_MW
    scaled_params.ATPase_Km_ATP = 1e-9
    scaled_params.ATPase_Km_ADP = glycolysis_params.ATPase_Km_ADP
    scaled_params.ATPase_Km_Phosphate = glycolysis_params.ATPase_Km_Phosphate

    # scaled_params.HK1_Conc = glycolysis_params.HK1_Conc
    # scaled_params.HK1_Vmax = glycolysis_params.HK1_Vmax
    # scaled_params.HK1_K_a_ADP = glycolysis_params.HK1_K_a_ADP
    # scaled_params.HK1_K_a_ATP = glycolysis_params.HK1_K_a_ATP
    # scaled_params.HK1_K_Glucose = glycolysis_params.HK1_K_Glucose
    # scaled_params.HK1_K_G6P = glycolysis_params.HK1_K_G6P
    # scaled_params.HK1_β_Glucose_ATP = glycolysis_params.HK1_β_Glucose_ATP

    prob = ODEProblem(glycolysis_ODEs, glycolysis_init_conc, tspan, scaled_params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(
        ensemble_prob,
        Rodas5P(),
        EnsembleSerial(),
        trajectories = n_Vmax_ATPase_values,
        abstol = 1e-13,
        reltol = 1e-6,
        save_everystep = false,
        save_start = false,
    )

    ratio_ATPprod_ATPase = 0.0
    for (i, sol) in enumerate(sim)
        if sol.retcode == ReturnCode.Success
            ratio_ATPprod_ATPase += Glycolysis.conc_to_rates(sol.u[end], scaled_params).ATPprod / ATPases[i]
        else
            ratio_ATPprod_ATPase += 0.0
        end
    end
    return ratio_ATPprod_ATPase / n_Vmax_ATPase_values
end

function parallel_output(p, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values)
    res = pmap(
        x -> glycolysis_output(x, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values),
        [p[:, i] for i in axes(p, 2)],
    )
    return permutedims(res)
end

##
fold_range = 10
n_Vmax_ATPase_values = 100
n_bootstrap = 10_000
res = pmap(
    x -> glycolysis_output(x, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values),
    [sqrt(fold_range) .^ (-2 .+ 4 * rand(length(glycolysis_params))) for i = 1:n_bootstrap],
)
CSV.write("$(Dates.format(now(),"mmddyy"))_hist_ATP_prod_AUC_$(n_bootstrap)runs_$(fold_range)x_range.csv", DataFrame(all_params = res))
##

# generate design matrices
n_samples = 10_000 #number of bootstrapped datasets to use
fold_range = 10
lb = ones(length(glycolysis_params)) / fold_range #lower bound for parameter
ub = ones(length(glycolysis_params)) * fold_range #upper bound for parameter
sampler = SobolSample()
A, B = QuasiMonteCarlo.generate_design_matrices(n_samples, lb, ub, sampler)

# do global sensitivity analysis
n_Vmax_ATPase_values = 100
sobol_sens = @time gsa(
    x -> parallel_output(x, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values),
    Sobol(),
    A,
    B;
    batch = true,
)

S1 = @LArray sobol_sens.S1 propertynames(glycolysis_params)
ST = @LArray sobol_sens.ST propertynames(glycolysis_params)


df = DataFrame()
push!(df, merge((Row = "S1",), convert(NamedTuple, S1)))
push!(df, merge((Row = "ST",), convert(NamedTuple, ST)))

CSV.write("$(Dates.format(now(),"mmddyy"))_ATPprod_AUC_gsa_sobol_$(n_samples)_$(fold_range)x_range.csv", df)

# Remove all the workers
for n in workers()
    rmprocs(n)
end
