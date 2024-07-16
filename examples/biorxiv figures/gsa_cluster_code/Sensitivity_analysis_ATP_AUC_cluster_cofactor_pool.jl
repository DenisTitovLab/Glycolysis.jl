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

# write a function outputting AUC of [ATP] estimated as the sum of [ATP] over log scale
@everywhere function glycolysis_output(
    cofactor_pool_scaler,
    glycolysis_params,
    glycolysis_init_conc,
    n_Vmax_ATPase_values,
)
    Adenine_pool_scaler, NAD_pool_scaler, Pi_scaler = cofactor_pool_scaler
    tspan = (0.0, 1e8)
    Glycolysis_Vmax = 2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc
    ATPases = 10 .^ range(log10(0.001), log10(1.0), n_Vmax_ATPase_values) .* Glycolysis_Vmax
    scaled_init_conc = deepcopy(glycolysis_init_conc)
    scaled_init_conc.ATP = Adenine_pool_scaler * glycolysis_init_conc.ATP
    scaled_init_conc.ADP = Adenine_pool_scaler * glycolysis_init_conc.ADP
    scaled_init_conc.AMP = Adenine_pool_scaler * glycolysis_init_conc.AMP

    scaled_init_conc.NAD = NAD_pool_scaler * glycolysis_init_conc.NAD
    scaled_init_conc.NADH = NAD_pool_scaler * glycolysis_init_conc.NADH

    scaled_init_conc.G6P = Pi_scaler * glycolysis_init_conc.G6P
    scaled_init_conc.F6P = Pi_scaler * glycolysis_init_conc.F6P
    scaled_init_conc.F16BP = Pi_scaler * glycolysis_init_conc.F16BP
    scaled_init_conc.GAP = Pi_scaler * glycolysis_init_conc.GAP
    scaled_init_conc.DHAP = Pi_scaler * glycolysis_init_conc.DHAP
    scaled_init_conc.BPG = Pi_scaler * glycolysis_init_conc.BPG
    scaled_init_conc.ThreePG = Pi_scaler * glycolysis_init_conc.ThreePG
    scaled_init_conc.TwoPG = Pi_scaler * glycolysis_init_conc.TwoPG
    scaled_init_conc.PEP = Pi_scaler * glycolysis_init_conc.PEP
    scaled_init_conc.Phosphate = Pi_scaler * glycolysis_init_conc.Phosphate

    Adenine_pool = scaled_init_conc.ATP + scaled_init_conc.ADP + scaled_init_conc.AMP
    prob = ODEProblem(glycolysis_ODEs, scaled_init_conc, tspan, glycolysis_params)
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

    sum_ATP = 0.0
    for sol in sim
        if sol.retcode == ReturnCode.Success
            sum_ATP += sol.u[end].ATP / Adenine_pool
        else
            sum_ATP += 0.0
        end
    end
    return sum_ATP / n_Vmax_ATPase_values
end

function parallel_output(p, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values)
    res = pmap(
        x -> glycolysis_output(x, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values),
        [p[:, i] for i in axes(p, 2)],
    )
    return permutedims(res)
end

##
# Calculate data for histogram of variance
fold_range = 3
n_Vmax_ATPase_values = 100
n_bootstrap = 10_000
res = @time pmap(
    x -> glycolysis_output(x, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values),
    [sqrt(fold_range) .^ (-2 .+ 4 * rand(3)) for i = 1:n_bootstrap],
)
CSV.write("$(Dates.format(now(),"mmddyy"))_hist_ATP_AUC_$(n_bootstrap)_runs_$(fold_range)x_range_pool_sizes.csv", DataFrame(all_params = res))
##

# generate design matrices
n_samples = 10_000 #number of bootstrapped datasets to use
n_params = 3
fold_range = 3
lb = ones(n_params) / fold_range #lower bound for parameter
ub = ones(n_params) * fold_range #upper bound for parameter
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

# save the results
S1 = @LArray sobol_sens.S1 (:Adenine_pool_scaler, :NAD_pool_scaler, :Pi_scaler)
ST = @LArray sobol_sens.ST (:Adenine_pool_scaler, :NAD_pool_scaler, :Pi_scaler)

df = DataFrame()
push!(df, merge((Row = "S1",), convert(NamedTuple, S1)))
push!(df, merge((Row = "ST",), convert(NamedTuple, ST)))

CSV.write(
    "$(Dates.format(now(),"mmddyy"))_ATP_AUC_gsa_sobol_cofactor_pool_$(n_samples)_$(fold_range)x_range.csv",
    df,
)

# Remove all the workers
for n in workers()
    rmprocs(n)
end
