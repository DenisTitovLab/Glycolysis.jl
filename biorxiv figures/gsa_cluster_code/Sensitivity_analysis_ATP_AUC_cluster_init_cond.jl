cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Distributed, ClusterManagers

# ENV["SLURM_JOB_CPUS_PER_NODE"] give number of cores in "40(x4),32,20" format
# "40(x2),32,20" can be parsed to 132 using code below to get total number of allocated cores

subs = Dict("x" => "*", "(" => "", ")" => "");
np = sum(eval(Meta.parse(replace(ENV["SLURM_JOB_CPUS_PER_NODE"], r"x|\(|\)" => s -> subs[s]))))
addprocs(SlurmManager(np); exeflags = "--project")

# using Distributed
# addprocs(8; exeflags = "--project")

@everywhere using Glycolysis
@everywhere using DifferentialEquations
using GlobalSensitivity, QuasiMonteCarlo, LabelledArrays, DataFrames, CSV, Dates

# write a function outputting AUC of [ATP] estimated as the sum of [ATP] over log scale
@everywhere function glycolysis_output(init_cond_scaler, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values)
    scaled_init_cond = deepcopy(glycolysis_init_conc)
    scaled_init_cond = init_cond_scaler .* scaled_init_cond
    tspan = (0.0, 1e8)
    Glycolysis_Vmax = 2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc
    ATPases = 10 .^ range(log10(0.001), log10(1.0), n_Vmax_ATPase_values) .* Glycolysis_Vmax

    Adenine_pool = scaled_init_cond.ATP + scaled_init_cond.ADP + scaled_init_cond.AMP
    prob = ODEProblem(glycolysis_ODEs, scaled_init_cond, tspan, glycolysis_params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(
        ensemble_prob,
        Rodas4(),
        EnsembleSerial(),
        trajectories = n_Vmax_ATPase_values,
        abstol = 1e-12,
        reltol = 1e-5,
        save_everystep = false,
        save_start = false,
    )

    sum_ATP = 0.0
    for sol in sim
        if sol.retcode == :Success
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

# ##
# using CairoMakie
# fold_range = 3
# n_Vmax_ATPase_values = 100
# res = @time pmap(
#     x -> glycolysis_output(x, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values),
#     [sqrt(fold_range) .^ (-2 .+ 4 * rand(length(glycolysis_init_conc))) for i = 1:10_000],
# ) 
# using CSV
# CSV.write("$(Dates.format(now(),"mmddyy"))_hist_ATP_AUC_10K_3x_init_cond.csv", DataFrame(all_params = res))
# hist(res, bins = Base.range(0.0, 1.0, 50))
# ##

# generate design matrices
n_bootstrap = 20_000 #number of bootstrapped datasets to use
lb = ones(length(glycolysis_init_conc)) / 3 #lower bound for parameter
ub = ones(length(glycolysis_init_conc)) * 3 #upper bound for parameter
sampler = SobolSample()
A, B = QuasiMonteCarlo.generate_design_matrices(n_bootstrap, lb, ub, sampler)

# do global sensitivity analysis
n_Vmax_ATPase_values = 100
sobol_sens = gsa(
    x -> parallel_output(x, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values),
    Sobol(),
    A,
    B;
    batch = true,
)

S1_labelled = @LArray sobol_sens.S1 propertynames(glycolysis_init_conc)
ST_labelled = @LArray sobol_sens.ST propertynames(glycolysis_init_conc)

df = DataFrame()
push!(df, merge((Row = "S1",), convert(NamedTuple, S1_labelled)))
push!(df, merge((Row = "ST",), convert(NamedTuple, ST_labelled)))

CSV.write(
    "/global/home/users/titov/gsa_cluster_code/$(Dates.format(now(),"mmddyy"))_ATP_AUC_gsa_sobol_init_cond_$(n_bootstrap)_3x_range.csv",
    df,
)

# Remove all the workers
for n in workers()
    rmprocs(n)
end

