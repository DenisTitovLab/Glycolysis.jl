cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.resolve()
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
@everywhere function glycolysis_output(cofactor_pool_scaler, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values)
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


# ##
# using CairoMakie
# fold_range = 3
# n_Vmax_ATPase_values = 100
# n_bootstrap = 10_000
# res = @time pmap(
#     x -> glycolysis_output(x, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values),
#     [sqrt(fold_range) .^ (-2 .+ 4 * rand(3)) for i = 1:n_bootstrap],
# ) 
# using CSV
# CSV.write("$(Dates.format(now(),"mmddyy"))_hist_ATP_AUC_$(n_bootstrap)_runs_3x_pool_sizes.csv", DataFrame(all_params = res))
# hist(res, bins = Base.range(0.0, 1.0, 50))
# ##

# generate design matrices
n_bootstrap = 20_000 #number of bootstrapped datasets to use
lb = ones(3) / 3 #lower bound for parameter
ub = ones(3) * 3 #upper bound for parameter
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

S1_labelled = @LArray sobol_sens.S1 (:Adenine_pool_scaler, :NAD_pool_scaler, :Pi_scaler)
ST_labelled = @LArray sobol_sens.ST (:Adenine_pool_scaler, :NAD_pool_scaler, :Pi_scaler)

df = DataFrame()
push!(df, merge((Row = "S1",), convert(NamedTuple, S1_labelled)))
push!(df, merge((Row = "ST",), convert(NamedTuple, ST_labelled)))

CSV.write(
    "/global/home/users/titov/gsa_cluster_code/$(Dates.format(now(),"mmddyy"))_ATP_AUC_gsa_sobol_cofactor_pool_$(n_bootstrap)_3x_range.csv",
    df,
)

# Remove all the workers
for n in workers()
    rmprocs(n)
end

