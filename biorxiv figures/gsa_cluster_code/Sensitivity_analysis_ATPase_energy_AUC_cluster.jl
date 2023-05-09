# cd(@__DIR__)
# using Pkg
# Pkg.activate(".")
# Pkg.resolve()
# Pkg.instantiate()

# using Distributed, ClusterManagers

# # ENV["SLURM_JOB_CPUS_PER_NODE"] give number of cores in "40(x4),32,20" format
# # "40(x2),32,20" can be parsed to 132 using code below to get total number of allocated cores

# subs = Dict("x" => "*", "(" => "", ")" => "");
# np = sum(eval(Meta.parse(replace(ENV["SLURM_JOB_CPUS_PER_NODE"], r"x|\(|\)" => s -> subs[s]))))
# addprocs(SlurmManager(np); exeflags = "--project")

using Distributed
addprocs(8; exeflags = "--project")

@everywhere using Glycolysis
@everywhere using DifferentialEquations
using GlobalSensitivity, QuasiMonteCarlo, LabelledArrays, DataFrames, CSV, Dates

# write a function outputting AUC of [ATP] estimated as the sum of [ATP] over log scale
@everywhere function glycolysis_output(params_scaler, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values)
    scaled_params = params_scaler .* glycolysis_params
    tspan = (0.0, 1e8)
    
    #Value below has to be calculated from 2 * glycolysis_params.HK1_Vmax * glycolysis_params.HK1_Conc to prevent the latter from influencing analysis through affect of ATPases rate
    # ATPases = 10 .^ range(log10(0.005), log10(0.5), n_Vmax_ATPase_values) .* 0.02920110192837465
    
    Glycolysis_Vmax = min(
        # 2 * scaled_params.GLUT_Vmax * scaled_params.GLUT_Conc,
        2 * scaled_params.HK1_Vmax * scaled_params.HK1_Conc,
        # 2 * scaled_params.GPI_Vmax * scaled_params.GPI_Conc,
        2 * scaled_params.PFKP_Vmax * scaled_params.PFKP_Conc,
        # 2 * scaled_params.ALDO_Vmax * scaled_params.ALDO_Conc,
        # 2 * scaled_params.TPI_Vmax * scaled_params.TPI_Conc,
        # scaled_params.GAPDH_Vmax * scaled_params.GAPDH_Conc,
        # scaled_params.PGK_Vmax * scaled_params.PGK_Conc,
        # scaled_params.PGM_Vmax * scaled_params.PGM_Conc,
        # scaled_params.ENO_Vmax * scaled_params.ENO_Conc,
        # scaled_params.PKM2_Vmax * scaled_params.PKM2_Conc,
        # scaled_params.LDH_Vmax * scaled_params.LDH_Conc,
        # scaled_params.MCT_Vmax * scaled_params.MCT_Conc,
    )
    ATPases = 10 .^ range(log10(0.001), log10(1.0), n_Vmax_ATPase_values) .* Glycolysis_Vmax

    #Keep glycolysis_params for ATPase and AK the same
    scaled_params.ATPase_Keq = glycolysis_params.ATPase_Keq
    scaled_params.AK_Km_ADP = glycolysis_params.AK_Km_ADP
    scaled_params.AK_Km_ATP = glycolysis_params.AK_Km_ATP
    scaled_params.AK_Km_AMP = glycolysis_params.AK_Km_AMP
    scaled_params.AK_Vmax = glycolysis_params.AK_Vmax
    scaled_params.AK_Keq = glycolysis_params.AK_Keq
    scaled_params.AK_MW = glycolysis_params.AK_MW
    scaled_params.ATPase_Km_ATP = glycolysis_params.ATPase_Km_ATP
    scaled_params.ATPase_Km_ADP = glycolysis_params.ATPase_Km_ADP
    scaled_params.ATPase_Km_Phosphate = glycolysis_params.ATPase_Km_Phosphate

    prob = ODEProblem(glycolysis_ODEs, glycolysis_init_conc, tspan, scaled_params)
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
    
    sum_ATP_Energy = 0.0
    for sol in sim
        if sol.retcode == ReturnCode.Success
            sum_ATP_Energy += -log(Glycolysis.conc_to_disequilibrium_ratios(sol.u[end], scaled_params).Q_Keq_ATPase)
        else
            sum_ATP_Energy += 0.0
        end
    end
    return sum_ATP_Energy / n_Vmax_ATPase_values
end

function parallel_output(p, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values)
    res = pmap(
        x -> glycolysis_output(x, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values),
        [p[:, i] for i in axes(p, 2)],
    )
    return permutedims(res)
end

##
using CairoMakie, Statistics
fold_range = 3
n_Vmax_ATPase_values = 100
n_bootstrap = 10_000
res = pmap(
    x -> glycolysis_output(x, glycolysis_params, glycolysis_init_conc, n_Vmax_ATPase_values),
    [sqrt(fold_range) .^ (-2 .+ 4 * rand(length(glycolysis_params))) for i = 1:n_bootstrap],
)

using CSV
CSV.write("$(Dates.format(now(),"mmddyy"))_hist_ATP_energy_AUC_$(n_bootstrap)runs_only_HK_PFK_Vmax.csv", DataFrame(all_params = res))
print(std(res)/mean(res))
hist(res, bins=50)
##

# generate design matrices
n_bootstrap = 20_000 #number of bootstrapped datasets to use
lb = ones(length(glycolysis_params)) / 3 #lower bound for parameter
ub = ones(length(glycolysis_params)) * 3 #upper bound for parameter
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

S1_labelled = @LArray sobol_sens.S1 propertynames(glycolysis_params)
ST_labelled = @LArray sobol_sens.ST propertynames(glycolysis_params)

df = DataFrame()
push!(df, merge((Row = "S1",), convert(NamedTuple, S1_labelled)))
push!(df, merge((Row = "ST",), convert(NamedTuple, ST_labelled)))

CSV.write(
    "/global/home/users/titov/gsa_cluster_code/$(Dates.format(now(),"mmddyy"))_ATPase_energy_AUC_gsa_sobol_$(n_bootstrap)_3x_range_HK_PFK_Vmax.csv",
    df,
)

# Remove all the workers
for n in workers()
    rmprocs(n)
end

