using Glycolysis, Revise
using DifferentialEquations
using DataFrames, DataFramesMeta, CSV, Distributions, StatsBase, Dates


#
#Calculate 13C tracing ensemble data at one ATPase value with bootstrap 95%CI

tspan = (0.0, 12_000.0)
ATPase_Vmax_frac = 0.15
C13_addition_time = [6_000.0]
n_repeats = 10_000

function affect!(integrator)
    #13C Glucose
    # #Try setting it to 1e-9 instead of zero to make solution more stable
    integrator.u.Glucose_media_12C = 0.0
    integrator.u.Glucose_media_13C = 25e-3

    # #13C Lactate 
    # integrator.u.Lactate_media_12C = 0.0
    # integrator.u.Lactate_media_13C = 5e-3
end

PresetTime_cb = PresetTimeCallback(C13_addition_time, affect!)
cb = PresetTimeCallback(C13_addition_time, affect!)

prob = ODEProblem(
    Glycolysis.glycolysis_13C_tracing_ODEs,
    Glycolysis.glycolysis_13C_tracing_init_conc,
    tspan,
    glycolysis_params,
)

# function prob_func(prob, i, repeat)
#     prob_copy = remake(
#         prob,
#         p = rand.(truncated.(Normal.(glycolysis_params, glycolysis_params_uncertainty); lower = 0.0)),
#     )
#     # Glycolysis_Vmax = min(
#     #     2 * prob_copy.p.HK1_Vmax * prob_copy.p.HK1_Conc,
#     #     2 * prob_copy.p.PFKP_Vmax * prob_copy.p.PFKP_Conc,
#     # )
#     # prob_copy.p.ATPase_Vmax = ATPase_Vmax_frac * Glycolysis_Vmax
#     prob_copy.p.ATPase_Vmax = ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax
#     # prob_copy.p.ATPase_Keq = glycolysis_params.ATPase_Keq
#     # prob_copy.p.AK_Km_ADP = glycolysis_params.AK_Km_ADP
#     # prob_copy.p.AK_Km_ATP = glycolysis_params.AK_Km_ATP
#     # prob_copy.p.AK_Km_AMP = glycolysis_params.AK_Km_AMP
#     # prob_copy.p.AK_Vmax = glycolysis_params.AK_Vmax
#     # prob_copy.p.AK_Keq = glycolysis_params.AK_Keq
#     # prob_copy.p.AK_MW = glycolysis_params.AK_MW
#     # prob_copy.p.ATPase_Km_ATP = glycolysis_params.ATPase_Km_ATP
#     # prob_copy.p.ATPase_Km_ADP = glycolysis_params.ATPase_Km_ADP
#     # prob_copy.p.ATPase_Km_Phosphate = glycolysis_params.ATPase_Km_Phosphate
#     return prob_copy
# end

Glycolysis.glycolysis_13C_tracing_init_conc.Citrate_12C = 1.0
Glycolysis.glycolysis_13C_tracing_init_conc.F26BP_12C = 1.0
Glycolysis.glycolysis_13C_tracing_init_conc.Phenylalanine_12C = 1.0
Glycolysis.glycolysis_13C_tracing_init_conc.Citrate_13C = 1.0
Glycolysis.glycolysis_13C_tracing_init_conc.F26BP_13C = 1.0
Glycolysis.glycolysis_13C_tracing_init_conc.Phenylalanine_13C = 1.0

function prob_func(prob, i, repeat)
    prob.p .= rand.(truncated.(Normal.(glycolysis_params, glycolysis_params_uncertainty); lower = 0.0))
    prob.u0 .=
        rand.(
            truncated.(
                Normal.(
                    Glycolysis.glycolysis_13C_tracing_init_conc,
                    Glycolysis.glycolysis_13C_tracing_init_conc_uncertainty,
                );
                lower = 0.0,
            )
        )
    # Glycolysis_Vmax = min(
    #     2 * prob.p.HK1_Vmax * prob.p.HK1_Conc,
    #     2 * prob.p.PFKP_Vmax * prob.p.PFKP_Conc,
    # )
    # prob.p.ATPase_Vmax = ATPase_Vmax_frac * Glycolysis_Vmax
    prob.p.ATPase_Vmax = ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax

    # below is due to bug in rand where it's hard to sample 0.0 distribution
    prob.u0.Lactate_media_12C = 0.0
    prob.u0.Lactate_media_13C = 0.0
    prob.u0.Citrate_12C = 0.0
    prob.u0.F26BP_12C = 0.0
    prob.u0.Phenylalanine_12C = 0.0
    prob.u0.Citrate_13C = 0.0
    prob.u0.F26BP_13C = 0.0
    prob.u0.Phenylalanine_13C = 0.0
    return prob
end

ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)

sim = solve(
    ensemble_prob,
    Rodas5P(),
    EnsembleThreads(),
    trajectories = n_repeats,
    abstol = 1e-12,
    reltol = 1e-5,
    callback = cb,
    saveat = [k for k = 5999:0.1:6030],
)

#
# Process 13C tracing ensemble data to calculate labeling fraction mean +/- 95%CI
Data_raw = DataFrame()
for (k, sol) in enumerate(sim)
    if sol.retcode == ReturnCode.Success# && all(sol.u[end] .>= -1e-16)
        # if sol.u[end].ATP > sol.u[end].AMP
        Glucose = sol.u[end].Glucose_12C + sol.u[end].Glucose_13C
        G6P = sol.u[end].G6P_12C + sol.u[end].G6P_13C
        F6P = sol.u[end].F6P_12C + sol.u[end].F6P_13C
        F16BP = sol.u[end].F16BP_12C + sol.u[end].F16BP_13C
        GAP = sol.u[end].GAP_12C + sol.u[end].GAP_13C
        BPG = sol.u[end].BPG_12C + sol.u[end].BPG_13C
        ThreePG = sol.u[end].ThreePG_12C + sol.u[end].ThreePG_13C
        PEP = sol.u[end].PEP_12C + sol.u[end].PEP_13C
        Pyruvate = sol.u[end].Pyruvate_12C + sol.u[end].Pyruvate_13C
        ATP = sol.u[end].ATP
        ADP = sol.u[end].ADP
        AMP = sol.u[end].AMP
        Phosphate = sol.u[end].Phosphate
        Citrate = sol.u[end].Citrate_12C + sol.u[end].Citrate_13C
        F26BP = sol.u[end].F26BP_12C + sol.u[end].F26BP_13C
        Phenylalanine = sol.u[end].Phenylalanine_12C + sol.u[end].Phenylalanine_13C
        ATPprod = (
            -Glycolysis.rate_HK1(Glucose, G6P, ATP, Phosphate, ADP, sol.prob.p) -
            Glycolysis.rate_PFKP(F6P, ATP, F16BP, ADP, Phosphate, Citrate, F26BP, sol.prob.p) +
            Glycolysis.rate_PGK(BPG, ADP, ATP, ThreePG, sol.prob.p) +
            Glycolysis.rate_PKM2(PEP, ADP, Pyruvate, ATP, F16BP, Phenylalanine, sol.prob.p) +
            Glycolysis.rate_AK(ATP, ADP, AMP, sol.prob.p)
        )
        if ATPprod > 0.99 * sol.prob.p.ATPase_Vmax
            append!(
                Data_raw,
                DataFrame([
                    merge((trajectory = k, time = (timepoint - 6000)), convert(NamedTuple, sol[i])) for
                    (i, timepoint) in enumerate(sol.t)
                ]),
            )
        end
    end
end
print("$(100 * length(unique(Data_raw.trajectory)) / n_repeats)% trajectories supported ATPase rate")

metabolite_names = replace.(names(Data_raw[:, Between(:Glucose_12C, :Lactate_12C)]), "_12C" => "")
Data_processed = DataFrame(
    [
        Data_raw[!, name*"_13C"] ./ (Data_raw[!, name*"_13C"] .+ Data_raw[!, name*"_12C"]) for
        name in metabolite_names
    ],
    metabolite_names,
)

Data_processed = hcat(Data_raw[!, Cols(:time, :trajectory)], Data_processed)
qlow(x) = percentile(x, 2.5)
qhigh(x) = percentile(x, 97.5)
Data = combine(
    groupby(Data_processed, :time),
    names(Data_processed) .=> mean,
    names(Data_processed) .=> qlow,
    names(Data_processed) .=> qhigh,
)
rename!(Data, replace.(names(Data), "_12C" => ""))

CSV.write(
    "/Users/Denis/Library/Mobile Documents/com~apple~CloudDocs/Research Projects/Glycolysis Model/JuliaGlycolysisModel/Results data and figures/$(Dates.format(now(),"mmddyy"))_13C_glucose_labeling_w_CI_$(ATPase_Vmax_frac).csv",
    Data,
);
# CSV.write(
#     "/Users/Denis/Library/Mobile Documents/com~apple~CloudDocs/Research Projects/Glycolysis Model/JuliaGlycolysisModel/Results data and figures/$(Dates.format(now(),"mmddyy"))_13C_lactate_labeling_w_CI_$(ATPase_Vmax_frac).csv",
#     Data,
# );

##
metabolite_names = propertynames(Data_processed[:, Between(:Glucose, :Lactate)])
fig = Figure()
ax = Axis(fig[1, 1])
for name in metabolite_names[1:3]
    lines!(ax, Data_processed.time, Data_processed[:, name])
end
fig

##
Data_raw[:, "Glucose_12C"]
name = "Glucose"
metabolite_names = replace.(names(Data_raw[:, Between(:Glucose_12C, :Lactate_12C)]), "_12C" => "")
DataFrame(
    [
        Data_raw[!, name*"_13C"] ./ (Data_raw[!, name*"_13C"] .+ Data_raw[!, name*"_12C"]) for
        name in metabolite_names
    ],
    metabolite_names,
)