using Glycolysis
using DifferentialEquations
using DataFrames, DataFramesMeta, CSV, Dates
using Distributions, StatsBase, ProgressMeter


#This code takes about 30 min to run for 10 ATPase values with 10,000 repeats each on a 4 core computer

tspan = (0.0, 1e8)
# ATPase_Vmax_frac_list = collect(0.01:0.005:0.15)
# ATPase_Vmax_frac_list = 10 .^ (range(log10(0.015), log10(0.15), 10))
ATPase_Vmax_frac_list = 10 .^ (range(log10(0.02), log10(0.2), 10))
n_repeats = 10_000
# ATPase_Vmax_frac_list = [0.1, 0.2]
glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
glycolysis_init_conc_copy.Citrate = 1.0
glycolysis_init_conc_copy.F26BP = 1.0
glycolysis_init_conc_copy.Phenylalanine = 1.0
# glycolysis_init_conc_copy.Lactate_media = 1.0
# glycolysis_init_conc_copy.Glucose_media = 5e-3

# glycolysis_params_copy.PKM2_K_a_F16BP *= 30
# glycolysis_params_copy.PKM2_L = 0.0
glycolysis_params_copy = deepcopy(glycolysis_params)
glycolysis_params_copy.CK_Vmax = 1.0
glycolysis_params_copy.NDPK_Vmax = 1.0

Data_Total = DataFrame()
Data_Free = DataFrame()
@showprogress for ATPase_Vmax_frac in ATPase_Vmax_frac_list
    local prob = ODEProblem(
        glycolysis_ODEs, glycolysis_init_conc_copy, tspan, glycolysis_params_copy)
    function prob_func(prob, i, repeat)
        prob_copy = remake(
            prob,
            p = rand.(truncated.(
                Normal.(glycolysis_params_copy, glycolysis_params_uncertainty); lower = 0.0)),
            u0 = rand.(
                truncated.(
                Normal.(
                    glycolysis_init_conc_copy, Glycolysis.glycolysis_init_conc_uncertainty);
                lower = 0.0
            )
            )
        )
        Glycolysis_Vmax = min(
            2 * prob_copy.p.HK1_Vmax * prob_copy.p.HK1_Conc,
            2 * prob_copy.p.PFKP_Vmax * prob_copy.p.PFKP_Conc
        )
        prob_copy.p.ATPase_Vmax = ATPase_Vmax_frac * Glycolysis_Vmax

        # below is due to bug in rand where it's hard to sample 0.0 distribution
        prob_copy.u0.Citrate = 0.0
        prob_copy.u0.F26BP = 0.0
        prob_copy.u0.Phenylalanine = 0.0
        prob_copy.u0.Lactate_media = 0.0
        prob_copy.p.CK_Vmax = 0.0
        prob_copy.p.NDPK_Vmax = 0.0
        #remove PKM2 regulation
        # prob_copy.p.PKM2_L = 0.0

        return prob_copy
    end

    local ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    local sim = solve(
        ensemble_prob,
        Rodas5P(),
        EnsembleThreads(),
        trajectories = n_repeats,
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false
    )
    append!(
        Data_Total,
        DataFrame([merge(
                       (ATPase_Vmax_frac = ATPase_Vmax_frac,),
                       convert(
                           NamedTuple, Glycolysis.free_to_total_conc(sol[end], sol.prob.p)),
                       (
                           ATPprod_ATPase_ratios = Glycolysis.conc_to_rates(
                           sol.u[end], sol.prob.p).ATPprod /
                                                   sol.prob.p.ATPase_Vmax,
                       )
                   ) for sol in sim if sol.retcode == ReturnCode.Success])
    )
    append!(
        Data_Free,
        DataFrame([merge(
                       (ATPase_Vmax_frac = ATPase_Vmax_frac,),
                       convert(NamedTuple, sol[end]),
                       (
                           ATPprod_ATPase_ratios = Glycolysis.conc_to_rates(
                           sol.u[end], sol.prob.p).ATPprod /
                                                   sol.prob.p.ATPase_Vmax,
                       )
                   ) for sol in sim if sol.retcode == ReturnCode.Success])
    )
end

CSV.write(
    "Results/$(Dates.format(now(),"mmddyy"))_Glycolysis_Total_Metabolite_Results_$(n_repeats)_reps_w_ATPase_range_$(Int(round(ATPase_Vmax_frac_list[1]*100, sigdigits =1)))_$(Int(round(ATPase_Vmax_frac_list[end]*100, sigdigits=1)))_percent_Lact_media_0_Glucose_media_25.csv",
    Data_Total
);
CSV.write(
    "Results/$(Dates.format(now(),"mmddyy"))_Glycolysis_Free_Metabolite_Results_$(n_repeats)_reps_w_ATPase_range_$(Int(round(ATPase_Vmax_frac_list[1]*100, sigdigits =1)))_$(Int(round(ATPase_Vmax_frac_list[end]*100, sigdigits=1)))_percent_Lact_media_0_Glucose_media_25.csv",
    Data_Total
);

qlow(x) = percentile(x, 2.5)
qhigh(x) = percentile(x, 97.5)

# Filter NaN from DataFrame
filter!(row -> !any(isnan(x) for x in row), Data_Total)

Processed_Data_Total = @chain Data_Total begin

    # Filter rows based on ATP prod ≈ ATPase Vmax
    @rsubset!(1.01>=:ATPprod_ATPase_ratios>=0.99)

    # Calculate mean, qlow and qhigh
    @by(:ATPase_Vmax_frac,
        :count=length(:ATP),
        $(propertynames(Data_Total) .=> median),
        $(propertynames(Data_Total) .=> qlow),
        $(propertynames(Data_Total) .=> qhigh))
end

CSV.write(
    "Results/$(Dates.format(now(),"mmddyy"))_Glycolysis_Processed_Total_Metabolite_Results_$(n_repeats)_reps_w_ATPase_range_$(Int(round(ATPase_Vmax_frac_list[1]*100, sigdigits =1)))_$(Int(round(ATPase_Vmax_frac_list[end]*100, sigdigits=1)))_percent_Lact_media_0_Glucose_media_25.csv",
    Processed_Data_Total
);


Processed_Data_Free = @chain Data_Free begin

    # Filter rows based on ATP prod ≈ ATPase Vmax
    @rsubset!(1.01>=:ATPprod_ATPase_ratios>=0.99)

    # Calculate mean, qlow and qhigh
    @by(:ATPase_Vmax_frac,
        :count=length(:ATP),
        $(propertynames(Data_Total) .=> median),
        $(propertynames(Data_Total) .=> qlow),
        $(propertynames(Data_Total) .=> qhigh))
end

CSV.write(
    "Results/$(Dates.format(now(),"mmddyy"))_Glycolysis_Processed_Free_Metabolite_Results_$(n_repeats)_reps_w_ATPase_range_$(Int(round(ATPase_Vmax_frac_list[1]*100, sigdigits =1)))_$(Int(round(ATPase_Vmax_frac_list[end]*100, sigdigits=1)))_percent_Lact_media_0_Glucose_media_25.csv",
    Processed_Data_Free
);

##
Data_Total = CSV.read(
    "Results/071224_Glycolysis_Total_Metabolite_Results_10000_reps_w_ATPase_range_2_20_percent_Lact_media_0_Glucose_media_25.csv",
    DataFrame
)

sum(Data_Total.ATP .> Data_Total.AMP) / nrow(Data_Total)
sum(Data_Total.ATPprod_ATPase_ratios .> 0.99) / nrow(Data_Total)

data_for_plot = combine(groupby(Data_Total, :ATPase_Vmax_frac),
    [:ATP, :AMP] => (atp, amp) -> sum(atp .> amp) / length(atp))
lines(data_for_plot.ATPase_Vmax_frac, data_for_plot.ATP_AMP_function)
current_figure()
