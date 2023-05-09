using Glycolysis
using OrdinaryDiffEq
using DataFrames, DataFramesMeta, CSV, Dates
using Distributions, StatsBase, ProgressMeter

tspan = (0.0, 1e8)
# ATPase_Vmax_frac_list = collect(0.01:0.005:0.15)
ATPase_Vmax_frac_list = 10 .^ (range(log10(0.015), log10(0.15), 10))
n_repeats = 10_000
# ATPase_Vmax_frac_list = [0.1, 0.2]

glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
glycolysis_init_conc_copy.Citrate = 1.0
glycolysis_init_conc_copy.F26BP = 1.0
glycolysis_init_conc_copy.Phenylalanine = 1.0

Data = DataFrame()
@showprogress for ATPase_Vmax_frac in ATPase_Vmax_frac_list
    local prob = ODEProblem(glycolysis_ODEs, glycolysis_init_conc_copy, tspan, glycolysis_params)
    function prob_func(prob, i, repeat)
        # prob_copy = remake(
        #     prob,
        #     p = rand.(truncated.(Normal.(glycolysis_params, glycolysis_params_uncertainty); lower = 0.0)),
        # )
        # prob_copy = remake(prob)
        prob.p .=
            rand.(truncated.(Normal.(glycolysis_params, glycolysis_params_uncertainty); lower = 0.0))
        prob.u0 .=
            rand.(
                truncated.(
                    Normal.(glycolysis_init_conc_copy, Glycolysis.glycolysis_init_conc_uncertainty);
                    lower = 0.0,
                )
            )
        Glycolysis_Vmax = min(
            2 * prob.p.HK1_Vmax * prob.p.HK1_Conc,
            2 * prob.p.PFKP_Vmax * prob.p.PFKP_Conc,
        )
        prob.p.ATPase_Vmax = ATPase_Vmax_frac * Glycolysis_Vmax
        # prob.p.ATPase_Vmax =
        #     ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax
        # prob.p.ATPase_Keq = glycolysis_params.ATPase_Keq
        # prob.p.AK_Km_ADP = glycolysis_params.AK_Km_ADP
        # prob.p.AK_Km_ATP = glycolysis_params.AK_Km_ATP
        # prob.p.AK_Km_AMP = glycolysis_params.AK_Km_AMP
        # prob.p.AK_Vmax = glycolysis_params.AK_Vmax
        # prob.p.AK_Keq = glycolysis_params.AK_Keq
        # prob.p.AK_MW = glycolysis_params.AK_MW
        # prob.p.ATPase_Km_ATP = glycolysis_params.ATPase_Km_ATP
        # prob.p.ATPase_Km_ADP = glycolysis_params.ATPase_Km_ADP
        # prob.p.ATPase_Km_Phosphate = glycolysis_params.ATPase_Km_Phosphate

        # below is due to bug in rand where it's hard to sample 0.0 distribution
        prob.u0.Citrate = 0.0
        prob.u0.F26BP = 0.0
        prob.u0.Phenylalanine = 0.0
        # prob.u0.Lactate_media = 0.0
        return prob
    end

    local ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    local sim = solve(
        ensemble_prob,
        Rodas5P(),
        EnsembleThreads(),
        trajectories = n_repeats,
        abstol = 1e-15,
        reltol = 1e-5,
        save_everystep = false,
    )
    append!(
        Data,
        DataFrame([
            merge(
                (ATPase_Vmax_frac = ATPase_Vmax_frac,),
                convert(NamedTuple, Glycolysis.free_to_total_conc(sol[end], sol.prob.p)),
                (
                    ATPprod_ATPase_ratios = Glycolysis.conc_to_rates(sol.u[end], sol.prob.p).ATPprod /
                                            sol.prob.p.ATPase_Vmax,
                ),
            ) for sol in sim if sol.retcode == ReturnCode.Success
        ]),
    )
    # append!(
    #     Data,
    #     DataFrame([
    #         merge(
    #             (ATPase_Vmax_frac = ATPase_Vmax_frac,),
    #             convert(NamedTuple, sol[end]),
    #             (
    #                 ATPprod_ATPase_ratios = Glycolysis.conc_to_rates(sol.u[end], sol.prob.p).ATPprod /
    #                                         sol.prob.p.ATPase_Vmax,
    #             ),
    #         ) for sol in sim if sol.retcode == ReturnCode.Success
    #     ]),
    # )
end

CSV.write(
    "/Users/Denis/Library/Mobile Documents/com~apple~CloudDocs/Research Projects/Glycolysis Model/JuliaGlycolysisModel/Results data and figures/$(Dates.format(now(),"mmddyy"))_Glycolysis_Results_$(n_repeats)_reps_w_CI_min_HK_PFK.csv",
    Data,
);

qlow(x) = percentile(x, 2.5)
qhigh(x) = percentile(x, 97.5)

Processed_Data = @chain Data begin

    # Filter rows based on ATP prod ≈ ATPase Vmax
    @rsubset!(:ATPprod_ATPase_ratios > 0.9)

    # Calculate ATP/ADP Pi and Energy charge
    # @rtransform!(
    #     :ATP_ADP_Pi = :ATP / (:ADP * :Phosphate),
    #     :EnergyCharge = (:ATP + 0.5 * :ADP) / (:ATP + :ADP + :AMP)
    # )

    # Calculate mean, qlow and qhigh
    @by(
        :ATPase_Vmax_frac,
        :count = length(:ATP),
        $(propertynames(Data) .=> median),
        $(propertynames(Data) .=> qlow),
        $(propertynames(Data) .=> qhigh)
    )
end


CSV.write(
    "/Users/Denis/Library/Mobile Documents/com~apple~CloudDocs/Research Projects/Glycolysis Model/JuliaGlycolysisModel/Results data and figures/$(Dates.format(now(),"mmddyy"))_Glycolysis_Processed_Results_$(n_repeats)_reps_w_CI_min_HK_PFK.csv",
    Processed_Data,
);

##
Data = CSV.read(
    "/Users/Denis/Library/Mobile Documents/com~apple~CloudDocs/Research Projects/Glycolysis Model/JuliaGlycolysisModel/Results data and figures/042023_Glycolysis_Results_10000_reps_w_CI_min_HK_PFK.csv",
    DataFrame,
)

sum(Data.ATP .> Data.AMP) / nrow(Data)
sum(Data.ATPprod_ATPase_ratios .> 0.99) / nrow(Data)

data_for_plot =
    combine(groupby(Data, :ATPase_Vmax_frac), [:ATP, :AMP] => (atp, amp) -> sum(atp .> amp) / length(atp))
lines(data_for_plot.ATPase_Vmax_frac, data_for_plot.ATP_AMP_function)
current_figure()