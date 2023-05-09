include("MinimalGlycolysisModel_w_feedback_w_constant_pi.jl")
using OrdinaryDiffEq, CairoMakie, LabelledArrays
using DataFrames, CSV, Dates
using FileIO

##
#Precalculate MC simulation of two enzyme model
using Distributed
addprocs(8; exeflags = "--project")
@everywhere include("MinimalGlycolysisModel_w_feedback_w_constant_pi.jl")
@everywhere using OrdinaryDiffEq

using DataFrames, CSV, Dates

@everywhere function minimal_model_output(params, initial_concentrations, n_Vmax_ATPase_values)
    tspan = (0.0, 1e8)
    adenine_pool_size = initial_concentrations.ATP + initial_concentrations.ADP
    pathway_Vmax = min(params.Enz1_Vmax, params.Enz2_Vmax)

    ATPases = 10 .^ range(log10(0.01), log10(1.0), n_Vmax_ATPase_values) .* pathway_Vmax
    prob = ODEProblem(minimal_glycolysis_ODEs, initial_concentrations, tspan, params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(
        ensemble_prob,
        Rodas4P2(),
        EnsembleThreads(),
        trajectories = n_Vmax_ATPase_values,
        abstol = 1e-15,
        reltol = 1e-8,
        save_everystep = false,
        save_start = false,
    )
    ATP_prod_eq_cons_frac = 0.0
    ATP_conc = 0.0
    ATP_energy = 0.0
    counter = 0
    for sol in sim
        if sol.retcode == ReturnCode.Success && all(sol.u[end] .> 0)
            if isapprox(
                minimal_glycolysis_conc_to_rates(sol.u[end], sol.prob.p).ATPprod,
                sol.prob.p.ATPase_Vmax;
                rtol = 0.01,
            )
                ATP_prod_eq_cons_frac += 1
            else
                ATP_prod_eq_cons_frac += 0
            end
            ATP_conc += sol.u[end].ATP / adenine_pool_size
            ATP_energy +=
                -log(minimal_glycolysis_conc_to_disequilibrium_ratios(sol.u[end], sol.prob.p).Q_Keq_ATPase)
            counter += 1
        end
    end
    ATP_prod_eq_cons_frac = ATP_prod_eq_cons_frac / counter
    ATP_conc = ATP_conc / counter
    ATP_energy = ATP_energy / counter
    Fract_success_sims = counter / n_Vmax_ATPase_values
    return (
        ATP_prod_eq_cons_frac = ATP_prod_eq_cons_frac,
        ATP_conc = ATP_conc,
        ATP_energy = ATP_energy,
        Fract_success_sims = Fract_success_sims,
        convert(NamedTuple, params)...,
    )
end

random_log_range(start_val, end_val) = 10^(start_val + (end_val - start_val) * rand())

function generate_random_params()
    rand_model_params = LVector(
        Enz1_L = random_log_range(-5, -1),
        Enz1_n = rand([1, 2, 3, 4]),
        # Enz1_L = 0.0,
        # Enz1_n = 1,
        Enz1_Vmax = random_log_range(-3, -1),
        Enz1_Keq = random_log_range(1, 4),
        Enz2_Vmax = random_log_range(-3, -1),
        Enz2_Keq = random_log_range(1, 4),
        ATPase_Vmax = 1e-4,
        ATPase_Keq = 1e3,
    )
    return rand_model_params
end

n_Vmax_ATPase_values = 100
n_repeats = 1_000
res = @time pmap(
    x -> minimal_model_output(x, minimal_model_initial_concentrations, n_Vmax_ATPase_values),
    [generate_random_params() for i = 1:n_repeats],
)

two_model_sims_res_df = DataFrame(res)
CSV.write(
    "$(Dates.format(now(),"mmddyy"))_hist_two_enzyme_model_fixed_Km_10mM_Pi_$(n_repeats)_repeats_fixed_Pi.csv",
    two_model_sims_res_df,
)

# Remove all the workers
for n in workers()
    rmprocs(n)
end

##
#Plot results

# Load simulation at a param range for histograms
two_model_sims_res_df = CSV.read("042423_hist_two_enzyme_model_fixed_Km_10mM_Pi_1000_repeats_fixed_Pi.csv", DataFrame)
# two_model_sims_res_df = CSV.read("$(Dates.format(now(),"mmddyy"))_hist_two_enzyme_model_fixed_Km_10mM_Pi_$(n_repeats)_repeats.csv", DataFrame)
filter!(df -> df.Fract_success_sims == 1.0, two_model_sims_res_df)

two_model_sims_res_df_no_allost =
    CSV.read("042423_hist_two_enzyme_model_fixed_Km_10mM_Pi_no_allost_10000_repeats_fixed_Pi.csv", DataFrame)
filter!(df -> df.Fract_success_sims == 1.0, two_model_sims_res_df_no_allost)

# Plot the results
size_inches = (6.5, 3)
size_pt = 72 .* size_inches
set_theme!(Theme(fontsize = 6, Axis = (
    xticksize = 1,
    yticksize = 1,
    # xticklabelsize = 5,
    # yticklabelsize = 5,
    yticklabelpad = 1,
    ylabelpadding = 3,
)))
fig = Figure(resolution = size_pt)


two_enzyme_schematic = load(
    "/Users/Denis/Library/Mobile Documents/com~apple~CloudDocs/My Articles/Glycolysis model paper/Figures/Two_enzyme_glycolysis_schematic_w_feedback_5mM_const_Pi.png",
)

ax_two_enzyme_schematic, im = image(
    fig[1, 1:3],
    rotr90(two_enzyme_schematic),
    axis = (aspect = DataAspect(), title = "Two-enzyme model with\nconstant [Phosphate]=5mM"),
)
ax_two_enzyme_schematic.alignmode = Mixed(top = -15, bottom = -5, left = -10)
ax_two_enzyme_schematic.width = 70
ax_two_enzyme_schematic.tellwidth = false

hidedecorations!(ax_two_enzyme_schematic)
hidespines!(ax_two_enzyme_schematic)

# Plot hist of ATP supply and demand match in Monte Carlo sims of two enzyme model
ax_hist = Axis(
    fig[1, 4:6],
    # limits = ((nothing), (0, 0.01)),
    # xscale = log10,
    title = "Ratio of ATP prod / ATPase of\ntwo enzyme model at a\nrange of model parameters",
    ylabel = "Fraction of simulations",
    xlabel = "ATP prod / ATPase",
    yticklabelrotation = π / 2,
)
complete_hist = stephist!(
    ax_hist,
    two_model_sims_res_df.ATP_prod_eq_cons_frac,
    bins = range(0, 1.01, 20),
    normalization = :probability,
    color = Makie.wong_colors()[1],
)
no_allo_hist = stephist!(
    ax_hist,
    two_model_sims_res_df_no_allost.ATP_prod_eq_cons_frac,
    bins = range(0, 1.01, 20),
    normalization = :probability,
    color = Makie.wong_colors()[6],
    linestyle = :dot,
)
axislegend(
    ax_hist,
    [complete_hist, no_allo_hist],
    ["+ allost.", "– allost."],
    position = :lt,
    # position = (1.075, 1.05),
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
    patchsize = (10, 5),
)

# Plot hist of ATP level in Monte Carlo sims of two enzyme model
ax_hist = Axis(
    fig[1, 7:9],
    limits = ((nothing), (0.0, 0.1)),
    # limits = ((nothing), (1e-3, 1)),
    # yscale = log10,
    title = "ATP levels of sims. of\ntwo-enzyme model at a\nrange of model parameters",
    ylabel = "Fraction of simulations",
    # xlabel = "[ATP] / ([ATP]+[ADP])",
    xlabel = L"\frac{[ATP]}{[ATP]+[ADP]}",
    yticklabelrotation = π / 2,
)
complete_hist = stephist!(
    ax_hist,
    two_model_sims_res_df.ATP_conc,
    bins = range(0, 1, 20),
    normalization = :probability,
    color = Makie.wong_colors()[1],
)
no_allo_hist = stephist!(
    ax_hist,
    two_model_sims_res_df_no_allost.ATP_conc,
    bins = range(0, 1, 20),
    normalization = :probability,
    color = Makie.wong_colors()[6],
    linestyle = :dot,
)
axislegend(
    ax_hist,
    [complete_hist, no_allo_hist],
    ["+ allost.", "– allost."],
    position = :rt,
    # position = (1.075, 1.05),
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
    patchsize = (10, 5),
)

# Plot hist of ATP energy in Monte Carlo sims of two enzyme model
ax_hist = Axis(
    fig[1, 10:12],
    # limits = ((nothing), (0, 0.01)),
    # xscale = log10,
    title = "Energy of ATP hydrolysis of\ntwo enzyme model at a\nrange of model parameters",
    ylabel = "Fraction of simulations",
    xlabel = rich("Energy of ATP hydrolysis, k", subscript("B"), "T"),
    yticklabelrotation = π / 2,
)
complete_hist = stephist!(
    ax_hist,
    two_model_sims_res_df.ATP_energy,
    bins = range(0, 1.1 * maximum(two_model_sims_res_df.ATP_energy), 20),
    normalization = :probability,
    color = Makie.wong_colors()[1],
)
no_allo_hist = stephist!(
    ax_hist,
    two_model_sims_res_df_no_allost.ATP_energy,
    bins = range(0, 1.1 * maximum(two_model_sims_res_df.ATP_energy), 20),
    normalization = :probability,
    color = Makie.wong_colors()[6],
    linestyle = :dot,
)
axislegend(
    ax_hist,
    [complete_hist, no_allo_hist],
    ["+ allost.", "– allost."],
    position = :lt,
    # position = (1.075, 1.05),
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
    patchsize = (10, 5),
)

colgap!(fig.layout, 10)
rowgap!(fig.layout, 5)


label_a = fig[1, 1, TopLeft()] = Label(fig, "A", fontsize = 12, halign = :right, padding = (0, 10, 10, 0))
label_b = fig[1, 4, TopLeft()] = Label(fig, "B", fontsize = 12, halign = :right, padding = (0, 10, 10, 0))
label_c = fig[1, 7, TopLeft()] = Label(fig, "C", fontsize = 12, halign = :right, padding = (0, 10, 10, 0))
label_d = fig[1, 10, TopLeft()] = Label(fig, "D", fontsize = 12, halign = :right, padding = (0, 10, 10, 0))
# label_h = fig[2, 4, TopLeft()] = Label(fig, "H", fontsize = 12, halign = :right, padding = (0, 10, 5, 0))

fig

# save("/Users/Denis/Library/Mobile Documents/com~apple~CloudDocs/Research Projects/Glycolysis Model/JuliaGlycolysisModel/Results data and figures/$(Dates.format(now(),"mmddyy"))_FigS7_two_enzyme_model_w_constant_pi.png", fig, px_per_unit = 4)
