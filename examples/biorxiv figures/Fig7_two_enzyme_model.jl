include("for_Fig7_TwoEnzymeModel.jl")
using OrdinaryDiffEq, CairoMakie, LabelledArrays
using DataFrames, CSV, Dates
using FileIO


function sol_at_ATPase_range(
    ODEs,
    minimal_model_params,
    minimal_model_initial_concentrations;
    n_Vmax_ATPase_values=1000,
)
    tspan = (0.0, 1e8)
    pathway_Vmax = min(minimal_model_params.Enz1_Vmax, minimal_model_params.Enz2_Vmax)
    # pathway_Vmax = minimal_model_params.Enz1_Vmax
    Adenine_pool = minimal_model_initial_concentrations.ATP + minimal_model_initial_concentrations.ADP
    ATPases = 10 .^ range(log10(0.001), log10(1.0), n_Vmax_ATPase_values) .* pathway_Vmax
    prob = ODEProblem(ODEs, minimal_model_initial_concentrations, tspan, minimal_model_params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
    sim = solve(
        ensemble_prob,
        Rodas5P(),
        EnsembleThreads(),
        trajectories=n_Vmax_ATPase_values,
        abstol=1e-18,
        reltol=1e-11,
        save_everystep=false,
        save_start=false,
    )
    ATP_energy =
        -log.([
            minimal_glycolysis_conc_to_disequilibrium_ratios(sol.u[end], sol.prob.p).Q_Keq_ATPase for
            sol in sim if sol.retcode == ReturnCode.Success
        ])
    ATP_conc = [sol.u[end].ATP for sol in sim if sol.retcode == ReturnCode.Success]
    X_conc = [sol.u[end].X for sol in sim if sol.retcode == ReturnCode.Success]
    ATPases_succ_sims = [ATPases[i] for (i, sol) in enumerate(sim) if sol.retcode == ReturnCode.Success]

    return (ATP_conc=ATP_conc, ATP_energy=ATP_energy, ATPase_Vmax=ATPases_succ_sims / pathway_Vmax, X_conc=X_conc)
end

# Calculate complete model results

optimal_minimal_model_params = @LVector Float64 propertynames(minimal_model_params)

optimal_minimal_model_params .=
    [0.000552043, 4, 0.014041193, 1453.992081, 0.022866635, 9567.335851, 0.0001, 1000]

complete_model = sol_at_ATPase_range(
    minimal_glycolysis_ODEs,
    optimal_minimal_model_params,
    minimal_model_initial_concentrations;
    n_Vmax_ATPase_values=100,
)

# Calculate no allostery model results
optimal_minimal_model_params_no_allostery = deepcopy(optimal_minimal_model_params)
optimal_minimal_model_params_no_allostery.Enz1_L = 0.0
optimal_minimal_model_params_no_allostery.Enz1_n = 1

model_no_allostery = sol_at_ATPase_range(
    minimal_glycolysis_ODEs,
    optimal_minimal_model_params_no_allostery,
    minimal_model_initial_concentrations;
    n_Vmax_ATPase_values=100,
)

model_no_allostery_fixed_Pi = sol_at_ATPase_range(
    minimal_glycolysis_ODEs_fixed_Pi,
    optimal_minimal_model_params,
    minimal_model_initial_concentrations;
    n_Vmax_ATPase_values=100,
)

model_no_allostery_Enz1_eq_ATPase = sol_at_ATPase_range(
    minimal_glycolysis_ODEs_Enz1_eq_ATPase,
    optimal_minimal_model_params,
    minimal_model_initial_concentrations;
    n_Vmax_ATPase_values=100,
)

# Calculate n effect on model results
n_value_model_params = deepcopy(optimal_minimal_model_params)
n_value_model = []
n_values = [4, 2, 1]
for n in n_values
    n_value_model_params.Enz1_n = n
    res = sol_at_ATPase_range(
        minimal_glycolysis_ODEs,
        n_value_model_params,
        minimal_model_initial_concentrations;
        n_Vmax_ATPase_values=100,
    )
    push!(n_value_model, res)
end

# Calculate effect of ENz 1 Keq on model results
Keq_value_model_params = deepcopy(optimal_minimal_model_params)
Keq_value_model_params = deepcopy(optimal_minimal_model_params)
Keq_value_model_params.Enz1_L = 0.0
Keq_value_model_params.Enz1_n = 1

Keq_value_model = []
Keq_values = [0.1, 10]
for Keq in Keq_values
    Keq_value_model_params.Enz1_Keq = Keq
    res = sol_at_ATPase_range(
        minimal_glycolysis_ODEs,
        Keq_value_model_params,
        minimal_model_initial_concentrations;
        n_Vmax_ATPase_values=100,
    )
    push!(Keq_value_model, res)
end



##
#Precalculate MC simulation of two enzyme model
#This code block takes ~ 1 hour to run 10,000 repeats on 8 core machine
using Distributed
addprocs(8; exeflags="--project")
@everywhere include("for_Fig7_TwoEnzymeModel.jl")
@everywhere using OrdinaryDiffEq

using DataFrames, CSV, Dates

@everywhere function minimal_model_output(params, minimal_ODEs, initial_concentrations, n_Vmax_ATPase_values)
    tspan = (0.0, 1e8)
    adenine_pool_size = initial_concentrations.ATP + initial_concentrations.ADP
    pathway_Vmax = min(params.Enz1_Vmax, params.Enz2_Vmax)

    ATPases = 10 .^ range(log10(0.01), log10(1.0), n_Vmax_ATPase_values) .* pathway_Vmax
    prob = ODEProblem(minimal_ODEs, initial_concentrations, tspan, params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
    sim = solve(
        ensemble_prob,
        Rodas5P(),
        EnsembleSerial(),
        trajectories=n_Vmax_ATPase_values,
        abstol=1e-15,
        reltol=1e-8,
        save_everystep=false,
        save_start=false,
        maxiters=1e4,
    )
    mean_ATP_prod_eq_cons_frac = 0.0
    mean_ATP_conc = 0.0
    mean_ATP_energy = 0.0
    mean_X_conc = 0.0
    counter = 0
    for sol in sim
        if sol.retcode == ReturnCode.Success && all(sol.u[end] .> 0)
            if isapprox(
                minimal_glycolysis_conc_to_rates(sol.u[end], sol.prob.p).ATPprod,
                sol.prob.p.ATPase_Vmax;
                rtol=0.01,
            )
                mean_ATP_prod_eq_cons_frac += 1
            else
                mean_ATP_prod_eq_cons_frac += 0
            end
            mean_ATP_conc += sol.u[end].ATP / adenine_pool_size
            mean_ATP_energy +=
                -log(minimal_glycolysis_conc_to_disequilibrium_ratios(sol.u[end], sol.prob.p).Q_Keq_ATPase)
            mean_X_conc += sol.u[end].X
            counter += 1
        end
    end
    mean_ATP_prod_eq_cons_frac = mean_ATP_prod_eq_cons_frac / counter
    mean_ATP_conc = mean_ATP_conc / counter
    mean_ATP_energy = mean_ATP_energy / counter
    mean_X_conc = mean_X_conc / counter
    Fract_success_sims = counter / n_Vmax_ATPase_values
    return (
        mean_ATP_prod_eq_cons_frac=mean_ATP_prod_eq_cons_frac,
        mean_ATP_conc=mean_ATP_conc,
        mean_ATP_energy=mean_ATP_energy,
        mean_X_conc=mean_X_conc,
        Fract_success_sims=Fract_success_sims,
        convert(NamedTuple, params)...,
    )
end

random_log_range(start_val, end_val) = 10^(start_val + (end_val - start_val) * rand())

function generate_random_params()
    rand_model_params = LVector(
        # Enz1_L = random_log_range(-5, -1), #uncomment to include allostery
        # Enz1_n = rand([1, 2, 3, 4]), #uncomment to include allostery
        Enz1_L=0.0, #uncomment to fix L=0 and remove allostery
        Enz1_n=1, #uncomment to fix n=1 and remove allostery
        Enz1_Vmax=random_log_range(-3, -1),
        Enz1_Keq=random_log_range(1, 4),
        Enz2_Vmax=random_log_range(-3, -1),
        Enz2_Keq=random_log_range(1, 4),
        ATPase_Vmax=1e-4,
        ATPase_Keq=1e3,
    )
    return rand_model_params
end
n_Vmax_ATPase_values = 100
n_repeats = 1_000
res = @time pmap(
    x -> minimal_model_output(x, minimal_glycolysis_ODEs, minimal_model_initial_concentrations, n_Vmax_ATPase_values),
    [generate_random_params() for i = 1:n_repeats],
)
vscodedisplay(two_model_sims_res_df)
two_model_sims_res_df = DataFrame(res)
# CSV.write(
#     "Results/$(Dates.format(now(),"mmddyy"))_hist_two_enzyme_model_$(n_repeats)_repeats_no_allost_const_Pi.csv",
#     two_model_sims_res_df,
# )
# CSV.write(
#     "Results/$(Dates.format(now(),"mmddyy"))_hist_two_enzyme_model_$(n_repeats)_repeats_w_allost_const_Pi.csv",
#     two_model_sims_res_df,
# )
# CSV.write(
#     "Results/$(Dates.format(now(),"mmddyy"))_hist_two_enzyme_model_$(n_repeats)_repeats_w_allost_variable_Pi.csv",
#     two_model_sims_res_df,
# )
CSV.write(
    "Results/$(Dates.format(now(),"mmddyy"))_hist_two_enzyme_model_$(n_repeats)_repeats_no_allost_variable_Pi.csv",
    two_model_sims_res_df,
)

# Remove all the workers
for n in workers()
    rmprocs(n)
end

##
#Plot results

# Load simulation at a param range for histograms
two_model_sims_res_df = CSV.read("Results/071524_hist_two_enzyme_model_10000_repeats_w_allost_variable_Pi.csv", DataFrame)
# two_model_sims_res_df = CSV.read("$(Dates.format(now(),"mmddyy"))_hist_two_enzyme_model_$(n_repeats)_repeats.csv", DataFrame)
filter!(df -> df.Fract_success_sims == 1.0, two_model_sims_res_df)

two_model_sims_res_df_no_allost_fixed_Pi =
    CSV.read("Results/071524_hist_two_enzyme_model_10000_repeats_no_allost_const_Pi.csv", DataFrame)
filter!(df -> df.Fract_success_sims == 1.0, two_model_sims_res_df_no_allost_fixed_Pi)
two_model_sims_res_df_no_allost =
    CSV.read("Results/071524_hist_two_enzyme_model_10000_repeats_no_allost_variable_Pi.csv", DataFrame)
filter!(df -> df.Fract_success_sims == 1.0, two_model_sims_res_df_no_allost)

# Plot the results
size_inches = (6.5, 5)
size_pt = 72 .* size_inches
set_theme!(Theme(fontsize=6, Axis=(
    xticksize=1,
    yticksize=1,
    # xticklabelsize = 5,
    # yticklabelsize = 5,
    yticklabelpad=1,
    ylabelpadding=3,
)))
fig = Figure(size=size_pt)


two_enzyme_schematic = load(
    "111523_Two_enzyme_glycolysis_schematic_w_feedback_w_pi_trap.png",
)

ax_two_enzyme_schematic, im = image(
    fig[1, 1:3],
    rotr90(two_enzyme_schematic),
    axis=(aspect=DataAspect(), title="Two-enzyme model"),
)
ax_two_enzyme_schematic.alignmode = Mixed(top=-15, bottom=-5, left=-10)
ax_two_enzyme_schematic.width = 80
ax_two_enzyme_schematic.tellwidth = false

hidedecorations!(ax_two_enzyme_schematic)
hidespines!(ax_two_enzyme_schematic)

# Plot hist of ATP supply and demand match in Monte Carlo sims of two enzyme model
ax_hist = Axis(
    fig[1, 4:6],
    # limits = ((nothing), (0, 0.01)),
    # xscale = log10,
    title="Ratio of ATP prod / ATPase of\ntwo enzyme model at a\nrange of model parameters",
    ylabel="Fraction of simulations",
    xlabel="ATP prod / ATPase",
    yticklabelrotation=π / 2,
)
complete_hist = stephist!(
    ax_hist,
    two_model_sims_res_df.mean_ATP_prod_eq_cons_frac,
    bins=range(0, 1.01, 20),
    normalization=:probability,
    color=Makie.wong_colors()[1],
)
no_allo_hist = stephist!(
    ax_hist,
    two_model_sims_res_df_no_allost.mean_ATP_prod_eq_cons_frac,
    bins=range(0, 1.01, 20),
    normalization=:probability,
    color=Makie.wong_colors()[6],
    linestyle=:dot,
)
no_allo_fixed_Pi_hist = stephist!(
    ax_hist,
    two_model_sims_res_df_no_allost_fixed_Pi.mean_ATP_prod_eq_cons_frac,
    bins=range(0, 1.01, 20),
    normalization=:probability,
    color=Makie.wong_colors()[3],
    linestyle=:dash,
)
axislegend(
    ax_hist,
    [complete_hist, no_allo_hist, no_allo_fixed_Pi_hist],
    ["+ allost.", "– allost.", "– allost.\n[Pi]=5mM"],
    position=(-0.3, 1.1),
    rowgap=1,
    framevisible=false,
    patchsize=(10, 5),
)

# Plot hist of ATP level in Monte Carlo sims of two enzyme model
ax_hist = Axis(
    fig[1, 7:9],
    limits=((nothing), (0.0, 0.1)),
    # limits = ((nothing), (1e-3, 1)),
    # yscale = log10,
    title="ATP levels of sims. of\ntwo-enzyme model at a\nrange of model parameters",
    ylabel="Fraction of simulations",
    xlabel="[ATP] / ([ATP]+[ADP])",
    # xlabel = L"\frac{[ATP]}{[ATP]+[ADP]}",
    yticklabelrotation=π / 2,
)
complete_hist = stephist!(
    ax_hist,
    two_model_sims_res_df.mean_ATP_conc,
    bins=range(0, 1, 20),
    normalization=:probability,
    color=Makie.wong_colors()[1],
)
no_allo_hist = stephist!(
    ax_hist,
    two_model_sims_res_df_no_allost.mean_ATP_conc,
    bins=range(0, 1, 20),
    normalization=:probability,
    color=Makie.wong_colors()[6],
    linestyle=:dot,
)
no_allo_fixed_Pi_hist = stephist!(
    ax_hist,
    two_model_sims_res_df_no_allost_fixed_Pi.mean_ATP_conc,
    bins=range(0, 1, 20),
    normalization=:probability,
    color=Makie.wong_colors()[3],
    linestyle=:dash,
)
axislegend(
    ax_hist,
    [complete_hist, no_allo_hist, no_allo_fixed_Pi_hist],
    ["+ allost.", "– allost.", "– allost.\n[Pi]=5mM"],
    position=(0.5, 1.1),
    rowgap=1,
    framevisible=false,
    patchsize=(10, 5),
)

# Plot hist of ATP energy in Monte Carlo sims of two enzyme model
ax_hist = Axis(
    fig[1, 10:12],
    # limits = ((nothing), (0, 0.01)),
    # xscale = log10,
    title="Energy of ATP hydrolysis of\ntwo enzyme model at a\nrange of model parameters",
    ylabel="Fraction of simulations",
    xlabel=rich("Energy of ATP hydrolysis, k", subscript("B"), "T"),
    yticklabelrotation=π / 2,
)
complete_hist = stephist!(
    ax_hist,
    two_model_sims_res_df.mean_ATP_energy,
    bins=range(0, 1.1 * maximum(two_model_sims_res_df.mean_ATP_energy), 20),
    normalization=:probability,
    color=Makie.wong_colors()[1],
)
no_allo_hist = stephist!(
    ax_hist,
    two_model_sims_res_df_no_allost.mean_ATP_energy,
    bins=range(0, 1.1 * maximum(two_model_sims_res_df.mean_ATP_energy), 20),
    normalization=:probability,
    color=Makie.wong_colors()[6],
    linestyle=:dot,
)
no_allo_fixed_Pi_hist = stephist!(
    ax_hist,
    two_model_sims_res_df_no_allost_fixed_Pi.mean_ATP_energy,
    bins=range(0, 1.1 * maximum(two_model_sims_res_df.mean_ATP_energy), 20),
    normalization=:probability,
    color=Makie.wong_colors()[3],
    linestyle=:dash,
)
axislegend(
    ax_hist,
    [complete_hist, no_allo_hist, no_allo_fixed_Pi_hist],
    ["+ allost.", "– allost.", "– allost.\n[Pi]=5mM"],
    position=(-0.3, 1.1),
    rowgap=1,
    framevisible=false,
    patchsize=(10, 5),
)


# Plot [ATP] maintenance depending on allostery
ax_ATPase_range = Axis(
    fig[2, 1:3],
    limits=((0.001, 1.2), (-0.5e-3, 11e-3)),
    xlabel="ATPase, % of pathway Vmax",
    ylabel="[ATP],mM",
    title="Effect of allostery on\ntwo-enzyme model",
    xscale=log10,
    # yscale = log10,
    xticks=[0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0],
    # xtickformat = xs -> ["$(Int(round(x*100)))%" for x in xs],
    xtickformat=xs -> ["$(100*x < 1.0 ? round(x*100, sigdigits=1) : Int(round(x*100)))" for x in xs],
    ytickformat=ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys],
    # yticklabelcolor = Makie.wong_colors()[1],
    # ylabelcolor = Makie.wong_colors()[1],
)
for (i, n) in enumerate(n_values)
    lines!(
        ax_ATPase_range,
        n_value_model[i].ATPase_Vmax,
        n_value_model[i].ATP_conc,
        color=Makie.wong_colors()[2+i],
        label="+ allost. n=$(n)",
    )
end
ATP_line_no_reg = lines!(
    ax_ATPase_range,
    model_no_allostery.ATPase_Vmax,
    model_no_allostery.ATP_conc,
    color=Makie.wong_colors()[6],
    # linestyle = :dot,
    label="– allost. L=0",
)
lines!(
    ax_ATPase_range,
    [0.0009, 1.1],
    repeat([minimal_model_initial_concentrations.ATP + minimal_model_initial_concentrations.ADP], 2),
    color=:grey,
    linestyle=:dash,
)
text!(
    ax_ATPase_range,
    0.011,
    1.01 * (minimal_model_initial_concentrations.ATP + minimal_model_initial_concentrations.ADP),
    text="Adenine pool size",
    align=(:left, :bottom),
    color=:grey,
)
axislegend(
    ax_ATPase_range,
    position=(1.4, -0.01),
    rowgap=0,
    framevisible=false,
    patchsize=(10, 5),
)

# Plot reversibility of constant Pi
ax_ATPase_range = Axis(
    fig[2, 4:6],
    limits=((0.001, 1.2), (-0.5e-3, 11e-3)),
    xlabel="ATPase, % of pathway Vmax",
    ylabel="[ATP],mM",
    title="Two-enzyme model with\nconstant phosphate",
    xscale=log10,
    # yscale = log10,
    xticks=[0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0],
    # xtickformat = xs -> ["$(Int(round(x*100)))%" for x in xs],
    xtickformat=xs -> ["$(100*x < 1.0 ? round(x*100, sigdigits=1) : Int(round(x*100)))" for x in xs],
    ytickformat=ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys],
    # yticklabelcolor = Makie.wong_colors()[1],
    # ylabelcolor = Makie.wong_colors()[1],
)

ATP_line_no_reg = lines!(
    ax_ATPase_range,
    model_no_allostery.ATPase_Vmax,
    model_no_allostery.ATP_conc,
    color=Makie.wong_colors()[6],
    # linestyle = :dot,
    label="– allost. L=0",
)
ATP_line_no_reg = lines!(
    ax_ATPase_range,
    model_no_allostery_fixed_Pi.ATPase_Vmax,
    model_no_allostery_fixed_Pi.ATP_conc,
    color=Makie.wong_colors()[3],
    # linestyle = :dot,
    label="– allost. L=0\nconst. [Pi]=5mM",
)
lines!(
    ax_ATPase_range,
    [0.0009, 1.1],
    repeat([minimal_model_initial_concentrations.ATP + minimal_model_initial_concentrations.ADP], 2),
    color=:grey,
    linestyle=:dash,
)
text!(
    ax_ATPase_range,
    0.011,
    1.01 * (minimal_model_initial_concentrations.ATP + minimal_model_initial_concentrations.ADP),
    text="Adenine pool size",
    align=(:left, :bottom),
    color=:grey,
)
axislegend(
    ax_ATPase_range,
    position=(0.0, -0.01),
    rowgap=3,
    framevisible=false,
    patchsize=(10, 5),
)

# Plot effect of Keq
ax_ATPase_range = Axis(
    fig[2, 7:9],
    limits=((0.001, 1.2), (-0.5e-3, 11e-3)),
    xlabel="ATPase, % of pathway Vmax",
    ylabel="[ATP],mM",
    title="Two-enzyme model with\ndifferent Keq of Enz. 1",
    xscale=log10,
    # yscale = log10,
    xticks=[0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0],
    # xtickformat = xs -> ["$(Int(round(x*100)))%" for x in xs],
    xtickformat=xs -> ["$(100*x < 1.0 ? round(x*100, sigdigits=1) : Int(round(x*100)))" for x in xs],
    ytickformat=ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys],
    # yticklabelcolor = Makie.wong_colors()[1],
    # ylabelcolor = Makie.wong_colors()[1],
)

for (i, Keq) in enumerate(Keq_values)
    lines!(
        ax_ATPase_range,
        Keq_value_model[i].ATPase_Vmax,
        Keq_value_model[i].ATP_conc,
        color=Makie.wong_colors()[2+i],
        label="- allost. L=0\nEnz. 1 Keq=$(Keq)",
    )
end
lines!(
    ax_ATPase_range,
    model_no_allostery.ATPase_Vmax,
    model_no_allostery.ATP_conc,
    color=Makie.wong_colors()[6],
    # linestyle = :dot,
    label="– allost. L=0\nEnz. 1 Keq=$(round(optimal_minimal_model_params.Enz1_Keq))",
)
lines!(
    ax_ATPase_range,
    [0.0009, 1.1],
    repeat([minimal_model_initial_concentrations.ATP + minimal_model_initial_concentrations.ADP], 2),
    color=:grey,
    linestyle=:dash,
)
text!(
    ax_ATPase_range,
    0.011,
    1.01 * (minimal_model_initial_concentrations.ATP + minimal_model_initial_concentrations.ADP),
    text="Adenine pool size",
    align=(:left, :bottom),
    color=:grey,
)
axislegend(
    ax_ATPase_range,
    position=(0.7, 0.3),
    rowgap=3,
    framevisible=false,
    patchsize=(10, 5),
)

# Plot effect of V(Enz1)=V(ATPase)
ax_ATPase_range = Axis(
    fig[2, 10:12],
    limits=((0.001, 1.2), (-0.5e-3, 11e-3)),
    xlabel="ATPase, % of pathway Vmax",
    ylabel="[ATP],mM",
    title="Two-enzyme model with\nV(Enz.1)=V(ATPase)",
    xscale=log10,
    # yscale = log10,
    xticks=[0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0],
    # xtickformat = xs -> ["$(Int(round(x*100)))%" for x in xs],
    xtickformat=xs -> ["$(100*x < 1.0 ? round(x*100, sigdigits=1) : Int(round(x*100)))" for x in xs],
    ytickformat=ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys],
    # yticklabelcolor = Makie.wong_colors()[1],
    # ylabelcolor = Makie.wong_colors()[1],
)
lines!(
    ax_ATPase_range,
    model_no_allostery.ATPase_Vmax,
    model_no_allostery.ATP_conc,
    color=Makie.wong_colors()[6],
    # linestyle = :dot,
    label="– allost. L=0",
)
lines!(
    ax_ATPase_range,
    model_no_allostery_Enz1_eq_ATPase.ATPase_Vmax,
    model_no_allostery_Enz1_eq_ATPase.ATP_conc,
    color=Makie.wong_colors()[3],
    # linestyle = :dot,
    label="– allost. L=0\nV(Enz. 1)=V(ATPase)",
)
lines!(
    ax_ATPase_range,
    [0.0009, 1.1],
    repeat([minimal_model_initial_concentrations.ATP + minimal_model_initial_concentrations.ADP], 2),
    color=:grey,
    linestyle=:dash,
)
text!(
    ax_ATPase_range,
    0.011,
    1.01 * (minimal_model_initial_concentrations.ATP + minimal_model_initial_concentrations.ADP),
    text="Adenine pool size",
    align=(:left, :bottom),
    color=:grey,
)
axislegend(
    ax_ATPase_range,
    position=(0.5, 0.1),
    rowgap=3,
    framevisible=false,
    patchsize=(10, 5),
)

colgap!(fig.layout, 10)
rowgap!(fig.layout, 5)


label_a = fig[1, 1, TopLeft()] = Label(fig, "A", fontsize=12, halign=:right, padding=(0, 10, 10, 0))
label_b = fig[1, 4, TopLeft()] = Label(fig, "B", fontsize=12, halign=:right, padding=(0, 10, 10, 0))
label_c = fig[1, 7, TopLeft()] = Label(fig, "C", fontsize=12, halign=:right, padding=(0, 10, 10, 0))
label_d = fig[1, 10, TopLeft()] = Label(fig, "D", fontsize=12, halign=:right, padding=(0, 10, 10, 0))
label_e = fig[2, 1, TopLeft()] = Label(fig, "E", fontsize=12, halign=:right, padding=(0, 10, 5, 0))
label_f = fig[2, 4, TopLeft()] = Label(fig, "F", fontsize=12, halign=:right, padding=(0, 10, 5, 0))
label_g = fig[2, 7, TopLeft()] = Label(fig, "G", fontsize=12, halign=:right, padding=(0, 10, 5, 0))
label_h = fig[2, 10, TopLeft()] = Label(fig, "H", fontsize=12, halign=:right, padding=(0, 10, 5, 0))

fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_Fig7_two_enzyme_model.png", fig, px_per_unit = 4)
