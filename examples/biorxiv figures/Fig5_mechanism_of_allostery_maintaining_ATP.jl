using Glycolysis
using OrdinaryDiffEq, DiffEqCallbacks, ProgressMeter, LabelledArrays
using CairoMakie, Dates, Printf, Statistics, StatsBase
using DataFrames, CSV
using FileIO

##
# Precalculate output of complete model and model without regulation
# This code takes ~20 minutes to run on an 8 core machine
using Distributed
addprocs()

@everywhere using OrdinaryDiffEq, Glycolysis

start = now()
@everywhere function glycolysis_ODEs_fixed_Pi(ds, s, params, t)
    ds.Glucose_media = 0.0
    ds.Glucose = Glycolysis.rate_GLUT(s, params) -
                 Glycolysis.rate_HK1(s, params)
    ds.G6P = Glycolysis.rate_HK1(s, params) -
             Glycolysis.rate_GPI(s, params)
    ds.F6P = (
        Glycolysis.rate_GPI(s, params) -
        Glycolysis.rate_PFKP(s, params)
    )
    ds.F16BP = (
        Glycolysis.rate_PFKP(s, params) -
        Glycolysis.rate_ALDO(s, params)
    )
    ds.GAP = (
        Glycolysis.rate_ALDO(s, params) + Glycolysis.rate_TPI(s, params) -
        Glycolysis.rate_GAPDH(s, params)
    )
    ds.DHAP = Glycolysis.rate_ALDO(s, params) - Glycolysis.rate_TPI(s, params)
    ds.BPG = Glycolysis.rate_GAPDH(s, params) -
             Glycolysis.rate_PGK(s, params)
    ds.ThreePG = Glycolysis.rate_PGK(s, params) -
                 Glycolysis.rate_PGM(s, params)
    ds.TwoPG = Glycolysis.rate_PGM(s, params) - Glycolysis.rate_ENO(s, params)
    ds.PEP = Glycolysis.rate_ENO(s, params) -
             Glycolysis.rate_PKM2(s, params)
    ds.Pyruvate = Glycolysis.rate_PKM2(s, params) -
                  Glycolysis.rate_LDH(s, params)
    ds.Lactate = Glycolysis.rate_LDH(s, params) -
                 Glycolysis.rate_MCT(s, params)
    ds.Lactate_media = 0.0
    ds.ATP = (
        -Glycolysis.rate_HK1(s, params) -
        Glycolysis.rate_PFKP(s, params) +
        Glycolysis.rate_PGK(s, params) +
        Glycolysis.rate_PKM2(s, params) -
        Glycolysis.rate_ATPase(s, params) +
        Glycolysis.rate_AK(s, params)
    )
    ds.ADP = (
        Glycolysis.rate_HK1(s, params) +
        Glycolysis.rate_PFKP(s, params) -
        Glycolysis.rate_PGK(s, params) -
        Glycolysis.rate_PKM2(s, params) +
        Glycolysis.rate_ATPase(s, params) -
        2 * Glycolysis.rate_AK(s, params)
    )
    ds.AMP = Glycolysis.rate_AK(s, params)
    # ds.Phosphate =
    #     Glycolysis.rate_ATPase(s, params) -
    #     Glycolysis.rate_GAPDH(s, params)
    ds.Phosphate = 0.0
    ds.NAD = Glycolysis.rate_LDH(s, params) -
             Glycolysis.rate_GAPDH(s, params)
    ds.NADH = Glycolysis.rate_GAPDH(s, params) -
              Glycolysis.rate_LDH(s, params)
    ds.F26BP = 0.0
    ds.Citrate = 0.0
    ds.Phenylalanine = 0.0
end

@everywhere function glycolysis_ODEs_Vmax_HK_PFK_eq_ATPase(ds, s, params, t)
    # 0.5 * params.ATPase_Vmax * (1 - s.G6P * s.ADP / (params.HK1_Keq * s.Glucose * s.ATP))
    # 0.5 * params.ATPase_Vmax * (1 - s.F16BP * s.ADP / (params.PFKP_Keq * s.F6P * s.ATP))
    # params.ATPase_Vmax * (1 - s.Phosphate * s.ADP / (params.ATPase_Keq * s.ATP))
    # Glycolysis.rate_HK1(s, params)
    # Glycolysis.rate_PFKP(s, params)
    # Glycolysis.rate_ATPase(s, params)
    ds.Glucose_media = 0.0
    ds.Glucose = Glycolysis.rate_GLUT(s, params) -
                 0.5 * Glycolysis.rate_ATPase(s, params)
    ds.G6P = 0.5 * Glycolysis.rate_ATPase(s, params) -
             Glycolysis.rate_GPI(s, params)
    ds.F6P = (
        Glycolysis.rate_GPI(s, params) -
        0.5 * Glycolysis.rate_ATPase(s, params)
    )
    ds.F16BP = (
        0.5 * Glycolysis.rate_ATPase(s, params) -
        Glycolysis.rate_ALDO(s, params)
    )
    ds.GAP = (
        Glycolysis.rate_ALDO(s, params) + Glycolysis.rate_TPI(s, params) -
        Glycolysis.rate_GAPDH(s, params)
    )
    ds.DHAP = Glycolysis.rate_ALDO(s, params) - Glycolysis.rate_TPI(s, params)
    ds.BPG = Glycolysis.rate_GAPDH(s, params) -
             Glycolysis.rate_PGK(s, params)
    ds.ThreePG = Glycolysis.rate_PGK(s, params) -
                 Glycolysis.rate_PGM(s, params)
    ds.TwoPG = Glycolysis.rate_PGM(s, params) - Glycolysis.rate_ENO(s, params)
    ds.PEP = Glycolysis.rate_ENO(s, params) -
             Glycolysis.rate_PKM2(s, params)
    ds.Pyruvate = Glycolysis.rate_PKM2(s, params) -
                  Glycolysis.rate_LDH(s, params)
    ds.Lactate = Glycolysis.rate_LDH(s, params) -
                 Glycolysis.rate_MCT(s, params)
    ds.Lactate_media = 0.0
    ds.ATP = (
        -0.5 * Glycolysis.rate_ATPase(s, params) -
        0.5 * Glycolysis.rate_ATPase(s, params) +
        Glycolysis.rate_PGK(s, params) +
        Glycolysis.rate_PKM2(s, params) -
        Glycolysis.rate_ATPase(s, params) +
        Glycolysis.rate_AK(s, params)
    )
    ds.ADP = (
        0.5 * Glycolysis.rate_ATPase(s, params) +
        0.5 * Glycolysis.rate_ATPase(s, params) -
        Glycolysis.rate_PGK(s, params) -
        Glycolysis.rate_PKM2(s, params) +
        Glycolysis.rate_ATPase(s, params) -
        2 * Glycolysis.rate_AK(s, params)
    )
    ds.AMP = Glycolysis.rate_AK(s, params)
    ds.Phosphate =
        Glycolysis.rate_ATPase(s, params) -
        Glycolysis.rate_GAPDH(s, params)
    ds.NAD = Glycolysis.rate_LDH(s, params) -
             Glycolysis.rate_GAPDH(s, params)
    ds.NADH = Glycolysis.rate_GAPDH(s, params) -
              Glycolysis.rate_LDH(s, params)
    ds.F26BP = 0.0
    ds.Citrate = 0.0
    ds.Phenylalanine = 0.0
end

function find_ATP_at_ATPase_range(
    ODE_func,
    params,
    init_conc;
    n_Vmax_ATPase_values=500,
    min_ATPase=0.001,
    max_ATPase=1.0
)
    tspan = (0.0, 1e8)
    pathway_Vmax = 2 * params.HK1_Vmax * params.HK1_Conc
    ATPases = 10 .^ range(log10(min_ATPase), log10(max_ATPase), n_Vmax_ATPase_values) .*
              pathway_Vmax

    prob = ODEProblem(ODE_func, init_conc, tspan, params)
    function prob_func(prob, i, repeat)
        prob.p.ATPase_Vmax = ATPases[i]
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
    @time sim = solve(
        ensemble_prob,
        RadauIIA9(autodiff=false),
        EnsembleDistributed(),
        trajectories=length(ATPases),
        abstol=1e-15,
        reltol=1e-8,
        save_everystep=false,
        save_start=false,
        maxiters=1e3
    )
    ATP_conc = [sol.u[end].ATP for sol in sim if sol.retcode == ReturnCode.Success]
    ATPase_Vmax = [ATPases[i]
                   for (i, sol) in enumerate(sim) if sol.retcode == ReturnCode.Success] ./
                  (2 * params.HK1_Vmax * params.HK1_Conc)
    ATP_energy = [-log(Glycolysis.conc_to_disequilibrium_ratios(
        sol.u[end], params).Q_Keq_ATPase)
                  for
                  sol in sim if sol.retcode == ReturnCode.Success]
    return (ATP_conc=ATP_conc, ATPase_Vmax=ATPase_Vmax, ATP_energy=ATP_energy, sim=sim)
end

glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
glycolysis_params.ATPase_Km_ATP = 1e-9
Complete_Model_Simulation_Data = find_ATP_at_ATPase_range(
    glycolysis_ODEs, BigFloat.(glycolysis_params), BigFloat.(glycolysis_init_conc_copy))

no_reg_params = deepcopy(glycolysis_params)
no_reg_params.HK1_K_a_G6P_cat = Inf
no_reg_params.HK1_K_i_G6P_reg = Inf
no_reg_params.HK1_K_a_Pi = Inf
no_reg_params.PFKP_L = 0.0
no_reg_params.GAPDH_L = 0.0
no_reg_params.PKM2_L = 0.0
no_reg_params.ATPase_Km_ATP = 1e-9

glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
No_Reg_Model_Simulation_Data = find_ATP_at_ATPase_range(
    glycolysis_ODEs, BigFloat.(no_reg_params), BigFloat.(glycolysis_init_conc_copy))

glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
glycolysis_init_conc_copy.Phosphate = 1e-3
Const_Pi_No_Reg_Model_Simulation_Data = find_ATP_at_ATPase_range(
    glycolysis_ODEs_fixed_Pi, BigFloat.(no_reg_params), BigFloat.(glycolysis_init_conc_copy))

no_reg_Keq_params = deepcopy(glycolysis_params)
no_reg_Keq_params.HK1_K_a_G6P_cat = Inf
no_reg_Keq_params.HK1_K_i_G6P_reg = Inf
no_reg_Keq_params.HK1_K_a_Pi = Inf
no_reg_Keq_params.PFKP_L = 0.0
no_reg_Keq_params.GAPDH_L = 0.0
no_reg_Keq_params.PKM2_L = 0.0
no_reg_Keq_params.ATPase_Km_ATP = 1e-9
no_reg_Keq_params.HK1_Keq = 0.05
no_reg_Keq_params.PFKP_Keq = 0.05
no_reg_Keq_params.PGK_Keq = 5.0
no_reg_Keq_params.GAPDH_Keq = 0.5

Keq_HK_PFK_GAPDH_PGK_No_Reg_Model_Simulation_Data = find_ATP_at_ATPase_range(
    glycolysis_ODEs,
    BigFloat.(no_reg_Keq_params),
    BigFloat.(glycolysis_init_conc)
)

# glycolysis_init_conc_copy = deepcopy(BigFloat.(glycolysis_init_conc))
# no_reg_params = deepcopy(BigFloat.(no_reg_params))
glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
no_reg_params = deepcopy(no_reg_params)
Vmax_HK_PFK_eq_ATPase_No_Reg_Model_Simulation_Data = find_ATP_at_ATPase_range(
    glycolysis_ODEs_Vmax_HK_PFK_eq_ATPase, BigFloat.(no_reg_params), BigFloat.(glycolysis_init_conc_copy))


elapsed_minutes = round(now() - start, Minute)
println("Duration: ", elapsed_minutes)
#free workers
rmprocs(workers())

##
# Precalculate dynamic change in [Metabolite] after allostery removal
# calculate results without const Pi
glycolysis_params_copy = deepcopy(glycolysis_params)
no_reg_params = deepcopy(glycolysis_params)
no_reg_params.HK1_K_a_G6P_cat = Inf
no_reg_params.HK1_K_i_G6P_reg = Inf
no_reg_params.HK1_K_a_Pi = Inf
no_reg_params.PFKP_L = 0.0
no_reg_params.GAPDH_L = 0.0
no_reg_params.PKM2_L = 0.0
no_reg_params.ATPase_Km_ATP = 1e-9

tspan = (0.0, 120.0)
callback_time1 = 60
Initial_ATPase_Vmax_frac = 0.10
glycolysis_params_copy.ATPase_Vmax = Initial_ATPase_Vmax_frac * 2 *
                                     glycolysis_params_copy.HK1_Conc *
                                     glycolysis_params_copy.HK1_Vmax
no_reg_params.ATPase_Vmax = Initial_ATPase_Vmax_frac * 2 * no_reg_params.HK1_Conc *
                            no_reg_params.HK1_Vmax
function affect1!(integrator)
    integrator.p .= no_reg_params
    # integrator.p.HK1_K_a_G6P_cat = Inf
    # integrator.p.HK1_K_i_G6P_reg = Inf
    # integrator.p.PFKP_L = 0.0
    # integrator.p.GAPDH_L = 0.0
    # integrator.p.PKM2_L = 0.0
    # integrator.p.ATPase_Km_ATP = 1e-9
end

PresetTime_cb1 = PresetTimeCallback(callback_time1, affect1!)
init_cond_prob = ODEProblem(
    glycolysis_ODEs, glycolysis_init_conc, (0, 1e8), glycolysis_params_copy)
init_cond_sol = solve(
    init_cond_prob, RadauIIA9(), abstol=1e-15, reltol=1e-8, save_everystep=false)
new_init_cond = init_cond_sol.u[end]
prob = ODEProblem(glycolysis_ODEs, BigFloat.(new_init_cond), tspan,
    BigFloat.(glycolysis_params_copy), callback=PresetTime_cb1)
sol = solve(
    prob,
    RadauIIA9(autodiff=false),
    abstol=1e-15,
    reltol=1e-8,
    saveat=[k for k in tspan[1]:((tspan[2]-tspan[1])/10_000):tspan[2]]
)

# calculate results with const Pi
glycolysis_init_conc_copy = deepcopy(glycolysis_init_conc)
glycolysis_init_conc_copy.Phosphate = 1e-3
glycolysis_params_copy = deepcopy(glycolysis_params)
no_reg_params = deepcopy(glycolysis_params)
no_reg_params.HK1_K_a_G6P_cat = Inf
no_reg_params.HK1_K_i_G6P_reg = Inf
no_reg_params.HK1_K_a_Pi = Inf
no_reg_params.PFKP_L = 0.0
no_reg_params.GAPDH_L = 0.0
no_reg_params.PKM2_L = 0.0
no_reg_params.ATPase_Km_ATP = 1e-9

glycolysis_params_copy.ATPase_Vmax = Initial_ATPase_Vmax_frac * 2 *
                                     glycolysis_params_copy.HK1_Conc *
                                     glycolysis_params_copy.HK1_Vmax
no_reg_params.ATPase_Vmax = Initial_ATPase_Vmax_frac * 2 * no_reg_params.HK1_Conc *
                            no_reg_params.HK1_Vmax
function affect1!(integrator)
    integrator.p .= no_reg_params
end
PresetTime_cb1 = PresetTimeCallback(callback_time1, affect1!)

init_cond_const_pi_prob = ODEProblem(
    glycolysis_ODEs_fixed_Pi, BigFloat.(glycolysis_init_conc_copy), (0, 1e8), BigFloat.(glycolysis_params_copy))
init_cond_const_pi_sol = solve(
    init_cond_const_pi_prob, RadauIIA9(autodiff=false), abstol=1e-15, reltol=1e-8, save_everystep=false)
new_init_const_pi_cond = BigFloat.(init_cond_const_pi_sol.u[end])
tspan = (0.0, 180.0)
prob_const_pi = ODEProblem(glycolysis_ODEs_fixed_Pi, BigFloat.(new_init_const_pi_cond), tspan,
    BigFloat.(glycolysis_params_copy), callback=PresetTime_cb1)

sol_const_pi = solve(
    prob_const_pi,
    RadauIIA9(autodiff=false),
    abstol=1e-15,
    reltol=1e-8,
    saveat=[k for k in tspan[1]:((tspan[2]-tspan[1])/10_000):tspan[2]]
)

##
# Plot the results
size_inches = (6.5, 5)
size_pt = 72 .* size_inches
set_theme!(Theme(fontsize=6,
    Axis=(
        xticksize=1,
        yticksize=1,
        # xticklabelsize = 6,
        # yticklabelsize = 6,
        yticklabelpad=1,
        ylabelpadding=3
    )))
fig = Figure(size=size_pt)

full_color = Makie.wong_colors()[1]
no_reg_color = Makie.wong_colors()[6]
full_color = :Black
no_reg_color = :Red
condition_no_reg_color = Makie.wong_colors()[3]
full_linestyle = :dot
no_reg_linestyle = :dot
condition_no_reg_linestyle = nothing
adenine_pool_color = :Grey
adenine_pool_linestyle = :dash

line_ATPase_width = 2
line_ATPase_color = :grey
# no_allo_line_models_style = :dot
# no_reg_linestyle = [0.5, 1.5, 2.5, 3.5]
no_allo_line_models_style = [0.5, 1, 1.5, 2] .* 2

# pi_trap = load("110623_Glycolysis_schematic_pi_trap.png")
pi_trap = load("031025_Glycolysis_schematic_Harden_Young.png")



ax_pi_trap, im = image(fig[1, 1:2], rotr90(pi_trap), axis=(aspect=DataAspect(),))
# ax_pi_trap.alignmode = Mixed(top=-20, bottom=-20, left=-25, right=-15)
# ax_pi_trap.alignmode = Mixed(left=-20)
ax_pi_trap.width = 250
# ax_pi_trap.tellwidth = false
hidedecorations!(ax_pi_trap)
hidespines!(ax_pi_trap)

# Plot changing in Pi, intermediate before GAPDH and after GAPDH

# Plot metabolite with switch off/on of allostery without Pi uptake
ax_allost_on_off = Axis(
    fig[1, 3:4],
    limits=((0, maximum(sol.t)), nothing),
    xlabel="Time, min",
    # ylabel = "[Metabolite],M",
    ylabel="[Metabolite], normalized to Time=0",
    title="[Metabolite] changes in response to removal\nof HK and PFK allostery without Pi uptake",
    yscale=log10    # yticks = ([1e-8, 1e-6, 1e-4, 1e-2], ["0.01µM", "1µM", "0.1mM", "10mM"]),    # yticks = ([0.001, 0.01, 0.1, 1, 10, 100, 1000]),
)
allost_on_off_lines = []
ATP_color = :Red
ATP_style = nothing
ADP_color = :Red
ADP_style = :dot
Phosphate_color = :Red
Phosphate_style = :dash
Upper_color = Makie.wong_colors()[1]
Upper_style = nothing
Lower_color = Makie.wong_colors()[3]
Lower_style = nothing
ATP_color = :Red
ATP_style = nothing
ADP_color = Makie.wong_colors()[1]
ADP_style = nothing
Phosphate_color = Makie.wong_colors()[3]
Phosphate_style = nothing
Upper_color = (:Grey, 0.5)
Upper_style = nothing
Lower_color = (:Grey, 0.5)
Lower_style = :dot
for metabolite in [
    :G6P, :F6P, :F16BP, :GAP, :DHAP, :BPG, :TwoPG, :ThreePG, :PEP, :ATP, :ADP, :Phosphate]
    if metabolite ∈ [:G6P, :F6P, :F16BP, :GAP, :DHAP]
        color = Upper_color
        linestyle = Upper_style
    elseif metabolite ∈ [:BPG, :TwoPG, :ThreePG, :PEP]
        color = Lower_color
        linestyle = Lower_style
    elseif metabolite ∈ [:ATP]
        color = ATP_color
        linestyle = ATP_style
    elseif metabolite ∈ [:ADP]
        color = ADP_color
        linestyle = ADP_style
    elseif metabolite ∈ [:Phosphate]
        color = Phosphate_color
        linestyle = Phosphate_style
    end
    allost_on_off_line = lines!(
        ax_allost_on_off,
        sol.t,
        # sol[metabolite, :] ./ sol[metabolite, 1],
        [sol.u[i][metabolite] for i in 1:length(sol.u)] ./ sol.u[1][metabolite],
        label=String(metabolite),
        color=color,
        linestyle=linestyle
    )
    push!(allost_on_off_lines, allost_on_off_line)
end
axislegend(
    ax_allost_on_off,
    [
        allost_on_off_lines[end-2],
        allost_on_off_lines[end-1],
        allost_on_off_lines[end],
        allost_on_off_lines[4],
        allost_on_off_lines[end-3]],
    ["ATP", "ADP", "Phosphate", "Metabolites\nbefore GAPDH", "Metabolites\nafter GAPDH"],
    position=:lb,
    # position = (1.075, 1.05),
    rowgap=1,
    patchlabelgap=2,
    framevisible=false,
    padding=(-2, -5, 0, -4),
    patchsize=(7.5, 5)
)

# Plot metabolite with switch off/on of allostery with Pi uptake
ax_allost_on_off_const_pi = Axis(
    fig[1, 5:6],
    limits=((0, maximum(sol_const_pi.t)), nothing),
    xlabel="Time, min",
    ylabel="[Metabolite],M",
    title="[Metabolite] changes in response to removal\nof HK and PFK allostery with Pi uptake",
    # yscale=log10    # yticks = ([1e-8, 1e-6, 1e-4, 1e-2], ["0.01µM", "1µM", "0.1mM", "10mM"]),    # yticks = ([0.001, 0.01, 0.1, 1, 10, 100, 1000]),
)
allost_on_off_const_pi_lines = []
for metabolite in [
    :GAP, :DHAP, :BPG, :TwoPG, :ThreePG, :PEP, :ATP, :ADP, :Phosphate, :F6P, :G6P, :F16BP]
    if metabolite ∈ [:F16BP]
        color = Makie.wong_colors()[1]
        linestyle = nothing
    elseif metabolite ∈ [:G6P]
        color = Makie.wong_colors()[2]
        linestyle = nothing
    else
        color = (:Grey, 0.5)
        linestyle = nothing
    end
    allost_on_off_const_pi_line = lines!(
        ax_allost_on_off_const_pi,
        sol_const_pi.t,
        # sol[metabolite, :] ./ sol[metabolite, 1],
        [sol_const_pi.u[i][metabolite] for i in 1:length(sol_const_pi.u)], # ./ sol_const_pi.u[1][metabolite],
        label=String(metabolite),
        color=color,
        linestyle=linestyle
    )
    push!(allost_on_off_const_pi_lines, allost_on_off_const_pi_line)
end
lines!(
    ax_allost_on_off_const_pi,
    [sol_const_pi.t[1], sol_const_pi.t[end]],
    repeat([0.2], 2),
    color=:Grey,
    linestyle=:dash
)
text!(
    ax_allost_on_off_const_pi,
    0.01 * sol_const_pi.t[end],
    1.1 * 0.2,
    text="Estimated level of\nintracellular molecules",
    align=(:left, :bottom),
    color=adenine_pool_color
)
axislegend(
    ax_allost_on_off_const_pi,
    [
        allost_on_off_const_pi_lines[end],
        allost_on_off_const_pi_lines[end-1],
        allost_on_off_const_pi_lines[1]],
    ["F16BP", "G6P", "All other metabolites"],
    position=:lt,
    # position = (1.075, 1.05),
    rowgap=1,
    patchlabelgap=2,
    framevisible=false,
    padding=(-2, -5, 0, -4),
    patchsize=(7.5, 5)
)

adenine_pool_size = glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP
# Plot [ATP] maintenance with constant Phosphate
ax_ATPase_range = Axis(
    fig[2, 1:2],
    limits=((0.001, 1.0), (-0.2e-3, 1.75 * adenine_pool_size)),
    xlabel="ATPase, % of pathway Vmax",
    ylabel="[ATP],mM",
    title="Metabolite concentration changes\nwith constant [Phosphate]",
    xscale=log10,
    # yscale = log10,
    xticks=([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.9],
        ["0.1%", "0.3%", "1%", "3%", "10%", "30%", "90%"]),
    ytickformat=ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys]
)
lines!(
    ax_ATPase_range,
    Complete_Model_Simulation_Data.ATPase_Vmax,
    Complete_Model_Simulation_Data.ATP_conc,
    color=full_color,
    linestyle=full_linestyle,
    label="Full model"
)
lines!(
    ax_ATPase_range,
    No_Reg_Model_Simulation_Data.ATPase_Vmax,
    No_Reg_Model_Simulation_Data.ATP_conc,
    color=no_reg_color,
    linestyle=no_reg_linestyle,
    label="Model w/o allost."
)
lines!(
    ax_ATPase_range,
    Const_Pi_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    Const_Pi_No_Reg_Model_Simulation_Data.ATP_conc,
    color=condition_no_reg_color,
    linestyle=condition_no_reg_linestyle,
    label="Model w/o allost.\n[Phosphate]=1mM"
)
lines!(
    ax_ATPase_range,
    [0.001, 1.0],
    repeat(
        [glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP],
        2),
    color=adenine_pool_color,
    linestyle=adenine_pool_linestyle
)
text!(
    ax_ATPase_range,
    0.0033,
    1.02 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text="Adenine pool size",
    align=(:left, :bottom),
    color=adenine_pool_color
)
axislegend(
    ax_ATPase_range,
    position=:lt,
    # position = (1.075, 1.05),
    rowgap=1,
    framevisible=false,
    padding=(-2, -5, 0, -4),
    patchsize=(10, 5)
)

# Plot [ATP] maintenance with different Keq HK PFK GAPDH PGK
ax_ATPase_range = Axis(
    fig[2, 3:4],
    limits=((0.001, 1.0), (-0.2e-3, 1.75 * adenine_pool_size)),
    # limits = ((0.0001, 1.0), (-0.2e-3, 14e-3)),
    xlabel="ATPase, % of pathway Vmax",
    ylabel="[ATP],mM",
    title="Maintaining ATP concentration with lower\nKeq for HK1, PFKP, GAPDH, and PGK reactions",
    xscale=log10,
    # yscale = log10,
    xticks=([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.9],
        ["0.1%", "0.3%", "1%", "3%", "10%", "30%", "90%"]),
    ytickformat=ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys]
)
lines!(
    ax_ATPase_range,
    Complete_Model_Simulation_Data.ATPase_Vmax,
    Complete_Model_Simulation_Data.ATP_conc,
    color=full_color,
    linestyle=full_linestyle,
    label="Full model"
)
lines!(
    ax_ATPase_range,
    No_Reg_Model_Simulation_Data.ATPase_Vmax,
    No_Reg_Model_Simulation_Data.ATP_conc,
    color=no_reg_color,
    linestyle=no_reg_linestyle,
    label="Model w/o allost."
)
lines!(
    ax_ATPase_range,
    Keq_HK_PFK_GAPDH_PGK_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    Keq_HK_PFK_GAPDH_PGK_No_Reg_Model_Simulation_Data.ATP_conc,
    color=condition_no_reg_color,
    linestyle=condition_no_reg_linestyle,
    label="Model w/o allost.\nKeq HK1=Keq PFKP=0.05\nKeq GAPDH=0.5, Keq PGK=5"
)
lines!(
    ax_ATPase_range,
    [0.001, 1.0],
    repeat(
        [glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP],
        2),
    color=adenine_pool_color,
    linestyle=adenine_pool_linestyle
)
text!(
    ax_ATPase_range,
    0.0033,
    1.02 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text="Adenine pool size",
    align=(:left, :bottom),
    color=adenine_pool_color
)
axislegend(
    ax_ATPase_range,
    position=:lt,
    # position = (1.075, 1.05),
    rowgap=1,
    framevisible=false,
    padding=(-2, -5, 0, -4),
    patchsize=(10, 5)
)
ax_ATPase_range.width = 120

# Plot [ATP] maintenance with HK and PFK rate sub'd to ATPase rate
ax_ATPase_range = Axis(
    fig[2, 5:6],
    limits=((0.001, 1.0), (-0.2e-3, 1.75 * adenine_pool_size)),
    xlabel="ATPase, % of pathway Vmax",
    ylabel="[ATP],mM",
    title="Maintaining ATP concentration\nwith V(HK1)=V(PFKP)=0.5V(ATPase)",
    xscale=log10,
    # yscale = log10,
    xticks=([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.9],
        ["0.1%", "0.3%", "1%", "3%", "10%", "30%", "90%"]),
    ytickformat=ys -> ["$(Int(round(y*1000, sigdigits=2)))" for y in ys]
)
lines!(
    ax_ATPase_range,
    Complete_Model_Simulation_Data.ATPase_Vmax,
    Complete_Model_Simulation_Data.ATP_conc,
    color=full_color,
    linestyle=full_linestyle,
    label="Full model"
)
lines!(
    ax_ATPase_range,
    No_Reg_Model_Simulation_Data.ATPase_Vmax,
    No_Reg_Model_Simulation_Data.ATP_conc,
    color=no_reg_color,
    linestyle=no_reg_linestyle,
    label="Model w/o allost."
)
lines!(
    ax_ATPase_range,
    Vmax_HK_PFK_eq_ATPase_No_Reg_Model_Simulation_Data.ATPase_Vmax,
    Vmax_HK_PFK_eq_ATPase_No_Reg_Model_Simulation_Data.ATP_conc,
    color=condition_no_reg_color,
    linestyle=condition_no_reg_linestyle,
    label="Model w/o allost.\nV(HK1)=V(PFKP)=0.5V(ATPase)"
)
lines!(
    ax_ATPase_range,
    [0.001, 1.0],
    repeat(
        [glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP],
        2),
    color=adenine_pool_color,
    linestyle=adenine_pool_linestyle
)
text!(
    ax_ATPase_range,
    0.0033,
    1.02 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text="Adenine pool size",
    align=(:left, :bottom),
    color=adenine_pool_color
)
axislegend(
    ax_ATPase_range,
    # position = :lt,
    position=(0.5, 1.0),
    rowgap=1,
    framevisible=false,
    padding=(-2, -5, 0, -4),
    patchsize=(10, 5)
)

# linkaxes!(ax_allost_on_off, ax_allost_on_off_const_pi)
colgap!(fig.layout, 10)
rowgap!(fig.layout, 10)

label_a = fig[1, 1, TopLeft()] = Label(
    fig, "A", fontsize=12, halign=:right, padding=(0, 5, 5, 0))
label_b = fig[1, 3, TopLeft()] = Label(
    fig, "B", fontsize=12, halign=:right, padding=(0, 15, 5, 0))
label_c = fig[1, 5, TopLeft()] = Label(
    fig, "C", fontsize=12, halign=:right, padding=(0, 10, 5, 0))
label_d = fig[2, 1, TopLeft()] = Label(
    fig, "D", fontsize=12, halign=:right, padding=(0, 10, 5, 0))
label_e = fig[2, 3, TopLeft()] = Label(
    fig, "E", fontsize=12, halign=:right, padding=(0, 15, 5, 0))
label_f = fig[2, 5, TopLeft()] = Label(
    fig, "F", fontsize=12, halign=:right, padding=(0, 5, 5, 0))

fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_Fig5_mechanism_of_ATP_maintenence_by_allost.png", fig, px_per_unit = 4)
