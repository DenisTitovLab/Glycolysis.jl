using Glycolysis
using DifferentialEquations
using CairoMakie, Dates, Printf

##
include("for_FigS2_Teusink2014Model.jl")
# Precalculate output of Teusink model
initial_concentrations = LVector(
    Glucose_media = 110e-3,
    Glucose = 0.087e-3,
    G6P = 3.085e-3,
    F6P = 0.752e-3,
    F16BP = 0.836e-3,
    GAP = 0.518e-3 * 0.045,
    DHAP = 0.518e-3 * (1 - 0.045),
    BPG = 0.111e-3,
    ThreePG = 0.825e-3,
    TwoPG = 0.138e-3,
    PEP = 0.140e-3,
    Pyruvate = 0.884e-3,
    Acetald = 0.047e-3,
    EtOH = 0.0e-3,
    ATP = 2.06e-3,
    ADP = 0.87e-3,
    AMP = 0.165e-3,
    NTP = 1.6e-3,
    NDP = 0.0,
    Phosphate = 10.0e-3,
    NAD = 1.5e-3,
    NADH = 0.044e-3,
    F26BP = 0.0,
)

tspan = (0.0, 240.0 / 4)
Initial_ATPase_Vmax_frac = 0.03
High_ATPase_Vmax_frac = Initial_ATPase_Vmax_frac * 2
Low_ATPase_Vmax_frac = Initial_ATPase_Vmax_frac / 2
ATPase_change_time1 = 60 / 4
ATPase_change_time2 = 120 / 4
ATPase_change_time3 = 180 / 4
#params.ATPase_Vmax = Initial_ATPase_Vmax_frac * params.MCT_Conc * params.MCT_Vmax
params.ATPase_Vmax = 2 * Initial_ATPase_Vmax_frac * params.PFKP_Vmax
params.ATPase_Km_ATP = 1e-9

function affect1!(integrator)
    integrator.p.ATPase_Vmax = 2 * High_ATPase_Vmax_frac * params.PFKP_Vmax
end
function affect2!(integrator)
    integrator.p.ATPase_Vmax = 2 * Low_ATPase_Vmax_frac * params.PFKP_Vmax
end
function affect3!(integrator)
    integrator.p.ATPase_Vmax = 2 * Initial_ATPase_Vmax_frac * params.PFKP_Vmax
end

PresetTime_cb1 = PresetTimeCallback(ATPase_change_time1, affect1!)
PresetTime_cb2 = PresetTimeCallback(ATPase_change_time2, affect2!)
PresetTime_cb3 = PresetTimeCallback(ATPase_change_time3, affect3!)
cb_set = CallbackSet(PresetTime_cb1, PresetTime_cb2, PresetTime_cb3)

init_cond_prob = ODEProblem(Glycolysis_ODEs, initial_concentrations, (0, 1e8), params)
init_cond_sol = solve(init_cond_prob, Rodas4P2(), abstol = 1e-15, reltol = 1e-9, save_everystep = false)
new_init_cond = init_cond_sol.u[end]
prob = ODEProblem(Glycolysis_ODEs, new_init_cond, tspan, params, callback = cb_set)
sol = solve(
    prob,
    Rodas4P2(),
    abstol = 1e-15,
    reltol = 1e-9,
    saveat = [k for k = tspan[1]:((tspan[2]-tspan[1])/10_000):tspan[2]],
)

teusink_model_timepoints = sol.t
teusink_model_ATPprod = [Conc_to_Rates(conc, params).ATPprod for conc in sol.u] / (2 * params.PFKP_Vmax)
teusink_model_ATPenergy = [Conc_to_MassAction(conc, params).Q_Keq_ATPase for conc in sol.u]
teusink_model_ATPase_rate = Initial_ATPase_Vmax_frac * ones(length(sol))
teusink_model_ATPase_rate[sol.t.>ATPase_change_time1.&&sol.t.<ATPase_change_time2] .= High_ATPase_Vmax_frac
teusink_model_ATPase_rate[sol.t.>ATPase_change_time2.&&sol.t.<ATPase_change_time3] .= Low_ATPase_Vmax_frac
teusink_model_ATP = [conc.ATP for conc in sol.u]

glycolysis_params.ATPase_Vmax =
    Initial_ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax
glycolysis_params.ATPase_Km_ATP = 1e-9

function affect1!(integrator)
    integrator.p.ATPase_Vmax =
        High_ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax
end
function affect2!(integrator)
    integrator.p.ATPase_Vmax =
        Low_ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax
end
function affect3!(integrator)
    integrator.p.ATPase_Vmax =
        Initial_ATPase_Vmax_frac * 2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax
end

PresetTime_cb1 = PresetTimeCallback(ATPase_change_time1, affect1!)
PresetTime_cb2 = PresetTimeCallback(ATPase_change_time2, affect2!)
PresetTime_cb3 = PresetTimeCallback(ATPase_change_time3, affect3!)
cb_set = CallbackSet(PresetTime_cb1, PresetTime_cb2, PresetTime_cb3)

init_cond_prob = ODEProblem(glycolysis_ODEs, glycolysis_init_conc, (0, 1e8), glycolysis_params)
init_cond_sol = solve(init_cond_prob, Rodas5P(), abstol = 1e-15, reltol = 1e-5, save_everystep = false)
new_init_cond = init_cond_sol.u[end]
prob = ODEProblem(glycolysis_ODEs, new_init_cond, tspan, glycolysis_params, callback = cb_set)
sol = solve(
    prob,
    Rodas5P(),
    abstol = 1e-15,
    reltol = 1e-5,
    saveat = [k for k = tspan[1]:((tspan[2]-tspan[1])/10_000):tspan[2]],
)
our_model_timepoints = sol.t
our_model_ATPprod =
    [Glycolysis.conc_to_rates(conc, glycolysis_params).ATPprod for conc in sol.u] / (2 * glycolysis_params.HK1_Conc * glycolysis_params.HK1_Vmax)
our_model_ATPenergy =
    [Glycolysis.conc_to_disequilibrium_ratios(conc, glycolysis_params).Q_Keq_ATPase for conc in sol.u]
our_model_ATPase_rate = Initial_ATPase_Vmax_frac * ones(length(sol))
our_model_ATPase_rate[sol.t.>ATPase_change_time1.&&sol.t.<ATPase_change_time2] .= High_ATPase_Vmax_frac
our_model_ATPase_rate[sol.t.>ATPase_change_time2.&&sol.t.<ATPase_change_time3] .= Low_ATPase_Vmax_frac
our_model_ATP = [conc.ATP for conc in sol.u]



##

# Plot the results
size_inches = (6.5, 3)
size_pt = 72 .* size_inches
set_theme!(
    Theme(
        fontsize = 6,
        Axis = (
            xticksize = 1,
            yticksize = 1,
            # xticklabelsize = 5,
            # yticklabelsize = 5,
            yticklabelpad = 1,
            ylabelpadding = 3,
        ),
    ),
)
fig = Figure(resolution = size_pt)

# full_color = :black
# teusink_color = :red

full_color = Makie.wong_colors()[1]
teusink_color = Makie.wong_colors()[6]
line_ATPase_width = 2
line_ATPase_color = :grey
teusink_linestyle = :dash
line_models_width = 4

# Plot ATP production
ax_ATP_prod = Axis(
    fig[1, 1],
    limits = (nothing, (0.0, 1.5 * maximum(our_model_ATPase_rate))),
    # limits = (nothing, (0.99, 1.01)),
    xlabel = "Time, min",
    ylabel = "ATP production rate, relative to glycolysis Vmax",
    title = "Matching ATP supply and demand",
)
ax_ATPase = Axis(
    fig[1, 1],
    limits = (nothing, (0.0, 1.5 * maximum(our_model_ATPase_rate))),
    ylabel = "ATPase rate, relative to glycolysis Vmax",
    yaxisposition = :right,
    ygridvisible = false,
)

hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)
our_model_ATP_prod_line = lines!(ax_ATP_prod, our_model_timepoints, our_model_ATPprod, color = full_color)
teusink_model_ATP_prod_line = lines!(
    ax_ATP_prod,
    teusink_model_timepoints,
    teusink_model_ATPprod,
    linestyle = teusink_linestyle,
    color = teusink_color,
)
ATPase_line = lines!(
    ax_ATPase,
    our_model_timepoints,
    our_model_ATPase_rate,
    linestyle = :dot,
    color = line_ATPase_color,
    linewidth = line_ATPase_width,
)

axislegend(
    ax_ATP_prod,
    [our_model_ATP_prod_line, teusink_model_ATP_prod_line, ATPase_line],
    ["Our Model", "van Heerden Model", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
    patchsize = (7.5, 7.5),
)

# Plot [ATP]
ax_ATP_conc = Axis(
    fig[1, 2],
    limits = (nothing, (2e-6, 15e-3)),
    xlabel = "Time, min",
    ylabel = "[ATP], mM",
    title = "Maintaining ATP concentration",
    yscale = log10,
    ytickformat = ys -> ["$(round(y*1000, sigdigits = 3))" for y in ys],
)
ax_ATPase = Axis(
    fig[1, 2],
    limits = (nothing, (0.0, 1.5 * maximum(our_model_ATPase_rate))),
    ylabel = "ATPase rate, relative to glycolysis Vmax",
    yaxisposition = :right,
    ygridvisible = false,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)
ATP_line = lines!(ax_ATP_conc, our_model_timepoints, our_model_ATP, color = full_color)
teusink_ATP_line = lines!(
    ax_ATP_conc,
    teusink_model_timepoints,
    teusink_model_ATP,
    linestyle = teusink_linestyle,
    color = teusink_color,
)
ATPase_line = lines!(
    ax_ATPase,
    our_model_timepoints,
    our_model_ATPase_rate,
    linestyle = :dot,
    color = line_ATPase_color,
    linewidth = line_ATPase_width,
)
lines!(
    ax_ATP_conc,
    [tspan[1], tspan[2]],
    repeat([glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP], 2),
    color = :grey,
    linestyle = :dash,
)
text!(
    ax_ATP_conc,
    0,
    1.05 * (glycolysis_init_conc.ATP + glycolysis_init_conc.ADP + glycolysis_init_conc.AMP),
    text = "Adenine pool size",
    align = (:left, :bottom),
    color = :grey,
)
axislegend(
    ax_ATP_conc,
    [ATP_line, teusink_ATP_line, ATPase_line],
    ["Our Model", "van Heerden Model", "ATPase rate"],
    # position = :rt,
    position = (1, 0.9),
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
    patchsize = (7.5, 7.5),
)

# Plot energy of ATPase reaction
ax_ATP_energy = Axis(
    fig[1, 3],
    limits = (nothing, (0, 30)),
    xlabel = "Time, min",
    ylabel = "Energy of ATPase reaction, kT",
    title = "Maintaining ATP energy",
    ytickformat = ys -> ["$(Int(round(y)))" for y in ys],
)
ax_ATPase = Axis(
    fig[1, 3],
    limits = (nothing, (0.0, 1.5 * maximum(our_model_ATPase_rate))),
    ylabel = "ATPase rate, relative to glycolysis Vmax",
    yaxisposition = :right,
    ygridvisible = false,
)
hidespines!(ax_ATPase)
hidexdecorations!(ax_ATPase)
ATPase_energy_line =
    lines!(ax_ATP_energy, our_model_timepoints, -log.(our_model_ATPenergy), color = full_color)
teusink_ATPase_energy_line = lines!(
    ax_ATP_energy,
    teusink_model_timepoints,
    -log.(teusink_model_ATPenergy),
    linestyle = teusink_linestyle,
    color = teusink_color,
)
ATPase_line = lines!(
    ax_ATPase,
    our_model_timepoints,
    our_model_ATPase_rate,
    linestyle = :dot,
    color = line_ATPase_color,
    linewidth = line_ATPase_width,
)
axislegend(
    ax_ATP_energy,
    [ATPase_energy_line, teusink_ATPase_energy_line, ATPase_line],
    ["Our Model", "van Heerden Model", "ATPase rate"],
    position = :rt,
    rowgap = 1,
    framevisible = false,
    padding = (-5, -5, 0, -8),
    patchsize = (7.5, 7.5),
)

colgap!(fig.layout, 7.5)

label_a = fig[1, 1, TopLeft()] = Label(fig, "A", fontsize = 12, halign = :right, padding = (0, 5, 5, 0))
label_b = fig[1, 2, TopLeft()] = Label(fig, "B", fontsize = 12, halign = :right, padding = (0, 5, 5, 0))
label_c = fig[1, 3, TopLeft()] = Label(fig, "C", fontsize = 12, halign = :right, padding = (0, 5, 5, 0))

fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_FigS2_Teusink_model_comparison.png", fig, px_per_unit = 4)


