using Glycolysis, CairoMakie, Dates

default_color = :black
activator_color = :green
max_ATP = 10e-3
max_F6P = 1e-3
max_G6P = 1e-3

set_theme!(
    Theme(
        fontsize = 20,
        Axis = (
            xlabelpadding = -20,
        ),
    ),
)
fig = Figure(; resolution = (900, 250))

ax_GPI = Axis(
    fig[1, 1],
    xlabel = "G6P, mM",
    ylabel = rich("Rate", superscript("GPI")),
    xtickformat = xs -> ["$(Int(round(x*1000, sigdigits=1)))" for x in xs],
    xticks = [0, max_G6P],
    xgridvisible = false,
    ygridvisible = false,
    yticksvisible =false,
    yticklabelsvisible =false,
)
lines!(
    0.0 .. max_G6P,
    x ->
        Glycolysis.rate_GPI(x, 0.0, glycolysis_params) /
        (glycolysis_params.GPI_Vmax * glycolysis_params.GPI_Conc),
    color = default_color,
)

ax_PFKP_F6P = Axis(
    fig[1, 2],
    xlabel = "F6P, mM",
    ylabel = rich("Rate", superscript("PFK")),
    xtickformat = xs -> ["$(Int(round(x*1000, sigdigits=1)))" for x in xs],
    xticks = [0, max_F6P],
    xgridvisible = false,
    ygridvisible = false,
    yticksvisible =false,
    yticklabelsvisible =false,
)
lines!(
    0.0 .. max_F6P,
    x ->
        Glycolysis.rate_PFKP(x, 5e-3, 0.0, 0.0, 0.0, 0.0, 0.0, glycolysis_params) /
        (glycolysis_params.PFKP_Vmax * glycolysis_params.PFKP_Conc),
    color = default_color,
)
lines!(
    0.0 .. max_F6P,
    x ->
        Glycolysis.rate_PFKP(x, 5e-3, 0.0, 0.5e-3, 1e-3, 0.0, 0.0, glycolysis_params) /
        (glycolysis_params.PFKP_Vmax * glycolysis_params.PFKP_Conc),
    color = activator_color,
)

ax_PFKP_ATP = Axis(
    fig[1, 3],
    xlabel = "ATP, mM",
    xtickformat = xs -> ["$(Int(round(x*1000, sigdigits=1)))" for x in xs],
    xticks = [0, max_ATP],
    xgridvisible = false,
    ygridvisible = false,
    yticksvisible =false,
    yticklabelsvisible =false,
)
lines!(
    0.0 .. max_ATP,
    x ->
        Glycolysis.rate_PFKP(0.5e-3, x, 0.0, 0.0, 0.0, 0.0, 0.0, glycolysis_params) /
        (glycolysis_params.PFKP_Vmax * glycolysis_params.PFKP_Conc),
    color = default_color,
)
lines!(
    0.0 .. max_ATP,
    x ->
        Glycolysis.rate_PFKP(0.5e-3, x, 0.0, 0.5e-3, 1e-3, 0.0, 0.0, glycolysis_params) /
        (glycolysis_params.PFKP_Vmax * glycolysis_params.PFKP_Conc),
    color = activator_color,
)
colgap!(fig.layout, 10)
fig

# uncomment the line below to save the plot
# save("Results/$(Dates.format(now(),"mmddyy"))_Fig1_GPI_PFK_plots.png", fig, px_per_unit = 4)

