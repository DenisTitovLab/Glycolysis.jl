# include("GlycolysisModel.jl")
using Glycolysis
using DataFrames, CSV, CairoMakie, PrettyTables, CairoMakie, Dates

# GSA_ATP_AUC_100K = CSV.read("110422_ATP_AUC_gsa_sobol_100000_3x_range.csv", DataFrame)
# GSA_ATP_AUC_20K_no_HK = CSV.read("110522_ATP_AUC_gsa_sobol_20000_3x_range_fixed_HK_Km_Vmax.csv", DataFrame)
# GSA_ATP_energy_AUC_20K = CSV.read("110422_ATPase_energy_AUC_gsa_sobol_20000_3x_range.csv", DataFrame)
# GSA_ATP_energy_AUC_20K_no_HK = CSV.read("110622_ATPase_energy_AUC_gsa_sobol_20000_3x_range_fixed_HK_Km_Vmax.csv", DataFrame)
# GSA_ATP_prod_AUC_20K = CSV.read("110422_ATPprod_AUC_gsa_sobol_20000_3x_range.csv", DataFrame)
# GSA_ATP_prod_AUC_20K_no_HK = CSV.read("110622_ATPprod_AUC_gsa_sobol_20000_3x_range_fixed_HK_Km_Vmax.csv", DataFrame)


# GSA_ATP_AUC_20K = CSV.read("111522_ATP_AUC_gsa_sobol_20000_3x_range_HK_PFK_Vmax.csv", DataFrame)
GSA_ATP_AUC_20K = CSV.read("050723_ATP_AUC_gsa_sobol_20000_3x_range_HK_PFK_Vmax.csv", DataFrame)
# GSA_ATP_energy_AUC_20K =
#     CSV.read("111522_ATPase_energy_AUC_gsa_sobol_20000_3x_range_HK_PFK_Vmax.csv", DataFrame)
GSA_ATP_energy_AUC_20K =
    CSV.read("050723_ATPase_energy_AUC_gsa_sobol_20000_3x_range_HK_PFK_Vmax.csv", DataFrame)
# GSA_ATP_prod_AUC_20K = CSV.read("111522_ATPprod_AUC_gsa_sobol_20000_3x_range_HK_PFK_Vmax.csv", DataFrame)
GSA_ATP_prod_AUC_20K = CSV.read("050723_ATPprod_AUC_gsa_sobol_20000_3x_range_HK_PFK_Vmax.csv", DataFrame)



vscodedisplay(GSA_ATP_AUC_20K)
vscodedisplay(GSA_ATP_energy_AUC_20K)
vscodedisplay(GSA_ATP_prod_AUC_20K)




scatterlines(collect(values(GSA_ATP_AUC_20K[1, 2:end])))
scatterlines(collect(values(GSA_ATP_AUC_20K[1, 2:end])))
scatterlines(collect(values(GSA_ATP_energy_AUC_20K[1, 2:end])))
scatterlines(collect(values(GSA_ATP_energy_AUC_20K[1, 2:end])))
scatterlines(collect(values(GSA_ATP_prod_AUC_20K[1, 2:end])))
scatterlines(collect(values(GSA_ATP_prod_AUC_20K[1, 2:end])))


##
HK1_Km_Vmax =
    [:HK1_K_Glucose, :HK1_K_a_ATP, :HK1_β_Glucose_ATP, :HK1_K_G6P, :HK1_K_a_ADP, :HK1_Conc, :HK1_Vmax]
HK1_reg = [:HK1_K_i_G6P_reg, :HK1_K_a_G6P_cat, :HK1_K_a_Pi]
PFKP_Km_Vmax = [:PFKP_K_a_F6P, :PFKP_K_ATP, :PFKP_K_F16BP, :PFKP_K_ADP, :PFKP_Conc, :PFKP_Vmax]
PFKP_reg = [:PFKP_K_i_ATP_reg, :PFKP_K_a_ADP_reg, :PFKP_K_Phosphate]
HK1_PFKP_reg = [HK1_reg; PFKP_reg]
PKM2_reg = [:PKM2_a_KdF16BP, :PKM2_i_KdF16BP]
GAPDH = [param for param in propertynames(glycolysis_params) if occursin("GAPDH", string(param))]
PGK = [param for param in propertynames(glycolysis_params) if occursin("PGK", string(param))]
Keqs = [param for param in propertynames(glycolysis_params) if occursin("Keq", string(param))]
Other = [
    param for param in propertynames(glycolysis_params) if param ∉ [HK1_Km_Vmax; HK1_reg; PFKP_Km_Vmax; PFKP_reg; GAPDH; PGK]
]

# params_group_names = [:All, :HK1_Km_Vmax, :HK1_reg, :PFKP_Km_Vmax, :PFKP_reg, :GAPDH_all, :Other]
# params_names =
#     [[name for name in propertynames(glycolysis_params)], HK1_Km_Vmax, HK1_reg, PFKP_Km_Vmax, PFKP_reg, GAPDH_all, Other]

params_group_names = [:All, :HK1_Km_Vmax, :HK1_PFKP_reg, :PFKP_Km_Vmax, :GAPDH, :PGK, :All_Other]
params_names =
    [[name for name in propertynames(glycolysis_params)], HK1_Km_Vmax, HK1_PFKP_reg, PFKP_Km_Vmax, GAPDH, PGK, Other]

Sensitivity_indexes = DataFrame(
    Parameters = params_group_names,
    S1 = [sum(Iterators.filter(!isnan, GSA_ATP_AUC_20K[1, params_name])) for params_name in params_names],
    ST = [sum(Iterators.filter(!isnan, GSA_ATP_AUC_20K[2, params_name])) for params_name in params_names],
    # S1_no_HK = [sum(GSA_ATP_AUC_20K_no_HK[1, params_name]) for params_name in params_names],
    # ST_no_HK = [sum(GSA_ATP_AUC_20K_no_HK[2, params_name]) for params_name in params_names],
)
Sensitivity_indexes = DataFrame(
    Parameters = replace.(String.(params_group_names), "_" => " "),
    S1 = [sum(Iterators.filter(!isnan, GSA_ATP_AUC_20K[1, params_name])) for params_name in params_names],
    ST = [sum(Iterators.filter(!isnan, GSA_ATP_AUC_20K[2, params_name])) for params_name in params_names],
    # S1_no_HK = [sum(GSA_ATP_AUC_20K_no_HK[1, params_name]) for params_name in params_names],
    # ST_no_HK = [sum(GSA_ATP_AUC_20K_no_HK[2, params_name]) for params_name in params_names],
)

formatter = (v, i, j) -> (j > 1) ? round(v, sigdigits = 2) : v
pretty_table(Sensitivity_indexes, crop = :none, formatters = formatter)
CSV.write("$(Dates.format(now(),"mmddyy"))_ATP_AUC_gsa_sobol_sens_ind_HK_PFK_Vmax.csv", Sensitivity_indexes)


##

Sensitivity_indexes = DataFrame(
    Parameters = params_group_names,
    S1 = [
        sum(Iterators.filter(!isnan, GSA_ATP_energy_AUC_20K[1, params_name])) for params_name in params_names
    ],
    ST = [
        sum(Iterators.filter(!isnan, GSA_ATP_energy_AUC_20K[2, params_name])) for params_name in params_names
    ],
    # S1_no_HK = [sum(GSA_ATP_energy_AUC_20K_no_HK[1, params_name]) for params_name in params_names],
    # ST_no_HK = [sum(GSA_ATP_energy_AUC_20K_no_HK[2, params_name]) for params_name in params_names],
)
Sensitivity_indexes = DataFrame(
    Parameters = replace.(String.(params_group_names), "_" => " "),
    S1 = [
        sum(Iterators.filter(!isnan, GSA_ATP_energy_AUC_20K[1, params_name])) for params_name in params_names
    ],
    ST = [
        sum(Iterators.filter(!isnan, GSA_ATP_energy_AUC_20K[2, params_name])) for params_name in params_names
    ],
    # S1_no_HK = [sum(GSA_ATP_energy_AUC_20K_no_HK[1, params_name]) for params_name in params_names],
    # ST_no_HK = [sum(GSA_ATP_energy_AUC_20K_no_HK[2, params_name]) for params_name in params_names],
)

formatter = (v, i, j) -> (j > 1) ? round(v, sigdigits = 2) : v
pretty_table(Sensitivity_indexes, crop = :none, formatters = formatter)
Sensitivity_indexes
CSV.write(
    "$(Dates.format(now(),"mmddyy"))_ATP_energy_AUC_gsa_sobol_sens_ind_HK_PFK_Vmax.csv",
    Sensitivity_indexes,
)

##

Sensitivity_indexes = DataFrame(
    Parameters = params_group_names,
    S1 = [
        sum(Iterators.filter(!isnan, GSA_ATP_prod_AUC_20K[1, params_name])) for params_name in params_names
    ],
    ST = [
        sum(Iterators.filter(!isnan, GSA_ATP_prod_AUC_20K[2, params_name])) for params_name in params_names
    ],
    # S1_no_HK = [sum(GSA_ATP_prod_AUC_20K_no_HK[1, params_name]) for params_name in params_names],
    # ST_no_HK = [sum(GSA_ATP_prod_AUC_20K_no_HK[2, params_name]) for params_name in params_names],
)
Sensitivity_indexes = DataFrame(
    Parameters = replace.(String.(params_group_names), "_" => " "),
    S1 = [
        sum(Iterators.filter(!isnan, GSA_ATP_prod_AUC_20K[1, params_name])) for params_name in params_names
    ],
    ST = [
        sum(Iterators.filter(!isnan, GSA_ATP_prod_AUC_20K[2, params_name])) for params_name in params_names
    ],
    # S1_no_HK = [sum(GSA_ATP_prod_AUC_20K_no_HK[1, params_name]) for params_name in params_names],
    # ST_no_HK = [sum(GSA_ATP_prod_AUC_20K_no_HK[2, params_name]) for params_name in params_names],
)

formatter = (v, i, j) -> (j > 1) ? round(v, sigdigits = 2) : v
pretty_table(Sensitivity_indexes, crop = :none, formatters = formatter)

CSV.write(
    "$(Dates.format(now(),"mmddyy"))_ATP_prod_AUC_gsa_sobol_sens_ind_HK_PFK_Vmax.csv",
    Sensitivity_indexes,
)


##
##



df = DataFrame(a = 1:5, b = 'a':'e', c = rand(5))
io = IOBuffer()
pretty_table(io, df, tf = tf_markdown)
str = String(io.data)

fig = Figure();
ax = Axis(fig[1, 1])
xlims!(ax, -1, 5);
ylims!(ax, 0, 2);
hidespines!(ax)
hidedecorations!(ax)
text!(ax, str; family = "monospace", pointsize = 8)
fig


