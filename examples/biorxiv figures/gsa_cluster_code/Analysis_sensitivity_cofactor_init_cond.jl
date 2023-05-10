# Define a MWC function for PFKFB kinase domain
function PFK2(F6P, ATP, Pi, PEP, ADP, F26BP, P::Array{Float64})
    # Assign the values in parameter array (P) to kinetic variables
    dE,
    Vmax_a,
    Vmax_i,
    VmaxR_a,
    KmF6P_a,
    KmF6P_i,
    KmATP_a,
    KmATP_i,
    KdP_a,
    KdP_i,
    KdPEP_a,
    KdPEP_i,
    KdADP_a,
    KdADP_i,
    KdF26BP_a,
    KdF26BP_i,
    VmaxR_a,
    VmaxR_i = P

    Vmax_a = 1.0

    # Calculate Keq for PFKFB1 Kinase using eQuibilibrator (37C, pH = 7.4, pMg = 1, I = 0.25 M)
    dG = 1.8 # +/- 5 kJ/mol, from eQuibilibrator
    R = 8.3144598 * 10^-3 # constant, kJ/mol*K
    T = 310.15 # 37C in Kelvin

    Keq = exp(-dG / (R*T))

    # Active PFKFB1 Kinase in Active Conformation
    VmaxR_a = (Vmax_a / Keq) * ((1 / (KmF6P_a * KmATP_a)) / (1 / (KdF26BP_a * KdADP_a)))

    cat_a =
        (1 + (F6P / KmF6P_a) + (PEP / KdPEP_a) + (F26BP / KdF26BP_a)) *
        (1 + (ATP / KmATP_a) + (ADP / KdADP_a))

    reg_a = 1 + (Pi / KdP_a)

    a = Vmax_a * (1 / (KmF6P_a * KmATP_a)) * (cat_a) * (reg_a)^2

    # Active PFKFB1 Kinase in Inactive Conformation

    VmaxR_i = (Vmax_i / Keq) * ((1 / (KmF6P_i * KmATP_i)) / (1 / (KdF26BP_i * KdADP_i)))

    cat_i =
        (1 + (F6P / KmF6P_i) + (PEP / KdPEP_i) + (F26BP / KdF26BP_i)) *
        (1 + (ATP / KmATP_i) + (ADP / KdADP_i))

    reg_i = 1 + (Pi / KdP_i)

    i = Vmax_i * (1 / (KmF6P_i * KmATP_i)) * (cat_i) * (reg_i)^2

    # Total PFKFB1 States for both Active & Inactive Confromation

    Za = (
        (
            1 +
            (ATP / KmATP_a) +
            (ADP / KdADP_a) +
            (F26BP / KdF26BP_a) +
            (F6P / KmF6P_a) +
            ((F6P * ATP) / (KmF6P_a * KmATP_a)) +
            ((F6P * ADP) / (KmF6P_a * KdADP_a)) +
            (PEP / KdPEP_a) +
            ((PEP * ATP) / (KdPEP_a * KmATP_a)) +
            ((PEP * ADP) / (KdPEP_a * KdADP_a)) +
            ((F26BP * ADP) / (KdF26BP_a * KdADP_a))
        )^2 * (1 + (Pi / KdP_a))^2
    )

    Zi = (
        (
            1 +
            (ATP / KmATP_i) +
            (ADP / KdADP_i) +
            (F26BP / KdF26BP_i) +
            (F6P / KmF6P_i) +
            ((F6P * ATP) / (KmF6P_i * KmATP_i)) +
            ((F6P * ADP) / (KmF6P_i * KdADP_i)) +
            (PEP / KdPEP_i) +
            ((PEP * ATP) / (KdPEP_i * KmATP_i)) +
            ((PEP * ADP) / (KdPEP_i * KdADP_i)) +
            ((F26BP * ADP) / (KdF26BP_i * KdADP_i))
        )^2 * (1 + (Pi / KdP_i))^2
    )

    # Return MWC

    rate = ((a + exp(-dE) * i) / (Za + exp(-dE) * Zi)) * ((F6P * ATP) - ((F26BP * ADP) / Keq))

    return rate
end



##
using Glycolysis, DataFrames, CSV, CairoMakie, PrettyTables, CairoMakie, Dates


GSA_ATP_init_cond_AUC_20K = CSV.read("112122_ATP_AUC_gsa_sobol_init_cond_20000_3x_range.csv", DataFrame)
GSA_ATP_AUC_cofactor_pools_20K =
    CSV.read("050723_ATP_AUC_gsa_sobol_cofactor_pool_20000_3x_range.csv", DataFrame)

vscodedisplay(GSA_ATP_init_cond_AUC_20K)
vscodedisplay(GSA_ATP_AUC_cofactor_pools_20K)
sum(GSA_ATP_init_cond_AUC_20K[1, 2:end])
sum(GSA_ATP_init_cond_AUC_20K[2, 2:end])
sum(GSA_ATP_AUC_cofactor_pools_20K[1, 2:end])
sum(GSA_ATP_AUC_cofactor_pools_20K[2, 2:end])


scatterlines(collect(values(GSA_ATP_AUC_100K[1, 2:end])))
scatterlines(collect(values(GSA_ATP_AUC_20K_no_HK[1, 2:end])))
scatterlines(collect(values(GSA_ATP_energy_AUC_20K[1, 2:end])))
scatterlines(collect(values(GSA_ATP_energy_AUC_20K_no_HK[1, 2:end])))
scatterlines(collect(values(GSA_ATP_prod_AUC_20K[1, 2:end])))
scatterlines(collect(values(GSA_ATP_prod_AUC_20K_no_HK[1, 2:end])))


##
HK1_Km_Vmax =
    [:HK1_K_Glucose, :HK1_K_a_ATP, :HK1_β_Glucose_ATP, :HK1_K_G6P, :HK1_K_a_ADP, :HK1_Conc, :HK1_Vmax]
HK1_reg = [:HK1_K_i_G6P_reg, :HK1_K_a_G6P_cat, :HK1_K_a_Pi]
PFKP_Km_Vmax = [:PFKP_K_a_F6P, :PFKP_K_ATP, :PFKP_K_F16BP, :PFKP_K_ADP, :PFKP_Conc, :PFKP_Vmax]
PFKP_reg = [:PFKP_K_i_ATP_reg, :PFKP_K_a_ADP_reg, :PFKP_K_Phosphate]
HK1_PFKP_reg = [HK1_reg; PFKP_reg]
PKM2_reg = [:PKM2_a_KdF16BP, :PKM2_i_KdF16BP]
GAPDH = [param for param in propertynames(params) if occursin("GAPDH", string(param))]
Keqs = [param for param in propertynames(params) if occursin("Keq", string(param))]
Other = [
    param for param in propertynames(params) if param ∉ [HK1_Km_Vmax; HK1_reg; PFKP_Km_Vmax; PFKP_reg; GAPDH]
]

# params_group_names = [:All, :HK1_Km_Vmax, :HK1_reg, :PFKP_Km_Vmax, :PFKP_reg, :GAPDH_all, :Other]
# params_names =
#     [[name for name in propertynames(params)], HK1_Km_Vmax, HK1_reg, PFKP_Km_Vmax, PFKP_reg, GAPDH_all, Other]

params_group_names = [:All, :HK1_Km_Vmax, :HK1_PFKP_reg, :PFKP_Km_Vmax, :GAPDH, :Other]
params_names =
    [[name for name in propertynames(params)], HK1_Km_Vmax, HK1_PFKP_reg, PFKP_Km_Vmax, GAPDH, Other]

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


##
using Distributions
Normal.(params, model_params_uncertainty)
res =
    truncated.(
        Normal.(convert(Vector{Float64}, params), convert(Vector{Float64}, model_params_uncertainty));
        lower = 0.0,
    )
params
rand(res[136])
counter = 0
for stuff in res
    counter += 1
    println(counter)
    rand(stuff)
end
res[1:2]
rand.(truncated.(Normal.(params, model_params_uncertainty); lower = 0.0))
[rand.(truncated.(Normal.(params, model_params_uncertainty); lower = 0.0)) for i = 1:number_boostraps]