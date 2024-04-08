#=
    This file contains functions that take M concentrations of glycolytic metabolites as input and
    and generate rates of glycolytic reactions in M/s
=#

"""
    rate_GLUT(metabs, params)

Calculate rate (M/s units) of GLUT transporter from concentrations (M units) of `Glucose_media` and `Glucose` according to the following equation:

```math
Rate = \frac{{V_{max} \cdot Conc}}{{K_{M}^{Glucose}}} \cdot \frac{{Glucose_{media} - \frac{1}{K_{eq}} \cdot Glucose}}{1 + \frac{Glucose_{media}}{K_{M}^{Glucose}} + \frac{Glucose}{K_{M}^{Glucose}}}
```

# Arguments
- `metabs`: LArray or NamedTuple and struct that contains fields Glucose_media, Glucose with corresponding metabolite concentrations of GLUT substrates and products. Glycolysis.jl exports LArray `glycolysis_init_cond` that contains estimates of cellular metabolite concentrations.
- `params`: LArray or NamedTuple and struct of kinetic parameters of GLUT. Glycolysis.jl exports LArray `glycolysis_params` that contains kinetic parameters of all glycolytic enzymes.

# Example
```julia-repl
julia> metabs = (Glucose_media = 25e-3, Glucose = 8e-3,)
julia> Glycolysis.rate_GLUT(metabs, glycolysis_params)
0.02267962645321136
```
"""
function rate_GLUT(metabs, params)
    Rate = (
        (params.GLUT_Vmax * params.GLUT_Conc / params.GLUT_Km_Glucose) *
        (metabs.Glucose_media - (1 / params.GLUT_Keq) * metabs.Glucose) /
        (1 + metabs.Glucose_media / params.GLUT_Km_Glucose + metabs.Glucose / params.GLUT_Km_Glucose)
    )
    return Rate
end

function rate_HK1(metabs, params)
    Z = (
        (
            1 +
            (metabs.Glucose / params.HK1_K_Glucose) +
            (metabs.ATP / params.HK1_K_a_ATP) +
            (metabs.G6P / params.HK1_K_G6P) +
            (metabs.G6P / params.HK1_K_a_G6P_cat) +
            (metabs.ADP / params.HK1_K_a_ADP) +
            (params.HK1_β_Glucose_ATP) * (metabs.Glucose / params.HK1_K_Glucose) * (metabs.ATP / params.HK1_K_a_ATP) +
            (metabs.Glucose / params.HK1_K_Glucose) * (metabs.ADP / params.HK1_K_a_ADP) +
            (metabs.Glucose / params.HK1_K_Glucose) * (metabs.G6P / params.HK1_K_a_G6P_cat) +
            (metabs.G6P / params.HK1_K_G6P) * (metabs.ADP / params.HK1_K_a_ADP) +
            (metabs.G6P / params.HK1_K_G6P) * (metabs.G6P / params.HK1_K_a_G6P_cat)
        ) * (1 + (metabs.Phosphate / params.HK1_K_a_Pi)) +
        (1 + (metabs.Glucose / params.HK1_K_Glucose) + (metabs.G6P / params.HK1_K_G6P)) * (metabs.G6P / params.HK1_K_i_G6P_reg)
    )
    Rate =
        (
            (
                params.HK1_Vmax *
                params.HK1_Conc *
                (params.HK1_β_Glucose_ATP) *
                (1 / params.HK1_K_Glucose) *
                (1 / params.HK1_K_a_ATP) *
                (1 + (metabs.Phosphate / params.HK1_K_a_Pi))
            ) * (metabs.Glucose * metabs.ATP - metabs.G6P * metabs.ADP / params.HK1_Keq)
        ) / Z
    return Rate
end

function rate_GPI(metabs, params)
    Rate = (
        (params.GPI_Vmax * params.GPI_Conc / params.GPI_Km_G6P) * (metabs.G6P - (1 / params.GPI_Keq) * metabs.F6P) /
        (1 + metabs.G6P / params.GPI_Km_G6P + metabs.F6P / params.GPI_Km_F6P)
    )
    return Rate
end

function rate_PFKP(metabs, params)

    Z_a_cat = (
        1 +
        (metabs.F6P / params.PFKP_K_a_F6P) +
        (metabs.ATP / params.PFKP_K_ATP) +
        (metabs.F16BP / params.PFKP_K_F16BP) +
        (metabs.ADP / params.PFKP_K_ADP) +
        (metabs.F6P / params.PFKP_K_a_F6P) * (metabs.ATP / params.PFKP_K_ATP) +
        (metabs.F16BP / params.PFKP_K_F16BP) * (metabs.ADP / params.PFKP_K_ADP)
    )
    Z_i_cat = (
        1 +
        (metabs.ATP / params.PFKP_K_ATP) +
        (metabs.F16BP / params.PFKP_K_F16BP) +
        (metabs.ADP / params.PFKP_K_ADP) +
        (metabs.F16BP / params.PFKP_K_F16BP) * (metabs.ADP / params.PFKP_K_ADP)
    )
    Z_a_reg = (
        (1 + metabs.Phosphate / params.PFKP_K_Phosphate) *
        (1 + metabs.ADP / params.PFKP_K_a_ADP_reg) *
        (1 + metabs.F26BP / params.PFKP_K_a_F26BP)
    )
    Z_i_reg = (
        (1 + metabs.ATP / params.PFKP_K_i_ATP_reg + metabs.Phosphate / params.PFKP_K_Phosphate) *
        (1 + metabs.F26BP / params.PFKP_K_i_F26BP) *
        (1 + metabs.Citrate / params.PFKP_K_i_Citrate)
    )

    Rate =
        params.PFKP_Vmax *
        params.PFKP_Conc *
        (metabs.F6P * metabs.ATP - metabs.F16BP * metabs.ADP / params.PFKP_Keq) *
        (1 / params.PFKP_K_a_F6P) *
        (1 / params.PFKP_K_ATP) *
        (Z_a_cat^3) *
        (Z_a_reg^4) / ((Z_a_cat^4) * (Z_a_reg^4) + params.PFKP_L * (Z_i_cat^4) * (Z_i_reg^4))

    return Rate
end

function rate_ALDO(metabs, params)
    Rate = (
        (params.ALDO_Vmax * params.ALDO_Conc / params.ALDO_Km_F16BP) * (
            (metabs.F16BP - (1 / params.ALDO_Keq) * (metabs.DHAP * metabs.GAP)) / (
                1 +
                metabs.GAP * metabs.DHAP / (params.ALDO_Kd_DHAP * params.ALDO_Km_GAP) +
                metabs.DHAP / params.ALDO_Kd_DHAP +
                metabs.F16BP * metabs.GAP / (params.ALDO_Ki_GAP * params.ALDO_Km_F16BP) +
                metabs.F16BP / params.ALDO_Km_F16BP +
                metabs.GAP * params.ALDO_Km_DHAP / (params.ALDO_Kd_DHAP * params.ALDO_Km_GAP)
            )
        )
    )
    return Rate
end

function rate_TPI(metabs, params)
    Rate = (
        (params.TPI_Vmax * params.TPI_Conc / params.TPI_Km_DHAP) * (metabs.DHAP - (1 / params.TPI_Keq) * metabs.GAP) /
        (1 + (metabs.DHAP / params.TPI_Km_DHAP) + (metabs.GAP / params.TPI_Km_GAP))
    )
    return Rate
end

function rate_GAPDH(metabs, params)
    Z_a =
        (
            1 +
            metabs.GAP / params.GAPDH_K_GAP * (1 + metabs.Phosphate / params.GAPDH_K_a_Phosphate) +
            metabs.BPG / params.GAPDH_K_BPG
        ) * (1 + metabs.NAD / params.GAPDH_K_a_NAD + metabs.NADH / params.GAPDH_K_a_NADH)

    Z_i =
        (1 + metabs.NAD / params.GAPDH_K_i_NAD) * (
            1 +
            metabs.GAP / params.GAPDH_K_GAP * (1 + metabs.Phosphate / params.GAPDH_K_i_Phosphate) +
            metabs.BPG / params.GAPDH_K_BPG
        ) +
        metabs.NADH / params.GAPDH_K_i_NADH * (
            1 +
            metabs.GAP / params.GAPDH_K_GAP * (1 + metabs.Phosphate / params.GAPDH_K_i_Phosphate) +
            metabs.BPG / (params.GAPDH_α_i_BPG * params.GAPDH_K_BPG)
        )

    Rate = (
        (
            params.GAPDH_Vmax * params.GAPDH_Conc /
            (params.GAPDH_K_GAP * params.GAPDH_K_a_NAD * params.GAPDH_K_a_Phosphate)
        ) *
        Z_a^3 *
        (metabs.GAP * metabs.NAD * metabs.Phosphate - (1 / params.GAPDH_Keq) * metabs.BPG * metabs.NADH) / (Z_a^4 + params.GAPDH_L * Z_i^4)
    )
    return Rate
end

function rate_PGK(metabs, params)
    Rate = (
        (params.PGK_Vmax * params.PGK_Conc / (params.PGK_α * params.PGK_K_BPG * params.PGK_K_ADP)) *
        (metabs.BPG * metabs.ADP - (1 / params.PGK_Keq) * (metabs.ThreePG * metabs.ATP)) / (
            1 +
            metabs.BPG / params.PGK_K_BPG +
            metabs.ADP / params.PGK_K_ADP +
            metabs.ThreePG / params.PGK_K_ThreePG +
            metabs.ATP / params.PGK_K_ATP +
            metabs.BPG * metabs.ADP / (params.PGK_α * params.PGK_K_BPG * params.PGK_K_ADP) +
            metabs.ThreePG * metabs.ATP / (params.PGK_β * params.PGK_K_ThreePG * params.PGK_K_ATP) +
            metabs.ThreePG * metabs.ADP / (params.PGK_γ * params.PGK_K_ThreePG * params.PGK_K_ADP)
        )
    )
    return Rate
end

function rate_PGM(metabs, params)
    Rate = (
        (params.PGM_Vmax * params.PGM_Conc / params.PGM_Km_ThreePG) *
        (metabs.ThreePG - (1 / params.PGM_Keq) * metabs.TwoPG) /
        (1 + metabs.ThreePG / params.PGM_Km_ThreePG + metabs.TwoPG / params.PGM_Km_TwoPG)
    )
    return Rate
end

function rate_ENO(metabs, params)
    Rate = (
        (params.ENO_Vmax * params.ENO_Conc / params.ENO_Km_TwoPG) * (metabs.TwoPG - (1 / params.ENO_Keq) * metabs.PEP) /
        (1 + metabs.TwoPG / params.ENO_Km_TwoPG + metabs.PEP / params.ENO_Km_PEP)
    )
    return Rate
end

function rate_PKM2(metabs, params)

    Z_a_cat = (
        1 +
        (metabs.PEP / params.PKM2_K_a_PEP) +
        (metabs.ATP / params.PKM2_K_ATP) +
        (metabs.ADP / params.PKM2_K_ADP) +
        (metabs.PEP / params.PKM2_K_a_PEP) * (metabs.ADP / params.PKM2_K_ADP) +
        (metabs.Pyruvate / params.PKM2_K_Pyruvate) * (metabs.ATP / params.PKM2_K_ATP) +
        (metabs.PEP / params.PKM2_K_a_PEP) * (metabs.ATP / params.PKM2_K_ATP) +
        (metabs.Pyruvate / params.PKM2_K_Pyruvate) * (metabs.ADP / params.PKM2_K_ADP)
    )
    Z_i_cat = (
        1 +
        (metabs.PEP / params.PKM2_K_i_PEP) +
        (metabs.ATP / params.PKM2_K_ATP) +
        (metabs.ADP / params.PKM2_K_ADP) +
        (metabs.PEP / params.PKM2_K_i_PEP) * (metabs.ADP / params.PKM2_K_ADP) +
        (metabs.Pyruvate / params.PKM2_K_Pyruvate) * (metabs.ATP / params.PKM2_K_ATP) +
        params.PKM2_β_i_PEP_ATP * (metabs.PEP / params.PKM2_K_i_PEP) * (metabs.ATP / params.PKM2_K_ATP) +
        (metabs.Pyruvate / params.PKM2_K_Pyruvate) * (metabs.ADP / params.PKM2_K_ADP)
    )
    Z_a_reg = ((1 + metabs.F16BP / params.PKM2_K_a_F16BP) * (1 + metabs.Phenylalanine / params.PKM2_K_a_Phenylalanine))
    Z_i_reg = (1 + metabs.Phenylalanine / params.PKM2_K_i_Phenylalanine)

    Rate =
        params.PKM2_Conc *
        (metabs.ADP * metabs.PEP - metabs.ATP * metabs.Pyruvate / params.PKM2_Keq) *
        (
            (params.PKM2_Vmax_a * (1.0 / params.PKM2_K_a_PEP) * (1.0 / params.PKM2_K_ADP)) *
            (Z_a_cat^3) *
            (Z_a_reg^4) +
            params.PKM2_L *
            (params.PKM2_Vmax_i * (1.0 / params.PKM2_K_i_PEP) * (1.0 / params.PKM2_K_ADP)) *
            (Z_i_cat^3) *
            (Z_i_reg^4)
        ) / ((Z_a_cat^4) * (Z_a_reg^4) + params.PKM2_L * (Z_i_cat^4) * (Z_i_reg^4))

    return Rate
end

function rate_LDH(metabs, params)
    Rate = (
        (params.LDH_Vmax * params.LDH_Conc / (params.LDH_Km_Pyruvate * params.LDH_Kd_NADH)) *
        (metabs.Pyruvate * metabs.NADH - (1 / params.LDH_Keq) * (metabs.Lactate * metabs.NAD)) / (
            1 +
            metabs.Pyruvate * params.LDH_Km_NADH / (params.LDH_Kd_NADH * params.LDH_Km_Pyruvate) +
            metabs.Lactate * params.LDH_Km_NAD / (params.LDH_Kd_NAD * params.LDH_Km_Lactate) +
            metabs.NADH / params.LDH_Kd_NADH +
            metabs.Lactate * metabs.NAD / (params.LDH_Kd_NAD * params.LDH_Km_Lactate) +
            metabs.Lactate * metabs.NADH * params.LDH_Km_NAD /
            (params.LDH_Kd_NAD * params.LDH_Kd_NADH * params.LDH_Km_Lactate) +
            metabs.Pyruvate * metabs.NADH / (params.LDH_Kd_NADH * params.LDH_Km_Pyruvate) +
            metabs.NAD / params.LDH_Kd_NAD +
            metabs.Pyruvate * metabs.NAD * params.LDH_Km_NADH /
            (params.LDH_Kd_NAD * params.LDH_Kd_NADH * params.LDH_Km_Pyruvate)
        )
    )
    return Rate
end

function rate_MCT(metabs, params)
    Rate = (
        (params.MCT_Vmax * params.MCT_Conc / params.MCT_Km_Lactate) *
        (metabs.Lactate - (1 / params.MCT_Keq) * metabs.Lactate_media) /
        (1 + metabs.Lactate / params.MCT_Km_Lactate + metabs.Lactate_media / params.MCT_Km_Lactate)
    )
    return Rate
end

function rate_AK(metabs, params)
    Rate = (
        (params.AK_Vmax / (params.AK_Km_ADP^2)) * (metabs.ADP^2 - (1 / params.AK_Keq) * (metabs.ATP * metabs.AMP)) / (
            (1 + metabs.ADP / params.AK_Km_ADP + metabs.ATP / params.AK_Km_ATP) *
            (1 + metabs.ADP / params.AK_Km_ADP + metabs.AMP / params.AK_Km_AMP)
        )
    )
    return Rate
end

function rate_NDPK(metabs, params)
    Rate = (
        (params.NDPK_Vmax / (params.NDPK_Km_ATP * params.NDPK_Km_NDP)) * (
            (metabs.ATP * metabs.NDP - (1 / params.NDPK_Keq) * (metabs.NTP * metabs.ADP)) / (
                (1 + metabs.ATP / params.NDPK_Km_ATP + metabs.ADP / params.NDPK_Km_ADP) *
                (1 + metabs.NTP / params.NDPK_Km_NTP + metabs.NDP / params.NDPK_Km_NDP)
            )
        )
    )
    return Rate
end

function rate_CK(metabs, params)
    Rate = (
        (params.CK_Vmax / (params.CK_Km_ATP * params.CK_Km_Creatine)) * (
            (metabs.ATP * metabs.Creatine - (1 / params.CK_Keq) * (metabs.Phosphocreatine * metabs.ADP)) / (
                (1 + metabs.ATP / params.CK_Km_ATP + metabs.ADP / params.CK_Km_ADP) *
                (1 + metabs.Phosphocreatine / params.CK_Km_Phosphocreatine + metabs.Creatine / params.CK_Km_Creatine)
            )
        )
    )
    return Rate
end

function rate_ATPase(metabs, params)
    Rate =
        (params.ATPase_Vmax / params.ATPase_Km_ATP) *
        (
            1 / (
                1 +
                metabs.ATP / params.ATPase_Km_ATP +
                metabs.ADP / params.ATPase_Km_ADP +
                metabs.Phosphate / params.ATPase_Km_Phosphate +
                (metabs.ADP / params.ATPase_Km_ADP) * (metabs.Phosphate / params.ATPase_Km_Phosphate)
            )
        ) *
        (metabs.ATP - metabs.ADP * metabs.Phosphate / params.ATPase_Keq)
    return Rate
end
