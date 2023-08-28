#=
    This file contains functions that take M concentrations of glycolytic metabolites as input and
    and generate rates of glycolytic reactions in M/s
=#

using LabelledArrays

function conc_to_rates(s, params)
    r = @LVector eltype(s) (
        :GLUT,
        :HK1,
        :GPI,
        :PFKP,
        :ALDO,
        :TPI,
        :GAPDH,
        :PGK,
        :PGM,
        :ENO,
        :PKM2,
        :LDH,
        :MCT,
        :ATPprod,
        :ATPase,
    )
    r.GLUT = rate_GLUT(s.Glucose_media, s.Glucose, params)
    r.HK1 = rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params)
    r.GPI = rate_GPI(s.G6P, s.F6P, params)
    r.PFKP = rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params)
    r.ALDO = rate_ALDO(s.F16BP, s.GAP, s.DHAP, params)
    r.TPI = rate_TPI(s.GAP, s.DHAP, params)
    r.GAPDH = rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    r.PGK = rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params)
    r.PGM = rate_PGM(s.ThreePG, s.TwoPG, params)
    r.ENO = rate_ENO(s.TwoPG, s.PEP, params)
    r.PKM2 = rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params)
    r.LDH = rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params)
    r.MCT = rate_MCT(s.Lactate, s.Lactate_media, params)
    r.ATPprod = (
        -rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) -
        rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) +
        rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) +
        rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) +
        rate_AK(s.ATP, s.ADP, s.AMP, params)  -
        rate_NDPK(s.NTP, s.NDP, s.ATP, s.ADP, params) -
        rate_CK(s.Phosphocreatine, s.Creatine, s.ATP, s.ADP, params)
    )
    r.ATPase = rate_ATPase(s.ATP, s.ADP, s.Phosphate, params)
    return r
end

"""
    rate_GLUT(Glucose_media, Glucose, glycolysis_params)

Calculate rate (M/s units) of GLUT transporter from concentrations (M units) of `Glucose_media` and `Glucose`

# Arguments
- `params::LArray`: kinetic parameters of GLUT. Glycolysis.jl exports LArray `glycolysis_params` that contains kinetic parameters of all glycolytic enzymes. 

# Example
```julia-repl
julia> Glucose_media, Glucose = 25e-3, 8e-3
(0.025, 0.008)
julia> Glycolysis.rate_GLUT(Glucose_media, Glucose, glycolysis_params)
0.02267962645321136
```
"""
function rate_GLUT(Glucose_media, Glucose, params)
    Rate = (
        (params.GLUT_Vmax * params.GLUT_Conc / params.GLUT_Km_Glucose) *
        (Glucose_media - (1 / params.GLUT_Keq) * Glucose) /
        (1 + Glucose_media / params.GLUT_Km_Glucose + Glucose / params.GLUT_Km_Glucose)
    )
    return Rate
end

function rate_HK1(Glucose, G6P, ATP, ADP, Phosphate, params)
    Z = (
        (
            1 +
            (Glucose / params.HK1_K_Glucose) +
            (ATP / params.HK1_K_a_ATP) +
            (G6P / params.HK1_K_G6P) +
            (G6P / params.HK1_K_a_G6P_cat) +
            (ADP / params.HK1_K_a_ADP) +
            (params.HK1_β_Glucose_ATP) * (Glucose / params.HK1_K_Glucose) * (ATP / params.HK1_K_a_ATP) +
            (Glucose / params.HK1_K_Glucose) * (ADP / params.HK1_K_a_ADP) +
            (Glucose / params.HK1_K_Glucose) * (G6P / params.HK1_K_a_G6P_cat) +
            (G6P / params.HK1_K_G6P) * (ADP / params.HK1_K_a_ADP) +
            (G6P / params.HK1_K_G6P) * (G6P / params.HK1_K_a_G6P_cat)
        ) * (1 + (Phosphate / params.HK1_K_a_Pi)) +
        (1 + (Glucose / params.HK1_K_Glucose) + (G6P / params.HK1_K_G6P)) * (G6P / params.HK1_K_i_G6P_reg)
    )
    Rate =
        (
            (
                params.HK1_Vmax *
                params.HK1_Conc *
                (params.HK1_β_Glucose_ATP) *
                (1 / params.HK1_K_Glucose) *
                (1 / params.HK1_K_a_ATP) *
                (1 + (Phosphate / params.HK1_K_a_Pi))
            ) * (Glucose * ATP - G6P * ADP / params.HK1_Keq)
        ) / Z
    return Rate
end

function rate_GPI(G6P, F6P, params)
    Rate = (
        (params.GPI_Vmax * params.GPI_Conc / params.GPI_Km_G6P) * (G6P - (1 / params.GPI_Keq) * F6P) /
        (1 + G6P / params.GPI_Km_G6P + F6P / params.GPI_Km_F6P)
    )
    return Rate
end

function rate_PFKP(F6P, ATP, F16BP, ADP, Phosphate, Citrate, F26BP, params)

    Z_a_cat = (
        1 +
        (F6P / params.PFKP_K_a_F6P) +
        (ATP / params.PFKP_K_ATP) +
        (F16BP / params.PFKP_K_F16BP) +
        (ADP / params.PFKP_K_ADP) +
        (F6P / params.PFKP_K_a_F6P) * (ATP / params.PFKP_K_ATP) +
        (F16BP / params.PFKP_K_F16BP) * (ADP / params.PFKP_K_ADP)
    )
    Z_i_cat = (
        1 +
        (ATP / params.PFKP_K_ATP) +
        (F16BP / params.PFKP_K_F16BP) +
        (ADP / params.PFKP_K_ADP) +
        (F16BP / params.PFKP_K_F16BP) * (ADP / params.PFKP_K_ADP)
    )
    Z_a_reg = (
        (1 + Phosphate / params.PFKP_K_Phosphate) *
        (1 + ADP / params.PFKP_K_a_ADP_reg) *
        (1 + F26BP / params.PFKP_K_a_F26BP)
    )
    Z_i_reg = (
        (1 + ATP / params.PFKP_K_i_ATP_reg + Phosphate / params.PFKP_K_Phosphate) *
        (1 + F26BP / params.PFKP_K_i_F26BP) *
        (1 + Citrate / params.PFKP_K_i_Citrate)
    )

    Rate =
        params.PFKP_Vmax *
        params.PFKP_Conc *
        (F6P * ATP - F16BP * ADP / params.PFKP_Keq) *
        (1 / params.PFKP_K_a_F6P) *
        (1 / params.PFKP_K_ATP) *
        (Z_a_cat^3) *
        (Z_a_reg^4) / ((Z_a_cat^4) * (Z_a_reg^4) + params.PFKP_L * (Z_i_cat^4) * (Z_i_reg^4))

    return Rate
end

function rate_ALDO(F16BP, GAP, DHAP, params)
    Rate = (
        (params.ALDO_Vmax * params.ALDO_Conc / params.ALDO_Km_F16BP) * (
            (F16BP - (1 / params.ALDO_Keq) * (DHAP * GAP)) / (
                1 +
                GAP * DHAP / (params.ALDO_Kd_DHAP * params.ALDO_Km_GAP) +
                DHAP / params.ALDO_Kd_DHAP +
                F16BP * GAP / (params.ALDO_Ki_GAP * params.ALDO_Km_F16BP) +
                F16BP / params.ALDO_Km_F16BP +
                GAP * params.ALDO_Km_DHAP / (params.ALDO_Kd_DHAP * params.ALDO_Km_GAP)
            )
        )
    )
    return Rate
end

function rate_TPI(GAP, DHAP, params)
    Rate = (
        (params.TPI_Vmax * params.TPI_Conc / params.TPI_Km_DHAP) * (DHAP - (1 / params.TPI_Keq) * GAP) /
        (1 + (DHAP / params.TPI_Km_DHAP) + (GAP / params.TPI_Km_GAP))
    )
    return Rate
end

function rate_GAPDH(GAP, NAD, Phosphate, BPG, NADH, params)
    Z_a =
        (
            1 +
            GAP / params.GAPDH_K_GAP * (1 + Phosphate / params.GAPDH_K_a_Phosphate) +
            BPG / params.GAPDH_K_BPG
        ) * (1 + NAD / params.GAPDH_K_a_NAD + NADH / params.GAPDH_K_a_NADH)

    Z_i =
        (1 + NAD / params.GAPDH_K_i_NAD) * (
            1 +
            GAP / params.GAPDH_K_GAP * (1 + Phosphate / params.GAPDH_K_i_Phosphate) +
            BPG / params.GAPDH_K_BPG
        ) +
        NADH / params.GAPDH_K_i_NADH * (
            1 +
            GAP / params.GAPDH_K_GAP * (1 + Phosphate / params.GAPDH_K_i_Phosphate) +
            BPG / (params.GAPDH_α_i_BPG * params.GAPDH_K_BPG)
        )

    Rate = (
        (
            params.GAPDH_Vmax * params.GAPDH_Conc /
            (params.GAPDH_K_GAP * params.GAPDH_K_a_NAD * params.GAPDH_K_a_Phosphate)
        ) *
        Z_a^3 *
        (GAP * NAD * Phosphate - (1 / params.GAPDH_Keq) * BPG * NADH) / (Z_a^4 + params.GAPDH_L * Z_i^4)
    )
    return Rate
end

function rate_PGK(BPG, ADP, ATP, ThreePG, params)
    Rate = (
        (params.PGK_Vmax * params.PGK_Conc / (params.PGK_α * params.PGK_K_BPG * params.PGK_K_ADP)) *
        (BPG * ADP - (1 / params.PGK_Keq) * (ThreePG * ATP)) / (
            1 +
            BPG / params.PGK_K_BPG +
            ADP / params.PGK_K_ADP +
            ThreePG / params.PGK_K_ThreePG +
            ATP / params.PGK_K_ATP +
            BPG * ADP / (params.PGK_α * params.PGK_K_BPG * params.PGK_K_ADP) +
            ThreePG * ATP / (params.PGK_β * params.PGK_K_ThreePG * params.PGK_K_ATP) +
            ThreePG * ADP / (params.PGK_γ * params.PGK_K_ThreePG * params.PGK_K_ADP)
        )
    )
    return Rate
end

function rate_PGM(ThreePG, TwoPG, params)
    Rate = (
        (params.PGM_Vmax * params.PGM_Conc / params.PGM_Km_ThreePG) *
        (ThreePG - (1 / params.PGM_Keq) * TwoPG) /
        (1 + ThreePG / params.PGM_Km_ThreePG + TwoPG / params.PGM_Km_TwoPG)
    )
    return Rate
end

function rate_ENO(TwoPG, PEP, params)
    Rate = (
        (params.ENO_Vmax * params.ENO_Conc / params.ENO_Km_TwoPG) * (TwoPG - (1 / params.ENO_Keq) * PEP) /
        (1 + TwoPG / params.ENO_Km_TwoPG + PEP / params.ENO_Km_PEP)
    )
    return Rate
end

# function rate_PKM2(PEP, ADP, F16BP, ATP, Pyruvate, params)
#     Z = (
#         ((1 + PEP / params.PKM2_a_KmPEP)^4) *
#         ((1 + ADP / params.PKM2_a_KmADP + ATP / params.PKM2_a_KdATP)^4) *
#         ((1 + F16BP / params.PKM2_a_KdF16BP)^4) +
#         params.PKM2_L *
#         ((1 + PEP / params.PKM2_i_KmPEP)^4) *
#         ((1 + ADP / params.PKM2_i_KmADP + ATP / params.PKM2_i_KdATP)^4) *
#         ((1 + F16BP / params.PKM2_i_KdF16BP)^4)
#     )
#     Pa = (
#         4 *
#         (1 / params.PKM2_a_KmPEP) *
#         (1 / params.PKM2_a_KmADP) *
#         ((1 + PEP / params.PKM2_a_KmPEP)^(4 - 1)) *
#         ((1 + ADP / params.PKM2_a_KmADP + ATP / params.PKM2_a_KdATP)^(4 - 1)) *
#         ((1 + F16BP / params.PKM2_a_KdF16BP)^4) / Z
#     )
#     Pi = (
#         4 *
#         params.PKM2_L *
#         (1 / params.PKM2_i_KmPEP) *
#         (1 / params.PKM2_i_KmADP) *
#         ((1 + PEP / params.PKM2_i_KmPEP)^(4 - 1)) *
#         ((1 + ADP / params.PKM2_i_KmADP + ATP / params.PKM2_i_KdATP)^(4 - 1)) *
#         ((1 + F16BP / params.PKM2_i_KdF16BP)^4) / Z
#     )
#     Rate = (
#         (1 / 4) *
#         (params.PKM2_a_Vmax * params.PKM2_Conc * Pa + params.PKM2_i_Vmax * params.PKM2_Conc * Pi) *
#         (ADP * PEP - ATP * Pyruvate / params.PKM2_Keq)
#     )
#     return Rate
# end

function rate_PKM2(PEP, ADP, Pyruvate, ATP, F16BP, Phenylalanine, params)

    Z_a_cat = (
        1 +
        (PEP / params.PKM2_K_a_PEP) +
        (ATP / params.PKM2_K_ATP) +
        (ADP / params.PKM2_K_ADP) +
        (PEP / params.PKM2_K_a_PEP) * (ADP / params.PKM2_K_ADP) +
        (Pyruvate / params.PKM2_K_Pyruvate) * (ATP / params.PKM2_K_ATP) +
        (PEP / params.PKM2_K_a_PEP) * (ATP / params.PKM2_K_ATP) +
        (Pyruvate / params.PKM2_K_Pyruvate) * (ADP / params.PKM2_K_ADP)
    )
    Z_i_cat = (
        1 +
        (PEP / params.PKM2_K_i_PEP) +
        (ATP / params.PKM2_K_ATP) +
        (ADP / params.PKM2_K_ADP) +
        (PEP / params.PKM2_K_i_PEP) * (ADP / params.PKM2_K_ADP) +
        (Pyruvate / params.PKM2_K_Pyruvate) * (ATP / params.PKM2_K_ATP) +
        params.PKM2_β_i_PEP_ATP * (PEP / params.PKM2_K_i_PEP) * (ATP / params.PKM2_K_ATP) +
        (Pyruvate / params.PKM2_K_Pyruvate) * (ADP / params.PKM2_K_ADP)
    )
    Z_a_reg = ((1 + F16BP / params.PKM2_K_a_F16BP) * (1 + Phenylalanine / params.PKM2_K_a_Phenylalanine))
    Z_i_reg = (1 + Phenylalanine / params.PKM2_K_i_Phenylalanine)

    Rate =
        params.PKM2_Conc *
        (ADP * PEP - ATP * Pyruvate / params.PKM2_Keq) *
        (
            (
                params.PKM2_Vmax_a * (1.0 / params.PKM2_K_a_PEP) * (1.0 / params.PKM2_K_ADP)
            ) *
            (Z_a_cat^3) *
            (Z_a_reg^4) +
            params.PKM2_L *
            (
                params.PKM2_Vmax_i * (1.0 / params.PKM2_K_i_PEP) * (1.0 / params.PKM2_K_ADP)
            ) *
            (Z_i_cat^3) *
            (Z_i_reg^4)
        ) / ((Z_a_cat^4) * (Z_a_reg^4) + params.PKM2_L * (Z_i_cat^4) * (Z_i_reg^4))

    return Rate
end

function rate_LDH(Pyruvate, NADH, NAD, Lactate, params)
    Rate = (
        (params.LDH_Vmax * params.LDH_Conc / (params.LDH_Km_Pyruvate * params.LDH_Kd_NADH)) *
        (Pyruvate * NADH - (1 / params.LDH_Keq) * (Lactate * NAD)) / (
            1 +
            Pyruvate * params.LDH_Km_NADH / (params.LDH_Kd_NADH * params.LDH_Km_Pyruvate) +
            Lactate * params.LDH_Km_NAD / (params.LDH_Kd_NAD * params.LDH_Km_Lactate) +
            NADH / params.LDH_Kd_NADH +
            Lactate * NAD / (params.LDH_Kd_NAD * params.LDH_Km_Lactate) +
            Lactate * NADH * params.LDH_Km_NAD /
            (params.LDH_Kd_NAD * params.LDH_Kd_NADH * params.LDH_Km_Lactate) +
            Pyruvate * NADH / (params.LDH_Kd_NADH * params.LDH_Km_Pyruvate) +
            NAD / params.LDH_Kd_NAD +
            Pyruvate * NAD * params.LDH_Km_NADH /
            (params.LDH_Kd_NAD * params.LDH_Kd_NADH * params.LDH_Km_Pyruvate)
        )
    )
    return Rate
end

function rate_MCT(Lactate, Lactate_media, params)
    Rate = (
        (params.MCT_Vmax * params.MCT_Conc / params.MCT_Km_Lactate) *
        (Lactate - (1 / params.MCT_Keq) * Lactate_media) /
        (1 + Lactate / params.MCT_Km_Lactate + Lactate_media / params.MCT_Km_Lactate)
    )
    return Rate
end

function rate_AK(ATP, ADP, AMP, params)
    Rate = (
        (params.AK_Vmax / (params.AK_Km_ADP^2)) * (ADP^2 - (1 / params.AK_Keq) * (ATP * AMP)) / (
            (1 + ADP / params.AK_Km_ADP + ATP / params.AK_Km_ATP) *
            (1 + ADP / params.AK_Km_ADP + AMP / params.AK_Km_AMP)
        )
    )
    return Rate
end

function rate_NDPK(NTP, NDP, ATP, ADP, params)
    Rate = (
        (params.NDPK_Vmax / (params.NDPK_Km_ATP * params.NDPK_Km_NDP)) * (
            (ATP * NDP - (1 / params.NDPK_Keq) * (NTP * ADP)) / (
                1 +
                ATP / params.NDPK_Km_ATP +
                ADP / params.NDPK_Km_ADP +
                NTP / params.NDPK_Km_NTP +
                NDP / params.NDPK_Km_NDP
            )
        )
    )
    return Rate
end

function rate_CK(Phosphocreatine, Creatine, ATP, ADP, params)
    Rate = (
        (params.CK_Vmax / (params.CK_Km_ATP * params.CK_Km_Creatine)) * (
            (ATP * Creatine - (1 / params.CK_Keq) * (Phosphocreatine * ADP)) / (
                1 +
                ATP / params.CK_Km_ATP +
                ADP / params.CK_Km_ADP +
                Phosphocreatine / params.CK_Km_Phosphocreatine +
                Creatine / params.CK_Km_Creatine
            )
        )
    )
    return Rate
end

function rate_ATPase(ATP, ADP, Phosphate, params)
    Rate =
        (params.ATPase_Vmax / params.ATPase_Km_ATP) *
        (
            1 / (
                1 +
                ATP / params.ATPase_Km_ATP +
                ADP / params.ATPase_Km_ADP +
                Phosphate / params.ATPase_Km_Phosphate +
                (ADP / params.ATPase_Km_ADP) * (Phosphate / params.ATPase_Km_Phosphate)
            )
        ) *
        (ATP - ADP * Phosphate / params.ATPase_Keq)
    return Rate
end