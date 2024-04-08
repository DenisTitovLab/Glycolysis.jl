#=
    This file contains functions that take M concentrations of glycolytic metabolites as input and
    and generate concentration of metabolites bound to corresponding enzymes in M
=#

function binding_GLUT(metabs, params)
    Glucose_bound = (
        (params.GLUT_Conc / params.GLUT_MW) *
        (metabs.Glucose_media / params.GLUT_Km_Glucose + metabs.Glucose / params.GLUT_Km_Glucose) /
        (1 + metabs.Glucose_media / params.GLUT_Km_Glucose + metabs.Glucose / params.GLUT_Km_Glucose)
    )
    return (Glucose = Glucose_bound,)
end

function binding_HK1(metabs, params)
    Z_a_cat = (
        1 +
        (metabs.Glucose / params.HK1_K_Glucose) +
        (metabs.ATP / params.HK1_K_a_ATP) +
        (metabs.G6P / params.HK1_K_G6P) +
        (metabs.G6P / params.HK1_K_a_G6P_cat) +
        (metabs.ADP / params.HK1_K_a_ADP) +
        (params.HK1_β_Glucose_ATP) *
        (metabs.Glucose / params.HK1_K_Glucose) *
        (metabs.ATP / params.HK1_K_a_ATP) +
        (metabs.Glucose / params.HK1_K_Glucose) * (metabs.ADP / params.HK1_K_a_ADP) +
        (metabs.Glucose / params.HK1_K_Glucose) * (metabs.G6P / params.HK1_K_a_G6P_cat) +
        (metabs.G6P / params.HK1_K_G6P) * (metabs.ADP / params.HK1_K_a_ADP) +
        (metabs.G6P / params.HK1_K_G6P) * (metabs.G6P / params.HK1_K_a_G6P_cat)
    )
    Z_i_cat = (1 + (metabs.Glucose / params.HK1_K_Glucose) + (metabs.G6P / params.HK1_K_G6P))
    Z_a_reg = (1 + (metabs.Phosphate / params.HK1_K_a_Pi))
    Z_i_reg = (metabs.G6P / params.HK1_K_i_G6P_reg)
    Z = Z_a_cat * Z_a_reg + Z_i_cat * Z_i_reg

    Glucose_bound =
        (
            (params.HK1_Conc / (params.HK1_MW)) * (
                (metabs.Glucose / params.HK1_K_Glucose) * Z_a_reg +
                (metabs.Glucose / params.HK1_K_Glucose) * Z_i_reg
            )
        ) / Z
    G6P_bound =
        (
            (params.HK1_Conc / (params.HK1_MW)) * (
                (
                    (metabs.G6P / params.HK1_K_G6P) +
                    (metabs.G6P / params.HK1_K_a_G6P_cat) +
                    (metabs.Glucose / params.HK1_K_Glucose) * (metabs.G6P / params.HK1_K_a_G6P_cat) +
                    (metabs.G6P / params.HK1_K_G6P) * (metabs.ADP / params.HK1_K_a_ADP) +
                    (metabs.G6P / params.HK1_K_G6P) * (metabs.G6P / params.HK1_K_a_G6P_cat)
                ) * Z_a_reg + Z_i_cat * (metabs.G6P / params.HK1_K_i_G6P_reg)
            )
        ) / Z
    ATP_bound =
        (
            (params.HK1_Conc / (params.HK1_MW)) * (
                (
                    (metabs.ATP / params.HK1_K_a_ATP) +
                    (params.HK1_β_Glucose_ATP) *
                    (metabs.Glucose / params.HK1_K_Glucose) *
                    (metabs.ATP / params.HK1_K_a_ATP)
                ) * Z_a_reg
            )
        ) / Z
    ADP_bound =
        (
            (params.HK1_Conc / (params.HK1_MW)) * (
                (
                    (metabs.ADP / params.HK1_K_a_ADP) +
                    (metabs.Glucose / params.HK1_K_Glucose) * (metabs.ADP / params.HK1_K_a_ADP) +
                    (metabs.G6P / params.HK1_K_G6P) * (metabs.ADP / params.HK1_K_a_ADP)
                ) * Z_a_reg
            )
        ) / Z
    Phosphate_bound =
        ((params.HK1_Conc / (params.HK1_MW)) * (Z_a_cat * (metabs.Phosphate / params.HK1_K_a_Pi))) / Z
    return (
        Glucose = Glucose_bound,
        G6P = G6P_bound,
        ATP = ATP_bound,
        ADP = ADP_bound,
        Phosphate = Phosphate_bound,
    )
end

function binding_GPI(metabs, params)
    G6P_bound = (
        (params.GPI_Conc / params.GPI_MW) * (metabs.G6P / params.GPI_Km_G6P) /
        (1 + metabs.G6P / params.GPI_Km_G6P + metabs.F6P / params.GPI_Km_F6P)
    )
    F6P_bound = (
        (params.GPI_Conc / params.GPI_MW) * (metabs.F6P / params.GPI_Km_F6P) /
        (1 + metabs.G6P / params.GPI_Km_G6P + metabs.F6P / params.GPI_Km_F6P)
    )
    return (G6P = G6P_bound, F6P = F6P_bound)
end

function binding_PFKP(metabs, params)

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
    Z = (Z_a_cat^4) * (Z_a_reg^4) + params.PFKP_L * (Z_i_cat^4) * (Z_i_reg^4)
    F6P_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) * (
                (metabs.F6P / params.PFKP_K_a_F6P) *
                (1 + (metabs.ATP / params.PFKP_K_ATP)) *
                (Z_a_cat^3) *
                (Z_a_reg^4)
            )
        ) / Z
    ATP_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) * (
                ((metabs.ATP / params.PFKP_K_ATP) * (1 + (metabs.F6P / params.PFKP_K_a_F6P))) *
                (Z_a_cat^3) *
                (Z_a_reg^4) +
                params.PFKP_L *
                (
                    Z_i_reg * (metabs.ATP / params.PFKP_K_ATP) +
                    Z_i_cat * (metabs.ATP / params.PFKP_K_i_ATP_reg)
                ) *
                (Z_i_cat^3) *
                (Z_i_reg^3)
            )
        ) / Z

    F16BP_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) * (
                (metabs.F16BP / params.PFKP_K_F16BP) *
                (1 + (metabs.ADP / params.PFKP_K_ADP)) *
                (Z_a_cat^3) *
                (Z_a_reg^4) +
                params.PFKP_L *
                (metabs.F16BP / params.PFKP_K_F16BP) *
                (1 + (metabs.ADP / params.PFKP_K_ADP)) *
                (Z_i_cat^3) *
                (Z_i_reg^4)
            )
        ) / Z

    ADP_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) * (
                (
                    Z_a_reg * (metabs.ADP / params.PFKP_K_ADP) * (1 + (metabs.F16BP / params.PFKP_K_F16BP)) +
                    Z_a_cat * (metabs.ADP / params.PFKP_K_a_ADP_reg)
                ) *
                (Z_a_cat^3) *
                (Z_a_reg^3) +
                params.PFKP_L *
                (metabs.ADP / params.PFKP_K_ADP) *
                (1 + (metabs.F16BP / params.PFKP_K_F16BP)) *
                (Z_i_cat^3) *
                (Z_i_reg^4)
            )
        ) / Z
    Phosphate_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) *
            (metabs.Phosphate / params.PFKP_K_Phosphate) *
            ((Z_a_cat^4) * (Z_a_reg^3) + params.PFKP_L * (Z_i_cat^4) * (Z_i_reg^3))
        ) / Z
    Citrate_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) *
            (params.PFKP_L * (metabs.Citrate / params.PFKP_K_i_Citrate) * (Z_i_cat^4) * (Z_i_reg^3))
        ) / Z
    F26BP_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) * (
                (metabs.F26BP / params.PFKP_K_a_F26BP) * (Z_a_cat^4) * (Z_a_reg^3) +
                params.PFKP_L * (metabs.F26BP / params.PFKP_K_i_F26BP) * (Z_i_cat^4) * (Z_i_reg^3)
            )
        ) / Z
    return (
        F6P = F6P_bound,
        ATP = ATP_bound,
        F16BP = F16BP_bound,
        ADP = ADP_bound,
        Phosphate = Phosphate_bound,
        Citrate = Citrate_bound,
        F26BP = F26BP_bound,
    )
end

function binding_ALDO(metabs, params)
    denom = (
        1 +
        metabs.GAP * metabs.DHAP / (params.ALDO_Kd_DHAP * params.ALDO_Km_GAP) +
        metabs.DHAP / params.ALDO_Kd_DHAP +
        metabs.F16BP * metabs.GAP / (params.ALDO_Ki_GAP * params.ALDO_Km_F16BP) +
        metabs.F16BP / params.ALDO_Km_F16BP +
        metabs.GAP * params.ALDO_Km_DHAP / (params.ALDO_Kd_DHAP * params.ALDO_Km_GAP)
    )
    F16BP_bound =
        0.5 * (
            (params.ALDO_Conc / params.ALDO_MW) * (
                metabs.GAP * metabs.DHAP / (params.ALDO_Kd_DHAP * params.ALDO_Km_GAP) +
                metabs.F16BP * metabs.GAP / (params.ALDO_Ki_GAP * params.ALDO_Km_F16BP) +
                metabs.F16BP * params.ALDO_Keq * params.ALDO_Km_DHAP /
                (params.ALDO_Kd_DHAP^2 * params.ALDO_Km_GAP)
            )
        ) / denom
    GAP_bound = 0.5 * F16BP_bound
    DHAP_bound =
        0.5 * F16BP_bound + ((params.ALDO_Conc / params.ALDO_MW) * (
            metabs.DHAP / params.ALDO_Kd_DHAP #+ metabs.F16BP / params.ALDO_Km_F16BP - metabs.F16BP * params.ALDO_Keq * params.ALDO_Km_DHAP / (params.ALDO_Kd_DHAP^2 * params.ALDO_Km_GAP)
        ) / denom)
    return (F16BP = F16BP_bound, GAP = GAP_bound, DHAP = DHAP_bound)
end

function binding_TPI(metabs, params)
    GAP_bound = (
        (params.TPI_Conc / params.TPI_MW) * (metabs.GAP / params.TPI_Km_GAP) /
        (1 + metabs.DHAP / params.TPI_Km_DHAP + metabs.GAP / params.TPI_Km_GAP)
    )
    DHAP_bound = (
        (params.TPI_Conc / params.TPI_MW) * (metabs.DHAP / params.TPI_Km_DHAP) /
        (1 + metabs.DHAP / params.TPI_Km_DHAP + metabs.GAP / params.TPI_Km_GAP)
    )
    return (GAP = GAP_bound, DHAP = DHAP_bound)
end

function binding_GAPDH(metabs, params)
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
    Z = (Z_a^4 + params.GAPDH_L * Z_i^4)
    GAP_bound =
        (
            (params.GAPDH_Conc / params.GAPDH_MW) *
            (metabs.GAP / params.GAPDH_K_GAP) *
            (
                (1 + metabs.NAD / params.GAPDH_K_a_NAD + metabs.NADH / params.GAPDH_K_a_NADH) *
                (1 + metabs.Phosphate / params.GAPDH_K_a_Phosphate) *
                Z_a^3 +
                params.GAPDH_L *
                (1 + metabs.NAD / params.GAPDH_K_i_NAD + metabs.NADH / params.GAPDH_K_i_NADH) *
                (1 + metabs.Phosphate / params.GAPDH_K_i_Phosphate) *
                Z_i^3
            )
        ) / Z
    NAD_bound =
        (
            (params.GAPDH_Conc / params.GAPDH_MW) * (
                (metabs.NAD / params.GAPDH_K_a_NAD) *
                (
                    1 +
                    metabs.GAP / params.GAPDH_K_GAP * (1 + metabs.Phosphate / params.GAPDH_K_a_Phosphate) +
                    metabs.BPG / params.GAPDH_K_BPG
                ) *
                Z_a^3 +
                params.GAPDH_L *
                (metabs.NAD / params.GAPDH_K_i_NAD) *
                (
                    1 +
                    metabs.GAP / params.GAPDH_K_GAP * (1 + metabs.Phosphate / params.GAPDH_K_i_Phosphate) +
                    metabs.BPG / params.GAPDH_K_BPG
                ) *
                Z_i^3
            )
        ) / Z
    Phosphate_bound =
        (
            (params.GAPDH_Conc / params.GAPDH_MW) *
            (metabs.Phosphate / params.GAPDH_K_i_Phosphate) *
            (metabs.GAP / params.GAPDH_K_GAP) *
            (
                (1 + metabs.NAD / params.GAPDH_K_a_NAD + metabs.NADH / params.GAPDH_K_a_NADH) * Z_a^3 +
                params.GAPDH_L *
                (1 + metabs.NAD / params.GAPDH_K_i_NAD + metabs.NADH / params.GAPDH_K_i_NADH) *
                Z_i^3
            )
        ) / Z
    BPG_bound =
        (
            (params.GAPDH_Conc / params.GAPDH_MW) *
            (metabs.BPG / params.GAPDH_K_BPG) *
            (
                (1 + metabs.NAD / params.GAPDH_K_a_NAD + metabs.NADH / params.GAPDH_K_a_NADH) *
                (1 + metabs.Phosphate / params.GAPDH_K_a_Phosphate) *
                Z_a^3 +
                params.GAPDH_L *
                (
                    1 +
                    metabs.NAD / params.GAPDH_K_i_NAD +
                    metabs.NADH / (params.GAPDH_α_i_BPG * params.GAPDH_K_i_NADH)
                ) *
                (1 + metabs.Phosphate / params.GAPDH_K_i_Phosphate) *
                Z_i^3
            )
        ) / Z
    NADH_bound =
        (
            (params.GAPDH_Conc / params.GAPDH_MW) * (
                (metabs.NADH / params.GAPDH_K_a_NADH) *
                (
                    1 +
                    metabs.GAP / params.GAPDH_K_GAP * (1 + metabs.Phosphate / params.GAPDH_K_a_Phosphate) +
                    metabs.BPG / params.GAPDH_K_BPG
                ) *
                Z_a^3 +
                params.GAPDH_L *
                (metabs.NADH / params.GAPDH_K_i_NADH) *
                (
                    1 +
                    metabs.GAP / params.GAPDH_K_GAP * (1 + metabs.Phosphate / params.GAPDH_K_i_Phosphate) +
                    metabs.BPG / (params.GAPDH_α_i_BPG * params.GAPDH_K_BPG)
                ) *
                Z_i^3
            )
        ) / Z
    return (GAP = GAP_bound, NAD = NAD_bound, Phosphate = Phosphate_bound, BPG = BPG_bound, NADH = NADH_bound)
end

function binding_PGK(metabs, params)
    denom = (
        1 +
        metabs.BPG / params.PGK_K_BPG +
        metabs.ADP / params.PGK_K_ADP +
        metabs.ThreePG / params.PGK_K_ThreePG +
        metabs.ATP / params.PGK_K_ATP +
        metabs.BPG * metabs.ADP / (params.PGK_α * params.PGK_K_BPG * params.PGK_K_ADP) +
        metabs.ThreePG * metabs.ATP / (params.PGK_β * params.PGK_K_ThreePG * params.PGK_K_ATP) +
        metabs.ThreePG * metabs.ADP / (params.PGK_γ * params.PGK_K_ThreePG * params.PGK_K_ADP)
    )
    BPG_bound =
        (
            (params.PGK_Conc / params.PGK_MW) * (
                metabs.BPG / params.PGK_K_BPG +
                metabs.BPG * metabs.ADP / (params.PGK_α * params.PGK_K_BPG * params.PGK_K_ADP)
            )
        ) / denom
    ADP_bound =
        (
            (params.PGK_Conc / params.PGK_MW) * (
                metabs.ADP / params.PGK_K_ADP +
                metabs.BPG * metabs.ADP / (params.PGK_α * params.PGK_K_BPG * params.PGK_K_ADP) +
                metabs.ThreePG * metabs.ADP / (params.PGK_γ * params.PGK_K_ThreePG * params.PGK_K_ADP)
            )
        ) / denom
    ThreePG_bound =
        (
            (params.PGK_Conc / params.PGK_MW) * (
                metabs.ThreePG / params.PGK_K_ThreePG +
                metabs.ThreePG * metabs.ATP / (params.PGK_β * params.PGK_K_ThreePG * params.PGK_K_ATP) +
                metabs.ThreePG * metabs.ADP / (params.PGK_γ * params.PGK_K_ThreePG * params.PGK_K_ADP)
            )
        ) / denom
    ATP_bound =
        (
            (params.PGK_Conc / params.PGK_MW) * (
                metabs.ATP / params.PGK_K_ATP +
                metabs.ThreePG * metabs.ATP / (params.PGK_β * params.PGK_K_ThreePG * params.PGK_K_ATP)
            )
        ) / denom
    return (BPG = BPG_bound, ADP = ADP_bound, ATP = ATP_bound, ThreePG = ThreePG_bound)
end

function binding_PGM(metabs, params)
    ThreePG_bound = (
        (params.PGM_Conc / params.PGM_MW) * (metabs.ThreePG / params.PGM_Km_ThreePG) /
        (1 + metabs.ThreePG / params.PGM_Km_ThreePG + metabs.TwoPG / params.PGM_Km_TwoPG)
    )
    TwoPG_bound = (
        (params.PGM_Conc / params.PGM_MW) * (metabs.TwoPG / params.PGM_Km_TwoPG) /
        (1 + metabs.ThreePG / params.PGM_Km_ThreePG + metabs.TwoPG / params.PGM_Km_TwoPG)
    )
    return (ThreePG = ThreePG_bound, TwoPG = TwoPG_bound)
end

function binding_ENO(metabs, params)
    TwoPG_bound = (
        (params.ENO_Conc / params.ENO_MW) * (metabs.TwoPG / params.ENO_Km_TwoPG) /
        (1 + metabs.TwoPG / params.ENO_Km_TwoPG + metabs.PEP / params.ENO_Km_PEP)
    )
    PEP_bound = (
        (params.ENO_Conc / params.ENO_MW) * (metabs.PEP / params.ENO_Km_PEP) /
        (1 + metabs.TwoPG / params.ENO_Km_TwoPG + metabs.PEP / params.ENO_Km_PEP)
    )
    return (TwoPG = TwoPG_bound, PEP = PEP_bound)
end

function binding_PKM2(metabs, params)
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
    Z_a_reg =
        ((1 + metabs.F16BP / params.PKM2_K_a_F16BP) * (1 + metabs.Phenylalanine / params.PKM2_K_a_Phenylalanine))
    Z_i_reg = (1 + metabs.Phenylalanine / params.PKM2_K_i_Phenylalanine)
    Z = ((Z_a_cat^4) * (Z_a_reg^4) + params.PKM2_L * (Z_i_cat^4) * (Z_i_reg^4))

    PEP_bound =
        (
            (params.PKM2_Conc / params.PKM2_MW) * (
                (metabs.PEP / params.PKM2_K_a_PEP) *
                (1 + (metabs.ADP / params.PKM2_K_ADP) + (metabs.ATP / params.PKM2_K_ATP)) *
                (Z_a_cat^3) *
                (Z_a_reg^4) +
                params.PKM2_L *
                (metabs.PEP / params.PKM2_K_i_PEP) *
                (
                    1 +
                    (metabs.ADP / params.PKM2_K_ADP) +
                    params.PKM2_β_i_PEP_ATP * (metabs.ATP / params.PKM2_K_ATP)
                ) *
                (Z_i_cat^3) *
                (Z_i_reg^4)
            )
        ) / Z
    ADP_bound =
        (
            (params.PKM2_Conc / params.PKM2_MW) * (
                (metabs.ADP / params.PKM2_K_ADP) *
                (1 + (metabs.PEP / params.PKM2_K_a_PEP) + (metabs.Pyruvate / params.PKM2_K_Pyruvate)) *
                (Z_a_cat^3) *
                (Z_a_reg^4) +
                params.PKM2_L *
                (metabs.ADP / params.PKM2_K_ADP) *
                (1 + (metabs.PEP / params.PKM2_K_i_PEP) + (metabs.Pyruvate / params.PKM2_K_Pyruvate)) *
                (Z_i_cat^3) *
                (Z_i_reg^4)
            )
        ) / Z
    Pyruvate_bound =
        (
            (params.PKM2_Conc / params.PKM2_MW) * (
                (metabs.PEP / params.PKM2_K_a_PEP) *
                (1 + (metabs.ADP / params.PKM2_K_ADP) + (metabs.ATP / params.PKM2_K_ATP)) *
                (Z_a_cat^3) *
                (Z_a_reg^4) +
                params.PKM2_L *
                (metabs.PEP / params.PKM2_K_i_PEP) *
                (1 + (metabs.ADP / params.PKM2_K_ADP) + (metabs.ATP / params.PKM2_K_ATP)) *
                (Z_i_cat^3) *
                (Z_i_reg^4)
            )
        ) / Z
    ATP_bound =
        (
            (params.PKM2_Conc / params.PKM2_MW) * (
                (metabs.ATP / params.PKM2_K_ATP) *
                (1 + (metabs.PEP / params.PKM2_K_a_PEP) + (metabs.Pyruvate / params.PKM2_K_Pyruvate)) *
                (Z_a_cat^3) *
                (Z_a_reg^4) +
                params.PKM2_L *
                (metabs.ATP / params.PKM2_K_ATP) *
                (
                    1 +
                    params.PKM2_β_i_PEP_ATP * (metabs.PEP / params.PKM2_K_i_PEP) +
                    (metabs.Pyruvate / params.PKM2_K_Pyruvate)
                ) *
                (Z_i_cat^3) *
                (Z_i_reg^4)
            )
        ) / Z
    F16BP_bound =
        (
            (params.PKM2_Conc / params.PKM2_MW) *
            ((metabs.F16BP / params.PKM2_K_a_F16BP) * (Z_a_cat^4) * (Z_a_reg^3))
        ) / Z
    Phenylalanine_bound =
        (
            (params.PKM2_Conc / params.PKM2_MW) * (
                ((metabs.Phenylalanine / params.PKM2_K_a_Phenylalanine) * (Z_a_cat^4) * (Z_a_reg^3)) +
                params.PKM2_L * ((metabs.Phenylalanine / params.PKM2_K_i_Phenylalanine) * (Z_i_cat^4) * (Z_i_reg^3))
            )
        ) / Z
    return (
        PEP = PEP_bound,
        ADP = ADP_bound,
        Pyruvate = Pyruvate_bound,
        ATP = ATP_bound,
        F16BP = F16BP_bound,
        Phenylalanine = Phenylalanine_bound,
    )
end
function binding_LDH(metabs, params)
    denom = (
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
    NAD_bound =
        (
            (params.LDH_Conc / params.LDH_MW) * (
                metabs.NAD / params.LDH_Kd_NAD +
                metabs.Pyruvate * metabs.NADH / (params.LDH_Kd_NADH * params.LDH_Km_Pyruvate) +
                metabs.Pyruvate * metabs.NAD * params.LDH_Km_NADH /
                (params.LDH_Kd_NAD * params.LDH_Kd_NADH * params.LDH_Km_Pyruvate)
            )
        ) / denom
    NADH_bound =
        (
            (params.LDH_Conc / params.LDH_MW) * (
                metabs.NADH / params.LDH_Kd_NADH +
                metabs.Lactate * metabs.NAD / (params.LDH_Kd_NAD * params.LDH_Km_Lactate) +
                metabs.Lactate * metabs.NADH * params.LDH_Km_NAD /
                (params.LDH_Kd_NAD * params.LDH_Kd_NADH * params.LDH_Km_Lactate)
            )
        ) / denom

    return (NAD = NAD_bound, NADH = NADH_bound)
end

function binding_MCT(metabs, params)
    Lactate_bound = (
        (params.MCT_Conc / params.MCT_MW) *
        (metabs.Lactate / params.MCT_Km_Lactate + metabs.Lactate_media / params.MCT_Km_Lactate) /
        (1 + metabs.Lactate / params.MCT_Km_Lactate + metabs.Lactate_media / params.MCT_Km_Lactate)
    )
    return (Lactate = Lactate_bound,)
end
