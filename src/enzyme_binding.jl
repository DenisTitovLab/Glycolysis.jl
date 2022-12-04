function free_to_total_conc(f, params)
    b = similar(f)
    b.Glucose_media = 0.0
    b.Glucose = (
        bindingGLUT(f.Glucose_media, f.Glucose, params).Glucose +
        bindingHK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).Glucose
    )
    b.G6P = (
        bindingHK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).G6P +
        bindingGPI(f.G6P, f.F6P, params).G6P
    )
    b.F6P = (
        bindingGPI(f.G6P, f.F6P, params).F6P +
        bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).F6P
    )
    b.F16BP = (
        bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).F16BP +
        bindingALDO(f.F16BP, f.GAP, f.DHAP, params).F16BP +
        bindingPKM2(f.PEP, f.ADP, f.F16BP, f.ATP, params).F16BP
    )
    b.GAP = (
        bindingALDO(f.F16BP, f.GAP, f.DHAP, params).GAP +
        bindingTPI(f.GAP, f.DHAP, params).GAP +
        bindingGAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).GAP
    )
    b.DHAP = (bindingALDO(f.F16BP, f.GAP, f.DHAP, params).DHAP + bindingTPI(f.GAP, f.DHAP, params).DHAP)
    b.BPG = (
        bindingGAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).BPG +
        bindingPGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).BPG
    )
    b.ThreePG = (
        bindingPGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).ThreePG +
        bindingPGM(f.ThreePG, f.TwoPG, params).ThreePG
    )
    b.TwoPG = (bindingPGM(f.ThreePG, f.TwoPG, params).TwoPG + bindingENO(f.TwoPG, f.PEP, params).TwoPG)
    b.PEP = (bindingENO(f.TwoPG, f.PEP, params).PEP + bindingPKM2(f.PEP, f.ADP, f.F16BP, f.ATP, params).PEP)
    b.Pyruvate = 0.0
    b.Lactate = bindingMCT(f.Lactate, f.Lactate_media, params).Lactate
    b.Lactate_media = 0.0
    b.ATP = (
        bindingHK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).ATP +
        bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).ATP +
        bindingPGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).ATP +
        bindingPKM2(f.PEP, f.ADP, f.F16BP, f.ATP, params).ATP
    )
    b.ADP = (
        bindingHK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).ADP +
        bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).ADP +
        bindingPGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).ADP +
        bindingPKM2(f.PEP, f.ADP, f.F16BP, f.ATP, params).ADP
    )
    b.AMP = 0
    b.Phosphate = (
        bindingGAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).Phosphate +
        bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).Phosphate +
        bindingHK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).Phosphate
    )
    b.NAD = (
        bindingLDH(f.Pyruvate, f.NADH, f.NAD, f.Lactate, params).NAD +
        bindingGAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).NAD
    )
    b.NADH = (
        bindingGAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).NADH +
        bindingLDH(f.Pyruvate, f.NADH, f.NAD, f.Lactate, params).NADH
    )
    b.F26BP = (bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).F26BP)
    b.Citrate = (bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).Citrate)
    return (b + f) .* cell_volume_correction
end

function free_to_bound_conc(f, params)
    b = similar(f)
    b.Glucose_media = 0.0
    b.Glucose = (
        bindingGLUT(f.Glucose_media, f.Glucose, params).Glucose +
        bindingHK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).Glucose
    )
    b.G6P = (
        bindingHK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).G6P +
        bindingGPI(f.G6P, f.F6P, params).G6P
    )
    b.F6P = (
        bindingGPI(f.G6P, f.F6P, params).F6P +
        bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).F6P
    )
    b.F16BP = (
        bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).F16BP +
        bindingALDO(f.F16BP, f.GAP, f.DHAP, params).F16BP +
        bindingPKM2(f.PEP, f.ADP, f.F16BP, f.ATP, params).F16BP
    )
    b.GAP = (
        bindingALDO(f.F16BP, f.GAP, f.DHAP, params).GAP +
        bindingTPI(f.GAP, f.DHAP, params).GAP +
        bindingGAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).GAP
    )
    b.DHAP = (ALDO(f.F16BP, f.GAP, f.DHAP, params).DHAP + TPI(f.GAP, f.DHAP, params).DHAP)
    b.BPG = (
        bindingGAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).BPG +
        bindingPGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).BPG
    )
    b.ThreePG = (
        bindingPGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).ThreePG +
        bindingPGM(f.ThreePG, f.TwoPG, params).ThreePG
    )
    b.TwoPG = (bindingPGM(f.ThreePG, f.TwoPG, params).TwoPG + bindingENO(f.TwoPG, f.PEP, params).TwoPG)
    b.PEP = (bindingENO(f.TwoPG, f.PEP, params).PEP + bindingPKM2(f.PEP, f.ADP, f.F16BP, f.ATP, params).PEP)
    b.Pyruvate = 0.0
    b.Lactate = bindingMCT(f.Lactate, f.Lactate_media, params).Lactate
    b.Lactate_media = 0.0
    b.ATP = (
        bindingHK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).ATP +
        bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).ATP +
        bindingPGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).ATP +
        bindingPKM2(f.PEP, f.ADP, f.F16BP, f.ATP, params).ATP
    )
    b.ADP = (
        bindingHK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).ADP +
        bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).ADP +
        bindingPGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).ADP +
        bindingPKM2(f.PEP, f.ADP, f.F16BP, f.ATP, params).ADP
    )
    b.AMP = 0.0
    b.Phosphate = (
        bindingGAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).Phosphate +
        bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).Phosphate +
        bindingHK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).Phosphate
    )
    b.NAD = (
        bindingLDH(f.Pyruvate, f.NADH, f.NAD, f.Lactate, params).NAD +
        bindingGAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).NAD
    )
    b.NADH = (
        bindingGAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).NADH +
        bindingLDH(f.Pyruvate, f.NADH, f.NAD, f.Lactate, params).NADH
    )
    b.F26BP = (bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).F26BP)
    b.Citrate = (bindingPFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).Citrate)
    return b
end

function bindingGLUT(Glucose_media, Glucose, params)
    Glucose_bound = (
        (params.GLUT_Conc / params.GLUT_MW) *
        (Glucose_media / params.GLUT_Km_Glucose + Glucose / params.GLUT_Km_Glucose) /
        (1 + Glucose_media / params.GLUT_Km_Glucose + Glucose / params.GLUT_Km_Glucose)
    )
    return (Glucose = Glucose_bound,)
end

function bindingHK1(Glucose, G6P, ATP, ADP, Phosphate, params)
    Z_a_cat = (
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
    )
    Z_i_cat = (1 + (Glucose / params.HK1_K_Glucose) + (G6P / params.HK1_K_G6P))
    Z_a_reg = (1 + (Phosphate / params.HK1_K_a_Pi))
    Z_i_reg = (G6P / params.HK1_K_i_G6P_reg)
    Z = Z_a_cat * Z_a_reg + Z_i_cat * Z_i_reg

    Glucose_bound =
        (
            (params.HK1_Conc / (params.HK1_MW)) *
            ((Glucose / params.HK1_K_Glucose) * Z_a_reg + (Glucose / params.HK1_K_Glucose) * Z_i_reg)
        ) / Z
    G6P_bound =
        (
            (params.HK1_Conc / (params.HK1_MW)) * (
                (
                    (G6P / params.HK1_K_G6P) +
                    (G6P / params.HK1_K_a_G6P_cat) +
                    (Glucose / params.HK1_K_Glucose) * (G6P / params.HK1_K_a_G6P_cat) +
                    (G6P / params.HK1_K_G6P) * (ADP / params.HK1_K_a_ADP) +
                    (G6P / params.HK1_K_G6P) * (G6P / params.HK1_K_a_G6P_cat)
                ) * Z_a_reg + Z_i_cat * (G6P / params.HK1_K_i_G6P_reg)
            )
        ) / Z
    ATP_bound =
        (
            (params.HK1_Conc / (params.HK1_MW)) * (
                (
                    (ATP / params.HK1_K_a_ATP) +
                    (params.HK1_β_Glucose_ATP) *
                    (Glucose / params.HK1_K_Glucose) *
                    (ATP / params.HK1_K_a_ATP)
                ) * Z_a_reg
            )
        ) / Z
    ADP_bound =
        (
            (params.HK1_Conc / (params.HK1_MW)) * (
                (
                    (ADP / params.HK1_K_a_ADP) +
                    (Glucose / params.HK1_K_Glucose) * (ADP / params.HK1_K_a_ADP) +
                    (G6P / params.HK1_K_G6P) * (ADP / params.HK1_K_a_ADP)
                ) * Z_a_reg
            )
        ) / Z
    Phosphate_bound = ((params.HK1_Conc / (params.HK1_MW)) * (Z_a_cat * (Phosphate / params.HK1_K_a_Pi))) / Z
    return (
        Glucose = Glucose_bound,
        G6P = G6P_bound,
        ATP = ATP_bound,
        ADP = ADP_bound,
        Phosphate = Phosphate_bound,
    )
end

function bindingGPI(G6P, F6P, params)
    G6P_bound = (
        (params.GPI_Conc / params.GPI_MW) * (G6P / params.GPI_Km_G6P) /
        (1 + G6P / params.GPI_Km_G6P + F6P / params.GPI_Km_F6P)
    )
    F6P_bound = (
        (params.GPI_Conc / params.GPI_MW) * (F6P / params.GPI_Km_F6P) /
        (1 + G6P / params.GPI_Km_G6P + F6P / params.GPI_Km_F6P)
    )
    return (G6P = G6P_bound, F6P = F6P_bound)
end

function bindingPFKP(F6P, ATP, F16BP, ADP, Phosphate, Citrate, F26BP, params)

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
    Z = (Z_a_cat^4) * (Z_a_reg^4) + params.PFKP_L * (Z_i_cat^4) * (Z_i_reg^4)
    F6P_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) *
            ((F6P / params.PFKP_K_a_F6P) * (1 + (ATP / params.PFKP_K_ATP)) * (Z_a_cat^3) * (Z_a_reg^4))
        ) / Z
    ATP_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) * (
                ((ATP / params.PFKP_K_ATP) * (1 + (F6P / params.PFKP_K_a_F6P))) * (Z_a_cat^3) * (Z_a_reg^4) +
                params.PFKP_L *
                (Z_i_reg * (ATP / params.PFKP_K_ATP) + Z_i_cat * (ATP / params.PFKP_K_i_ATP_reg)) *
                (Z_i_cat^3) *
                (Z_i_reg^3)
            )
        ) / Z

    F16BP_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) * (
                (F16BP / params.PFKP_K_F16BP) * (1 + (ADP / params.PFKP_K_ADP)) * (Z_a_cat^3) * (Z_a_reg^4) +
                params.PFKP_L *
                (F16BP / params.PFKP_K_F16BP) *
                (1 + (ADP / params.PFKP_K_ADP)) *
                (Z_i_cat^3) *
                (Z_i_reg^4)
            )
        ) / Z

    ADP_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) * (
                (
                    Z_a_reg * (ADP / params.PFKP_K_ADP) * (1 + (F16BP / params.PFKP_K_F16BP)) +
                    Z_a_cat * (ADP / params.PFKP_K_a_ADP_reg)
                ) *
                (Z_a_cat^3) *
                (Z_a_reg^3) +
                params.PFKP_L *
                (ADP / params.PFKP_K_ADP) *
                (1 + (F16BP / params.PFKP_K_F16BP)) *
                (Z_i_cat^3) *
                (Z_i_reg^4)
            )
        ) / Z
    Phosphate_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) *
            (Phosphate / params.PFKP_K_Phosphate) *
            ((Z_a_cat^4) * (Z_a_reg^3) + params.PFKP_L * (Z_i_cat^4) * (Z_i_reg^3))
        ) / Z
    Citrate_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) *
            (params.PFKP_L * (Citrate / params.PFKP_K_i_Citrate) * (Z_i_cat^4) * (Z_i_reg^3))
        ) / Z
    F26BP_bound =
        (
            (params.PFKP_Conc / (params.PFKP_MW)) * (
                (F26BP / params.PFKP_K_a_F26BP) * (Z_a_cat^4) * (Z_a_reg^3) +
                params.PFKP_L * (F26BP / params.PFKP_K_i_F26BP) * (Z_i_cat^4) * (Z_i_reg^3)
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

function bindingALDO(F16BP, GAP, DHAP, params)
    denom = (
        1 +
        GAP * DHAP / (params.ALDO_Kd_DHAP * params.ALDO_Km_GAP) +
        DHAP / params.ALDO_Kd_DHAP +
        F16BP * GAP / (params.ALDO_Ki_GAP * params.ALDO_Km_F16BP) +
        F16BP / params.ALDO_Km_F16BP +
        GAP * params.ALDO_Km_DHAP / (params.ALDO_Kd_DHAP * params.ALDO_Km_GAP)
    )
    F16BP_bound =
        0.5 * (
            (params.ALDO_Conc / params.ALDO_MW) * (
                GAP * DHAP / (params.ALDO_Kd_DHAP * params.ALDO_Km_GAP) +
                F16BP * GAP / (params.ALDO_Ki_GAP * params.ALDO_Km_F16BP) +
                F16BP * params.ALDO_Keq * params.ALDO_Km_DHAP / (params.ALDO_Kd_DHAP^2 * params.ALDO_Km_GAP)
            )
        ) / denom
    GAP_bound = 0.5 * F16BP_bound
    DHAP_bound =
        0.5 * F16BP_bound + ((params.ALDO_Conc / params.ALDO_MW) * (
            DHAP / params.ALDO_Kd_DHAP #+ F16BP / params.ALDO_Km_F16BP - F16BP * params.ALDO_Keq * params.ALDO_Km_DHAP / (params.ALDO_Kd_DHAP^2 * params.ALDO_Km_GAP)
        ) / denom)
    return (F16BP = F16BP_bound, GAP = GAP_bound, DHAP = DHAP_bound)
end

function bindingTPI(GAP, DHAP, params)
    GAP_bound = (
        (params.TPI_Conc / params.TPI_MW) * (GAP / params.TPI_Km_GAP) /
        (1 + DHAP / params.TPI_Km_DHAP + GAP / params.TPI_Km_GAP)
    )
    DHAP_bound = (
        (params.TPI_Conc / params.TPI_MW) * (DHAP / params.TPI_Km_DHAP) /
        (1 + DHAP / params.TPI_Km_DHAP + GAP / params.TPI_Km_GAP)
    )
    return (GAP = GAP_bound, DHAP = DHAP_bound)
end

function bindingGAPDH(GAP, NAD, Phosphate, BPG, NADH, params)
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
    Z = (Z_a^4 + params.GAPDH_L * Z_i^4)
    GAP_bound =
        (
            (params.GAPDH_Conc / params.GAPDH_MW) *
            (GAP / params.GAPDH_K_GAP) *
            (
                (1 + NAD / params.GAPDH_K_a_NAD + NADH / params.GAPDH_K_a_NADH) *
                (1 + Phosphate / params.GAPDH_K_a_Phosphate) *
                Z_a^3 +
                params.GAPDH_L *
                (1 + NAD / params.GAPDH_K_i_NAD + NADH / params.GAPDH_K_i_NADH) *
                (1 + Phosphate / params.GAPDH_K_i_Phosphate) *
                Z_i^3
            )
        ) / Z
    NAD_bound =
        (
            (params.GAPDH_Conc / params.GAPDH_MW) * (
                (NAD / params.GAPDH_K_a_NAD) *
                (
                    1 +
                    GAP / params.GAPDH_K_GAP * (1 + Phosphate / params.GAPDH_K_a_Phosphate) +
                    BPG / params.GAPDH_K_BPG
                ) *
                Z_a^3 +
                params.GAPDH_L *
                (NAD / params.GAPDH_K_i_NAD) *
                (
                    1 +
                    GAP / params.GAPDH_K_GAP * (1 + Phosphate / params.GAPDH_K_i_Phosphate) +
                    BPG / params.GAPDH_K_BPG
                ) *
                Z_i^3
            )
        ) / Z
    Phosphate_bound =
        (
            (params.GAPDH_Conc / params.GAPDH_MW) *
            (Phosphate / params.GAPDH_K_i_Phosphate) *
            (GAP / params.GAPDH_K_GAP) *
            (
                (1 + NAD / params.GAPDH_K_a_NAD + NADH / params.GAPDH_K_a_NADH) * Z_a^3 +
                params.GAPDH_L * (1 + NAD / params.GAPDH_K_i_NAD + NADH / params.GAPDH_K_i_NADH) * Z_i^3
            )
        ) / Z
    BPG_bound =
        (
            (params.GAPDH_Conc / params.GAPDH_MW) *
            (BPG / params.GAPDH_K_BPG) *
            (
                (1 + NAD / params.GAPDH_K_a_NAD + NADH / params.GAPDH_K_a_NADH) *
                (1 + Phosphate / params.GAPDH_K_a_Phosphate) *
                Z_a^3 +
                params.GAPDH_L *
                (1 + NAD / params.GAPDH_K_i_NAD + NADH / (params.GAPDH_α_i_BPG * params.GAPDH_K_i_NADH)) *
                (1 + Phosphate / params.GAPDH_K_i_Phosphate) *
                Z_i^3
            )
        ) / Z
    NADH_bound =
        (
            (params.GAPDH_Conc / params.GAPDH_MW) * (
                (NADH / params.GAPDH_K_a_NADH) *
                (
                    1 +
                    GAP / params.GAPDH_K_GAP * (1 + Phosphate / params.GAPDH_K_a_Phosphate) +
                    BPG / params.GAPDH_K_BPG
                ) *
                Z_a^3 +
                params.GAPDH_L *
                (NADH / params.GAPDH_K_i_NADH) *
                (
                    1 +
                    GAP / params.GAPDH_K_GAP * (1 + Phosphate / params.GAPDH_K_i_Phosphate) +
                    BPG / (params.GAPDH_α_i_BPG * params.GAPDH_K_BPG)
                ) *
                Z_i^3
            )
        ) / Z
    return (GAP = GAP_bound, NAD = NAD_bound, Phosphate = Phosphate_bound, BPG = BPG_bound, NADH = NADH_bound)
end

function bindingPGK(BPG, ADP, ATP, ThreePG, params)
    denom = (
        1 +
        BPG / params.PGK_K_BPG +
        ADP / params.PGK_K_ADP +
        ThreePG / params.PGK_K_ThreePG +
        ATP / params.PGK_K_ATP +
        BPG * ADP / (params.PGK_α * params.PGK_K_BPG * params.PGK_K_ADP) +
        ThreePG * ATP / (params.PGK_β * params.PGK_K_ThreePG * params.PGK_K_ATP) +
        ThreePG * ADP / (params.PGK_γ * params.PGK_K_ThreePG * params.PGK_K_ADP)
    )
    BPG_bound =
        (
            (params.PGK_Conc / params.PGK_MW) *
            (BPG / params.PGK_K_BPG + BPG * ADP / (params.PGK_α * params.PGK_K_BPG * params.PGK_K_ADP))
        ) / denom
    ADP_bound =
        (
            (params.PGK_Conc / params.PGK_MW) * (
                ADP / params.PGK_K_ADP +
                BPG * ADP / (params.PGK_α * params.PGK_K_BPG * params.PGK_K_ADP) +
                ThreePG * ADP / (params.PGK_γ * params.PGK_K_ThreePG * params.PGK_K_ADP)
            )
        ) / denom
    ThreePG_bound =
        (
            (params.PGK_Conc / params.PGK_MW) * (
                ThreePG / params.PGK_K_ThreePG +
                ThreePG * ATP / (params.PGK_β * params.PGK_K_ThreePG * params.PGK_K_ATP) +
                ThreePG * ADP / (params.PGK_γ * params.PGK_K_ThreePG * params.PGK_K_ADP)
            )
        ) / denom
    ATP_bound =
        (
            (params.PGK_Conc / params.PGK_MW) * (
                ATP / params.PGK_K_ATP +
                ThreePG * ATP / (params.PGK_β * params.PGK_K_ThreePG * params.PGK_K_ATP)
            )
        ) / denom
    return (BPG = BPG_bound, ADP = ADP_bound, ATP = ATP_bound, ThreePG = ThreePG_bound)
end

function bindingPGM(ThreePG, TwoPG, params)
    ThreePG_bound = (
        (params.PGM_Conc / params.PGM_MW) * (ThreePG / params.PGM_Km_ThreePG) /
        (1 + ThreePG / params.PGM_Km_ThreePG + TwoPG / params.PGM_Km_TwoPG)
    )
    TwoPG_bound = (
        (params.PGM_Conc / params.PGM_MW) * (TwoPG / params.PGM_Km_TwoPG) /
        (1 + ThreePG / params.PGM_Km_ThreePG + TwoPG / params.PGM_Km_TwoPG)
    )
    return (ThreePG = ThreePG_bound, TwoPG = TwoPG_bound)
end

function bindingENO(TwoPG, PEP, params)
    TwoPG_bound = (
        (params.ENO_Conc / params.ENO_MW) * (TwoPG / params.ENO_Km_TwoPG) /
        (1 + TwoPG / params.ENO_Km_TwoPG + PEP / params.ENO_Km_PEP)
    )
    PEP_bound = (
        (params.ENO_Conc / params.ENO_MW) * (PEP / params.ENO_Km_PEP) /
        (1 + TwoPG / params.ENO_Km_TwoPG + PEP / params.ENO_Km_PEP)
    )
    return (TwoPG = TwoPG_bound, PEP = PEP_bound)
end

function bindingPKM2(PEP, ADP, F16BP, ATP, params)
    Z = (
        ((1 + PEP / params.PKM2_a_KmPEP)^params.PKM2_n) *
        ((1 + ADP / params.PKM2_a_KmADP + ATP / params.PKM2_a_KdATP)^params.PKM2_n) *
        ((1 + F16BP / params.PKM2_a_KdF16BP)^params.PKM2_n) +
        params.PKM2_L *
        ((1 + PEP / params.PKM2_i_KmPEP)^params.PKM2_n) *
        ((1 + ADP / params.PKM2_i_KmADP + ATP / params.PKM2_i_KdATP)^params.PKM2_n) *
        ((1 + F16BP / params.PKM2_i_KdF16BP)^params.PKM2_n)
    )
    Pa = (
        ((1 + PEP / params.PKM2_a_KmPEP)^params.PKM2_n) *
        ((1 + ADP / params.PKM2_a_KmADP + ATP / params.PKM2_a_KdATP)^params.PKM2_n) *
        ((1 + F16BP / params.PKM2_a_KdF16BP)^params.PKM2_n) / Z
    )
    Pi = (
        params.PKM2_L *
        ((1 + PEP / params.PKM2_i_KmPEP)^params.PKM2_n) *
        ((1 + ADP / params.PKM2_i_KmADP + ATP / params.PKM2_i_KdATP)^params.PKM2_n) *
        ((1 + F16BP / params.PKM2_i_KdF16BP)^params.PKM2_n) / Z
    )
    PEP_bound = (
        (params.PKM2_Conc / params.PKM2_MW) * (
            Pa * (PEP / params.PKM2_a_KmPEP) / (1 + PEP / params.PKM2_a_KmPEP) +
            Pi * (PEP / params.PKM2_i_KmPEP) / (1 + PEP / params.PKM2_i_KmPEP)
        )
    )
    ADP_bound = (
        (params.PKM2_Conc / params.PKM2_MW) * (
            Pa * (ADP / params.PKM2_a_KmADP) / (1 + ADP / params.PKM2_a_KmADP + ATP / params.PKM2_a_KdATP) +
            Pi * (ADP / params.PKM2_i_KmADP) / (1 + ADP / params.PKM2_i_KmADP + ATP / params.PKM2_i_KdATP)
        )
    )
    ATP_bound = (
        (params.PKM2_Conc / params.PKM2_MW) * (
            Pa * (ATP / params.PKM2_a_KdATP) / (1 + ADP / params.PKM2_a_KmADP + ATP / params.PKM2_a_KdATP) +
            Pi * (ATP / params.PKM2_i_KdATP) / (1 + ADP / params.PKM2_i_KmADP + ATP / params.PKM2_i_KdATP)
        )
    )
    F16BP_bound = (
        (params.PKM2_Conc / params.PKM2_MW) * (
            Pa * (F16BP / params.PKM2_a_KdF16BP) / (1 + F16BP / params.PKM2_a_KdF16BP) +
            Pi * (F16BP / params.PKM2_i_KdF16BP) / (1 + F16BP / params.PKM2_i_KdF16BP)
        )
    )
    return (PEP = PEP_bound, ADP = ADP_bound, F16BP = F16BP_bound, ATP = ATP_bound)
end

function bindingLDH(Pyruvate, NADH, NAD, Lactate, params)
    denom = (
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
    NAD_bound =
        (
            (params.LDH_Conc / params.LDH_MW) * (
                NAD / params.LDH_Kd_NAD +
                Pyruvate * NADH / (params.LDH_Kd_NADH * params.LDH_Km_Pyruvate) +
                Pyruvate * NAD * params.LDH_Km_NADH /
                (params.LDH_Kd_NAD * params.LDH_Kd_NADH * params.LDH_Km_Pyruvate)
            )
        ) / denom
    NADH_bound =
        (
            (params.LDH_Conc / params.LDH_MW) * (
                NADH / params.LDH_Kd_NADH +
                Lactate * NAD / (params.LDH_Kd_NAD * params.LDH_Km_Lactate) +
                Lactate * NADH * params.LDH_Km_NAD /
                (params.LDH_Kd_NAD * params.LDH_Kd_NADH * params.LDH_Km_Lactate)
            )
        ) / denom

    return (NAD = NAD_bound, NADH = NADH_bound)
end

function bindingMCT(Lactate, Lactate_media, params)
    Lactate_bound = (
        (params.MCT_Conc / params.MCT_MW) *
        (Lactate / params.MCT_Km_Lactate + Lactate_media / params.MCT_Km_Lactate) /
        (1 + Lactate / params.MCT_Km_Lactate + Lactate_media / params.MCT_Km_Lactate)
    )
    return (Lactate = Lactate_bound,)
end