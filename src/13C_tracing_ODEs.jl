function glycolysis_13C_tracing_ODEs(ds, s, params, t)
    Glucose_media = s.Glucose_media_12C + s.Glucose_media_13C
    Glucose = s.Glucose_12C + s.Glucose_13C
    G6P = s.G6P_12C + s.G6P_13C
    F6P = s.F6P_12C + s.F6P_13C
    F16BP = s.F16BP_12C + s.F16BP_13C
    GAP = s.GAP_12C + s.GAP_13C
    DHAP = s.DHAP_12C + s.DHAP_13C
    BPG = s.BPG_12C + s.BPG_13C
    ThreePG = s.ThreePG_12C + s.ThreePG_13C
    TwoPG = s.TwoPG_12C + s.TwoPG_13C
    PEP = s.PEP_12C + s.PEP_13C
    Pyruvate = s.Pyruvate_12C + s.Pyruvate_13C
    Lactate = s.Lactate_12C + s.Lactate_13C
    Lactate_media = s.Lactate_media_12C + s.Lactate_media_13C
    ATP = s.ATP
    ADP = s.ADP
    AMP = s.AMP
    Phosphate = s.Phosphate
    NAD = s.NAD
    NADH = s.NADH
    Citrate = s.Citrate_12C + s.Citrate_13C
    F26BP = s.F26BP_12C + s.F26BP_13C

    ds.Glucose_media_12C = 0
    ds.Glucose_12C =
        rate_13C_GLUT(s.Glucose_media_12C, s.Glucose_12C, Glucose_media, Glucose, params) -
        rate_13C_HK1(s.Glucose_12C, s.G6P_12C, Glucose, G6P, ATP, Phosphate, ADP, params)
    ds.G6P_12C =
        rate_13C_HK1(s.Glucose_12C, s.G6P_12C, Glucose, G6P, ATP, Phosphate, ADP, params) -
        rate_13C_GPI(s.G6P_12C, s.F6P_12C, G6P, F6P, params)
    ds.F6P_12C = (
        rate_13C_GPI(s.G6P_12C, s.F6P_12C, G6P, F6P, params) -
        rate_13C_PFKP(s.F6P_12C, s.F16BP_12C, F6P, ATP, F16BP, ADP, Phosphate, Citrate, F26BP, params)
    )
    ds.F16BP_12C = (
        rate_13C_PFKP(s.F6P_12C, s.F16BP_12C, F6P, ATP, F16BP, ADP, Phosphate, Citrate, F26BP, params) -
        rate_13C_ALDO(s.F16BP_12C, s.GAP_12C, s.DHAP_12C, F16BP, GAP, DHAP, params)
    )
    ds.GAP_12C = (
        rate_13C_ALDO(s.F16BP_12C, s.GAP_12C, s.DHAP_12C, F16BP, GAP, DHAP, params) +
        rate_13C_TPI(s.GAP_12C, s.DHAP_12C, GAP, DHAP, params) -
        rate_13C_GAPDH(s.GAP_12C, s.BPG_12C, GAP, NAD, Phosphate, BPG, NADH, params)
    )
    ds.DHAP_12C =
        rate_13C_ALDO(s.F16BP_12C, s.GAP_12C, s.DHAP_12C, F16BP, GAP, DHAP, params) -
        rate_13C_TPI(s.GAP_12C, s.DHAP_12C, GAP, DHAP, params)
    ds.BPG_12C =
        rate_13C_GAPDH(s.GAP_12C, s.BPG_12C, GAP, NAD, Phosphate, BPG, NADH, params) -
        rate_13C_PGK(s.BPG_12C, s.ThreePG_12C, BPG, ADP, ATP, ThreePG, params)
    ds.ThreePG_12C =
        rate_13C_PGK(s.BPG_12C, s.ThreePG_12C, BPG, ADP, ATP, ThreePG, params) -
        rate_13C_PGM(s.ThreePG_12C, s.TwoPG_12C, ThreePG, TwoPG, params)
    ds.TwoPG_12C =
        rate_13C_PGM(s.ThreePG_12C, s.TwoPG_12C, ThreePG, TwoPG, params) -
        rate_13C_ENO(s.TwoPG_12C, s.PEP_12C, TwoPG, PEP, params)
    ds.PEP_12C =
        rate_13C_ENO(s.TwoPG_12C, s.PEP_12C, TwoPG, PEP, params) -
        rate_13C_PKM2(s.PEP_12C, s.Pyruvate_12C, PEP, ADP, F16BP, ATP, Pyruvate, params)
    ds.Pyruvate_12C =
        rate_13C_PKM2(s.PEP_12C, s.Pyruvate_12C, PEP, ADP, F16BP, ATP, Pyruvate, params) -
        rate_13C_LDH(s.Pyruvate_12C, s.Lactate_12C, Pyruvate, NADH, NAD, Lactate, params)
    ds.Lactate_12C =
        rate_13C_LDH(s.Pyruvate_12C, s.Lactate_12C, Pyruvate, NADH, NAD, Lactate, params) -
        rate_13C_MCT(s.Lactate_12C, s.Lactate_media_12C, Lactate, Lactate_media, params)
    ds.Lactate_media_12C = 0
    ds.F26BP_12C = 0
    ds.Citrate_12C = 0

    ds.Glucose_media_13C = 0
    ds.Glucose_13C =
        rate_13C_GLUT(s.Glucose_media_13C, s.Glucose_13C, Glucose_media, Glucose, params) -
        rate_13C_HK1(s.Glucose_13C, s.G6P_13C, Glucose, G6P, ATP, Phosphate, ADP, params)
    ds.G6P_13C =
        rate_13C_HK1(s.Glucose_13C, s.G6P_13C, Glucose, G6P, ATP, Phosphate, ADP, params) -
        rate_13C_GPI(s.G6P_13C, s.F6P_13C, G6P, F6P, params)
    ds.F6P_13C = (
        rate_13C_GPI(s.G6P_13C, s.F6P_13C, G6P, F6P, params) -
        rate_13C_PFKP(s.F6P_13C, s.F16BP_13C, F6P, ATP, F16BP, ADP, Phosphate, Citrate, F26BP, params)
    )
    ds.F16BP_13C = (
        rate_13C_PFKP(s.F6P_13C, s.F16BP_13C, F6P, ATP, F16BP, ADP, Phosphate, Citrate, F26BP, params) -
        rate_13C_ALDO(s.F16BP_13C, s.GAP_13C, s.DHAP_13C, F16BP, GAP, DHAP, params)
    )
    ds.GAP_13C = (
        rate_13C_ALDO(s.F16BP_13C, s.GAP_13C, s.DHAP_13C, F16BP, GAP, DHAP, params) +
        rate_13C_TPI(s.GAP_13C, s.DHAP_13C, GAP, DHAP, params) -
        rate_13C_GAPDH(s.GAP_13C, s.BPG_13C, GAP, NAD, Phosphate, BPG, NADH, params)
    )
    ds.DHAP_13C =
        rate_13C_ALDO(s.F16BP_13C, s.GAP_13C, s.DHAP_13C, F16BP, GAP, DHAP, params) -
        rate_13C_TPI(s.GAP_13C, s.DHAP_13C, GAP, DHAP, params)
    ds.BPG_13C =
        rate_13C_GAPDH(s.GAP_13C, s.BPG_13C, GAP, NAD, Phosphate, BPG, NADH, params) -
        rate_13C_PGK(s.BPG_13C, s.ThreePG_13C, BPG, ADP, ATP, ThreePG, params)
    ds.ThreePG_13C =
        rate_13C_PGK(s.BPG_13C, s.ThreePG_13C, BPG, ADP, ATP, ThreePG, params) -
        rate_13C_PGM(s.ThreePG_13C, s.TwoPG_13C, ThreePG, TwoPG, params)
    ds.TwoPG_13C =
        rate_13C_PGM(s.ThreePG_13C, s.TwoPG_13C, ThreePG, TwoPG, params) -
        rate_13C_ENO(s.TwoPG_13C, s.PEP_13C, TwoPG, PEP, params)
    ds.PEP_13C =
        rate_13C_ENO(s.TwoPG_13C, s.PEP_13C, TwoPG, PEP, params) -
        rate_13C_PKM2(s.PEP_13C, s.Pyruvate_13C, PEP, ADP, F16BP, ATP, Pyruvate, params)
    ds.Pyruvate_13C =
        rate_13C_PKM2(s.PEP_13C, s.Pyruvate_13C, PEP, ADP, F16BP, ATP, Pyruvate, params) -
        rate_13C_LDH(s.Pyruvate_13C, s.Lactate_13C, Pyruvate, NADH, NAD, Lactate, params)
    ds.Lactate_13C =
        rate_13C_LDH(s.Pyruvate_13C, s.Lactate_13C, Pyruvate, NADH, NAD, Lactate, params) -
        rate_13C_MCT(s.Lactate_13C, s.Lactate_media_13C, Lactate, Lactate_media, params)
    ds.Lactate_media_13C = 0
    ds.F26BP_13C = 0
    ds.Citrate_13C = 0

    ds.ATP = (
        -rate_HK1(Glucose, G6P, ATP, ADP, Phosphate, params) -
        rate_PFKP(F6P, ATP, F16BP, ADP, Phosphate, Citrate, F26BP, params) +
        rate_PGK(BPG, ADP, ATP, ThreePG, params) +
        rate_PKM2(PEP, ADP, F16BP, ATP, Pyruvate, params) - rate_ATPase(ATP, ADP, Phosphate, params) +
        rate_AK(ATP, ADP, AMP, params)
    )
    ds.ADP = (
        rate_HK1(Glucose, G6P, ATP, ADP, Phosphate, params) +
        rate_PFKP(F6P, ATP, F16BP, ADP, Phosphate, Citrate, F26BP, params) -
        rate_PGK(BPG, ADP, ATP, ThreePG, params) - rate_PKM2(PEP, ADP, F16BP, ATP, Pyruvate, params) +
        rate_ATPase(ATP, ADP, Phosphate, params) - 2 * rate_AK(ATP, ADP, AMP, params)
    )
    ds.AMP = rate_AK(ATP, ADP, AMP, params)
    ds.Phosphate =
        rate_ATPase(ATP, ADP, Phosphate, params) - rate_GAPDH(GAP, NAD, Phosphate, BPG, NADH, params)
    ds.NAD =
        rate_LDH(Pyruvate, NADH, NAD, Lactate, params) - rate_GAPDH(GAP, NAD, Phosphate, BPG, NADH, params)
    ds.NADH =
        rate_GAPDH(GAP, NAD, Phosphate, BPG, NADH, params) - rate_LDH(Pyruvate, NADH, NAD, Lactate, params)
end