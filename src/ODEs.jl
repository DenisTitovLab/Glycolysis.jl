#Glycolysis Model ODE system

function glycolysis_ODEs(ds, s, params, t)
    ds.Glucose_media = 0.0
    ds.Glucose =
        rate_GLUT(s.Glucose_media, s.Glucose, params) -
        rate_HK1(s.Glucose, s.G6P, s.ATP, s.ADP, s.Phosphate, params)
    ds.G6P = rate_HK1(s.Glucose, s.G6P, s.ATP, s.ADP, s.Phosphate, params) - rate_GPI(s.G6P, s.F6P, params)
    ds.F6P = (
        rate_GPI(s.G6P, s.F6P, params) -
        rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params)
    )
    ds.F16BP = (
        rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) -
        rate_ALDO(s.F16BP, s.GAP, s.DHAP, params)
    )
    ds.GAP = (
        rate_ALDO(s.F16BP, s.GAP, s.DHAP, params) + rate_TPI(s.GAP, s.DHAP, params) -
        rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    )
    ds.DHAP = rate_ALDO(s.F16BP, s.GAP, s.DHAP, params) - rate_TPI(s.GAP, s.DHAP, params)
    ds.BPG =
        rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params) -
        rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params)
    ds.ThreePG = rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) - rate_PGM(s.ThreePG, s.TwoPG, params)
    ds.TwoPG = rate_PGM(s.ThreePG, s.TwoPG, params) - rate_ENO(s.TwoPG, s.PEP, params)
    ds.PEP = rate_ENO(s.TwoPG, s.PEP, params) - rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params)
    ds.Pyruvate =
        rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) -
        rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params)
    ds.Lactate =
        rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params) - rate_MCT(s.Lactate, s.Lactate_media, params)
    ds.Lactate_media = 0.0
    ds.ATP = (
        -rate_HK1(s.Glucose, s.G6P, s.ATP, s.ADP, s.Phosphate, params) -
        rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) +
        rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) +
        rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) -
        rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) + rate_AK(s.ATP, s.ADP, s.AMP, params) -
        rate_NDPK(s.NTP, s.NDP, s.ATP, s.ADP, params) -
        rate_CK(s.Phosphocreatine, s.Creatine, s.ATP, s.ADP, params)
    )
    ds.ADP = (
        rate_HK1(s.Glucose, s.G6P, s.ATP, s.ADP, s.Phosphate, params) +
        rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) -
        rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) -
        rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) +
        rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) - 2 * rate_AK(s.ATP, s.ADP, s.AMP, params) +
        rate_NDPK(s.NTP, s.NDP, s.ATP, s.ADP, params) +
        rate_CK(s.Phosphocreatine, s.Creatine, s.ATP, s.ADP, params)
    )
    ds.AMP = rate_AK(s.ATP, s.ADP, s.AMP, params)
    ds.Phosphate =
        rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) -
        rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    #ds.Phosphate = 0
    ds.NTP = rate_NDPK(s.NTP, s.NDP, s.ATP, s.ADP, params)
    ds.NDP = -rate_NDPK(s.NTP, s.NDP, s.ATP, s.ADP, params)
    ds.Phosphocreatine = rate_CK(s.Phosphocreatine, s.Creatine, s.ATP, s.ADP, params)
    ds.Creatine = -rate_CK(s.Phosphocreatine, s.Creatine, s.ATP, s.ADP, params)
    ds.NAD =
        rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params) -
        rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    ds.NADH =
        rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params) -
        rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params)
    ds.F26BP = 0.0
    ds.Citrate = 0.0
    ds.Phenylalanine = 0.0
end
