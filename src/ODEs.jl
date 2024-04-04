#Glycolysis Model ODE system

function glycolysis_ODEs(ds, s, params, t)
    ds.Glucose_media = 0.0
    ds.Glucose = rate_GLUT(s, params) - rate_HK1(s, params)
    ds.G6P = rate_HK1(s, params) - rate_GPI(s, params)
    ds.F6P = (rate_GPI(s, params) - rate_PFKP(s, params))
    ds.F16BP = (rate_PFKP(s, params) - rate_ALDO(s, params))
    ds.GAP = (rate_ALDO(s, params) + rate_TPI(s, params) - rate_GAPDH(s, params))
    ds.DHAP = rate_ALDO(s, params) - rate_TPI(s, params)
    ds.BPG = rate_GAPDH(s, params) - rate_PGK(s, params)
    ds.ThreePG = rate_PGK(s, params) - rate_PGM(s, params)
    ds.TwoPG = rate_PGM(s, params) - rate_ENO(s, params)
    ds.PEP = rate_ENO(s, params) - rate_PKM2(s, params)
    ds.Pyruvate = rate_PKM2(s, params) - rate_LDH(s, params)
    ds.Lactate = rate_LDH(s, params) - rate_MCT(s, params)
    ds.Lactate_media = 0.0
    ds.ATP = (
        -rate_HK1(s, params) - rate_PFKP(s, params) + rate_PGK(s, params) + rate_PKM2(s, params) -
        rate_ATPase(s, params) + rate_AK(s, params) - rate_NDPK(s, params) - rate_CK(s, params)
    )
    ds.ADP = (
        rate_HK1(s, params) + rate_PFKP(s, params) - rate_PGK(s, params) - rate_PKM2(s, params) +
        rate_ATPase(s, params) - 2 * rate_AK(s, params) +
        rate_NDPK(s, params) +
        rate_CK(s, params)
    )
    ds.AMP = rate_AK(s, params)
    ds.Phosphate = rate_ATPase(s, params) - rate_GAPDH(s, params)
    ds.NTP = rate_NDPK(s, params)
    ds.NDP = -rate_NDPK(s, params)
    ds.Phosphocreatine = rate_CK(s, params)
    ds.Creatine = -rate_CK(s, params)
    ds.NAD = rate_LDH(s, params) - rate_GAPDH(s, params)
    ds.NADH = rate_GAPDH(s, params) - rate_LDH(s, params)
    ds.F26BP = 0.0
    ds.Citrate = 0.0
    ds.Phenylalanine = 0.0
end
