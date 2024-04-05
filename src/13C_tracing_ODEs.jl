function glycolysis_13C_tracing_ODEs(ds, s, params, t)
    total_conc_larray = LArrays((
        Glucose_media = s.Glucose_media_12C + s.Glucose_media_13C,
        Glucose = s.Glucose_12C + s.Glucose_13C,
        G6P = s.G6P_12C + s.G6P_13C,
        F6P = s.F6P_12C + s.F6P_13C,
        F16BP = s.F16BP_12C + s.F16BP_13C,
        GAP = s.GAP_12C + s.GAP_13C,
        DHAP = s.DHAP_12C + s.DHAP_13C,
        BPG = s.BPG_12C + s.BPG_13C,
        ThreePG = s.ThreePG_12C + s.ThreePG_13C,
        TwoPG = s.TwoPG_12C + s.TwoPG_13C,
        PEP = s.PEP_12C + s.PEP_13C,
        Pyruvate = s.Pyruvate_12C + s.Pyruvate_13C,
        Lactate = s.Lactate_12C + s.Lactate_13C,
        Lactate_media = s.Lactate_media_12C + s.Lactate_media_13C,
        ATP = s.ATP,
        ADP = s.ADP,
        AMP = s.AMP,
        Phosphate = s.Phosphate,
        NAD = s.NAD,
        NADH = s.NADH,
        Citrate = s.Citrate_12C + s.Citrate_13C,
        F26BP = s.F26BP_12C + s.F26BP_13C,
        Phenylalanine = s.Phenylalanine_12C + s.Phenylalanine_13C,
    ))

    ds.Glucose_media_12C = 0.0
    ds.Glucose_12C = rate_isotope_tracing_GLUT(s, params, :C12) - rate_isotope_tracing_HK1(s, params, :C12)
    ds.G6P_12C = rate_isotope_tracing_HK1(s, params, :C12) - rate_isotope_tracing_GPI(s, params, :C12)
    ds.F6P_12C = (rate_isotope_tracing_GPI(s, params, :C12) - rate_isotope_tracing_PFKP(s, params, :C12))
    ds.F16BP_12C = (rate_isotope_tracing_PFKP(s, params, :C12) - rate_isotope_tracing_ALDO(s, params, :C12))
    ds.GAP_12C = (
        rate_isotope_tracing_ALDO(s, params, :C12) + rate_isotope_tracing_TPI(s, params, :C12) -
        rate_isotope_tracing_GAPDH(s, params, :C12)
    )
    ds.DHAP_12C = rate_isotope_tracing_ALDO(s, params, :C12) - rate_isotope_tracing_TPI(s, params, :C12)
    ds.BPG_12C = rate_isotope_tracing_GAPDH(s, params, :C12) - rate_isotope_tracing_PGK(s, params, :C12)
    ds.ThreePG_12C = rate_isotope_tracing_PGK(s, params, :C12) - rate_isotope_tracing_PGM(s, params, :C12)
    ds.TwoPG_12C = rate_isotope_tracing_PGM(s, params, :C12) - rate_isotope_tracing_ENO(s, params, :C12)
    ds.PEP_12C = rate_isotope_tracing_ENO(s, params, :C12) - rate_isotope_tracing_PKM2(s, params, :C12)
    ds.Pyruvate_12C = rate_isotope_tracing_PKM2(s, params, :C12) - rate_isotope_tracing_LDH(s, params, :C12)
    ds.Lactate_12C = rate_isotope_tracing_LDH(s, params, :C12) - rate_isotope_tracing_MCT(s, params, :C12)
    ds.Lactate_media_12C = 0.0
    ds.F26BP_12C = 0.0
    ds.Citrate_12C = 0.0
    ds.Phenylalanine_12C = 0.0

    ds.Glucose_media_13C = 0.0
    ds.Glucose_13C = rate_isotope_tracing_GLUT(s, params, :C13) - rate_isotope_tracing_HK1(s, params, :C13)
    ds.G6P_13C = rate_isotope_tracing_HK1(s, params, :C13) - rate_isotope_tracing_GPI(s, params, :C13)
    ds.F6P_13C = (rate_isotope_tracing_GPI(s, params, :C13) - rate_isotope_tracing_PFKP(s, params, :C13))
    ds.F16BP_13C = (rate_isotope_tracing_PFKP(s, params, :C13) - rate_isotope_tracing_ALDO(s, params, :C13))
    ds.GAP_13C = (
        rate_isotope_tracing_ALDO(s, params, :C13) + rate_isotope_tracing_TPI(s, params, :C13) -
        rate_isotope_tracing_GAPDH(s, params, :C13)
    )
    ds.DHAP_13C = rate_isotope_tracing_ALDO(s, params, :C13) - rate_isotope_tracing_TPI(s, params, :C13)
    ds.BPG_13C = rate_isotope_tracing_GAPDH(s, params, :C13) - rate_isotope_tracing_PGK(s, params, :C13)
    ds.ThreePG_13C = rate_isotope_tracing_PGK(s, params, :C13) - rate_isotope_tracing_PGM(s, params, :C13)
    ds.TwoPG_13C = rate_isotope_tracing_PGM(s, params, :C13) - rate_isotope_tracing_ENO(s, params, :C13)
    ds.PEP_13C = rate_isotope_tracing_ENO(s, params, :C13) - rate_isotope_tracing_PKM2(s, params, :C13)
    ds.Pyruvate_13C = rate_isotope_tracing_PKM2(s, params, :C13) - rate_isotope_tracing_LDH(s, params, :C13)
    ds.Lactate_13C = rate_isotope_tracing_LDH(s, params, :C13) - rate_isotope_tracing_MCT(s, params, :C13)
    ds.Lactate_media_13C = 0.0
    ds.F26BP_13C = 0.0
    ds.Citrate_13C = 0.0
    ds.Phenylalanine_13C = 0.0

    ds.ATP = (
        -rate_HK1(total_conc_larray, params) - rate_PFKP(total_conc_larray, params) +
        rate_PGK(total_conc_larray, params) +
        rate_PKM2(total_conc_larray, params) - rate_ATPase(total_conc_larray, params) +
        rate_AK(total_conc_larray, params)
    )
    ds.ADP = (
        rate_HK1(total_conc_larray, params) + rate_PFKP(total_conc_larray, params) -
        rate_PGK(total_conc_larray, params) - rate_PKM2(total_conc_larray, params) +
        rate_ATPase(total_conc_larray, params) - 2 * rate_AK(total_conc_larray, params)
    )
    ds.AMP = rate_AK(total_conc_larray, params)
    ds.Phosphate = rate_ATPase(total_conc_larray, params) - rate_GAPDH(total_conc_larray, params)
    ds.NAD = rate_LDH(total_conc_larray, params) - rate_GAPDH(total_conc_larray, params)
    ds.NADH = rate_GAPDH(total_conc_larray, params) - rate_LDH(total_conc_larray, params)
end
