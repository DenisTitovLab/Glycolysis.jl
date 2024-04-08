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
    r.GLUT = rate_GLUT(s, params)
    r.HK1 = rate_HK1(s, params)
    r.GPI = rate_GPI(s, params)
    r.PFKP = rate_PFKP(s, params)
    r.ALDO = rate_ALDO(s, params)
    r.TPI = rate_TPI(s, params)
    r.GAPDH = rate_GAPDH(s, params)
    r.PGK = rate_PGK(s, params)
    r.PGM = rate_PGM(s, params)
    r.ENO = rate_ENO(s, params)
    r.PKM2 = rate_PKM2(s, params)
    r.LDH = rate_LDH(s, params)
    r.MCT = rate_MCT(s, params)
    r.ATPprod = (
        -rate_HK1(s, params) - rate_PFKP(s, params) +
        rate_PGK(s, params) +
        rate_PKM2(s, params) +
        rate_AK(s, params) - rate_NDPK(s, params) - rate_CK(s, params)
    )
    r.ATPase = rate_ATPase(s, params)
    return r
end

function conc_to_disequilibrium_ratios(s, params)
    k = @LVector eltype(s) (
        :Q_Keq_GLUT,
        :Q_Keq_HK1,
        :Q_Keq_GPI,
        :Q_Keq_PFKP,
        :Q_Keq_ALDO,
        :Q_Keq_TPI,
        :Q_Keq_GAPDH,
        :Q_Keq_PGK,
        :Q_Keq_PGM,
        :Q_Keq_ENO,
        :Q_Keq_PKM2,
        :Q_Keq_LDH,
        :Q_Keq_MCT,
        :Q_Keq_ATPase,
    )
    k.Q_Keq_GLUT = (s.Glucose / s.Glucose_media) / params.GLUT_Keq
    k.Q_Keq_HK1 = ((s.G6P * s.ADP) / (s.Glucose * s.ATP)) / params.HK1_Keq
    k.Q_Keq_GPI = (s.F6P / s.G6P) / params.GPI_Keq
    k.Q_Keq_PFKP = ((s.F16BP * s.ADP) / (s.F6P * s.ATP)) / params.PFKP_Keq
    k.Q_Keq_ALDO = ((s.GAP * s.DHAP) / s.F16BP) / params.ALDO_Keq
    k.Q_Keq_TPI = (s.GAP / s.DHAP) / params.TPI_Keq
    k.Q_Keq_GAPDH = ((s.BPG * s.NADH) / (s.GAP * s.NAD * s.Phosphate)) / params.GAPDH_Keq
    k.Q_Keq_PGK = (s.ATP * s.ThreePG) / (s.BPG * s.ADP) / params.PGK_Keq
    k.Q_Keq_PGM = (s.TwoPG / s.ThreePG) / params.PGM_Keq
    k.Q_Keq_ENO = (s.PEP / s.TwoPG) / params.ENO_Keq
    k.Q_Keq_PKM2 = ((s.ATP * s.Pyruvate) / (s.PEP * s.ADP)) / params.PKM2_Keq
    k.Q_Keq_LDH = ((s.NAD * s.Lactate) / (s.Pyruvate * s.NADH)) / params.LDH_Keq
    k.Q_Keq_MCT = (s.Lactate_media / s.Lactate) / params.MCT_Keq
    k.Q_Keq_ATPase = ((s.Phosphate * s.ADP) / s.ATP) / params.ATPase_Keq
    return k
end

function free_to_total_conc(f, params)
    b = similar(f)
    b.Glucose_media = 0.0
    b.Glucose = (binding_GLUT(f, params).Glucose + binding_HK1(f, params).Glucose)
    b.G6P = (binding_HK1(f, params).G6P + binding_GPI(f, params).G6P)
    b.F6P = (binding_GPI(f, params).F6P + binding_PFKP(f, params).F6P)
    b.F16BP = (binding_PFKP(f, params).F16BP + binding_ALDO(f, params).F16BP + binding_PKM2(f, params).F16BP)
    b.GAP = (binding_ALDO(f, params).GAP + binding_TPI(f, params).GAP + binding_GAPDH(f, params).GAP)
    b.DHAP = (binding_ALDO(f, params).DHAP + binding_TPI(f, params).DHAP)
    b.BPG = (binding_GAPDH(f, params).BPG + binding_PGK(f, params).BPG)
    b.ThreePG = (binding_PGK(f, params).ThreePG + binding_PGM(f, params).ThreePG)
    b.TwoPG = (binding_PGM(f, params).TwoPG + binding_ENO(f, params).TwoPG)
    b.PEP = (binding_ENO(f, params).PEP + binding_PKM2(f, params).PEP)
    b.Pyruvate = binding_PKM2(f, params).Pyruvate
    b.Lactate = binding_MCT(f, params).Lactate
    b.Lactate_media = 0.0
    b.ATP = (
        binding_HK1(f, params).ATP +
        binding_PFKP(f, params).ATP +
        binding_PGK(f, params).ATP +
        binding_PKM2(f, params).ATP
    )
    b.ADP = (
        binding_HK1(f, params).ADP +
        binding_PFKP(f, params).ADP +
        binding_PGK(f, params).ADP +
        binding_PKM2(f, params).ADP
    )
    b.AMP = 0
    b.Phosphate = (
        binding_GAPDH(f, params).Phosphate +
        binding_PFKP(f, params).Phosphate +
        binding_HK1(f, params).Phosphate
    )
    b.NAD = (binding_LDH(f, params).NAD + binding_GAPDH(f, params).NAD)
    b.NADH = (binding_GAPDH(f, params).NADH + binding_LDH(f, params).NADH)
    b.F26BP = (binding_PFKP(f, params).F26BP)
    b.Citrate = (binding_PFKP(f, params).Citrate)
    b.Phenylalanine = (binding_PKM2(f, params).Phenylalanine)
    total = (b + f)
    total.Glucose_media
    total.Lactate_media
    return total
end

function free_to_bound_conc(f, params)
    b = similar(f)
    b.Glucose_media = 0.0
    b.Glucose = (binding_GLUT(f, params).Glucose + binding_HK1(f, params).Glucose)
    b.G6P = (binding_HK1(f, params).G6P + binding_GPI(f, params).G6P)
    b.F6P = (binding_GPI(f, params).F6P + binding_PFKP(f, params).F6P)
    b.F16BP = (binding_PFKP(f, params).F16BP + binding_ALDO(f, params).F16BP + binding_PKM2(f, params).F16BP)
    b.GAP = (binding_ALDO(f, params).GAP + binding_TPI(f, params).GAP + binding_GAPDH(f, params).GAP)
    b.DHAP = (ALDO(f, params).DHAP + TPI(f, params).DHAP)
    b.BPG = (binding_GAPDH(f, params).BPG + binding_PGK(f, params).BPG)
    b.ThreePG = (binding_PGK(f, params).ThreePG + binding_PGM(f, params).ThreePG)
    b.TwoPG = (binding_PGM(f, params).TwoPG + binding_ENO(f, params).TwoPG)
    b.PEP = (binding_ENO(f, params).PEP + binding_PKM2(f, params).PEP)
    b.Pyruvate = binding_PKM2(f, params).Pyruvate
    b.Lactate = binding_MCT(f, params).Lactate
    b.Lactate_media = 0.0
    b.ATP = (
        binding_HK1(f, params).ATP +
        binding_PFKP(f, params).ATP +
        binding_PGK(f, params).ATP +
        binding_PKM2(f, params).ATP
    )
    b.ADP = (
        binding_HK1(f, params).ADP +
        binding_PFKP(f, params).ADP +
        binding_PGK(f, params).ADP +
        binding_PKM2(f, params).ADP
    )
    b.AMP = 0.0
    b.Phosphate = (
        binding_GAPDH(f, params).Phosphate +
        binding_PFKP(f, params).Phosphate +
        binding_HK1(f, params).Phosphate
    )
    b.NAD = (binding_LDH(f, params).NAD + binding_GAPDH(f, params).NAD)
    b.NADH = (binding_GAPDH(f, params).NADH + binding_LDH(f, params).NADH)
    b.F26BP = (binding_PFKP(f, params).F26BP)
    b.Citrate = (binding_PFKP(f, params).Citrate)
    b.Phenylalanine = (binding_PKM2(f, params).Phenylalanine)
    return b
end
