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
        rate_AK(s.ATP, s.ADP, s.AMP, params)
    )
    r.ATPase = rate_ATPase(s.ATP, s.ADP, s.Phosphate, params)
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
    b.Glucose = (
        binding_GLUT(f.Glucose_media, f.Glucose, params).Glucose +
        binding_HK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).Glucose
    )
    b.G6P = (
        binding_HK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).G6P +
        binding_GPI(f.G6P, f.F6P, params).G6P
    )
    b.F6P = (
        binding_GPI(f.G6P, f.F6P, params).F6P +
        binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).F6P
    )
    b.F16BP = (
        binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).F16BP +
        binding_ALDO(f.F16BP, f.GAP, f.DHAP, params).F16BP +
        binding_PKM2(f.PEP, f.ADP, f.Pyruvate, f.ATP, f.F16BP, f.Phenylalanine, params).F16BP
    )
    b.GAP = (
        binding_ALDO(f.F16BP, f.GAP, f.DHAP, params).GAP +
        binding_TPI(f.GAP, f.DHAP, params).GAP +
        binding_GAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).GAP
    )
    b.DHAP = (binding_ALDO(f.F16BP, f.GAP, f.DHAP, params).DHAP + binding_TPI(f.GAP, f.DHAP, params).DHAP)
    b.BPG = (
        binding_GAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).BPG +
        binding_PGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).BPG
    )
    b.ThreePG = (
        binding_PGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).ThreePG +
        binding_PGM(f.ThreePG, f.TwoPG, params).ThreePG
    )
    b.TwoPG = (binding_PGM(f.ThreePG, f.TwoPG, params).TwoPG + binding_ENO(f.TwoPG, f.PEP, params).TwoPG)
    b.PEP = (
        binding_ENO(f.TwoPG, f.PEP, params).PEP +
        binding_PKM2(f.PEP, f.ADP, f.Pyruvate, f.ATP, f.F16BP, f.Phenylalanine, params).PEP
    )
    b.Pyruvate = binding_PKM2(f.PEP, f.ADP, f.Pyruvate, f.ATP, f.F16BP, f.Phenylalanine, params).Pyruvate
    b.Lactate = binding_MCT(f.Lactate, f.Lactate_media, params).Lactate
    b.Lactate_media = 0.0
    b.ATP = (
        binding_HK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).ATP +
        binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).ATP +
        binding_PGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).ATP +
        binding_PKM2(f.PEP, f.ADP, f.Pyruvate, f.ATP, f.F16BP, f.Phenylalanine, params).ATP
    )
    b.ADP = (
        binding_HK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).ADP +
        binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).ADP +
        binding_PGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).ADP +
        binding_PKM2(f.PEP, f.ADP, f.Pyruvate, f.ATP, f.F16BP, f.Phenylalanine, params).ADP
    )
    b.AMP = 0
    b.Phosphate = (
        binding_GAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).Phosphate +
        binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).Phosphate +
        binding_HK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).Phosphate
    )
    b.NAD = (
        binding_LDH(f.Pyruvate, f.NADH, f.NAD, f.Lactate, params).NAD +
        binding_GAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).NAD
    )
    b.NADH = (
        binding_GAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).NADH +
        binding_LDH(f.Pyruvate, f.NADH, f.NAD, f.Lactate, params).NADH
    )
    b.F26BP = (binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).F26BP)
    b.Citrate = (binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).Citrate)
    b.Phenylalanine =
        (binding_PKM2(f.PEP, f.ADP, f.Pyruvate, f.ATP, f.F16BP, f.Phenylalanine, params).Phenylalanine)
    total = (b + f) .* cell_volume_correction
    total.Glucose_media /= cell_volume_correction
    total.Lactate_media /= cell_volume_correction
    return total
end

function free_to_bound_conc(f, params)
    b = similar(f)
    b.Glucose_media = 0.0
    b.Glucose = (
        binding_GLUT(f.Glucose_media, f.Glucose, params).Glucose +
        binding_HK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).Glucose
    )
    b.G6P = (
        binding_HK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).G6P +
        binding_GPI(f.G6P, f.F6P, params).G6P
    )
    b.F6P = (
        binding_GPI(f.G6P, f.F6P, params).F6P +
        binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).F6P
    )
    b.F16BP = (
        binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).F16BP +
        binding_ALDO(f.F16BP, f.GAP, f.DHAP, params).F16BP +
        binding_PKM2(f.PEP, f.ADP, f.Pyruvate, f.ATP, f.F16BP, f.Phenylalanine, params).F16BP
    )
    b.GAP = (
        binding_ALDO(f.F16BP, f.GAP, f.DHAP, params).GAP +
        binding_TPI(f.GAP, f.DHAP, params).GAP +
        binding_GAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).GAP
    )
    b.DHAP = (ALDO(f.F16BP, f.GAP, f.DHAP, params).DHAP + TPI(f.GAP, f.DHAP, params).DHAP)
    b.BPG = (
        binding_GAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).BPG +
        binding_PGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).BPG
    )
    b.ThreePG = (
        binding_PGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).ThreePG +
        binding_PGM(f.ThreePG, f.TwoPG, params).ThreePG
    )
    b.TwoPG = (binding_PGM(f.ThreePG, f.TwoPG, params).TwoPG + binding_ENO(f.TwoPG, f.PEP, params).TwoPG)
    b.PEP = (
        binding_ENO(f.TwoPG, f.PEP, params).PEP +
        binding_PKM2(f.PEP, f.ADP, f.Pyruvate, f.ATP, f.F16BP, f.Phenylalanine, params).PEP
    )
    b.Pyruvate = binding_PKM2(f.PEP, f.ADP, f.Pyruvate, f.ATP, f.F16BP, f.Phenylalanine, params).Pyruvate
    b.Lactate = binding_MCT(f.Lactate, f.Lactate_media, params).Lactate
    b.Lactate_media = 0.0
    b.ATP = (
        binding_HK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).ATP +
        binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).ATP +
        binding_PGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).ATP +
        binding_PKM2(f.PEP, f.ADP, f.Pyruvate, f.ATP, f.F16BP, f.Phenylalanine, params).ATP
    )
    b.ADP = (
        binding_HK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).ADP +
        binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).ADP +
        binding_PGK(f.BPG, f.ADP, f.ATP, f.ThreePG, params).ADP +
        binding_PKM2(f.PEP, f.ADP, f.Pyruvate, f.ATP, f.F16BP, f.Phenylalanine, params).ADP
    )
    b.AMP = 0.0
    b.Phosphate = (
        binding_GAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).Phosphate +
        binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).Phosphate +
        binding_HK1(f.Glucose, f.G6P, f.ATP, f.ADP, f.Phosphate, params).Phosphate
    )
    b.NAD = (
        binding_LDH(f.Pyruvate, f.NADH, f.NAD, f.Lactate, params).NAD +
        binding_GAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).NAD
    )
    b.NADH = (
        binding_GAPDH(f.GAP, f.NAD, f.Phosphate, f.BPG, f.NADH, params).NADH +
        binding_LDH(f.Pyruvate, f.NADH, f.NAD, f.Lactate, params).NADH
    )
    b.F26BP = (binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).F26BP)
    b.Citrate = (binding_PFKP(f.F6P, f.ATP, f.F16BP, f.ADP, f.Phosphate, f.Citrate, f.F26BP, params).Citrate)
    b.Phenylalanine =
        (binding_PKM2(f.PEP, f.ADP, f.Pyruvate, f.ATP, f.F16BP, f.Phenylalanine, params).Phenylalanine)
    return b
end