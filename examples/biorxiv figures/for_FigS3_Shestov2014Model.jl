using LabelledArrays, Measurements

# Define cellular parameters and kinetic+thermodynamic constants of enzymes
Vmax_normalizer = 1e-3

water_fraction_cell_volume = 0.66
cytosol_fraction_cell_volume = 0.66
cell_volume_correction = water_fraction_cell_volume * cytosol_fraction_cell_volume
cell_protein_density = 0.2

shestov_glycolysis_params_w_uncertainty = LVector(
    GLUT_Km_glc = 2.1e-3,
    GLUT_Vmax_tr = Vmax_normalizer * 100.0,
    GLUT_Keq = 1.0 ± 0.0,
    GLUT_MW = 54084.0 / 1000,
    HK1_L = 1.0,
    HK1_Km_f = 7.5e-3,
    HK1_Km_r = 2.94e-13,
    HK1_Ki = 0.2e-3,
    HK1_Vmax_f = Vmax_normalizer * 176.0,
    HK1_Keq = exp(-(-19.22)*1000/(8.314*310.15)),
    # HK1_Keq = 2700.0 ± 800.0,
    HK1_MW = 102486.0 / 1000,
    GPI_Km_f = 0.25e-3,
    GPI_Km_r = 7.73e-5,
    GPI_Vmax_f = Vmax_normalizer * 858.0,
    GPI_Keq = exp(-(2.78)*1000/(8.314*310.15)),
    # GPI_Keq = 0.36 ± 0.11,
    GPI_MW = 63147.0 / 1000,
    PFKP_Km_f = 0.23e-3,
    PFKP_Km_r = 1.82e-13,
    PFKP_Ka = 1e-6,
    PFKP_Vmax_f = Vmax_normalizer * 321.0,
    PFKP_Keq = exp(-(-15.62)*1000/(8.314*310.15)),
    # PFKP_Keq = 760.0 ± 380.0,
    PFKP_MW = 85596.0 / 1000,
    ALDO_Km_f = 0.16e-3,
    ALDO_Km_r = 8.29e-5,
    ALDO_Vmax_f = Vmax_normalizer * 321.0,
    ALDO_Keq = exp(-(24.64)*1000/(8.314*310.15)),
    # ALDO_Keq = 2.7e-6 ± 1.2e-6,
    ALDO_MW = 39420.0 / 1000,
    TPI_Km_f = 0.04e-3,
    TPI_Km_r = 0.002e-3,
    TPI_Vmax_f = Vmax_normalizer * 859.0,
    TPI_Keq = exp(-(7.57)*1000/(8.314*310.15)),
    # TPI_Keq = 0.0045 ± 0.0024,
    TPI_MW = 26669.0 / 1000,
    GAPDH_Km_f = 4.392e-6,
    GAPDH_Km_r = 1e-14,
    GAPDH_Km_nad = 0.55e-3,
    GAPDH_Km_nadh = 0.001e-3,
    GAPDH_Km_pi = 4e-3,
    GAPDH_Km_bpg = 1e-4,
    GAPDH_Km_gap = 2e-6,
    GAPDH_Vmax_f = Vmax_normalizer * 781.0,
    GAPDH_Keq = exp(-(2.6)*1000/(8.314*310.15)),
    # GAPDH_Keq = 16 ± 5,
    GAPDH_MW = 36053.0 / 1000,
    PGK_Km_f = 1.18e-6,
    PGK_Km_r = 1.5e-3,
    PGK_Vmax_f = Vmax_normalizer * 221.0,
    PGK_Keq = exp(-(-21.6)*1000/(8.314*310.15)),
    # PGK_Keq = 2000 ± 700,
    PGK_MW = 44615.0 / 1000,
    PGM_Km_f = 0.5e-3,
    PGM_Km_r = 0.03e-3,
    PGM_Vmax_f = Vmax_normalizer * 527.895,
    PGM_Keq = exp(-(6.35)*1000/(8.314*310.15)),
    # PGM_Keq = 0.18 ± 0.05,
    PGM_MW = 28804.0 / 1000,
    ENO_Km_f = 0.03e-3,
    ENO_Km_r = 0.15e-3,
    ENO_Vmax_f = Vmax_normalizer * 1340.089,
    ENO_Keq = exp(-(-4.47)*1000/(8.314*310.15)),
    # ENO_Keq = 4.4 ± 1.0,
    ENO_MW = 47169.0 / 1000,
    PKM2_Vmax_f = Vmax_normalizer * 211.525,
    PKM2_Km_f = 1.77e-6,
    PKM2_Km_r = 0.13e-3,
    PKM2_Ka = 0.5e-3,
    PKM2_Keq = exp(-(-27.18)*1000/(8.314*310.15)),
    # PKM2_Keq = 20000.0 ± 6000.0,
    PKM2_MW = 57937.0 / 1000,
    LDH_Km_Pyruvate = 0.5e-3,
    LDH_Km_NADH = 0.001e-3,
    LDH_Km_Lactate = 5e-3,
    LDH_Km_NAD = 0.549e-3,
    LDH_Km_f = 5e-7,
    LDH_Km_r = 2.745e-3,
    LDH_Vmax_f = Vmax_normalizer * 434.0,
    LDH_Keq = exp(-(-23.9)*1000/(8.314*310.15)),
    # LDH_Keq = 13000 ± 5000,
    LDH_MW = 36689.0 / 1000,
    MCT_Km_Lactate = 3e-3,
    MCT_Vmax_ltr = Vmax_normalizer * 60.0,
    MCT_Keq = 1.0,
    MCT_MW = 53944.0 / 1000,
    AK_Km_f = 5e-3,
    AK_Km_r = 2e-3,
    AK_Vmax_f = Vmax_normalizer * 2000.0,
    AK_Keq = exp(-(0.0)*1000/(8.314*310.15)),
    # AK_Keq = 0.48,
    AK_MW = 21635.0 / 1000,
    ATPase_Km_f = 3e-3,
    ATPase_Km_r = 4.71e-12,
    ATPase_Keq = exp(-(-32.42)*1000/(8.314*310.15)),
    # ATPase_Keq = 83000.0,
    ATPase_Vmax_f = Vmax_normalizer * 390.0,
)
shestov_glycolysis_params = Measurements.value.(shestov_glycolysis_params_w_uncertainty)

shestov_glycolysis_init_conc = LVector(
    Glucose_media = 5e-3,
    Glucose = 2.5e-3,
    G6P = 0.25e-3,
    F6P = 7.73e-5,
    F16BP = 1.55e-4,
    GAP = 2e-6,
    DHAP = 4.14e-5,
    BPG = 1e-4,
    ThreePG = 0.5e-3,
    TwoPG = 3e-5,
    PEP = 0.15e-3,
    Pyruvate = 5e-4,
    Lactate = 5e-3,
    Lactate_media = 5e-3,
    ATP = 3e-3,
    ADP = 1.18e-5,
    AMP = 4.62e-8,
    Phosphate = 4e-3,
    NAD = 5.49e-4,
    NADH = 1e-6,
    F26BP = 0.0,
    Citrate = 0.0,
    Phenylalanine = 0.0,
)

# Warning for when all params are not Float64 (i.e. 1.0 instead 1). ODE solution speed slows down 10x.

if eltype(shestov_glycolysis_params) != Float64
    @warn "One or more of the kinetic or thermodynamic constants in params are not Float64 (e.g., 1 instead of 1.0). This will slow down Glycolysis Model solution by 10x but it will still work. Fix this to get better performance."
end

#Rate equations for glycolytic enzymes
function rate_GLUT(Glucose_media, Glucose, params)
    Rate = (
        (params.GLUT_Vmax_tr) * (Glucose_media / params.GLUT_Km_glc - Glucose / params.GLUT_Km_glc) /
        (1 + Glucose_media / params.GLUT_Km_glc + Glucose / params.GLUT_Km_glc)
    )
    return Rate
end

function rate_HK1(Glucose, G6P, ATP, ADP, Phosphate, params)
    Rate =
        (
            ((params.HK1_Vmax_f / params.HK1_Km_f) * (Glucose * ATP - G6P * ADP / params.HK1_Keq)) /
            (1 + Glucose * ATP / params.HK1_Km_f + G6P * ADP / params.HK1_Km_r)
        ) * (params.HK1_Ki / (params.HK1_Ki + params.HK1_L * G6P))
    return Rate
end

function rate_GPI(G6P, F6P, params)
    Rate = (
        (params.GPI_Vmax_f / params.GPI_Km_f) * (G6P - (1 / params.GPI_Keq) * F6P) /
        (1 + G6P / params.GPI_Km_f + F6P / params.GPI_Km_r)
    )
    return Rate
end

function rate_PFKP(F6P, ATP, F16BP, ADP, Phosphate, Citrate, F26BP, params)
    H_plus = 1e-7
    Rate =
        (
            ((params.PFKP_Vmax_f / params.PFKP_Km_f) * (F6P * ATP - F16BP * ADP * H_plus / params.PFKP_Keq)) /
            (1 + F6P * ATP / params.PFKP_Km_f + F16BP * ADP * H_plus / params.PFKP_Km_r)
        ) * (ADP / (params.PFKP_Ka + ADP))
    return Rate
end

function rate_ALDO(F16BP, GAP, DHAP, params)
    Rate = (
        (params.ALDO_Vmax_f / params.ALDO_Km_f) * (
            (F16BP - (1 / params.ALDO_Keq) * (DHAP * GAP)) /
            (1 + GAP * DHAP / (params.ALDO_Km_r) + F16BP / params.ALDO_Km_f)
        )
    )
    return Rate
end

function rate_TPI(GAP, DHAP, params)
    Rate = (
        (params.TPI_Vmax_f / params.TPI_Km_f) * (DHAP - (1 / params.TPI_Keq) * GAP) /
        (1 + (DHAP / params.TPI_Km_f) + (GAP / params.TPI_Km_r))
    )
    return Rate
end

function rate_GAPDH(GAP, NAD, Phosphate, BPG, NADH, params)
    H_plus = 1e-7
    Rate = (
        (params.GAPDH_Vmax_f / (params.GAPDH_Km_f)) *
        (GAP * NAD * Phosphate - (1 / params.GAPDH_Keq) * BPG * NADH * H_plus) / (
            1 +
            GAP / params.GAPDH_Km_gap +
            Phosphate / params.GAPDH_Km_pi +
            GAP * NAD * Phosphate / params.GAPDH_Km_f +
            BPG / params.GAPDH_Km_bpg +
            NADH / params.GAPDH_Km_nadh +
            BPG * NADH * H_plus / (params.GAPDH_Km_r)
        )
    )
    return Rate
end

function rate_PGK(BPG, ADP, ATP, ThreePG, params)
    Rate = (
        (params.PGK_Vmax_f / (params.PGK_Km_f)) * (BPG * ADP - (1 / params.PGK_Keq) * (ThreePG * ATP)) /
        (1 + BPG * ADP / (params.PGK_Km_f) + ThreePG * ATP / (params.PGK_Km_r))
    )
    return Rate
end

function rate_PGM(ThreePG, TwoPG, params)
    Rate = (
        (params.PGM_Vmax_f / params.PGM_Km_f) * (ThreePG - (1 / params.PGM_Keq) * TwoPG) /
        (1 + ThreePG / params.PGM_Km_f + TwoPG / params.PGM_Km_r)
    )
    return Rate
end

function rate_ENO(TwoPG, PEP, params)
    Rate = (
        (params.ENO_Vmax_f / params.ENO_Km_f) * (TwoPG - (1 / params.ENO_Keq) * PEP) /
        (1 + TwoPG / params.ENO_Km_f + PEP / params.ENO_Km_r)
    )
    return Rate
end

function rate_PKM2(PEP, ADP, Pyruvate, ATP, F16BP, Phenylalanine, params)

    Rate =
        (params.PKM2_Vmax_f / params.PKM2_Km_f) * (ADP * PEP - ATP * Pyruvate / params.PKM2_Keq) /
        (1 + ADP * PEP / params.PKM2_Km_f + ATP * Pyruvate / params.PKM2_Km_r) *
        (F16BP / (params.PKM2_Ka + F16BP))
    return Rate
end

function rate_LDH(Pyruvate, NADH, NAD, Lactate, params)
    Rate = (
        (params.LDH_Vmax_f / (params.LDH_Km_f)) * (Pyruvate * NADH - (1 / params.LDH_Keq) * (Lactate * NAD)) / (
            1 +
            Pyruvate / params.LDH_Km_Pyruvate +
            Lactate / params.LDH_Km_Lactate +
            NADH / params.LDH_Km_NADH +
            NAD / params.LDH_Km_NAD +
            Lactate * NAD / params.LDH_Km_r +
            Pyruvate * NADH / params.LDH_Km_f
        )
    )
    return Rate
end

function rate_MCT(Lactate, Lactate_media, params)
    Rate = (
        (params.MCT_Vmax_ltr / params.MCT_Km_Lactate) * (Lactate - (1 / params.MCT_Keq) * Lactate_media) /
        (1 + Lactate / params.MCT_Km_Lactate + Lactate_media / params.MCT_Km_Lactate)
    )
    return Rate
end

function rate_AK(ATP, ADP, AMP, params)
    Rate = (
        (params.AK_Vmax_f / (params.AK_Km_r)) * (ADP^2 - (1 / params.AK_Keq) * (ATP * AMP)) /
        (1 + ADP^2 / params.AK_Km_r + ATP * AMP / params.AK_Km_f)
    )
    return Rate
end

function rate_ATPase(ATP, ADP, Phosphate, params)
    H_plus = 1e-7
    Rate =
        (params.ATPase_Vmax_f / params.ATPase_Km_f) *
        (
            1 / (
                1 +
                ATP / params.ATPase_Km_f +
                ADP * Phosphate * H_plus / params.ATPase_Km_r
            )
        ) *
        (ATP - ADP * Phosphate * H_plus / params.ATPase_Keq)
    return Rate
end


#Glycolysis Model ODE system
function shestov_glycolysis_ODEs(ds, s, params, t)
    ds.Glucose_media = 0.0
    ds.Glucose =
        rate_GLUT(s.Glucose_media, s.Glucose, params) -
        rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params)
    ds.G6P = rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) - rate_GPI(s.G6P, s.F6P, params)
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
    ds.PEP =
        rate_ENO(s.TwoPG, s.PEP, params) -
        rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params)
    ds.Pyruvate =
        rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) -
        rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params)
    ds.Lactate =
        rate_LDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params) - rate_MCT(s.Lactate, s.Lactate_media, params)
    ds.Lactate_media = 0.0
    ds.ATP = (
        -rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) -
        rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) +
        rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) +
        rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) -
        rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) + rate_AK(s.ATP, s.ADP, s.AMP, params)
    )
    ds.ADP = (
        rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) +
        rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) -
        rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) -
        rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) +
        rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) - 2 * rate_AK(s.ATP, s.ADP, s.AMP, params)
    )
    ds.AMP = rate_AK(s.ATP, s.ADP, s.AMP, params)
    ds.Phosphate =
        rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) -
        rate_GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    #ds.Phosphate = 0
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

function shestov_conc_to_rates(s, params)
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

function shestov_conc_to_disequilibrium_ratios(s, params)
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