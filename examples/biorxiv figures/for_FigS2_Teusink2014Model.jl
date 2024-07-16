using LabelledArrays

# Define cellular parameters and kinetic+thermodynamic constants of enzymes
Vmax_normalizer = 20

params = LVector(
    GLUT_Km_Glucose_media = 1.19e-3,
    GLUT_Km_Glucose = 1.19e-3,
    GLUT_Ki = 0.91,
    GLUT_Vmax = (198 / 97.264) * 0.36 / Vmax_normalizer,
    GLUT_Keq = 1.0,
    GLUT_MW = 54084.0 / 1000,
    HK1_KmGlucose = 0.08e-3,
    HK1_KmATP = 0.15e-3,
    HK1_KmG6P = 30e-3,
    HK1_KmADP = 0.23e-3,
    HK1_KiG6P = 0.07e-3,
    HK1_Vmax = 0.84 / Vmax_normalizer,
    HK1_Keq = 3.8e3,
    HK1_MW = 102486.0 / 1000,
    GPI_Km_G6P = 1.4e-3,
    GPI_Km_F6P = 0.3e-3,
    GPI_Vmax = 1.26 / Vmax_normalizer,
    GPI_Keq = 0.314,
    GPI_MW = 50000.0 / 1000,
    PFKP_KmATP = 0.71e-3,
    PFKP_cATP = 3.0,
    PFKP_KmF6P = 0.1e-3,
    PFKP_KdATP = 0.65e-3,
    PFKP_KdAMP = 0.0995e-3,
    PFKP_KdF26BP = 6.82e-7,
    PFKP_KdF16BP = 0.111e-3,
    PFKP_Ci_ATP = 100.0,
    PFKP_Ci_AMP = 0.0845,
    PFKP_Ci_F26BP = 0.0174,
    PFKP_Ci_F16BP = 0.397,
    PFKP_L0 = 0.66,
    PFKP_gR = 5.12,
    PFKP_gT = 1.0,
    PFKP_Vmax = 0.68 / Vmax_normalizer,
    PFKP_Keq = 8e2,
    PFKP_MW = 85596.0 / 1000,
    ALDO_K_F16BP = 0.3e-3,
    ALDO_K_DHAP = 2.4e-3,
    ALDO_K_GAP = 2e-3,
    ALDO_Ki_GAP = 10e-3,
    ALDO_Vmax = 1.19 / Vmax_normalizer,
    ALDO_Keq = 0.069,
    ALDO_MW = 39420.0 / 1000,
    TPI_Km_DHAP = 0.7e-3,
    TPI_Km_GAP = 0.37e-3,
    TPI_Vmax = 8.4 / Vmax_normalizer,
    TPI_Keq = 0.045,
    TPI_MW = 26669.0 / 1000,
    GAPDH_Km_GAP = 0.21e-3,
    GAPDH_Km_NAD = 0.09e-3,
    GAPDH_Km_BPG = 9.8e-6,
    GAPDH_Km_NADH = 0.06e-3,
    GAPDH_Km_Phosphate = 1e-3,
    GAPDH_Vmax_for = 4.4 / Vmax_normalizer,
    GAPDH_Vmax_rev = 24.3 / Vmax_normalizer,
    GAPDH_Keq = 0.0056,
    GAPDH_MW = 36053.0 / 1000,
    PGK_Km_BPG = 3e-6,
    PGK_Km_ADP = 0.2e-3,
    PGK_Km_ThreePG = 0.53e-3,
    PGK_Km_ATP = 0.3e-3,
    PGK_Vmax = 4.8 / Vmax_normalizer,
    PGK_Keq = 3.2e3,
    PGK_MW = 44615.0 / 1000,
    PGM_Km_ThreePG = 1.2e-3,
    PGM_Km_TwoPG = 0.1e-3,
    PGM_Vmax = 9.4 / Vmax_normalizer,
    PGM_Keq = 0.19,
    PGM_MW = 28804.0 / 1000,
    ENO_Km_TwoPG = 0.04e-3,
    ENO_Km_PEP = 0.5e-3,
    ENO_Vmax = 1.35 / Vmax_normalizer,
    ENO_Keq = 6.7,
    ENO_MW = 47169.0 / 1000,
    PKM2_L0 = 60_000.0,
    PKM2_KmPEP = 0.19e-3,
    PKM2_KmADP = 0.3e-3,
    PKM2_KmPyr = 0.14e-3,
    PKM2_K_F16BP = 0.2e-3,
    PKM2_K_ATP = 9.3e-3,
    PKM2_Vmax = 4.05 / Vmax_normalizer,
    PKM2_n = 4.0,
    PKM2_Keq = 6.5e3,
    PKM2_MW = 57937.0 / 1000,
    PDC_K_Pyruvate = 4.33e-3,
    PDC_Vmax = (1062.58 / 174.19) * 0.65 / Vmax_normalizer,
    PDC_n = 1.9,
    ADH_K_Acetald = 1.11e-3,
    ADH_K_NADH = 0.11e-3,
    ADH_K_EtOH = 17e-3,
    ADH_K_NAD = 0.17e-3,
    ADH_Ki_Acetald = 1.1e-3,
    ADH_Ki_NADH = 0.031e-3,
    ADH_Ki_EtOH = 90e-3,
    ADH_Ki_NAD = 0.92e-3,
    ADH_Vmax_for = 3.0 / Vmax_normalizer,
    ADH_Vmax_rev = 3.0 / Vmax_normalizer,
    ADH_Keq = 2000.0,
    ADH_MW = 36689.0 / 1000,
    k_Pi_exchange = 0.0,
    AK_Km_ADP = 0.12e-3,
    AK_Km_ATP = 0.11e-3,
    AK_Km_AMP = 0.09e-3,
    AK_Vmax = 100.0 / Vmax_normalizer,
    AK_Keq = 0.48,
    AK_MW = 21635.0 / 1000,
    NDPK_Km_ATP = 2e-3,
    NDPK_Km_ADP = 0.1e-3,
    NDPK_Km_NTP = 0.5e-3,
    NDPK_Km_NDP = 0.2e-3,
    NDPK_Vmax = 100.0 / Vmax_normalizer,
    NDPK_Keq = 2.16,
    NDPK_MW = 17149.0 / 1000,
    ATPase_Km_ATP = 1e-6,
    ATPase_Keq = 83000.0,
    ATPase_Vmax = 0.1,
);


# Warning for when all params are not Float64 (i.e. 1.0 instead 1). ODE solution speed slows down 10x.

if eltype(params) != Float64
    @warn "One or more of the kinetic or thermodynamic constants in params are not Float64 (e.g., 1 instead of 1.0). This will slow down Glycolysis Model solution by 10x but it will still work. Fix this to get better performance."
end

#Rate equations for glycolytic enzymes
function GLUT(Glucose_media, Glucose, params)
    Rate = (
        (params.GLUT_Vmax / params.GLUT_Km_Glucose_media) * (Glucose_media - Glucose) / (
            1 +
            Glucose_media / params.GLUT_Km_Glucose_media +
            Glucose / params.GLUT_Km_Glucose +
            params.GLUT_Ki *
            (Glucose_media / params.GLUT_Km_Glucose_media) *
            (Glucose / params.GLUT_Km_Glucose)
        )
    )
    return Rate
end

function HK1(Glucose, G6P, ATP, ADP, params)
    Rate = (
        params.HK1_Vmax / (params.HK1_KmGlucose * params.HK1_KmATP) *
        (Glucose * ATP - G6P * ADP / params.HK1_Keq) / (
            (1 + ATP / params.HK1_KmATP + ADP / params.HK1_KmADP) *
            (1 + Glucose / params.HK1_KmGlucose + G6P / params.HK1_KmG6P + G6P / params.HK1_KiG6P)
        )
    )
    return Rate
end

function GPI(G6P, F6P, params)
    Rate = (
        (params.GPI_Vmax / params.GPI_Km_G6P) * (G6P - (1 / params.GPI_Keq) * F6P) /
        (1 + G6P / params.GPI_Km_G6P + F6P / params.GPI_Km_F6P)
    )
    return Rate
end

function PFKP(F6P, ATP, AMP, F26BP, F16BP, params)
    λ1 = F6P / params.PFKP_KmF6P
    λ2 = ATP / params.PFKP_KmATP
    R = 1 + λ1 * λ2 + params.PFKP_gR * λ1 * λ2
    T = 1 + params.PFKP_cATP * λ2
    L =
        params.PFKP_L0 *
        ((1 + params.PFKP_Ci_ATP * ATP / params.PFKP_KdATP) / (1 + ATP / params.PFKP_KdATP))^2 *
        ((1 + params.PFKP_Ci_AMP * AMP / params.PFKP_KdAMP) / (1 + AMP / params.PFKP_KdAMP))^2 *
        (
            (
                1 +
                params.PFKP_Ci_F26BP * F26BP / params.PFKP_KdF26BP +
                params.PFKP_Ci_F16BP * F16BP / params.PFKP_KdF16BP
            ) / (1 + F26BP / params.PFKP_KdF26BP + F16BP / params.PFKP_KdF16BP)
        )
    Rate = params.PFKP_Vmax * params.PFKP_gR * λ1 * λ2 * R / (R^2 + L * T^2)
    return Rate
end

function ALDO(F16BP, GAP, DHAP, params)
    Rate = (
        (params.ALDO_Vmax / params.ALDO_K_F16BP) * (
            (F16BP - (1 / params.ALDO_Keq) * (DHAP * GAP)) / (
                1 +
                F16BP / params.ALDO_K_F16BP +
                GAP / params.ALDO_K_GAP +
                DHAP / params.ALDO_K_DHAP +
                (GAP / params.ALDO_K_GAP) * (DHAP / params.ALDO_K_DHAP) +
                (GAP / params.ALDO_Ki_GAP) * (F16BP / params.ALDO_K_F16BP)
            )
        )
    )
    return Rate
end

function TPI(GAP, DHAP, params)
    Rate = (
        (params.TPI_Vmax / params.TPI_Km_DHAP) * (DHAP - (1 / params.TPI_Keq) * GAP) /
        (1 + (DHAP / params.TPI_Km_DHAP) + (GAP / params.TPI_Km_GAP))
    )
    return Rate
end

function GAPDH(GAP, NAD, Phosphate, BPG, NADH, params)
    Rate = (
        (
            params.GAPDH_Vmax_for * GAP * NAD * Phosphate /
            (params.GAPDH_Km_GAP * params.GAPDH_Km_NAD * params.GAPDH_Km_Phosphate) -
            params.GAPDH_Vmax_rev * BPG * NADH / (params.GAPDH_Km_BPG * params.GAPDH_Km_NADH)
        ) / (
            (1 + Phosphate / params.GAPDH_Km_Phosphate) *
            (1 + NAD / params.GAPDH_Km_NAD + NADH / params.GAPDH_Km_NADH) *
            (1 + GAP / params.GAPDH_Km_GAP + BPG / params.GAPDH_Km_BPG)
        )
    )
    return Rate
end

function PGK(BPG, ADP, ATP, ThreePG, params)
    Rate = (
        (params.PGK_Vmax / (params.PGK_Km_BPG * params.PGK_Km_ADP)) *
        (BPG * ADP - (1 / params.PGK_Keq) * (ThreePG * ATP)) / (
            (1 + BPG / params.PGK_Km_BPG + ThreePG / params.PGK_Km_ThreePG) *
            (1 + ADP / params.PGK_Km_ADP + ATP / params.PGK_Km_ATP)
        )
    )
    return Rate
end

function PGM(ThreePG, TwoPG, params)
    Rate = (
        (params.PGM_Vmax / params.PGM_Km_ThreePG) * (ThreePG - (1 / params.PGM_Keq) * TwoPG) /
        (1 + ThreePG / params.PGM_Km_ThreePG + TwoPG / params.PGM_Km_TwoPG)
    )
    return Rate
end

function ENO(TwoPG, PEP, params)
    Rate = (
        (params.ENO_Vmax / params.ENO_Km_TwoPG) * (TwoPG - (1 / params.ENO_Keq) * PEP) /
        (1 + TwoPG / params.ENO_Km_TwoPG + PEP / params.ENO_Km_PEP)
    )
    return Rate
end

function PKM2(PEP, ADP, F16BP, ATP, params)
    Rate = (
        params.PKM2_Vmax *
        (PEP / params.PKM2_KmPEP) *
        (ADP / ADP + params.PKM2_KmADP) *
        (1 + PEP / params.PKM2_KmPEP)^(params.PKM2_n - 1) / ((
            params.PKM2_L0 *
            ((1 + ATP / params.PKM2_K_ATP) / (1 + F16BP / params.PKM2_K_F16BP))^params.PKM2_n +
            (1 + PEP / params.PKM2_KmPEP)^params.PKM2_n
        ))
    )
    return Rate
end

function PDC(Pyruvate, params)
    Rate = (
        params.PDC_Vmax * (Pyruvate / params.PDC_K_Pyruvate)^params.PDC_n /
        (1 + (Pyruvate / params.PDC_K_Pyruvate)^params.PDC_n)
    )
    return Rate
end

function ADH(Acetald, NADH, NAD, EtOH, params)
    Rate = (
        params.ADH_Vmax_for * Acetald * NADH / (params.ADH_Ki_NADH * params.ADH_K_Acetald) -
        params.ADH_Vmax_rev * EtOH * NAD / (params.ADH_Ki_EtOH * params.ADH_K_NAD) / (
            1 +
            EtOH / params.ADH_Ki_EtOH +
            params.ADH_K_EtOH * NAD / (params.ADH_K_NAD * params.ADH_Ki_EtOH) +
            params.ADH_K_NADH * Acetald / (params.ADH_K_Acetald * params.ADH_Ki_NADH) +
            NADH / params.ADH_Ki_NADH +
            EtOH * NAD / (params.ADH_Ki_EtOH * params.ADH_K_NAD) +
            params.ADH_K_NADH * EtOH * Acetald /
            (params.ADH_Ki_EtOH * params.ADH_K_Acetald * params.ADH_Ki_NADH) +
            params.ADH_K_EtOH * NAD * NADH / (params.ADH_Ki_EtOH * params.ADH_K_NAD * params.ADH_Ki_NADH) +
            Acetald * NADH / (params.ADH_Ki_NADH * params.ADH_K_Acetald) +
            EtOH * NAD * Acetald / (params.ADH_Ki_EtOH * params.ADH_K_NAD * params.ADH_Ki_Acetald) +
            NAD * Acetald * NADH / (params.ADH_Ki_EtOH * params.ADH_K_NAD * params.ADH_Ki_NADH)
        )
    )
    return Rate
end

function NDPK(NTP, NDP, ATP, ADP, params)
    Rate = (
        (params.NDPK_Vmax / (params.NDPK_Km_ATP * params.NDPK_Km_NDP)) * (
            (ATP * NDP - (1 / params.NDPK_Keq) * (NTP * ADP)) / (
                1 +
                ATP / params.NDPK_Km_ADP +
                ADP / params.NDPK_Km_ADP +
                NTP / params.NDPK_Km_NTP +
                NDP / params.NDPK_Km_NDP
            )
        )
    )
    return Rate
end

function AK(ATP, ADP, AMP, params)
    Rate = (
        (params.AK_Vmax / (params.AK_Km_ADP^2)) * (ADP^2 - (1 / params.AK_Keq) * (ATP * AMP)) / (
            (1 + ADP / params.AK_Km_ADP + ATP / params.AK_Km_ATP) *
            (1 + ADP / params.AK_Km_ADP + AMP / params.AK_Km_AMP)
        )
    )
    return Rate
end

# Edit this equation to include ADP/Km_ADP and Pi/Km_Pi in the denominator to make it more accurate
function rate_ATPase(ATP, ADP, Phosphate, params)
    Rate =
        (params.ATPase_Vmax / params.ATPase_Km_ATP) *
        (1 / (1 + ATP / params.ATPase_Km_ATP)) *
        (ATP - ADP * Phosphate / params.ATPase_Keq)
    return Rate
end


#Glycolysis Model ODE system
function Glycolysis_ODEs(ds, s, params, t)
    ds.Glucose_media = 0
    ds.Glucose = GLUT(s.Glucose_media, s.Glucose, params) - HK1(s.Glucose, s.G6P, s.ATP, s.ADP, params)
    ds.G6P = HK1(s.Glucose, s.G6P, s.ATP, s.ADP, params) - GPI(s.G6P, s.F6P, params)
    ds.F6P = (GPI(s.G6P, s.F6P, params) - PFKP(s.F6P, s.ATP, s.AMP, s.F26BP, s.F16BP, params))
    ds.F16BP = (PFKP(s.F6P, s.ATP, s.AMP, s.F26BP, s.F16BP, params) - ALDO(s.F16BP, s.GAP, s.DHAP, params))
    ds.GAP = (
        ALDO(s.F16BP, s.GAP, s.DHAP, params) + TPI(s.GAP, s.DHAP, params) -
        GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    )
    ds.DHAP = ALDO(s.F16BP, s.GAP, s.DHAP, params) - TPI(s.GAP, s.DHAP, params)
    ds.BPG =
        GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params) - PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params)
    ds.ThreePG = PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) - PGM(s.ThreePG, s.TwoPG, params)
    ds.TwoPG = PGM(s.ThreePG, s.TwoPG, params) - ENO(s.TwoPG, s.PEP, params)
    ds.PEP = ENO(s.TwoPG, s.PEP, params) - PKM2(s.PEP, s.ADP, s.F16BP, s.ATP, params)
    ds.Pyruvate = PKM2(s.PEP, s.ADP, s.F16BP, s.ATP, params) - PDC(s.Pyruvate, params)
    ds.Acetald = PDC(s.Pyruvate, params) - ADH(s.Acetald, s.NADH, s.NAD, s.EtOH, params)
    ds.EtOH = 0
    ds.ATP = (
        -HK1(s.Glucose, s.G6P, s.ATP, s.ADP, params) - PFKP(s.F6P, s.ATP, s.AMP, s.F26BP, s.F16BP, params) +
        PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) +
        PKM2(s.PEP, s.ADP, s.F16BP, s.ATP, params) - rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) +
        AK(s.ATP, s.ADP, s.AMP, params) - NDPK(s.NTP, s.NDP, s.ATP, s.ADP, params)
    )
    ds.ADP = (
        HK1(s.Glucose, s.G6P, s.ATP, s.ADP, params) + PFKP(s.F6P, s.ATP, s.AMP, s.F26BP, s.F16BP, params) -
        PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) - PKM2(s.PEP, s.ADP, s.F16BP, s.ATP, params) +
        rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) - 2 * AK(s.ATP, s.ADP, s.AMP, params) +
        NDPK(s.NTP, s.NDP, s.ATP, s.ADP, params)
    )
    ds.AMP = AK(s.ATP, s.ADP, s.AMP, params)
    ds.NTP = NDPK(s.NTP, s.NDP, s.ATP, s.ADP, params)
    ds.NDP = -NDPK(s.NTP, s.NDP, s.ATP, s.ADP, params)
    ds.Phosphate =
        rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) - GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params) +
        params.k_Pi_exchange * (10e-3 - s.Phosphate)
    ds.NAD =
        ADH(s.Acetald, s.NADH, s.NAD, s.EtOH, params) -
        GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    ds.NADH =
        GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params) -
        ADH(s.Acetald, s.NADH, s.NAD, s.EtOH, params)
    ds.F26BP = 0
end

function Conc_to_Rates(s, params)
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
        :ADH,
        :PDC,
        :ATPprod,
        :ATPase,
    )
    r.GLUT = GLUT(s.Glucose_media, s.Glucose, params)
    r.HK1 = HK1(s.Glucose, s.G6P, s.ATP, s.ADP, params)
    r.GPI = GPI(s.G6P, s.F6P, params)
    r.PFKP = PFKP(s.F6P, s.ATP, s.AMP, s.F26BP, s.F16BP, params)
    r.ALDO = ALDO(s.F16BP, s.GAP, s.DHAP, params)
    r.TPI = TPI(s.GAP, s.DHAP, params)
    r.GAPDH = GAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    r.PGK = PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params)
    r.PGM = PGM(s.ThreePG, s.TwoPG, params)
    r.ENO = ENO(s.TwoPG, s.PEP, params)
    r.PKM2 = PKM2(s.PEP, s.ADP, s.F16BP, s.ATP, params)
    r.PDC = PDC(s.Pyruvate, params)
    r.ADH = ADH(s.Acetald, s.NADH, s.NAD, s.EtOH, params)
    r.ATPprod = (
        -HK1(s.Glucose, s.G6P, s.ATP, s.ADP, params) - PFKP(s.F6P, s.ATP, s.AMP, s.F26BP, s.F16BP, params) +
        PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) +
        PKM2(s.PEP, s.ADP, s.F16BP, s.ATP, params) +
        AK(s.ATP, s.ADP, s.AMP, params) - NDPK(s.NTP, s.NDP, s.ATP, s.ADP, params)
    )
    r.ATPase = rate_ATPase(s.ATP, s.ADP, s.Phosphate, params)
    return r
end

function Conc_to_MassAction(s, params)
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
        :Q_Keq_ADH,
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
    k.Q_Keq_ADH = ((s.NAD * s.EtOH) / (s.Acetald * s.NADH)) / params.ADH_Keq
    k.Q_Keq_ATPase = ((s.Phosphate * s.ADP) / s.ATP) / params.ATPase_Keq
    return k
end;