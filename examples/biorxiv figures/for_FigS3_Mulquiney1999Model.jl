using LabelledArrays, Measurements

# Define cellular parameters and kinetic+thermodynamic constants of enzymes

# Adjust Vmax so that slowest enzyme has similar rate to our model
Vmax_normalizer = 3000

# water_fraction_cell_volume = 0.66
# cytosol_fraction_cell_volume = 0.66
# cell_volume_correction = water_fraction_cell_volume * cytosol_fraction_cell_volume
# cell_protein_density = 0.2
cell_volume_correction = 1.0
mulquiney_glycolysis_init_conc_w_uncertainty = LVector(
    Glucose_media = 5e-3,
    Glucose = (0.00304524 ± 0.000494698) / cell_volume_correction,
    G6P = (0.0390e-3) / cell_volume_correction,
    F6P = (0.0130e-3) / cell_volume_correction,
    F16BP = (0.00270e-3) / cell_volume_correction,
    GAP = (0.00570e-3) / cell_volume_correction,
    DHAP = (0.0170e-3) / cell_volume_correction,
    BPG = (0.000700e-3) / cell_volume_correction,
    ThreePG = (0.069e-3) / cell_volume_correction,
    TwoPG = (0.010e-3) / cell_volume_correction,
    PEP = (0.017e-3) / cell_volume_correction,
    Pyruvate = (0.085e-3) / cell_volume_correction,
    Lactate = (1.40e-3) / cell_volume_correction,
    Lactate_media = 1.83e-3,
    ATP = (2.1e-3) / cell_volume_correction,
    ADP = (0.31e-3) / cell_volume_correction,
    AMP = (0.03e-3) / cell_volume_correction,
    Phosphate = (1.00e-3) / cell_volume_correction,
    NAD = (0.0599e-3) / cell_volume_correction,
    NADH = (0.000245e-3) / cell_volume_correction,
    F26BP = 0.0 / cell_volume_correction,
    Citrate = 0.0 / cell_volume_correction,
    Phenylalanine = 0.0 / cell_volume_correction,
)
mulquiney_glycolysis_init_conc = Measurements.value.(mulquiney_glycolysis_init_conc_w_uncertainty)


mulquiney_glycolysis_params_w_uncertainty = LVector(
    GLUT_k1 = Vmax_normalizer * 1.0,
    GLUT_Keq = 1.0 ± 0.0,
    GLUT_MW = 54084.0 / 1000,
    HK1_Km_Glucose = 47e-6,
    HK1_Km_ATP = 1e-3,
    HK1_Km_G6P = 47e-6,
    HK1_Km_ADP = 1e-3,
    HK1_Ki_G6P = 10e-6,
    HK1_e0 = Vmax_normalizer * 24e-9,
    HK1_kcat_f = 180.0,
    HK1_kcat_r = 1.16,
    HK1_Keq = 155.0,
    HK1_MW = 102486.0 / 1000,
    GPI_Km_G6P = 181e-6,
    GPI_Km_F6P = 71e-6,
    GPI_e0 = Vmax_normalizer * 218e-9,
    GPI_kcat_f = 1760.0,
    GPI_kcat_r = 1470.0,
    GPI_Keq = 0.36 ± 0.11,
    GPI_MW = 63147.0 / 1000,
    PFKP_L = 1.0,
    PFKP_pKa = 7.05,
    PFKP_Km_R_F6P = 0.075e-3,
    PFKP_Km_R_ATP = 0.068e-3,
    PFKP_Km_R_F16BP = 0.50e-3,
    PFKP_Km_R_ADP = 0.54e-3,
    PFKP_K_R_Phosphate = 30e-3,
    PFKP_K_T_ATP = 100e-6,
    PFKP_K_T_Mg = 4e-3,
    PFKP_K_R_AMP = 0.3e-3,
    PFKP_e0 = Vmax_normalizer * 1.1e-7,
    PFKP_kcat_f = 822.0,
    PFKP_kcat_r = 36.0,
    PFKP_Keq = 760.0 ± 380.0,
    PFKP_MW = 85596.0 / 1000,
    ALDO_Km_F16BP = 7.1e-6,
    ALDO_Km_DHAP = 35e-6,
    ALDO_Ki_DHAP = 11e-6,
    ALDO_Km_GAP = 189e-6,
    ALDO_Ki_F16BP = 19.8e-6,
    ALDO_e0 = Vmax_normalizer * 37e-6,
    ALDO_kcat_f = 68.0,
    ALDO_kcat_r = 234.0,
    ALDO_Keq = 2.7e-6 ± 1.2e-6,
    ALDO_MW = 39420.0 / 1000,
    TPI_Km_DHAP = 162.4e-6 ± 80e-6,
    TPI_Km_GAP = 446e-6,
    TPI_e0 = Vmax_normalizer * 1.14e-6,
    TPI_Keq = 0.0045 ± 0.0024,
    TPI_kcat_f = 14560.0,
    TPI_kcat_r = 1280.0,
    TPI_MW = 26669.0 / 1000,
    GAPDH_Km_NAD = 45e-6,
    GAPDH_Ki_NAD = 45e-6,
    GAPDH_Km_Phosphate = 3.16e-3,
    GAPDH_Ki_Phosphate = 3.16e-3,
    GAPDH_Km_GAP = 95e-6,
    GAPDH_Ki_GAP = 1.59e-19,
    GAPDH_Ki_sharp_GAP = 0.031e-3,
    GAPDH_Km_NADH = 3.3e-6,
    GAPDH_Ki_NADH = 10e-6,
    GAPDH_Km_BPG = 0.671e-6,
    GAPDH_Ki_BPG = 1.52e-21,
    GAPDH_Ki_sharp_BPG = 1e-6,
    GAPDH_kcat_f = 232.0,
    GAPDH_kcat_r = 171.0,
    GAPDH_e0 = Vmax_normalizer * 7.66e-6,
    GAPDH_Keq = 16 ± 5,
    GAPDH_MW = 36053.0 / 1000,
    PGK_Km_BPG = 2e-6,
    PGK_Ki_BPG = 1.6e-3,
    PGK_Km_ADP = 100e-6,
    PGK_Ki_ADP = 80e-6,
    PGK_Km_ThreePG = 1.1e-3,
    PGK_Ki_ThreePG = 0.205e-3,
    PGK_Km_ATP = 1e-3,
    PGK_Ki_ATP = 0.186e-3,
    PGK_e0 = Vmax_normalizer * 2.74e-6,
    PGK_kcat_f = 2290.0,
    PGK_kcat_r = 917.0,
    PGK_Vmax = 3500.0 ± 900,
    PGK_Keq = 2000 ± 700,
    PGK_MW = 44615.0 / 1000,
    PGM_Km_ThreePG = 168e-6,
    PGM_Km_TwoPG = 26.6e-6,
    PGM_e0 = Vmax_normalizer * 410e-9,
    PGM_kcat_f = 795.0,
    PGM_kcat_r = 714.0,
    PGM_Keq = 0.18 ± 0.05,
    PGM_MW = 28804.0 / 1000,
    ENO_Km_TwoPG = 140e-6,
    ENO_Km_PEP = 110.5e-6,
    ENO_Km_Mg = 0.5e-3,
    ENO_e0 = Vmax_normalizer * 0.22e-6,
    ENO_kcat_f = 190.0,
    ENO_kcat_r = 50.0,
    ENO_Keq = 4.4 ± 1.0,
    ENO_MW = 47169.0 / 1000,
    PKM2_L = 1.0,
    PKM2_pKa = 6.8,
    PKM2_Km_R_PEP = 0.225e-3,
    PKM2_Km_R_ADP = 0.474e-3,
    PKM2_Km_R_Pyruvate = 2e-3,
    PKM2_Km_R_ATP = 3e-3,
    PKM2_K_T_ATP = 3.39e-3,
    PKM2_K_R_F16BP = 5e-6,
    PKM2_e0 = Vmax_normalizer * 87e-9,
    PKM2_kcat_f = 1386.0,
    PKM2_kcat_r = 3.26,
    PKM2_Keq = 20000.0 ± 6000.0,
    PKM2_MW = 57937.0 / 1000,
    LDH_Km_Pyruvate = 137e-6,
    LDH_Ki_Pyruvate = 228e-6,
    LDH_Ki_sharp_Pyruvate = 101e-6,
    LDH_Km_NADH = 8.44e-6,
    LDH_Ki_NADH = 2.45e-6,
    LDH_Km_Lactate = 1070e-6,
    LDH_Ki_Lactate = 7.33e-3,
    LDH_Km_NAD = 107e-6,
    LDH_Ki_NAD = 503e-6,
    LDH_e0 = Vmax_normalizer * 3.43e-6,
    LDH_kcat_f = 458.0,
    LDH_kcat_r = 115.0,
    LDH_Keq = 13000 ± 5000,
    LDH_MW = 36689.0 / 1000,
    MCT_k1 = Vmax_normalizer * 5.06e-3,
    MCT_Keq = 1.0,
    MCT_MW = 53944.0 / 1000,
    AK_Km_ADP = 0.0026e-3,
    AK_Km_ATP = 0.066e-3,
    AK_Km_AMP = 0.11e-3,
    AK_Km_MgADP = 0.11e-3,
    AK_e0 = Vmax_normalizer * 0.97e-6,
    AK_kcat_f = 2.08e3,
    AK_kcat_r = 3.8e3,
    AK_Keq = 0.48,
    AK_MW = 21635.0 / 1000,
    ATPase_Km_ATP = 1e-6,
    ATPase_Km_ADP = 1e-3,
    ATPase_Km_Phosphate = 1e-3,
    ATPase_Keq = 83000.0,
    ATPase_Vmax = 0.0002,
)
mulquiney_glycolysis_params = Measurements.value.(mulquiney_glycolysis_params_w_uncertainty)


# Warning for when all params are not Float64 (i.e. 1.0 instead 1). ODE solution speed slows down 10x.

if eltype(mulquiney_glycolysis_params) != Float64
    @warn "One or more of the kinetic or thermodynamic constants in params are not Float64 (e.g., 1 instead of 1.0). This will slow down Glycolysis Model solution by 10x but it will still work. Fix this to get better performance."
end

#Rate equations for glycolytic enzymes
function rate_GLUT(Glucose_media, Glucose, params)
    Rate = params.GLUT_k1 * (Glucose_media - Glucose / params.GLUT_Keq)
    return Rate
end

function rate_HK1(Glucose, G6P, ATP, ADP, Phosphate, params)
    pH = 7.2
    kcat_pH_adjustement = 1.0 / (1.0 + 10^(-pH) / 10^(-7.02) + 10^(-9.0) / 10^(-pH))
    Z = (
        1 +
        (Glucose / params.HK1_Km_Glucose) +
        (ATP / params.HK1_Km_ATP) +
        (G6P / params.HK1_Km_G6P) +
        (ADP / params.HK1_Km_ADP) +
        (Glucose / params.HK1_Km_Glucose) * (ATP / params.HK1_Km_ATP) +
        (G6P / params.HK1_Km_G6P) * (ADP / params.HK1_Km_ADP) +
        (Glucose / params.HK1_Km_Glucose) * (G6P / params.HK1_Ki_G6P)
    )
    Rate =
        (
            params.HK1_e0 * (
                kcat_pH_adjustement *
                params.HK1_kcat_f *
                (Glucose / params.HK1_Km_Glucose) *
                (ATP / params.HK1_Km_ATP) -
                kcat_pH_adjustement *
                params.HK1_kcat_r *
                (G6P / params.HK1_Km_G6P) *
                (ADP / params.HK1_Km_ADP)
            )
        ) / Z
    return Rate
end

function rate_GPI(G6P, F6P, params)
    Rate = (
        params.GPI_e0 *
        (params.GPI_kcat_f * G6P / params.GPI_Km_G6P - params.GPI_kcat_r * F6P / params.GPI_Km_F6P) /
        (1 + G6P / params.GPI_Km_G6P + F6P / params.GPI_Km_F6P)
    )
    return Rate
end

function rate_PFKP(F6P, ATP, F16BP, ADP, Phosphate, AMP, params)
    pH = 7.2
    Mg = 3.4e-3
    L_star =
        params.PFKP_L * (
            (10^(-pH) / 10^(-params.PFKP_pKa))^4 *
            (1 + ATP / params.PFKP_K_T_ATP)^4 *
            (1 + Mg / params.PFKP_K_T_Mg)^4
        ) / (
            (1 + F6P / params.PFKP_Km_R_F6P + F16BP / params.PFKP_Km_R_F16BP)^4 *
            (1 + AMP / params.PFKP_K_R_AMP)^4 *
            (1 + Phosphate / params.PFKP_K_R_Phosphate)^4
        )

    Z_R = (
        1 +
        (F6P / params.PFKP_Km_R_F6P) +
        (ATP / params.PFKP_Km_R_ATP) +
        (F16BP / params.PFKP_Km_R_F16BP) +
        (ADP / params.PFKP_Km_R_ADP) +
        (F6P / params.PFKP_Km_R_F6P) * (ATP / params.PFKP_Km_R_ATP) +
        (F16BP / params.PFKP_Km_R_F16BP) * (ADP / params.PFKP_Km_R_ADP)
    )

    Rate =
        (1 / (1 + L_star)) * (
            params.PFKP_e0 * (
                params.PFKP_kcat_f * (F6P / params.PFKP_Km_R_F6P) * (ATP / params.PFKP_Km_R_ATP) -
                params.PFKP_kcat_r * (F16BP / params.PFKP_Km_R_F16BP) * (ADP / params.PFKP_Km_R_ADP)
            )
        ) / Z_R

    return Rate
end

function rate_TPI(GAP, DHAP, params)
    Rate = (
        params.TPI_e0 * (
            params.TPI_kcat_f * (DHAP / params.TPI_Km_DHAP) - params.TPI_kcat_r * (GAP / params.TPI_Km_GAP)
        ) / (1 + (DHAP / params.TPI_Km_DHAP) + (GAP / params.TPI_Km_GAP))
    )
    return Rate
end

function rate_ALDO(F16BP, GAP, DHAP, params)
    Rate = (
        (params.ALDO_e0) * (
            (
                params.ALDO_kcat_f * F16BP / params.ALDO_Km_F16BP -
                params.ALDO_kcat_r * GAP * DHAP / (params.ALDO_Km_GAP * params.ALDO_Ki_DHAP)
            ) / (
                1 +
                F16BP / params.ALDO_Km_F16BP +
                GAP * params.ALDO_Km_DHAP / (params.ALDO_Km_GAP * params.ALDO_Ki_DHAP) +
                GAP * DHAP / (params.ALDO_Km_GAP * params.ALDO_Ki_DHAP) +
                DHAP / params.ALDO_Ki_DHAP +
                F16BP * GAP * params.ALDO_Km_DHAP /
                (params.ALDO_Ki_F16BP * params.ALDO_Km_GAP * params.ALDO_Ki_DHAP)
            )
        )
    )
    return Rate
end

function rate_GAPDH(GAP, NAD, Phosphate, BPG, NADH, params)
    pH = 7.2
    kcat_pH_adjustement = 1.0 / (1.0 + 10^(-pH) / 10^(-7.5) + 10^(-10.0) / 10^(-pH))
    Rate =
        params.GAPDH_e0 * (
            kcat_pH_adjustement *
            params.GAPDH_kcat_f *
            (NAD / params.GAPDH_Km_NAD) *
            (Phosphate / params.GAPDH_Ki_Phosphate) *
            (GAP / params.GAPDH_Ki_GAP) -
            kcat_pH_adjustement *
            params.GAPDH_kcat_r *
            (BPG / params.GAPDH_Ki_GAP) *
            (NADH / params.GAPDH_Km_NADH) *
            (10^-(pH))
        ) / (
            (
                GAP / params.GAPDH_Ki_GAP +
                BPG / params.GAPDH_Ki_BPG +
                Phosphate * GAP / (params.GAPDH_Ki_GAP * params.GAPDH_Ki_Phosphate)
            ) * (1 + GAP / params.GAPDH_Ki_sharp_GAP) +
            params.GAPDH_Km_GAP * NAD * Phosphate /
            (params.GAPDH_Km_NAD * params.GAPDH_Ki_Phosphate * params.GAPDH_Ki_GAP) +
            NAD * GAP / (params.GAPDH_Ki_GAP * params.GAPDH_Ki_NAD) +
            (NAD / params.GAPDH_Km_NAD) *
            (Phosphate / params.GAPDH_Ki_Phosphate) *
            (GAP / params.GAPDH_Ki_GAP) +
            (BPG / params.GAPDH_Ki_GAP) * (NADH / params.GAPDH_Km_NADH) * 10^(-pH) +
            NAD * BPG / (params.GAPDH_Ki_BPG * params.GAPDH_Ki_NAD) +
            params.GAPDH_Km_GAP * NAD * Phosphate * BPG / (
                params.GAPDH_Ki_GAP *
                params.GAPDH_Ki_Phosphate *
                params.GAPDH_Km_NAD *
                params.GAPDH_Ki_sharp_BPG
            ) +
            10^(-pH) * (
                params.GAPDH_Km_BPG * NADH / (params.GAPDH_Ki_BPG * params.GAPDH_Km_NADH) +
                params.GAPDH_Km_BPG * NADH * Phosphate /
                (params.GAPDH_Ki_BPG * params.GAPDH_Km_NADH * params.GAPDH_Ki_Phosphate) +
                GAP * NADH / (params.GAPDH_Ki_NADH * params.GAPDH_Ki_GAP) +
                BPG * NADH / (params.GAPDH_Ki_BPG * params.GAPDH_Ki_NADH) +
                Phosphate * GAP * NADH /
                (params.GAPDH_Ki_NADH * params.GAPDH_Ki_GAP * params.GAPDH_Ki_Phosphate) +
                params.GAPDH_Km_BPG * Phosphate * BPG * NADH / (
                    params.GAPDH_Km_NADH *
                    params.GAPDH_Ki_BPG *
                    params.GAPDH_Ki_Phosphate *
                    params.GAPDH_Ki_sharp_BPG
                )
            )
        )
    return Rate
end

function rate_PGK(BPG, ADP, ATP, ThreePG, params)
    Z = (
        1 +
        (BPG / params.PGK_Ki_BPG) +
        (ATP / params.PGK_Ki_ATP) +
        (ThreePG / params.PGK_Ki_ThreePG) +
        (ADP / params.PGK_Ki_ADP) +
        (BPG / params.PGK_Km_BPG) * (ADP / params.PGK_Km_ADP) +
        (ThreePG / params.PGK_Km_ThreePG) * (ATP / params.PGK_Ki_ATP)
    )
    Rate =
        (
            params.PGK_e0 * (
                params.PGK_kcat_f * (BPG / params.PGK_Km_BPG) * (ADP / params.PGK_Ki_ADP) -
                params.PGK_kcat_r * (ThreePG / params.PGK_Km_ThreePG) * (ATP / params.PGK_Ki_ATP)
            )
        ) / Z
    return Rate
end

function rate_PGM(ThreePG, TwoPG, params)
    Rate = (
        params.PGM_e0 * (
            params.PGM_kcat_f * ThreePG / params.PGM_Km_ThreePG -
            params.PGM_kcat_r * TwoPG / params.PGM_Km_TwoPG
        ) / (1 + ThreePG / params.PGM_Km_ThreePG + TwoPG / params.PGM_Km_TwoPG)
    )
    return Rate
end

function rate_ENO(TwoPG, PEP, params)
    Mg = 3.4e-3
    Rate = (
        params.ENO_e0 * (
            params.ENO_kcat_f * (TwoPG / params.ENO_Km_TwoPG) * (Mg / params.ENO_Km_Mg) -
            params.ENO_kcat_r * (PEP / params.ENO_Km_PEP) * (Mg / params.ENO_Km_Mg)
        ) / (
            1 +
            TwoPG / params.ENO_Km_TwoPG +
            PEP / params.ENO_Km_PEP +
            (TwoPG / params.ENO_Km_TwoPG) * (Mg / params.ENO_Km_Mg) +
            (PEP / params.ENO_Km_PEP) * (Mg / params.ENO_Km_Mg) +
            (Mg / params.ENO_Km_Mg)
        )
    )
    return Rate
end

function rate_PKM2(PEP, ADP, Pyruvate, ATP, F16BP, Phenylalanine, params)
    pH = 7.2
    L_star =
        params.PKM2_L * ((10^(-pH) / 10^(-params.PKM2_pKa))^4 * (1 + ATP / params.PKM2_K_T_ATP)^4) / (
            (1 + PEP / params.PKM2_Km_R_PEP + Pyruvate / params.PKM2_Km_R_Pyruvate)^4 *
            (1 + F16BP / params.PKM2_K_R_F16BP)^4
        )

    Z_R = (
        1 +
        (PEP / params.PKM2_Km_R_PEP) +
        (ATP / params.PKM2_Km_R_ATP) +
        (Pyruvate / params.PKM2_Km_R_Pyruvate) +
        (ADP / params.PKM2_Km_R_ADP) +
        (PEP / params.PKM2_Km_R_PEP) * (ADP / params.PKM2_Km_R_ATP) +
        (Pyruvate / params.PKM2_Km_R_Pyruvate) * (ATP / params.PKM2_Km_R_ADP)
    )

    Rate =
        (1 / (1 + L_star)) * (
            params.PKM2_e0 * (
                params.PKM2_kcat_f * (PEP / params.PKM2_Km_R_PEP) * (ADP / params.PKM2_Km_R_ADP) -
                params.PKM2_kcat_r * (Pyruvate / params.PKM2_Km_R_Pyruvate) * (ATP / params.PKM2_Km_R_ATP)
            )
        ) / Z_R

    return Rate
end

function rate_LDH(Pyruvate, NADH, NAD, Lactate, params)
    Rate = (
        params.LDH_e0 * (
            params.LDH_kcat_f * Pyruvate * NADH / (params.LDH_Ki_NADH * params.LDH_Km_Pyruvate) -
            params.LDH_kcat_r * (Lactate * NAD / (params.LDH_Ki_NAD * params.LDH_Km_Lactate))
        ) / (
            (
                1 +
                params.LDH_Km_NADH * Pyruvate / (params.LDH_Ki_NADH * params.LDH_Km_Pyruvate) +
                params.LDH_Km_NAD * Lactate / (params.LDH_Ki_NAD * params.LDH_Km_Lactate)
            ) * (1 + Pyruvate / params.LDH_Ki_sharp_Pyruvate) +
            NADH / (params.LDH_Ki_NADH) +
            NAD / (params.LDH_Ki_NAD) +
            Pyruvate * NAD / (params.LDH_Ki_NAD * params.LDH_Km_Pyruvate) +
            params.LDH_Km_NAD * NADH * Lactate /
            (params.LDH_Ki_NADH * params.LDH_Km_Lactate * params.LDH_Ki_NAD) +
            params.LDH_Km_NADH * Pyruvate * NAD /
            (params.LDH_Ki_NADH * params.LDH_Km_Pyruvate * params.LDH_Ki_NAD) +
            Lactate * NAD / (params.LDH_Ki_NAD * params.LDH_Km_Lactate) +
            NADH * Pyruvate * Lactate /
            (params.LDH_Ki_NADH * params.LDH_Km_Pyruvate * params.LDH_Ki_Lactate) +
            Pyruvate * Lactate * NAD / (params.LDH_Ki_Pyruvate * params.LDH_Km_Lactate * params.LDH_Ki_NAD)
        )
    )
    return Rate
end

function rate_MCT(Lactate, Lactate_media, params)
    Rate = params.MCT_k1 * (Lactate - Lactate_media / params.MCT_Keq)
    return Rate
end

function rate_AK(ATP, ADP, AMP, params)
    Mg = 3.4e-3
    K_Mg = 3568 # M^-1 average from https://doi.org/10.1021/bi00317a026
    MgADP = (1 / K_Mg + ADP + Mg - sqrt((1 / K_Mg + ADP + Mg)^2 - 4 * ADP * Mg)) / 2
    freeADP = ADP - MgADP
    Rate = (
        (params.AK_e0) * (
            params.AK_kcat_f * (freeADP / params.AK_Km_ADP) * (MgADP / params.AK_Km_MgADP) -
            params.AK_kcat_r * (ATP / params.AK_Km_ATP) * (AMP / params.AK_Km_AMP)
        ) / (
            1 +
            freeADP / params.AK_Km_ADP +
            MgADP / params.AK_Km_MgADP +
            ATP / params.AK_Km_ATP +
            AMP / params.AK_Km_AMP +
            (freeADP / params.AK_Km_ADP) * (MgADP / params.AK_Km_MgADP) +
            (ATP / params.AK_Km_ATP) * (AMP / params.AK_Km_AMP)
        )
    )
    return Rate
end

function rate_ATPase(ATP, ADP, Phosphate, params)
    Rate =
    # ((0.01 * params.HK1_e0 * params.HK1_kcat_f + params.ATPase_Vmax) / params.ATPase_Km_ATP) *
        (params.ATPase_Vmax / params.ATPase_Km_ATP) *
        (
            1 / (
                1 +
                ATP / params.ATPase_Km_ATP +
                ADP / params.ATPase_Km_ADP +
                Phosphate / params.ATPase_Km_Phosphate +
                (ADP / params.ATPase_Km_ADP) * (Phosphate / params.ATPase_Km_Phosphate)
            )
        ) *
        (ATP - ADP * Phosphate / params.ATPase_Keq)
    return Rate
end


#Glycolysis Model ODE system
function mulquiney_glycolysis_ODEs(ds, s, params, t)
    ds.Glucose_media = 0.0
    ds.Glucose =
        rate_GLUT(s.Glucose_media, s.Glucose, params) -
        rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params)
    ds.G6P = rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) - rate_GPI(s.G6P, s.F6P, params)
    ds.F6P =
        (rate_GPI(s.G6P, s.F6P, params) - rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.AMP, params))
    ds.F16BP = (
        rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.AMP, params) -
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
        rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.AMP, params) +
        rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) +
        rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) -
        rate_ATPase(s.ATP, s.ADP, s.Phosphate, params) + rate_AK(s.ATP, s.ADP, s.AMP, params)
    )
    ds.ADP = (
        rate_HK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) +
        rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.AMP, params) -
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

function mulquiney_conc_to_rates(s, params)
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
    r.PFKP = rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.AMP, params)
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
        rate_PFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.AMP, params) +
        rate_PGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) +
        rate_PKM2(s.PEP, s.ADP, s.Pyruvate, s.ATP, s.F16BP, s.Phenylalanine, params) +
        rate_AK(s.ATP, s.ADP, s.AMP, params)
    )
    r.ATPase = rate_ATPase(s.ATP, s.ADP, s.Phosphate, params)
    return r
end

function mulquiney_conc_to_disequilibrium_ratios(s, params)
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