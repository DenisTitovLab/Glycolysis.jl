using LabelledArrays, Measurements

#=
This file contains all the parameters needed to run Glycolysis.jl including
kinetic constants, thermodynamic constants, enzyme concentrations, initial
conditions for ODE simulations (i.e., estimates of intracellular concentration of metabolites),
and various constants for unit convertions.

Parameters are estimated from experimental data and many are in the form of Mean ± SEM to
allow for propagation of uncertainty.

=#

###########################################################################################
# Various constants for unit convertions
###########################################################################################
"Fraction of intracellular volume occupied by water with macromolecules making up the rest"
water_fraction_cell_volume = 0.66

"Fraction of intracellular volume occupied by cytosol"
cytosol_fraction_cell_volume = 0.66

"""
Correction for cytosolic volume and its water content to better estimate concentration of
cytosolic molecules from total cellular concentration and vice versa
"""
cell_volume_correction = water_fraction_cell_volume * cytosol_fraction_cell_volume

"Intracellular protein density, mg protein per µl of cell volume"
cell_protein_density = 0.2

###########################################################################################
# Initial conditions for ODE simulations
###########################################################################################

"""
Mean ± SEM of intracellular concentrations of glycolytic intermediates in mammalian cells (M units).
Values are corrected for protein concentration and cytosol fraction to better represent
metabolite concentration in cytosol.
"""
glycolysis_init_conc_w_uncertainty = LVector(
    Glucose_media = 25e-3,
    Glucose = (4.1e-3 ± 8e-4) / cell_volume_correction,
    G6P = (1.7e-4 ± 3e-5) / cell_volume_correction,
    F6P = (9.8e-5 ± 1.0e-5) / cell_volume_correction,
    F16BP = (5.5e-4 ± 1.0e-4) / cell_volume_correction,
    GAP = (9.6e-5 ± 2.0e-5) / cell_volume_correction,
    DHAP = (8.2e-4 ± 1.9e-4) / cell_volume_correction,
    BPG = (1.1e-6 ± 3e-7) / cell_volume_correction,
    ThreePG = (1.6e-3 ± 6e-4) / cell_volume_correction,
    TwoPG = (4.9e-4 ± 2e-4) / cell_volume_correction,
    PEP = (9.7e-5 ± 1.9e-5) / cell_volume_correction,
    Pyruvate = (6.1e-4 ± 2.3e-4) / cell_volume_correction,
    Lactate = (9.6e-3 ± 2.0e-3) / cell_volume_correction,
    Lactate_media = 2e-3,
    ATP = (3.6e-3 ± 5e-4) / cell_volume_correction,
    ADP = (7.3e-4 ± 1.5e-4) / cell_volume_correction,
    AMP = (1.6e-4 ± 6e-5) / cell_volume_correction,
    Phosphate = (4.0e-3 ± 9e-4) / cell_volume_correction,
    NTP = (0.00186 ± 0.00022) / cell_volume_correction,
    NDP = (0.000354 ± 0.000058) / cell_volume_correction,
    Phosphocreatine = 0.003 / cell_volume_correction,
    Creatine = 0.0003 / cell_volume_correction,
    NAD = (9.2e-4 ± 4.1e-4) / cell_volume_correction,
    NADH = (8.4e-5 ± 3.6e-5) / cell_volume_correction,
    F26BP = 0.0 / cell_volume_correction,
    Citrate = 0.0 / cell_volume_correction,
    Phenylalanine = 0.0 / cell_volume_correction,
)

"""
Intracellular concentrations of glycolytic intermediates in mammalian cells (M units).
`glycolysis_init_conc` can be used as initial condition for simulating `glycolysis_ODEs` glycolysis ODEs model.
Values are corrected for protein concentration and cytosol fraction to better represent
metabolite concentration in cytosol.
"""
glycolysis_init_conc = Measurements.value.(glycolysis_init_conc_w_uncertainty)

"""
SEM of intracellular concentrations of glycolytic intermediates in mammalian cells (M units).
Values are corrected for protein concentration and cytosol fraction to better represent
metabolite concentration in cytosol.
"""
glycolysis_init_conc_uncertainty = Measurements.uncertainty.(glycolysis_init_conc_w_uncertainty)

tracing_13C_multiplier = 1e-6
"""
Mean ± SEM of intracellular concentrations of 12C and 13C glycolytic intermediates in mammalian cells (M units).
Values are corrected for protein concentration and cytosol fraction to better represent
metabolite concentration in cytosol.
"""
glycolysis_13C_tracing_init_conc_w_uncertainty = LVector(
    Glucose_media_12C = 25e-3,
    Glucose_12C = (4.1e-3 ± 8e-4) / cell_volume_correction,
    G6P_12C = (1.7e-4 ± 3e-5) / cell_volume_correction,
    F6P_12C = (9.8e-5 ± 1.0e-5) / cell_volume_correction,
    F16BP_12C = (5.5e-4 ± 1.0e-4) / cell_volume_correction,
    GAP_12C = (9.6e-5 ± 2.0e-5) / cell_volume_correction,
    DHAP_12C = (8.2e-4 ± 1.9e-4) / cell_volume_correction,
    BPG_12C = (1.1e-6 ± 3e-7) / cell_volume_correction,
    ThreePG_12C = (1.6e-3 ± 6e-4) / cell_volume_correction,
    TwoPG_12C = (4.9e-4 ± 2e-4) / cell_volume_correction,
    PEP_12C = (9.7e-5 ± 1.9e-5) / cell_volume_correction,
    Pyruvate_12C = (6.1e-4 ± 2.3e-4) / cell_volume_correction,
    Lactate_12C = (9.6e-3 ± 2.0e-3) / cell_volume_correction,
    Lactate_media_12C = 2e-3,
    F26BP_12C = 0.0 / cell_volume_correction,
    Citrate_12C = 0.0 / cell_volume_correction,
    Phenylalanine_12C = 0.0 / cell_volume_correction,
    Glucose_media_13C = tracing_13C_multiplier * 25e-3,
    Glucose_13C = tracing_13C_multiplier * (4.1e-3 ± 8e-4) / cell_volume_correction,
    G6P_13C = tracing_13C_multiplier * (1.7e-4 ± 3e-5) / cell_volume_correction,
    F6P_13C = tracing_13C_multiplier * (9.8e-5 ± 1.0e-5) / cell_volume_correction,
    F16BP_13C = tracing_13C_multiplier * (5.5e-4 ± 1.0e-4) / cell_volume_correction,
    GAP_13C = tracing_13C_multiplier * (9.6e-5 ± 2.0e-5) / cell_volume_correction,
    DHAP_13C = tracing_13C_multiplier * (8.2e-4 ± 1.9e-4) / cell_volume_correction,
    BPG_13C = tracing_13C_multiplier * (1.1e-6 ± 3e-7) / cell_volume_correction,
    ThreePG_13C = tracing_13C_multiplier * (1.6e-3 ± 6e-4) / cell_volume_correction,
    TwoPG_13C = tracing_13C_multiplier * (4.9e-4 ± 2e-4) / cell_volume_correction,
    PEP_13C = tracing_13C_multiplier * (9.7e-5 ± 1.9e-5) / cell_volume_correction,
    Pyruvate_13C = tracing_13C_multiplier * (6.1e-4 ± 2.3e-4) / cell_volume_correction,
    Lactate_13C = tracing_13C_multiplier * (9.6e-3 ± 2.0e-3) / cell_volume_correction,
    Lactate_media_13C = tracing_13C_multiplier * 2e-3,
    F26BP_13C = tracing_13C_multiplier * 0.0 / cell_volume_correction,
    Citrate_13C = tracing_13C_multiplier * 0.0 / cell_volume_correction,
    Phenylalanine_13C = tracing_13C_multiplier * 0.0 / cell_volume_correction,
    ATP = (3.6e-3 ± 5e-4) / cell_volume_correction,
    ADP = (7.3e-4 ± 1.5e-4) / cell_volume_correction,
    AMP = (1.6e-4 ± 6e-5) / cell_volume_correction,
    Phosphate = (4.0e-3 ± 9e-4) / cell_volume_correction,
    NAD = (9.2e-4 ± 4.1e-4) / cell_volume_correction,
    NADH = (8.4e-5 ± 3.6e-5) / cell_volume_correction,
)

"""
Intracellular concentrations of 12C and 13C glycolytic intermediates in mammalian cells (M units).
`glycolysis_init_conc` can be used as initial condition for simulating `glycolysis_13C_tracing_ODEs` glycolysis ODEs model.
Values are corrected for protein concentration and cytosol fraction to better represent
metabolite concentration in cytosol.
"""
glycolysis_13C_tracing_init_conc = Measurements.value.(glycolysis_13C_tracing_init_conc_w_uncertainty)

"""
SEM of intracellular concentrations of 12C and 13C glycolytic intermediates in mammalian cells (M units).
Values are corrected for protein concentration and cytosol fraction to better represent
metabolite concentration in cytosol.
"""
glycolysis_13C_tracing_init_conc_uncertainty = Measurements.uncertainty.(glycolysis_13C_tracing_init_conc_w_uncertainty)


#= Values of parameters with uncertainty
Parameters have the following units:
K - M
Vmax - µmol/min per mg protein or U per mg
Conc - mg protein/ul of cytosolic volume (converted from mg protein/mg proteome)
Conc⋅Vmax - M/min
Keq - either unitless or concentration
MW - g / mole /1000
=#

###########################################################################################
# Kinetic and thermodynamic parameters for enzyme rate equations
###########################################################################################

"""
Mean ± SEM of kinetic and thermodynamic parameters for enzyme rate equations.
Use `propertynames(glycolysis_params_w_uncertainty)` for list of all parameters given
as "EnzymeGeneName_ParameterName".

    Parameters have the following units:
    K - M
    β - unitless
    Vmax - µmol/min per mg protein or U per mg
    Keq - depends on reaction stoichiometry
    Conc - mg protein/ul of cytosolic volume (converted from mg protein/mg proteome)
    Conc⋅Vmax - M/min
    Keq - either unitless or concentration
    MW - g / mole /1000
"""
glycolysis_params_w_uncertainty = LVector(
    GLUT_Km_Glucose = 20e-3 ± 4e-3,
    GLUT_Conc = (cell_protein_density / cell_volume_correction) * (1.4e-4 ± 4e-5),
    GLUT_Vmax = 1100.0 ± 300.0,
    GLUT_Keq = 1.0 ± 0.0,
    GLUT_MW = 54084.0 / 1000,
    HK1_K_Glucose = 4.9e-8 ± 0.4e-8,
    HK1_K_a_ATP = 9.3e-7 ± 1e-7,
    HK1_β_Glucose_ATP = 0.0016 ± 0.0002,
    HK1_K_G6P = 480e-6 ± 130e-6,
    HK1_K_a_ADP = 390e-6 ± 50e-6,
    HK1_K_a_G6P_cat = 33e-6 ± 8e-6,
    HK1_K_i_G6P_reg = 6.9e-6 ± 0.9e-6,
    HK1_K_a_Pi = 460e-6 ± 120e-6,
    HK1_Conc = (cell_protein_density / cell_volume_correction) * (6.0e-4 ± 1.2e-4),
    HK1_Vmax = 53.0 ± 5.0,
    HK1_Keq = 2700.0 ± 800.0,
    HK1_MW = 102486.0 / 1000,
    GPI_Km_G6P = 330e-6 ± 100e-6,
    GPI_Km_F6P = 70e-6 ± 20e-6,
    GPI_Conc = (cell_protein_density / cell_volume_correction) * (1.4e-3 ± 1e-4),
    GPI_Vmax = 790.0 ± 190.0,
    GPI_Keq = 0.36 ± 0.11,
    GPI_MW = 63147.0 / 1000,
    PFKP_L = 6.9 ± 4.0,
    PFKP_K_a_F6P = 1.7e-3 ± 0.5e-3,
    PFKP_K_ATP = 250e-6 ± 130e-6,
    PFKP_K_F16BP = 1.5e-3 ± 0.2e-3,
    PFKP_K_ADP = 360e-6 ± 190e-6,
    PFKP_K_Phosphate = 1.1e-3 ± 0.3e-3,
    PFKP_K_a_ADP_reg = 340e-6 ± 100e-6,
    PFKP_K_i_ATP_reg = 1.1e-3 ± 0.3e-3,
    PFKP_K_a_F26BP = 3.4e-7 ± 1.2e-7,
    PFKP_K_i_F26BP = 1.1e-6 ± 0.3e-6,
    PFKP_K_i_Citrate = 3.6e-3 ± 0.3e-3,
    PFKP_Conc = (cell_protein_density / cell_volume_correction) * (8.6e-4 ± 1.2e-4),
    PFKP_Vmax = 80 ± 22,
    PFKP_Keq = 760.0 ± 380.0,
    PFKP_MW = 85596.0 / 1000,
    ALDO_Km_F16BP = 17e-6 ± 5e-6,
    ALDO_Km_DHAP = 430e-6 ± 240e-6,
    ALDO_Kd_DHAP = 12e-6 ± 5e-6,
    ALDO_Km_GAP = 17e-6 ± 14e-6,
    ALDO_Ki_GAP = 1.8e-6 ± 0.7e-6,
    ALDO_Conc = (cell_protein_density / cell_volume_correction) * (3.0e-3 ± 0.3e-3),
    ALDO_Vmax = 29 ± 8,
    ALDO_Keq = 2.7e-6 ± 1.2e-6,
    ALDO_MW = 39420.0 / 1000,
    TPI_Km_DHAP = 470e-6 ± 80e-6,
    TPI_Km_GAP = 11e-6 ± 2e-6,
    TPI_Conc = (cell_protein_density / cell_volume_correction) * (3.0e-3 ± 0.2e-3),
    TPI_Vmax = 6600.0 ± 800.0,
    TPI_Keq = 0.0045 ± 0.0024,
    TPI_MW = 26669.0 / 1000,
    GAPDH_L = 2.9 ± 1.7,
    GAPDH_K_GAP = 1.7e-6 ± 0.1e-6,
    GAPDH_K_a_NAD = 78e-6 ± 9e-6,
    GAPDH_K_i_NAD = 130e-6 ± 34e-6,
    GAPDH_K_a_Phosphate = 1.5e-3 ± 0.4e-3,
    GAPDH_K_i_Phosphate = 5.5e-3 ± 3.8e-3,
    GAPDH_K_BPG = 0.82e-6 ± 0.14e-6,
    GAPDH_K_a_NADH = 12e-6 ± 7e-6,
    GAPDH_K_i_NADH = 2.5e-6 ± 1.1e-6,
    GAPDH_β_i_BPG = 0.22 ± 0.04,
    GAPDH_Conc = (cell_protein_density / cell_volume_correction) * (5.6e-3 ± 0.5e-3),
    GAPDH_Vmax = 130 ± 20,
    # GAPDH_Keq = 15 ± 6,
    GAPDH_Keq = 68 ± 20,
    GAPDH_MW = 36053.0 / 1000,
    PGK_K_BPG = 3e-6 ± 0.7e-6,
    PGK_K_ADP = 42e-6 ± 10e-6,
    PGK_K_ThreePG = 660e-6 ± 220e-6,
    PGK_K_ATP = 580e-6 ± 230e-6,
    PGK_α = 2.1 ± 0.7,
    PGK_β = 0.22 ± 0.11,
    PGK_γ = 2.0 ± 0.8,
    PGK_Conc = (cell_protein_density / cell_volume_correction) * (3.9e-3 ± 0.4e-3),
    PGK_Vmax = 3500.0 ± 900,
    PGK_Keq = 2000 ± 700,
    PGK_MW = 44615.0 / 1000,
    PGM_Km_ThreePG = 270e-6 ± 60e-6,
    PGM_Km_TwoPG = 19e-6 ± 7e-6,
    PGM_Conc = (cell_protein_density / cell_volume_correction) * (2.3e-3 ± 0.2e-3),
    PGM_Vmax = 3000 ± 900,
    PGM_Keq = 0.18 ± 0.05,
    PGM_MW = 28804.0 / 1000,
    ENO_Km_TwoPG = 40e-6 ± 9e-6,
    ENO_Km_PEP = 160e-6 ± 30e-6,
    ENO_Conc = (cell_protein_density / cell_volume_correction) * (9.6e-3 ± 0.5e-3),
    ENO_Vmax = 71 ± 9,
    ENO_Keq = 4.4 ± 1.0,
    ENO_MW = 47169.0 / 1000,
    PKM2_L = 26.0 ± 7.0,
    PKM2_Vmax_a = 540.0 ± 100.0,
    PKM2_Vmax_i = 86.0 ± 41.0,
    PKM2_K_a_PEP = 150e-6 ± 20e-6,
    PKM2_K_ADP = 320e-6 ± 50e-6,
    PKM2_K_Pyruvate = 1800e-6 ± 400e-6,
    PKM2_K_ATP = 940e-6 ± 250e-6,
    PKM2_K_a_F16BP = 1.4e-6 ± 0.4e-6,
    PKM2_K_a_Phenylalanine = 2100e-6 ± 600e-6,
    PKM2_K_i_PEP = 2400e-6 ± 600e-6,
    PKM2_K_i_Phenylalanine = 200e-6 ± 30e-6,
    PKM2_β_i_PEP_ATP = 12.0 ± 5.0,
    PKM2_Keq = 20000.0 ± 6000.0,
    PKM2_Conc = (cell_protein_density / cell_volume_correction) * (7.5e-3 ± 1.0e-3),
    PKM2_MW = 57937.0 / 1000,
    LDH_Km_Pyruvate = 140e-6 ± 20e-6,
    LDH_Km_NADH = 20e-6 ± 2e-6,
    LDH_Km_Lactate = 6.6e-3 ± 0.7e-3,
    LDH_Km_NAD = 220e-6 ± 30e-6,
    LDH_Kd_NADH = 5.7e-6 ± 0.4e-6,
    LDH_Kd_NAD = 810e-6 ± 80e-6,
    LDH_Conc = (cell_protein_density / cell_volume_correction) * (5.1e-3 ± 0.5e-3),
    LDH_Vmax = 320 ± 40,
    LDH_Keq = 13000 ± 5000,
    LDH_MW = 36689.0 / 1000,
    MCT_Km_Lactate = 5.5e-3 ± 0.9e-3,
    MCT_Conc = (cell_protein_density / cell_volume_correction) * (1.5e-4 ± 0.3e-4),
    MCT_Vmax = 470.0 ± 230,
    MCT_Keq = 1.0,
    MCT_MW = 53944.0 / 1000,
    AK_Km_ADP = 0.12e-3,
    AK_Km_ATP = 0.11e-3,
    AK_Km_AMP = 0.09e-3,
    AK_Vmax = 0.03,
    AK_Keq = 0.48,
    AK_MW = 21635.0 / 1000,
    NDPK_Km_ATP = 2e-3,
    NDPK_Km_ADP = 0.1e-3,
    NDPK_Km_NTP = 0.5e-3,
    NDPK_Km_NDP = 0.2e-3,
    NDPK_Vmax = 0.0,
    NDPK_Keq = 2.16, #Average of 2.57 (UDP), 2.71 (GDP) and 1.21 (CTP)
    NDPK_MW = 17149.0 / 1000,
    CK_Km_ATP = 0.16e-3,
    CK_Km_ADP = 0.11e-3,
    CK_Km_Phosphocreatine = 0.9e-3,
    CK_Km_Creatine = 2.3e-3,
    CK_Vmax = 0.0,
    CK_Keq = 0.00598 ± 0.00093,
    CK_MW = 17149.0 / 1000,
    ATPase_Km_ATP = 1e-6,
    ATPase_Km_ADP = 1e-3,
    ATPase_Km_Phosphate = 1e-3,
    ATPase_Keq = 83000.0,
    ATPase_Vmax = 0.0002,
)

"""
LArray containing mean of kinetic and thermodynamic parameters for enzyme rate equations.
Use `propertynames(glycolysis_params)` for list of all parameters given
as "EnzymeGeneName_ParameterName".

    Parameters have the following units:
    K - M
    β - unitless
    Vmax - µmol/min per mg protein or U per mg
    Keq - depends on reaction stoichiometry
    Conc - mg protein/ul of cytosolic volume (converted from mg protein/mg proteome)
    Conc⋅Vmax - M/min
    Keq - either unitless or concentration
    MW - g / mole /1000
"""
glycolysis_params = Measurements.value.(glycolysis_params_w_uncertainty)

"""
SEM of kinetic and thermodynamic parameters for enzyme rate equations.
Use `propertynames(glycolysis_params_uncertainty)` for list of all parameters given
as "EnzymeGeneName_ParameterName".

    Parameters have the following units:
    K - M
    β - unitless
    Vmax - µmol/min per mg protein or U per mg
    Keq - depends on reaction stoichiometry
    Conc - mg protein/ul of cytosolic volume (converted from mg protein/mg proteome)
    Conc⋅Vmax - M/min
    Keq - either unitless or concentration
    MW - g / mole /1000
"""
glycolysis_params_uncertainty = Measurements.uncertainty.(glycolysis_params_w_uncertainty)
