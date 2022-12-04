using LabelledArrays, Measurements
#=
To do list

- add Km for products for PFKP and PKM2
- rename params to something else as it conflicts with too many other library params

=#

# Define cellular parameters and kinetic+thermodynamic constants of enzymes
water_fraction_cell_volume = 0.66
cytosol_fraction_cell_volume = 0.66
cell_volume_correction = water_fraction_cell_volume * cytosol_fraction_cell_volume
cell_protein_density = 0.2 #mg/ul

# Define initial concentrations of metabolites for model simulations
glycolysis_init_conc = LVector(
    Glucose_media = 25e-3,
    Glucose = 6.12e-3 / cell_volume_correction,
    G6P = 1.73e-4 / cell_volume_correction,
    F6P = 8.17e-5 / cell_volume_correction,
    F16BP = 8.75e-4 / cell_volume_correction,
    GAP = 1.12e-4 / cell_volume_correction,
    DHAP = 7.19e-4 / cell_volume_correction,
    BPG = 1.3e-6 / cell_volume_correction,
    ThreePG = 2.21e-4 / cell_volume_correction,
    TwoPG = 2.01e-5 / cell_volume_correction,
    PEP = 4.29e-5 / cell_volume_correction,
    Pyruvate = 6.77e-4 / cell_volume_correction,
    Lactate = 2.74e-3 / cell_volume_correction,
    Lactate_media = 1e-3 / cell_volume_correction,
    ATP = 3.32e-3 / cell_volume_correction,
    ADP = 3.96e-4 / cell_volume_correction,
    AMP = 9.41e-5 / cell_volume_correction,
    Phosphate = 1e-3 / cell_volume_correction,
    NAD = 0.346e-3 / cell_volume_correction,
    NADH = 3.58e-5 / cell_volume_correction,
    F26BP = 0.0 / cell_volume_correction,
    Citrate = 0.0 / cell_volume_correction,
)

tracing_13C_divider = 1e-3
glycolysis_13C_tracing_init_conc = LVector(
    Glucose_media_12C = 25e-3,
    Glucose_12C = 6.12e-3 / cell_volume_correction,
    G6P_12C = 1.73e-4 / cell_volume_correction,
    F6P_12C = 8.17e-5 / cell_volume_correction,
    F16BP_12C = 8.75e-4 / cell_volume_correction,
    GAP_12C = 1.12e-4 / cell_volume_correction,
    DHAP_12C = 7.19e-4 / cell_volume_correction,
    BPG_12C = 1.3e-6 / cell_volume_correction,
    ThreePG_12C = 2.21e-4 / cell_volume_correction,
    TwoPG_12C = 2.01e-5 / cell_volume_correction,
    PEP_12C = 4.29e-5 / cell_volume_correction,
    Pyruvate_12C = 6.77e-4 / cell_volume_correction,
    Lactate_12C = 2.74e-3 / cell_volume_correction,
    Lactate_media_12C = 1e-3 / cell_volume_correction,
    F26BP_12C = 0.0 / cell_volume_correction,
    Citrate_12C = 0.0 / cell_volume_correction,
    
    Glucose_media_13C = tracing_13C_divider * 25e-3,
    Glucose_13C = tracing_13C_divider * 6.12e-3 / cell_volume_correction,
    G6P_13C = tracing_13C_divider * 1.73e-4 / cell_volume_correction,
    F6P_13C = tracing_13C_divider * 8.17e-5 / cell_volume_correction,
    F16BP_13C = tracing_13C_divider * 8.75e-4 / cell_volume_correction,
    GAP_13C = tracing_13C_divider * 1.12e-4 / cell_volume_correction,
    DHAP_13C = tracing_13C_divider * 7.19e-4 / cell_volume_correction,
    BPG_13C = tracing_13C_divider * 1.3e-6 / cell_volume_correction,
    ThreePG_13C = tracing_13C_divider * 2.21e-4 / cell_volume_correction,
    TwoPG_13C = tracing_13C_divider * 2.01e-5 / cell_volume_correction,
    PEP_13C = tracing_13C_divider * 4.29e-5 / cell_volume_correction,
    Pyruvate_13C = tracing_13C_divider * 6.77e-4 / cell_volume_correction,
    Lactate_13C = tracing_13C_divider * 2.74e-3 / cell_volume_correction,
    Lactate_media_13C = tracing_13C_divider * 1e-3 / cell_volume_correction,
    F26BP_13C = tracing_13C_divider * 0.0 / cell_volume_correction,
    Citrate_13C = tracing_13C_divider * 0.0 / cell_volume_correction,

    ATP = 3.32e-3 / cell_volume_correction,
    ADP = 3.96e-4 / cell_volume_correction,
    AMP = 9.41e-5 / cell_volume_correction,
    Phosphate = 1e-3 / cell_volume_correction,
    NAD = 0.346e-3 / cell_volume_correction,
    NADH = 3.58e-5 / cell_volume_correction,
)


#= Values of parameters with uncertainty
Parameters have the following units:
K - M
Vmax - umole/min per mg protein or U per mg
Conc - mg protein/ul of cytosolic volume (converted from mg protein/mg proteome)
Conc⋅Vmax - M/min
Keq - either unitless or concentration
MW - g / mole /1000
=#

glycolysis_params_w_uncertainty = LVector(
    # Intermediate_cons_frac = 0.0,
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
    HK1_n = 1.0, # To be removed once bound equations are edited
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
    PFKP_Conc = (cell_protein_density / cell_volume_correction) * (8.6e-4 ± 1.2e-4), #Enzyme abundance in mg enzyme per ul cellular protein
    PFKP_Vmax = 80 ± 22, #Max forward rate of active PFKP form enzyme, umole/min per mg enzyme
    PFKP_Keq = 760.0 ± 380.0, #Equilibrium constant of PFK reaction at pH=7.4 and Ionic Strength=0.1M
    PFKP_MW = 85596.0 / 1000, #Molecular weight of enzyme, mg/umol
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
    GAPDH_L = 1.6 ± 0.3,
    GAPDH_K_GAP = 1.7e-6 ± 0.1e-6,
    GAPDH_K_a_NAD = 83e-6 ± 5e-6,
    GAPDH_K_i_NAD = 180e-6 ± 20e-6,
    GAPDH_K_a_Phosphate = 1.9e-3 ± 0.1e-3,
    GAPDH_K_i_Phosphate = 10e-3 ± 4e-3,
    GAPDH_K_BPG = 0.78e-6 ± 0.1e-6,
    GAPDH_K_a_NADH = 7.8e-6 ± 1.2e-6,
    GAPDH_K_i_NADH = 1.8e-6 ± 0.3e-6,
    GAPDH_α_i_BPG = 4.5 ± 0.7,
    GAPDH_Conc = (cell_protein_density / cell_volume_correction) * (5.6e-3 ± 0.5e-3),
    GAPDH_Vmax = 130 ± 20,
    GAPDH_Keq = 16 ± 5,
    GAPDH_MW = 36053.0 / 1000,

    # Below GAPDH constants to be removed after adjusting bound equations
    GAPDH_Km_Phosphate = 7.1e-4, #Km, M
    GAPDH_Km_NAD = 0.047e-3, #Km, M
    GAPDH_Km_NADH = 0.01e-3, #Km, M
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
    ENO_Conc = (cell_protein_density / cell_volume_correction) * (9.6e-3 ± 0.5e-3), #Enzyme abundance in mg enzyme per ul cellular protein
    ENO_Vmax = 71 ± 9,
    ENO_Keq = 4.4 ± 1.0,
    ENO_MW = 47169.0 / 1000,
    PKM2_L = exp(-(-9.06346122e+00)),
    #PKM2_L = 0.0,
    PKM2_a_KmPEP = 7.81331145e-05,  #Km of PKM2 active form for PEP, M
    PKM2_a_KmADP = 9.10888386e-05,  #Km of PKM2 active form for ADP, M
    PKM2_a_KdF16BP = 10e-6,  #Kd of PKM2 active form for allosteric activator F16BP, M
    #PKM2_a_KdF16BP = Inf,  #Kd of PKM2 active form for allosteric activator F16BP, M
    PKM2_a_KdATP = 4.13743762e-04,  #Kd of PKM2 active form for inhibitor ATP, M
    PKM2_i_KmPEP = 2.04925488e-03,  #Km of PKM2 inactive form for PEP, M
    PKM2_i_KmADP = 1.41515745e-04,  #Km of PKM2 inactive form for ADP, M
    PKM2_i_KdF16BP = 1000.0,  #Kd of PKM2 inactive form for allosteric activator F16BP, M
    #PKM2_i_KdF16BP = Inf,  #Kd of PKM2 inactive form for allosteric activator F16BP, M
    PKM2_i_KdATP = 4.26282043e-04,  #Kd of PKM2 active form for inhibitor ATP, M
    PKM2_Conc = (cell_protein_density / cell_volume_correction) * (7.5e-3 ± 1.0e-3), #Enzyme abundance in mg enzyme per ul cellular protein
    PKM2_a_Vmax = 220.0,  #Max forward rate of active PKM2 form enzyme, M/min OR umole/min per ul of cell volume
    PKM2_i_Vmax = 220.0,  #Max forward rate of inactive PKM2 form enzyme, M/min OR umole/min per ul of cell volume
    PKM2_Vmax = 220.0,  #Max forward rate of active PKM2 form enzyme, M/min OR umole/min per ul of cell volume
    PKM2_n = 4.0,  #number of subunits/oligomeric state of PKM2 in solution
    PKM2_Keq = 31000.0,  #Equilibrium constant of reaction at pH=7.4 and Ionic Strength=0.1M
    PKM2_MW = 57937.0 / 1000, #Molecular weight of enzyme, mg/umol
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
    MCT_Conc = (cell_protein_density / cell_volume_correction) * (1.5e-4 ± 0.3e-4), #Enzyme abundance in mg enzyme per ul cellular protein
    MCT_Vmax = 470.0 ± 230, #Max forward rate of enzyme, umole/min per mg of enzyme
    MCT_Keq = 1.0, #Equilibrium constant of reaction at pH=7.4 and Ionic Strength=0.1M
    MCT_MW = 53944.0 / 1000, #Molecular weight of enzyme, mg/umol
    AK_Km_ADP = 0.12e-3, #Km, M
    AK_Km_ATP = 0.11e-3, #Km, M
    AK_Km_AMP = 0.09e-3, #Km, M
    AK_Vmax = 0.03, #M/min OR umole/min per ul of cell volume. Calculated based on estimate of 10 umole/mg (7-200 range) per hour of dryg weight in tissues
    AK_Keq = 0.48, #Equilibrium constant of enzyme reaction
    AK_MW = 21635.0 / 1000, #Molecular weight of enzyme, mg/umol
    ATPase_Km_ATP = 1e-6, #Km, M
    ATPase_Km_ADP = 1e-3, #Km, M
    ATPase_Km_Phosphate = 1e-3, #Km, M
    ATPase_Keq = 83000.0, #Equilibrium constant of enzyme reaction
    ATPase_Vmax = 0.0002, #Km, M
)

glycolysis_params = Measurements.value.(glycolysis_params_w_uncertainty)
glycolysis_params_uncertainty = Measurements.uncertainty.(glycolysis_params_w_uncertainty)
