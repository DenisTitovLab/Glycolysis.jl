using LabelledArrays

minimal_model_initial_concentrations =
    LVector(Glu = 25e-3, X = 5e-3, Lac = 1e-3, ATP = 9.9e-3, ADP = 0.1e-3, Phosphate = 5e-3)
minimal_model_params = LVector(
    Enz1_L = 1e-3,
    Enz1_n = 4,
    Enz1_Vmax = 1e-3,
    Enz1_Keq = 1e3,
    Enz2_Vmax = 1e0,
    Enz2_Keq = 1e3,
    ATPase_Vmax = 1e-4,
    ATPase_Keq = 1e3,
)

function rate_Enz1(Glu, ATP, X, ADP, Phosphate, params)
    Γ = (X * ADP) / (Glu * ATP)
    Rate =
        (params.Enz1_Vmax / (1 + (params.Enz1_L * ATP / (ADP * Phosphate))^params.Enz1_n)) *
        (1 - Γ / params.Enz1_Keq)
    return Rate
end

function rate_Enz2(X, Phosphate, ADP, Lac, ATP, params)
    Γ = (Lac * ATP^2) / (X * Phosphate * ADP^2)
    Rate = params.Enz2_Vmax * (1 - Γ / params.Enz2_Keq)
    return Rate
end

function minimal_model_rateATPase(ATP, ADP, Phosphate, params)
    Γ = (ADP * Phosphate) / (ATP)
    Rate = params.ATPase_Vmax * (1 - Γ / params.ATPase_Keq)
    return Rate
end

function minimal_glycolysis_ODEs(ds, s, params, t)
    ds.Glu = 0.0
    ds.X = (
        rate_Enz1(s.Glu, s.ATP, s.X, s.ADP, s.Phosphate, params) -
        rate_Enz2(s.X, s.Phosphate, s.ADP, s.Lac, s.ATP, params)
    )
    ds.Lac = 0
    ds.ATP = (
        -rate_Enz1(s.Glu, s.ATP, s.X, s.ADP, s.Phosphate, params) +
        2 * rate_Enz2(s.X, s.Phosphate, s.ADP, s.Lac, s.ATP, params) -
        minimal_model_rateATPase(s.ATP, s.ADP, s.Phosphate, params)
    )
    ds.ADP = (
        rate_Enz1(s.Glu, s.ATP, s.X, s.ADP, s.Phosphate, params) -
        2 * rate_Enz2(s.X, s.Phosphate, s.ADP, s.Lac, s.ATP, params) +
        minimal_model_rateATPase(s.ATP, s.ADP, s.Phosphate, params)
    )
    ds.Phosphate =
        minimal_model_rateATPase(s.ATP, s.ADP, s.Phosphate, params) -
        rate_Enz2(s.X, s.Phosphate, s.ADP, s.Lac, s.ATP, params)
end

function minimal_glycolysis_conc_to_disequilibrium_ratios(s, params)
    k = @LVector eltype(s) (:Q_Keq_Enz1, :Q_Keq_Enz2, :Q_Keq_ATPase)
    k.Q_Keq_Enz1 = ((s.X * s.ADP) / (s.Glu * s.ATP)) / params.Enz1_Keq
    k.Q_Keq_Enz2 = ((s.Lac * s.ATP^2) / (s.X * s.Phosphate * s.ADP^2)) / params.Enz2_Keq
    k.Q_Keq_ATPase = ((s.Phosphate * s.ADP) / s.ATP) / params.ATPase_Keq
    return k
end

function minimal_glycolysis_conc_to_rates(s, params)
    r = @LVector eltype(s) (:Enz1, :Enz2, :ATPprod, :ATPase)
    r.Enz1 = rate_Enz1(s.Glu, s.ATP, s.X, s.ADP, s.Phosphate, params)
    r.Enz2 = rate_Enz2(s.X, s.Phosphate, s.ADP, s.Lac, s.ATP, params)
    r.ATPprod = (
        -rate_Enz1(s.Glu, s.ATP, s.X, s.ADP, s.Phosphate, params) +
        2 * rate_Enz2(s.X, s.Phosphate, s.ADP, s.Lac, s.ATP, params)
    )
    r.ATPase = minimal_model_rateATPase(s.ATP, s.ADP, s.Phosphate, params)
    return r
end