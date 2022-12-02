#Glycolysis Model ODE system

function glycolysis_ODEs(ds, s, params, t)
    ds.Glucose_media = 0
    ds.Glucose =
        rateGLUT(s.Glucose_media, s.Glucose, params) -
        rateHK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params)
    ds.G6P = rateHK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) - rateGPI(s.G6P, s.F6P, params)
    ds.F6P = (
        rateGPI(s.G6P, s.F6P, params) -
        ratePFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params)
    )
    ds.F16BP = (
        ratePFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) -
        rateALDO(s.F16BP, s.GAP, s.DHAP, params)
    )
    ds.GAP = (
        rateALDO(s.F16BP, s.GAP, s.DHAP, params) + rateTPI(s.GAP, s.DHAP, params) -
        rateGAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    )
    ds.DHAP = rateALDO(s.F16BP, s.GAP, s.DHAP, params) - rateTPI(s.GAP, s.DHAP, params)
    ds.BPG =
        rateGAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params) -
        ratePGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params)
    ds.ThreePG = ratePGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) - ratePGM(s.ThreePG, s.TwoPG, params)
    ds.TwoPG = ratePGM(s.ThreePG, s.TwoPG, params) - rateENO(s.TwoPG, s.PEP, params)
    ds.PEP = rateENO(s.TwoPG, s.PEP, params) - ratePKM2(s.PEP, s.ADP, s.F16BP, s.ATP, s.Pyruvate, params)
    ds.Pyruvate =
        ratePKM2(s.PEP, s.ADP, s.F16BP, s.ATP, s.Pyruvate, params) -
        rateLDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params)
    ds.Lactate =
        rateLDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params) - rateMCT(s.Lactate, s.Lactate_media, params)
    ds.Lactate_media = 0
    ds.ATP = (
        -rateHK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) -
        ratePFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) +
        ratePGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) +
        ratePKM2(s.PEP, s.ADP, s.F16BP, s.ATP, s.Pyruvate, params) -
        rateATPase(s.ATP, s.ADP, s.Phosphate, params) + rateAK(s.ATP, s.ADP, s.AMP, params)
    )
    ds.ADP = (
        rateHK1(s.Glucose, s.G6P, s.ATP, s.Phosphate, s.ADP, params) +
        ratePFKP(s.F6P, s.ATP, s.F16BP, s.ADP, s.Phosphate, s.Citrate, s.F26BP, params) -
        ratePGK(s.BPG, s.ADP, s.ATP, s.ThreePG, params) -
        ratePKM2(s.PEP, s.ADP, s.F16BP, s.ATP, s.Pyruvate, params) +
        rateATPase(s.ATP, s.ADP, s.Phosphate, params) - 2 * rateAK(s.ATP, s.ADP, s.AMP, params)
    )
    ds.AMP = rateAK(s.ATP, s.ADP, s.AMP, params)
    ds.Phosphate =
        rateATPase(s.ATP, s.ADP, s.Phosphate, params) -
        rateGAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    #ds.Phosphate = 0

    ds.NAD =
        rateLDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params) -
        rateGAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params)
    ds.NADH =
        rateGAPDH(s.GAP, s.NAD, s.Phosphate, s.BPG, s.NADH, params) -
        rateLDH(s.Pyruvate, s.NADH, s.NAD, s.Lactate, params)
    ds.F26BP = 0
    ds.Citrate = 0
end