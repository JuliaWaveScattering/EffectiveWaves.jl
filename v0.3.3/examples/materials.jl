# The material below are explicity exported with the package EffectiveWaves.jl, so no need to include this file.

## possible materials

# note that many of these materials below are automatically loaded, see "src/EffectiveWaves.jl".

# http://www.rfcafe.com/references/general/velocity-sound-media.htm

  Brick        = Acoustic(3; ρ=1.80*1000, c = 3650.0)
  IronArmco    = Acoustic(3; ρ=7.85*1000, c = 5960.0)
  LeadAnnealed = Acoustic(3; ρ=11.4*1000, c = 2160.)
  RubberGum    = Acoustic(3; ρ=0.95*1000, c = 1550.0)
  FusedSilica  = Acoustic(3; ρ=2.2*1000,  c = 5760.0)
  GlassPyrex   = Acoustic(3; ρ=2.32*1000, c = 5640.0)
  ClayRock     = Acoustic(3; ρ=2.2*1000,  c = 3480.0)

# Fluids # http://www.rshydro.co.uk/sound-speeds/
  WaterDistilled= Acoustic(3; ρ=0.998*1000, c = 1496.0)
  Glycerol      = Acoustic(3; ρ=1.26*1000,  c = 1904.0)
  Hexadecane    = Acoustic(3; ρ=0.773*1000, c = 1338.0)
  Acetone       = Acoustic(3; ρ=0.791*1000, c = 1174.0)
  Benzene       = Acoustic(3; ρ=0.87*1000,  c = 1295.0)
  Nitrobenzene  = Acoustic(3; ρ=1.204*1000, c = 1415.0)
  OliveOil      = Acoustic(3; ρ=0.912*1000, c = 1431.0)
  SodiumNitrate = Acoustic(3; ρ=1.8*1000,   c = 1760.)
# gases
  AirDry       = Acoustic(3; ρ=1.293,     c = 331.4)

# non-porous minerals, i.e. add air-bubbles to be more realistic
# https://www.spec2000.net/05-mineralprops.htm
  ρ=1.58*1000;
  Clay         = Acoustic(3; ρ=ρ, c = p_speed(ρ,1.5*10^9,1.4*10^9))
  ρ=2.65*1000;
  SilicaQuartz = Acoustic(3; ρ=ρ, c = p_speed(ρ,36.6*10^9,45.0*10^9))
  ρ=2.71*1000;
  Calcite = Acoustic(3; ρ=ρ, c = p_speed(ρ,70.8*10^9,32.0*10^9)) # approximately limestone
  LimeStone    = Acoustic(3; ρ=2460.0,  c = 4855.0) # 4% porosity, base level http://www.sciencedirect.com/science/article/pii/S1365160915300101
