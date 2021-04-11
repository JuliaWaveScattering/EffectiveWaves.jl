discrete_system(ω::AbstractFloat, source::Source, material::Material; kws...) = discrete_system(ω::AbstractFloat, source::Source, material::Material, setupsymmetry(source,material); kws...)
