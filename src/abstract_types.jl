## A the top of the type tree we have:

# """
#     PairCorrelationType

# A type used to specify what type of pair correlation is to be used. This is like a tag, or a option, to specify which pair correlation is wanted.
# """
# abstract type PairCorrelationType end

# """
#     PairCorrelation

# A type used to store a pair-correlation. This represents the calculated pair-correlation.
# """
# abstract type PairCorrelation end

"""
    Microstructure{Dim}

Every [`Material`](@ref) has a shape and has a microstructure. Currently [`ParticulateMicrostructure`](@ref) is the only example microstructure.
"""
abstract type Microstructure{Dim} end

abstract type AbstractWaveMode end

abstract type AbstractRegularWaveMode <: AbstractWaveMode end
