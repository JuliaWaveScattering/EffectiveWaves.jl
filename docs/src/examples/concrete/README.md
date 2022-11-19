
## Code to generate the figures below, which compare different approximations for the effective wavenumber in a 2D concrete. If you are new to Julia, please open the ipynb files above, which show how to run the code in your browser on [JuliaBox](https://www.juliabox.com/).

The examples below only run for Julia 0.6 and the package tag v0.1.0. To obtain the version of this package to run these example run in julia:
```julia
Pkg.clone("https://github.com/arturgower/EffectiveWaves.jl.git") 
Pkg.checkout("EffectiveWaves.jl", "master@v0.1.0")
using EffectiveWaves
```

[Code: concrete_species.ipynb](concrete_species.ipynb)
![Compare effective wavenumber for 2D concrete](media/compare_concrete.png)
![Compare effective wavenumber for 2D concrete](media/compare_concrete_zoom.png)

[Code: concrete_species-large_freq.ipynb](concrete_species_large-freq.ipynb)
![Compare effective wavenumber for 2D concrete](media/compare_concrete_large-w.png)

[Code: concrete_species_volfrac.ipynb](concrete_species_volfrac.ipynb)
![Compare effective wavenumber for 2D concrete](media/compare_concrete_volfrac.png)
