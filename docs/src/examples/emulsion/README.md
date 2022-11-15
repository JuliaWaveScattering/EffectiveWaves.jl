
## Code to generate the figures below, which compare different approximations for the effective wavenumber in a 2D emulsion. If you are new to Julia, please open the ipynb files above, which show how to run the code in your browser on [JuliaBox](https://www.juliabox.com/).

The examples below only run for Julia 0.6 and the package tag v0.1.0. To obtain the version of this package to run these example run in julia:
```julia
Pkg.clone("https://github.com/arturgower/EffectiveWaves.jl.git") 
Pkg.checkout("EffectiveWaves.jl", "master@v0.1.0")
using EffectiveWaves
```

[Code: fluid_species.ipynb](fluid_species.ipynb)
![Compare effective wavenumber for 2D emulsion](media/compare_fluid_small.png)
![Compare effective wavenumber for 2D emulsion](media/compare_fluid.png)

[Code: fluid_species_large-freq.ipynb](fluid_species_large-freq.ipynb)
![Compare effective wavenumber for 2D emulsion](media/compare_fluid_large-w.png)

[Code: fluid_species_volfrac.ipynb](fluid_species_volfrac.ipynb)
![Compare effective wavenumber for 2D emulsion](media/compare_fluid_volfrac.png)
