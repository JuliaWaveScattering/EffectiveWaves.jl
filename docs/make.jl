using Documenter, EffectiveWaves

# cp("../example","example")

makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    sitename="EffectiveWaves.jl",
    authors = "Artur L. Gower",
    source= "src",
    modules=[EffectiveWaves],
    pages=[
        "Home" =>"index.md",
        "Manual" => [
            "base.md"
        ]
        ,
        "Examples" => [
            "examples/vary_two_species/README.md",
            "examples/many_wavenumbers/README.md",
            "examples/matched_method/README.md"
        ]
    ]
)

deploydocs(
    branch = "gh-pages",
    latest = "master",
    julia = "1.0",
    osname = "linux",
    target = "build",
    repo = "github.com/arturgower/EffectiveWaves.jl.git"
)
