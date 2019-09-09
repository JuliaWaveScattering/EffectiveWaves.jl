using Documenter, EffectiveWaves

makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    debug = true,
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
    repo = "github.com/arturgower/EffectiveWaves.jl.git",
    # forcepush = true,
    branch = "gh-pages",
    latest = "master",
    target = "build"
)
