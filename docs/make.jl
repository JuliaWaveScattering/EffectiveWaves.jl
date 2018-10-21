using Documenter, EffectiveWaves

cp("../example","./example")

makedocs(
    format=:html,
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
            "examples/intro/README.md",
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
    make = nothing,
    deps = nothing,
    repo = "github.com/arturgower/EffectiveWaves.jl.git"
)
