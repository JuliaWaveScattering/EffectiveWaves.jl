using Documenter, EffectiveWaves

makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    debug = true,
    sitename ="EffectiveWaves.jl",
    authors = "Artur L. Gower",
    source = "src",
    # modules = [EffectiveWaves, MultipleScattering],
    pages = [
        "Home" =>"index.md",
        "Manual" => [
            "manual/introduction.md",
            "manual/background.md",
            "manual/wavenumbers.md",
            "manual/reflection.md",
            "manual/sphere.md"
        ],
        "Library" => "library/library.md"
        ,
        "Examples" => [
            "examples/vary_two_species/README.md",
            "examples/equivalent-symmetries/README.md",
            "examples/matched_method/README.md"
        ]
    ]
)

deploydocs(
    branch = "gh-pages",
    target = "build",
    versions = ["stable" => "v^", "v#.#.#"],
    # forcepush = true,
    repo = "github.com/arturgower/EffectiveWaves.jl.git"
)
