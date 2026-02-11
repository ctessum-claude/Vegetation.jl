using Documenter
using Vegetation

DocMeta.setdocmeta!(Vegetation, :DocTestSetup, :(using Vegetation); recursive = true)

makedocs(;
    modules = [Vegetation],
    authors = "EarthSciML authors and contributors",
    sitename = "Vegetation.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://earthsciml.github.io/Vegetation.jl",
        repolink = "https://github.com/earthsciml/Vegetation.jl",
        size_threshold = 500 * 2^10,
    ),
    pages = [
        "Home" => "index.md",
        "Scheller & Mladenoff (2004)" => [
            "LANDIS Biomass Module" => "landis_biomass.md",
        ],
        "Stage (1973)" => [
            "Prognosis Model" => "stage_prognosis.md",
        ],
        "API" => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/earthsciml/Vegetation.jl",
)
