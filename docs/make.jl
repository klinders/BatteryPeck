using Documenter, BatteryToolkit, DocumenterCitations, DocumenterVitepress

bib = CitationBibliography(
    joinpath(@__DIR__, "src/assets", "BatteryToolkit.bib");
    style=:numeric
)

makedocs(
    modules=[BatteryToolkit],
    authors="Koen Linders",
    sitename="BatteryToolkit.jl";
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/klinders/BatteryToolkit",
        devbranch = "main",
        devurl = "dev",
    ),
    plugins=[bib],
    pages=[
        "Home" => "index.md",
        "Getting Started" => [
            "Quick Start" => "guide/quickstart.md",
            "Parameter Sets" => "guide/parameters.md",
            "Experiments" => "guide/experiments.md",
            "Examples" => "guide/examples.md",
        ],
        "Models" => [
            "SPMe Cell Model" => "models/spme.md",
            "Pack Models" => "models/pack-models.md",
            "Side Reactions" => "models/side-reactions.md",
        ],
        "API Reference" => [
            "Parameters" => "api/parameters.md",
            "Cell Models" => "api/cellmodels.md",
            "Pack Models" => "api/packmodels.md",
            "Experiments" => "api/experiments.md",
            "Solvers" => "api/solvers.md",
        ],
        "Advanced" => [
            "Finite Volume Method" => "Finite Volume Method/index.md",
        ],
        "References" => "references.md"
    ]
)

deploydocs(
    repo = "github.com/klinders/BatteryToolkit.git",
    devbranch = "dev"
)