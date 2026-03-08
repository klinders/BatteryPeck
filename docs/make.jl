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
        devbranch = "main", # or master, trunk, ...
        devurl = "dev",
    ),
    plugins=[bib],
    pages=[
        "Home" => "index.md",
        "FVM" => "Finite Volume Method/index.md",
        "References"=>"references.md"
    ]
)