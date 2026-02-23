using Documenter, BatteryToolkit, DocumenterCitations

bib = CitationBibliography(
    joinpath(@__DIR__, "src/assets", "BatteryToolkit.bib");
    style=:numeric
)

makedocs(sitename="BatteryToolkit.jl";  plugins=[bib])