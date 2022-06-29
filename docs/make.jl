using Documenter, SNMRForward, Literate

fmscript = joinpath("..", "scripts", "forward_model.jl")
outdir = joinpath(@__DIR__, "src", "generated")
Literate.markdown(fmscript, outdir; execute=true, documenter=true)

makedocs(sitename="SNMRForward",
   pages = ["index.md",
            "Examples" => ["generated/forward_model.md"],
            "API docs" => ["forward_modelling.md", "constants.md"]])


