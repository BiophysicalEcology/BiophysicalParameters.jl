# Compile a biophysical parameters database from raw traits.build data and query it.
#
# Prerequisites:
#   ENV["HEATBUDGETDB_PATH"]       = path to the root of the heat budget databases
#   ENV["BIOPHYSICAL_PARAMS_PATH"] = directory where parameters.arrow will be written
#
# The example builds parameters for two taxa (one squamate, one marsupial),
# then loads them back and inspects point estimates and fitted distributions.

using BiophysicalParameters
using RData            # activates the RData extension for loading .rds files
using Distributions    # for inspecting reconstructed distribution objects
using Statistics       # for mean(dist)

# ── 1. Build the compiled database ────────────────────────────────────────────

taxa = ["Pogona vitticeps", "Myrmecobius fasciatus"]

taxon_classes = Dict(
    "Pogona vitticeps"    => "Squamate",
    "Myrmecobius fasciatus" => "EutherianMammal",
)

result = compile_parameters_database(
    taxa;
    taxon_classes,
)

println("Compiled $(nrow(result)) parameter rows for $(length(taxa)) taxa.")
println(first(result, 5))

# ── 2. Query a ParameterRecord ────────────────────────────────────────────────

record = load_parameters("Pogona vitticeps")

println("\nParameterRecord for: ", record.taxon_name)
println("  Struct keys: ", collect(keys(record.parameters)))

# Radiation: per-individual solar absorptivity
dorsal = record[:RadiationParameters][:body_absorptivity_dorsal]
println("\nDorsal absorptivity:")
println("  mean = ", round(dorsal.mean; digits=3))
println("  sd   = ", round(dorsal.sd; digits=4))
println("  n    = ", dorsal.n)
println("  95% CI: [",
    round(dorsal.ci_lower_95; digits=3), ", ",
    round(dorsal.ci_upper_95; digits=3), "]")
println("  distribution family: ", dorsal.distribution_family)

# Reconstruct and inspect the fitted Beta distribution
dist = distribution(dorsal)
println("  distribution: ", dist)
println("  distribution mean: ", round(mean(dist); digits=3))
println("  distribution std:  ", round(std(dist); digits=4))
@assert dist isa Beta

# Convenience access via ParameterRecord
dist2 = distribution(record, :RadiationParameters, :body_absorptivity_dorsal)
@assert dist2 isa Beta

# ── 3. Raw Arrow access via CompiledParametersDB ──────────────────────────────

db = CompiledParametersDB()   # reads path from ENV["BIOPHYSICAL_PARAMS_PATH"]
df = gettraits(db; taxon="Pogona vitticeps")
println("\nRaw Arrow rows for Pogona vitticeps: ", nrow(df))
println(df[:, [:parameter_struct, :parameter_field, :mean, :n, :distribution_family]])

# ── 4. Numbat (marsupial endotherm) — metabolic parameters ───────────────────

numbat = load_parameters("Myrmecobius fasciatus")
if haskey(numbat.parameters, :MetabolismParameters)
    qten = numbat[:MetabolismParameters][:q10]
    println("\nNumbat Q10: ", round(qten.mean; digits=2),
            " (provenance: ", qten.provenance, ")")
end
