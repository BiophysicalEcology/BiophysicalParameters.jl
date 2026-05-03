# Lizard solar absorptivity example
#
# Computes per-individual solar-weighted absorptivity for Pogona vitticeps from
# raw spectral reflectance data (Smith et al. 2016), compares against the
# pre-computed values in the radiationDB parameters export, and builds a
# RadiationParameters struct with provenance tracking.

using BiophysicalParameters
using CSV
using DataFrames
using Statistics

const RADIATION_DB = joinpath(
    homedir(), "Dropbox", "Current Research Projects",
    "trait_database", "heat_budget_databases", "radiationDB",
)

# ── 1. Load raw spectral data (Smith et al. 2016) ───────────────────────────
spectral_data = CSV.read(
    joinpath(RADIATION_DB, "data", "Smith_etal_2016", "data.csv"),
    DataFrame,
)
println("Raw spectral data: $(nrow(spectral_data)) rows, $(length(unique(spectral_data.Individual_ID))) individuals")
println("Wavelength range: $(minimum(spectral_data.wave_length_nm))–$(maximum(spectral_data.wave_length_nm)) nm")
println("Body temperatures: $(sort(unique(spectral_data.temp_C))) °C")
println("Regions: $(unique(spectral_data.region_of_measurement))")
println()

# ── 2. Compute per-individual solar-weighted absorptivities ──────────────────
#
# Smith et al. measured each lizard at two body temperatures (15 °C and 40 °C).
# Bearded dragons change colour with temperature — dorsal absorptivity is higher
# at 15 °C (dark) than at 40 °C (pale), while ventral is stable.
println("── Per-individual absorptivity ──────────────────────────────────────")
individual_absorptivities = per_individual_absorptivity(spectral_data)
println(individual_absorptivities)
println()

# Summary by region and temperature
summary = combine(
    groupby(individual_absorptivities, [:region_of_measurement, :temperature_celsius]),
    :absorptivity => mean => :mean_absorptivity,
    :absorptivity => std  => :sd_absorptivity,
    nrow          => :n,
)
sort!(summary, [:region_of_measurement, :temperature_celsius])
println("── Summary by region and temperature ────────────────────────────────")
println(summary)
println()

# ── 3. Compare with pre-computed parameters export ───────────────────────────
precomputed = CSV.read(
    joinpath(RADIATION_DB, "export", "data", "parameters",
             "tableB_lizard_morphology_radiative_properties_DB.csv"),
    DataFrame,
)
# The parameters export uses "lizard01"… IDs mapped to M3…M14 in Smith et al.
# Compare means: our computation vs. the R pipeline output
for region in ["dorsal", "ventral"]
    ours = mean(filter(r -> r.region_of_measurement == region, individual_absorptivities).absorptivity)
    theirs = mean(filter(r -> r.region_of_measurement == region, precomputed).abs)
    println("$region: ours = $(round(ours; digits=4))   R pipeline = $(round(theirs; digits=4))   diff = $(round(ours - theirs; digits=4))")
end
println()

# ── 4. Build RadiationParameters with provenance tracking ────────────────────
#
# radiation_parameters detects the :spectral format automatically, computes
# per-individual absorptivities internally, and averages by region.
# Geometry areas are not available here so they will default with warnings.
println("── RadiationParameters for Pogona vitticeps ─────────────────────────")
radiation_params, provenance = radiation_parameters(
    "Pogona vitticeps",
    spectral_data,
    nothing,              # no geometry data — areas will default
)
println()
println("body_absorptivity_dorsal:  $(round(radiation_params.body_absorptivity_dorsal;  digits=4))")
println("body_absorptivity_ventral: $(round(radiation_params.body_absorptivity_ventral; digits=4))")
println("body_emissivity_dorsal:    $(radiation_params.body_emissivity_dorsal)")
println()
data_gaps(provenance; name="Pogona vitticeps — radiation")
