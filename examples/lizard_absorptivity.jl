# Lizard solar absorptivity example
#
# Computes per-individual solar-weighted absorptivity for Pogona vitticeps from
# raw spectral reflectance data (Smith et al. 2016), compares against the
# pre-computed values in the radiationDB parameters export, and builds a
# RadiationParameters struct with provenance tracking.
#
# Data source: traits.build .rds database loaded via RData.jl.
# join_contexts merges the method context (region_of_measurement) onto the traits
# table. pivot_traits_build_wide pivots to one row per wavelength measurement.

using BiophysicalParameters
using CSV
using DataFrames
using RData
using Statistics

const RADIATION_DB = joinpath(
    homedir(), "Dropbox", "Current Research Projects",
    "trait_database", "heat_budget_databases", "radiationDB",
)

# ── 1. Load traits.build database from .rds ──────────────────────────────────
raw_db      = RData.load(joinpath(RADIATION_DB, "export", "data", "current_DB",
                                   "Morphology_radiative_properties.rds"))
traits_long = join_contexts(DataFrame(raw_db["traits"]), DataFrame(raw_db["contexts"]))

println("Morphology_radiative_properties database: $(nrow(traits_long)) trait rows")
println("Taxa: $(unique(traits_long.taxon_name))")
println("Traits: $(unique(traits_long.trait_name))")
println()

# ── 2. Filter and pivot for Pogona vitticeps ──────────────────────────────────
#
# id_cols: each unique combination of (observation_id, repeat_measurements_id,
# region_of_measurement) is one wavelength measurement for one individual.
# repeat_measurements_id indexes the wavelength steps within each individual+region.

pogona_long = filter(r -> r.taxon_name == "Pogona vitticeps", traits_long)
pogona_wide = pivot_traits_build_wide(
    pogona_long,
    [:taxon_name, :observation_id, :repeat_measurements_id, :entity_type,
     :region_of_measurement],
)

println("Pogona vitticeps: $(nrow(pogona_wide)) rows")
println("Columns: $(names(pogona_wide))")
println("Wavelength range: $(minimum(skipmissing(pogona_wide[!, "wave_length(nm)"])))–$(maximum(skipmissing(pogona_wide[!, "wave_length(nm)"]))) nm")
println("Body temperatures: $(sort(unique(skipmissing(pogona_wide[!, "temperature_body(Cel)"])))) °C")
println("Regions: $(unique(pogona_wide.region_of_measurement))")
println()

# Rename traits.build column names to what per_individual_absorptivity expects.
# observation_id → Individual_ID so the function groups per individual correctly.
rename!(pogona_wide,
    "wave_length(nm)"               => "wave_length_nm",
    "temperature_body(Cel)"         => "temp_C",
    "reflectance({dimensionless})"  => "reflectance",
    "observation_id"                => "Individual_ID",
)

# ── 3. Compute per-individual solar-weighted absorptivities ──────────────────
#
# Smith et al. measured each lizard at two body temperatures (15 °C and 40 °C).
# Bearded dragons change colour with temperature — dorsal absorptivity is higher
# at 15 °C (dark) than at 40 °C (pale), while ventral is stable.
println("── Per-individual absorptivity ──────────────────────────────────────")
individual_absorptivities = per_individual_absorptivity(pogona_wide)
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

# ── 4. Compare with pre-computed parameters export ───────────────────────────
precomputed = CSV.read(
    joinpath(RADIATION_DB, "export", "data", "parameters",
             "tableB_lizard_morphology_radiative_properties_DB.csv"),
    DataFrame,
)
for region in ["dorsal", "ventral"]
    ours   = mean(filter(r -> r.region_of_measurement == region, individual_absorptivities).absorptivity)
    theirs = mean(filter(r -> r.region_of_measurement == region, precomputed).abs)
    println("$region: ours = $(round(ours; digits=4))   R pipeline = $(round(theirs; digits=4))   diff = $(round(ours - theirs; digits=4))")
end
println()

# ── 5. Build RadiationParameters with provenance tracking ────────────────────
#
# radiation_parameters detects the :spectral format automatically, computes
# per-individual absorptivities internally, and averages by region.
# Geometry areas are not available here so they will default with warnings.
println("── RadiationParameters for Pogona vitticeps ─────────────────────────")
radiation_params, provenance = radiation_parameters(
    "Pogona vitticeps",
    pogona_wide,
    nothing,              # no geometry data — areas will default
)
println()
println("body_absorptivity_dorsal:  $(round(radiation_params.body_absorptivity_dorsal;  digits=4))")
println("body_absorptivity_ventral: $(round(radiation_params.body_absorptivity_ventral; digits=4))")
println("body_emissivity_dorsal:    $(radiation_params.body_emissivity_dorsal)")
println()
data_gaps(provenance; name="Pogona vitticeps — radiation")
