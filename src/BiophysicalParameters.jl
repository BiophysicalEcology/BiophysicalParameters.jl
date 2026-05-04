module BiophysicalParameters

using Arrow
using CSV
using DataFrames
using Dates
using Distributions
using HeatExchange
using Interpolations
using JSON3
using SolarRadiation
using Statistics
using TraitDataSources
using Unitful

# ── Exports ───────────────────────────────────────────────────────────────────

# Provenance
export DataSource, Measured, Allometric, PhylogeneticMean, TaxonDefault, GlobalDefault
export FieldProvenance, ParameterProvenance
export completeness_score, data_gaps

# Model trait resolution
export resolve_insulation_type

# Parameter builders — each returns (params, provenance)
export radiation_parameters
export per_individual_absorptivity
export radiation_data_format
export location_solar_spectrum
export respiration_parameters
export metabolism_parameters
export respiration_data_format
export endotherm_default_core_temperature
export ENDOTHERM_DEFAULT_CORE_TEMPERATURES
# export shape_parameters          # Phase 1 step 4
# export insulation_parameters     # Phase 1 step 5
# export respiration_parameters    # Phase 1 step 6
# export metabolism_parameters     # Phase 1 step 6
# export convection_parameters     # Phase 1 step 4
# export conduction_parameters_external  # Phase 1 step 4
# export conduction_parameters_internal  # Phase 1 step 4
# export evaporation_parameters    # Phase 1 step 8
# export hydraulic_parameters      # Phase 1 step 4

# Allometry fitting
# export fit_allometry             # Phase 1 step 7

# Inverse fitting
# export fit_skin_wetness          # Phase 1 step 8
# export fit_absorptivity          # Phase 1 step 8

# Organism builder
# export build_organism            # Phase 1 step 9

# Solar absorptivity utility (public for testing/reuse)
export solar_weighted_absorptivity

# traits.build integration
# join_contexts lives in TraitDataSources and is re-exported here for backwards compatibility
export join_contexts, pivot_traits_build_wide

# TraitDataSources re-exports — `using BiophysicalParameters` is sufficient for TraitDataSources API
export TraitDataSource, TraitDomain
export GeometryDomain, InsulationDomain, RadiationDomain, RespirationDomain
export HeatBudgetDB, TraitsBuildFile, CompiledParametersDB
export gettraits, traitnames, traitpath

# Compiled parameters database
export ParameterSummary, ParameterRecord
export compile_parameters_database, parameters_table, load_parameters
export distribution

# ── Includes ──────────────────────────────────────────────────────────────────

include("provenance.jl")
include("model_traits.jl")
include("radiation.jl")
include("respiration.jl")
include("traits_build.jl")
include("compiled_db.jl")

end
