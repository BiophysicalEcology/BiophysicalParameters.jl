# Respiration and metabolism parameters from respirometry data
#
# Demonstrates loading a traits.build database (.rds) via RData.jl, joining
# context variables (Ta, MR_estimate_type), pivoting to wide format, and feeding
# the result into respiration_parameters / metabolism_parameters.
#
# Two species:
#
# 1. Pogona vitticeps (Wild et al. 2023) — individual-level open-circuit respirometry
#    at a single temperature (33°C). Provides RQ from CO₂/O₂ ratios and metabolic
#    heat flow per individual. Ectotherm: ambient temperature used as core temperature.
#
# 2. Myrmecobius fasciatus — numbat (Cooper & Withers 2002) — population-level
#    metabolic rate across five ambient temperatures (15–32.5°C). Body temperature
#    (T_b) is included in the data, digitised from Fig. 1 of the paper: T_b is nearly
#    constant at ~33.7–34.1°C across 15–30°C ambient, rising to ~35.7°C at 32.5°C.
#    Q10 is regressed against T_b (body temperature drives metabolic rate in endotherms).

using BiophysicalAllometry
using BiophysicalParameters
using DataFrames
using HeatExchange
using RData          # activates TraitDataSources RData extension
using Statistics
using Unitful

celsius(t) = round(ustrip(uconvert(u"°C", t)); digits=2)
watts(w)   = round(ustrip(uconvert(u"W",  w)); digits=5)

# ── Load traits.build database ────────────────────────────────────────────────
# ENV["HEATBUDGETDB_PATH"] must point to the heat budget databases root directory.
# HeatBudgetDB{RespirationDomain} resolves to:
#   $HEATBUDGETDB_PATH/respirationDB/export/data/current_DB/Physiology_respirometry.rds

traits_long = gettraits(HeatBudgetDB{RespirationDomain}())

println("Physiology_respirometry database: $(nrow(traits_long)) trait rows")
println("Taxa: $(unique(traits_long.taxon_name))")
println("Traits: $(unique(traits_long.trait_name))")
println()

# ── 1. Pogona vitticeps — open-circuit respirometry (Wild et al. 2023) ────────
#
# Individual-level data: id columns are taxon_name, entity_type, observation_id,
# repeat_measurements_id, Ta. Mirrors the R Analysis.R pivot logic.

pogona_long = filter(r -> r.taxon_name == "Pogona vitticeps", traits_long)
pogona_raw = pivot_traits_build_wide(
    pogona_long,
    [:taxon_name, :entity_type, :observation_id, :repeat_measurements_id, :Ta],
)

println("Pogona vitticeps: $(nrow(pogona_raw)) rows, $(length(unique(pogona_raw.observation_id))) individuals")
println("Columns: $(names(pogona_raw))")
println("Data format: $(respiration_data_format(pogona_raw))")
println()

# Per-individual RQ: CO₂ produced / O₂ consumed
co2_col = findfirst(c -> occursin("change_co2", lowercase(string(c))), names(pogona_raw))
o2_col  = findfirst(c -> occursin("change_o2",  lowercase(string(c))), names(pogona_raw))
pogona_raw[!, :respiratory_quotient] = pogona_raw[!, names(pogona_raw)[co2_col]] ./
                                        pogona_raw[!, names(pogona_raw)[o2_col]]
println("── Per-individual respiratory quotients ─────────────────────────────")
println("  Mean RQ:   $(round(mean(skipmissing(pogona_raw.respiratory_quotient)); digits=3))")
println("  Median RQ: $(round(median(skipmissing(pogona_raw.respiratory_quotient)); digits=3))")
println("  SD RQ:     $(round(std(skipmissing(pogona_raw.respiratory_quotient)); digits=3))")
println()

pogona_resp_params, pogona_resp_prov = respiration_parameters(
    "Pogona vitticeps",
    pogona_raw,
)
println("── RespirationParameters — Pogona vitticeps ──────────────────────────")
println("  respiratory_quotient:         $(round(pogona_resp_params.respiratory_quotient; digits=3))")
println("  oxygen_extraction_efficiency: $(pogona_resp_params.oxygen_extraction_efficiency)  (default)")
println()
data_gaps(pogona_resp_prov; name="Pogona vitticeps — respiration")
println()

pogona_metab_params, pogona_metab_prov = metabolism_parameters(
    "Pogona vitticeps",
    pogona_raw;
    respiratory_quotient = pogona_resp_params.respiratory_quotient,
)
pogona_core_tc = celsius(pogona_metab_params.core_temperature)
pogona_heat_w  = watts(pogona_metab_params.metabolic_heat_flow)
println("── MetabolismParameters — Pogona vitticeps ───────────────────────────")
println("  core_temperature:    $pogona_core_tc °C  (= Ta, ectotherm)")
println("  metabolic_heat_flow: $pogona_heat_w W  (STP-corrected)")
println("  q10:                 $(pogona_metab_params.q10)  (GlobalDefault — single temperature)")
println()
data_gaps(pogona_metab_prov; name="Pogona vitticeps — metabolism")
println()

# ── 2. Myrmecobius fasciatus — metabolic rate vs temperature (Cooper & Withers 2002) ──
#
# Population-level data: id columns are taxon_name, replicates, Ta, MR_estimate_type.
# Tb (body temperature) is a trait in the database; after pivoting it appears as the
# column temperature_body_resp(Cel), which the builders detect automatically.

numbat_long = filter(r -> r.taxon_name == "Myrmecobius fasciatus", traits_long)
numbat_raw = pivot_traits_build_wide(
    numbat_long,
    [:taxon_name, :replicates, :Ta, :MR_estimate_type],
)

println("Cooper & Withers 2002: $(nrow(numbat_raw)) rows")
println("Columns: $(names(numbat_raw))")
println("Data format: $(respiration_data_format(numbat_raw))")
println()
println(numbat_raw)
println()

numbat_resp_params, _ = respiration_parameters("Myrmecobius fasciatus", numbat_raw)
println("── RespirationParameters — Myrmecobius fasciatus ─────────────────────")
println("  respiratory_quotient: $(numbat_resp_params.respiratory_quotient)  (GlobalDefault — no CO₂ data)")
println()

numbat_metab_params, numbat_metab_prov = metabolism_parameters(
    "Myrmecobius fasciatus",
    numbat_raw;
    estimate_type = "min",
    thermal_group = :marsupial,
)
numbat_core_tc = celsius(numbat_metab_params.core_temperature)
numbat_heat_w  = watts(numbat_metab_params.metabolic_heat_flow)
println("── MetabolismParameters — Myrmecobius fasciatus (estimate_type = min) ─")
println("  core_temperature:    $numbat_core_tc °C  (mean T_b from data)")
println("  metabolic_heat_flow: $numbat_heat_w W  (at mean T_b, STP-corrected)")
println("  q10:                 $(round(numbat_metab_params.q10; digits=3))  (fitted from T_b vs MR_min)")
println()
data_gaps(numbat_metab_prov; name="Myrmecobius fasciatus — metabolism (min)")
println()

numbat_metab_mean, _ = metabolism_parameters(
    "Myrmecobius fasciatus",
    numbat_raw;
    estimate_type = "mean",
    thermal_group = :marsupial,
)
numbat_heat_mean_w = watts(numbat_metab_mean.metabolic_heat_flow)
println("── MetabolismParameters — Myrmecobius fasciatus (estimate_type = mean) ─")
println("  metabolic_heat_flow: $numbat_heat_mean_w W")
println("  q10:                 $(round(numbat_metab_mean.q10; digits=3))")
println()

# ── 3. Endotherm defaults when no T_b data is available ───────────────────────
println("── Endotherm default core temperatures ──────────────────────────────")
for group in (:monotreme, :marsupial, :eutherian, :bird)
    t_default = endotherm_default_core_temperature(group)
    t_celsius = round(ustrip(uconvert(u"°C", t_default)); digits=1)
    println("  :$group → $t_celsius °C")
end
println()

# ── 4. Effect of O₂ conversion model on metabolic heat flow ──────────────────
pogona_metab_kleiber1961, _ = metabolism_parameters(
    "Pogona vitticeps",
    pogona_raw;
    respiratory_quotient     = pogona_resp_params.respiratory_quotient,
    oxygen_joules_conversion = HeatExchange.Kleiber1961(),
)
kleiber1961_heat_w = watts(pogona_metab_kleiber1961.metabolic_heat_flow)
println("── Pogona vitticeps: O₂ conversion model comparison ─────────────────")
println("  Typical()    (20.1 J/mL):   $pogona_heat_w W")
println("  Kleiber1961  (RQ-dependent): $kleiber1961_heat_w W  (RQ = $(round(pogona_resp_params.respiratory_quotient; digits=3)))")
println()

# ── 5. Measured vs AndrewsPough2 allometric prediction — Pogona vitticeps ────
mass_col_idx           = findfirst(c -> occursin("mass_in", lowercase(string(c))), names(pogona_raw))
mean_pogona_mass_grams = mean(skipmissing(pogona_raw[!, names(pogona_raw)[mass_col_idx]]))
pogona_mass            = mean_pogona_mass_grams * u"g"

andrews_pough_standard_w = watts(allometric(
    StandardMetabolicRate(), Squamate(), pogona_mass, pogona_metab_params.core_temperature,
))
andrews_pough_resting_w = watts(allometric(
    StandardMetabolicRate(), Squamate(), pogona_mass, pogona_metab_params.core_temperature;
    metabolic_state = 1.0,
))
ratio_to_standard = round(pogona_heat_w / andrews_pough_standard_w; digits=2)
println("── Pogona vitticeps: measured vs AndrewsPough2 allometric prediction ─")
println("  Mean body mass:              $(round(mean_pogona_mass_grams; digits=1)) g")
println("  Core temperature:            $pogona_core_tc °C")
println("  Measured (respirometry):     $pogona_heat_w W  (STP-corrected, Typical)")
println("  AndrewsPough2 standard:      $andrews_pough_standard_w W  (metabolic_state = 0)")
println("  AndrewsPough2 resting:       $andrews_pough_resting_w W  (metabolic_state = 1)")
println("  Measured / standard:         $ratio_to_standard")
println()

# ── 6. Measured vs marsupial BMR prediction — Myrmecobius fasciatus ──────────
mass_col_numbat_idx = findfirst(c -> occursin("mass_in", lowercase(string(c))), names(numbat_raw))
mean_numbat_mass_g  = first(skipmissing(numbat_raw[!, names(numbat_raw)[mass_col_numbat_idx]]))
numbat_mass         = mean_numbat_mass_g * u"g"

marsupial_bmr_w  = watts(allometric(BasalMetabolicRate(), Marsupial(), numbat_mass))
eutherian_bmr_w  = watts(allometric(BasalMetabolicRate(), EutherianMammal(), numbat_mass))
ratio_numbat_min  = round(numbat_heat_w       / marsupial_bmr_w; digits=2)
ratio_numbat_mean = round(numbat_heat_mean_w  / marsupial_bmr_w; digits=2)
println("── Myrmecobius fasciatus: measured vs allometric BMR prediction ──────")
println("  Mean body mass:                    $(round(mean_numbat_mass_g; digits=1)) g")
println("  Core temperature:                  $numbat_core_tc °C  (mean T_b)")
println("  Marsupial BMR (Dawson & Hulbert):   $marsupial_bmr_w W")
println("  Eutherian BMR (Kleiber / Schmidt-Nielsen): $eutherian_bmr_w W")
println("  Measured min  (respirometry):      $numbat_heat_w W")
println("  Measured mean (respirometry):      $numbat_heat_mean_w W")
println("  Measured min  / marsupial BMR:     $ratio_numbat_min")
println("  Measured mean / marsupial BMR:     $ratio_numbat_mean")
