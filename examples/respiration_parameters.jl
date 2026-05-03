# Respiration and metabolism parameters from respirometry data
#
# Demonstrates respiration_parameters and metabolism_parameters for two species:
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

using BiophysicalParameters
using CSV
using DataFrames
using HeatExchange
using Statistics
using Unitful

const RESPIRATION_DB = joinpath(
    homedir(), "Dropbox", "Current Research Projects",
    "trait_database", "heat_budget_databases", "respirationDB",
)

celsius(t) = round(ustrip(uconvert(u"°C", t)); digits=2)
watts(w)   = round(ustrip(uconvert(u"W",  w)); digits=5)

# ── 1. Pogona vitticeps — open-circuit respirometry (Wild et al. 2023) ────────

pogona_raw = CSV.read(
    joinpath(RESPIRATION_DB, "data", "Wild_etal_2023", "data.csv"),
    DataFrame,
)
println("Wild et al. 2023: $(nrow(pogona_raw)) rows, $(length(unique(pogona_raw.Individual_ID))) individuals")
println("Columns: $(names(pogona_raw))")
println("Data format: $(respiration_data_format(pogona_raw))")
println()

# Per-individual RQ: CO₂ produced / O₂ consumed
pogona_raw.respiratory_quotient = pogona_raw.change_CO2 ./ pogona_raw.change_O2
println("── Per-individual respiratory quotients ─────────────────────────────")
println("  Mean RQ:   $(round(mean(pogona_raw.respiratory_quotient); digits=3))")
println("  Median RQ: $(round(median(pogona_raw.respiratory_quotient); digits=3))")
println("  SD RQ:     $(round(std(pogona_raw.respiratory_quotient); digits=3))")
println()

# Respiration parameters — RQ measured, oxygen_extraction_efficiency defaults
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

# Metabolism parameters — single temperature (33°C), Q10 not fittable.
# Ectotherm: no thermal_group, ambient Ta used as core_temperature.
# STP correction applied using bp column (hPa) and Ta = 33°C → factor ≈ 0.89.
pogona_metab_params, pogona_metab_prov = metabolism_parameters(
    "Pogona vitticeps",
    pogona_raw;
    respiratory_quotient = pogona_resp_params.respiratory_quotient,
)
pogona_core_tc  = celsius(pogona_metab_params.core_temperature)
pogona_heat_w   = watts(pogona_metab_params.metabolic_heat_flow)
println("── MetabolismParameters — Pogona vitticeps ───────────────────────────")
println("  core_temperature:    $pogona_core_tc °C  (= Ta, ectotherm)")
println("  metabolic_heat_flow: $pogona_heat_w W  (STP-corrected)")
println("  q10:                 $(pogona_metab_params.q10)  (GlobalDefault — single temperature)")
println()
data_gaps(pogona_metab_prov; name="Pogona vitticeps — metabolism")
println()

# ── 2. Myrmecobius fasciatus — metabolic rate vs temperature (Cooper & Withers 2002) ──

numbat_raw = CSV.read(
    joinpath(RESPIRATION_DB, "data", "Cooper_Withers_2002", "data.csv"),
    DataFrame,
)
println("Cooper & Withers 2002: $(nrow(numbat_raw)) rows")
println("Columns: $(names(numbat_raw))")
println("Data format: $(respiration_data_format(numbat_raw))")
println()
println(numbat_raw)
println()

# Respiration parameters — no CO₂ data in Cooper dataset, RQ defaults
numbat_resp_params, numbat_resp_prov = respiration_parameters(
    "Myrmecobius fasciatus",
    numbat_raw,
)
println("── RespirationParameters — Myrmecobius fasciatus ─────────────────────")
println("  respiratory_quotient: $(numbat_resp_params.respiratory_quotient)  (GlobalDefault — no CO₂ data)")
println()

# Metabolism parameters — minimum metabolic rate (BMR), T_b from data.csv,
# Q10 regressed against T_b. thermal_group = :marsupial documents the taxon class;
# T_b is already in the data so the group default is not needed for core_temperature,
# but it is recorded in provenance and would be used if T_b were absent.
# STP correction uses Ta (chamber gas temperature), not T_b.
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

# Note: the Q10 is fitted against T_b. Since T_b is nearly constant (~33.7–34.1°C)
# across Ta = 15–30°C and only rises at Ta = 32.5°C, the regression captures the
# metabolic increase near the upper thermal limit rather than a classic Q10 relationship.
# For endotherms, metabolic rate is primarily driven by the thermal gradient (T_b - T_a)
# for heat loss, not by T_b itself changing — this is a structural limitation of using
# Q10 as a parameter for homeotherms.

# For comparison: estimate_type = "mean"
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
pogona_metab_kleiber, _ = metabolism_parameters(
    "Pogona vitticeps",
    pogona_raw;
    respiratory_quotient     = pogona_resp_params.respiratory_quotient,
    oxygen_joules_conversion = HeatExchange.Kleiber1961(),
)
kleiber_heat_w = watts(pogona_metab_kleiber.metabolic_heat_flow)
println("── Pogona vitticeps: O₂ conversion model comparison ─────────────────")
println("  Typical()    (20.1 J/mL):   $pogona_heat_w W")
println("  Kleiber1961  (RQ-dependent): $kleiber_heat_w W  (RQ = $(round(pogona_resp_params.respiratory_quotient; digits=3)))")
