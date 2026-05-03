# ─── Format detection ─────────────────────────────────────────────────────────

"""
    respiration_data_format(df) → Symbol

Detect the format of a respiration DataFrame:

- `:open_circuit` — individual-level flow-through respirometry with O₂ and CO₂ change
  rates and chamber flow rate (e.g. Wild et al. 2023). Provides RQ and metabolic heat flow.
- `:metabolic_rate` — O₂ consumption rate, possibly at multiple ambient temperatures
  (e.g. Cooper & Withers 2002). Provides metabolic heat flow; Q10 if ≥ 3 temperatures.
"""
function respiration_data_format(df::DataFrame)
    column_names_lower = lowercase.(string.(names(df)))
    has_oxygen_change       = any(c -> occursin("change_o2",  c), column_names_lower)
    has_carbon_dioxide_change = any(c -> occursin("change_co2", c), column_names_lower)
    has_chamber_flow_rate   = any(c -> occursin("flow_rate",   c), column_names_lower)
    has_metabolic_rate_o2   = any(c -> occursin("mr_o2", c) || occursin("mr(o2", c), column_names_lower)
    if has_oxygen_change && has_carbon_dioxide_change && has_chamber_flow_rate
        return :open_circuit
    elseif has_metabolic_rate_o2 || has_oxygen_change
        return :metabolic_rate
    else
        error("Cannot determine respiration data format from columns: " * join(names(df), ", "))
    end
end

# ─── Column finders ───────────────────────────────────────────────────────────

function _oxygen_rate_column(df)
    index = findfirst(c -> occursin("change_o2", lowercase(string(c))), names(df))
    index === nothing && return nothing
    Symbol(names(df)[index])
end

function _carbon_dioxide_rate_column(df)
    index = findfirst(c -> occursin("change_co2", lowercase(string(c))), names(df))
    index === nothing && return nothing
    Symbol(names(df)[index])
end

function _metabolic_rate_o2_column(df)
    index = findfirst(
        c -> let s = lowercase(string(c)); occursin("mr_o2", s) || occursin("mr(o2", s); end,
        names(df),
    )
    index === nothing && return nothing
    Symbol(names(df)[index])
end

# Returns true when the O₂ column name indicates mass-specific units (mL O₂ g⁻¹ h⁻¹).
# Such columns must be multiplied by body mass before energy conversion.
function _is_mass_specific_o2_column(column_name)
    name_lower = lowercase(string(column_name))
    occursin("_g_h", name_lower) || occursin("/(g", name_lower) || occursin("/g/h", name_lower)
end

function _mass_column(df)
    candidate = _detect_column(
        df,
        [:Mean_Resp_mass, :mass_in, :mass_in_resp, :body_mass, :mass],
        nothing,
    )
    candidate !== nothing && return candidate
    index = findfirst(c -> occursin("mass", lowercase(string(c))), names(df))
    index === nothing ? nothing : Symbol(names(df)[index])
end

# Body temperature (T_b) takes priority over ambient temperature (T_a) for
# core_temperature and Q10 fitting. For endotherms, add a T_b column to the
# DataFrame or pass reference_temperature directly.
function _body_temperature_column(df)
    # Check known symbol names first, then fall back to pattern matching for export-format
    # column names like "temperature_body_resp(Cel)"
    candidate = _detect_column(df, [:T_b, :Tb, :body_temperature, :body_temperature_celsius], nothing)
    candidate !== nothing && return candidate
    index = findfirst(c -> occursin("temperature_body", lowercase(string(c))), names(df))
    index === nothing ? nothing : Symbol(names(df)[index])
end

function _ambient_temperature_column(df)
    _detect_column(df, [:Ta, :ta, :ambient_temperature, :ambient_temperature_celsius], nothing)
end

function _estimate_type_column(df)
    _detect_column(df, [:MR_estimate_type, :mr_estimate_type], nothing)
end

# Infer O₂ volume time-unit from column name: mL/hr for Cooper-style, mL/min otherwise.
# Handles both total (ml_h) and mass-specific (_g_h) per-hour variants.
function _infer_oxygen_time_unit(column_name)
    name_lower = lowercase(string(column_name))
    is_per_hour = (
        occursin("_g_h", name_lower) ||
        occursin("ml_h", name_lower) ||
        occursin("/h)", name_lower)  ||
        occursin("{o2}/h", name_lower)
    )
    is_per_hour ? u"mL/hr" : u"mL/minute"
end

function _barometric_pressure_column(df)
    _detect_column(df, [:bp, :barometric_pressure, :atmospheric_pressure], nothing) |>
    col -> if col === nothing
        n = findfirst(c -> occursin("pa_chamber", lowercase(string(c))), names(df))
        n === nothing ? nothing : Symbol(names(df)[n])
    else
        col
    end
end

# Convert a measured O₂ volume rate from experimental conditions to standard conditions
# (STP: 0°C, 101325 Pa). The correction mirrors the R formula in
# Physiology_respirometry_database_functions.R:
#   CorVol = Vol * (P_exp / P_STP) * (T_STP / T_exp)
# Note: the R code uses 101.325 in the denominator, which is correct only when bp is in
# kPa — the Wild et al. 2023 data records bp in hPa (labeled Pa in the original metadata,
# corrected to hPa). Always attach Unitful units so mismatches are caught at runtime.
#
# If no barometric pressure is available the correction still applies the temperature term,
# assuming sea-level pressure (P_exp ≈ P_STP, pressure ratio ≈ 1).
function _stp_correction_factor(
    experimental_temperature::Unitful.Temperature,
    barometric_pressure::Union{Unitful.Pressure, Nothing},
)
    T_STP = 273.15u"K"
    P_STP = 101325.0u"Pa"
    T_experimental = uconvert(u"K", experimental_temperature)
    P_experimental = isnothing(barometric_pressure) ? P_STP : uconvert(u"Pa", barometric_pressure)
    return ustrip(NoUnits, (P_experimental / P_STP) * (T_STP / T_experimental))
end

# Read barometric pressure from a DataFrame column, attaching Unitful units.
# Performs a sanity check: values in 900–1100 with a "Pa" annotation are almost
# certainly hPa (standard atmospheric pressure in hPa is ~1013). A warning is
# emitted if this mismatch is detected, so the unit error can be traced to its source.
function _read_barometric_pressure(df, pressure_column)
    isnothing(pressure_column) && return nothing
    raw_values   = Float64.(df[!, pressure_column])
    column_label = string(pressure_column)
    mean_value   = mean(raw_values)
    # Unit annotation in export-format column names like "Pa_chamber_resp(Pa)"
    annotated_pa = occursin("(pa)", lowercase(column_label)) && !occursin("(hpa)", lowercase(column_label))
    if annotated_pa && 900.0 < mean_value < 1100.0
        @warn "Column $column_label is annotated as Pa but values (mean=$mean_value) are " *
              "consistent with hPa. Using hPa. Check metadata.yml unit_in for this trait."
        return mean_value * u"hPa"
    end
    # Plain column names (e.g. :bp) with values in the hPa range
    if !occursin("(", column_label) && 900.0 < mean_value < 1100.0
        return mean_value * u"hPa"
    end
    return mean_value * u"Pa"
end

# ─── Q10 fitting ─────────────────────────────────────────────────────────────

# Fit Q10 via log₁₀-linear regression: log₁₀(MR) = a + b·T → Q10 = 10^(10b).
function _fit_q10(temperatures_celsius::AbstractVector, metabolic_rates::AbstractVector)
    temperatures      = Float64.(temperatures_celsius)
    log_rates         = log10.(Float64.(metabolic_rates))
    n_observations    = length(temperatures)
    mean_temperature  = sum(temperatures) / n_observations
    mean_log_rate     = sum(log_rates)    / n_observations
    temperature_sensitivity = (
        sum((temperatures .- mean_temperature) .* (log_rates .- mean_log_rate)) /
        sum((temperatures .- mean_temperature) .^ 2)
    )
    return 10.0 ^ (10.0 * temperature_sensitivity)
end

# ─── Temperature unit helpers ─────────────────────────────────────────────────

_to_kelvin(temperature::Real)                = uconvert(u"K", Float64(temperature) * u"°C")
_to_kelvin(temperature::Unitful.Temperature) = uconvert(u"K", temperature)

_to_celsius_value(temperature::Real)                = Float64(temperature)
_to_celsius_value(temperature::Unitful.Temperature) = ustrip(u"°C", uconvert(u"°C", temperature))

# ─── Default constants ────────────────────────────────────────────────────────

const DEFAULT_RESPIRATORY_QUOTIENT         = 0.8
const DEFAULT_OXYGEN_EXTRACTION_EFFICIENCY = 0.2
const DEFAULT_CORE_TEMPERATURE             = uconvert(u"K", 37.0u"°C")
const DEFAULT_METABOLIC_HEAT_FLOW          = 0.0u"W"
const DEFAULT_Q10                          = 2.0

# Default core temperatures by endotherm group (TaxonDefault provenance).
# Monotremes have the lowest mammalian body temperatures; marsupials are intermediate;
# eutherians typically ~37°C. Birds vary substantially by order — the value below is a
# broad avian average; passerines average ~42°C, ratites ~38°C, penguins ~38°C.
# Source: Withers (1992) Comparative Animal Physiology; McNab (2002) The Physiological
# Ecology of Vertebrates.
const ENDOTHERM_DEFAULT_CORE_TEMPERATURES = (
    monotreme = uconvert(u"K", 32.0u"°C"),
    marsupial = uconvert(u"K", 35.5u"°C"),
    eutherian = uconvert(u"K", 37.0u"°C"),
    bird      = uconvert(u"K", 41.0u"°C"),
)

"""
    endotherm_default_core_temperature(thermal_group::Symbol) → Quantity{K}

Return the typical core temperature for an endotherm taxonomic group.

Valid groups: `:monotreme` (~32°C), `:marsupial` (~35.5°C), `:eutherian` (~37°C),
`:bird` (~41°C). Raises an error for unknown groups.
"""
function endotherm_default_core_temperature(thermal_group::Symbol)
    haskey(ENDOTHERM_DEFAULT_CORE_TEMPERATURES, thermal_group) ||
        error("Unknown thermal group :$thermal_group. Valid groups: " *
              join(keys(ENDOTHERM_DEFAULT_CORE_TEMPERATURES), ", "))
    ENDOTHERM_DEFAULT_CORE_TEMPERATURES[thermal_group]
end

# ─── Builder: RespirationParameters ──────────────────────────────────────────

"""
    respiration_parameters(taxon, df_respiration; kwargs...)
    → (RespirationParameters, ParameterProvenance)

Build a `RespirationParameters` struct from respirometry data.

## Parameters estimated from data

| Parameter                      | Requires                                  |
|--------------------------------|-------------------------------------------|
| `respiratory_quotient`         | open-circuit data with CO₂ and O₂ columns |
| `oxygen_extraction_efficiency` | ventilation rate (not in standard chamber data) |

`oxygen_extraction_efficiency` is the physiological fraction of inspired O₂ extracted by
the lungs per breath. It cannot be derived from flow-through chamber measurements (which
give whole-animal O₂ consumption, not per-breath extraction). The NicheMapR default of
0.2 (20%) is used when no ventilation rate data is available.

## Data formats (detected automatically via `respiration_data_format`)
- `:open_circuit` — Wild et al. 2023 style: `change_O2`, `change_CO2`, `flow_rate`, `Ta`
- `:metabolic_rate` — Cooper & Withers 2002 style: `MR_O2_ml_h`, `Ta`, `MR_estimate_type`

## Keyword arguments
- `estimate_type` — which metabolic rate rows to use when a `MR_estimate_type` column
  is present. Default `"min"` (minimum rate ≈ BMR). Pass `"mean"` for mean rate.
- `O2conversion` — `OxygenJoulesConversion` model (default: `HeatExchange.Typical()`).
"""
function respiration_parameters(
    taxon,
    df_respiration;
    estimate_type::Union{String,Nothing}                      = nothing,
    oxygen_joules_conversion::HeatExchange.OxygenJoulesConversion = HeatExchange.Typical(),
)
    provenance = Dict{Symbol,FieldProvenance}()

    if isnothing(df_respiration) || isempty(df_respiration)
        for field in (:respiratory_quotient, :oxygen_extraction_efficiency, :pant,
                      :exhaled_temperature_offset, :exhaled_relative_humidity, :mouth_fraction)
            provenance[field] = FieldProvenance(GlobalDefault)
        end
        return HeatExchange.RespirationParameters(), NamedTuple(provenance)
    end

    df        = _filter_by_taxon(df_respiration, taxon)
    data_format = respiration_data_format(df)
    taxon_str = string(taxon)

    # ── Respiratory quotient ─────────────────────────────────────────────────
    respiratory_quotient = if data_format == :open_circuit
        oxygen_column          = _oxygen_rate_column(df)
        carbon_dioxide_column  = _carbon_dioxide_rate_column(df)
        rq_values = Float64.(df[!, carbon_dioxide_column]) ./ Float64.(df[!, oxygen_column])
        valid_rq_values = filter(isfinite, rq_values)
        if !isempty(valid_rq_values)
            provenance[:respiratory_quotient] = FieldProvenance(Measured, length(valid_rq_values), taxon_str)
            mean(valid_rq_values)
        else
            _warn_fallback(:respiratory_quotient, GlobalDefault, "all CO₂/O₂ ratios were non-finite")
            provenance[:respiratory_quotient] = FieldProvenance(GlobalDefault)
            DEFAULT_RESPIRATORY_QUOTIENT
        end
    else
        _warn_fallback(:respiratory_quotient, GlobalDefault,
            "no CO₂ data in supplied DataFrame; using default RQ = $DEFAULT_RESPIRATORY_QUOTIENT")
        provenance[:respiratory_quotient] = FieldProvenance(GlobalDefault)
        DEFAULT_RESPIRATORY_QUOTIENT
    end

    # ── Oxygen extraction efficiency — requires ventilation rate ─────────────
    _warn_fallback(:oxygen_extraction_efficiency, GlobalDefault,
        "physiological O₂ extraction efficiency requires ventilation rate data; " *
        "using NicheMapR default of $DEFAULT_OXYGEN_EXTRACTION_EFFICIENCY")
    provenance[:oxygen_extraction_efficiency] = FieldProvenance(GlobalDefault)

    # ── Remaining fields — literature defaults ───────────────────────────────
    for field in (:pant, :exhaled_temperature_offset, :exhaled_relative_humidity, :mouth_fraction)
        provenance[field] = FieldProvenance(GlobalDefault)
    end

    params = HeatExchange.RespirationParameters(; respiratory_quotient)
    return params, NamedTuple(provenance)
end

# ─── Builder: MetabolismParameters ───────────────────────────────────────────

"""
    metabolism_parameters(taxon, df_respiration; kwargs...)
    → (MetabolismParameters, ParameterProvenance)

Build a `MetabolismParameters` struct from respirometry data.

## Parameters estimated from data

| Parameter             | Source                                      | Data required           |
|-----------------------|---------------------------------------------|-------------------------|
| `metabolic_heat_flow` | O₂ consumed → W via `O2_to_Joules`          | any respirometry format |
| `core_temperature`    | measured body temperature, then ambient `Ta` | `T_b` or `Ta` column   |
| `q10`                 | log-linear regression of MR vs temperature   | ≥ 3 distinct temperatures |

`core_temperature` is taken from a `T_b`/`Tb`/`body_temperature` column if present,
falling back to `Ta`. For endotherms, body temperature differs from ambient temperature
and the paper containing the metabolic data often reports both (e.g. Cooper & Withers 2002
Fig. 1 shows T_b ≈ 33.7–34.1°C across 15–30°C ambient, rising to ~35.7°C at 32.5°C).
Add a `T_b` column to the DataFrame or supply `reference_temperature` directly.

Q10 is regressed against body temperature when a `T_b` column is present, and against
ambient temperature otherwise. For endotherms with near-constant T_b across the measured
temperature range, fewer than 3 distinct T_b values will prevent fitting.

## Keyword arguments
- `estimate_type` — filter rows by `MR_estimate_type` column when present.
  Default `"min"` (minimum rate ≈ standard/basal metabolic rate).
- `reference_temperature` — override `core_temperature`. Accepts `Real` (°C) or
  `Unitful.Temperature`. When supplied and a temperature column exists, the metabolic
  heat flow is taken from the row(s) nearest this temperature.
- `oxygen_joules_conversion` — `OxygenJoulesConversion` model (default: `HeatExchange.Typical()`).
  Use `HeatExchange.Kleiber1961()` for RQ-dependent energy equivalents.
- `respiratory_quotient` — used in the O₂→W conversion (default: $DEFAULT_RESPIRATORY_QUOTIENT).
  Pass the value from `respiration_parameters` for consistency.
- `thermal_group` — endotherm taxonomic group for `core_temperature` fallback when neither
  a `T_b` column nor `reference_temperature` is available. Valid values: `:monotreme`,
  `:marsupial`, `:eutherian`, `:bird`. When supplied, the group default temperature is used
  with `TaxonDefault` provenance rather than falling through to ambient temperature.
  Ectotherms should leave this `nothing` (ambient temperature is the appropriate fallback).
"""
function metabolism_parameters(
    taxon,
    df_respiration;
    estimate_type::Union{String,Nothing}                          = nothing,
    reference_temperature                                          = nothing,
    oxygen_joules_conversion::HeatExchange.OxygenJoulesConversion = HeatExchange.Typical(),
    respiratory_quotient::Real                                     = DEFAULT_RESPIRATORY_QUOTIENT,
    thermal_group::Union{Symbol,Nothing}                           = nothing,
)
    provenance = Dict{Symbol,FieldProvenance}()

    if isnothing(df_respiration) || isempty(df_respiration)
        for field in (:core_temperature, :metabolic_heat_flow, :q10)
            provenance[field] = FieldProvenance(GlobalDefault)
        end
        return HeatExchange.MetabolismParameters(), NamedTuple(provenance)
    end

    df        = _filter_by_taxon(df_respiration, taxon)
    data_format = respiration_data_format(df)
    taxon_str = string(taxon)

    # ── Filter by estimate type ──────────────────────────────────────────────
    estimate_type_col = _estimate_type_column(df)
    if !isnothing(estimate_type_col)
        chosen_estimate_type = something(estimate_type, "min")
        filtered_df = filter(r -> r[estimate_type_col] == chosen_estimate_type, df)
        df = isempty(filtered_df) ? df : filtered_df
    end

    # ── Core temperature ─────────────────────────────────────────────────────
    # Priority: reference_temperature > T_b column > thermal_group default > Ta column.
    # Ambient temperature (Ta) is appropriate for ectotherms but not endotherms — supply
    # thermal_group for endotherms without T_b measurements.
    body_temperature_col    = _body_temperature_column(df)
    ambient_temperature_col = _ambient_temperature_column(df)
    measurement_temp_col    = something(body_temperature_col, ambient_temperature_col, nothing)

    core_temperature, core_temperature_provenance = if !isnothing(reference_temperature)
        _to_kelvin(reference_temperature), FieldProvenance(Measured, 0, taxon_str)
    elseif !isnothing(body_temperature_col)
        mean_tb_celsius = mean(Float64.(df[!, body_temperature_col]))
        uconvert(u"K", mean_tb_celsius * u"°C"), FieldProvenance(Measured, nrow(df), taxon_str)
    elseif !isnothing(thermal_group)
        endotherm_default_core_temperature(thermal_group),
        FieldProvenance(TaxonDefault, 0, string(thermal_group))
    elseif !isnothing(ambient_temperature_col)
        mean_ta_celsius = mean(Float64.(df[!, ambient_temperature_col]))
        uconvert(u"K", mean_ta_celsius * u"°C"), FieldProvenance(Measured, nrow(df), taxon_str)
    else
        _warn_fallback(:core_temperature, GlobalDefault,
            "no T_b column, no thermal_group, and no Ta column found; " *
            "for endotherms pass thermal_group = :marsupial / :eutherian / :bird / :monotreme")
        DEFAULT_CORE_TEMPERATURE, FieldProvenance(GlobalDefault)
    end
    provenance[:core_temperature] = core_temperature_provenance

    reference_celsius = _to_celsius_value(core_temperature)

    # ── Metabolic heat flow ──────────────────────────────────────────────────
    metabolic_heat_flow, metabolic_heat_flow_provenance = _metabolic_heat_flow_from_data(
        df, data_format, taxon_str,
        measurement_temp_col, _to_celsius_value(core_temperature), reference_temperature,
        ambient_temperature_col,
        oxygen_joules_conversion, respiratory_quotient,
    )
    provenance[:metabolic_heat_flow] = metabolic_heat_flow_provenance

    # ── Q10 — only fittable from ≥ 3 distinct temperatures ───────────────────
    fitted_q10, q10_provenance = _fit_q10_from_data(
        df, data_format, taxon_str, ambient_temperature_col,
    )
    provenance[:q10] = q10_provenance

    params = HeatExchange.MetabolismParameters(;
        core_temperature    = core_temperature,
        metabolic_heat_flow = metabolic_heat_flow,
        q10                 = fitted_q10,
    )
    return params, NamedTuple(provenance)
end

# ─── Internal helpers ─────────────────────────────────────────────────────────

function _filter_by_taxon(df, taxon)
    taxon_column = _detect_column(df, [:taxon_name, :Species], nothing)
    isnothing(taxon_column) && return df
    taxon_str = string(taxon)
    filtered = filter(r -> r[taxon_column] == taxon_str, df)
    isempty(filtered) ? df : filtered
end

function _metabolic_heat_flow_from_data(
    df, data_format, taxon_str,
    measurement_temperature_column, reference_celsius, reference_temperature,
    ambient_temperature_column,
    oxygen_joules_conversion, respiratory_quotient,
)
    # Narrow to rows near the reference temperature when one is provided.
    # Row selection uses the measurement temperature (T_b or Ta); the STP correction
    # always uses ambient temperature (Ta) since that is the gas temperature in the chamber.
    selected_df = if !isnothing(reference_temperature) && !isnothing(measurement_temperature_column)
        nearest_index = argmin(abs.(Float64.(df[!, measurement_temperature_column]) .- reference_celsius))
        df[nearest_index:nearest_index, :]
    else
        df
    end

    oxygen_column = data_format == :open_circuit ?
        _oxygen_rate_column(selected_df) : _metabolic_rate_o2_column(selected_df)
    isnothing(oxygen_column) && begin
        _warn_fallback(:metabolic_heat_flow, GlobalDefault, "no O₂ rate column found")
        return DEFAULT_METABOLIC_HEAT_FLOW, FieldProvenance(GlobalDefault)
    end

    oxygen_time_unit        = _infer_oxygen_time_unit(oxygen_column)
    mean_oxygen_rate        = mean(Float64.(selected_df[!, oxygen_column]))

    # When the O₂ column is mass-specific (mL O₂ g⁻¹ h⁻¹), multiply by mean body mass to
    # obtain total O₂ consumption (mL O₂ h⁻¹) before applying the energy conversion.
    mean_oxygen_consumption = if _is_mass_specific_o2_column(oxygen_column)
        mass_col = _mass_column(selected_df)
        if !isnothing(mass_col)
            mean_mass_grams = mean(Float64.(selected_df[!, mass_col]))
            (mean_oxygen_rate * mean_mass_grams) * oxygen_time_unit
        else
            @warn "Mass-specific O₂ column $(oxygen_column) detected but no mass column found; " *
                  "values will be treated as total (not per-gram). Add a mass column to correct this."
            mean_oxygen_rate * oxygen_time_unit
        end
    else
        mean_oxygen_rate * oxygen_time_unit
    end

    # STP correction: O₂ volumes at experimental conditions → standard conditions
    # (0°C, 101325 Pa). Factor = (P_exp/P_STP) × (T_STP/T_exp). At 33°C ~0.89.
    # Always use ambient temperature (Ta) for T_exp — it is the gas temperature in the
    # chamber, not the animal's body temperature.
    barometric_pressure_column = _barometric_pressure_column(selected_df)
    barometric_pressure        = _read_barometric_pressure(selected_df, barometric_pressure_column)
    chamber_temperature_celsius = if !isnothing(ambient_temperature_column)
        mean(Float64.(selected_df[!, ambient_temperature_column]))
    else
        ustrip(u"°C", uconvert(u"°C", DEFAULT_CORE_TEMPERATURE))
    end
    stp_correction = _stp_correction_factor(chamber_temperature_celsius * u"°C", barometric_pressure)
    mean_oxygen_consumption_stp = mean_oxygen_consumption * stp_correction

    metabolic_heat_flow = uconvert(u"W",
        HeatExchange.O2_to_Joules(oxygen_joules_conversion, mean_oxygen_consumption_stp, respiratory_quotient))
    return metabolic_heat_flow, FieldProvenance(Measured, nrow(selected_df), taxon_str)
end

function _fit_q10_from_data(df, data_format, taxon_str, ambient_temperature_col)
    # Q10 should be regressed against body temperature when available.
    # For endotherms with near-constant T_b (e.g. numbats ~34°C across 15–30°C ambient),
    # fewer than 3 distinct T_b values will prevent fitting; add a T_b column to supply the data.
    body_temperature_col = _body_temperature_column(df)
    temperature_column   = something(body_temperature_col, ambient_temperature_col, nothing)
    temperature_label    = !isnothing(body_temperature_col) ? "T_b" : "Ta"

    isnothing(temperature_column) && begin
        _warn_fallback(:q10, GlobalDefault, "no temperature column; Q10 cannot be fitted")
        return DEFAULT_Q10, FieldProvenance(GlobalDefault)
    end
    temperatures = Float64.(df[!, temperature_column])
    n_distinct_temperatures = length(unique(temperatures))
    n_distinct_temperatures < 3 && begin
        _warn_fallback(:q10, GlobalDefault,
            "fewer than 3 distinct $temperature_label values ($n_distinct_temperatures found); " *
            "Q10 cannot be fitted — for endotherms add a T_b column to the DataFrame")
        return DEFAULT_Q10, FieldProvenance(GlobalDefault)
    end
    oxygen_column = data_format == :open_circuit ?
        _oxygen_rate_column(df) : _metabolic_rate_o2_column(df)
    isnothing(oxygen_column) && begin
        _warn_fallback(:q10, GlobalDefault, "no O₂ rate column for Q10 regression")
        return DEFAULT_Q10, FieldProvenance(GlobalDefault)
    end
    metabolic_rates = Float64.(df[!, oxygen_column])
    fitted_q10      = _fit_q10(temperatures, metabolic_rates)
    return fitted_q10, FieldProvenance(Measured, nrow(df), taxon_str)
end
