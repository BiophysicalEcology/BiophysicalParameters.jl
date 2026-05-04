
# ── ParameterSummary ──────────────────────────────────────────────────────────

"""
    ParameterSummary

Summary statistics and a fitted parametric distribution for one parameter of one taxon.
The distribution can be reconstructed as a `Distributions.jl` object via `distribution(s)`.

Distribution families and their `distribution_params` JSON keys:

| `distribution_family` | JSON keys |
|---|---|
| `"Beta"` | `alpha`, `beta` |
| `"LogNormal"` | `mu`, `sigma` |
| `"Normal"` | `mu`, `sigma` |
| `"PointMass"` | `value` (single default with no observed uncertainty) |
"""
struct ParameterSummary
    mean                ::Float64
    sd                  ::Union{Float64, Missing}
    n                   ::Int
    ci_lower_95         ::Union{Float64, Missing}
    ci_upper_95         ::Union{Float64, Missing}
    distribution_family ::String
    distribution_params ::String               # JSON
    provenance          ::String
    source_reference    ::String
end

function ParameterSummary(default_value::Float64;
                           provenance::String  = "GlobalDefault",
                           source::String      = "HeatExchange.jl default")
    ParameterSummary(
        default_value, missing, 0, missing, missing,
        "PointMass", JSON3.write(Dict("value" => default_value)),
        provenance, source,
    )
end

function ParameterSummary(values::AbstractVector{<:Real};
                           bounds::Tuple{Float64,Float64} = (-Inf, Inf),
                           provenance::String             = "Measured",
                           source::String                 = "")
    vals = Float64.(filter(x -> !ismissing(x) && isfinite(x), values))
    isempty(vals) && return ParameterSummary(NaN; provenance="NoData", source)
    n = length(vals)
    m = mean(vals)
    s = n >= 2 ? std(vals)  : missing
    ci_lo, ci_hi = if n >= 3
        se = std(vals) / sqrt(n)
        t  = quantile(TDist(n - 1), 0.975)
        (m - t * se, m + t * se)
    else
        (missing, missing)
    end
    family, params = _fit_distribution(vals, bounds)
    ParameterSummary(m, s, n, ci_lo, ci_hi, family, JSON3.write(params), provenance, source)
end

"""
    distribution(s::ParameterSummary) → Distribution

Reconstruct a `Distributions.jl` object from the stored distribution family and params.
"""
function distribution(s::ParameterSummary)
    p = JSON3.read(s.distribution_params, Dict{String, Float64})
    if s.distribution_family == "Beta"
        return Beta(p["alpha"], p["beta"])
    elseif s.distribution_family == "LogNormal"
        return LogNormal(p["mu"], p["sigma"])
    elseif s.distribution_family == "Normal"
        return Normal(p["mu"], p["sigma"])
    elseif s.distribution_family == "PointMass"
        return Dirac(p["value"])
    else
        error("Unknown distribution family: $(s.distribution_family)")
    end
end

# ── Distribution fitting ──────────────────────────────────────────────────────

function _fit_distribution(values::AbstractVector{Float64},
                            bounds::Tuple{Float64,Float64})
    n = length(values)
    n < 2 && return ("PointMass", Dict("value" => first(values)))

    lo, hi = bounds
    if lo == 0.0 && hi == 1.0
        clamped = clamp.(values, 1e-6, 1.0 - 1e-6)
        d = Distributions.fit_mle(Beta, clamped)
        return ("Beta", Dict("alpha" => d.α, "beta" => d.β))
    elseif lo == 0.0 && isinf(hi)
        pos = max.(values, 1e-10)
        d = Distributions.fit_mle(LogNormal, pos)
        return ("LogNormal", Dict("mu" => d.μ, "sigma" => d.σ))
    else
        d = Distributions.fit_mle(Normal, values)
        return ("Normal", Dict("mu" => d.μ, "sigma" => d.σ))
    end
end

# ── ParameterRecord ───────────────────────────────────────────────────────────

"""
    ParameterRecord

All available parameter summaries for one taxon, keyed by HeatExchange.jl parameter
struct name (as Symbol) and field name (as Symbol).

Retrieve a summary: `record[:RadiationParameters][:body_absorptivity_dorsal]`
Get a distribution: `distribution(record[:RadiationParameters][:body_absorptivity_dorsal])`
"""
struct ParameterRecord
    taxon_name    ::String
    taxon_class   ::String
    parameters    ::Dict{Symbol, Dict{Symbol, ParameterSummary}}
    date_compiled ::Date
end

Base.getindex(r::ParameterRecord, struct_name::Symbol) =
    get(r.parameters, struct_name, Dict{Symbol, ParameterSummary}())

function distribution(r::ParameterRecord, struct_name::Symbol, field::Symbol)
    d = get(r.parameters, struct_name, nothing)
    d === nothing && error("No parameters for struct $struct_name in record for $(r.taxon_name)")
    s = get(d, field, nothing)
    s === nothing && error("No parameter $field in $struct_name for $(r.taxon_name)")
    distribution(s)
end

# ── Arrow table helpers ───────────────────────────────────────────────────────

function _summary_to_row(taxon::String, taxon_class::String,
                          param_struct::String, field::String, unit::String,
                          s::ParameterSummary)
    (
        taxon_name          = taxon,
        taxon_class         = taxon_class,
        parameter_struct    = param_struct,
        parameter_field     = field,
        unit                = unit,
        mean                = s.mean,
        sd                  = s.sd,
        n                   = Int32(s.n),
        ci_lower_95         = s.ci_lower_95,
        ci_upper_95         = s.ci_upper_95,
        distribution_family = s.distribution_family,
        distribution_params = s.distribution_params,
        provenance          = s.provenance,
        source_reference    = s.source_reference,
        date_compiled       = string(today()),
    )
end

# ── Path resolution ───────────────────────────────────────────────────────────

function _resolve_db_path()
    path = get(ENV, "BIOPHYSICAL_PARAMS_PATH", nothing)
    path !== nothing && return path
    error(
        "Biophysical parameters database path not set.\n" *
        "Set ENV[\"BIOPHYSICAL_PARAMS_PATH\"] to a directory where the compiled " *
        "parameters database can be stored and read.\n" *
        "Example: ENV[\"BIOPHYSICAL_PARAMS_PATH\"] = \"/home/user/biophys_params\""
    )
end

_params_arrow_path(db_path::String) = joinpath(db_path, "parameters.arrow")

# ── Per-domain compile helpers ────────────────────────────────────────────────

# Radiation: per-individual absorptivities from spectral reflectance.
function _compile_radiation_params(taxon::String, taxon_class::String)
    rows = NamedTuple[]
    try
        long = gettraits(HeatBudgetDB{RadiationDomain}(); taxon)
        isempty(long) && return rows

        id_cols = [:taxon_name, :observation_id, :repeat_measurements_id,
                   :entity_type, :region_of_measurement]
        id_cols_present = filter(c -> hasproperty(long, c), id_cols)
        wide = pivot_traits_build_wide(long, id_cols_present)

        renames = Pair{String,String}[]
        for (old, new) in [
            ("wave_length(nm)" => "wave_length_nm"),
            ("temperature_body(Cel)" => "temp_C"),
            ("reflectance({dimensionless})" => "reflectance"),
            ("observation_id" => "Individual_ID"),
        ]
            old in names(wide) && push!(renames, old => new)
        end
        isempty(renames) || rename!(wide, renames...)

        hasproperty(wide, :reflectance) || return rows

        indiv = per_individual_absorptivity(wide)

        for (region, field) in [("dorsal",  "body_absorptivity_dorsal"),
                                 ("ventral", "body_absorptivity_ventral")]
            sub = filter(r -> r.region_of_measurement == region, indiv)
            isempty(sub) && continue
            vals = collect(skipmissing(sub.absorptivity))
            isempty(vals) && continue
            s = ParameterSummary(vals; bounds=(0.0, 1.0), source="radiationDB")
            push!(rows, _summary_to_row(taxon, taxon_class, "RadiationParameters", field, "{dimensionless}", s))
        end
    catch e
        @debug "No radiation absorptivity data for $taxon" exception=e
    end
    return rows
end

# Respiration: per-individual respiratory quotient and metabolic rate.
function _compile_respiration_params(taxon::String, taxon_class::String)
    rows = NamedTuple[]
    try
        long = gettraits(HeatBudgetDB{RespirationDomain}(); taxon)
        isempty(long) && return rows

        id_cols = _respiration_id_cols(long)
        wide = pivot_traits_build_wide(long, id_cols)

        fmt = respiration_data_format(wide)

        # ── Respiratory quotient ─────────────────────────────────────────────
        if fmt == :open_circuit
            o2_col  = _oxygen_rate_column(wide)
            co2_col = _carbon_dioxide_rate_column(wide)
            if o2_col !== nothing && co2_col !== nothing
                rq_vals = filter(isfinite, Float64.(wide[!, co2_col]) ./ Float64.(wide[!, o2_col]))
                if !isempty(rq_vals)
                    s = ParameterSummary(rq_vals; bounds=(0.4, 1.3), source="respirationDB")
                    push!(rows, _summary_to_row(taxon, taxon_class, "RespirationParameters",
                                                "respiratory_quotient", "{dimensionless}", s))
                end
            end
        end

        # ── Metabolic heat flow per row → LogNormal summary ──────────────────
        heat_vals = _per_row_metabolic_heat_watts(wide)
        if !isempty(heat_vals)
            s = ParameterSummary(heat_vals; bounds=(0.0, Inf), source="respirationDB")
            push!(rows, _summary_to_row(taxon, taxon_class, "MetabolismParameters",
                                        "metabolic_heat_flow", "W", s))
        end

        # ── Q10 — only if multiple temperatures ──────────────────────────────
        q10_val = _estimate_q10_from_wide(wide)
        if q10_val !== nothing
            s = ParameterSummary(q10_val; provenance="GlobalDefault", source="respirationDB fit")
            push!(rows, _summary_to_row(taxon, taxon_class, "MetabolismParameters",
                                        "q10", "{dimensionless}", s))
        end

    catch e
        @debug "No respiration data for $taxon" exception=e
    end
    return rows
end

function _respiration_id_cols(long::DataFrame)
    base = [:taxon_name, :observation_id, :entity_type]
    extras = [:repeat_measurements_id, :Ta, :MR_estimate_type, :replicates]
    vcat(base, filter(c -> hasproperty(long, c), extras))
end

# Convert each row of a wide respirometry DataFrame to metabolic heat in watts.
function _per_row_metabolic_heat_watts(wide::DataFrame)
    o2_col   = _metabolic_rate_o2_column(wide)
    o2_col === nothing && (o2_col = _oxygen_rate_column(wide))
    o2_col === nothing && return Float64[]

    mass_col = _mass_column(wide)
    time_unit = _infer_oxygen_time_unit(o2_col)

    vals = Float64[]
    for row in eachrow(wide)
        v = row[o2_col]
        (ismissing(v) || !isfinite(v)) && continue
        o2_ml = Float64(v)
        if _is_mass_specific_o2_column(o2_col)
            mass_col === nothing && continue
            mass = row[mass_col]
            (ismissing(mass) || !isfinite(mass)) && continue
            o2_ml *= Float64(mass)
        end
        # Convert to mL/s then to W (20.1 J/mL O₂)
        o2_ml_s = time_unit == u"mL/hr" ? o2_ml / 3600.0 : o2_ml / 60.0
        push!(vals, o2_ml_s * 20.1)
    end
    return vals
end

function _estimate_q10_from_wide(wide::DataFrame)
    temp_col = _body_temperature_column(wide)
    temp_col === nothing && (temp_col = _ambient_temperature_column(wide))
    temp_col === nothing && return nothing

    o2_col = _metabolic_rate_o2_column(wide)
    o2_col === nothing && (o2_col = _oxygen_rate_column(wide))
    o2_col === nothing && return nothing

    temps = Float64[]
    rates = Float64[]
    for row in eachrow(wide)
        t = row[temp_col]; r = row[o2_col]
        (ismissing(t) || ismissing(r) || !isfinite(t) || !isfinite(r) || r <= 0) && continue
        push!(temps, Float64(t))
        push!(rates, Float64(r))
    end

    length(unique(temps)) < 3 && return nothing
    return _fit_q10(temps, rates)
end

# ── compile_parameters_database ──────────────────────────────────────────────

"""
    compile_parameters_database(taxa; output_path, append) → DataFrame

Build a compiled parameters database for each taxon in `taxa` by loading raw
traits.build data from `ENV["HEATBUDGETDB_PATH"]` and running per-domain
parameter estimators.

Each row of the returned (and written) DataFrame represents one parameter field
for one taxon, with columns for mean, SD, n, 95% CI, parametric distribution
family, and distribution parameters stored as JSON.

Writes `parameters.arrow` to `output_path` (default: `ENV["BIOPHYSICAL_PARAMS_PATH"]`).

# Arguments
- `taxa` — vector of taxon name strings (must match `taxon_name` in traits.build databases)
- `output_path` — directory to write `parameters.arrow`; defaults to `ENV["BIOPHYSICAL_PARAMS_PATH"]`
- `append` — if `true`, merge with any existing `parameters.arrow`, replacing rows for the
  specified taxa; if `false` (default), overwrites
- `taxon_classes` — optional `Dict{String,String}` mapping taxon name → class label
  (e.g. `"Squamate"`, `"EutherianMammal"`). Used for the `taxon_class` column.
"""
function compile_parameters_database(
    taxa::AbstractVector{<:AbstractString};
    output_path::String               = _resolve_db_path(),
    append::Bool                      = false,
    taxon_classes::Dict{String,String} = Dict{String,String}(),
)
    all_rows = NamedTuple[]

    for taxon in taxa
        taxon_class = get(taxon_classes, string(taxon), "Unknown")
        t = string(taxon)
        append!(all_rows, _compile_radiation_params(t, taxon_class))
        append!(all_rows, _compile_respiration_params(t, taxon_class))
    end

    if isempty(all_rows)
        @warn "No parameter data found for any of the specified taxa."
        return DataFrame()
    end

    result = DataFrame(all_rows)

    arrow_path = _params_arrow_path(output_path)
    if append && isfile(arrow_path)
        existing = DataFrame(Arrow.Table(arrow_path))
        taxa_set = Set(string.(taxa))
        existing = filter(r -> !(r.taxon_name in taxa_set), existing)
        result   = vcat(existing, result)
    end

    mkpath(output_path)
    Arrow.write(arrow_path, result)
    return result
end

# ── Query: parameters_table ───────────────────────────────────────────────────

"""
    parameters_table(taxon; db_path) → DataFrame

Return the raw compiled-parameters rows for `taxon` from the Arrow database.
Each row is one parameter field with mean, SD, n, CI, and distribution columns.

Set `ENV["BIOPHYSICAL_PARAMS_PATH"]` or pass `db_path` explicitly.
"""
function parameters_table(taxon::AbstractString;
                           db_path::String = _resolve_db_path())
    arrow_path = _params_arrow_path(db_path)
    isfile(arrow_path) || error(
        "Parameters database not found at $arrow_path.\n" *
        "Run compile_parameters_database([\"$taxon\"]) to build it first."
    )
    df = DataFrame(Arrow.Table(arrow_path))
    filter(r -> r.taxon_name == string(taxon), df)
end

# ── Query: load_parameters ────────────────────────────────────────────────────

"""
    load_parameters(taxon; db_path) → ParameterRecord

Load all compiled parameter summaries for `taxon` from the Arrow database,
returning a `ParameterRecord` struct.  Each parameter summary includes mean, SD, n,
95% CI, and a parametric distribution that can be retrieved via `distribution(record, ...)`.

# Example
```julia
record = load_parameters("Pogona vitticeps")
record[:RadiationParameters][:body_absorptivity_dorsal].mean
distribution(record, :RadiationParameters, :body_absorptivity_dorsal)  # → Beta(…)
```
"""
function load_parameters(taxon::AbstractString;
                          db_path::String = _resolve_db_path())
    df = parameters_table(taxon; db_path)
    isempty(df) && error("No parameters found for taxon \"$taxon\" in database at $db_path")

    taxon_class = first(df.taxon_class)
    params = Dict{Symbol, Dict{Symbol, ParameterSummary}}()

    for row in eachrow(df)
        struct_key = Symbol(row.parameter_struct)
        field_key  = Symbol(row.parameter_field)
        get!(params, struct_key, Dict{Symbol, ParameterSummary}())
        params[struct_key][field_key] = ParameterSummary(
            row.mean,
            row.sd,
            Int(row.n),
            row.ci_lower_95,
            row.ci_upper_95,
            row.distribution_family,
            row.distribution_params,
            row.provenance,
            row.source_reference,
        )
    end

    ParameterRecord(string(taxon), taxon_class, params, today())
end
