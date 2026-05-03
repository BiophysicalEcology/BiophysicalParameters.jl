"""
    DataSource

Records how a parameter value was obtained, from most to least preferred.
"""
@enum DataSource begin
    Measured         # from the DataFrame supplied by the user
    Allometric       # PowerLaw prediction at lowest available taxon rank
    PhylogeneticMean # mean of nearest relatives present in the database
    TaxonDefault     # hardcoded literature value for the class/order/family
    GlobalDefault    # Param() default embedded in HeatExchange.jl struct definition
end

"""
    FieldProvenance

Records the provenance of a single parameter field.

# Fields
- `source` — how the value was obtained
- `n_obs` — number of individual observations behind the value (0 for defaults)
- `taxon_used` — taxon whose data was actually used (may differ from the query taxon when
  falling back to a relative)
"""
struct FieldProvenance
    source::DataSource
    n_obs::Int
    taxon_used::String
end

FieldProvenance(source::DataSource) = FieldProvenance(source, 0, "")

"""
    ParameterProvenance

A `NamedTuple` mapping parameter struct field names to their `FieldProvenance`.
Returned alongside each parameter struct by every builder function.
"""
const ParameterProvenance = NamedTuple

"""
    completeness_score(prov) → Float64

Fraction of fields in `prov` whose source is `Measured` (0.0–1.0).
"""
function completeness_score(prov::NamedTuple)
    fields = values(prov)
    isempty(fields) && return 1.0
    sum(f.source == Measured for f in fields) / length(fields)
end

"""
    data_gaps(prov; name="")

Print a structured table of parameter sources. Fields not yet `Measured` are
highlighted to guide future data collection.
"""
function data_gaps(prov::NamedTuple; name::String="")
    header = isempty(name) ? "Parameter provenance" : "Parameter provenance — $name"
    println(header)
    println("─" ^ length(header))
    for (k, v) in pairs(prov)
        marker = v.source == Measured ? "✓" : "!"
        taxon_note = isempty(v.taxon_used) ? "" : " ($(v.taxon_used))"
        nobs_note = v.n_obs > 0 ? " [n=$(v.n_obs)]" : ""
        println("  $marker  $(rpad(string(k), 30))  $(v.source)$(taxon_note)$(nobs_note)")
    end
    score = completeness_score(prov)
    println("─" ^ length(header))
    println("  Completeness: $(round(100score; digits=1))%  ($(sum(v.source == Measured for v in values(prov)))/$(length(prov)) fields measured)")
end

# Internal helper used by all builders to emit a warning when falling back
function _warn_fallback(field::Symbol, source::DataSource, reason::String)
    @warn "Using $source for $field: $reason"
end
