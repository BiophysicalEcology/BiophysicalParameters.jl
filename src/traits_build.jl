
# ── traits.build integration ───────────────────────────────────────────────────
#
# join_contexts has moved to TraitDataSources.jl and is re-exported from
# BiophysicalParameters for backwards compatibility.
#
# pivot_traits_build_wide stays here: it is a domain-aware BiophysicalParameters
# utility that callers invoke after gettraits() with domain-specific id_cols.

"""
    pivot_traits_build_wide(traits_long, id_cols) → DataFrame

Pivot a traits.build long-format DataFrame to wide format with `trait_name(unit)`
column headers, suitable for use with BiophysicalParameters.jl parameter builders.

Equivalent to the R pipeline:
```R
df %>%
  mutate(trait_name_with_unit = paste(trait_name, "(", unit, ")", sep = "")) %>%
  pivot_wider(names_from = trait_name_with_unit, values_from = value)
```

Value columns are converted from String to Float64 (missing where conversion fails).

# Arguments
- `traits_long`: long-format DataFrame, typically the output of `join_contexts`
- `id_cols`: column names (Symbols or Strings) that uniquely identify each row
  together with the trait name. Typically includes taxon_name, observation_id,
  repeat_measurements_id, and any context columns (Ta, MR_estimate_type).
"""
function pivot_traits_build_wide(traits_long::DataFrame, id_cols)
    df = select(traits_long, vcat(string.(id_cols), ["trait_name", "value", "unit"]))
    trait_name_with_unit = df.trait_name .* "(" .* df.unit .* ")"
    df[!, :trait_name_with_unit] = trait_name_with_unit

    wide = unstack(df, string.(id_cols), :trait_name_with_unit, :value)

    # Convert all string columns to Float64 where all non-missing values parse successfully.
    # This handles both trait value columns and numeric id columns such as Ta and replicates.
    # Non-numeric id columns (taxon_name, MR_estimate_type, etc.) are left as strings.
    for col in names(wide)
        col_values = wide[!, col]
        eltype(col_values) <: Union{Missing, AbstractString} || continue
        parsed = [ismissing(v) ? missing : tryparse(Float64, string(v)) for v in col_values]
        all(x -> ismissing(x) || isa(x, Float64), parsed) || continue
        wide[!, col] = parsed
    end
    return wide
end
