
# ── traits.build integration ───────────────────────────────────────────────────
#
# Functions for working with traits.build database objects loaded from .rds files.
# The traits.build format stores context variables (treatment temperature, metabolic
# rate estimate type, etc.) in a separate `contexts` table. They must be joined onto
# the `traits` table before the data can be used by parameter builders.
#
# This is the Julia equivalent of the R function:
#   join_context_properties(include_description = FALSE, format = "many_columns")

"""
    join_contexts(traits, contexts) → DataFrame

Join context variables from a traits.build `contexts` table onto the long-format
`traits` table, adding one column per context property (e.g. `Ta`, `MR_estimate_type`).

This is the Julia equivalent of `traits.build::join_context_properties()` used in R
analysis scripts. Context variables (ambient temperature, metabolic rate estimate type,
etc.) are stored separately in `contexts` and must be joined before the data can be
filtered or pivoted.

# Arguments
- `traits`: the `traits` slot from a traits.build .rds file (loaded via RData.jl)
- `contexts`: the `contexts` slot from the same file

# Returns
A copy of `traits` with additional columns for each unique `context_property` in the
contexts table. Values are strings; convert to numeric types after pivoting to wide
format as needed.
"""
function join_contexts(traits::DataFrame, contexts::DataFrame)
    isempty(contexts) && return copy(traits)

    result = copy(traits)
    for property_name in unique(contexts.context_property)
        property_rows = filter(r -> r.context_property == property_name, contexts)
        # All rows for a single context property share the same link_id column name
        link_column = Symbol(first(property_rows.link_id))

        # Context IDs are scoped per dataset — must key lookup by (dataset_id, link_vals)
        # to avoid collisions when multiple datasets share the same numeric context IDs.
        lookup = Dict{Tuple{String,String}, String}(
            (r.dataset_id, r.link_vals) => r.value
            for r in eachrow(property_rows)
        )

        result[!, Symbol(property_name)] = [
            (ismissing(row[link_column]) || !haskey(lookup, (row.dataset_id, row[link_column]))) ?
                missing : lookup[(row.dataset_id, row[link_column])]
            for row in eachrow(result)
        ]
    end
    return result
end

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
