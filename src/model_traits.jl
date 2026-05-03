"""
    resolve_insulation_type(; df_insulation=nothing) → AbstractInsulation

Infer the appropriate insulation type from whether insulation data is available.
Returns `Naked()` when `df_insulation` is `nothing` or empty, otherwise `FibrousLayer()`.

Taxonomic model traits (thermal strategy, activity period) cannot be inferred from
data and must be supplied explicitly to `build_organism`.
"""
function resolve_insulation_type(; df_insulation=nothing)
    if isnothing(df_insulation) || isempty(df_insulation)
        return HeatExchange.Naked()
    else
        return HeatExchange.FibrousLayer()
    end
end
