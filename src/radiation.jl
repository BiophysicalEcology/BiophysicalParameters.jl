"""
    solar_weighted_absorptivity(reflectance, measurement_wavelengths; kwargs...)
    → Float64

Compute solar-weighted absorptivity from a measured spectral reflectance curve.

Integrates:  absorptivity = 1 − ∫ Iλ(λ) ρ(λ) dλ / ∫ Iλ(λ) dλ

Units are kept throughout: wavelengths in nm, irradiance in W m⁻² nm⁻¹.
The absorptivity result is dimensionless.

# Arguments
- `reflectance` — measured reflectance values (0–1), dimensionless
- `measurement_wavelengths` — wavelengths at which reflectance was measured.
  Accepts `Unitful.Length` quantities (any unit) or bare `Real` values (assumed nm).

# Keyword arguments
- `solar_wavelengths` — wavelength grid (default: `SolarRadiation.DEFAULT_WAVELENGTHS`, nm)
- `solar_spectral_irradiance` — spectral irradiance (default: `SolarRadiation.DEFAULT_SOLAR_SPECTRAL_IRRADIANCE`, W m⁻² nm⁻¹)
"""
function solar_weighted_absorptivity(
    reflectance::AbstractVector{<:Real},
    measurement_wavelengths;
    solar_wavelengths         = SolarRadiation.DEFAULT_WAVELENGTHS,
    solar_spectral_irradiance = SolarRadiation.DEFAULT_SOLAR_SPECTRAL_IRRADIANCE,
)
    λ_solar = solar_wavelengths
    λ_meas  = uconvert.(u"nm", _to_nm.(measurement_wavelengths))

    # InterpolationsUnitfulExt handles Unitful length grids directly
    interpolated_reflectance = linear_interpolation(λ_meas, Float64.(reflectance); extrapolation_bc=Flat())
    ρ_on_solar_grid = interpolated_reflectance.(λ_solar)

    # Integrate: Iλ [W/m²/nm] × dλ [nm] → W/m²; ratio is dimensionless
    total_irradiance      = _trapz(solar_spectral_irradiance,                     λ_solar)
    weighted_reflectance  = _trapz(solar_spectral_irradiance .* ρ_on_solar_grid,  λ_solar)

    return 1.0 - ustrip(NoUnits, weighted_reflectance / total_irradiance)
end

_to_nm(x::Unitful.Length) = uconvert(u"nm", x)
_to_nm(x::Real)           = Float64(x) * u"nm"

function _trapz(y::AbstractVector, x::AbstractVector)
    n  = length(x)
    total = zero(first(y) * (x[2] - x[1]))
    for i in 2:n
        total += (y[i] + y[i-1]) * (x[i] - x[i-1]) / 2
    end
    return total
end

# ─── Location-specific solar spectrum ────────────────────────────────────────

"""
    location_solar_spectrum(latitude, elevation; kwargs...)
    → (wavelengths::Vector{Quantity}, spectral_irradiance::Vector{Quantity})

Compute a clear-sky solar spectrum representative of a given location, using
SolarRadiation.jl. Returns the global (direct + diffuse) spectral irradiance at
solar noon on the summer solstice.

The spectrum varies with latitude (solar zenith angle changes the atmospheric path
length, shifting the UV:NIR ratio) and elevation (less atmosphere above means more
UV reaching the surface).

# Arguments
- `latitude` — latitude in degrees (positive north). Accepts `Real` or `Unitful.°`.
- `elevation` — elevation above sea level. Accepts `Real` (assumed m) or `Unitful.Length`.

# Keyword arguments
- `scattered_uv` — use diffuse UV scattering model (default: `false`; `true` is more accurate but slower and the difference is small)
- `solar_model` — `SolarProblem` instance; overrides `scattered_uv` when supplied explicitly
- `albedo` — surface reflectance (default: 0.15, typical soil/vegetation)
- `slope` — terrain slope (default: 0°, flat)
- `aspect` — terrain aspect (default: 0°)

# Notes
The solstice day used is day 355 for the Southern Hemisphere (December) and day 172
for the Northern Hemisphere (June), maximising solar elevation for each.
Atmospheric pressure is derived from elevation using the standard atmosphere formula.
Horizon angles are set to 0° (open flat terrain).
"""
function location_solar_spectrum(
    latitude,
    elevation;
    scattered_uv::Bool         = false,
    solar_model::SolarProblem  = SolarProblem(; scattered_uv),
    albedo::Real               = 0.15,
    slope                      = 0.0u"°",
    aspect                     = 0.0u"°",
)
    lat_deg  = _to_degrees(latitude)
    elev_m   = _to_metres(elevation)

    # Standard atmosphere pressure at elevation: P = P₀·exp(−h/8500 m)
    atmospheric_pressure = 101325.0u"Pa" * exp(-ustrip(u"m", elev_m) / 8500.0)

    # Summer solstice: day 355 (SH, Dec 21) for southern latitudes, day 172 (NH, Jun 21) otherwise
    solstice_day = ustrip(u"°", lat_deg) < 0 ? 355.0 : 172.0

    terrain = SolarTerrain(;
        elevation            = elev_m,
        latitude             = lat_deg,
        longitude            = 0.0u"°",       # noon spectrum is longitude-independent
        slope                = _to_degrees(slope),
        aspect               = _to_degrees(aspect),
        albedo               = Float64(albedo),
        atmospheric_pressure = atmospheric_pressure,
        horizon_angles       = zeros(36) .* u"°",  # flat open terrain
    )

    # Run for a single day at half-hour resolution around midday
    hours  = collect(10.0:0.5:14.0)
    result = solar_radiation(solar_model; solar_terrain=terrain, days=[solstice_day], hours)

    # Pick the timestep with the highest total irradiance (solar noon)
    peak_step = argmax(result.global_horizontal)
    spectrum  = result.global_spectra[peak_step, :]

    return solar_model.wavelengths, spectrum
end

_to_degrees(x::Real)           = Float64(x) * u"°"
_to_degrees(x::Unitful.DimensionlessQuantity) = uconvert(u"°", x)
_to_degrees(x)                 = uconvert(u"°", x)
_to_metres(x::Real)            = Float64(x) * u"m"
_to_metres(x::Unitful.Length)  = uconvert(u"m", x)

# ─── Per-individual absorptivity from raw spectral data ──────────────────────

"""
    per_individual_absorptivity(df_spectral; kwargs...) → DataFrame

Compute solar-weighted absorptivity for each individual in a raw spectral DataFrame.

Handles the case where the same individual was measured at multiple temperatures (e.g.
Smith et al. 2016: 15 °C and 40 °C), returning one row per unique
`(individual_id, region_of_measurement, temperature)` combination.

# Required columns
- `reflectance` — spectral reflectance (0–1)
- one wavelength column: `wave_length_nm`, `wavelength_nm`, or `wavelength`

# Optional columns (used for grouping if present)
- `individual_id` or `Individual_ID`
- `region_of_measurement`
- `temp_C` or `temperature_body`
- `taxon_name` or `Species`

# Returns
`DataFrame` with columns: `taxon_name`, `individual_id`, `region_of_measurement`,
`temperature_celsius`, `absorptivity`, `reflectance_integrated`
"""
function per_individual_absorptivity(
    df_spectral::DataFrame;
    solar_wavelengths         = SolarRadiation.DEFAULT_WAVELENGTHS,
    solar_spectral_irradiance = SolarRadiation.DEFAULT_SOLAR_SPECTRAL_IRRADIANCE,
)
    wavelength_col = _detect_wavelength_column(df_spectral)
    taxon_col      = _detect_column(df_spectral, [:taxon_name, :Species], nothing)
    id_col         = _detect_column(df_spectral, [:individual_id, :Individual_ID], nothing)
    region_col     = _detect_column(df_spectral, [:region_of_measurement], nothing)
    temp_col       = _detect_column(df_spectral, [:temp_C, :temperature_body], nothing)

    group_cols = Symbol[]
    taxon_col  !== nothing && push!(group_cols, taxon_col)
    id_col     !== nothing && push!(group_cols, id_col)
    region_col !== nothing && push!(group_cols, region_col)
    temp_col   !== nothing && push!(group_cols, temp_col)

    rows = NamedTuple[]
    for grp in groupby(df_spectral, group_cols)
        key = first(grp)
        sorted = sort(grp, wavelength_col)

        λ = Float64.(sorted[!, wavelength_col]) .* u"nm"
        ρ = Float64.(sorted.reflectance)
        abs_val = solar_weighted_absorptivity(ρ, λ; solar_wavelengths, solar_spectral_irradiance)

        push!(rows, (
            taxon_name            = taxon_col  !== nothing ? key[taxon_col]  : missing,
            individual_id         = id_col     !== nothing ? key[id_col]     : missing,
            region_of_measurement = region_col !== nothing ? key[region_col] : missing,
            temperature_celsius   = temp_col   !== nothing ? key[temp_col]   : missing,
            absorptivity          = abs_val,
            reflectance_integrated = 1.0 - abs_val,
        ))
    end
    return DataFrame(rows)
end

function _detect_wavelength_column(df)
    for col in [:wave_length_nm, :wavelength_nm, :wavelength]
        hasproperty(df, col) && return col
    end
    error("No wavelength column found. Expected one of: wave_length_nm, wavelength_nm, wavelength")
end

function _detect_column(df, candidates, fallback)
    for col in candidates
        hasproperty(df, col) && return col
    end
    return fallback
end

# ─── Detect which data format is present ─────────────────────────────────────

"""
    radiation_data_format(df) → Symbol

Detect the format of a radiation DataFrame:
- `:spectral` — has wavelength + reflectance columns; requires integration
- `:integrated` — has pre-computed absorptivity column (`abs`, `absortivity`, or `absorptivity`)
- `:integrated_reflectance` — has reflectance only (no wavelength); absorptivity = 1 - reflectance
"""
function radiation_data_format(df::DataFrame)
    has_wavelength = any(hasproperty(df, col) for col in [:wave_length_nm, :wavelength_nm, :wavelength])
    has_abs        = any(hasproperty(df, col) for col in [:abs, :absortivity, :absorptivity])
    has_rho        = hasproperty(df, :reflectance) || hasproperty(df, :rho)

    if has_wavelength
        return :spectral
    elseif has_abs
        return :integrated
    elseif has_rho
        return :integrated_reflectance
    else
        error("Cannot determine radiation data format from columns: $(names(df))")
    end
end

# ─── Default values ───────────────────────────────────────────────────────────

const DEFAULT_ABSORPTIVITY_DORSAL  = 0.85
const DEFAULT_ABSORPTIVITY_VENTRAL = 0.85
const DEFAULT_EMISSIVITY_DORSAL    = 0.95
const DEFAULT_EMISSIVITY_VENTRAL   = 0.95
const DEFAULT_SKY_VIEW_FACTOR      = 0.5
const DEFAULT_GROUND_VIEW_FACTOR   = 0.5
const DEFAULT_VENTRAL_FRACTION     = 0.5

# ─── Builder ─────────────────────────────────────────────────────────────────

"""
    radiation_parameters(taxon, df_radiation, df_geometry; kwargs...)
    → (RadiationParameters, ParameterProvenance)

Build a `RadiationParameters` struct from radiation and geometry DataFrames.

## Absorptivity (estimation trait)

Accepts three radiation data formats, detected automatically:

1. **Raw spectral** (has `wave_length_nm`/`wavelength_nm` column): calls
   `per_individual_absorptivity` to solar-weight each spectral curve, then averages
   across individuals.

2. **Integrated scalar** (has `abs`/`absortivity`/`absorptivity` column, no wavelength):
   reads the pre-computed absorptivity directly and averages across individuals.

3. **Integrated reflectance** (has `reflectance`/`rho` column, no wavelength):
   computes `absorptivity = 1 − mean(reflectance)`.

In all formats the `region_of_measurement` column (values `"dorsal"` / `"ventral"`) is
used to split dorsal and ventral absorptivity.

## Areas and fractions (parameter traits from geometry data)
`total_area`, `silhouette_area`, `conduction_area` (m²) and `ventral_fraction` (dimensionless).

## View factors
Always default — they are habitat geometry parameters, not organism traits.

## Keyword arguments
- `solar_wavelengths`, `solar_spectral_irradiance` — override the default solar spectrum
  (passed through to `solar_weighted_absorptivity`; only used for spectral format data)
"""
function radiation_parameters(
    taxon,
    df_radiation,
    df_geometry;
    solar_wavelengths         = SolarRadiation.DEFAULT_WAVELENGTHS,
    solar_spectral_irradiance = SolarRadiation.DEFAULT_SOLAR_SPECTRAL_IRRADIANCE,
)
    provenance = Dict{Symbol,FieldProvenance}()

    # ── Absorptivities ───────────────────────────────────────────────────────
    absorptivity_dorsal, absorptivity_ventral = _absorptivities_from_data(
        taxon, df_radiation, solar_wavelengths, solar_spectral_irradiance, provenance,
    )

    # ── Emissivities ─────────────────────────────────────────────────────────
    emissivity_dorsal  = _scalar_from_trait_data(
        df_radiation, "emissivity_dorsal", DEFAULT_EMISSIVITY_DORSAL,
        :body_emissivity_dorsal, provenance,
    )
    emissivity_ventral = _scalar_from_trait_data(
        df_radiation, "emissivity_ventral", DEFAULT_EMISSIVITY_VENTRAL,
        :body_emissivity_ventral, provenance,
    )

    # ── Surface areas (m²) ───────────────────────────────────────────────────
    total_area      = _area_from_geometry_data(df_geometry, "total_area",      0.0u"m^2", :total_area,      provenance)
    silhouette_area = _area_from_geometry_data(df_geometry, "silhouette_area", 0.0u"m^2", :silhouette_area, provenance)
    conduction_area = _area_from_geometry_data(df_geometry, "conduction_area", 0.0u"m^2", :conduction_area, provenance)

    # ── Ventral fraction ─────────────────────────────────────────────────────
    ventral_fraction = _scalar_from_trait_data(
        df_geometry, "ventral_fraction", DEFAULT_VENTRAL_FRACTION, :ventral_fraction, provenance,
    )

    # ── View factors — habitat geometry, always defaulted ───────────────────
    for field in (:sky_view_factor, :ground_view_factor, :vegetation_view_factor, :bush_view_factor)
        provenance[field] = FieldProvenance(GlobalDefault)
        _warn_fallback(field, GlobalDefault, "view factors reflect habitat geometry, not organism traits")
    end

    params = HeatExchange.RadiationParameters(;
        body_absorptivity_dorsal  = absorptivity_dorsal,
        body_absorptivity_ventral = absorptivity_ventral,
        body_emissivity_dorsal    = emissivity_dorsal,
        body_emissivity_ventral   = emissivity_ventral,
        total_area                = total_area,
        silhouette_area           = silhouette_area,
        conduction_area           = conduction_area,
        sky_view_factor           = DEFAULT_SKY_VIEW_FACTOR,
        ground_view_factor        = DEFAULT_GROUND_VIEW_FACTOR,
        vegetation_view_factor    = 0.0,
        bush_view_factor          = 0.0,
        ventral_fraction          = ventral_fraction,
    )

    return params, NamedTuple(provenance)
end

# ─── Internal helpers ─────────────────────────────────────────────────────────

function _absorptivities_from_data(taxon, df_radiation, solar_wavelengths, solar_spectral_irradiance, provenance)
    if isnothing(df_radiation) || isempty(df_radiation)
        for field in (:body_absorptivity_dorsal, :body_absorptivity_ventral)
            _warn_fallback(field, GlobalDefault, "no radiation data supplied")
            provenance[field] = FieldProvenance(GlobalDefault)
        end
        return DEFAULT_ABSORPTIVITY_DORSAL, DEFAULT_ABSORPTIVITY_VENTRAL
    end

    fmt = radiation_data_format(df_radiation)
    taxon_str = string(taxon)

    # Filter to this taxon if a taxon column is present
    taxon_col = _detect_column(df_radiation, [:taxon_name, :Species], nothing)
    df = if taxon_col !== nothing
        filter(r -> r[taxon_col] == taxon_str, df_radiation)
    else
        df_radiation
    end

    if isempty(df)
        for field in (:body_absorptivity_dorsal, :body_absorptivity_ventral)
            _warn_fallback(field, TaxonDefault, "no radiation data found for $taxon_str; using defaults")
            provenance[field] = FieldProvenance(TaxonDefault, 0, taxon_str)
        end
        return DEFAULT_ABSORPTIVITY_DORSAL, DEFAULT_ABSORPTIVITY_VENTRAL
    end

    if fmt == :spectral
        # Compute per-individual absorptivities, then average by region
        per_individual = per_individual_absorptivity(df; solar_wavelengths, solar_spectral_irradiance)
        abs_dorsal  = _mean_absorptivity_by_region(per_individual, "dorsal",  :body_absorptivity_dorsal,  taxon_str, provenance)
        abs_ventral = _mean_absorptivity_by_region(per_individual, "ventral", :body_absorptivity_ventral, taxon_str, provenance)

    elseif fmt == :integrated
        abs_col = _detect_column(df, [:abs, :absortivity, :absorptivity], nothing)
        abs_dorsal  = _mean_column_by_region(df, abs_col, "dorsal",  :body_absorptivity_dorsal,  taxon_str, provenance)
        abs_ventral = _mean_column_by_region(df, abs_col, "ventral", :body_absorptivity_ventral, taxon_str, provenance)

    else  # :integrated_reflectance
        rho_col = _detect_column(df, [:reflectance, :rho], nothing)
        abs_dorsal  = 1.0 - _mean_column_by_region(df, rho_col, "dorsal",  :body_absorptivity_dorsal,  taxon_str, provenance)
        abs_ventral = 1.0 - _mean_column_by_region(df, rho_col, "ventral", :body_absorptivity_ventral, taxon_str, provenance)
    end

    return abs_dorsal, abs_ventral
end

function _mean_absorptivity_by_region(per_individual_df, region, field, taxon_str, provenance)
    rows = filter(r -> !ismissing(r.region_of_measurement) && r.region_of_measurement == region, per_individual_df)
    if isempty(rows)
        default = region == "dorsal" ? DEFAULT_ABSORPTIVITY_DORSAL : DEFAULT_ABSORPTIVITY_VENTRAL
        _warn_fallback(field, TaxonDefault, "no $region spectral data for $taxon_str; using default $default")
        provenance[field] = FieldProvenance(TaxonDefault, 0, taxon_str)
        return default
    end
    provenance[field] = FieldProvenance(Measured, nrow(rows), taxon_str)
    return mean(rows.absorptivity)
end

function _mean_column_by_region(df, value_col, region, field, taxon_str, provenance)
    region_col = _detect_column(df, [:region_of_measurement], nothing)
    rows = if region_col !== nothing
        filter(r -> !ismissing(r[region_col]) && r[region_col] == region, df)
    else
        df
    end
    if isempty(rows)
        default = region == "dorsal" ? DEFAULT_ABSORPTIVITY_DORSAL : DEFAULT_ABSORPTIVITY_VENTRAL
        _warn_fallback(field, TaxonDefault, "no $region data for $taxon_str; using default $default")
        provenance[field] = FieldProvenance(TaxonDefault, 0, taxon_str)
        return default
    end
    provenance[field] = FieldProvenance(Measured, nrow(rows), taxon_str)
    return mean(Float64.(rows[!, value_col]))
end

function _scalar_from_trait_data(df, trait_name, default_value, field, provenance)
    if isnothing(df) || isempty(df)
        provenance[field] = FieldProvenance(GlobalDefault)
        return default_value
    end
    hasproperty(df, :trait_name) || begin
        provenance[field] = FieldProvenance(GlobalDefault)
        return default_value
    end
    rows = filter(r -> r.trait_name == trait_name, df)
    if isempty(rows)
        _warn_fallback(field, GlobalDefault, "no $trait_name row found; using default $default_value")
        provenance[field] = FieldProvenance(GlobalDefault)
        return default_value
    end
    provenance[field] = FieldProvenance(Measured, nrow(rows), "")
    return mean(rows.value)
end

function _area_from_geometry_data(df_geometry, trait_name, default_value, field, provenance)
    if isnothing(df_geometry) || isempty(df_geometry)
        _warn_fallback(field, GlobalDefault, "no geometry data supplied")
        provenance[field] = FieldProvenance(GlobalDefault)
        return default_value
    end
    hasproperty(df_geometry, :trait_name) || begin
        _warn_fallback(field, GlobalDefault, "geometry DataFrame has no trait_name column")
        provenance[field] = FieldProvenance(GlobalDefault)
        return default_value
    end
    rows = filter(r -> r.trait_name == trait_name, df_geometry)
    if isempty(rows)
        _warn_fallback(field, GlobalDefault, "no $trait_name row found; using default $default_value")
        provenance[field] = FieldProvenance(GlobalDefault)
        return default_value
    end
    provenance[field] = FieldProvenance(Measured, nrow(rows), "")
    return mean(rows.value) * u"m^2"
end
