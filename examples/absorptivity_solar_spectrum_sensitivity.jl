# Absorptivity sensitivity to solar spectrum: elevation and latitude
#
# Demonstrates how the solar-weighted absorptivity of Pogona vitticeps
# varies when the solar spectrum is computed for specific locations rather
# than using the default (standard atmosphere) spectrum.
#
# The mechanism: at high elevation and low latitude, the shorter atmospheric
# path length shifts the spectrum toward shorter wavelengths (more UV/blue
# relative to NIR). Since lizard skin absorbs strongly in the UV, absorptivity
# increases with elevation and proximity to the equator.

using BiophysicalParameters
using CSV
using DataFrames
using Statistics
using Unitful

const RADIATION_DB = joinpath(
    homedir(), "Dropbox", "Current Research Projects",
    "trait_database", "heat_budget_databases", "radiationDB",
)

# ── 1. Build mean spectral reflectance curve for Pogona vitticeps ─────────────
spectral_data = CSV.read(
    joinpath(RADIATION_DB, "data", "Smith_etal_2016", "data.csv"),
    DataFrame,
)

# Average reflectance across all individuals and temperatures, separately for
# dorsal and ventral, to get a representative population-mean spectral curve.
function mean_spectral_reflectance(df, region)
    subset = filter(r -> r.region_of_measurement == region, df)
    combine(groupby(subset, :wave_length_nm), :reflectance => mean => :reflectance)
end

dorsal_spectrum  = mean_spectral_reflectance(spectral_data, "dorsal")
ventral_spectrum = mean_spectral_reflectance(spectral_data, "ventral")

println("Mean spectral curves built from $(length(unique(spectral_data.Individual_ID))) individuals")
println("Wavelength range: $(minimum(dorsal_spectrum.wave_length_nm))–$(maximum(dorsal_spectrum.wave_length_nm)) nm")
println()

# ── 2. Absorptivity with default spectrum (standard atmosphere) ───────────────
abs_dorsal_default  = solar_weighted_absorptivity(dorsal_spectrum.reflectance,  dorsal_spectrum.wave_length_nm)
abs_ventral_default = solar_weighted_absorptivity(ventral_spectrum.reflectance, ventral_spectrum.wave_length_nm)
println("── Default spectrum (standard atmosphere) ───────────────────────────")
println("  Dorsal absorptivity:  $(round(abs_dorsal_default;  digits=4))")
println("  Ventral absorptivity: $(round(abs_ventral_default; digits=4))")
println()

# ── 3. Sensitivity to elevation ───────────────────────────────────────────────
# Walpeup collection latitude, varying elevation 0–4000 m.
collection_latitude = -35.14   # degrees

elevations = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000] .* u"m"

println("── Sensitivity to elevation (latitude = $(collection_latitude)°) ────────────")
println("  elevation_m  abs_dorsal  abs_ventral  Δ_dorsal_vs_sea_level")
for elev in elevations
    lam, Ilam = location_solar_spectrum(collection_latitude, ustrip(u"m", elev), scattered_uv=false)
    abs_d = solar_weighted_absorptivity(dorsal_spectrum.reflectance,  dorsal_spectrum.wave_length_nm;
                                         solar_wavelengths=lam, solar_spectral_irradiance=Ilam)
    abs_v = solar_weighted_absorptivity(ventral_spectrum.reflectance, ventral_spectrum.wave_length_nm;
                                         solar_wavelengths=lam, solar_spectral_irradiance=Ilam)
    delta = abs_d - abs_dorsal_default
    println("  $(lpad(elev, 11))  $(round(abs_d; digits=4))      $(round(abs_v; digits=4))      $(round(delta; digits=4))")
end
println()

# ── 4. Sensitivity to latitude ────────────────────────────────────────────────
# Sea level, varying latitude from equator to high southern latitudes.
latitudes_deg = [-70, -65, -60, -55, -50, -45, -40, -35, -25, -15, -5, 5, 15, 25, 35, 40, 45, 50, 55, 60, 65, 70]

println("── Sensitivity to latitude (elevation = 0 m) ────────────────────────")
println("  latitude_deg  abs_dorsal  abs_ventral  Δ_dorsal_vs_default")
for lat in latitudes_deg
    lam, Ilam = location_solar_spectrum(lat, 0.0)
    abs_d = solar_weighted_absorptivity(dorsal_spectrum.reflectance,  dorsal_spectrum.wave_length_nm;
                                         solar_wavelengths=lam, solar_spectral_irradiance=Ilam)
    abs_v = solar_weighted_absorptivity(ventral_spectrum.reflectance, ventral_spectrum.wave_length_nm;
                                         solar_wavelengths=lam, solar_spectral_irradiance=Ilam)
    delta = abs_d - abs_dorsal_default
    println("  $(lpad(lat, 12))  $(round(abs_d; digits=4))      $(round(abs_v; digits=4))      $(round(delta; digits=4))")
end
println()

# ── 5. Joint sensitivity: elevation × latitude ────────────────────────────────
elevations_grid = [0, 1000, 2000, 3000, 4000] .* u"m"
latitudes_grid  = [-35, -20, -5, 10, 25]

println("── Joint sensitivity: dorsal absorptivity (elevation × latitude) ────")
header = rpad("elev\\lat", 10) * join(lpad.(latitudes_grid, 8), "")
println(header)
for elev in elevations_grid
    row = rpad(string(ustrip(u"m", elev)) * "m", 10)
    for lat in latitudes_grid
        lam, Ilam = location_solar_spectrum(lat, elev)
        abs_d = solar_weighted_absorptivity(dorsal_spectrum.reflectance, dorsal_spectrum.wave_length_nm;
                                             solar_wavelengths=lam, solar_spectral_irradiance=Ilam)
        row *= lpad(round(abs_d; digits=3), 8)
    end
    println(row)
end
println()
println("Default (standard atmosphere): $(round(abs_dorsal_default; digits=3))")
