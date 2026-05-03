using BiophysicalParameters
using SolarRadiation
using Test

@testset "solar_weighted_absorptivity" begin
    # Flat reflectance of 1.0 → absorptivity = 0.0
    λ = [300.0, 700.0, 1500.0, 4000.0]
    @test solar_weighted_absorptivity(ones(4), λ) ≈ 0.0 atol=1e-10

    # Flat reflectance of 0.0 → absorptivity = 1.0
    @test solar_weighted_absorptivity(zeros(4), λ) ≈ 1.0 atol=1e-10

    # Flat reflectance of 0.3 → absorptivity ≈ 0.7
    @test solar_weighted_absorptivity(fill(0.3, 4), λ) ≈ 0.7 atol=1e-10
end

@testset "radiation_parameters — no data falls back to defaults" begin
    params, prov = @test_logs (:warn,) (:warn,) (:warn,) (:warn,) (:warn,) (:warn,) (:warn,) (:warn,) match_mode=:any radiation_parameters(
        nothing, nothing, nothing
    )
    @test completeness_score(prov) == 0.0
    @test params isa BiophysicalParameters.HeatExchange.RadiationParameters
end

@testset "radiation_parameters — measured reflectance" begin
    using DataFrames
    # Minimal df_radiation with flat 30% reflectance on dorsal and ventral
    wavelengths = [400.0, 700.0, 1000.0, 2000.0]
    rows = vcat(
        DataFrame(taxon_name="Pogona vitticeps", trait_name="reflectance_dorsal",
                  wavelength_nm=wavelengths, value=fill(0.3, 4)),
        DataFrame(taxon_name="Pogona vitticeps", trait_name="reflectance_ventral",
                  wavelength_nm=wavelengths, value=fill(0.5, 4)),
        DataFrame(taxon_name="Pogona vitticeps", trait_name="emissivity_dorsal",
                  wavelength_nm=missing, value=[0.95]),
        DataFrame(taxon_name="Pogona vitticeps", trait_name="emissivity_ventral",
                  wavelength_nm=missing, value=[0.95]),
    )

    params, prov = radiation_parameters("Pogona vitticeps", rows, nothing)

    @test prov.body_absorptivity_dorsal.source == Measured
    @test prov.body_absorptivity_ventral.source == Measured
    @test params.body_absorptivity_dorsal ≈ 0.7 atol=0.01
    @test params.body_absorptivity_ventral ≈ 0.5 atol=0.01
end
