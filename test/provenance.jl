using BiophysicalParameters
using Test

@testset "completeness_score" begin
    prov_all_measured = (
        a = FieldProvenance(Measured, 10, "Pogona vitticeps"),
        b = FieldProvenance(Measured, 5, "Pogona vitticeps"),
    )
    @test completeness_score(prov_all_measured) == 1.0

    prov_none_measured = (
        a = FieldProvenance(GlobalDefault),
        b = FieldProvenance(TaxonDefault, 0, "Squamata"),
    )
    @test completeness_score(prov_none_measured) == 0.0

    prov_mixed = (
        a = FieldProvenance(Measured, 10, "Pogona vitticeps"),
        b = FieldProvenance(GlobalDefault),
        c = FieldProvenance(GlobalDefault),
        d = FieldProvenance(GlobalDefault),
    )
    @test completeness_score(prov_mixed) == 0.25
end

@testset "data_gaps prints without error" begin
    prov = (
        body_absorptivity_dorsal = FieldProvenance(Measured, 12, "Pogona vitticeps"),
        sky_view_factor          = FieldProvenance(GlobalDefault),
    )
    @test_nowarn data_gaps(prov; name="test")
end
