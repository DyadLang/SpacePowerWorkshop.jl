using SpacePowerWorkshop
using DyadInterface
using Test


@testset "Battery energy never goes below 0.5" begin
    res = SpacePowerWorkshop.SolarPanelOrbitalTransient()

    stored_energy = res.sol(res.time_rel; idxs = res.spec.model.converter.stored_energy)

    @test all(>(0.5), stored_energy)
    # if that's false, try `findall(<(0.5), stored_energy[:])` to see where it's failing
    @test all(<(5), stored_energy)
end
