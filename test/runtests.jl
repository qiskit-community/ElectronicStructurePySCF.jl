using ElectronicStructurePySCF
using ElectronicStructure: Atom, Geometry, MolecularData, MolecularSpec, to_pyscf
using Test

using Artifacts
using LazyArtifacts
import JLD2

function get_molecular_data()
    rootpath = artifact"molecular_data"
    filename = joinpath(rootpath, "molecular_data_examples.jld2")
    data = JLD2.load(filename)
    return data
end

@testset "ElectronicStructurePySCF.jl" begin
    atom = Atom(:Li, (0.0, 0.0, 1.4))
    @test to_pyscf(atom) == "Li 0.0 0.0 1.4"
    geom = Geometry(Atom(:H, (0., 0., 0.)), Atom(:H, (0., 0., 0.7414)))
    @test to_pyscf(geom) == "H 0.0 0.0 0.0;H 0.0 0.0 0.7414"
end

@testset "compare data file" begin
    data_top = get_molecular_data()
    h2 = data_top["data"]["H2_sto-3g"]
    geom = Geometry(Atom(:H, (0., 0., 0.)), Atom(:H, (0., 0., 0.7414)))
    ms = MolecularSpec(geometry=geom)
    pymd = to_pyscf(ms)
    @test pymd.atom == "H 0.0 0.0 0.0;H 0.0 0.0 0.7414"
    @test pymd.atom == to_pyscf(geom)
    @test pymd.charge == 0
    @test pymd.multiplicity == 1
end
