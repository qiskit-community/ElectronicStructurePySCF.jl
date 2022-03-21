module ElectronicStructurePySCF

using ElectronicStructure: Atom, Geometry, MolecularData, MolecularSpec
import ElectronicStructure

import PyCall
using PyCall: PyObject

export PyMol, PySCF
export hartree_fock!

## In order to make this a global const it is defined outside the
## try/catch block. See the `PyCall` documentation for this idiom.
const _PYSCF = PyCall.PyNULL()

const electronic_structure = PyCall.PyNULL()

function __init__()
    try
        copy!(_PYSCF, PyCall.pyimport("pyscf"))

        pypath = dirname(@__FILE__)
        if !(pypath in PyCall.pyimport("sys")."path")
            println("Adding ", pypath)
            pushfirst!(PyCall.PyVector(PyCall.pyimport("sys")."path"), pypath)
        end
        copy!(electronic_structure, PyCall.pyimport("electronic_structure"))
    catch
        error("Unable to load python module pyscf.")
    end
end


function ElectronicStructure.MolecularSpec(mol_spec::PyObject)
    class = PyCall.pytypeof(mol_spec).__name__
    if class != "MolecularSpec"
            throw(ErrorException("`MolecularSpec` expecting python class `MolecularSpec`., got ", class, "."))
    end
    # Following is how PyCall returns a Python Geometry object
    geom = Geometry([Atom(Symbol(atom[1]), atom[2]) for atom in mol_spec.geometry[1]]...)
    return MolecularSpec(geom, mol_spec.multiplicity, mol_spec.charge, mol_spec.basis)
end


function ElectronicStructure.MolecularData(mol_data::PyObject)
    class = PyCall.pytypeof(mol_data).__name__
    if class != "MolecularData"
        throw(ErrorException("`MolecularData` expecting python class `MolecularData`., got " * class * "."))
    end
    jspec = MolecularSpec(mol_data.spec)
    return MolecularData(jspec, mol_data.nuclear_repulsion,
                         mol_data.one_body_integrals, mol_data.two_body_integrals)
end


include("wrap_all.jl")


end # module ElectronicStructurePySCF
