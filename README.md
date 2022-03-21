# ElectronicStructurePySCF

[![Build Status](https://github.com/jlapeyre/ElectronicStructurePySCF.jl/workflows/CI/badge.svg)](https://github.com/jlapeyre/ElectronicStructurePySCF.jl/actions)
[![Coverage](https://codecov.io/gh/jlapeyre/ElectronicStructurePySCF.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jlapeyre/ElectronicStructurePySCF.jl)


This package includes types and methods to integrate [`pyscf`](https://github.com/pyscf/pyscf) with
[ElectronicStructure.jl](https://github.ibm.com/John-Lapeyre/ElectronicStructure.jl).


You can try this, to help work around some (but by no means all) problems with loading PyCall:
```julia
ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python")
```
Note this must be done before loading any Julia package that depends on PyCall, including `ElectronicStructurePySCF.jl`.
