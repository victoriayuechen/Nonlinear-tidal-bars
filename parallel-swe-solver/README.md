# Parallel SWE Solver 

This implements our solvers in a parallel manner using [GridapDistributed](https://github.com/gridap/GridapDistributed.jl) and [GridapPETSc](https://github.com/carlodev/GridapPETSc.jl) (by carlodev). This requires a working MPI library and [MPI.jl](https://github.com/JuliaParallel/MPI.jl). 


The solvers correspond to the sequential ones in the `swe-solver/` by name, only with `_parallel` appended at the end.  