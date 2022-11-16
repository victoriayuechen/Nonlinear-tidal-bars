using Pkg
Pkg.activate(".")

using Revise
using Serialization
using Profile, ProfileView # to profile the code, compute tree rep. 
using Plots 


includet("src/params.jl")                 # Params(...)
includet("src/allocate_arrays.jl")        # allocatArrays(...)
includet("src/grid.jl")                   # Grid(...)
includet("src/initial_bottomheight.jl")   # h₀_coscos(...) etc.

includet("src/write_output.jl")           # make_output_dir(...), output_grid(...)
includet("src/swec.jl")                   # spinup(...), SWEC(...)

"""
run a hydrodynamical simulation
"""
function swec_run()
    #make output folder
    output_path = make_output_dir(test=true)

    #make struct with default parameters + changed ones
    Pars = Params(
        B                = 1e3,
        φ                = 50.0,
        # t_end            = 3 * 44700.0, # integrate for 3 tidal cycles
        output_hydro     = true,
        output_residual  = true,
        output_path      = output_path,
    );

    # allocate arrays
    @unpack nx, ny, U = Pars;
    workingArrays = allocateArrays(nx,ny);  # initial conditions ζ,u,v,c are zeros

    # make grid
    myGrid = Grid(Pars);
    output_grid(Pars,myGrid);
    
    # initial bottom height
    # h₀ = h₀_bump(0.1,200,Pars,myGrid);
    h₀ = h₀_coscos(0.1,1,1,Pars,myGrid);
    # h₀ = h₀_random(0.1,Pars);

    workingArrays.h .= h₀;             # put h₀ in workingArrays
    # serialize(output_path *"h",h₀)     # save h₀ to file
    
    # This is the part where the actual simulation happens 
    spinup(workingArrays,Pars,2;verbose=true) # spin up for 2 tidal cycles (i.e., run swec twice without saving)
    ProfileView.@profview SWEC(workingArrays,Pars);                 # integrate ζ,u,v,c in time.
end
#rm("output/test",recursive=true)

sol = swec_run();
