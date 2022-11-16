using Revise

includet("src/params.jl")                 # Params(...)
includet("src/allocate_arrays.jl")        # allocatArrays(...)
includet("src/grid.jl")                   # Grid(...)
includet("src/initial_bottomheight.jl")   # h₀_coscos(...) etc.

includet("src/write_output.jl")           # make_output_dir(...), output_grid(...)
includet("src/morphodynamics.jl")         # bed_evolution(...) spinup(...) (morphodynamics includes swec)

function morpho()
    #make output folder
    output_path = make_output_dir(test=true)

    #make struct with default parameters
    Pars = Params(
        B                  = 1e3,
        L                  = 10e3,
        φ                  = 50.0,
        u_orbital          = 0.0,
        u_critical         = 0.0,
        concentration_on   = false,
        t_end_in_years     = 20.0,
        Δτ_out             = 10*6.048,
        # t_end_in_years     = 2000.0,
        output_path = output_path
    );

    @unpack nr_Δτ_steps, τ_end, ε, sec_in_year = Pars
    last_timestep = τ_end/(ε * sec_in_year)

    # allocate arrays
    @unpack nx, ny = Pars;
    workingArrays = allocateArrays(nx,ny);  # initial conditions ζ,u,v,c, h are zeros

    # make grid
    myGrid = Grid(Pars);
    output_grid(Pars,myGrid);

    # make initial bottom pattern
    h₀ = h₀_random(0.1,Pars);
    # h₀ = h₀_coscos(0.1,1,2,Pars,myGrid);
    workingArrays.h .= h₀;

    sol_spin_up = spinup(workingArrays,Pars,1);
    
    bed_evolution(workingArrays,Pars);
end
sol = morpho();
