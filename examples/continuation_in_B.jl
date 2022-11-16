
#----------------------------------------------------------------------------------------------------
# This script is a messy not optimized and will not work locally in its current state (I only ran it on the server).
# At some point I had the problem that B was not correctly updated. The first increase was ok but the second was way to big. 
# I didnt investigate this further.
# So take this script maybe only as inspiration/a proof of concept.
# In addition, for a pseudo-arc length continuation in cd, one would not use the for loop below, but something like:
# optcont = ContinuationPar(some options)
# br, _ = continuation(f, jac, out, par_af, (@lens _.cd), optcont, filename = output_path * "cd")
# See docu of BifurcationKit.jl for more information.     
#----------------------------------------------------------------------------------------------------

# activate virtual environment on the server
using Pkg
Pkg.activate("path to model folder")
Pkg.instantiate()

using Revise, Dates, Serialization
using Interpolations
using ForwardDiff:  jacobian
using LinearAlgebra, BifurcationKit, Setfield
const BK = BifurcationKit

includet("params.jl")                 # Params(...)
includet("allocate_arrays.jl")        # allocatArrays(...)
includet("grid.jl")                   # Grid(...)
includet("initial_bottomheight.jl")   # h₀_coscos(...) etc.

includet("morphodynamics.jl")
includet("write_output.jl")
includet("read_output.jl");

"""
Make a bottom pattern h, NΔy wider, with N ∈ ℤ, by streching and interpolating.
returns the new h and the new y grid
see docu of Interpolations package for details.
"""
function add_NΔy(h,x,y,N)
    Δy = y[2]-y[1]
    Ny = length(y)
    B_old = y[end]-y[1] + Δy
    B_new = B_old + N*Δy
    B_new_cgrid = B_new - Δy
    y_shifted_scaled = LinRange(0,B_new_cgrid,Ny)
    x_range = x[1]:x[2]-x[1]:x[end]
    itp = interpolate(h, BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp,x_range,y_shifted_scaled) 
    y_new_shifted = LinRange(0,B_new_cgrid,Ny+N) 
    h_new = [sitp(x,y) for x in x, y in y_new_shifted]
    y_new = y_new_shifted .+ Δy/2 # physical domain
    return h_new, y_new
end

""" previous_out is folder name and N is number of Δy to add (positive or negative)"""
function find_equilibrium(start_run,N)
    
    Δxₑ, Δyₑ, = 100.0, 50.0  # ₑ stands for external (used in the other functions as external/global-ish param)
    nxₑ = 103
    Lₑ  = 10e3
    

    #read initial bottom pattern
    path_output_folder_start_run = "path to model folder/output/"* start_run * "/arrays/"
    fid = h5open(path_output_folder_start_run * "continuation_B.h5","cw")
    fid_B = fid["700.0"]
    out = Array(fid_B["out_1"])
    close(fid)

    output_path = make_output_dir(
        out_path="path to model folder/output/",
        scr_path="path to model folder/",
        description_long = "
        Continuation start from  B = 700, tol = 1e-5
        lati=0, uc=0, uw=0, 
        ∇ζ_term = off, c_top_bottom_dep_on_D  = false,
        rest default.",
        logbook     = true,
        logbookFile = "path to model folder/logbook.csv",
        description_logbook = "Continuation B = 700, tol = 1e-5, N = 2"
        )
    
    # for loop to widen/narrowing the model domain 50 times.
    for _ in 1:50

        ny_ic = Int(size(out,1)/(nxₑ-2) + 1) # with 'extra' gridpoint
        B_ic  = Δyₑ * (ny_ic - 1)
        equi = reshape(out,(nxₑ-2,ny_ic-1))
        x_grid = collect(range(0,     stop = Lₑ,         length = nxₑ-2))   # physcial domain
        y_ic   = collect(range(Δyₑ/2, stop = B_ic-Δyₑ/2, length = ny_ic-1)) # physcial domain
        
        if !iszero(N)
            h₀, y_grid = add_NΔy(equi,x_grid,y_ic,N) #scale and interpolate to wider/narrower C grid
            h₀ = vec(h₀)
            
            ny_new = (ny_ic) + N #including the 'extra' gridpoint
            @assert(length(y_grid) == ny_new - 1, "ny_new wrong value")
            B_new  = Δyₑ * (ny_new - 1)
        else
            h₀ = out
            y_grid = y_ic
            ny_new = ny_ic
            B_new = B_ic
        end
        
        @show B_new
          
        """ calculate rhs of dhdτ = f(h) """
        function f(h,p)
            # convert vector h into Matrix h̃ (with ghost points)
            h̃ = zeros(eltype(h),nxₑ,ny_new) # Matrix h̃ (with ghost points)
            h̃[2:nxₑ-1,1:ny_new-1] .= reshape(h,(nxₑ-2,ny_new-1))
            periodic_bc!(h̃)

            Pars = Params(
                B                      = B_new,
                L                      = Lₑ,
                φ                      = 0.0,
                u_critical             = 0.0,     # avoids heaviside function
                u_orbital              = 0.0,    
                c_diffusion_∇ζ_term_on = false,
                c_top_bottom_dep_on_D  = false,   # reduces number of exp() calculations
                C∞_mode                = true,    # needed when u_critical≠0 or linear_friction = false
            );
            
            # allocate arrays
            @unpack nx, ny = Pars;  
            @assert((nx == nxₑ) & (ny == ny_new))
            A = allocateArrays(eltype(h̃),nx,ny);  # initial conditions ζ,u,v,c, h are zeros
            
            # make grid
            myGrid = Grid(Pars);
            
            #initialize h and calculate gradient
            A.h .= h̃;
            ∇h!(A,Pars)

            @unpack Δt_num, t_span, save_moments_t = Pars
            p_swec           = (A,Pars,true)
            savedvalues_swec = SavedValues(Float64, Array{eltype(h),2})
            cb_swec          = SavingCallback(save_div_q⃗, saveat=save_moments_t, savedvalues_swec)
            prob_swec        = ODEProblem(rhs_swec!,A.ic,t_span,p_swec,save_everystep=false)
            integrator_swec  = init(prob_swec,RK4(),adaptive=false,dt=Δt_num,callback=cb_swec)

            #integrate swec for 2 tidalcycles, first spinup, than real one
            solve!(integrator_swec)
            #update ic with last values of ζ, u, v, c
            A.ic .= integrator_swec.sol[:,:,:,end]
            # integrate again
            reinit!(integrator_swec, A.ic)
            solve!(integrator_swec)
            
            #calculate average divergence of total sediment transport
            mean_sediment_transport!(A,savedvalues_swec.saveval,Pars)
            
            #calculate rhs
            @unpack p̂,ε = Pars
            dhdτ = vec( -1.0/(1.0-p̂) .* A.avg_div_q⃗[2:nx-1,1:ny-1] ./ ε ) # make vector of dhdτ on on physical domain
            
            #replace equation of last gridpoint with integral condition (conservation of sediment)
            dhdτ[end] = zero(eltype(dhdτ))
            for i in eachindex(h)
                dhdτ[end] +=  h[i]  
            end
            
            #force first gridpoint to be zero to fix phase of pattern.
            dhdτ[1] = h[1]

            return dhdτ
        end

        """ calculate jacobian """
        @inline function jac(h,p)
            f̃= z -> f(z,p)
            jacobian(f̃,h)
        end

        """ callback to save output Newton """
        function cb(x,f,jac,res,it,itl,optN; kwargs...)
            name = string(it)
            namejac = string(it-1)

            # Note that in Newton.jl (BifurcationKit.jl) the callback is called multiple times in 1 newton iteration step to determine the flag.
            # Hence, the necessity to check if datasets exists (and a simple check for it = 0 is not enough).
            h5open(output_path * "continuation_B.h5","cw") do fid
                if !exists(fid, string(B_new))
                    fid_B = g_create(fid,string(B_new))
                else
                    fid_B = fid[string(B_new)]
                end
                if !exists(fid_B,"x")
                    fid_B["x"]  = x_grid
                end
                if !exists(fid_B,"y")
                    fid_B["y"]  = collect(y_grid)
                end
                if !exists(fid_B,"out_$name")
                    fid_B["out_$name"]  = x
                end
                if !exists(fid_B,"res_$name")
                    fid_B["res_$name"]  = res
                end
                if !exists(fid_B,"jac_$namejac") && it > 0
                    fid_B["jac_$namejac"]  = jac
                end
            end
            return true
        end

        #Newton iterations
        par_af = (p1= 0.0, p2 = 0.0) # not used
        optnewton = NewtonPar(tol = 1e-5, verbose = false, maxIter= 10)
        out, hist,flag, nr_its = BK.newton(f, jac, h₀, par_af , optnewton, normN = x -> norm(x, Inf64),callback=cb)
        
        #calculate jacobian of last iterate
        #This chould (should?) be changed to a jacobian of f(h) without the dirichlet and integral condition (and mass_matrix below removed)
        jac_end = jac(out,par_af)
        
        #calculate eigenvalues and eigenvectors of last jacobian
        N,M = size(jac_end)
        @assert N == M
        mass_matrix = Matrix(Diagonal([zeros(1); ones(N-2) ; zeros(1)]))
        F = eigen(Array(jac_end),mass_matrix)
        p = sortperm(F.values, by = x-> real(x), rev = true)
        eigval = F.values[p]
        eigvec = F.vectors[:,p]

        # save jacobian and its spectral decomposition
        h5open(output_path * "continuation_B.h5","cw") do fid
            fid_B = fid[string(B_new)]
            fid_B["jac_end"]    = jac_end
            fid_B["eigval_end"] = eigval
            fid_B["eigvec_end"] = eigvec
        end
    end

end

find_equilibrium("2020_10_01__09",2);
