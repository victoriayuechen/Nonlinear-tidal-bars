using Test, Revise

includet("params.jl")
includet("initialize.jl")
includet("utils.jl")
includet("write_output.jl")
includet("swec.jl")
includet("morphodynamics.jl")

function periodic_bc_test(array)
    boel=true
    for i in 1:ny
        boel = boel & (array[1,:] == array[end-1,:]) & (array[end,:] == array[2,:])
    end
    return boel
end

@testset "periodic boundary conditions" begin
    Pars = Params()
    @unpack nx, ny = Pars;   
    A = allocateArrays(nx,ny);  # initial conditions ζ,u,v,c, h are zeros
    myGrid = Grid(Pars);
    # A.h .= h₀_coscos(0.1,1,1,Pars,myGrid);
    A.h .= h₀_random(0.1,Pars);
    @unpack Δt_num, t_span, save_moments_t = Pars;
    p_swec           = (A,Pars,false);
    savedvalues_swec = SavedValues(Float64, Array{eltype(A.h),2});
    cb_swec          = SavingCallback(save_div_q⃗, saveat=save_moments_t, savedvalues_swec);
    prob_swec        = ODEProblem(rhs_swec!,A.ic,t_span,p_swec,save_everystep=true);
    integrator_swec  = init(prob_swec,RK4(),dt=Δt_num,callback=cb_swec);
    solve!(integrator_swec); 

    @test periodic_bc_test(A.h)
    @test periodic_bc_test(A.continuityᵤ)
    @test periodic_bc_test(A.continuityᵥ)
    @test periodic_bc_test(A.frictionᵤ)
    @test periodic_bc_test(A.frictionᵥ)
    @test periodic_bc_test(A.pressgradᵤ)
    @test periodic_bc_test(A.pressgradᵥ)
    @test periodic_bc_test(A.advectionᵤ)
    @test periodic_bc_test(A.advectionᵥ)
    @test periodic_bc_test(A.coriolisᵤ)
    @test periodic_bc_test(A.coriolisᵥ)
    @test periodic_bc_test(A.D)
    @test periodic_bc_test(A.Dᵤ)
    @test periodic_bc_test(A.Dᵥ)
    @test periodic_bc_test(A.Dᵤu)
    @test periodic_bc_test(A.Dᵥv)
    @test periodic_bc_test(A.uᵥ)
    @test periodic_bc_test(A.vᵤ)
    @test periodic_bc_test(A.Ûᵤ)
    @test periodic_bc_test(A.Ûᵥ)
    @test periodic_bc_test(A.Û)
    @test periodic_bc_test(A.u_eff²)
    @test periodic_bc_test(A.u_effᵤ²)
    @test periodic_bc_test(A.u_effᵥ²)
    @test periodic_bc_test(A.Hu)
    @test periodic_bc_test(A.Huᵤ)
    @test periodic_bc_test(A.Huᵥ)
    @test periodic_bc_test(A.dudx)
    @test periodic_bc_test(A.dudy)
    @test periodic_bc_test(A.dvdx)
    @test periodic_bc_test(A.dvdy)
    @test periodic_bc_test(A.dcdx)
    @test periodic_bc_test(A.dcdy)
    @test periodic_bc_test(A.dζdx)
    @test periodic_bc_test(A.dζdy)
    @test periodic_bc_test(A.dhdx)
    @test periodic_bc_test(A.dhdy)
    @test periodic_bc_test(A.cᵤ)
    @test periodic_bc_test(A.cᵥ)
    @test periodic_bc_test(A.c_bottomᵤ)
    @test periodic_bc_test(A.c_bottomᵥ)
    @test periodic_bc_test(A.c_bottom)
    @test periodic_bc_test(A.c_topᵤ)
    @test periodic_bc_test(A.c_topᵥ)
    @test periodic_bc_test(A.erosion)
    @test periodic_bc_test(A.deposition)
    @test periodic_bc_test(A.qs_advectionᵤ)
    @test periodic_bc_test(A.qs_advectionᵥ)
    @test periodic_bc_test(A.qs_diffᵤ)
    @test periodic_bc_test(A.qs_diffᵥ)
    @test periodic_bc_test(A.qsᵤ)
    @test periodic_bc_test(A.qsᵥ)
    @test periodic_bc_test(A.dqsᵤdx)
    @test periodic_bc_test(A.dqsᵥdy)
    @test periodic_bc_test(A.div_q⃗s)
    @test periodic_bc_test(A.qb_advectionᵤ)
    @test periodic_bc_test(A.qb_advectionᵥ)
    @test periodic_bc_test(A.qb_diffᵤ)
    @test periodic_bc_test(A.qb_diffᵥ)
    @test periodic_bc_test(A.qbᵤ)
    @test periodic_bc_test(A.qbᵥ)
    @test periodic_bc_test(A.dqbᵤdx)
    @test periodic_bc_test(A.dqbᵥdy)
    @test periodic_bc_test(A.div_q⃗b)
    @test periodic_bc_test(A.avg_div_q⃗)
end

@testset "zero v at walls" begin
    Pars = Params()
    @unpack nx, ny = Pars;   
    A = allocateArrays(nx,ny);  # initial conditions ζ,u,v,c, h are zeros
    myGrid = Grid(Pars);
    A.h .= h₀_coscos(0.1,1,1,Pars,myGrid);
    # A.h .= h₀_random(0.1,Pars,myGrid);
    @unpack Δt_num, t_span, save_moments_t = Pars;
    p_swec           = (A,Pars,false);
    savedvalues_swec = SavedValues(Float64, Array{eltype(A.h),2});
    cb_swec          = SavingCallback(save_div_q⃗, saveat=save_moments_t, savedvalues_swec);
    prob_swec        = ODEProblem(rhs_swec!,A.ic,t_span,p_swec,save_everystep=true);
    integrator_swec  = init(prob_swec,RK4(),dt=Δt_num,callback=cb_swec);
    solve!(integrator_swec); 
    vB = @view integrator_swec.sol[:,end,3,:]
    v0 = @view integrator_swec.sol[:,1,3,:]
    @test vB == zeros(size(vB))
    @test v0 == zeros(size(v0))
end

@testset "lateral sediment transport is zero at walls" begin
    Pars = Params(output_residual=true)
    @unpack nx, ny = Pars;   
    A = allocateArrays(nx,ny);  # initial conditions ζ,u,v,c, h are zeros
    myGrid = Grid(Pars);
    A.h .= h₀_coscos(0.1,1,1,Pars,myGrid);
    # A.h .= h₀_random(0.1,Pars,myGrid);
    @unpack Δt_num, t_span, save_moments_t = Pars;
    p_swec           = (A,Pars,false);
    savedvalues_swec = SavedValues(Float64, Tuple{Array{Float64,2},Array{Float64,2}});
    cb_swec          = SavingCallback(save_func_swec, saveat=save_moments_t, savedvalues_swec);
    prob_swec        = ODEProblem(rhs_swec!,A.ic,t_span,p_swec,save_everystep=true);
    integrator_swec  = init(prob_swec,RK4(),dt=Δt_num,callback=cb_swec);
    solve!(integrator_swec); 
    
   
    for tn in 1:length(savedvalues_swec.saveval)
        qvB = @view savedvalues_swec.saveval[tn][2][:,end]
        qv0 = @view savedvalues_swec.saveval[tn][2][:,1]
        @test qvB == zeros(size(qvB))
        @test qv0 == zeros(size(qv0))
    end

end
