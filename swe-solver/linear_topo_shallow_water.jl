using Pkg
Pkg.activate(".")

using Gridap
using SparseMatricesCSR
using SparseArrays
using WriteVTK
using LinearAlgebra
#Solves linear shallow water equations on a 2d plane

function uh(u₀,h₀,F₀,X,Y,dΩ)
    a((u,h,u2),(w,ϕ,w2)) = ∫(w⋅u + ϕ*h + w2⋅u2)dΩ
    b((w,ϕ,w2)) = ∫(w⋅u₀ + ϕ*h₀ + w2⋅F₀)dΩ
    solve(AffineFEOperator(a,b,X,Y))
end


function compute_mass_flux!(F,dΩ,V,RTMMchol,u)
    b(v) = ∫(v⋅u)dΩ
    Gridap.FESpaces.assemble_vector!(b, get_free_dof_values(F), V)
    ldiv!(RTMMchol,get_free_dof_values(F))
end


function new_vtk_step(Ω,file,_cellfields)
    n = size(_cellfields)[1]
    createvtk(Ω,
              file,
              cellfields=_cellfields,
              nsubcells=n)
end

function setup_and_factorize_mass_matrices(dΩ,V,Q,U,P)
    amm(a,b) = ∫(a⋅b)dΩ

    RTMM = assemble_matrix(amm,U,V)

    RTMMchol = lu(RTMM)

    RTMM,RTMMchol
end

clone_fe_function(space,f)=FEFunction(space,copy(get_free_dof_values(f)))

function perp(u,n)
    n×u
 end
 const ⟂ = perp

function Gridap.get_free_dof_values(functions...)
    map(get_free_dof_values,functions)
  end


function Shallow_water_theta_newton(
        order,degree,h₀,u₀,θ,topography,
        linear_solver::Gridap.Algebra.LinearSolver=Gridap.Algebra.BackslashSolver(),
        sparse_matrix_type::Type{<:AbstractSparseMatrix}=SparseMatrixCSC{Float64,Int})
    #Create model
    B = 1.0
    L = 1.0
    nlrtol = 1.0e-02
    dx = 1
    dy = 1
    domain = (0,B,0,L)
    
    partition = (100,100)
    dir = "swe-solver/output_linear_swe"
    model = CartesianDiscreteModel(domain,partition;isperiodic=(true,true))


    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,2,5])
    add_tag_from_tags!(labels,"left",[7])
    add_tag_from_tags!(labels,"right",[8])
    add_tag_from_tags!(labels,"top",[3,4,6])
    add_tag_from_tags!(labels,"inside",[9])
    nuemanntaggs = ["left","right"]

    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    Γ = BoundaryTriangulation(model,tags=nuemanntaggs)
    dΓ = Measure(Γ,degree)
    dω = Measure(Ω,degree,ReferenceDomain())


    

    reffe_rt = ReferenceFE(raviart_thomas,Float64,order)
    V = FESpace(model,reffe_rt;conformity=:HDiv,dirichlet_tags = "boundary")
    U = TrialFESpace(V)

    reffe_lgn = ReferenceFE(lagrangian,Float64,order)
    Q = FESpace(model,reffe_lgn;conformity=:L2)
    P = TrialFESpace(Q)
    
    RTMM,RTMMchol = setup_and_factorize_mass_matrices(dΩ,V,Q,U,P)

    # reffe_lgn = ReferenceFE(lagrangian,Float64,order+1)
    # S = FESpace(model,reffe_lgn;conformity=:H1)
    # R = TrailFESpace(S)

    X = MultiFieldFESpace([V,Q,V])#∇u, ∇h
    Y = MultiFieldFESpace([U,P,U])
    #Create initial solutions
    a1(u,v) = ∫(v⋅u)dΩ
    l1(v) = ∫(v⋅u₀)dΩ
    un = solve(AffineFEOperator(a1,l1,U,V))

    a2(u,v) = ∫(v*u)dΩ
    l2(v) = ∫(v*h₀)dΩ
    hn = solve(AffineFEOperator(a2,l2,P,Q))

    b = interpolate_everywhere(topography,P)
    unv,hnv = get_free_dof_values(un,hn)
    F₀ = clone_fe_function(V,un)
    compute_mass_flux!(F₀,dΩ,V,RTMMchol,un*hn)
    uhn = uh(un,hn,F₀,X,Y,dΩ)
    du, du, F = uhn


    E = [0 -1; 1 0 ]

    
    function run_simulation(pvd=nothing)
        N = 100
        dt = 100
        g = 9.81
        for step=1:N
            function residual((du,dh,F),(w,ϕ,w2))
                one_θ = (1.0-θ)
                u_star = un + one_θ*du
                h_star = hn + one_θ*dh
                h_η = hn + b + one_θ*dh

                ∫((1.0/dt)*w⋅(du) + (∇⋅(w))*(g*h_η)
                    + (1.0/dt)*ϕ*(dh))dΩ + ∫(ϕ*(DIV(F)))dω +
                ∫(w2⋅(F - u_star*h_star))dΩ
            end

            # function residual((u,h,F),(w,ϕ,w2))

            #     #+ α*dt*f*(w⋅(E×u))
            #     #+ α*dt*f*(w⋅(E×un))
            #     ∫((w⋅u) - g*α*dt*DIV(w)*h - w⋅un - α*dt*g*DIV(w)*hn
            #     + ϕ*h + α*dt*(b)*ϕ*DIV(u) -ϕ*hn + α*dt*(b)*ϕ*DIV(u) + w2⋅(F-h*u))dΩ
            # end
            
            dY = get_fe_basis(Y)
            residualdudhF = residual(uhn,dY)
            r = assemble_vector(residualdudhF,Y)
            assem = SparseMatrixAssembler(sparse_matrix_type,Vector{Float64},X,Y)
            op = FEOperator(residual,X,Y,assem)
            nls = NLSolver(linear_solver;
                            show_trace = true,
                            method=:newton,
                            ftol = nlrtol*norm(r,Inf),
                            xtol = 1.0e-02)
            solver = FESolver(nls)
            solve!(uhn,solver,op)


            unv .= unv + get_free_dof_values(du)
            hnv .= unv + get_free_dof_values(dh)
            println("$(step)/$(N) complete")
            #add = output results
            pvd[dt*Float64(step)] = createvtk(Ω,joinpath(dir,"nswe_cells_linear_shallow_water_$(dt*step)"*".vtu"),cellfields=["u"=>un,"h"=>hn,"F" =>F])
        end
    end
    if isdir(dir)
        pvdfile = joinpath(dir,"nswe_cells_linear_shallow_water"*".pvd")
        paraview_collection(run_simulation,pvdfile)
    else
        mkdir(dir)
        pvdfile = joinpath(dir,"nswe_cells_linear_shallow_water"*".pvd")
        paraview_collection(run_simulation,pvdfile)
    end
end

function h₀((x,y))
    h = sin(2*π*x)*sin(2*π*y)
    h
end

function topography((x,y))
    p = 0.5*sin(π*x)*sin(π*y)
    p
end

function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

Shallow_water_theta_newton(1,3,h₀,u₀,0.9,topography)