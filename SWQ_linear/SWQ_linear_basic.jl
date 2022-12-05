# using Pkg
# Pkg.add("Gridap")
# Pkg.add("SparseMatricesCSR")
# Pkg.add("SparseArrays")
# Pkg.add("WriteVTK")
# Pkg.add("ProgressBars")
# Pkg.activate(".")

using Gridap
using SparseMatricesCSR
using SparseArrays
using WriteVTK
using ProgressBars

#Solves linear shallow water equations on a 2d plane

function uh(u₀,h₀,X,Y,dΩ)
    a((u,h),(w,ϕ)) = ∫(w⋅u + ϕ*h)dΩ
    b((w,ϕ)) = ∫(w⋅u₀ + ϕ*h₀)dΩ
    solve(AffineFEOperator(a,b,X,Y))
end

function new_vtk_step(Ω,file,_cellfields)
    n = size(_cellfields)[1]
    createvtk(Ω,
              file,
              cellfields=_cellfields,
              nsubcells=n)
end

function perp(u,n)
    n×u
 end
 const ⟂ = perp

function Gridap.get_free_dof_values(functions...)
    map(get_free_dof_values,functions)
  end


function Shallow_water_theta_newton(
        order,degree,h₀,u₀,α,
        linear_solver::Gridap.Algebra.LinearSolver=Gridap.Algebra.BackslashSolver(),
        sparse_matrix_type::Type{<:AbstractSparseMatrix}=SparseMatrixCSC{Float64,Int})
    #Create model
    B = 1
    L = 1
    dx = 1
    dy = 1
    domain = (0,B,0,L)
    
    partition = (100,100)
    dir = "./RESULTS/"
    model = CartesianDiscreteModel(domain,partition)


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
    

    reffe_rt = ReferenceFE(raviart_thomas,Float64,order)
    V = FESpace(model,reffe_rt;conformity=:HDiv,dirichlet_tags = "boundary")
    U = TrialFESpace(V,VectorValue(0.0,0.0))

    reffe_lgn = ReferenceFE(lagrangian,Float64,order)
    Q = FESpace(model,reffe_lgn;conformity=:L2)
    P = TrialFESpace(Q)

    # reffe_lgn = ReferenceFE(lagrangian,Float64,order+1)
    # S = FESpace(model,reffe_lgn;conformity=:H1)
    # R = TrailFESpace(S)

    X = MultiFieldFESpace([V,Q])#∇u, ∇h
    Y = MultiFieldFESpace([U,P])

    E = [0 -1; 1 0]
    #Create initial solutions
    a1(u,v) = ∫(v⋅u)dΩ
    l1(v) = ∫(v⋅u₀)dΩ
    un = solve(AffineFEOperator(a1,l1,U,V))

    a2(u,v) = ∫(v*u)dΩ
    l2(v) = ∫(v*h₀)dΩ
    hn = solve(AffineFEOperator(a2,l2,P,Q))

    b = 2.0
    unv,hnv = get_free_dof_values(un,hn)
    E = [0 -1; 1 0 ]
    uhn = uh(un,hn,X,Y,dΩ)
    un, hn = uhn

    function run_simulation(pvd=nothing)
        N = 100
        dt = 100
        g = 9.81
        for step = ProgressBar(1:N)
            function residual((u,h),(w,ϕ))
                #+ α*dt*f*(w⋅(E×u))
                #+ α*dt*f*(w⋅(E×un))
                ∫((w⋅u) - g*α*dt*DIV(w)*h - w⋅un - α*dt*g*DIV(w)*hn
                + ϕ*h + α*dt*(b)*ϕ*DIV(u) -ϕ*hn + α*dt*(b)*ϕ*DIV(u))dΩ
            end

            assem = SparseMatrixAssembler(sparse_matrix_type,Vector{Float64},X,Y)
            op = FEOperator(residual,X,Y,assem)
            nls = NLSolver(linear_solver)
            solver = FESolver(nls)
            solve!(uhn,solver,op)


            unv .= get_free_dof_values(un)
            hnv .= get_free_dof_values(hn)
            # println("$(step)/$(N) complete")
            #add = output results
            pvd[dt*Float64(step)] = createvtk(Ω,joinpath(dir,"SWQ_linear_$(dt*step)"*".vtu"),cellfields=["u"=>un,"h"=>hn])
        end
    end
    if isdir(dir)
        pvdfile = joinpath(dir,"SWQ_linear"*".pvd")
        paraview_collection(run_simulation,pvdfile)
    else
        mkdir(dir)
        pvdfile = joinpath(dir,"SWQ_linear"*".pvd")
        paraview_collection(run_simulation,pvdfile)
    end
end

function h₀((x,y))
    if (x < 0.25 && y < 0.25)
        h = 0.0
    elseif (x > 0.75 && y > 0.75)
        h = 0.1
    else
        # h = sin(2*π*x)*sin(2*π*y)
        h = 0.5
    end
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

Shallow_water_theta_newton(1,3,h₀,u₀,0.5)