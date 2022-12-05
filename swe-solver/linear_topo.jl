using Pkg
Pkg.activate(".")

using Gridap
using SparseMatricesCSR
using SparseArrays
using WriteVTK

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
        order,degree,h₀,u₀,α,topography,
        linear_solver::Gridap.Algebra.LinearSolver=Gridap.Algebra.BackslashSolver(),
        sparse_matrix_type::Type{<:AbstractSparseMatrix}=SparseMatrixCSC{Float64,Int})
    #Create model
    B = 1
    L = 1
    dx = 
    dy = 1
    domain = (0,B,0,L)
    
    partition = (100,100)
    dir = "swe-solver/output_linear_swe"
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
    U = TransientTrialFESpace(V)

    reffe_lgn = ReferenceFE(lagrangian,Float64,order)
    Q = FESpace(model,reffe_lgn;conformity=:L2)
    P = TransientTrialFESpace(Q)

    # reffe_lgn = ReferenceFE(lagrangian,Float64,order+1)
    # S = FESpace(model,reffe_lgn;conformity=:H1)
    # R = TrailFESpace(S)

    Y = MultiFieldFESpace([V,Q])#∇u, ∇h
    X = TransientMultiFieldFESpace([U,P])

    E = [0 -1; 1 0]
    #Create initial solutions
    a1(u,v) = ∫(v⋅u)dΩ
    l1(v) = ∫(v⋅u₀)dΩ
    un = solve(AffineFEOperator(a1,l1,U,V))

    a2(u,v) = ∫(v*u)dΩ
    l2(v) = ∫(v*h₀)dΩ
    hn = solve(AffineFEOperator(a2,l2,P,Q))

    b = interpolate_everywhere(topography,P(0.0))
    unv,hnv = get_free_dof_values(un,hn)
    E = [0 -1; 1 0 ]
    uhn = uh(un,hn,X,Y,dΩ)
    un, hn = uhn

    g = 9.81
    res(t,(u,h),(w,ϕ)) = ∫(∂t(u)⋅w -g*DIV(w)*h + ∂t(h)*ϕ -∇(ϕ)⋅((h+b)*u))dΩ
    jac(t,(u,h),(du,dh),(w,ϕ)) = ∫(-g*DIV(w)*dh -∇(ϕ)⋅((dh)*u + (h+b)*du))dΩ
    jac_t(t,(u,h),(dut,dht),(w,ϕ)) = ∫(dut⋅w + dht*ϕ)dΩ
    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = 
    ode_solver = ThetaMethod(nls,1,0.5)
    x = solve(ode_solver,op,uhn,0.0,100)
    dir = "swe-solver/1d-topo-output"
    if isdir(dir)
        output_file = paraview_collection("swe-solver/1d-topo-output")do pvd
            for (x,t) in x
                u,h = x
                pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"h"=>h])
                println("done $t/10")
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection("swe-solver/1d-topo-output") do pvd
            for (x,t) in x
                u,h = x
                pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"h"=>h])
                println("done $t/10")
            end
        end
    end
end

function h₀((x,y))
    h = 3 + 0.5*sin(π*x)*sin(π*y)
    h
end

function topography((x,y))
    p = 1.0
    p
end

function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

Shallow_water_theta_newton(1,3,h₀,u₀,topography,0.5)