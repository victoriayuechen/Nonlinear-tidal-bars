using Pkg
Pkg.activate(".")

using Gridap
using SparseMatricesCSR
using SparseArrays
using WriteVTK
using LinearAlgebra
using LineSearches: BackTracking

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

clone_fe_function(space,f)=FEFunction(space,copy(get_free_dof_values(f)))

function setup_and_factorize_mass_matrices(dΩ,V,Q,U,P)
    amm(a,b) = ∫(a⋅b)dΩ

    RTMM = assemble_matrix(amm,U,V)

    RTMMchol = lu(RTMM)

    RTMM,RTMMchol
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
        order,degree,h₀,u₀,topography,
        linear_solver::Gridap.Algebra.LinearSolver=Gridap.Algebra.BackslashSolver(),
        sparse_matrix_type::Type{<:AbstractSparseMatrix}=SparseMatrixCSC{Float64,Int})
    #Create model
    B = 10
    L = 10
    dx = 1
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
    nuemanntaggs = ["left","right","top","bottom"]

    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    dω = Measure(Ω,degree,ReferenceDomain())
    Γ = BoundaryTriangulation(model,tags="boundary")
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)
    

    

    reffe_rt = ReferenceFE(raviart_thomas,Float64,order)
    V = TestFESpace(model,reffe_rt;conformity=:HDiv,dirichlet_tags="boundary")
    U = TransientTrialFESpace(V)

    reffe_lgn = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lgn;conformity=:L2)
    P = TransientTrialFESpace(Q)


    RTMM,RTMMchol = setup_and_factorize_mass_matrices(dΩ,V,Q,U,P)

    Y = MultiFieldFESpace([V,Q,V])#∇u, ∇h
    X = TransientMultiFieldFESpace([U,P,U])

    E = [0 -1; 1 0]
    #Create initial solutions
    a1(u,v) = ∫(v⋅u)dΩ
    l1(v) = ∫(v⋅u₀)dΩ
    un = solve(AffineFEOperator(a1,l1,U,V))

    a2(u,v) = ∫(v*u)dΩ
    l2(v) = ∫(v*h₀)dΩ
    hn = solve(AffineFEOperator(a2,l2,P,Q))

    a3(u,v) = ∫(v*u)dΩ
    l3(v) = ∫(v*topography)dΩ
    b = solve(AffineFEOperator(a3,l3,P,Q))
    unv,hnv = get_free_dof_values(un,hn)
    F₀ = clone_fe_function(V,un)
    compute_mass_flux!(F₀,dΩ,V,RTMMchol,un*hn)
    
    coriolis((x,y)) = [0 -1;1 0]
    uhn = uh(un,hn,F₀,X,Y,dΩ)
    un, hn,F = uhn
    A = [0 -1; 1 0]
    forcfunc(t) = VectorValue(0.5*π*cos(π*t),0)  

    g = 9.81
    res(t,(u,h,F),(w,ϕ,w2)) = ∫(∂t(u)⋅w -g*(∇⋅(w))*(b+h)  + ∂t(h)*ϕ  + w2⋅(F - u*h))dΩ + ∫(ϕ*DIV(F))dω + ∫(g*(h+b)*(w⋅nΓ) )dΓ
    jac(t,(u,h,F),(du,dh,dF),(w,ϕ,w2)) = ∫(-g*(∇⋅(w))*dh  + w2⋅(dF -du*h -u*dh))dΩ + ∫(ϕ*DIV(dF))dω  + ∫(g*dh*(w⋅nΓ))dΓ
    jac_t(t,(u,h),(dut,dht),(w,ϕ)) = ∫(dut⋅w + dht*ϕ)dΩ


    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=true,linesearch=BackTracking())
    Tend = 5
    ode_solver = ThetaMethod(nls,0.1,0.5)
    x = solve(ode_solver,op,uhn,0.0,Tend)
    dir = "swe-solver/1d-topo-output_zero"
    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"1d-topo-output"))do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"1d-topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn+b),"b"=>b])
            for (x,t) in x
                u,h,F = x
                pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"h"=>(h+b),"b"=>b])
                println("done $t/$Tend")
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"1d-topo-output")) do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"1d-topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn+b),"b"=>b])
            for (x,t) in x
                u,h,F = x
                pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"h"=>(h+b),"b"=>b])
                println("done $t/$Tend")
            end
        end
    end
end

function h₀((x,y))
    h = -topography((x,y)) +  1  + 1*exp(-10*(x-2.5)^2 -10*(y-2.5)^2)
    h
end

function topography((x,y))
    p = 0.8*exp(-5*(x-5)^2)
    p
end

function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

Shallow_water_theta_newton(1,3,h₀,u₀,topography)