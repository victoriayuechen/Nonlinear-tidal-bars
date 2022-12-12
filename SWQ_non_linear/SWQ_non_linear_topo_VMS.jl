# using Pkg
# Pkg.add("Gridap")
# Pkg.add("SparseMatricesCSR")
# Pkg.add("SparseArrays")
# Pkg.add("WriteVTK")
# Pkg.add("LineSearches")
# Pkg.add("LinearAlgebra")
# Pkg.activate(".")

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
        order,degree,h₀,u₀,topography,dt,Tend,
        linear_solver::Gridap.Algebra.LinearSolver=Gridap.Algebra.BackslashSolver(),
        sparse_matrix_type::Type{<:AbstractSparseMatrix}=SparseMatrixCSC{Float64,Int})
    #Create model
    B = 10
    L = 10
    dx = 1
    dy = 1
    g = 9.80655
    H = 0.5
    h = 0
    v = 10^(-6)

    #Stabilization parameters
    c₁ = 12
    c₂ = 2
    c₃ = 1

    domain = (0,B,0,L)
    
    partition = (100,100)
    dir = "./RESULTS/"
    model = CartesianDiscreteModel(domain,partition;isperiodic=(false,false))


    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,2,5])
    add_tag_from_tags!(labels,"left",[7])
    add_tag_from_tags!(labels,"right",[8])
    add_tag_from_tags!(labels,"top",[3,4,6])
    add_tag_from_tags!(labels,"inside",[9])
    DC = ["right","left"]

    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    dω = Measure(Ω,degree,ReferenceDomain())
    Γ = BoundaryTriangulation(model,tags=DC)
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)
    

    

    reffe_rt = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    V = TestFESpace(model,reffe_rt)
    U = TransientTrialFESpace(V)

    reffe_lgn = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lgn)
    P = TransientTrialFESpace(Q)


    RTMM,RTMMchol = setup_and_factorize_mass_matrices(dΩ,V,Q,U,P)

    Y = MultiFieldFESpace([V,Q,V])#∇u, ∇h
    X = TransientMultiFieldFESpace([U,P,U])
    E = [0 -1; 1 0]
    #Create initial solutions
    a1(u,v) = ∫(v⋅u)dΩ
    l1(v) = ∫(v⋅u₀)dΩ
    uₙ = solve(AffineFEOperator(a1,l1,U,V))

    a2(u,v) = ∫(v*u)dΩ
    l2(v) = ∫(v*h₀)dΩ
    ζₙ = solve(AffineFEOperator(a2,l2,P,Q))

    a3(u,v) = ∫(v*u)dΩ
    l3(v) = ∫(v*topography)dΩ
    b = solve(AffineFEOperator(a3,l3,P,Q))


    unv,hnv = get_free_dof_values(uₙ,ζₙ)
    F₀ = clone_fe_function(V,uₙ)
    compute_mass_flux!(F₀,dΩ,V,RTMMchol,uₙ*ζₙ)
    
    coriolis((x,y)) = [0 -1;1 0]
    uhn = uh(uₙ,ζₙ,F₀,X,Y,dΩ)
    uₙ, ζₙ, F = uhn
    A = [0 -1; 1 0]
    forcfunc(t) = VectorValue(0.5,0)  


    perp(u) = VectorValue(-u[2],u[1])
    norm(u) = sqrt(u⋅u)
    dnorm(u,du) = u ⋅ du / norm(u)
    I = [1,1]
    Rζ(u,ζ,b) = ∂t(ζ) + (ζ + b) * (∇ ⋅ (u)) + ∇(ζ)'⋅u
    Rᵤ(u,ζ,b) = ∂t(u) + ∇(u)'⋅u + cD * norm∘(u) * u / (ζ+b) + g * (∇(ζ)) #To be added, forcing function Fₚ ; + f*perp∘(u) : coriolis neglected
    dRζ(u,ζ,du,dζ,b) = dζ * (∇⋅(u)) + (ζ+b)*(∇⋅(du)) + du ⋅ (∇(ζ)) + u ⋅ (∇(ζ))
    dRᵤ(u,ζ,du,dζ,b) = ∇(u)'⋅du + ∇(du)'⋅u + cD * dnorm∘(u,du) * u / (ζ+b) + cD * norm∘(u) * du / (ζ+b) + cD * norm∘(u) * u * dζ /(ζ+b)*(ζ+b) +  g * ∇(dζ) # + f*perp∘(du) : coriolis neglected
    Lζ(v,w) = (h-H) * (∇⋅v) - (∇(w))⋅uₙ
    Lᵤ(v,w) = -∇(v)⋅uₙ - g * ∇(w) # - f * E * I : coriolis neglected
    τᵤ(a,ζ) = 1.0 / (c₁*ν/(Δxₒ.^2) + c₂*a/Δxₒ + c₃*cD*g*a/(ζ+1e-14))
    τζ(a,ζ) = Δxₒ.^2/(c₁*τᵤ(a,ζ))
    dτᵤdu(a,ζ,da) = - τᵤ(a,ζ)*τᵤ(a,ζ) * (c₂/Δxₒ + c₃*cD*g/(ζ+1e-14))*da
    dτᵤdζ(a,ζ,dζ) = τᵤ(a,ζ)*τᵤ(a,ζ) * c₃*cD*g*a/(ζ*ζ+1e-14)*dζ
    dτζdu(a,ζ,da) = τζ(a,ζ)/τᵤ(a,ζ)*dτᵤdu(a,ζ,da)
    dτζdζ(a,ζ,dζ) = τζ(a,ζ)/τᵤ(a,ζ)*dτᵤdζ(a,ζ,dζ)

    res(t, (u, ζ), (v, w)) = ∫(∂t(ζ)*w - (ζ + H - h)*u⋅(∇(w)) + (∂t(u) + (u⋅∇)*u + cD * norm∘(u)*u/(ζ + H - h))⋅v -g*ζ*(∇⋅v)    # Remember to add forcing function Fₚ
     - Rζ∘(u, ζ, (H-h))*(τζ∘(norm∘(u), ζ)*Lζ∘(v, w))
     - Rᵤ∘(u, ζ, (H-h))⋅(τᵤ∘(norm∘(u), ζ)*Lᵤ∘(v, w)))dΩ + ∫(g*ζ*v⋅nΓ)dΓ

    jac(t, (u, ζ), (du, dζ), (v, w)) = ∫(((ζ + H - h)*du + dζ)⋅∇(w) + ((du⋅∇)*u + (u⋅∇)*du + cD * dnorm∘(u, du) * u / (ζ+H-h) + cD*norm∘(u)/(ζ+H-h) * du + cD*norm∘(u)*u*dζ / (ζ+H-h)^2)⋅v - g*dζ*(∇⋅v)
     - dRζ∘(u,ζ,du,dζ,(H-h)) * τζ∘(norm∘(u), ζ)*Lζ∘(v, w) + Rζ∘(u, ζ, (H-h))*((dτζdu∘(norm∘(u), ζ, dnorm∘(u, du)) + dτζdζ∘(norm∘(u), ζ, dζ)) * Lζ∘(v, w))
     - dRᵤ∘(u,ζ,du,dζ,(H-h))⋅(τᵤ∘(norm∘(u), ζ)*Lᵤ∘(v, w)) + Rᵤ∘(u, ζ, (H-h))⋅((dτζdu∘(norm∘(u), ζ, dnorm∘(u, du)) + dτζdζ∘(norm∘(u), ζ, dζ)) * Lᵤ∘(v, w)))dΩ + ∫(g*dζ*v⋅nΓ)dΓ
    
    jac_t(t, (u, ζ), (dut, dζt), (v, w)) = ∫(dζt*w +dut⋅v - dζt*τζ∘(norm∘(u), ζ)* Lζ∘(v, w) - dut⋅(τᵤ∘(norm∘(u), ζ)*L∘(v, w)))dΩ

    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=true,linesearch=BackTracking())
    ode_solver = ThetaMethod(nls,dt,0.5)
    x = solve(ode_solver,op,uhn,0.0,Tend)

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
    h = -topography((x,y)) +  1  #+ 1*exp(-10*(x-2.5)^2 -10*(y-2.5)^2)
    h
end

function topography((x,y))
    p = 0.0 #0.8*exp(-5*(x-5)^2)
    p
end

function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end



#First a simple one to speed up the code --> tip from Oriol's meeting
Shallow_water_theta_newton(1,3,h₀,u₀,topography,0.1,0.2)

#Measurement 
Shallow_water_theta_newton(1,3,h₀,u₀,topography,5,60*60*24)