using Pkg
Pkg.activate(".")

using Gridap
using SparseMatricesCSR
using SparseArrays
using WriteVTK
using LinearAlgebra
using LineSearches: BackTracking
using GridapGmsh

#Solves linear shallow water equations on a 2d plane

function uh(u₀,h₀,F₀,q₀,X,Y,dΩ)
    a((u,h,r,u2),(w,ϕ,ϕ2,w2)) = ∫(w⋅u + ϕ*h + w2⋅u2 + ϕ2*r)dΩ
    b((w,ϕ,ϕ2,w2)) = ∫(w⋅u₀ + ϕ*h₀ + w2⋅F₀ + q₀*ϕ2)dΩ
    solve(AffineFEOperator(a,b,X,Y))
end

function compute_mass_flux!(F,dΩ,V,RTMMchol,u)
    b(v) = ∫(v⋅u)dΩ
    Gridap.FESpaces.assemble_vector!(b, get_free_dof_values(F), V)
    ldiv!(RTMMchol,get_free_dof_values(F))
end

function compute_potential_vorticity!(q,H1h,H1hchol,dΩ,R,S,h,u,f)
    a(r,s) = ∫(s*h*r)dΩ
    c(s)   = ∫(perp∘(∇(s))⋅(u) + s*f)dΩ
    Gridap.FESpaces.assemble_matrix_and_vector!(a, c, H1h, get_free_dof_values(q), R, S)
    lu!(H1hchol, H1h)
    ldiv!(H1hchol, get_free_dof_values(q))
  end

clone_fe_function(space,f)=FEFunction(space,copy(get_free_dof_values(f)))

function setup_and_factorize_mass_matrices(dΩ, R, S, U, V, P, Q)
    amm(a,b) = ∫(a⋅b)dΩ
    H1MM = assemble_matrix(amm, R, S)
    RTMM = assemble_matrix(amm, U, V)
    L2MM = assemble_matrix(amm, P, Q)
    H1MMchol = lu(H1MM)
    RTMMchol = lu(RTMM)
    L2MMchol = lu(L2MM)
  
    H1MM, RTMM, L2MM, H1MMchol, RTMMchol, L2MMchol
  end

function new_vtk_step(Ω,file,_cellfields)
    n = size(_cellfields)[1]
    createvtk(Ω,
              file,
              cellfields=_cellfields,
              nsubcells=n)
end

perp(u) = VectorValue(-u[2],u[1])
coriolis(u) = 0.5*perp(u)

function Gridap.get_free_dof_values(functions...)
    map(get_free_dof_values,functions)
  end


function Shallow_water_theta_newton(
        order,degree,h₀,u₀,topography,f,
        linear_solver::Gridap.Algebra.LinearSolver=Gridap.Algebra.BackslashSolver(),
        sparse_matrix_type::Type{<:AbstractSparseMatrix}=SparseMatrixCSC{Float64,Int})
    #Create model
    B = 10
    L = 10
    dx = 1
    dy = 1
    latitude = 52
    η = 7.29e-5
    #f = 0.05#2*η*sin(latitude*(π/180))
    g = 9.81

    #Domain properties
    domain = (0,B,0,L)
    partition = (50,50)

    model = CartesianDiscreteModel(domain,partition;isperiodic=(false,false))

    #Make labels
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,2,5])
    add_tag_from_tags!(labels,"left",[7])
    add_tag_from_tags!(labels,"right",[8])
    add_tag_from_tags!(labels,"top",[3,4,6])
    add_tag_from_tags!(labels,"inside",[9])
    DC = ["right","left"]

    DC = ["left","right"]

    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    dω = Measure(Ω,degree,ReferenceDomain())
    Γ = BoundaryTriangulation(model,tags=DC)
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)
    
    

    reffe_rt = ReferenceFE(raviart_thomas,Float64,order)
    V = TestFESpace(model,reffe_rt,dirichlet_tags="boundary")
    U = TrialFESpace(V)

    reffe_lgn = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lgn)
    P = TrialFESpace(Q)

    reffe_lgn = ReferenceFE(lagrangian, Float64, order+1)
    S = FESpace(model, reffe_lgn)
    R = TrialFESpace(S)


    H1MM, RTMM, L2MM, H1MMchol, RTMMchol, L2MMchol = setup_and_factorize_mass_matrices(dΩ,R,S,U,V,P,Q)

    Y = MultiFieldFESpace([V,Q,S,V])
    X = MultiFieldFESpace([U,P,R,U])

    E = [0 -1; 1 0]
    #Create initial solutions
    a1(u,v) = ∫(v⋅u)dΩ
    l1(v) = ∫(v⋅u₀)dΩ
    un = solve(AffineFEOperator(a1,l1,U,V))

    a2(u,v) = ∫(v*u)dΩ
    l2(v) = ∫(v*h₀)dΩ
    hn = solve(AffineFEOperator(a2,l2,P,Q))


    b = interpolate(topography,P)

    a4(u,v)=∫(v*u)dΩ
    l4(v)=∫(v*f)dΩ
    fn=solve(AffineFEOperator(a4,l4,R,S))



    q₀ = clone_fe_function(P,fn)
    compute_potential_vorticity!(q₀,H1MM,H1MMchol,dΩ,R,S,hn,un,fn)


    unv,hnv = get_free_dof_values(un,hn)
    F₀ = clone_fe_function(V,un)
    compute_mass_flux!(F₀,dΩ,V,RTMMchol,un*hn)
    
    
    uhn = uh(un,hn,F₀,q₀,X,Y,dΩ)
    un,hn,q,F= uhn

    #perp(u) = VectorValue(-u[2],u[1])

    forcfunc(t) = VectorValue(0.5*cos((1/10)*π*t),0.0)  

    N = 100
    θ = 0.5
    dt = 0.1
    τ=0.5*dt
    g = 9.81
    nlrtol = 1.0e-8
    dir = "swe-solver/1d-topo-output_zero"

    for step=1:N
        function res((Δu,Δh,qvort,F),(v,q,s,v2))
            one_m_θ = (1-θ)
            uiΔu    = un     + one_m_θ*Δu
            hiΔh    = hn     + one_m_θ*Δh
            hbiΔh   = hn + b + one_m_θ*Δh 
            ∫((1.0/dt)*v⋅(Δu)-(∇⋅(v))*(g*hbiΔh + 0.5*uiΔu⋅uiΔu) + ((qvort-τ*(uiΔu⋅∇(qvort)))*(perp∘(F)))⋅v + (1.0/dt)*q*(Δh))dΩ + ∫(q*(DIV(F)))dω + ∫(s*qvort*hiΔh +  perp∘(∇(s))⋅uiΔu - s*fn + v2⋅(F-hiΔh*uiΔu))dΩ   

        end
        
        function jac((Δu,Δh,qvort,F),(du,dh,dq,dF),(v,q,s,v2))
            one_m_θ = (1-θ)
            uiΔu  = un + one_m_θ*Δu
            hiΔh  = hn + one_m_θ*Δh
            uidu  = one_m_θ*du
            hidh  = one_m_θ*dh
            ∫((1.0/dt)*v⋅du +  (dq  - τ*(uiΔu⋅∇(dq)+uidu⋅∇(qvort)))*(perp∘(F))⋅v +  (qvort - τ*(uiΔu⋅∇(qvort)))*(perp∘(dF))⋅v -  (∇⋅(v))*(g*hidh +uiΔu⋅uidu) + (1.0/dt)*q*dh)dΩ + ∫(q*(DIV(dF)))dω  + ∫(s*(qvort*hidh+dq*hiΔh) + perp∘(∇(s))⋅uidu + v2⋅(dF-hiΔh*uidu-hidh*uiΔu))dΩ 
        end
        dY = get_fe_basis(Y)
        residualuhn=res(uhn,dY)
        r=assemble_vector(residualuhn,Y)
        assem = SparseMatrixAssembler(sparse_matrix_type,Vector{Float64},X,Y)
        op = FEOperator(res,jac,X,Y,assem)
        nls=NLSolver(show_trace=true)
        solver=FESolver(nls)
        solve!(uhn,solver,op)
        
        unv .=unv+get_free_dof_values(un)
        hnv .= hnv+get_free_dof_values(hn)
        pvd[dt*Float64(step)] = createvtk(Ω,joinpath(dir,"nswe_cells_linear_shallow_water_$(dt*step)"*".vtu"),cellfields=["u"=>un,"h"=>hn])
    end
end

function h₀((x,y))
    h = -topography((x,y)) +  1  + 0.1*exp(-100*(x-5)^2 -100*(y-5)^2)
    h
end

function topography((x,y))
    p = 0.0
    p
end

function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

function f((x,y))
    f = 0.05
    f
end

Shallow_water_theta_newton(1,3,h₀,u₀,topography,f)