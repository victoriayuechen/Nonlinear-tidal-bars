module MyNonlinearAPVM
using Pkg
Pkg.activate(".")

using Gridap
using SparseMatricesCSR
using SparseArrays
using GridapPardiso
using WriteVTK
using LinearAlgebra
using LineSearches: BackTracking
using GridapGmsh
using Gridap.TensorValues: meas
#Solves nonlinear shallow water equations on a 2d plane
#Makes use of the Anticipated Potential Vorticity Method (APVM) produced by McRae and Cotter: https://doi.org/10.1002/qj.2291 
#Implementation in Gridap for this 2D domain also uses similar code produced by the Gridap GeoSciences github: https://github.com/gridapapps/GridapGeosciences.jl

export APVM_run

function uD(u₀,D₀,F₀,q₀,X,Y,dΩ)
    a((u,D,r,u2),(w,ϕ,ϕ2,w2)) = ∫(w⋅u + ϕ*D + w2⋅u2 + ϕ2*r)dΩ
    b((w,ϕ,ϕ2,w2)) = ∫(w⋅u₀ + ϕ*D₀ + w2⋅F₀ + q₀*ϕ2)dΩ
    solve(AffineFEOperator(a,b,X(0.0),Y))
end

function compute_mass_flux(F,U,V,dΩ)
    a(DU,w) = ∫(DU⋅w)dΩ
    b(w) = ∫(w⋅F)dΩ
    solve(AffineFEOperator(a,b,U(0.0),V))
end

function compute_potential_vorticity(D,f,u,R,S,dΩ)
    a(r,s) = ∫(s*D*r)dΩ
    c(s)   = ∫(perp∘(∇(s))⋅(u) + s*f)dΩ
    solve(AffineFEOperator(a,c,R(0.0),S))
end

clone_fe_function(space,f)=FEFunction(space,copy(get_free_dof_values(f)))


function new_vtk_step(Ω,file,_cellfields)
    n = size(_cellfields)[1]
    createvtk(Ω,
              file,
              cellfields=_cellfields,
              nsubcells=n)
end

#Perpendicular operator
perp(u) = VectorValue(-u[2],u[1])

function Gridap.get_free_dof_values(functions...)
    map(get_free_dof_values,functions)
end


function APVM_run(order,degree,D₀,u₀,topography,forcefunc,dir,periodic::Bool,Tend,dt,model,DC)
    #Parameters
    latitude = 50 #Latitude of the model being analysed
    η = 7.29e-5
    f = 2*η*sin(latitude*(π/180))
    cd = 0.0025
    g = 9.81
    T0 = 0.0
    τ = dt*0.5
    


    #Make triangulations and boundaries
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    dω = Measure(Ω,degree,ReferenceDomain())
    if periodic
        Γ = BoundaryTriangulation(model,tags=DC)
    else
        Γ = BoundaryTriangulation(model,tags=DC)
    end
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)

    udc(x,t::Real) = VectorValue(0.0,0.0)
    udc(t::Real) = x -> udc(x,t)


    #Make FE spaces
    if periodic
        reffe_rt = ReferenceFE(raviart_thomas,Float64,order)#
        V = TestFESpace(model,reffe_rt;conformity=:HDiv,dirichlet_tags=DC)#
        U = TransientTrialFESpace(V,[udc,udc])
    else
        reffe_rt = ReferenceFE(raviart_thomas,Float64,order)
        V = TestFESpace(model,reffe_rt;conformity=:HDiv,dirichlet_tags=DC)#
        U = TransientTrialFESpace(V,[udc,udc])
    end

    reffe_lgn = ReferenceFE(lagrangian,Float64,order-1)
    Q = TestFESpace(model,reffe_lgn;conformity=:L2)#
    P = TransientTrialFESpace(Q)

    reffe_lgn = ReferenceFE(lagrangian, Float64, order)
    S = TestFESpace(model, reffe_lgn;conformity=:H1)#
    R = TransientTrialFESpace(S)


    Y = MultiFieldFESpace([V,Q,S,V])
    X = TransientMultiFieldFESpace([U,P,R,U])

    #Create initial solutions
    #u, velocity
    un = interpolate_everywhere(u₀,U(0.0))

    #h, fluid depth
    Dn = interpolate_everywhere(D₀,P(0.0))

    #b, topography
    h = interpolate_everywhere(topography,Q(0.0))

    #Constant FE space with the coriolis parameter, possibly can be removed.
    fn = interpolate_everywhere(f,S(0.0))

    #Build and compute initial potential vorticity
    q₀ = compute_potential_vorticity(D₀,f,u₀,R,S,dΩ)

    #Build and compute initial mass flux
    F₀ = compute_mass_flux(un*Dn,U,V,dΩ)

    
    #Solve whole coupled system with the initial values
    uDn = uD(un,Dn,F₀,q₀,X,Y,dΩ)
    un,Dn,q,F= uDn

    #Forcing function on u(t)
    forcfunc(t) = x -> forcefunc(x,t)

    #Norm operator
    norm(u) = meas∘(u) + 1e-14
    dnorm(u,du) = u ⋅ du / norm(u)

    #Define residual, jacobian in space and time
    res(t,(u,D,q,F),(w,ϕ,ϕ2,w2)) = ∫(∂t(u)⋅w + ((q-τ*(u⋅∇(q)))*(perp∘(F)))⋅w - (∇⋅(w))*(g*(D+h) + 0.5*(u⋅u)) + ∂t(D)*ϕ + w2⋅(F-u*D) + ϕ2*(q*D) +  (perp∘(∇(ϕ2)))⋅u - (ϕ2*fn)- forcfunc(t)⋅w +  ϕ*(∇⋅(F)) + ((cd*(norm(u))*u)/(D+1e-14))⋅w)dΩ + ∫((g*(D+h) + 0.5*(u⋅u))*(w⋅nΓ) + nΓ⋅(perp∘(u))*ϕ2)dΓ
    jac(t,(u,D,q,F),(du,dD,dq,dF),(w,ϕ,ϕ2,w2)) = ∫(((dq - τ*(u⋅∇(dq) + du⋅∇(q)))*(perp∘(F)))⋅w + ((q-τ*(u⋅∇(q)))*(perp∘(dF)))⋅w - (∇⋅(w))*(g*(dD) + (du⋅u)) + w2⋅(dF - du*D - dD*u) + ϕ2*(q*dD) + ϕ2*(dq*D) + du⋅(perp∘(∇(ϕ2))) +  ϕ*(∇⋅(dF))  + (cd*(dnorm(u,du)*u)/(D+h+1e-14) + cd*(norm(u)/(D+h+1e-14)*du) - cd*(norm(u)*u*dD)/((D+h+1e-14)*(D+h+1e-14)))⋅w )dΩ+ ∫((g*(dD) + (du⋅u))*(w⋅nΓ) + nΓ⋅(perp∘(du))*ϕ2)dΓ#
    jac_t(t,(u,D),(dut,dDt),(w,ϕ)) = ∫(dut⋅w + dDt*ϕ)dΩ

    #Define operators and solvers
    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=true,method =:newton,linesearch = BackTracking())
    ode_solver = ThetaMethod(nls,dt,0.5)
    x = solve(ode_solver,op,uDn,T0,Tend)

    probe = [Point(2500,500),Point(5000,500)]
    

    #Output results
    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"nonlinear_topo"))do pvd
            pvd[0.0] = createvtk(Ω,joinpath(dir,"nonlinear_topo0.0.vtu"),cellfields=["u"=>un,"D"=>(Dn+h),"h"=>h])
            for (x,t) in x
                u,D,F = x
                println(u(probe))
                println(D(probe))
                pvd[t] = createvtk(Ω,joinpath(dir,"nonlinear_topo$t.vtu"),cellfields=["u"=>u,"D"=>(D+h),"h"=>h])
                println("done $t/$Tend")
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"nonlinear_topo")) do pvd
            pvd[0.0] = createvtk(Ω,joinpath(dir,"nonlinear_topo0.0.vtu"),cellfields=["u"=>un,"D"=>(Dn+h),"h"=>h])
            for (x,t) in x
                u,D,F = x
                pvd[t] = createvtk(Ω,joinpath(dir,"nonlinear_topo$t.vtu"),cellfields=["u"=>u,"D"=>(D+h),"h"=>h])
                println("done $t/$Tend")
            end
        end
    end
end
end



