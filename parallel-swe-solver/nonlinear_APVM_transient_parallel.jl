using Pkg
Pkg.activate(".")

using Gridap
using SparseMatricesCSR
using SparseArrays
using LinearAlgebra
using LineSearches: BackTracking
using Gridap.TensorValues: meas
using GridapDistributed
using PartitionedArrays
using WriteVTK
using LineSearches: BackTracking
using GridapPETSc
using GridapPETSc: PETSC

B = 100
L = 100

#Variable functions to be used to setup model, used for local tests
function D₀((x,y))
    Dout = -topography((x,y)) +  0.5 + 0.01*exp(-0.01(x-50)^2 - 0.01*(y-30)^2)
    Dout
end

function topography((x,y))
    p =  0.1 * cos(π/B * x) * cos(2π/L * y)
    p
end

function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

function forcefunc((x,y),t)
    U_start =  0.5 
    η = 7.29e-5
    latitude = 52
    cD = 0.0025
    f = 2*η*sin(latitude*(π/180))
    σ = 2*pi/44700
    f = VectorValue(0.0,0.0)#VectorValue(-f*U_start*cos(σ*t),-σ*U_start*sin(σ*t)+cD*abs(U_start*cos(σ*t))*U_start*cos(σ*t)) 
    f
end


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

#Perpendicular operator
perp(u) = VectorValue(-u[2],u[1])

function Gridap.get_free_dof_values(functions...)
    map(get_free_dof_values,functions)
end


options = """ 
    -snes_type newtonls
    -snes_linesearch_type basic
    -snes_linesearch_damping 1.0
    -snes_rtol 1.0e-8
    -snes_atol 0.0
    -snes_monitor
    -snes_converged_reason
    -mm_ksp_type cg
    -mm_ksp_monitor
    -mm_ksp_rtol 1.0e-10
    -mm_pc_type jacobi
"""

function apvm_swe(parts)
    GridapPETSc.with(args=split(options)) do
        mesh_partition = (50,50)
        domain = (0,B,0,L)
        model = CartesianDiscreteModel(domain,mesh_partition;isperiodic=(false,true)) 
        labels = get_face_labeling(model)
        add_tag_from_tags!(labels,"bottom",[1,2,5])
        add_tag_from_tags!(labels,"left",[7])
        add_tag_from_tags!(labels,"right",[8])
        add_tag_from_tags!(labels,"top",[3,4,6])
        add_tag_from_tags!(labels,"inside",[9])
        DC = ["left","right"]
        Tend = 30
        order = 1 
        degree = 4 
        dir = "output2"
        periodic=true
        # Parameters
        latitude = 52 #Latitude of the model being analysed
        η = 7.29e-5
        f = 2*η*sin(latitude*(π/180))
        cd = 0.0025
        g = 9.81
        T0 = 0.0
        dt = 0.5
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
        uda(x,t::Real) = VectorValue(1*sin(π*(1/3600)*t),0.0)
        udc(t::Real) = x -> udc(x,t)
        uda(t::Real) = x -> uda(x,t)

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

        reffe_lgn = ReferenceFE(lagrangian,Float64,order)
        Q = TestFESpace(model,reffe_lgn;conformity=:L2)#
        P = TransientTrialFESpace(Q)

        reffe_lgn = ReferenceFE(lagrangian, Float64, order+1)
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
        res(t,(u,D,q,F),(w,ϕ,ϕ2,w2)) = ∫(∂t(u)⋅w + ((q-τ*(u⋅∇(q)))*(perp∘(F)))⋅w - (∇⋅(w))*(g*(D+h) + 0.5*(u⋅u)) + ∂t(D)*ϕ + w2⋅(F-u*D) + ϕ2*(q*D) +  (perp∘(∇(ϕ2)))⋅u - (ϕ2*fn) + ϕ*(∇⋅(F))- forcfunc(t)⋅w + ((cd*(norm(u))*u)/(D+1e-14))⋅w)dΩ + ∫((g*(D+h) + 0.5*(u⋅u))*(w⋅nΓ) + nΓ⋅(perp∘(u))*ϕ2)dΓ 
        jac(t,(u,D,q,F),(du,dD,dq,dF),(w,ϕ,ϕ2,w2)) = ∫(((dq - τ*(u⋅∇(dq) + du⋅∇(q)))*(perp∘(F)))⋅w + ((q-τ*(u⋅∇(q)))*(perp∘(dF)))⋅w - (∇⋅(w))*(g*(dD) + (du⋅u)) + w2⋅(dF - du*D - dD*u) + ϕ2*(q*dD) + ϕ2*(dq*D) + du⋅(perp∘(∇(ϕ2))) + ϕ*(∇⋅(dF)) + (cd*(dnorm(u,du)*u)/(D+h+1e-14) + cd*(norm(u)/(D+h+1e-14)*du) - cd*(norm(u)*u*dD)/((D+h+1e-14)*(D+h+1e-14)))⋅w )dΩ+ ∫((g*(dD) + (du⋅u))*(w⋅nΓ) + nΓ⋅(perp∘(du))*ϕ2)dΓ
        jac_t(t,(u,D),(dut,dDt),(w,ϕ)) = ∫(dut⋅w + dDt*ϕ)dΩ

        #Define operators and solvers
        op = TransientFEOperator(res,jac,jac_t,X,Y)
        nls =NLSolver(show_trace=true,method =:newton,linesearch=BackTracking())
        ode_solver = ThetaMethod(nls,dt,0.5)
        x = solve(ode_solver,op,uDn,T0,Tend)

        #Output results
        if isdir(dir)
            createpvd(joinpath(dir,"non_linear_topo")) do pvd
                for (x,t) in x
                    u,D,F = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"nonlinear_topo$t.vtu"),cellfields=["u"=>u,"D"=>(D+h),"h"=>h])
                end
            end
        else
            mkdir(dir)
            createpvd(joinpath(dir,"non_linear_topo")) do pvd
                for (x,t) in x
                    u,D,F = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"nonlinear_topo$t.vtu"),cellfields=["u"=>u,"D"=>(D+h),"h"=>h])                
                end
            end
        end
    end
end


partition = (2, 2)
prun(apvm_swe, mpi, partition) # use `sequential` instead of `mpi` for debugging 