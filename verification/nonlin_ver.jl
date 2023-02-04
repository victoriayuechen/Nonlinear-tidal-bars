using Pkg
Pkg.activate(".")

using Gridap
using SparseMatricesCSR
using SparseArrays
using WriteVTK
using LinearAlgebra
using LineSearches: BackTracking
using GridapGmsh
using Gridap.TensorValues: meas

U_start =  0.5 
η = 7.29e-5
latitude = 50
cD = 0.0025
f = 2*η*sin(latitude*(π/180))
T = 10
σ = 2*pi/T
H = 0.5

dt = 0.1


B = 2 #y
L = 10 #x
partition = (80,16)
domain = (0,L,0,B)
#model = GmshDiscreteModel("swe-solver/meshes/1x1periodic01.msh")
model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false)) 
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"bottom",[1,2,5])
add_tag_from_tags!(labels,"left",[7])
add_tag_from_tags!(labels,"right",[8])
add_tag_from_tags!(labels,"top",[3,4,6])
add_tag_from_tags!(labels,"inside",[9])
DC = ["bottom","top"]



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

function D₀((x,y))
    Dout = -topography((x,y)) +  H #+ 0.01*exp(-0.01(x-0.5)^2 - 0.01*(y-5)^2)
    Dout
end

function topography((x,y))
    p =  0.0#0.1 * cos(π/B * y) * cos(2π/L * x)
    p
end

function u₀((x,y))
    u = VectorValue(1.0,0.0)
    u
end

function forcfunc((x,y),t)
    U_start =  0.5 
    η = 7.29e-5
    latitude = 50
    cD = 0.0025
    f = 2*η*sin(latitude*(π/180))
    σ = 2*pi/T
    H = 3
    f = VectorValue(0.0,0.0)#VectorValue(σ*U_start*cos(σ*t)+cD/H*abs(U_start*sin(σ*t))*U_start*sin(σ*t),f*U_start*sin(σ*t))
    f
end



outputdir = "output_swe"
dir = "verification/output_nnl"
# if !isdir(outputdir)
#     mkdir(outputdir)
# end

# if !isdir(joinpath(outputdir,dir))
#     mkdir(joinpath(outputdir,dir))
# end


function APVM_run(order,degree,D₀,u₀,topography,forcefunc,dir,periodic::Bool,T,dt,model,DC)
    #Parameters
    latitude = 50 #Latitude of the model being analysed
    η = 7.29e-5
    f = 2*η*sin(latitude*(π/180))
    #f = 0
    cd = 0.0025
    g = 9.81
    T0 = 0.0
    τ = dt*0.5
    cd = cd
    T = 10
    dt = 0.1
    


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

    #Forcing function on u(t)
    forcfunc(t) = x -> forcefunc(x,t)

    # manufactured solutions

    H = 0.5
    u0 = 0.0#1.0
    top = 0.0

    A3 = (2*π/T)
    A1 = π/(1*L)
    A2 = π/(B)

    B3 = (2*π/T)
    B1 = π/(1*L)
    B2 = π/B

    C3 = (2*π/T)
    C1 = π/(1*L)
    C2 = π/B

    Sh = 0.05
    Su1 = 0.1
    Su2 = 0.1

    um(x,t::Real) = VectorValue(u0+Su1*sin(B3*t)*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2)), Su2*sin(C3*t)*cos(C1*(x[1]-L/2))*sin(C2*x[2]))
    um(t::Real) = x -> um(x,t)

    ζm(x,t::Real) = Sh*sin(A2*x[2])*cos(A1*(x[1]-L/2))*sin(A3*t)
    ζm(t::Real) = x -> ζm(x,t)

    Dm(x,t::Real) = Sh*sin(A2*x[2])*cos(A1*(x[1]-L/2))*sin(A3*t) + H - top
    Dm(t::Real) = x -> Dm(x,t)

    fh(x,t::Real) = C2*Su2*cos(C1*(x[1]-L/2))*cos(C2*x[2])*(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)*sin(C3*t)+A2*Sh*Su2*cos(A1*(x[1]-L/2))*cos(C1*(x[1]-L/2))*cos(A2*x[2])*sin(C2*x[2])*sin(A3*t)*sin(C3*t)-A1*Sh*sin(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)*(Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t)+u0)-B1*Su1*sin(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)*sin(B3*t)+A3*Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*cos(A3*t)
    fh(t::Real) = x -> fh(x,t)

    fu(x,t::Real) = VectorValue(-Su2*cos(C1*(x[1]-L/2))*sin(C2*x[2])*(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)*sin(C3*t)*((-C1*Su2*sin(C1*(x[1]-L/2))*sin(C2*x[2])*sin(C3*t)+
    B2*Su1*cos(B1*(x[1]-L/2))*sin(B2*(x[2]-B/2))*sin(B3*t)+f)/(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)-
    τ*((Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*
    sin(B3*t)+u0)*((A1*Sh*sin(A1*(x[1]-L/2))*
    sin(A2*x[2])*sin(A3*t)*
    (-C1*Su2*sin(C1*(x[1]-L/2))*
    sin(C2*x[2])*sin(C3*t)+
    B2*Su1*cos(B1*(x[1]-L/2))*sin(B2*(x[2]-B/2))*sin(B3*t)+f))/(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)^2+
    (-C1^2*Su2*cos(C1*(x[1]-L/2))*sin(C2*x[2])*sin(C3*t)-B1*B2*Su1*sin(B1*(x[1]-L/2))*sin(B2*(x[2]-B/2))*sin(B3*t))/(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H))+
    Su2*cos(C1*(x[1]-L/2))*sin(C2*x[2])*sin(C3*t)*((B2^2*Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t)-C1*C2*Su2*sin(C1*(x[1]-L/2))*cos(C2*x[2])*sin(C3*t))/(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)-
    (A2*Sh*cos(A1*(x[1]-L/2))*cos(A2*x[2])*sin(A3*t)*(-C1*Su2*sin(C1*(x[1]-L/2))*sin(C2*x[2])*sin(C3*t)+
    B2*Su1*cos(B1*(x[1]-L/2))*sin(B2*(x[2]-B/2))*sin(B3*t)+f))/(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)^2)))+
    (cD*(Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t)+u0)*sqrt(Su2^2*cos(C1*(x[1]-L/2))^2*sin(C2*x[2])^2*sin(C3*t)^2+
    (Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t)+u0)^2))/(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)+
    0.5*(-2*C1*Su2^2*cos(C1*(x[1]-L/2))*sin(C1*(x[1]-L/2))*sin(C2*x[2])^2*sin(C3*t)^2-
    2*B1*Su1*sin(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t)*(Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t)+u0))+
    B3*Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*cos(B3*t)-A1*Sh*g*sin(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t),
    (Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)*(Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t)+u0)*((-C1*Su2*sin(C1*(x[1]-L/2))*sin(C2*x[2])*sin(C3*t)+B2*Su1*cos(B1*(x[1]-L/2))*sin(B2*(x[2]-B/2))*sin(B3*t)+f)/(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)-τ*((Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t)+u0)*((A1*Sh*sin(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)*(-C1*Su2*sin(C1*(x[1]-L/2))*sin(C2*x[2])*sin(C3*t)+B2*Su1*cos(B1*(x[1]-L/2))*sin(B2*(x[2]-B/2))*sin(B3*t)+f))/(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)^2+(-C1^2*Su2*cos(C1*(x[1]-L/2))*sin(C2*x[2])*sin(C3*t)-B1*B2*Su1*sin(B1*(x[1]-L/2))*sin(B2*(x[2]-B/2))*sin(B3*t))/(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H))+Su2*cos(C1*(x[1]-L/2))*sin(C2*x[2])*sin(C3*t)*((B2^2*Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t)-C1*C2*Su2*sin(C1*(x[1]-L/2))*cos(C2*x[2])*sin(C3*t))/(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)-(A2*Sh*cos(A1*(x[1]-L/2))*cos(A2*x[2])*sin(A3*t)*(-C1*Su2*sin(C1*(x[1]-L/2))*sin(C2*x[2])*sin(C3*t)+B2*Su1*cos(B1*(x[1]-L/2))*sin(B2*(x[2]-B/2))*sin(B3*t)+f))/(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)^2)))+(Su2*cD*cos(C1*(x[1]-L/2))*sin(C2*x[2])*sin(C3*t)*sqrt(Su2^2*cos(C1*(x[1]-L/2))^2*sin(C2*x[2])^2*sin(C3*t)^2+(Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t)+u0)^2))/(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)+0.5*(2*C2*Su2^2*cos(C1*(x[1]-L/2))^2*cos(C2*x[2])*sin(C2*x[2])*sin(C3*t)^2-2*B2*Su1*cos(B1*(x[1]-L/2))*sin(B2*(x[2]-B/2))*sin(B3*t)*(Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t)+u0))+C3*Su2*cos(C1*(x[1]-L/2))*sin(C2*x[2])*cos(C3*t)+A2*Sh*g*cos(A1*(x[1]-L/2))*cos(A2*x[2])*sin(A3*t))
    fu(t::Real) = x -> fu(x,t)

    fD(x,t::Real) = Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H
    fD(t::Real) = x -> fD(x,t)

    fq(x,t::Real) = (-C1*Su2*sin(C1*(x[1]-L/2))*sin(C2*x[2])*sin(C3*t)+B2*Su1*cos(B1*(x[1]-L/2))*sin(B2*(x[2]-B/2))*sin(B3*t)+f)/(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)
    fq(t::Real) = x -> fq(x,t)

    fF(x,t::Real) = VectorValue((Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)*(Su1*cos(B1*(x[1]-L/2))*cos(B2*(x[2]-B/2))*sin(B3*t)+u0),
    Su2*cos(C1*(x[1]-L/2))*sin(C2*x[2])*(Sh*cos(A1*(x[1]-L/2))*sin(A2*x[2])*sin(A3*t)-top+H)*sin(C3*t))
    fF(t::Real) = x -> fF(x,t)



    #Make FE spaces
    if periodic
        reffe_rt = ReferenceFE(raviart_thomas, Float64,order)#
        V = TestFESpace(model,reffe_rt;conformity=:HDiv,dirichlet_tags=DC)#
        #U = TransientTrialFESpace(V,[udc,udc])
        U = TransientTrialFESpace(V,[um,um])
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
    #un = interpolate_everywhere(u₀,U(0.0))
    un = interpolate_everywhere(um(0.0),U(0.0))

    #h, fluid depth
    #Dn = interpolate_everywhere(D₀,P(0.0))
    Dn = interpolate_everywhere(Dm(0.0),P(0.0))

    #b, topography
    h = interpolate_everywhere(topography,Q(0.0))

    #Constant FE space with the coriolis parameter, possibly can be removed.
    fn = interpolate_everywhere(f,S(0.0))

    #Build and compute initial potential vorticity
    #q₀ = compute_potential_vorticity(D₀,f,u₀,R,S,dΩ)
    q₀ = compute_potential_vorticity(Dm(0.0),f,um(0.0),R,S,dΩ)

    #Build and compute initial mass flux
    F₀ = compute_mass_flux(un*Dn,U,V,dΩ)

    
    #Solve whole coupled system with the initial values
    uDn = uD(un,Dn,F₀,q₀,X,Y,dΩ)
    un,Dn,q,F= uDn

    #Norm operator
    norm(u) = meas∘(u) + 1e-14
    dnorm(u,du) = u ⋅ du / norm(u)



    #Define residual, jacobian in space and time
    res(t,(u,D,q,F),(w,ϕ,ϕ2,w2)) = ∫(∂t(u)⋅w + ((q-τ*(u⋅∇(q)))*(perp∘(F)))⋅w - (∇⋅(w))*(g*(D+h) + 0.5*(u⋅u)) + ∂t(D)*ϕ + w2⋅(F-u*D) + ϕ2*(q*D) +  (perp∘(∇(ϕ2)))⋅u - (ϕ2*fn) + ϕ*(∇⋅(F))- forcfunc(t)⋅w + ((cd*(norm(u))*u)/(D+1e-14))⋅w - fu(t)⋅w - fh(t)*ϕ)dΩ + ∫((g*(D+h) + 0.5*(u⋅u))*(w⋅nΓ) + nΓ⋅(perp∘(u))*ϕ2)dΓ #  - fq(t)*ϕ2 - fF(t)⋅w2
    jac(t,(u,D,q,F),(du,dD,dq,dF),(w,ϕ,ϕ2,w2)) = ∫(((dq - τ*(u⋅∇(dq) + du⋅∇(q)))*(perp∘(F)))⋅w + ((q-τ*(u⋅∇(q)))*(perp∘(dF)))⋅w - (∇⋅(w))*(g*(dD) + (du⋅u)) + w2⋅(dF - du*D - dD*u) + ϕ2*(q*dD) + ϕ2*(dq*D) + du⋅(perp∘(∇(ϕ2))) + ϕ*(∇⋅(dF)) + (cd*(dnorm(u,du)*u)/(D+h+1e-14) + cd*(norm(u)/(D+h+1e-14)*du) - cd*(norm(u)*u*dD)/((D+h+1e-14)*(D+h+1e-14)))⋅w )dΩ+ ∫((g*(dD) + (du⋅u))*(w⋅nΓ) + nΓ⋅(perp∘(du))*ϕ2)dΓ
    jac_t(t,(u,D),(dut,dDt),(w,ϕ)) = ∫(dut⋅w + dDt*ϕ)dΩ

    #Define operators and solvers
    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=true,method =:trust_region)
    ode_solver = ThetaMethod(nls,dt,0.5)
    xs = solve(ode_solver,op,uDn,T0,T)

    #Output results
    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"nonl_new"))do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"nonl_new.vtu"),cellfields=["u"=>un,"D"=>(Dn+h),"h"=>h])
            for (x,t) in xs
                u,D,F = x
                eh = (D - Dm(t))
                eu = (u - um(t))
                pvd[t] = createvtk(Ω,joinpath(dir,"nonl_new$t.vtu"),cellfields=["u"=>u,"D"=>(D+h),"h"=>h, "eu"=>eu, "eh"=>eh, "mu"=>um(t), "mζ"=>ζm(t), "fu"=>fu(t), "fh"=>fh(t)])
                println("done $t/$T")
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"nonl_new")) do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"nonlinear_topo0.0.vtu"),cellfields=["u"=>un,"D"=>(Dn+h),"h"=>h])
            for (x,t) in x
                u,D,F = x
                eh = (D - ζm(t))
                eu = (u - um(t))
                pvd[t] = createvtk(Ω,joinpath(dir,"nonl_new$t.vtu"),cellfields=["u"=>u,"D"=>(D+h),"h"=>h, "eu"=>eu, "eh"=>eh, "mu"=>um(t), "mh"=>ζm(t)])
                println("done $t/$T")
            end
        end
    end
end

order = 1
degree = 4
T = 10
dt = 0.1

APVM_run(order,degree,D₀,u₀,topography,forcfunc,dir,true,10,1,model,DC)