
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
f = 0# 2*η*sin(latitude*(π/180))
σ = 2*pi/Tend
H = 3.0


W = 1
L = 1
#partition = (50,50)
domain = (0,L,0,W)
model = GmshDiscreteModel("swe-solver/meshes/1x1periodic01.msh")
# model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false)) 
# labels = get_face_labeling(model)
# add_tag_from_tags!(labels,"bottom",[1,2,5])
# add_tag_from_tags!(labels,"left",[7])
# add_tag_from_tags!(labels,"right",[8])
# add_tag_from_tags!(labels,"top",[3,4,6])
# add_tag_from_tags!(labels,"inside",[9])
DC = ["bottom","top"]
Tend = 10
dt = 1

A = (π/W)
B = (π/L)
C = ((1 * π)/(2 * Tend))

D = (π/W)
E = (π/(1 * L))
F = A
G = B

top = 0.0
xref = W / 2
cd = cD
g = 9.81

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
    σ = 2*pi/Tend
    H = 3
    f = VectorValue(0.0,0.0)#VectorValue(σ*U_start*cos(σ*t)+cD/H*abs(U_start*sin(σ*t))*U_start*sin(σ*t),f*U_start*sin(σ*t))
    f
end

#Manufactured functions


function fu_((x,y),t)
    # fu = VectorValue((cd*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)*sqrt((0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)^2+(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)^2))/(0.05*sin(F*x)*sin(G*y)*sin(C*t)-top+2*H)- 
    # 0.1*E*sin(E*(x-xref))*sin(D*y)*sin(C*t)*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)+ 
    # 0.1*D*cos(E*(x-xref))*cos(D*y)*sin(C*t)*(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)- 
    # f*(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)+ 
    # 0.05*F*g*cos(F*x)*sin(G*y)*sin(C*t) + 
    # 0.1*C*cos(E*(x-xref))*sin(D*y)*cos(C*t),
    # (cd*(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)*sqrt((0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)^2+(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)^2))/(0.05*sin(F*x)*sin(G*y)*sin(C*t)-top+2*H)+ 
    # 0.1*A*cos(A*x)*sin(B*y)*sin(C*t)*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)+ 
    # f*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)+ 
    # 0.1*B*sin(A*x)*cos(B*y)*sin(C*t)*(0.1*sin(A*x)*sin(B*y)*sin(C*z)+1.0)+ 
    # 0.05*G*g*sin(F*x)*cos(G*y)*sin(C*t)+0.1*C*sin(A*x)*sin(B*y)*cos(C*t))
    #fu = VectorValue((cd*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)*sqrt((0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)^2+(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)^2))/(0.05*sin(F*x)*sin(G*y)*sin(C*t)-top+2*H)-0.1*E*sin(E*(x-xref))*sin(D*y)*sin(C*t)*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)+0.1*D*cos(E*(x-xref))*cos(D*y)*sin(C*t)*(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)-f*(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)+0.05*F*g*cos(F*x)*sin(G*y)*sin(C*t)+0.1*C*cos(E*(x-xref))*sin(D*y)*cos(C*t),(cd*(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)*sqrt((0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)^2+(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)^2))/(0.05*sin(F*x)*sin(G*y)*sin(C*t)-top+2*H)+0.1*A*cos(A*x)*sin(B*y)*sin(C*t)*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)+f*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)+0.1*B*sin(A*x)*cos(B*y)*sin(C*t)*(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)+0.05*G*g*sin(F*x)*cos(G*y)*sin(C*t)+0.1*C*sin(A*x)*sin(B*y)*cos(C*t))
    fu = VectorValue((cd*(0.1*sin(A*y)*sin(C*t)+1.0)*abs(0.1*sin(A*y)*sin(C*t)+1.0))/(0.05*x*sin(C*t)-top+2*H)+0.05*g*sin(C*t)+0.1*C*sin(A*y)*cos(C*t), f*(0.1*sin(A*y)*sin(C*t)+1.0))
    fu
end

# function fu_((x,y),t)
#     fu = VectorValue(sin(x)*sin(t) + y*t, x*t + y*t)
#     fu
# end

function fh_((x,y),t)
    # fh = -0.1*E*sin(E*(x-xref))*sin(D*y)*sin(C*t)*(0.05*sin(F*x)*sin(G*y)*sin(C*t)-top+2*H)+
    # 0.1*B*sin(A*x)*cos(B*y)*sin(C*t)*(0.05*sin(F*x)*sin(G*y)*sin(C*t)-top+2*H)+
    # 0.05*F*cos(F*x)*sin(G*y)*sin(C*t)*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)+
    # 0.05*G*sin(F*x)*cos(G*y)*sin(C*t)*(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)+
    # 0.05*C*sin(F*x)*sin(G*y)*cos(C*t)
    fh = 0.05*sin(C*t)*(0.1*sin(A*y)*sin(C*t)+1.0)+0.05*C*x*cos(C*t)
    fh
end

function mu_((x,y),t)
    # u = VectorValue(0.1*sin(D*y)*cos(E*(x-xref))*sin(C*t)+1.0,
    # 0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)
    u = VectorValue(1.0 + 0.1 * sin(A*y) * sin(C * t), 0.0)
    u    
end

function mh_((x,y),t)
    #h = 0.05*sin(F*x)*sin(G*y)*sin(C*t)+H
    h = H + 0.05*x*sin(C * t)
    h
end


outputdir = "output_swe"
dir = "NL_SWE_APVM_test_gc"
if !isdir(outputdir)
    mkdir(outputdir)
end

if !isdir(joinpath(outputdir,dir))
    mkdir(joinpath(outputdir,dir))
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

    #Forcing function on u(t)
    forcfunc(t) = x -> forcefunc(x,t)

    mu(t) = x -> mu_(x,t)

    mh(t) = x -> mh_(x,t)

    fu(t) = x -> fu_(x,t)

    fh(t) = x -> fh_(x,t)

    mu_b(x,t::Real) = 
    # VectorValue(0.1*sin(D*x[2])*cos(E*(x[1]-xref))*sin(C*t)+1.0,#0.1*sin(D*x[2])*cos(E*(x[1]-xref))*
    # 0.1*sin(A*x[1])*sin(B*x[2])*sin(C*t)+1.0)#0.1*sin(A*x[1])*sin(B*x[2])*
    VectorValue(0.1 * sin(A * x[2]) * sin(C * t) + 1.0, 0.0)
    mu_b(t::Real) = x -> mu_b(x,t)

    #Make FE spaces
    if periodic
        reffe_rt = ReferenceFE(raviart_thomas, Float64,order)#
        V = TestFESpace(model,reffe_rt;conformity=:HDiv,dirichlet_tags=DC)#
        #U = TransientTrialFESpace(V,[udc,udc])
        U = TransientTrialFESpace(V,[mu_b,mu_b])
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

    #  #Manufactured functions
    #  W = 1
    #  L = 10
 
    #  η = 7.29e-5
    #  latitude = 50
    #  cD = 0.0025
    #  f = 2*η*sin(latitude*(π/180))
    #  σ = 2*pi/Tend
    #  H = 3.0
 
    #  A = (π/W)
    #  B = (π/L)
    #  C = ((4 * π)/(1 * Tend))
 
    #  D = (π/W)
    #  E = (π/(7 * L))
    #  F = A
    #  G = B
 
    #  top = 0.0
    #  xref = W / 2
    #  cd = cD
 
    #  fu(x,t::Real) = VectorValue((cd*(0.1*cos(E*(x[1]-xref))*sin(D*x[2])*sin(C*t)+1.0)*sqrt((0.1*cos(E*(x[1]-xref))*sin(D*x[2])*sin(C*t)+1.0)^2+(0.1*sin(A*x[1])*sin(B*x[2])*sin(C*t)+1.0)^2))/(0.05*sin(F*x[1])*sin(G*x[2])*sin(C*t)-top+2*H)- 
    #  0.1*E*sin(E*(x[1]-xref))*sin(D*x[2])*sin(C*t)*(0.1*cos(E*(x[1]-xref))*sin(D*x[2])*sin(C*t)+1.0)+ 
    #  0.1*D*cos(E*(x[1]-xref))*cos(D*x[2])*sin(C*t)*(0.1*sin(A*x[1])*sin(B*x[2])*sin(C*t)+1.0)- 
    #  f*(0.1*sin(A*x[1])*sin(B*x[2])*sin(C*t)+1.0)+ 
    #  0.05*F*g*cos(F*x[1])*sin(G*x[2])*sin(C*t) + 
    #  0.1*C*cos(E*(x[1]-xref))*sin(D*x[2])*cos(C*t),
    #  (cd*(0.1*sin(A*x[1])*sin(B*x[2])*sin(C*t)+1.0)*sqrt((0.1*cos(E*(x[1]-xref))*sin(D*x[2])*sin(C*t)+1.0)^2+(0.1*sin(A*x[1])*sin(B*x[2])*sin(C*t)+1.0)^2))/(0.05*sin(F*x[1])*sin(G*x[2])*sin(C*t)-top+2*H)+ 
    #  0.1*A*cos(A*x[1])*sin(B*x[2])*sin(C*t)*(0.1*cos(E*(x[1]-xref))*sin(D*x[2])*sin(C*t)+1.0)+ 
    #  f*(0.1*cos(E*(x[1]-xref))*sin(D*x[2])*sin(C*t)+1.0)+ 
    #  0.1*B*sin(A*x[1])*cos(B*x[2])*sin(C*t)*(0.1*sin(A*x[1])*sin(B*x[2])*sin(C*z)+1.0)+ 
    #  0.05*G*g*sin(F*x[1])*cos(G*x[2])*sin(C*t)+0.1*C*sin(A*x[1])*sin(B*x[2])*cos(C*t))
    #  fu(t::Real) = x -> fu(x,t)
 
    #  fh(x,t::Real) = -0.1*E*sin(E*(x[1]-xref))*sin(D*x[2])*sin(C*t)*(0.05*sin(F*x[1])*sin(G*x[2])*sin(C*t)-top+2*H)+
    #  0.1*B*sin(A*x[1])*cos(B*x[2])*sin(C*t)*(0.05*sin(F*x[1])*sin(G*x[2])*sin(C*t)-top+2*H)+
    #  0.05*F*cos(F*x[1])*sin(G*x[2])*sin(C*t)*(0.1*cos(E*(x[1]-xref))*sin(D*x[2])*sin(C*t)+1.0)+
    #  0.05*G*sin(F*x[1])*cos(G*x[2])*sin(C*t)*(0.1*sin(A*x[1])*sin(B*x[2])*sin(C*t)+1.0)+
    #  0.05*C*sin(F*x[1])*sin(G*x[2])*cos(C*t)
    #  fh(t::Real) = x -> fh(x,t)
 
    #  mu(x,t::Real) = VectorValue(sin(C*t)+1.0,#0.1*sin(D*x[2])*cos(E*(x[1]-xref))*
    #  sin(C*t)+1.0)#0.1*sin(A*x[1])*sin(B*x[2])*
    #  mu(t::Real) = x -> mu(x,t)
 
    #  mh(x,t::Real) = 0.05*sin(F*x[1])*sin(G*x[2])*sin(C*t)+H
    #  mh(t::Real) = x -> mh(x,t)


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

    #Norm operator
    norm(u) = meas∘(u) + 1e-14
    dnorm(u,du) = u ⋅ du / norm(u)


    #Define residual, jacobian in space and time
    res(t,(u,D,q,F),(w,ϕ,ϕ2,w2)) = ∫(∂t(u)⋅w + ((q-τ*(u⋅∇(q)))*(perp∘(F)))⋅w - (∇⋅(w))*(g*(D+h) + 0.5*(u⋅u)) + ∂t(D)*ϕ + w2⋅(F-u*D) + ϕ2*(q*D) +  (perp∘(∇(ϕ2)))⋅u - (ϕ2*fn) + ϕ*(∇⋅(F))- forcfunc(t)⋅w + ((cd*(norm(u))*u)/(D+1e-14))⋅w - fu(t)⋅w - fh(t)*ϕ)dΩ + ∫((g*(D+h) + 0.5*(u⋅u))*(w⋅nΓ) + nΓ⋅(perp∘(u))*ϕ2)dΓ 
    jac(t,(u,D,q,F),(du,dD,dq,dF),(w,ϕ,ϕ2,w2)) = ∫(((dq - τ*(u⋅∇(dq) + du⋅∇(q)))*(perp∘(F)))⋅w + ((q-τ*(u⋅∇(q)))*(perp∘(dF)))⋅w - (∇⋅(w))*(g*(dD) + (du⋅u)) + w2⋅(dF - du*D - dD*u) + ϕ2*(q*dD) + ϕ2*(dq*D) + du⋅(perp∘(∇(ϕ2))) + ϕ*(∇⋅(dF)) + (cd*(dnorm(u,du)*u)/(D+h+1e-14) + cd*(norm(u)/(D+h+1e-14)*du) - cd*(norm(u)*u*dD)/((D+h+1e-14)*(D+h+1e-14)))⋅w )dΩ+ ∫((g*(dD) + (du⋅u))*(w⋅nΓ) + nΓ⋅(perp∘(du))*ϕ2)dΓ
    jac_t(t,(u,D),(dut,dDt),(w,ϕ)) = ∫(dut⋅w + dDt*ϕ)dΩ

    #Define operators and solvers
    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=true,method =:trust_region)
    ode_solver = ThetaMethod(nls,dt,0.5)
    xs = solve(ode_solver,op,uDn,T0,Tend)

    #Output results
    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"nonlinear_topo"))do pvd
            pvd[0.0] = createvtk(Ω,joinpath(dir,"nonlinear_topo0.0.vtu"),cellfields=["u"=>un,"D"=>(Dn+h),"h"=>h])
            for (x,t) in xs
                u,D,F = x
                eh = (D - mh(t))
                eu = (u - mu(t))
                pvd[t] = createvtk(Ω,joinpath(dir,"nonlinear_topo$t.vtu"),cellfields=["u"=>u,"D"=>(D+h),"h"=>h, "eu"=>eu, "eh"=>eh, "mu"=>mu(t), "mh"=>mh(t)])
                println("done $t/$Tend")
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"nonlinear_topo")) do pvd
            pvd[0.0] = createvtk(Ω,joinpath(dir,"nonlinear_topo0.0.vtu"),cellfields=["u"=>un,"D"=>(Dn+h),"h"=>h])
            for (x,t) in x
                u,D,F = x
                pvd[t] = createvtk(Ω,joinpath(dir,"nonlinear_topo$t.vtu"),cellfields=["u"=>u,"D"=>(D+h),"h"=>h, "mu"=>mu(t), "mh"=>mh(t)])
                println("done $t/$Tend")
            end
        end
    end
end

order = 1
degree = 4
APVM_run(order,degree,D₀,u₀,topography,forcfunc,dir,true,Tend,dt,model,DC)