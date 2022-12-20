#SWQ_Semi_Linear_Without_F
using Pkg
Pkg.activate(".")

using Gridap
using SparseMatricesCSR
using SparseArrays
using WriteVTK
using LinearAlgebra
using LineSearches: BackTracking
using Gridap.TensorValues: meas

function Shallow_water_equations_newton_solver(
        order,degree,Tend,dt,B,L,x_points,y_points)
        dir = "./RESULTSSL"

    #Parameters
    latitude = 52
    η = 7.29e-5
    f = 2*η*sin(latitude*(π/180))
    g = 9.81
    H = 3 #Constant layer depth at rest
    U_start = 0.5
    σ = 2*pi/44700
    cD = 0.0025

    #Make model
    domain = (0,B,0,L)
    partition = (x_points,y_points)

    # Generate the model
    model = CartesianDiscreteModel(domain,partition;isperiodic=(false,true))

    #Make labels
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,2,5])
    add_tag_from_tags!(labels,"left",[7])
    add_tag_from_tags!(labels,"right",[8])
    add_tag_from_tags!(labels,"top",[3,4,6])
    add_tag_from_tags!(labels,"inside",[9])
    DC = ["left","right"]
    
    #Define triangulation
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    dω = Measure(Ω,degree,ReferenceDomain())
    Γ = BoundaryTriangulation(model,tags=DC)
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)

    v_boundary(x) = VectorValue(0.0,0.0)
    reffe_rt = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    V = TestFESpace(model,reffe_rt,dirichlet_tags=DC,dirichlet_masks=[(true,false),(true,false)])
    U = TransientTrialFESpace(V)#v_boundary) 
    
    reffe_lgn = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lgn)
    P = TransientTrialFESpace(Q) 


   # RTMM,RTMMchol = setup_and_factorize_mass_matrices(dΩ,V,Q,U,P)

    Y = MultiFieldFESpace([V,Q])#∇u, ∇h
    X = TransientMultiFieldFESpace([U,P])

    #Initial solutions
    #ζ₀(x,t) = 0.05*exp(-0.01*(x[1]-50)^2 - 0.01*(x[2]-25)^2)
    ζ₀(x,t) = 0.0
    ζ₀(t::Real) = x->ζ₀(x,t)
    u₀(x,t) = VectorValue(0.0,0.5)
    #u₀(x,t) = VectorValue(0.0,U_start)
    u₀(t::Real) = x->u₀(x,t)
    uhn = interpolate_everywhere([u₀(0.0),ζ₀(0.0)],X(0.0))

    #We use ordinary Hepkema structure
    #h₀(x,t) = 0.4*exp(-0.01*(x[2]-50)^2)
    #h₀(x,t) = 0.0
    
    #h₀(x,t) = 0.8*(exp(-0.001*(x[2]-L/4)^2)+exp(-0.001*(x[2]-L/3)^2)+exp(-0.001*(x[2]-L/2)^2)+exp(-0.001*(x[2]-L/3*2)^2)+exp(-0.001*(x[2]-L/4*3)^2))
    h₀(x,t) = exp(-0.000002*(x[1])^2 - 0.000002*(x[2]-L*0.8)^2) + exp(-0.000002*(x[1])^2 - 0.000002*(x[2]-L*0.4)^2) + exp(-0.000002*(x[1]-B)^2 - 0.000002*(x[2]-L*0.6)^2) + exp(-0.000002*(x[1]-B)^2 - 0.000002*(x[2]-L*0.2)^2)
    h₀(t::Real) = x->h₀(x,t)
    h = interpolate_everywhere(h₀(0.0),P(0.0))
    coriolis(u) = f*VectorValue(u[2],-u[1])
    forcfunc(t) = VectorValue(-f*U_start*cos(σ*t),-σ*U_start*sin(σ*t)+cD*abs(U_start*cos(σ*t))*U_start*cos(σ*t)/H)

    #Stabilization Parameters
    α = 1e-6
    ν = 1.0e-6 #3.225


#Added a Stabilization term inside the residual following E. Hanert et al. / Ocean Modelling 5 (2002) 17–35
    res(t,(u,ζ),(w,ϕ)) = ∫(∂t(u)⋅w - g*(∇⋅(w))*ζ+ ∇(u)'⋅u⋅w - forcfunc(t)⋅w + (coriolis∘u)⋅w + cD* (meas∘u) * u⋅w/(ζ+H-h)+∂t(ζ)*ϕ-(ζ+H-h)*u ⋅(∇(ϕ)) + α*(∇(ζ)⋅∇(ϕ)) + ν * (∇⋅u)*(∇⋅w)  )dΩ + ∫(g*(h)*(w⋅nΓ))dΓ 
    jac(t,(u,ζ),(du,dζ),(w,ϕ)) = ∫(-g*(∇⋅(w))*dζ + ∇(du)'⋅u⋅w + ∇(u)'⋅du⋅w + (coriolis∘du)⋅w +cD*(meas∘u)*u*dζ⋅w/((ζ+H-h)*(ζ+H-h))+cD*u⋅du*u⋅w/((meas∘(u+1e-14))*(ζ+H-h)) + cD*(meas∘u)*du⋅w/(ζ+H-h) - (dζ*u+ζ*u) ⋅(∇(ϕ)) +α*(∇(dζ)⋅∇(ϕ)) + ν * (∇⋅du)*(∇⋅w) )dΩ + ∫(g*dζ*(w⋅nΓ))dΓ
    jac_t(t,(u,ζ),(dut,dζt),(w,ϕ)) = ∫(dut⋅w + dζt*ϕ)dΩ
    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=true,linesearch=BackTracking())
    #Tend = 20*5
    ode_solver = ThetaMethod(nls,dt,0.5)
    x = solve(ode_solver,op,uhn,0.0,Tend)
    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"1d-topo-output"))do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"1d-topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn+b),"b"=>b])
            for (x,t) in x
                u,ζ = x
                pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>ζ,"h"=>h])
                println("done $t/$Tend")
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"1d-topo-output")) do pvd
            #pvd[0.0] = createvtk(Ω,joinpath(dir,"1d-topo0.0.vtu"),cellfields=["u"=>un,"h"=>(hn+b),"b"=>b])
            for (x,t) in x
                u,ζ = x
                pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>ζ,"h"=>h])
                println("done $t/$Tend")
            end
        end
    end
end

order = 1
degree = 3
Tend = 500
dt = 4
B = 1000
L = 10000
x_points = 20
y_points = 100

#Initialization to speed up second time function is called
Shallow_water_equations_newton_solver(order,degree,2*dt,dt,B,L,x_points,y_points)

#Real function
Shallow_water_equations_newton_solver(order,degree,Tend,dt,B,L,x_points,y_points)