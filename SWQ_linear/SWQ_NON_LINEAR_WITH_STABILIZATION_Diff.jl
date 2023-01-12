#= SWQ_Non_Linear_Stabilization with simple diffusion added
In this code the original structure of Hepkema is used regarding equations.
The shallow water equations are solved using depth-averaged equations. The diffusion-terms Δu and Δζ with respectively ν and α as Stabilization Parameters
x and y are flipped compared to Hepkema
=#

#Add Pkg
using Pkg
Pkg.activate(".")

#Add packages
using Gridap
using SparseMatricesCSR
using SparseArrays
using WriteVTK
using LinearAlgebra
using LineSearches: BackTracking
using Gridap.TensorValues: meas
using CSV, DelimitedFiles


#function to solve the equations
function Shallow_water_equations_newton_solver(
        order,degree,Tend,dt,B,L,x_points,y_points,dir,show_result,latitude,ζ₀,u₀,h₀)

    #''''''''''''''Parameters from Hepkema''''''''''''''##
    η = 7.29e-5                             #angular speed of Earth rotation        (s^(-1))
    f = 2*η*sin(latitude*(π/180))           #coriolis parameter                     (s^(-1))
    g = 9.81                                #Gravitational constant                 (ms^(-2))
    H = 3                                   #Constant layer depth at rest           (m)
    U_start = 0.5                           #Background current amplitude           (ms^(-1))
    σ = 2*pi/44700                          #Tidal frequency                        (s^(-1))
    cD = 0.0025                             #Drag coefficient                       ()

    ##''''''''''''''Stabilization Parameters''''''''''''''##
    α = 1e-6                                #Based on ν
    ν = 1e-6                                #From Anna Louka

    ##''''''''''''''Functions''''''''''''''##
    coriolis(u) = f*VectorValue(u[2],-u[1])                                                                             #coriolis
    forcfunc(t) = VectorValue(-f*U_start*cos(σ*t),-σ*U_start*sin(σ*t)+cD*abs(U_start*cos(σ*t))*U_start*cos(σ*t))      #Fₚ from Hepkema

    #''''''''''''''Make model''''''''''''''##
    domain = (0,B,0,L)                      #Original x and y flipped
    partition = (x_points,y_points)         #Partition

    ##''''''''''''''Generate the model''''''''''''''##
    model = CartesianDiscreteModel(domain,partition;isperiodic=(false,true))  #Periodic on y-axis

    ##''''''''''''''Make labels''''''''''''''##
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"bottom",[1,2,5])
    add_tag_from_tags!(labels,"left",[7])
    add_tag_from_tags!(labels,"right",[8])
    add_tag_from_tags!(labels,"top",[3,4,6])
    add_tag_from_tags!(labels,"inside",[9])
    DC = ["left","right"]
    
    ##''''''''''''''Define triangulation''''''''''''''##
    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)
    Γ = BoundaryTriangulation(model,tags=DC)
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)

    ##''''''''''''''Define function spaces''''''''''''''##
    v_boundary(x) = VectorValue(0.0,0.0)
    reffe_rt = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    V = TestFESpace(model,reffe_rt,dirichlet_tags=DC,dirichlet_masks=[(true,false),(true,false)])  #Zero at the x-boundaries
    U = TransientTrialFESpace(V)#v_boundary) 
    
    reffe_lgn = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lgn)
    P = TransientTrialFESpace(Q) 

    Y = MultiFieldFESpace([V,Q])  #u, ζ
    X = TransientMultiFieldFESpace([U,P])

    ##''''''''''''''Initial solution of ζ''''''''''''''##
    #ζ₀(x,t) = 0.05*exp(-0.01*(x[1]-50)^2 - 0.01*(x[2]-25)^2)           #Druppel
    

    ##''''''''''''''Initial solution of u''''''''''''''##
    
    #Interpolate initial solution
    uhn = interpolate_everywhere([u₀(0.0),ζ₀(0.0)],X(0.0))


    ##''''''''''''''Define bottom topography''''''''''''''##
    #h₀(x,t) = 0.4*exp(-0.01*(x[2]-50)^2)                               #Druppel
    #h₀(x,t) = 0.0                                                      #No bottom topography
    #Different vertical bottoms:
    #h₀(x,t) = 0.8*(exp(-0.001*(x[2]-L/4)^2)+exp(-0.001*(x[2]-L/3)^2)+exp(-0.001*(x[2]-L/2)^2)+exp(-0.001*(x[2]-L/3*2)^2)+exp(-0.001*(x[2]-L/4*3)^2))
    #River bottom
    #h₀(x,t) = exp(-0.000002*(x[1])^2 - 0.000002*(x[2]-L*0.8)^2) + exp(-0.000002*(x[1])^2 - 0.000002*(x[2]-L*0.4)^2) + exp(-0.000002*(x[1]-B)^2 - 0.000002*(x[2]-L*0.6)^2) + exp(-0.000002*(x[1]-B)^2 - 0.000002*(x[2]-L*0.2)^2)
    h = interpolate_everywhere(h₀(0.0),P(0.0))
    

    ##''''''''''Different functions with their derivitave''''''''##
    #Velocity change 
    func_ut(t,(u,ζ),(w,ϕ)) =  ∂t(u)⋅w

    #Convection
    func_con(t,(u,ζ),(w,ϕ)) = ∇(u)'⋅u⋅w
    dfunc_con(t,(u,ζ),(du,dζ),(w,ϕ)) = ∇(du)'⋅u⋅w + ∇(u)'⋅du⋅w

    #Coriolis
    func_cor(t,(u,ζ),(w,ϕ)) = (coriolis∘u)⋅w
    dfunc_cor(t,(u,ζ),(du,dζ),(w,ϕ)) = (coriolis∘du)⋅w

    #Drag coefficient term
    func_cD(t,(u,ζ),(w,ϕ)) = cD* (meas∘u) * u⋅w/(ζ+H-h)
    dfunc_cD(t,(u,ζ),(du,dζ),(w,ϕ)) = cD/(ζ+H-h)*w⋅(-(meas∘u)*u*dζ/(ζ+H-h)+u⋅du*u/((meas∘(u+1e-14))) + (meas∘u)*du)

    #Gravitational (without boundary)
    func_g(t,(u,ζ),(w,ϕ)) = - g*(∇⋅(w))*ζ
    dfunc_g(t,(u,ζ),(du,dζ),(w,ϕ)) = -g*(∇⋅(w))*dζ

    #Forcing function
    func_Fₚ(t,(u,ζ),(w,ϕ)) = - forcfunc(t)⋅w

    #ζ change
    func_ζt(t,(u,ζ),(w,ϕ)) = ∂t(ζ)*ϕ

    #ζ+Velocity function
    func_h(t,(u,ζ),(w,ϕ)) = -(ζ+H-h)*u ⋅(∇(ϕ))
    dfunc_h(t,(u,ζ),(du,dζ),(w,ϕ)) = - (dζ*u+(ζ+H-h)*u) ⋅(∇(ϕ))

    #Boundary
    func_boun(t,(u,ζ),(w,ϕ)) = g*(ζ)*(w⋅nΓ)
    dfunc_boun(t,(u,ζ),(du,dζ),(w,ϕ)) = g*(dζ)*(w⋅nΓ)

    #Stabilization function ζ
    func_stabζ(t,(u,ζ),(w,ϕ)) = α*(∇(ζ)⋅∇(ϕ))                                               
    dfunc_stabζ(t,(u,ζ),(du,dζ),(w,ϕ)) = α*(∇(dζ)⋅∇(ϕ))                                   


    #Stabilization function u
    func_stabu(t,(u,ζ),(w,ϕ)) = ν * (∇⋅u)*(∇⋅w)                          
    dfunc_stabu(t,(u,ζ),(du,dζ),(w,ϕ)) = ν * (∇⋅du)*(∇⋅w)



    ##''''''''''''''Residual''''''''''''''##
    #Added a Stabilization term inside the residual following E. Hanert et al. / Ocean Modelling 5 (2002) 17–35
    res(t,(u,ζ),(w,ϕ)) = ∫(func_ut(t,(u,ζ),(w,ϕ)) +                     #Velocity change
        func_con(t,(u,ζ),(w,ϕ)) +                                       #Convection
        func_cor(t,(u,ζ),(w,ϕ)) +                                       #Coriolis
        func_cD(t,(u,ζ),(w,ϕ))  +                                       #Drag coefficient term
        func_g(t,(u,ζ),(w,ϕ))   +                                       #Gravitational
        func_Fₚ(t,(u,ζ),(w,ϕ))  +                                       #Forcing function
        func_ζt(t,(u,ζ),(w,ϕ))  +                                       #ζ change
        func_h(t,(u,ζ),(w,ϕ))   +                                       #ζ+Velocity function
        func_stabζ(t,(u,ζ),(w,ϕ)) +                                     #Stabilization ζ
        func_stabu(t,(u,ζ),(w,ϕ)))dΩ +                                  #Stabilization u
        ∫(func_boun(t,(u,ζ),(w,ϕ)))dΓ                                   #Boundary

    ##''''''''''''''Jacobian''''''''''''''##
    jac(t,(u,ζ),(du,dζ),(w,ϕ)) = ∫(dfunc_con(t,(u,ζ),(du,dζ),(w,ϕ)) +   #Convection
        dfunc_cor(t,(u,ζ),(du,dζ),(w,ϕ)) +                              #Coriolis
        dfunc_cD(t,(u,ζ),(du,dζ),(w,ϕ))  +                              #Drag coefficient term
        dfunc_g(t,(u,ζ),(du,dζ),(w,ϕ))   +                              #Gravitational
        dfunc_h(t,(u,ζ),(du,dζ),(w,ϕ))   +                              #ζ+Velocity function
        dfunc_stabζ(t,(u,ζ),(du,dζ),(w,ϕ)) +                            #Stabilization ζ
        dfunc_stabu(t,(u,ζ),(du,dζ),(w,ϕ)))dΩ +                         #Stabilization u
        ∫(dfunc_boun(t,(u,ζ),(du,dζ),(w,ϕ)))dΓ                          #Boundary

    ##''''''''''''''Jacobian for t''''''''''''''##
    jac_t(t,(u,ζ),(dut,dζt),(w,ϕ)) = ∫(dut⋅w + dζt*ϕ)dΩ

    ##''''''''''''''Solver with ThetaMethod''''''''''''''##
    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=show_result,linesearch=BackTracking())
    ode_solver = ThetaMethod(nls,dt,0.5)
    x = solve(ode_solver,op,uhn,0.0,Tend)  

    
    #''''''''''''''Saving''''''''''''''##
    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"1d-topo-output"))do pvd
            for (x,t) in x
                if t%(dt*10) ==0
                    u,ζ = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>ζ,"h"=>h])
                    println("Saved")
                end
                println("done $t/$Tend")
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"1d-topo-output")) do pvd
            for (x,t) in x
                if t%(dt*10) ==0
                    u,ζ = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>ζ,"h"=>h])
                    println("Saved")
                end
                println("done $t/$Tend")
            end
        end
    end
end

    #=
##''''''''''''''Shallow_water_equations_newton_solver input''''''''''''''##
    order                   Order of the spaces
    degree                  Degree for the boundaries
    Tend                    The total time
    dt                      Time step 
    B                       Channel width
    L                       Channel length
    x_points                Amount of points in the x-direction
    y_points                Amount of points in the y-direction
    dir                     Name of the direction the files get savid in
    show_result             True or false (print show_trace)
    latitude                Latitude   
    ζ₀                      Initial solution for ζ as a function of t
    u₀                      Initial solution for u₀ as a function of t
    h₀                      Given topography
    =#

##''''''''''''''Variable parameters for initial setup''''''''''''''##
order = 1                                               #Order of the spaces
degree = 3                                              #Degree for the boundaries
Tend = 44700                                            #The total time                              (s)
dt = 5                                                  #Time step                                   (s)
B = 1000                                                #Channel width                               (m)
L = 10000                                               #Channel length                              (m)
x_points = 20                                           #Amount of points in the x-direction -> dx = 50 m 
y_points = 100                                          #Amount of points in the y-direction -> dy = 100 m 
dir = "./setup_direction"                               #Name of the direction the files get savid in
show_results = false                                    #Show the intermediate steps
latitude = 52                                           #latitude                               (°)
ζ₀(x,t) = 0.0                                           #Hepkema original initial solution for ζ
ζ₀(t::Real) = x->ζ₀(x,t)
u₀(x,t) = VectorValue(0.0,0.0)                          #Hepkema original initial solution for u
u₀(t::Real) = x->u₀(x,t)
h₀(x,t) = 0.1 * cos(π/B * x[1]) * cos(2π/L * x[2])      #Hepkema original function for h
h₀(t::Real) = x->h₀(x,t)


##''''''''''''''Initialization to speed up second time function is called''''''''''''''##
Shallow_water_equations_newton_solver(order,degree,2*dt,dt,B,L,x_points,y_points,dir,show_results,latitude,ζ₀,u₀,h₀)



##''''''''''''''Variable parameters for the desired problem''''''''''''''##
order = 1                                               #Order of the spaces
degree = 3                                              #Degree for the boundaries
Tend = 44700                                            #The total time                              (s)
dt = 5                                                  #Time step                                   (s)
B = 1000                                                #Channel width                               (m)
L = 10000                                               #Channel length                              (m)
x_points = 20                                           #Amount of points in the x-direction -> dx = 50 m 
y_points = 100                                          #Amount of points in the y-direction -> dy = 100 m
dir = "./RESULTSDiffusion10"                            #Name of the direction the files get savid in
show_results = false                                    #Show the intermediate results
latitude = 52                                           #latitude                               (°)
ζ₀(x,t) = 0.0                                           #Hepkema original initial solution for ζ
ζ₀(t::Real) = x->ζ₀(x,t)
u₀(x,t) = VectorValue(0.0,0.0)                          #Hepkema original initial solution for u
u₀(t::Real) = x->u₀(x,t)
h₀(x,t) = 0.1 * cos(π/B * x[1]) * cos(2π/L * x[2])      #Hepkema original function for h
h₀(t::Real) = x->h₀(x,t)


##''''''''''''''Real function''''''''''''''##
Shallow_water_equations_newton_solver(order,degree,Tend,dt,B,L,x_points,y_points,dir,show_results,latitude,ζ₀,u₀,h₀)




