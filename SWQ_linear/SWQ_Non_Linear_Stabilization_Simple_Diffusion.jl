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
using CSV, DataFrames


#function to solve the equations
function Shallow_water_equations_newton_solver(
        order,degree,Tend,dt,B,L,x_points,y_points,dir)

    #Parameters from Hepkema
    latitude = 52                           #latitude                               (°)
    η = 7.29e-5                             #angular speed of Earth rotation        (s^(-1))
    f = 2*η*sin(latitude*(π/180))           #coriolis parameter                     (s^(-1))
    g = 9.81                                #Gravitational constant                 (ms^(-2))
    H = 3                                   #Constant layer depth at rest           (m)
    U_start = 0.5                           #Background current amplitude           (ms^(-1))
    σ = 2*pi/44700                          #Tidal frequency                        (s^(-1))
    cD = 0.0025                             #Drag coefficient                       ()

    #Stabilization Parameters
    α = 1e-6                                #Based on ν
    ν = 1e-6                                #From Anna Louka

    #Functions 
    coriolis(u) = f*VectorValue(u[2],-u[1])                                                                             #coriolis
    forcfunc(t) = VectorValue(-f*U_start*cos(σ*t),-σ*U_start*sin(σ*t)+cD*abs(U_start*cos(σ*t))*U_start*cos(σ*t)/H)      #Fₚ from Hepkema

    #Make model
    domain = (0,B,0,L)                      #Original x and y flipped
    partition = (x_points,y_points)         #Partition

    # Generate the model
    model = CartesianDiscreteModel(domain,partition;isperiodic=(false,true))  #Periodic on y-axis

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
    Γ = BoundaryTriangulation(model,tags=DC)
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,degree)

    #Define function spaces
    v_boundary(x) = VectorValue(0.0,0.0)
    reffe_rt = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    V = TestFESpace(model,reffe_rt,dirichlet_tags=DC,dirichlet_masks=[(true,false),(true,false)])  #Zero at the x-boundaries
    U = TransientTrialFESpace(V)#v_boundary) 
    
    reffe_lgn = ReferenceFE(lagrangian,Float64,order)
    Q = TestFESpace(model,reffe_lgn)
    P = TransientTrialFESpace(Q) 

    Y = MultiFieldFESpace([V,Q])  #u, ζ
    X = TransientMultiFieldFESpace([U,P])

    #Initial solution of ζ
    #ζ₀(x,t) = 0.05*exp(-0.01*(x[1]-50)^2 - 0.01*(x[2]-25)^2)           #Druppel
    ζ₀(x,t) = 0.0                                                       #Hepkema original initial solution
    ζ₀(t::Real) = x->ζ₀(x,t)

    #Initial solution of u
    #u₀(x,t) = VectorValue(0.0,0.0)                                     #Hepkema original initial solution
    u₀(x,t) = VectorValue(0.0,U_start)                                  #Maybe better initial solution
    u₀(t::Real) = x->u₀(x,t)

    #Interpolate initial solution
    uhn = interpolate_everywhere([u₀(0.0),ζ₀(0.0)],X(0.0))


    #Define bottom topography

    #h₀(x,t) = 0.4*exp(-0.01*(x[2]-50)^2)                               #Druppel
    #h₀(x,t) = 0.0                                                      #No bottom topography
    #Different vertical bottoms:
    #h₀(x,t) = 0.8*(exp(-0.001*(x[2]-L/4)^2)+exp(-0.001*(x[2]-L/3)^2)+exp(-0.001*(x[2]-L/2)^2)+exp(-0.001*(x[2]-L/3*2)^2)+exp(-0.001*(x[2]-L/4*3)^2))
    #River bottom
    h₀(x,t) = exp(-0.000002*(x[1])^2 - 0.000002*(x[2]-L*0.8)^2) + exp(-0.000002*(x[1])^2 - 0.000002*(x[2]-L*0.4)^2) + exp(-0.000002*(x[1]-B)^2 - 0.000002*(x[2]-L*0.6)^2) + exp(-0.000002*(x[1]-B)^2 - 0.000002*(x[2]-L*0.2)^2)
    h₀(t::Real) = x->h₀(x,t)
    h = interpolate_everywhere(h₀(0.0),P(0.0))
    




    #Residual 
    #Added a Stabilization term inside the residual following E. Hanert et al. / Ocean Modelling 5 (2002) 17–35
    res(t,(u,ζ),(w,ϕ)) = ∫(∂t(u)⋅w - g*(∇⋅(w))*ζ+ ∇(u)'⋅u⋅w - forcfunc(t)⋅w + (coriolis∘u)⋅w + cD* (meas∘u) * u⋅w/(ζ+H-h)+∂t(ζ)*ϕ-(ζ+H-h)*u ⋅(∇(ϕ)) + α*(∇(ζ)⋅∇(ϕ)) + ν * (∇⋅u)*(∇⋅w)  )dΩ + ∫(g*(ζ)*(w⋅nΓ))dΓ 

    #Jacobian
    jac(t,(u,ζ),(du,dζ),(w,ϕ)) = ∫(-g*(∇⋅(w))*dζ + ∇(du)'⋅u⋅w + ∇(u)'⋅du⋅w + (coriolis∘du)⋅w +cD*(meas∘u)*u*dζ⋅w/((ζ+H-h)*(ζ+H-h))+cD*u⋅du*u⋅w/((meas∘(u+1e-14))*(ζ+H-h)) + cD*(meas∘u)*du⋅w/(ζ+H-h) - (dζ*u+(ζ+H-h)*u) ⋅(∇(ϕ)) +α*(∇(dζ)⋅∇(ϕ)) + ν * (∇⋅du)*(∇⋅w) )dΩ + ∫(g*dζ*(w⋅nΓ))dΓ

    #Jacobian for t
    jac_t(t,(u,ζ),(dut,dζt),(w,ϕ)) = ∫(dut⋅w + dζt*ϕ)dΩ



    #Solver with ThetaMethod
    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=true,linesearch=BackTracking())
    ode_solver = ThetaMethod(nls,dt,0.5)
    x = solve(ode_solver,op,uhn,0.0,Tend)  
    
    x_hep_test = []
    y_hep_test = []
    i_grid_test = 0
    j_grid_test = 12.5
    while i_grid_test <= 10000
        append!(x_hep_test, i_grid_test)
        i_grid_test += 50
    end
    while j_grid_test <= 975
        append!(y_hep_test, j_grid_test)
        j_grid_test += 12.5
    end
    probe = [Point(i, j) for i in x_hep_test, j in y_hep_test]
        
    lDa = zeros(Float64, 1, length(probe))
    prbDa = DataFrame(lDa, :auto)

    #Saving
    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"1d-topo-output"))do pvd
            u,ζ = uhn
            for (x,t) in x
                u,ζ = x
                Di = interpolate_everywhere(ζ, P(0.0))
                push!(prbDa, Di.(probe)) 
                pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>ζ,"h"=>h])
                println("done $t/$Tend")
            end
            prbDa = Matrix(prbDa)
            writedlm("Diffusion_zeta.csv", prbDa, ',')
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"1d-topo-output")) do pvd
            u,ζ = uhn
            for (x,t) in x
                u,ζ = x
                Di = interpolate_everywhere(ζ, P(0.0))
                push!(prbDa, Di.(probe)) 
                pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>ζ,"h"=>h])
                println("done $t/$Tend")
            end
            prbDa = Matrix(prbDa)
            writedlm("Diffusion_zeta.csv", prbDa, ',')
        end
    end
end

#Variable parameters
order = 1               #Order of the spaces
degree = 3              #Degree for the boundaries
Tend = 500              #The total time
dt = 5                  #Time step
B = 1000                #Channel width                               (m)
L = 10000               #Channel length                              (m)
x_points = 20           #Amount of points in the x-direction -> dx = 50 m
y_points = 100          #Amount of points in the y-direction -> dy = 100 m
dir = "./RESULTSSL"     #Name of the direction the files get savid in

#=
Shallow_water_equations_newton_solver input:
order                   Order of the spaces
degree                  Degree for the boundaries
Tend                    The total time
dt                      Time step 
B                       Channel width
L                       Channel length
x_points                Amount of points in the x-direction
y_points                Amount of points in the y-direction
dir                     Name of the direction the files get savid in
=#

#Initialization to speed up second time function is called
Shallow_water_equations_newton_solver(order,degree,2*dt,dt,B,L,x_points,y_points,dir)

#Real function
Shallow_water_equations_newton_solver(order,degree,Tend,dt,B,L,x_points,y_points,dir)


