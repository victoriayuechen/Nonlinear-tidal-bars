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

##''''''''''Different functions with their derivitave''''''''##
#Velocity change 
function func_ut(t,(u,ζ),(w,ϕ)) 
    ut = ∂t(u)⋅w
    return ut
end

#Convection
function func_con(t,(u,ζ),(w,ϕ)) 
    con = ∇(u)'⋅u⋅w
    return con
end

function dfunc_con(t,(u,ζ),(du,dζ),(w,ϕ)) 
    dcon = ∇(du)'⋅u⋅w + ∇(u)'⋅du⋅w
    return dcon
end

#Coriolis
function func_cor(t,(u,ζ),(w,ϕ),f) 
    cor = (coriolis(u,f))⋅w
    return cor
end
function dfunc_cor(t,(u,ζ),(du,dζ),(w,ϕ),f) 
    dcor = (coriolis(du,f))⋅w
    return dcor
end


#Drag coefficient term
function func_cD(t,(u,ζ),(w,ϕ),cD,H,h) 
    fu_cD = cD* (meas∘u) * u⋅w/(ζ+H-h)
    return fu_cD
end
function dfunc_cD(t,(u,ζ),(du,dζ),(w,ϕ),cD,H,h) 
    dfu_cD =  cD/(ζ+H-h)*w⋅(-(meas∘u)*u*dζ/(ζ+H-h)+u⋅du*u/((meas∘(u+1e-14))) + (meas∘u)*du)
    return dfu_cD
end

#Gravitational (without boundary)
function func_g(t,(u,ζ),(w,ϕ),g) 
    fu_g = - g*(∇⋅(w))*ζ
    return fu_g
end
function dfunc_g(t,(u,ζ),(du,dζ),(w,ϕ),g) 
    dfu_g =  -g*(∇⋅(w))*dζ
    return dfu_g
end

#Forcing function
function func_Fₚ(t,(u,ζ),(w,ϕ),f,U_start,σ) 
    Fₚ = - forcfunc(t,f,U_start,σ)⋅w
    return Fₚ
end

#ζ change
function func_ζt(t,(u,ζ),(w,ϕ)) 
    ζt = ∂t(ζ)*ϕ
    return ζt
end
#ζ+Velocity function
function func_h(t,(u,ζ),(w,ϕ),H,h) 
    fu_h =  -(ζ+H-h)*u ⋅(∇(ϕ))
    return fu_h
end
function dfunc_h(t,(u,ζ),(du,dζ),(w,ϕ),H,h) 
    dfu_h = - (dζ*u+(ζ+H-h)*u) ⋅(∇(ϕ))
    return dfu_h
end

#Boundary
function func_boun(t,(u,ζ),(w,ϕ),g,nΓ) 
    boun = g*(ζ)*(w⋅nΓ)
    return boun
end
function dfunc_boun(t,(u,ζ),(du,dζ),(w,ϕ),g,nΓ) 
    dboun = g*(dζ)*(w⋅nΓ)
    return dboun
end

#Stabilization function ζ
function func_stabζ(t,(u,ζ),(w,ϕ),α) 
    stabζ = α*(∇(ζ)⋅∇(ϕ)) 
    return stabζ
end                                              
function dfunc_stabζ(t,(u,ζ),(du,dζ),(w,ϕ),α) 
    dstabζ = α*(∇(dζ)⋅∇(ϕ)) 
    return dstabζ
end                                  


#Stabilization function u
function func_stabu(t,(u,ζ),(w,ϕ),ν) 
    stabu = ν * (∇⋅u)*(∇⋅w)            
    return stabu
end

function dfunc_stabu(t,(u,ζ),(du,dζ),(w,ϕ),ν) 
    dstabu = ν * (∇⋅du)*(∇⋅w)
    return dstabu
end

#Coriolis
function coriolis(u,f) 
    corio = 0.0 #f*VectorValue(u[2],-u[1]) 
    return corio                                                                            
end

#forcfunc
function forcfunc(t,f,U_start,σ) 
    forc = 0.0 #VectorValue(-f*U_start*cos(σ*t),-σ*U_start*sin(σ*t)+cD*abs(U_start*cos(σ*t))*U_start*cos(σ*t))      #Fₚ from Hepkema
    return forc
end


function Make_model(B,L,x_points,y_points,order,degree)
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
    return Ω, dΩ, nΓ, dΓ, Y, X, P
end

##''''''''''''''''Define initial solution''''''''''''''''##
function init_sol(u₀,ζ₀,X)
    uhn = interpolate_everywhere([u₀(0.0),ζ₀(0.0)],X(0.0))
    return uhn
end

function find_h(h₀)
    h = interpolate_everywhere(h₀(0.0),P(0.0))
    return h
end


##''''''''''''''''Define solver''''''''''''''''##
function solver(Initial_conditions, latitude, Y, X, dt, Tstart, Tend, dΩ, dΓ, nΓ, h, show_result)




    res(t,(u,ζ),(w,ϕ)) = ∫(func_ut(t,(u,ζ),(w,ϕ)) +                     #Velocity change
        func_con(t,(u,ζ),(w,ϕ)) +                                       #Convection
        func_cor(t,(u,ζ),(w,ϕ),f) +                                       #Coriolis
        func_cD(t,(u,ζ),(w,ϕ),cD,H,h)  +                                       #Drag coefficient term
        func_g(t,(u,ζ),(w,ϕ),g)   +                                       #Gravitational
        func_Fₚ(t,(u,ζ),(w,ϕ),f,U_start,σ)  +                                       #Forcing function
        func_ζt(t,(u,ζ),(w,ϕ))  +                                       #ζ change
        func_h(t,(u,ζ),(w,ϕ),H,h)   +                                       #ζ+Velocity function
        func_stabζ(t,(u,ζ),(w,ϕ),α) +                                     #Stabilization ζ
        func_stabu(t,(u,ζ),(w,ϕ)),ν)dΩ +                                  #Stabilization u
        ∫(func_boun(t,(u,ζ),(w,ϕ),g,nΓ))dΓ                                   #Boundary

    ##''''''''''''''Jacobian''''''''''''''##
    jac(t,(u,ζ),(du,dζ),(w,ϕ)) = ∫(dfunc_con(t,(u,ζ),(du,dζ),(w,ϕ)) +   #Convection
        dfunc_cor(t,(u,ζ),(du,dζ),(w,ϕ),f) +                              #Coriolis
        dfunc_cD(t,(u,ζ),(du,dζ),(w,ϕ),cD,H,h)  +                              #Drag coefficient term
        dfunc_g(t,(u,ζ),(du,dζ),(w,ϕ),g)   +                              #Gravitational
        dfunc_h(t,(u,ζ),(du,dζ),(w,ϕ),H,h)   +                              #ζ+Velocity function
        dfunc_stabζ(t,(u,ζ),(du,dζ),(w,ϕ),α) +                            #Stabilization ζ
        dfunc_stabu(t,(u,ζ),(du,dζ),(w,ϕ)),ν)dΩ +                         #Stabilization u
        ∫(dfunc_boun(t,(u,ζ),(du,dζ),(w,ϕ),g,nΓ))dΓ                          #Boundary

    ##''''''''''''''Jacobian for t''''''''''''''##
    jac_t(t,(u,ζ),(dut,dζt),(w,ϕ)) = ∫(dut⋅w + dζt*ϕ)dΩ

    ##''''''''''''''Solver with ThetaMethod''''''''''''''##
    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=show_result,method= :newton,linesearch=BackTracking())
    ode_solver = ThetaMethod(nls,dt,0.5)
    x = solve(ode_solver,op,Initial_conditions,Tstart,Tend)
    return x
end

##''''''''''''''Save the beauty''''''''''''''##
function writing_output(dir, x, Ω)
    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"1d-topo-output"))do pvd
            for (x,t) in x
                if t%(500) ==0
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
                if t%(500) ==0
                    u,ζ = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>ζ,"h"=>h])
                    println("Saved")
                end
                println("done $t/$Tend")
            end
        end
    end
end

B = 1000
L = 10000
x_points = 20
y_points = 100
order = 1
degree = 3
global ν = 1e-6                                #From Anna Louka

@time Ω, dΩ, nΓ, dΓ, Y, X, P = Make_model(B,L,x_points,y_points,order,degree)

function u₀(t) 
    u = VectorValue(0.0,0.0)
    return u
end

function ζ₀(t)
    ζ = 0.0
    return ζ
end

h₀(x,t) = 0.1 * cos(π/B * x[1]) * cos(2π/L * x[2])      #Hepkema original function for h
h₀(t::Real) = x->h₀(x,t)
h = find_h(h₀)
# u₀(t) = VectorValue(0.0,0.0) evenly fast
# ζ₀(t) = 0.0

@time Initial_solutions = init_sol(u₀,ζ₀,X)

latitude = 52
dt = 5
Tstart = 0
Tend = 50
@time x = solver(Initial_solutions, latitude, Y, X, dt, Tstart, Tend, dΩ, dΓ, nΓ, h, true)

@time writing_output(dir, x, Ω)