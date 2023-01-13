using Gridap
using LineSearches: BackTracking
using Gridap.TensorValues: meas
using Parameters

includet("SWQ_DIFFUSION_FUNCTIONS.jl")
includet("SWQ_Parameters.jl")

##''''''''''''''''Define solver''''''''''''''''##
function solver(Initial_conditions, Y, X, dt, Tstart, Tend, dΩ, dΓ, show_result,Param)
    res(t,(u,ζ),(w,ϕ)) = ∫(func_ut(t,(u,ζ),(w,ϕ)) +                     #Velocity change
        func_con(t,(u,ζ),(w,ϕ)) +                                       #Convection
        func_cor(t,(u,ζ),(w,ϕ),Param) +                                       #Coriolis
        func_cD(t,(u,ζ),(w,ϕ),Param)  +                                       #Drag coefficient term
        func_g(t,(u,ζ),(w,ϕ),Param)   +                                       #Gravitational
        func_Fₚ(t,(u,ζ),(w,ϕ),Param)  +                                       #Forcing function
        func_ζt(t,(u,ζ),(w,ϕ))  +                                       #ζ change
        func_h(t,(u,ζ),(w,ϕ),Param)   +                                       #ζ+Velocity function
        func_stabζ(t,(u,ζ),(w,ϕ),Param) +                                     #Stabilization ζ
        func_stabu(t,(u,ζ),(w,ϕ),Param))dΩ +                                  #Stabilization u
        ∫(func_boun(t,(u,ζ),(w,ϕ),Param))dΓ                                   #Boundary

##''''''''''''''Jacobian''''''''''''''##
    jac(t,(u,ζ),(du,dζ),(w,ϕ)) = ∫(dfunc_con(t,(u,ζ),(du,dζ),(w,ϕ)) +   #Convection
        dfunc_cor(t,(u,ζ),(du,dζ),(w,ϕ),Param) +                              #Coriolis
        dfunc_cD(t,(u,ζ),(du,dζ),(w,ϕ),Param)  +                              #Drag coefficient term
        dfunc_g(t,(u,ζ),(du,dζ),(w,ϕ),Param)   +                              #Gravitational
        dfunc_h(t,(u,ζ),(du,dζ),(w,ϕ),Param)   +                              #ζ+Velocity function
        dfunc_stabζ(t,(u,ζ),(du,dζ),(w,ϕ),Param) +                            #Stabilization ζ
        dfunc_stabu(t,(u,ζ),(du,dζ),(w,ϕ),Param))dΩ +                         #Stabilization u
        ∫(dfunc_boun(t,(u,ζ),(du,dζ),(w,ϕ),Param))dΓ


    ##''''''''''''''Jacobian for t''''''''''''''##
    jac_t(t,(u,ζ),(dut,dζt),(w,ϕ)) = ∫(dut⋅w + dζt*ϕ)dΩ

    ##''''''''''''''Solver with ThetaMethod''''''''''''''##
    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=show_result,method= :newton,linesearch=BackTracking())
    ode_solver = ThetaMethod(nls,dt,0.5)
    x = solve(ode_solver,op,Initial_conditions,Tstart,Tend)
    return x
end