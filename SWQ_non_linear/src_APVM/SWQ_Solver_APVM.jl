using Gridap
# using LineSearches: BackTracking
using LineSearches
using Gridap.TensorValues: meas
using Parameters

# includet("SWQ_DIFFUSION_FUNCTIONS.jl")
includet("SWQ_Parameters_APVM.jl")
function forcfunc(t,Param) 
    @unpack f, U_start, σ, cD, H = Param
    # forc = VectorValue(-f*U_start*cos(σ*t),-σ*U_start*sin(σ*t)+cD/H*abs(U_start*cos(σ*t))*U_start*cos(σ*t))      #Fₚ from Hepkema
    forc = VectorValue(σ*U_start*cos(σ*t)+cD/H*abs(U_start*sin(σ*t))*U_start*sin(σ*t),f*U_start*sin(σ*t))
    return forc
end

##''''''''''''''''Define solver''''''''''''''''##
function solver(Initial_conditions, Y, X, dt, Tstart, Tend, theta, dΩ, dΓ, show_result,Param)
    #Forcing function on u(t)

    #Norm operator
    norm(u) = meas∘(u) + 1e-14
    dnorm(u,du) = u ⋅ du / norm(u)

    #Parameters
    cd = 0.0025
    g = 9.81
    τ = dt*0.5

    #Define residual, jacobian in space and time
    res(t,(u,D,q,F),(w,ϕ,ϕ2,w2)) = ∫(∂t(u)⋅w + ((q-τ*(u⋅∇(q)))*(perp∘(F)))⋅w - (∇⋅(w))*(g*(D+h) + 0.5*(u⋅u)) + ∂t(D)*ϕ + w2⋅(F-u*D) + ϕ2*(q*D) +  (perp∘(∇(ϕ2)))⋅u - (ϕ2*fn) + ϕ*(∇⋅(F))- forcfunc(t,Param)⋅w + ((cd*(norm(u))*u)/(D+1e-14))⋅w)dΩ + ∫((g*(D+h) + 0.5*(u⋅u))*(w⋅nΓ) + nΓ⋅(perp∘(u))*ϕ2)dΓ 
    jac(t,(u,D,q,F),(du,dD,dq,dF),(w,ϕ,ϕ2,w2)) = ∫(((dq - τ*(u⋅∇(dq) + du⋅∇(q)))*(perp∘(F)))⋅w + ((q-τ*(u⋅∇(q)))*(perp∘(dF)))⋅w - (∇⋅(w))*(g*(dD) + (du⋅u)) + w2⋅(dF - du*D - dD*u) + ϕ2*(q*dD) + ϕ2*(dq*D) + du⋅(perp∘(∇(ϕ2))) + ϕ*(∇⋅(dF)) + (cd*(dnorm(u,du)*u)/(D+h+1e-14) + cd*(norm(u)/(D+h+1e-14)*du) - cd*(norm(u)*u*dD)/((D+h+1e-14)*(D+h+1e-14)))⋅w )dΩ+ ∫((g*(dD) + (du⋅u))*(w⋅nΓ) + nΓ⋅(perp∘(du))*ϕ2)dΓ
    jac_t(t,(u,D),(dut,dDt),(w,ϕ)) = ∫(dut⋅w + dDt*ϕ)dΩ

    #Define operators and solvers
    op = TransientFEOperator(res,jac,jac_t,X,Y)
    nls = NLSolver(show_trace=show_result,method =:trust_region)
    ode_solver = ThetaMethod(nls,dt,theta)
    x = solve(ode_solver,op,Initial_conditions,Tstart,Tend)
    return x
end