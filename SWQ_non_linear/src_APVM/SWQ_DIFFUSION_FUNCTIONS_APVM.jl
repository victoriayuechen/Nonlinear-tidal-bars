using Gridap
using Gridap.TensorValues: meas
using Parameters

includet("SWQ_Parameters_APVM.jl")


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
function func_cor(t,(u,ζ),(w,ϕ),Param) 
    @unpack f = Param
    cor = f*(coriolis∘(u))⋅w
    return cor
end
function dfunc_cor(t,(u,ζ),(du,dζ),(w,ϕ),Param) 
    @unpack f = Param
    dcor = f*(coriolis∘(du))⋅w
    return dcor
end


#Drag coefficient term
function func_cD(t,(u,ζ),(w,ϕ),Param) 
    @unpack cD, H = Param
    fu_cD = cD* (meas∘u) * u⋅w/(ζ+H-h)
    return fu_cD
end
function dfunc_cD(t,(u,ζ),(du,dζ),(w,ϕ),Param) 
    @unpack cD, H = Param
    dfu_cD =  cD/(ζ+H-h)*w⋅(-(meas∘u)*u*dζ/(ζ+H-h)+u⋅du*u/((meas∘(u+1e-14))) + (meas∘u)*du)
    return dfu_cD
end

#Gravitational (without boundary)
function func_g(t,(u,ζ),(w,ϕ),Param) 
    @unpack g = Param
    fu_g = - g*(∇⋅(w))*ζ
    return fu_g
end
function dfunc_g(t,(u,ζ),(du,dζ),(w,ϕ),Param) 
    @unpack g = Param
    dfu_g =  -g*(∇⋅(w))*dζ
    return dfu_g
end

#Forcing function
function func_Fₚ(t,(u,ζ),(w,ϕ),Param) 
    Fₚ = - forcfunc(t,Param)⋅w
    return Fₚ
end

#ζ change
function func_ζt(t,(u,ζ),(w,ϕ)) 
    ζt = ∂t(ζ)*ϕ
    return ζt
end
#ζ+Velocity function
function func_h(t,(u,ζ),(w,ϕ),Param) 
    @unpack H = Param
    fu_h =  -(ζ+H-h)*u ⋅(∇(ϕ))
    return fu_h
end
function dfunc_h(t,(u,ζ),(du,dζ),(w,ϕ),Param) 
    @unpack H = Param
    dfu_h = - (dζ*u+(ζ+H-h)*u) ⋅(∇(ϕ))
    return dfu_h
end

#Boundary
function func_boun(t,(u,ζ),(w,ϕ),Param) 
    @unpack g = Param
    boun = g*(ζ)*(w⋅nΓ)
    return boun
end
function dfunc_boun(t,(u,ζ),(du,dζ),(w,ϕ),Param) 
    @unpack g = Param
    dboun = g*(dζ)*(w⋅nΓ)
    return dboun
end

#Stabilization function ζ
function func_stabζ(t,(u,ζ),(w,ϕ),Param) 
    @unpack α = Param
    stabζ = α*(∇(ζ)⋅∇(ϕ)) 
    return stabζ
end                                              
function dfunc_stabζ(t,(u,ζ),(du,dζ),(w,ϕ),Param)
    @unpack α = Param 
    dstabζ = α*(∇(dζ)⋅∇(ϕ)) 
    return dstabζ
end                                  


#Stabilization function u
function func_stabu(t,(u,ζ),(w,ϕ),Param) 
    @unpack ν = Param
    stabu = ν * (∇⋅u)*(∇⋅w)            
    return stabu
end

function dfunc_stabu(t,(u,ζ),(du,dζ),(w,ϕ),Param) 
    @unpack ν = Param
    dstabu = ν * (∇⋅du)*(∇⋅w)
    return dstabu
end

#Coriolis
function coriolis(u) 
    corio = VectorValue(-u[2],u[1]) 
    return corio                                                                            
end

#forcfunc
function forcfunc(t,Param) 
    @unpack f, U_start, σ, cD, H = Param
    # forc = VectorValue(-f*U_start*cos(σ*t),-σ*U_start*sin(σ*t)+cD/H*abs(U_start*cos(σ*t))*U_start*cos(σ*t))      #Fₚ from Hepkema
    forc = VectorValue(σ*U_start*cos(σ*t)+cD/H*abs(U_start*sin(σ*t))*U_start*sin(σ*t),f*U_start*sin(σ*t))
    return forc
end