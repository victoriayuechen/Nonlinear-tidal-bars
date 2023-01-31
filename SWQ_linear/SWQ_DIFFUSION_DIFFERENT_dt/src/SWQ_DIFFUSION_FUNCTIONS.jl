using Gridap
using Gridap.TensorValues: meas
using Parameters

includet("SWQ_Parameters.jl")


##''''''''''Different functions with their derivitave''''''''##
#Velocity change 
function func_ut(t,(u,ζ),(w,ϕ)) 
    ut = ∂t(u)⋅w
    return ut
end

#Convection
function func_con(t,(u,ζ),(w,ϕ),Param) 
    @unpack convection_on = Param
    if convection_on
        con = ∇(u)'⋅u⋅w
    else
        con = 0.0
    end
    return con
end

function dfunc_con(t,(u,ζ),(du,dζ),(w,ϕ),Param) 
    @unpack convection_on = Param
    if convection_on
        dcon = ∇(du)'⋅u⋅w + ∇(u)'⋅du⋅w
    else
        dcon = 0.0
    end
    return dcon
end

#Coriolis
function func_cor(t,(u,ζ),(w,ϕ),Param) 
    @unpack f, coriolis_on = Param
    if coriolis_on
        cor = f*(coriolis∘(u))⋅w
    else
        cor = 0.0
    end
    return cor
end
function dfunc_cor(t,(u,ζ),(du,dζ),(w,ϕ),Param) 
    @unpack f, coriolis_on = Param
    if coriolis_on
        dcor = f*(coriolis∘(du))⋅w
    else
        dcor = 0.0
    end
    return dcor
end


#Drag coefficient term
function func_cD(t,(u,ζ),(w,ϕ),Param, h) 
    @unpack cD, H, U_start, friction_on, linear_friction_on = Param
    if friction_on
        if linear_friction_on
            r = 8/3π * cD * U_start
            fu_cD = r*u⋅w/(ζ+H-h)
        else
            fu_cD = cD* (meas∘u) * u⋅w/(ζ+H-h)
        end
    else
        fu_cD = 0.0
    end
    return fu_cD
end
function dfunc_cD(t,(u,ζ),(du,dζ),(w,ϕ),Param, h) 
    @unpack cD, H, U_start, friction_on, linear_friction_on = Param
    if friction_on
        if linear_friction_on
            r = 8/3π * cD * U_start
            dfu_cD = r*du⋅w/(ζ+H-h)-r*u⋅w/((ζ+H-h)*(ζ+H-h))*dζ
        else
            dfu_cD =  cD/(ζ+H-h)*w⋅(-(meas∘u)*u*dζ/(ζ+H-h)+u⋅du*u/((meas∘(u+1e-14))) + (meas∘u)*du)
        end
    else
        dfu_cD = 0.0
    end    
    return dfu_cD
end

#Gravitational (without boundary)
function func_g(t,(u,ζ),(w,ϕ),Param) 
    @unpack g, gravitional_on = Param
    if gravitional_on
        fu_g = - g*(∇⋅(w))*ζ
    else
        fu_g = 0.0
    end        
    return fu_g
end
function dfunc_g(t,(u,ζ),(du,dζ),(w,ϕ),Param) 
    @unpack g, gravitional_on = Param
    if gravitional_on
        dfu_g =  -g*(∇⋅(w))*dζ
    else
        dfu_g = 0.0
    end
    return dfu_g
end

#Forcing function
function func_Fₚ(t,(u,ζ),(w,ϕ),Param) 
    @unpack forcing_on = Param
    if forcing_on
        Fₚ = - forcfunc(t,Param)⋅w
    else
        Fₚ = 0.0
    end
    return Fₚ
end

#ζ change
function func_ζt(t,(u,ζ),(w,ϕ)) 
    ζt = ∂t(ζ)*ϕ
    return ζt
end

#ζ+Velocity function
function func_h(t,(u,ζ),(w,ϕ),Param, h) 
    @unpack H, momentum_on, linear_on = Param
    if momentum_on
        if linear_on
            fu_h =  -(ζ+H-h)*u ⋅(∇(ϕ))
        else
            fu_h = - H * u ⋅(∇(ϕ))
        end
    else
        fu_h = 0.0
    end
    return fu_h
end
function dfunc_h(t,(u,ζ),(du,dζ),(w,ϕ),Param, h) 
    @unpack H, momentum_on, linear_on = Param
    if momentum_on
        if linear_on
            dfu_h = - H * du ⋅(∇(ϕ))
        else
            dfu_h = - (dζ*u+(ζ+H-h)*u) ⋅(∇(ϕ))
        end
    else
        dfu_h = 0.0
    end 
    return dfu_h
end

#Boundary
function func_boun(t,(u,ζ),(w,ϕ),Param,nΓ) 
    @unpack g, boundary_on = Param
    if boundary_on
        boun = g*(ζ)*(w⋅nΓ)
    else
        boun = 0.0
    end
    return boun
end
function dfunc_boun(t,(u,ζ),(du,dζ),(w,ϕ),Param,nΓ) 
    @unpack g, boundary_on = Param
    if boundary_on
        dboun = g*(dζ)*(w⋅nΓ)
    else
        dboun = 0.0
    end
    return dboun
end

#Stabilization function ζ
function func_stabζ(t,(u,ζ),(w,ϕ),Param) 
    @unpack α, stabilization_ζ = Param
    if stabilization_ζ
        stabζ = α*(∇(ζ)⋅∇(ϕ)) 
    else
        stabζ = 0.0
    end
    return stabζ
end                                              
function dfunc_stabζ(t,(u,ζ),(du,dζ),(w,ϕ),Param)
    @unpack α, stabilization_ζ = Param 
    if stabilization_ζ
        dstabζ = α*(∇(dζ)⋅∇(ϕ)) 
    else
        dstabζ = 0.0
    end
    return dstabζ
end                                  


#Stabilization function u
function func_stabu(t,(u,ζ),(w,ϕ),Param) 
    @unpack ν, stabilization_u = Param
    if stabilization_u
        stabu = ν * (∇⋅u)*(∇⋅w)   
    else
        stabu = 0.0  
    end       
    return stabu
end

function dfunc_stabu(t,(u,ζ),(du,dζ),(w,ϕ),Param) 
    @unpack ν, stabilization_u = Param
    if stabilization_u
        dstabu = ν * (∇⋅du)*(∇⋅w)
    else
        dstabu = 0.0
    end
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
    forc = VectorValue(σ*U_start*cos(σ*t)+cD/H*abs(U_start*sin(σ*t))*U_start*sin(σ*t),f*U_start*sin(σ*t))
    return forc
end