using Gridap
using Parameters
includet("SWQ_Parameters_APVM.jl")

##''''''''''''''''Define initial solution''''''''''''''''##
function uD(u₀,D₀,F₀,q₀,X,Y,dΩ)
    a((u,D,r,u2),(w,ϕ,ϕ2,w2)) = ∫(w⋅u + ϕ*D + w2⋅u2 + ϕ2*r)dΩ
    b((w,ϕ,ϕ2,w2)) = ∫(w⋅u₀ + ϕ*D₀ + w2⋅F₀ + q₀*ϕ2)dΩ
    solve(AffineFEOperator(a,b,X(0.0),Y))
end

function compute_potential_vorticity(D,f,u,R,S,dΩ)
    a(r,s) = ∫(s*D*r)dΩ
    c(s)   = ∫(perp∘(∇(s))⋅(u) + s*f)dΩ
    solve(AffineFEOperator(a,c,R(0.0),S))
end
function compute_mass_flux(F,U,V,dΩ)
    a(DU,w) = ∫(DU⋅w)dΩ
    b(w) = ∫(w⋅F)dΩ
    solve(AffineFEOperator(a,b,U(0.0),V))
end

function init_sol(X,P,R,S,U,V,Y,dΩ,Param)
    @unpack u_y,u_x,ζₙ,f,H = Param
    u₀ = VectorValue(u_x,u_y)
    un = interpolate_everywhere(u₀,U(0.0))

    find_h(S,Param)
    D₀ = ζₙ + H - h
    Dn = interpolate_everywhere(D₀,P(0.0))

    # fn = interpolate_everywhere(f,S(0.0))
    F₀ = compute_mass_flux(un*Dn,U,V,dΩ)
    q₀ = compute_potential_vorticity(D₀,f,u₀,R,S,dΩ)
    uhn = uD(un,Dn,F₀,q₀,X,Y,dΩ)
    global fn = interpolate_everywhere(f,S(0.0))


    return uhn
end

function find_h(S,Param)
    @unpack B, L = Param
    h₀(x,t) = 0.1 * cos(π/B * x[2]) * cos(2π/L * x[1])      #Hepkema original function for h
    # h₀(x,t) = 0.0
    h₀(t::Real) = x->h₀(x,t)
    global h = interpolate_everywhere(h₀(0.0),S(0.0))
end
