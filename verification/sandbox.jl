using Plots


S1 = π*cos(π*t) - 2*H*π*cos(π*x)*sin(π*x)*sin(π*y)*exp(t) + 2*H*π*cos(π*x)*sin(π*x)*sin(π*y) - H*π*cos(π*x)*sin(π*y)*exp(t)*cos(π*y) + H*π*cos(π*x)*sin(π*y)*cos(π*y) - H*π*cos(π*x)*sin(π*y)*sin(π*t)*exp(t) + H*π*cos(π*x)*sin(π*y)*sin(π*t) - 2*H*π*cos(π*y)*sin(π*x)*sin(π*y)*exp(t) + 2*H*π*cos(π*y)*sin(π*x)*sin(π*y) - H*sin(π*x)*π*cos(π*y)*exp(t)*cos(π*x) + H*sin(π*x)*π*cos(π*y)*cos(π*x) - H*sin(π*x)*π*cos(π*y)*sin(π*t)*exp(t) + H*sin(π*x)*π*cos(π*y)*sin(π*t)

S2 = -π*cos(π*t)*sin(π*x)*sin(π*y)*exp(t) + π*cos(π*t)*sin(π*x)*sin(π*y) - sin(π*x)*sin(π*y)*exp(t)*cos(π*x) - sin(π*x)*sin(π*y)^2*exp(t) - sin(π*x)*sin(π*y)*exp(t)*sin(π*t) - f*sin(π*x)^2*sin(π*y)*exp(t) + f*sin(π*x)^2*sin(π*y) - f*sin(π*x)*sin(π*y)*exp(t)*cos(π*y) + f*sin(π*x)*sin(π*y)*cos(π*y) - f*sin(π*x)*sin(π*y)*sin(π*t)*exp(t) + f*sin(π*x)*sin(π*y)*sin(π*t) + g*π*cos(π*y)

S3 = π*cos(π*t) - 2*H*π*cos(π*x)*sin(π*x)*sin(π*y)*exp(t) + 2*H*π*cos(π*x)*sin(π*x)*sin(π*y) - H*π*cos(π*x)*sin(π*y)*exp(t)*cos(π*y) + H*π*cos(π*x)*sin(π*y)*cos(π*y) - H*π*cos(π*x)*sin(π*y)*sin(π*t)*exp(t) + H*π*cos(π*x)*sin(π*y)*sin(π*t) - 2*H*π*cos(π*y)*sin(π*x)*sin(π*y)*exp(t) + 2*H*π*cos(π*y)*sin(π*x)*sin(π*y) - H*sin(π*x)*π*cos(π*y)*exp(t)*cos(π*x) + H*sin(π*x)*π*cos(π*y)*cos(π*x) - H*sin(π*x)*π*cos(π*y)*sin(π*t)*exp(t) + H*sin(π*x)*π*cos(π*y)*sin(π*t)

H = 0.5 + sin(π*x) + sin(π*y) + sin(π*t)

U = (sin(π*x) + cos(π*y) + sin(π*t))*sin(π*x)*sin(π*y)*(-exp(t) + 1)

V = (cos(π*x) + sin(π*y) + sin(π*t))*sin(π*x)*sin(π*y)*(-exp(t) + 1)

u_m(x,t::Real) = VectorValue()
u_m(t::Real) = x -> g(x,t)

function u1_res((x, y, t))
    H_const = 0.5
    
    u_res = VectorValue(π*cos(π*t) - 2*H_const*π*cos(π*x)*sin(π*x)*sin(π*y)*exp(t) + 2*H_const*π*cos(π*x)*sin(π*x)*sin(π*y) - H_const*π*cos(π*x)*sin(π*y)*exp(t)*cos(π*y) + H_const*π*cos(π*x)*sin(π*y)*cos(π*y) - H_const*π*cos(π*x)*sin(π*y)*sin(π*t)*exp(t) + H_const*π*cos(π*x)*sin(π*y)*sin(π*t) - 2*H_const*π*cos(π*y)*sin(π*x)*sin(π*y)*exp(t) + 2*H_const*π*cos(π*y)*sin(π*x)*sin(π*y) - H_const*sin(π*x)*π*cos(π*y)*exp(t)*cos(π*x) + H_const*sin(π*x)*π*cos(π*y)*cos(π*x) - H_const*sin(π*x)*π*cos(π*y)*sin(π*t)*exp(t) + H_const*sin(π*x)*π*cos(π*y)*sin(π*t), S2 = -π*cos(π*t)*sin(π*x)*sin(π*y)*exp(t) + π*cos(π*t)*sin(π*x)*sin(π*y) - sin(π*x)*sin(π*y)*exp(t)*cos(π*x) - sin(π*x)*sin(π*y)^2*exp(t) - sin(π*x)*sin(π*y)*exp(t)*sin(π*t) - f*sin(π*x)^2*sin(π*y)*exp(t) + f*sin(π*x)^2*sin(π*y) - f*sin(π*x)*sin(π*y)*exp(t)*cos(π*y) + f*sin(π*x)*sin(π*y)*cos(π*y) - f*sin(π*x)*sin(π*y)*sin(π*t)*exp(t) + f*sin(π*x)*sin(π*y)*sin(π*t) + g*π*cos(π*y)

    )
end