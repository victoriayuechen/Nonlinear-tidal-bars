using Gridap
function perp(u)
p = VectorValue(-u[2],u[1])
p
end
function ζ0((x,y))
h =0.01*exp(-0.1*(x-30)^2 -0.1*(y-50)^2)
h
end
function u0((x,y))
u = VectorValue(0.0,0.0)
u
end
function forcefunc((x,y),t)
f = VectorValue(0.0,0.0)
f
end


Tend = 20
dt = 0.1
B = 100
L = 100
partition = (100,100)
domain = (0,L,0,B)
order = 1
degree = 4
latitude = 50
η = 7.29e-5
f = 2*η*sin(latitude*(π/180))
g = 9.81
H = 0.5
model = CartesianDiscreteModel(domain,partition;isperiodic=(false,false))
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"bottom",[1,2,5])
add_tag_from_tags!(labels,"left",[7])
add_tag_from_tags!(labels,"right",[8])
add_tag_from_tags!(labels,"top",[3,4,6])
add_tag_from_tags!(labels,"inside",[9])
DC = ["left","right","top","bottom"]


Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
reffe_rt = ReferenceFE(raviart_thomas,Float64,order)
V = TestFESpace(model,reffe_rt,dirichlet_tags=DC,conformity=:HDiv)
udc(x,t::Real) = VectorValue(0.0,0.0)
udc(t::Real) = x -> udc(x,t)
U = TransientTrialFESpace(V,[udc,udc,udc,udc])
reffe_lg = ReferenceFE(lagrangian,Float64,order)
Q = TestFESpace(model,reffe_lg,conformity=:L2)
P = TransientTrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = TransientMultiFieldFESpace([U,P])

x0 = interpolate_everywhere([u0,ζ0],X(0.0))
un,ζn = x0
forcefunc_x(t) = x -> forcefunc(x,t)

res(t,(u,ζ),(w,φ)) = ∫(∂t(u)·w + w·(f*(perp∘(u)))- (∇·(w))*g*ζ + ∂t(ζ)*φ + φ*H*(∇·(u)) - forcefunc_x(t)·w)dΩ
jac(t,(u,ζ),(du,dζ),(w,φ)) = ∫(w·(f*(perp∘(du))) - (∇·(w))*g*dζ + φ*(H*(∇·(du))))dΩ
jac_t(t,(u,ζ),(dut,dζt),(w,φ)) = ∫(dut·w + dζt*φ)dΩ

op = TransientFEOperator(res,jac,jac_t,X,Y)
nls = LUSolver()
ode_solver = ThetaMethod(nls,dt,0.5)
x = solve(ode_solver,op,x0,0.0,Tend)

createpvd("linear_SWE") do pvd
    for (x,t) in x
    u,ζ = x
    pvd[t] = createvtk(Ω,"linear_SWE_$t.vtu",cellfields = ["u"=>u,"ζ"=>(ζ)])
    println("Done $t/Tend")
    end
end