using Pkg
Pkg.activate(".")
#create a new local environment, using my custom version of GridapPETSc
#] add https://github.com/carlodev/GridapPETSc.jl

using Gridap
using GridapDistributed
using PartitionedArrays
using WriteVTK
using LineSearches: BackTracking
using GridapPETSc
using GridapPETSc: PETSC

function perp(u)
    p = VectorValue(-u[2],u[1])
    p
end

function uζ(u₀,ζ₀,X,Y,dΩ)
    a((u,ζ),(w,ϕ)) = ∫(w⋅u +ϕ*ζ)dΩ
    b((w,ϕ)) = ∫(w⋅u₀ + ϕ*ζ₀)dΩ
    solve(AffineFEOperator(a,b,X,Y))
end

function ζ₀((x,y))
    h =0.05*exp(-0.01*(x-25)^2 -0.01*(y-25)^2)
    h
end

function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

function forcefunc((x,y),t)
    f = VectorValue(0.0,0.0*0.5*cos(π*(1/5)*t))
    f
end

order = 1
degree = 4

options = """ 
    -snes_type newtonls
    -snes_linesearch_type basic
    -snes_linesearch_damping 1.0
    -snes_rtol 1.0e-8
    -snes_atol 0.0
    -snes_monitor
    -snes_converged_reason
    -mm_ksp_type cg
    -mm_ksp_monitor
    -mm_ksp_rtol 1.0e-10
    -mm_pc_type jacobi
"""

function linear_swe(parts)
    GridapPETSc.with(args=split(options)) do
        #Parameters
        B = 100 
        L = 100 
        latitude = 52
        η = 7.29e-5
        f = 2*η*sin(latitude*(π/180))
        g = 9.81
        H = 0.5 #Constant layer depth at rest

        # Make model
        domain = (0,B,0,L)
        mesh_partition = (100,100)
        # Generate the model
        model = CartesianDiscreteModel(parts,domain,mesh_partition;isperiodic=(false,true))

        # Make labels
        labels = get_face_labeling(model)
        add_tag_from_tags!(labels,"bottom",[1,2,5])
        add_tag_from_tags!(labels,"left",[7])
        add_tag_from_tags!(labels,"right",[8])
        add_tag_from_tags!(labels,"top",[3,4,6])
        add_tag_from_tags!(labels,"inside",[9])
        DC = ["left","right"]
        dir = "output3/"
        Ω = Triangulation(model)
        dΩ = Measure(Ω,degree)
        dω = Measure(Ω,degree,ReferenceDomain())
        Γ = BoundaryTriangulation(model,tags=["top", "bottom"])
        nΓ = get_normal_vector(Γ)
        dΓ = Measure(Γ,degree)

        # Make reference spaces
        reffe_rt = ReferenceFE(raviart_thomas,Float64,order)
        V = TestFESpace(model,reffe_rt,conformity=:HDiv,dirichlet_tags=DC)
        udc(x,t::Real) = VectorValue(0.0,0.0)
        udc = x -> udc(x,t)
        U = TransientTrialFESpace(V)

        reffe_lg = ReferenceFE(lagrangian,Float64,order)
        Q = TestFESpace(model,reffe_lg,conformity=:L2)
        P = TransientTrialFESpace(Q)

        Y = MultiFieldFESpace([V,Q])
        X = MultiFieldFESpace([U,P])

        a1(u,v) = ∫(v⋅u)dΩ 
        l1(v) = ∫(v⋅u₀)dΩ
        un = solve(AffineFEOperator(a1,l1,U,V))

        a2(u,v) = ∫(v*u)dΩ
        l2(v) = ∫(v*ζ₀)dΩ
        ζn = solve(AffineFEOperator(a2,l2,P,Q))

        uζn = uζ(un,ζn,X,Y,dΩ)
        un,ζn = uζn

        forcefunc_x(t) = x -> forcefunc(x,t)

        res(t,(u,ζ),(w,ϕ)) = ∫(∂t(u)⋅w + w⋅((perp∘(f*u.cellfield))) - (∇⋅(w))*g*ζ + ∂t(ζ)*ϕ + ϕ*H*(∇⋅(u)) - forcefunc_x(t)⋅w)dΩ
        jac(t,(u,ζ),(du,dζ),(w,ϕ)) = ∫(w⋅(f*(perp∘(du))) - (∇⋅(w))*g*dζ + ϕ*(H*(∇⋅(du))))dΩ
        jac_t(t,(u,ζ),(dut,dζt),(w,ϕ)) = ∫(dut⋅w + dζt*ϕ)dΩ

        op = TransientFEOperator(res,jac,jac_t,X,Y)
        nls = PETScNonlinearSolver()
        Tend = 30
        ode_solver = ThetaMethod(nls,1.0,0.5)
        x = solve(ode_solver,op,uζn,0.0,Tend)

        if isdir(dir)
            createpvd(parts,joinpath(dir,"third_linear_topo")) do pvd
                for (x,t) in x
                    u,ζ = x
                    pvd[t] = createvtk(Ω, "third_linear_topo_$t"*".vtu", cellfields=["u"=>u,"ζ+H "=>ζ+H])
                end
            end
        else
            mkdir(dir)
            createpvd(parts,joinpath(dir,"linear_topo")) do pvd
                for (x,t) in x
                    u,ζ = x
                    pvd[t] = createvtk(Ω, "third_linear_topo_$t"*".vtu", cellfields=["u"=>u,"ζ+H "=>ζ+H])
                end
            end
        end
    end
end

partition = (2, 2)
prun(linear_swe, mpi, partition) # use `sequential` instead of `mpi` for debugging 