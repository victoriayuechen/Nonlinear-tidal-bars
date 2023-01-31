using Pkg
Pkg.activate(".")

using Gridap
using WriteVTK
using LineSearches: BackTracking
using GridapDistributed
using PartitionedArrays
using GridapPETSc
using GridapPETSc: PETSC

function perp(u)
    p = VectorValue(-u[2],u[1])
    p
end

function uh(u₀,h₀,X,Y,dΩ)
    a((u,h),(w,ϕ)) = ∫(w⋅u +ϕ*h)dΩ
    b((w,ϕ)) = ∫(w⋅u₀ + ϕ*h₀)dΩ
    solve(AffineFEOperator(a,b,X,Y))
end

function h₀((x,y))
    h = 0.05*exp(-0.01*(x-50)^2 -0.01*(y-25)^2)
    h
end

function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

function topography((x,y))
    b = 0.4*exp(-0.01*(y-50)^2)
    b
end

order = 0 
degree = 3 

options = """ 
    -mm_ksp_type cg
    -mm_ksp_monitor
    -mm_ksp_rtol 1.0e-10
    -mm_pc_type jacobi
    -log_view
"""

function linear_swe_topo(parts)
    GridapPETSc.with(args=split(options)) do
        # Parameters
        B = 100 #Channel width
        L = 100 #Channel Length
        latitude = 52
        η = 7.29e-5
        f = 2*η*sin(latitude*(π/180))
        g = 9.81
        H = 0.5 # Constant layer depth at rest

        # Make model
        domain = (0,B,0,L)
        mesh_partition = (5,5)
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
        dir = "output"
        Ω = Triangulation(model)
        dΩ = Measure(Ω,degree)
        dω = Measure(Ω,degree,ReferenceDomain())
        Γ = BoundaryTriangulation(model,tags=DC)
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
        l2(v) = ∫(v*h₀)dΩ
        hn = solve(AffineFEOperator(a2,l2,P,Q))

        a3(u,v) = ∫(v*u)dΩ
        l3(v) = ∫(v*topography)dΩ
        b = solve(AffineFEOperator(a3,l3,P,Q))

        uhn = uh(un,hn,X,Y,dΩ)
        un,hn = uhn
        #forcefunc(t) = VectorValue(0.0,0.0*0.5*cos(π*(1/5)*t))

        res(t,(u,h),(w,ϕ)) = ∫(∂t(u)⋅w + w⋅((perp∘(f*u.cellfield))) - (∇⋅(w))*g*h + ∂t(h)*ϕ + ϕ*(H-b)*(∇⋅(u))- ∇(b)⋅u)dΩ
        jac(t,(u,h),(du,dh),(w,ϕ)) = ∫(w⋅(f*(perp∘(du))) - (∇⋅(w))*g*dh + ϕ*((H-b)*(∇⋅(du))) - ∇(b)⋅u)dΩ
        jac_t(t,(u,h),(dut,dht),(w,ϕ)) = ∫(dut⋅w + dht*ϕ)dΩ

        op = TransientFEOperator(res,jac,jac_t,X,Y)
        nls = PETScLinearSolver()
        Tend = 2
        dt = 1.0
        ode_solver = ThetaMethod(nls, dt, 0.5)
        x = solve(ode_solver,op,uhn,0.0,Tend)

         #Output results
        if isdir(dir)
            createpvd(parts, joinpath(dir,"linear_topo")) do pvd
                for (x,t) in x
                    u,h = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"linear_topo$t.vtu"),cellfields=["u"=>u,"h"=>(h)])
                    println("done $t/$Tend")
                end
            end
        else
            mkdir(dir)
            createpvd(parts, joinpath(dir,"linear_topo")) do pvd
                for (x,t) in x
                    u,h = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"linear_topo$t.vtu"),cellfields=["u"=>u,"h"=>(h)])
                    println("done $t/$Tend")
                end
            end
        end

      

    end
end

partition = (2, 1)
prun(linear_swe_topo, mpi, partition) # use `sequential` instead of `mpi` for debugging 