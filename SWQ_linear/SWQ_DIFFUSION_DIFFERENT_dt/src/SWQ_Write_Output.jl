using WriteVTK
using Gridap


##''''''''''''''Save the beauty''''''''''''''##
function writing_output(dir, x, Ω, Tend)
    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"1d-topo-output"))do pvd
            for (x,t) in x
                if t==Tend
                    global xn = x
                end
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
                if t==Tend
                    global xn = x
                end

            end
        end
    end
end