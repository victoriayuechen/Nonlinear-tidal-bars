using WriteVTK
using Gridap


##''''''''''''''Save the beauty''''''''''''''##
function writing_output(dir, x, Ω, Tend)
    if isdir(dir)
        output_file = paraview_collection(joinpath(dir,"1d-topo-output"))do pvd
            for (x,t) in x
                # if t==Tend
                #     global xn = x 
                # end
                # if t==5000
                #     global x5000 = x
                #     println("done $t/$Tend")
                # end
                # if t==10000
                #     global x10000 = x
                #     println("done $t/$Tend")
                # end
                # if t==15000
                #     global x15000 = x
                #     println("done $t/$Tend")
                # end
                # if t==20000
                #     global x20000 = x
                #     println("done $t/$Tend")
                # end
                # if t==25000
                #     global x25000 = x
                #     println("done $t/$Tend")
                # end
                # if t==30000
                #     global x30000 = x
                #     println("done $t/$Tend")
                # end
                # if t==35000
                #     global x35000 = x
                #     println("done $t/$Tend")
                # end
                # if t==40000
                #     global x40000 = x
                #     println("done $t/$Tend")
                # end
                # if t==45000
                #     global x45000 = x
                # end

                if t%(400) ==0
                    u,ζ = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>ζ,"h"=>h])
                    println("Saved")
                    println("done $t/$Tend")
                end
                
            end
        end
    else
        mkdir(dir)
        output_file = paraview_collection(joinpath(dir,"1d-topo-output")) do pvd
            for (x,t) in x
                if t%(400) ==0
                    u,ζ = x
                    pvd[t] = createvtk(Ω,joinpath(dir,"1d-topo$t.vtu"),cellfields=["u"=>u,"ζ"=>ζ,"h"=>h])
                    println("Saved")
                    println("done $t/$Tend")
                end

            end
        end
    end
end