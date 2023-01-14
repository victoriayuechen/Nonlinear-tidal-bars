using Plots
using WriteVTK
using HDF5
# using VTK
gr();

using HDF5, Revise
includet("src/read_output.jl");

# paths
array_path = "output/test/arrays/"
plots_path = mkpath("output/test/plots/") * "/"

# read output
x,y,xᵤ,xᵥ,yᵤ,yᵥ = read_grid(array_path); # read grid
t, ζ, u, v, c, = read_swec(array_path);  # read swec 
ζ_avg, u_avg, v_avg, c_avg, qu_avg, qv_avg = read_residual_swec(array_path); # read residuals
t_mor,h = read_bottomheight(array_path); # read bottom height


# show(ζ)
# show(ndims(t))
# show(size(t))
# show(t)
# show(typeof(file))

# file = h5open(array_path*"swec.h5")
# show(typeof(file))
# data = read(file, "zeta")
# plot(data)
# savefig("file.pvd")









# # plot residual currents
# step_y=2
# step_x=5
# x_sm = x[1:step_x:end];
# y_sm = y[1:step_y:end];
# u_avg_sm = u_avg[1:step_x:end,1:step_y:end];
# v_avg_sm = v_avg[1:step_x:end,1:step_y:end];
# u_avg_sm .= u_avg_sm ./ maximum(u_avg_sm) .* 5e2;
# v_avg_sm .= v_avg_sm ./ maximum(v_avg_sm) .* 5e1;
# X = repeat(x_sm,1,length(y_sm));
# Y = repeat(y_sm',length(x_sm),1);
# contour(x,y,h[:,:,1]',fill=true,clims=(-0.2,0.2))
# quiver!(vec(X),vec(Y), quiver = (vec(u_avg_sm),vec(v_avg_sm)), c=:black )


# # plot ζ,u,v,c versus time at one location.
# index_x = 10
# index_y = 10
# str_x = string(x[index_x]/1000)
# str_y = string(y[index_y]/1000)
# loc = "x = " * str_x * " km   " * "y = " * str_y * " km"
# pζ = plot(t,ζ[index_x,index_y,:]                 , ylabel="ζ",leg=false,c=:black,lw=2,title=loc)
# pu = plot(t,u[index_x,index_y,:]                 , ylabel="u",leg=false,ylim = (-0.8,0.8)  ,c=:black,lw=2)
# pv = plot(t,v[index_x,index_y,:]                 , ylabel="v",leg=false,ylim = (-0.15,0.15),c=:black,lw=2)
# pc = plot(t,c[index_x,index_y,:],xlabel = "t [h]", ylabel="c",leg=false,ylim = (0,7e-5)     ,c=:black,lw=2)
# plot(pζ,pu,pv,pc,layout = (4, 1),size=(600,800))


# # animations
# animᵤ = @animate for tn in 1:size(u,3)
#     heatmap(xᵤ,yᵤ,u[:,:,tn]',clims=(-1,1))
# end
# gif(animᵤ,  plots_path * "u.gif",    fps = 10)

# animᵥ = @animate for tn in 1:size(v,3)
#     contour(xᵥ,yᵥ,v[:,:,tn]',fill=true,clims=(-0.01,0.01))
# end
# gif(animᵥ,  plots_path * "v.gif",    fps = 10)

# anim_ζ = @animate for tn in 1:size(ζ,3)
#     contour(x,y,ζ[:,:,tn]',fill=true,clims=(-0.2,0.2))
# end
# gif(anim_ζ, plots_path * "zeta.gif", fps = 10)

# anim_c = @animate for tn in 1:5:size(c,3)
#     contour(x,y,c[:,:,tn]',fill=true,clims=(0,1e-3))
# end
# gif(anim_c, plots_path * "c.gif",    fps = 10)


