using Gridap

x_hep_test = []
y_hep_test = []
i_grid_test = 0
j_grid_test = 12.5
while i_grid_test <= 10000
    append!(x_hep_test, i_grid_test)
    global i_grid_test += 50
end
while j_grid_test <= 975
    append!(y_hep_test, j_grid_test)
    global j_grid_test += 12.5
end
probe = [Point(i, j) for i in x_hep_test, j in y_hep_test]

show(x_hep_test)
show(y_hep_test)