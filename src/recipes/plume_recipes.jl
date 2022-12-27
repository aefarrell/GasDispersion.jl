@recipe function f(p::Plume; nx=100, ny=100, height=2)

    xlims --> (-1,1000)
    ylims --> (-50,50)
    seriestype := :contour
    seriescolor --> :thermal
    fill --> true
    xlabel --> "Downwind distance, m"
    ylabel --> "Crosswind distance, m"

    # define the x part of the grid
    xmin, xmax = plotattributes[:xlims]
    xs = range(xmin, xmax; length=nx)

    # define the y part of the grid
    ymin, ymax = plotattributes[:ylims]
    ys = range(ymin, ymax;length=ny)

    xs,ys,(x,y)->p(x,y,height)

end
