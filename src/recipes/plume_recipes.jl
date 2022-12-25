@recipe function f(p::Plume; nx::Integer=100, ny::Integer=100, height::Number=2)

    xlims --> (-1,1000)
    ylims --> (-50,50)
    seriestype := :contour
    seriescolor --> :thermal
    fill --> true

    # define the x part of the grid
    xmin, xmax = plotattributes[:xlims]
    xs = range(start=xmin,stop=xmax,length=nx)

    # define the y part of the grid
    ymin, ymax = plotattributes[:ylims]
    ys = range(start=ymin,stop=ymax,length=ny)

    xs,ys,(x,y)->p(x,y,height)

end
