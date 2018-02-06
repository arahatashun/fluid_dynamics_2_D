#=
 Incompressible Navier-Stokes 2D Flow Solver
 Author: Shun Arahata
=#
# using PlotlyJS
# using Plotly
using PyPlot
# using Plots
#FLOW CONDITIONS-----------------------------------
print("Input an Reynolds Number: ")
const RE =  parse(Float64, readline())# Reynolds Number
const CFL = 0.2 #  CFL Number

# SOR Pamameters
const OMEGAP = 1.00
const MAXITP = 10000
const ERRORP = 0.0001

# No. of Time Steps
const NLAST = 6000 # Steps
const NLP = 100

# set x-grid parameters
const MX = 401 # grid number of x
const I1 = 96 # left edge
const I2 = 106 # right edge
const DX = 1.0 / (I2 - I1)
const ICENT = (I1 + I2) / 2

# set y-grid parameters
const MY = 201 # grid number of y
const J1 = 96 # top edge
const J2 = 106 # bottom edge
const DY = 1.0 / (J2 - J1)
const JCENT = (J1 + J2) / 2

# set time step size
const DT = CFL * min(DX, DY)
#-----------------------------------------------------
# set grid
@inbounds function setgrd(X,Y)
    for i in 1:MX, j in 1:MY
        X[i] = DX * (i-ICENT)
        Y[j] = DY * (j-JCENT)
    end
end

function slvflw()
    # print conditions
    println("***  Comp. conditions,")
    println("       CFL = ", CFL)
    println("       DT  = ", DT)
    println("       ",NLAST,"Time Steps to go ...")
    println(" ")
    println(">> 2D Incompressible Flow Solver")
    println("   Re",RE)
    println("   No. of Grid Points", MX, MY)
    println("   CFL:",CFL,", DT:", DT,", Steps:",NLAST)
    # initial condition (uniform flow)
    nbegin  = 0
    time = 0.0
    u = ones(MX, MY)
    v = zeros(u)
    p = zeros(u)
    X = zeros(MX)
    Y = zeros(MY)
    setgrd(X, Y)
    bcforp(p)
    bcforv(u, v)
    # time marching
    println("Step /Res(p) at itr. /CD/CL/Cp1/Cp2")
    println("   Step Res(p) CD CL Cp1 Cp2")
    nsteps = 0
    cd_list = zeros(NLAST, 1)
    cl_list = zeros(cd_list)
    cp1_list = zeros(cd_list)
    cp2_list = zeros(cd_list)
    #anim = @animate
    for n in 1:NLAST
        nstep = n + nbegin
        time += DT
        # solve poisson for p
        resp, itrp = poiseq(p, u, v)
        bcforp(p)
        # update u,v
        veloeq(p, u, v)
        bcforv(u, v)

        # calculate CL & CD
        cd = 0.0
        for j in J1+1:J2
            cpfore = (2p[I1, j] + 2p[I1, j+1])/2
            cpback = (2p[I2, j] + 2p[I2, j+1])/2
            cd += (cpfore - cpback)*DY
        end
        cd_list[n] = cd
        cl = 0.0
        for i in I1+1:I2
            cpbtm = (2p[i, J1] + 2p[i+1, J1])/2
            cptop = (2p[i, J2] + 2p[i+1, J2])/2
            cl += (cpbtm - cptop)*DX
        end
        cl_list[n] = cl
        # heatmap(X, Y, transpose(v),clim = (-1, 1), aspect_ratio=1, c=:viridis, title = "$(time) sec" )
        # monitor by NLP steps
        cp1 = 2p[I2 + I2 - I1, J1]
        cp2 = 2p[I2 + I2 - I1, J2]
        if n % NLP == 0
            println("   ",nstep,",", resp,",", itrp,"," ,cd,"," ,cl,"," ,cp1,"," ,cp2)
        end
        cp1_list[n] = cp1
        cp2_list[n] = cp2
    end # every 100
    # gif(anim, "karman_Re=$(RE)gif", fps = 10)
    # write final results
    # return p, u, v
    return cd_list, cl_list, cp1_list, cp2_list
end

# boundary condition for velocity pressure
@inbounds function bcforp(p::Array{Float64, 2})
    for j in 1:MY
        # inflow condition i = 1
        p[1, j] = 0.0
        # dowmstream condition i = MX
        p[MX, j] = 0.0
    end
    for i in 1:MX
        # bottom condition j = 1
        p[i, 1] = 0.0
        # bottom condition j = MY
        p[i, MY] = 0.0
    end
    # wall condition
    # define four point
    p[I1, J1] = p[I1 - 1, J1 - 1]
    p[I1, J2] = p[I1 - 1, J2 + 1]
    p[I2, J1] = p[I2 + 1, J1 - 1]
    p[I2, J2] = p[I2 + 1, J2 + 1]
    # four sides
    for j in J1+1:J2-1
        p[I1, j] = p[I1 - 1, j]
        p[I2, j] = p[I2 + 1, j]
    end
    for i in I1+1:I2-1
        p[i, J1] = p[i, J1 - 1]
        p[i, J2] = p[i, J2 + 1]
    end
end


# boundary condition for velocity
@inbounds function bcforv(u::Array{Float64, 2}, v::Array{Float64, 2})
    for j in 1:MY
        # inflow condition i =1
        u[1, j] = 1.0
        v[1, j] = 0.0
        # downstream condition i = MX
        u[MX, j] = 2.0 * u[MX-1, j] - u[MX-2, j]
        v[MX, j] = 2.0 * v[MX-1, j] - v[MX-2, j]
    end
    for i in 1:MX
        # bottom condition j = 1
        # bottom condition j = MY
        u[i, 1] = 2.0 * u[i, 2] - u[i, 3]
        v[i, 1] = 2.0 * v[i, 2] - v[i, 3]
        u[i, MY] = 2.0 * u[i, MY - 1] - u[i, MY - 2]
        v[i, MY] = 2.0 * v[i, MY - 1] - v[i, MY - 2]
    end
    # wall condition
    for i in I1:I2, j in J1:J2
        u[i, j] = 0.0
        v[i, j] = 0.0
    end
end


# poisson equation
@inbounds function poiseq(p::Array{Float64, 2}, u::Array{Float64, 2}, v::Array{Float64, 2})
    # compute RHS
    rhs = zeros(p)
    res = 0.0
    itrp = 0
    for i in 2:MX-1, j in 2:MY-1
        if i<I1 || i>I2 || j<J1 || j>J2
            ux = (u[i + 1, j] - u[i - 1, j]) / (2DX)
            uy = (u[i, j + 1] - u[i, j - 1]) / (2DY)
            vx = (v[i + 1, j] - v[i - 1, j]) / (2DX)
            vy = (v[i, j + 1] - v[i, j - 1]) / (2DY)
            rhs[i, j] = (ux + vy) / DT - (ux^2 + 2uy * vx + vy^2)
        end
    end
    #iterations
    relaxation_flag = true
    for itr in 0:MAXITP
        res = 0.0
        # relaxation
        for i in 2:MX-1, j in 2:MY-1
            if i<I1 || i>I2 || j<J1 || j>J2
                dp = (p[i + 1, j] + p[i - 1, j]) / (DX^2) + (p[i, j + 1] + p[i, j - 1]) / (DY^2) -rhs[i, j]
                dp = dp / (2.0 / (DX^2) + 2.0 /(DY^2)) - p[i, j]
                res += dp^2
                p[i, j] += OMEGAP * dp
            end
        end
        # set boundary condition
        bcforp(p)
        res = sqrt(res/(MX*MY))
        itrp = itr
        if res < ERRORP
            relaxation_flag = false
            break
        end
    end
    relaxation_flag&&println("relaxation error")
    return res, itrp
end

# Kawamura scheme
@inbounds function veloeq(p::Array{Float64, 2}, u::Array{Float64, 2}, v::Array{Float64, 2})
    urhs = zeros(p)
    vrhs = zeros(p)
    # pressure gradient
    for i in 2:MX-1, j in 2:MY-1
        if i<I1 || i>I2 || j<J1 || j>J2
            urhs[i, j] = - (p[i + 1, j] - p[i - 1, j])/(2DX)
            vrhs[i ,j] = - (p[i, j + 1] - p[i, j - 1])/(2DY)
        end
    end
    # viscous term
    for i in 2:MX-1, j in 2:MY-1
        if i<I1 || i>I2 || j<J1 || j>J2
            urhs[i, j] +=
            (u[i + 1,j] - 2u[i, j] + u[i - 1, j]) / (RE * DX^2)
            + (u[i, j + 1] - 2u[i, j] + u[i, j - 1]) / (RE * DY^2)
            vrhs[i, j] +=
            (v[i + 1, j] - 2v[i, j] + v[i - 1, j]) / (RE * DX^2)
            + (v[i, j + 1] - 2v[i, j] + v[i, j - 1]) / (RE * DY^2)
        end
    end
    # advection term in x direction
    for j in J1+1:J2-1
        u[I1+1, j] = 2.0 * u[I1, j] - u[I1-1, j]
        u[I2-1, j] = 2.0 * u[I2, j] - u[I2+1, j]
        v[I1+1, j] = 2.0 * v[I1, j] - v[I1-1, j]
        v[I2-1, j] = 2.0 * v[I2, j] - v[I2+1, j]
    end
    for i in 3:MX-2, j in 3:MY-2
        if i<I1 || i>I2 || j<J1 || j>J2
            urhs[i, j] -= u[i, j] * (-u[i+2, j] + 8.0 * (u[i+1, j] - u[i-1, j]) + u[i-2, j]) / (12.0*DX) +
            abs(u[i, j]) * (u[i+2, j] - 4.0 * u[i+1, j] + 6.0 * u[i, j] - 4.0 * u[i-1, j] + u[i-2, j]) / (4.0*DX)
            vrhs[i, j] -= u[i, j] * (-v[i+2, j] + 8.0 * (v[i+1, j] - v[i-1, j]) + v[i-2, j]) / (12.0*DX) +
            abs(u[i, j]) * (v[i+2, j] - 4.0 * v[i+1, j] + 6.0 * v[i, j] - 4.0 * v[i-1, j] + v[i-2, j]) / (4.0*DX)
        end
    end
    # advection term in y direction
    for i in I1+2:I2-1
        u[i, J1 + 1] = 2.0 * u[i, J1] - u[i, J1 - 1]
        u[i, J2 - 1] = 2.0 * u[i, J2] - u[i, J2 + 1]
        v[i, J1 + 1] = 2.0 * v[i, J1] - v[i, J1 - 1]
        v[i, J2 - 1] = 2.0 * v[i, J2] - v[i, J2 + 1]
    end
    for i in 3:MX-2,j in 3:MY-2
        if i<I1 || i>I2 || j<J1 || j>J2
            urhs[i, j] -= v[i, j] * (-u[i, j+2] + 8.0 * (u[i, j+1] - u[i, j-1]) + u[i, j-2]) / (12.0*DY) +
            abs(v[i, j]) * (u[i, j+2] - 4.0 * u[i, j+1] + 6.0 * u[i, j] - 4.0 * u[i, j-1] + u[i, j-2]) / (4.0*DY)
            vrhs[i, j] -= v[i, j] * (-v[i, j+2] + 8.0 * (v[i, j+1] - v[i, j-1]) + v[i, j-2]) / (12.0*DY) +
            abs(v[i, j]) * (v[i, j+2] - 4.0 * v[i, j+1] + 6.0 * v[i, j] - 4.0 * v[i, j-1] + v[i, j-2]) / (4.0*DY)
        end
    end
    #update
    for i in 2:MX-1, j in 2:MY-1
        if i<I1 || i>I2 || j<J1 || j>J2
            u[i, j] += DT * urhs[i, j]
            v[i, j] += DT * vrhs[i, j]
        end
    end
end


function plot_matplot_all(p, u, v)
    X = zeros(MX)
    Y = zeros(MY)
    setgrd(X, Y)
    fig = figure()
    ax = fig[:add_subplot](311)
    cp = ax[:contour](X, Y, transpose(p), 10)
    ax[:clabel](cp, inline=1, fontsize=10)
    ylabel("Y")
    title("Contour Plot of Pressure")
    ax = fig[:add_subplot](312)
    cp = ax[:contour](X, Y, transpose(u), 10)
    ax[:clabel](cp, inline=1, fontsize=10)
    ylabel("Y")
    title("Contour Plot of U")
    ax = fig[:add_subplot](313)
    cp = ax[:contour](X, Y, transpose(v), 10)
    ax[:clabel](cp, inline=1, fontsize=10)
    xlabel("X")
    ylabel("Y")
    title("Contour Plot of V")
    PyPlot.plt[:savefig]("contour_plot.png",dpi=300)
end

function plot_byplotly(p, u, v)

    X = zeros(MX)
    Y = zeros(MY)
    setgrd(X, Y)
    # trace = contour(x=Y, y=X, z=p, colorscale="Jet", contours_coloring="lines")
    trace = surface(x=Y, y=X, z=p, colorscale="Jet", contours_coloring="lines")
    layout = Layout(title="P in a Contour Plot")
    my_plot = PlotlyJS.Plot(trace,layout)
    remote_plot = post(my_plot)
end


function plot_matplot(p, u, v)
    X = zeros(MX)
    Y = zeros(MY)
    setgrd(X, Y)
    fig = figure()
    ax = fig[:add_subplot](111)
    cp = ax[:contour](X, Y, transpose(p), 10)
    ax[:clabel](cp, inline=1, fontsize=10)
    ylabel("Y")
    xlabel("X")
    title("Contour Plot of Pressure $(DT*NLAST) sec")
    PyPlot.plt[:savefig]("contour_plot.pgf")
end

function plot_coefficient(cd, cl, cp1, cp2)
    fig = figure()
    ax = fig[:add_subplot](111)
    start_index = Int(NLAST*3/4)
    time_list = linspace(start_index*DT,DT * NLAST, NLAST - start_index+1)
    ax[:plot](time_list, cd[start_index:NLAST], label="CD")
    ax[:plot](time_list, cl[start_index:NLAST], label="CL")
    ax[:plot](time_list, cp1[start_index:NLAST], label="Cp1")
    ax[:plot](time_list, cp2[start_index:NLAST], label="Cp2")
    xlabel("time")
    cd_average = sum(cd)/length(cd)
    title("Re = $(RE), CD = $(cd_average) ")
    legend(loc = 1)
    PyPlot.plt[:savefig]("$(trunc(Int,RE)).pgf")
    PyPlot.plt[:show]()
    # plot_fft(cl[start_index:NLAST])
end

function plot_fft(signal)
    am = fft(signal)
    fig = figure()
    ax = fig[:add_subplot](111)
    plotnum = div(1.0,DT)
    freq = linspace(0, 1, plotnum)
    am = abs(am)[1:plotnum]
    ax[:plot](freq, am)
    main_freq = freq[indmax(am)]
    xlabel("Hz")
    ylabel("Amptitude")
    title("Re = $(RE), $(main_freq)Hz")
    PyPlot.plt[:savefig]("fft_$(trunc(Int,RE)).pgf")
    PyPlot.plt[:show]()
end

function main()
    # solve flow
    # p, u, v =  slvflw()
    # plot_matplot(p, u, v)
    # plot_byplotly(p, u, v)
    cd, cl, cp1, cp2  = slvflw()
    plot_coefficient(cd, cl, cp1, cp2)
end

@time main()
