#=
 Incompressible Navier-Stokes 2D Flow Solver
 Author: Shun Arahata
=#
# using PyPlot
#FLOW CONDITIONS-----------------------------------
const RE = 70.0 # Reynolds Number
const CFL = 0.2 #  CFL Number

# SOR Pamameters
const OMEGAP = 1.00
const MAXITP = 100
const ERRORP = 0.0001

# No. of Time Steps
const NLAST = 5000 # Steps
const NLP = 10

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
        X[i, j] = DX * (i-ICENT)
        Y[i, j] = DY * (j-JCENT)
    end
end

@inbounds function slvflw()
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
    bcforp(p)
    bcforv(u, v)
    # time marching
    println("Step /Res(p) at itr. /CD/CL/Cp1/Cp2")
    println("   Step Res(p) CD CL Cp1 Cp2")
    nsteps = 0
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
        cl = 0.0
        for i in I1+1:I2
            cpbtm = (2p[i, J1] + 2p[i+1, J1])/2
            cptop = (2p[i, J2] + 2p[i+1, J2])/2
            cl += (cpbtm - cptop)*DX
        end
        # monitor by NLP steps
        if n % NLP == 0
            cp1 = 2p[I2 + I2 - I1, J1]
            cp2 = 2p[I2 + I2 - I1, J2]
            println("   ",nstep,",", resp,",", itrp,"," ,cd,"," ,cl,"," ,cp1,"," ,cp2)
        end
    end
    # write final results
    return v
end

# boundary condition for velocity pressure
@inbounds function bcforp(p::Array{Float64, 2})
    @simd for j in 1:MY
        # inflow condition i = 1
        p[1, j] = 0.0
        # dowmstream condition i = MX
        p[MX, j] = 0.0
    end
    @simd for i in 1:MX
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
        u[MX, j] = 2u[MX-1, j] - u[MX-2, j]
        v[MX, j] = 2v[MX-1, j] - v[MX-2, j]
    end
    for i in 1:MX
        # bottom condition j = 1
        # bottom condition j = MY
        u[i, 1] = 2u[i, 2] - u[i, 3]
        v[i, 1] = 2v[i, 2] - v[i, 3]
        u[i, MY] = 2u[i, MY - 1] - u[i, MY - 2]
        v[i, MY] = 2v[i, MY - 1] - v[i, MY - 2]
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
            break
        end
    end
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
        u[I1+1, j] = 2u[I1, j] - u[I1-1, j]
        u[I2-1, j] = 2u[I2, j] - u[I2+1, j]
        v[I1+1, j] = 2v[I1, j] - v[I1-1, j]
        v[I2-1, j] = 2v[I2, j] - v[I2+1, j]
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
        u[i, J1 + 1] = 2u[i, J1] - u[i, J1 - 1]
        u[i, J2 - 1] = 2u[i, J2] - u[i, J2 + 1]
        v[i, J1 + 1] = 2v[i, J1] - v[i, J1 - 1]
        v[i, J2 - 1] = 2v[i, J2] - v[i, J2 + 1]
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

function main()
    # make Array
    X = zeros(MX, MY)
    Y = zeros(MX, MY)
    # set grd
    setgrd(X,Y)
    # solve flow
    foo =  slvflw()
    #=
    fig = figure()
    ax = fig[:add_subplot](111)
    img = ax[:imshow](transpose(foo))
    PyPlot.plt[:show]()
    =#
end

@time main()
