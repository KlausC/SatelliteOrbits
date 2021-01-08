
#=
struct SatelliteModel{E<:Epoch,X,Y}
    origin::CelestialBody
    time::E
    pos::X
    vel::V
end

struct GravitationData
    gp
end

function acceleration(state::AbstractVector, m::SatelliteModel, gd::GravitationData)
    pos = state[1:3]
    gp = gd.gp
    - gp / norm(pos)^3 * pos
end

function gen_rhs(m::SatelliteModel, gd::GravitationalData)
    function rhs!(deriv::AbstractVector, state::AbstractVector)
        deriv[1:3] = state[4:6]
        deriv[4:6] = acceleration(state, m, gd)
    end
end
=#

const G = AstroBase.gravitational_const()[1] # unit km^3/kg/s^2
const spkpath = abspath(joinpath(homedir(), "dev", "JPLEphemeris", "deps"))
const spk = SPK(joinpath(spkpath, "de440.bsp"))

function test1()
    gpe = grav_param(earth)
    h0 = mean_radius(earth) + 400
    g0 = gpe / h0^2
    v0 = sqrt(g0 * h0)
    me = gpe / G

    z = zeros(SVector{3})
    r = SVector(h0, 0.0, 0.0)
    v = SVector(0.0, v0, 0.0)
    mass = 1.25
    body1 = MassBody(r, v, mass)
    terra = MassBody(z, z, me)

    system = GravitationalSystem([terra, body1], G)

    tspan = (0.0, 2*3600.0)
    simulation = NBodySimulation(system, tspan)

    sim_result = run_simulation(simulation)

    animate(sim_result, "eart1.gif", fps=2)

    sim_result
end

struct SolarSystemParameters <: PotentialParameters
    spk
    t0
    u0
    v0
end

function SolarSystemParameters(spk, t0)
    st = state(spk, t0, earth)
    SolarSystemParameters(spk, t0, st[1], st[2])
end

import NBodySimulator.get_accelerating_function

function get_accelerating_function(p::SolarSystemParameters, simulation::NBodySimulation)
    gpe = grav_param(earth)
    gpm = grav_param(moon)
    gps = grav_param(sun)
    gpj = grav_param(jupiter)
    spk = p.spk
    t0 = p.t0
    u0 = p.u0
    v0 = p.v0
    r0 = mean_radius(earth) + 300
    tt = 600

    function ssaccel(dv, u, v, t, i)
        #println("ssaccel(::$(typeof(dv)), ::$((u))")

        if i == 1
            dv .= 0
        else
            tj = t0 + t*seconds
            u = u[:,i]
            v = v[:,i]
            pe, qe = SVector{3}.(state(spk, tj, earth, earth))
            pm = SVector{3}(position(spk, tj, earth_barycenter, moon))
            ps = SVector{3}(position(spk, tj, earth_barycenter, sun))
            ue = u  .- pe
            um = ue .- pm
            us = ue .- ps
            ve = v  .- qe
            println("time=$tj $ue $u")
            r = norm(ue)
            s = norm(v)
            gg = -gpe / r^3
            gm = -gpm / norm(um)^3
            gs = -gps / norm(us)^3
            dp = gm * um + gs * us
            #println("ssaccel gearth     $(gg * u)")
            #println("ssaccel g moon/sun $dp")
            #=
            if r <= r0
                dp .-= ve * (r0/r - 1) / tt
            end
            =#
            dv .= gg * ue + 0 * dp
        end
        nothing
    end
end

function test2(;hrs=1hours, speed=1, angle=0)
    t0 = TDBEpoch("2020-12-23T21:00:00")
    p = SolarSystemParameters(spk, t0)
    u0, v0 = p.u0, p.v0

    gpe = grav_param(earth)
    h0 = mean_radius(earth) + 400
    vx = sqrt(gpe / h0)

    z = zero(SVector{3})
    s, c = sincosd(angle)
    v = vx * speed
    re, ve = state(spk, t0, earth, earth)
    r = SVector{3}(re + [s * h0, c * h0, 0.0])
    v = SVector{3}(ve + [c*v, - s*v, 0.0])
    mass = 1.25
    body1 = MassBody(r, v, mass)
    bearth = MassBody(z, z, 0.0)

    pdict = Dict(:custom_potential_params => p)

    system = PotentialNBodySystem([bearth, body1], pdict)
    tspan = (0.0, toseconds(hrs))
    simulation = NBodySimulation(system, tspan)

    sim_result = run_simulation(simulation, rtol=1e-15, alg_hints=[:stiff])

    #animate(sim_result, "eart2.gif", fps=2)

    sim_result
end

function deviation(spk, tj, body, gpb, u)
    ub = (position(spk, tj, earth, body))
    us = ub .- u
    gs = gpb / norm(us)^3 * us
    gb = gpb / norm(ub)^3 * ub
    gs - gb
end

accel(spk, t0, dt=600) = accel(spk, t0, ssb, earth, dt)

function accel1(spk, t0, orig, body, dt)
    tm = t0 - dt*seconds
    tp = t0 + dt*seconds
    sp = velocity(spk, tp, orig, body)
    sm = velocity(spk, tm, orig, body)
    (sp - sm) ./ 2dt
end

function accel(spk, t0, orig, body, dt=600)
    g0 = accel1(spk, t0, orig, body, dt)
    g1 = accel1(spk, t0, orig, body, 2dt)
    (g0 - g1) / 3 + g0
end

toseconds(hrs) = seconds(hrs).Δt

struct SolarSystemP
    m::Float64
    spk::SPK
    t0::TDBEpoch
    inertial::CelestialBody # inertialsystem center (e.g. sun or earth_barycenter)
    origin::CelestialBody # origin of coordinates (e.g. earth)
    bodies
    gravis
    SolarSystemP(spk::SPK, t0::TDBEpoch, inertial::CelestialBody, origin::CelestialBody, bodies) =
        new(1.0, spk, t0, inertial, origin, bodies, grav_param.(bodies))
end

#=
function L(x, v, param, t)
    Q = q[1:3]
    P = p[1:3]
    tj = param.t0 + t * seconds
    pos = x - param.position.(Ref(spk), Ref(tj), Ref(param.origin), param.bodies)
    param.m / 2 * (v + velocity(spk, tj, param.inertial, param.origin))^2 + param.m * sum(param.gravis ./ pos)
end

p(v, t) = ∂L/∂v = m * (v + velocity(t, origin))  => v(p, t) = p / m - velocity(t, origin)

H(p, q) = sum(v(p, t) * p) - L(q, v(p, t), t) = p^2 / m - p * velocity(origin, t) - p^2 /2m + m * sum(gravis(body) / (q + pos(orig) - pos(body)))

=#

function H(p, q, param)
    t = q[4]    # artifical coordinate to formulate H independent of time argument
    E = p[4]    # corresponding energy term which keeps H constant
    Q = q[1:3]  # coordinates relative to param.origin
    P = p[1:3]  # linear momentum relative to param.inertial
    tj = param.t0 + t * seconds
    #distances between satellite and center of bodies
    pos= [Q .- position(spk, tj, param.origin, param.bodies[i]) for i = 1:length(param.bodies)]
    # kinetec energy
    K = norm(P)^2 / 2param.m - dot(P, velocity(spk, tj, param.inertial, param.origin))
    # gravitational potential
    V = - param.m * sum(param.gravis ./ norm.(pos))
    K + V + E
end

# ∂H/∂p (velocity relative to param.origin)
function dq(qdot, p, q, param::SolarSystemP, t)
    t = q[4]
    P = p[1:3]
    tj = param.t0 + t * seconds
    qdot .= ([P / param.m - velocity(spk, tj, param.inertial, param.origin); 1.0])
end

# -∂H/∂q (force relative to param.initial)
function dp(pdot, p, q, param::SolarSystemP, t)
    t = q[4]
    E = p[4]
    Q = q[1:3]
    P = p[1:3]
    tj = param.t0 + t * seconds
    n = length(param.bodies)

    pos = [Q .- position(spk, tj, param.origin, param.bodies[i]) for i = 1:n]
    dp13 = -param.m * sum( param.gravis[i] / norm(pos[i])^3 * pos[i] for i in 1:n)
    dp4 = dot(P, accel(spk, tj, param.inertial, param.origin))
    dp4 += param.m * sum(param.gravis[i] / norm(pos[i])^3 * dot(pos[i], accel(spk, tj, param.origin, param.bodies[i])) for i = 1:n)
    pdot .= ([dp13; dp4])
end

function startvalues(h, v, dir, rot)
    pos = normalize(dir) * h
    vel = normalize(cross(rot, dir)) * v
    vel, pos
end

function startvalues(ssp::SolarSystemP, t::Period, tv::Period, h, vfactor)
    spk = ssp.spk
    tj = ssp.t0 + t
    ts = toseconds(t)
    b = ssp.bodies
    if length(ssp.bodies) >= 2
        body1 = b[1]
        body2 = length(b) >= 2 ? b[2] : body1 != ssp.origin ? ssp.origin : ssp.inertial
        dir = position(spk, tj + tv, body1, body2)
        rot = cross(velocity(spk, tj + tv, body1, body2), position(spk, tj + tv, body1, body2))
    else
        rot = [0.0, 0.0, 1.0]
        dir = [1.0, 0.0, 0.0]
    end
    v = sqrt(ssp.gravis[1] / norm(h)) * vfactor
    vel, pos = startvalues(h, v, dir, rot)
    vel .+= velocity(spk, tj, ssp.inertial, body1)
    pos .+= position(spk, tj, ssp.origin, body1)
    E = 0.0 # H([pos; ts], [vel; 0.0], ssp)
    ([vel * ssp.m; E]), ([pos; ts])
end

function problem(ssp::SolarSystemP, tspan, tv, h, v)
    ts = toseconds(tspan[1]), toseconds(tspan[2])
    p0, q0 = startvalues(ssp, tspan[1], tv, h, v)
    HamiltonianProblem((dp, dq), p0, q0, ts, ssp)
end

function cb_moonapproach()
    ContinuousCallback(moonposvel, on_approach)
end

function moonstate(u, t, param)
    spk = param.spk
    tj = param.t0 + t*seconds
    V = u[1:3] / param.m
    R = u[5:7]
    pos = R .- position(spk, tj, param.origin, moon)
    vel = V .- velocity(spk, tj, param.origin, moon)
    vel, pos
end

function moonposvel(u, t, integrator)
    V, R = moonstate(u, t, integrator.p)
    dot(V, R)
end

function on_approach(integrator)
    u, t, param  = integrator.u, integrator.t, integrator.p
    V, R = moonstate(u, t, param)
    println("closest to moon: $(norm(R)) km t=$t s R=$R km")
end

function plotorbit(t0::TDBEpoch, tv::Period, n::Real, v, h=6731.0; dt = 300, id=2)
    ssp = SolarSystemP(spk, t0, earth_barycenter, earth, [earth, moon])
    prob = problem(ssp, (0days, n*days), tv, -h, v)
    sol = solve(prob, KahanLi8(), dt=dt, callback=cb_moonapproach())
    x = [sol(i*3600*12)[5] for i = 0:max(2n,1)]
    y = [sol(i*3600*12)[6] for i = 0:max(2n,1)]
    xaxis, yaxis = axis(ssp, tv, id)
    println(xaxis, " ", yaxis)
    p = plot!(sol, vars=(5,6), ratio=1.0)
    scatter!(p, x, y; xaxis, yaxis, zcolor = (1:length(x)).*100)
    display(p)
    sol
end

function axis(ssp, tv, i=2)
    spk = ssp.spk
    t0 = ssp.t0
    tj = t0 + tv
    if 1 <= i <= length(ssp.bodies)
        r = JPLEphemeris.mean_radius(ssp.bodies[i])
        pos = position(spk, tj, ssp.origin, ssp.bodies[i])
    else
        r = 4e5
        pos = zeros(3)
    end
    [pos[1] - r, pos[1] + r] .+ [5000, 5000]
    [pos[2] - r, pos[2] + r] .+ [5000, 5000]
end

function plotmoon(t0, n)
    n = n * 12
    p = plot(legend=false, aspect_ratio=1.0)
    pos = [position(spk, t0 + (i/12)*days, earth, moon) for i = 0:n]
    x, y = getindex.(pos, 1), getindex.(pos, 2)
    plot!(p, x, y, xaxis=[-4e5,0], yaxis=[-4e5,1e5])
    scatter!(p, x[1:6:end], y[1:6:end], zcolor = 1:6:length(x))
    display(p)
end

const t0 = TDBEpoch("2021-01-03T21:00:00")

radius(sol, t) = begin st = sol(t); hypot(st[5], st[6], st[7]) end

function eccentricity(t; t0=t0, v=1.0)
    sol=plotorbit(t0, t, 0.075, v, dt=60)
    a = /(extrema(radius.(Ref(sol), sol.t))...)
    (1 - a) / (1 + a)
end

function ap(t; kwargs...)
    ec = eccentricity(t; kwargs...)
    2*ec / (1 + ec) * 6731
end

plotap(t0=t0; kwargs...) = plot([ap(i/2*days; kwargs...) for i = 0:56], legend=nothing)

function relpos(sol, t, ori=earth)
    spk = sol.prob.p.spk
    orig = sol.prob.p.origin
    tj = sol.prob.p.t0 + t*seconds
    sol(t)[5:7] - position(spk, tj, orig, ori)
end
function relvel(sol, t, ori=earth)
    spk = sol.prob.p.spk
    orig = sol.prob.p.inertial
    m = sol.prob.p.m
    tj = sol.prob.p.t0 + t*seconds
    sol(t)[1:3] / m - velocity(spk, tj, orig, ori)
end


const interesting_parameters = [
    (4.5days, -1.40090),
    (4.5days, +1.40250),
    (4.0days, +1.40220), # earth, moon
    (4.0days, +1.40230), # with sun
    (4.5days, +1.40175),
]
