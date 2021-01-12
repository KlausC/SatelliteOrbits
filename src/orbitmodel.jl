const G = AstroBase.gravitational_const()[1] # unit km^3/kg/s^2
const spkpath = abspath(joinpath(homedir(), "dev", "JPLEphemeris", "deps"))
const spk = SPK(joinpath(spkpath, "de440.bsp"))

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
    h = norm(pos[1])
    dp13 = -param.m * sum( param.gravis[i] / norm(pos[i])^3 * pos[i] for i in 1:n)
    dp4 = dot(P, accel(spk, tj, param.inertial, param.origin))
    dp4 += param.m * sum(param.gravis[i] / norm(pos[i])^3 * dot(pos[i], accel(spk, tj, param.origin, param.bodies[i])) for i = 1:n)
    if h < 10000.0
        v = P / param.m - velocity(spk, tj, param.inertial, param.bodies[1])
        b = brake(h, norm(v))
        println("brake = $b h = $h")
        dp13 -= v / norm(v) * b
    end
    pdot .= ([dp13; dp4])
end

"""
    brake(h, v)
    h in km distance to center of earth
    v in km/s relative to earth
    result in km/s^2
"""
function brake(h, v)
    c = 1e-8 * exp(-(h - 7000.0)/100)
    c * norm(v)^2
end

function startvalues(h, v, dir, rot)
    pos = normalize(dir) * h
    vel = normalize(cross(rot, dir)) * v
    vel, pos
end

function speed_v1(b::CelestialBody, h)
    G = grav_param(b)
    sqrt(G / norm(h))
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
    v = speed_v1(ssp.bodies[1], h) * vfactor
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
    vmoon = velocity(spk, tj, param.origin, moon)
    pos = R .- position(spk, tj, param.origin, moon)
    vel = V .- vmoon
    vel, pos, vmoon
end

function moonposvel(u, t, integrator)
    V, R = moonstate(u, t, integrator.p)
    dot(V, R)
end

function on_approach(integrator)
    u, t, param  = integrator.u, integrator.t, integrator.p
    V, R, vmoon = moonstate(u, t, param)
    if norm(R) < 5000.0
        println("closest to moon: $(norm(R)) km t=$t s R=$R km")
        v0 = norm(V)
        v1 = speed_v1(moon, R)
        vnew = V .* v1 / v0 .+ vmoon
        integrator.u = ArrayPartition([vnew; 0.0], u.x[2])
        println("reduced rel. speed from $v0 to $v1")
    end
end

function cb_closeearth(h)
    DiscreteCallback(earth_nearer(h), at_earth_close)
end
function earth_nearer(h)
    function (u, t, integrator)
        ssp = integrator.p
        spk = ssp.spk
        tj = integrator.p.t0 + t * seconds
        v = u[5:7] .- position(spk, tj, ssp.origin, earth)
        norm(v) < h
    end
end
function at_earth_close(integrator)
    println("crash! after $(integrator.t) s")
    terminate!(integrator)
end

function plotorbit(t0::TDBEpoch, tv::Period, n::Real, v, h=6731.0; dt = 300, id=1)
    ssp = SolarSystemP(spk, t0, ssb, earth, [earth, moon, sun])
    prob = problem(ssp, (0days, n*days), tv, -h, v)

    cb = CallbackSet(cb_moonapproach(), cb_closeearth(6331.0))
    sol = solve(prob, KahanLi8(), dt=dt, callback=cb)
    x = [sol(i*3600*12)[5] for i = 0:max(2n,1)]
    y = [sol(i*3600*12)[6] for i = 0:max(2n,1)]
    xaxis, yaxis = axis(ssp, tv, id)
    # println(xaxis, " ", yaxis)
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
    R = 4000.0
    ([pos[1] - r, pos[1] + r] .+ [-R, R], [pos[2] - r, pos[2] + r] .+ [-R, R])
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
