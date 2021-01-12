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
