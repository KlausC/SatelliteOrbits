
#=
struct SatelliteModel{E<:Epoch,X,Y}
    origin::Body
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

    z = zero(SVector{3})
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
end

import NBodySimulator.get_accelerating_function

function get_accelerating_function(p::SolarSystemParameters, simulation::NBodySimulation)
    gpe = grav_param(earth)
    gpm = grav_param(moon)
    gps = grav_param(sun)
    gpj = grav_param(jupiter)
    spk = p.spk
    t0 = TDBEpoch("2020-12-23T21:00:00")
    function ssaccel(dv, u, v, t, i)
        #println("ssaccel(::$(typeof(dv)), ::$((u))")

        if i == 1
            dv .= 0
        else
            tj = t0 + t*seconds
            pm = SVector{3}(position(spk, tj, earth, moon))
            ps = SVector{3}(position(spk, tj, earth, sun))
            u = u[:,i]
            um = u .- pm
            us = u .- ps
            gg = -gpe / norm(u)^3
            gm = -gpm / norm(um)^3
            gs = -gps / norm(us)^3
            dve = accel2(spk, tj, 600)
            dp = gm * um + gs * us
            println("ssaccel gearth $(gg * u)")
            println("ssaccel gms    $dp")
            println("ssaccel earth  $dve")
            println("ssaccel delta  $(dp - dve)")
            dv .= gg * u + dp - dve
            dm = deviation(spk, tj, moon, gpm, u)
            ds = deviation(spk, tj, sun, gps, u)
            dj = deviation(spk, tj, jupiter_barycenter, gpj, u)
            println("ssaccel dmoon  $(dm)")
            println("ssaccel dsun   $(ds)")
            println("ssaccel djupi  $(dj)")
            println("ssaccel dall   $(dm + ds + dj)")
            dv .= dj + ds + dm + gg * u
        end
        nothing
    end
end

function deviation(spk, tj, body, gpb, u)
    ub = SVector{3}(position(spk, tj, earth, body))
    us = ub .- u
    gs = gpb / norm(us)^3 * us
    gb = gpb / norm(ub)^3 * ub
    gs - gb
end

function accel(spk, t0, dt)
    tm = t0 - dt*seconds
    tp = t0 + dt*seconds
    sp = velocity(spk, tp, earth)
    sm = velocity(spk, tm, earth)
    (sp - sm) ./ 2dt
end

function accel2(spk, t0, dt)
    g0 = accel(spk, t0, dt)
    g1 = accel(spk, t0, 2dt)
    (g0 - g1) / 3 + g0
end

toseconds(hrs) = seconds(hrs).Δt

function test2(;hrs=1hours, speed=1, angle=0)
    gpe = grav_param(earth)
    h0 = mean_radius(earth) + 400
    v0 = sqrt(gpe / h0)

    z = zero(SVector{3})
    s, c = sincosd(angle)
    v = v0 * speed
    r = SVector(s * h0, c * h0, 0.0)
    v = SVector(c*v, -s*v, 0.0)
    mass = 1.25
    body1 = MassBody(r, v, mass)
    bearth = MassBody(z, z, 0.0)

    parameters = SolarSystemParameters(spk)
    pdict = Dict(:custom_potential_params => parameters)

    system = PotentialNBodySystem([bearth, body1], pdict)
    tspan = (0.0, toseconds(hrs))
    simulation = NBodySimulation(system, tspan)

    sim_result = run_simulation(simulation)

    #animate(sim_result, "eart2.gif", fps=2)

    sim_result



end
