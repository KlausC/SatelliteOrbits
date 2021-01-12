
import SatelliteOrbits.AU

tj = 2451545.0 # julian date J2000 
@testset "kepler parameters $p, $dt, 1" for dt = (0, 100), p = 1:9
    (a, e, E, I, ω, Ω) = kepler_parameters(p, tj + dt, OrbitElements1())
    a0, e0, I0, L0, ω0, Ω0 = OrbitElements1().t1[p,1:6]
    a1, e1, I1, L1, ω1, Ω1 = OrbitElements1().t1[p,7:12]

    a0 *= AU
    a1 *= AU
    dt = dt / 36525
    
    a0 += a1 * dt
    e0 += e1 * dt
    I0 += I1 * dt
    L0 += L1 * dt
    ω0 += ω1 * dt
    Ω0 += Ω1 * dt

    @test a ≈ a0
    @test e ≈ e0
    @test I ≈ I0
    @test Ω ≈ Ω0
    @test ω ≈ ω0 - Ω0
    es = 180 / pi * e
    M = E - sind(E) * es
    M0 = L0 - ω0
    M0 = mod(M0 + 180, 360) - 180
    @test M ≈ M0
end

@testset "kepler parameters $p, $dt, 2" for dt = (0, 100), p = 1:9
    (a, e, E, I, ω, Ω) = kepler_parameters(p, tj + dt, OrbitElements2())
    a0, e0, I0, L0, ω0, Ω0 = OrbitElements2().t1[p,1:6]
    a1, e1, I1, L1, ω1, Ω1 = OrbitElements2().t1[p,7:12]

    a0 *= AU
    a1 *= AU
    dt = dt / 36525
    
    a0 += a1 * dt
    e0 += e1 * dt
    I0 += I1 * dt
    L0 += L1 * dt
    ω0 += ω1 * dt
    Ω0 += Ω1 * dt

    @test a ≈ a0
    @test e ≈ e0
    @test I ≈ I0
    @test Ω ≈ Ω0
    @test ω ≈ ω0 - Ω0
    es = 180 / pi * e
    M = E - sind(E) * es
    M0 = L0 - ω0
    if p > 4
        b, c, s, f = OrbitElements2().t2[p-4, 1:4]
        dM = b * dt^2 + c * cosd(f * dt) + s * sind(f * dt)
        M0 += dM
    end
    M0 = mod(M0 + 180, 360) - 180
    @test M ≈ M0
end
