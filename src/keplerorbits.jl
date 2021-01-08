
export kepler_ecliptic, kepler_icrf, eccanomaly, kepler_parameters

"""
 Planete positions derived from orbit elements and Kepler's calculations.
 Valid for times between 3000 BC and 3000 AD
 Accuracy limited.
 Data and procedures from https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
"""
const OrbitElements2a = [
#   a          e          I          L            ω₋          Ω
    0.38709843 0.20563661 7.00559432 252.25166724 77.45771895 48.33961819  0.00000000 0.00002123 -0.00590158 149472.67486623 0.15940013 -0.12214182
    0.72332102 0.00676399 3.39777545 181.97970850 131.76755713 76.67261496  -0.00000026 -0.00005107 0.00043494 58517.81560260 0.05679648 -0.27274174
    1.00000018 0.01673163 -0.00054346 100.46691572 102.93005885 -5.11260389  -0.00000003 -0.00003661 -0.01337178 35999.37306329 0.31795260 -0.24123856
    1.52371243 0.09336511 1.85181869 -4.56813164 -23.91744784 49.71320984  0.00000097 0.00009149 -0.00724757 19140.29934243 0.45223625 -0.26852431
    5.20248019 0.04853590 1.29861416 34.33479152 14.27495244 100.29282654  -0.00002864 0.00018026 -0.00322699 3034.90371757 0.18199196 0.13024619
    9.54149883 0.05550825 2.49424102 50.07571329 92.86136063 113.63998702  -0.00003065 -0.00032044 0.00451969 1222.11494724 0.54179478 -0.25015002
    19.18797948 0.04685740 0.77298127 314.20276625 172.43404441 73.96250215  -0.00020455 -0.00001550 -0.00180155 428.49512595 0.09266985 0.05739699
    30.06952752 0.00895439 1.77005520 304.22289287 46.68158724 131.78635853  0.00006447 0.00000818 0.00022400 218.46515314 0.01009938 -0.00606302
    39.48686035 0.24885238 17.14104260 238.96535011 224.09702598 110.30167986  0.00449751 0.00006016 0.00000501 145.18042903 -0.00968827 -0.00809981
]

const OrbitElements2b = [
#     b          c           s           f
    -0.00012452 0.06064060 -0.35635438 38.35125000
    0.00025899 -0.13434469 0.87320147 38.35125000
    0.00058331 -0.97731848 0.17689245 7.67025000
    -0.00041348 0.68346318 -0.10162547 7.67025000
    -0.01262724 0 0 0
]

const days_per_century = 36525.0
const julian_days_epoch2000 = 2451545.0
const obliquity2000 = 23.43928

"""
    a, E, e, I, ω, Ω = kepler(plant_no, time in days Julian Ephemerids Date)
    input data:
    1    a - major half axis (au)
    2    e - eccentricity (1)
    3    I - inclination (deg)
    4    L - mean longitude (deg)
    5    ω₀ - longitude of perihelon modified (deg)
    6    Ω - longitude of ascending node (deg)

    return values
    a major half axis
    e eccentricity
    E eccentric anomaly
    I inclination
    ω longitude of perihelion
    Ω longitude of ascending node
"""
function kepler_parameters(n::Integer, time::Float64)
    T = (time - julian_days_epoch2000) / days_per_century
    N = 6
    J = 4
    oelements = [OrbitElements2a[n,i+N] * T + OrbitElements2a[n,i] for i = 1:N]
    a, e, I, L, ωq, Ω = oelements
    ω = ωq - Ω
    M = L - ωq
    J = 4
    if n > J
        b, c, s, f = OrbitElements2b[n-J,:]
        sf, cf = sincosd(f * T) # check if f is in degrees/century
        M += T^2 * b + cf * c + sf * s
    end
    while M <= -180.0
        M += 360.0
    end
    while M > 180
        M -= 360.0
    end
    es = e * (180.0 / pi)
    E = eccanomaly(M, es)

    return a, e, E, I, ω, Ω
end

function kepler_ecliptic(n::Integer, time::Float64)
    a, e, E, I, ω, Ω = kepler_parameters(n, time)

    # 4. heliocentric coordinates in plnat's orbital plane
    s, c = sincosd(E)
    xs = (c - e) * a
    ys = sqrt(1 - e^2) * s * a

    # 5. coordinates in J2000 ecliptic plane - x-axis to equinox
    so, co = sincosd(ω)
    sO, cO = sincosd(Ω)
    sI, cI = sincosd(I)
    xecl = (co * cO - so * sO * cI) * xs + (-so * cO - co * sO * cI) * ys
    yecl = (co * sO + so * cO * cI) * xs + (-so * sO + co * cO * cI) * ys
    zecl = so * sI * xs + co * sI * ys

    return xecl, yecl, zecl
end

function kepler_icrf(n::Integer, time::Float64)
    xecl, yecl, zecl = kepler_ecliptic(n, time)

    # 6. coordinates in international celestial reference frame J2000
    se, ce = sincosd(obliquity2000)
    xeq = xecl
    yeq = ce * yecl - se * zecl
    zeq = se * yecl + ce * zecl
    
    return xeq, yeq, zeq
end

"""
    eccanomaly(M, es)
    eccentric anomaly `E` from mean anomaly `M` and eccentricity `es`

    `M = E - es * sin(E)`

    all arguments in degrees
"""
function eccanomaly(M, es)
    e = es * (pi / 180.0)
    dE = 1.0
    tol = 1e-6
    E = M + es * sind(M)
    while abs(dE) > tol
        sE, cE = sincosd(E)
        dE = (M - (E - es * sE)) / (1.0 - e * cE)
        E += dE
    end
    return E
end