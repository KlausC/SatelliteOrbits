
export kepler_ecliptic, kepler_icrf, eccanomaly, kepler_parameters, OrbitElements1, OrbitElements2
export kepler_ecliptic_d, kepler_icrf_d, eccanomaly_d, kepler_parameters_d

"""
 Planete positions derived from orbit elements and Kepler's calculations.
 Valid for times between 3000 BC and 3000 AD
 Accuracy limited.
 Data and procedures from https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
"""
const ORBIT1 = [
#    a          e          I            L           ω₋          Ω
    0.38709927 0.20563593 7.00497902 252.25032350 77.45779628 48.33076593  0.00000037 0.00001906 -0.00594749 149472.67411175 0.16047689 -0.12534081
    0.72333566 0.00677672 3.39467605 181.97909950 131.60246718 76.67984255  0.00000390 -0.00004107 -0.00078890 58517.81538729 0.00268329 -0.27769418
    1.00000261 0.01671123 -0.00001531 100.46457166 102.93768193 0.0 0.00000562 -0.00004392 -0.01294668 35999.37244981 0.32327364 0.0
    1.52371034 0.09339410 1.84969142 -4.55343205 -23.94362959 49.55953891  0.00001847 0.00007882 -0.00813131 19140.30268499 0.44441088 -0.29257343
    5.20288700 0.04838624 1.30439695 34.39644051 14.72847983 100.47390909  -0.00011607 -0.00013253 -0.00183714 3034.74612775 0.21252668 0.20469106
    9.53667594 0.05386179 2.48599187 49.95424423 92.59887831 113.66242448  -0.00125060 -0.00050991 0.00193609 1222.49362201 -0.41897216 -0.28867794
    19.18916464 0.04725744 0.77263783 313.23810451 170.95427630 74.01692503  -0.00196176 -0.00004397 -0.00242939 428.48202785 0.40805281 0.04240589
    30.06992276 0.00859048 1.77004347 -55.12002969 44.96476227 131.78422574  0.00026291 0.00005105 0.00035372 218.45945325 -0.32241464 -0.00508664
    39.48211675 0.24882730 17.14001206 238.92903833 224.06891629 110.30393684  -0.00031596 0.00005170 0.00004818 145.20780515 -0.04062942 -0.01183482
]
const ORBIT2A = [
#    a          e          I            L           ω₋          Ω
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

const ORBIT2B = [
#     b          c           s           f
    -0.00012452  0.06064060 -0.35635438 38.35125000
     0.00025899 -0.13434469  0.87320147 38.35125000
     0.00058331 -0.97731848  0.17689245 7.67025000
    -0.00041348  0.68346318 -0.10162547 7.67025000
    -0.01262724  0.0         0.0        0.0
]

abstract type OrbitElements end
struct OrbitElements1 <: OrbitElements
    t1
    OrbitElements1() = new(ORBIT1)
end
struct OrbitElements2 <: OrbitElements
    t1
    t2
    OrbitElements2() = new(ORBIT2A, ORBIT2B)
end

const days_per_century = 36525.0 # days / 100 years
const julian_days_epoch2000 = 2451545.0 # days
const obliquity2000 = 23.43928 # deg

const pi180 = pi / 180.0 # 1/deg
const AU = 149597870.700 # km
const seconds_per_day = 86400 # s/day
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
function kepler_parameters_d(n::Integer, time::Real, oe::OrbitElements)
    oe1 = oe.t1
    oe2 = oe isa OrbitElements2 ? oe.t2 : zeros(0, 0)
    T = (time - julian_days_epoch2000) / days_per_century
    T_d = inv(days_per_century)
    k = 6
    J = size(oe1, 1) - size(oe2, 1)
    oelements_d = oe1[n,k+1:2k] * T_d
    oelements = [oe1[n,i+k] * T + oe1[n,i] for i = 1:k]
    a_d, e_d, I_d, L_d, ωq_d, Ω_d = oelements_d
    a, e, I, L, ωq, Ω = oelements
    a *= AU
    a_d *= AU
    ω = ωq - Ω
    ω_d = ωq_d - Ω_d
    M = L - ωq
    M_d = L_d   - ωq_d
    if oe isa OrbitElements2 && n > J
        b, c, s, f = oe2[n-J,:]
        fpi8 = pi180 * f * T_d
        sf, cf = sincosd(f * T) # assuming f is in degrees/century
        sf_d, cf_d = fpi8 * cf, -fpi8 * sf
        M += T^2 * b + cf * c + sf * s
        M_d += 2T * T_d * b + cf_d * c + sf_d * s 
    end
    M = mod(M + 180.0, 360.0) - 180.0
    es = e / pi180
    E, E_d = eccanomaly_d((M, M_d), (e, e_d))

    return [a, e, E, I, ω, Ω], [a_d, e_d, E_d, I_d, ω_d, Ω_d]
end
function kepler_parameters(n::Integer, time::Real, oe::OrbitElements)
oe1 = oe.t1
oe2 = oe isa OrbitElements2 ? oe.t2 : zeros(0, 0)
T = (time - julian_days_epoch2000) / days_per_century
k = 6
J = size(oe1, 1) - size(oe2, 1)
oelements = [oe1[n,i+k] * T + oe1[n,i] for i = 1:k]
a, e, I, L, ωq, Ω = oelements
a *= AU
ω = ωq - Ω
M = L - ωq
if oe isa OrbitElements2 && n > J
    b, c, s, f = oe2[n-J,:]
    sf, cf = sincosd(f * T) # assuming f is in degrees/century
    M += T^2 * b + cf * c + sf * s
end
M = mod(M + 180.0, 360.0) - 180.0
es = e / pi180
E = eccanomaly(M, e)

return [a, e, E, I, ω, Ω]
end

function kepler_ecliptic_d(n::Integer, time::Real, oe::OrbitElements)
    (a, e, E, I, ω, Ω), (a_d, e_d, E_d, I_d, ω_d, Ω_d) = kepler_parameters_d(n, time, oe)

    # 4. heliocentric coordinates in planet's orbital plane
    s, c = sincosd(E); hE = pi180 * E_d
    s_d, c_d = c * hE, -s * hE
    xs = (c - e) * a
    xs_d = (c - e) * a_d + (c_d - e_d) * a
    sq = sqrt(1 - e^2)
    sq_d = - e * e_d / sq
    ys = sq * s * a
    ys_d = sq_d * s * a + sq * s_d * a + sq * s * a_d

    # 5. coordinates in J2000 ecliptic plane - x-axis to equinox
    so, co = sincosd(ω); hω = pi180 * ω_d
    so_d, co_d = co * hω , -so * hω
    sO, cO = sincosd(Ω); hΩ = pi180 * Ω_d
    sO_d, cO_d = cO * hΩ, -sO * hΩ
    sI, cI = sincosd(I); hI = pi180 * I_d
    sI_d, cI_d = cI * hI, -sI * hI
    sosO = so * sO; sosO_d = so_d * sO + so * sO_d
    cocO = co * cO; cocO_d = co_d * cO + co * cO_d 
    socO = so * cO; socO_d = so_d * cO + so * cO_d
    cosO = co * sO; cosO_d = co_d * sO + co * sO_d
    sosOcI = sosO * cI; sosOcI_d = sosO_d * cI + sosO * cI_d
    socOcI = socO * cI; socOcI_d = socO_d * cI + socO * cI_d
    cosOcI = cosO * cI; cosOcI_d = cosO_d * cI + cosO * cI_d
    cocOcI = cocO * cI; cocOcI_d = cocO_d * cI + cocO * cI_d
    sosI = so * sI; sosI_d = so_d * sI + so * sI_d
    cosI = co * sI; cosI_d = co_d * sI + co * sI_d
    mxx =  cocO - sosOcI; mxx_d =  cocO_d - sosOcI_d
    mxy = -socO - cosOcI; mxy_d = -socO_d - cosOcI_d
    myx =  cosO + socOcI; myx_d =  cosO_d + socOcI_d
    myy = -sosO + cocOcI; myy_d = -sosO_d + cocOcI_d

    xecl = mxx * xs + mxy * ys
    yecl = myx * xs + myy * ys
    zecl = sosI * xs + cosI * ys

    xecl_d = mxx_d * xs + mxx * xs_d + mxy_d * ys + mxy * ys_d
    yecl_d = myx_d * xs + myx * xs_d + myy_d * ys + myy * ys_d
    zecl_d = sosI_d * xs + sosI * xs_d + cosI_d * ys + cosI * ys_d

    return [xecl, yecl, zecl], [xecl_d, yecl_d, zecl_d]
end

function kepler_ecliptic(n::Integer, time::Real, oe::OrbitElements)
    a, e, E, I, ω, Ω = kepler_parameters(n, time, oe)

    # 4. heliocentric coordinates in planet's orbital plane
    s, c = sincosd(E)
    xs = (c - e) * a
    sq = sqrt(1 - e^2)
    ys = sq * s * a

    # 5. coordinates in J2000 ecliptic plane - x-axis to equinox
    so, co = sincosd(ω)
    sO, cO = sincosd(Ω)
    sI, cI = sincosd(I)
    sosO = so * sO
    cocO = co * cO
    socO = so * cO
    cosO = co * sO
    sosOcI = sosO * cI
    socOcI = socO * cI
    cosOcI = cosO * cI
    cocOcI = cocO * cI
    sosI = so * sI
    cosI = co * sI
    mxx =  cocO - sosOcI
    mxy = -socO - cosOcI
    myx =  cosO + socOcI
    myy = -sosO + cocOcI

    xecl = mxx * xs + mxy * ys
    yecl = myx * xs + myy * ys
    zecl = sosI * xs + cosI * ys

    return [xecl, yecl, zecl]
end

function kepler_icrf_d(n::Integer, time::Real, oe::OrbitElements)
    (xecl, yecl, zecl), (xecl_d, yecl_d, zecl_d) = kepler_ecliptic_d(n, time, oe)

    # 6. coordinates in international celestial reference frame J2000
    se, ce = sincosd(obliquity2000)
    xeq = xecl; xeq_d = xecl_d
    yeq = ce * yecl - se * zecl; yeq_d = ce * yecl_d - se * zecl_d
    zeq = se * yecl + ce * zecl; zeq_d = se * yecl_d + ce * zecl_d
    
    return [xeq, yeq, zeq], [xeq_d, yeq_d, zeq_d]
end

function kepler_icrf(n::Integer, time::Real, oe::OrbitElements)
    xecl, yecl, zecl= kepler_ecliptic(n, time, oe)

    # 6. coordinates in international celestial reference frame J2000
    se, ce = sincosd(obliquity2000)
    xeq = xecl
    yeq = ce * yecl - se * zecl
    zeq = se * yecl + ce * zecl
    
    return SVector{3}(xeq, yeq, zeq)
end

"""
    eccanomaly(M, e)
    eccentric anomaly `E` from mean anomaly `M` and eccentricity `e`

    `M = E - e * 180/π * sin(E)`

    M and E in degrees, e dimensionless
"""
function eccanomaly_d((M, M_d), (e, e_d))
    es = e / pi180
    dE = 1.0
    tol = 1e-6
    E = M + es * sind(M)
    while abs(dE) > tol
        sE, cE = sincosd(E)
        dE = (M - (E - es * sE)) / (1.0 - e * cE)
        E += dE
    end
    sE, cE = sincosd(E)
    E_d = (M_d + e_d * (180.0 / pi) * sE) / (1.0 - e * cE)
    return E, E_d
end
function eccanomaly(M, e)
    es = e / pi180
    tol = 1e-12
    dE = 1.0
    E = M + es * sind(M)
    while abs(dE) > tol
        sE, cE = sincosd(E)
        dE = (M - (E - es * sE)) / (1.0 - e * cE)
        E += dE
    end
    return E
end
