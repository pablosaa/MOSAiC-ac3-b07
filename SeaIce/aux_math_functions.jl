# *************************************************************
# ***** MATH FUNCTIONS 
"""
faktor(x,y) computes the ratio of x/y
"""
faktor(x, y) = x/y

"""
To compute uncertainties from function f = x*y or f=x/y 
δf(x, y, δx, δy, f) with x ± δx, y ± δy and f::Function
"""
δf(x, y, δx, δy, f) = f(x,y)*sqrt( (δx/x)^2 + (δy/y)^2 )

"""
Convert input to SIC scaled to reference SIC given by factor fₖ,
truncating the output to 100%
"""
fix(x, fₖ) = min(fₖ*x, 100)


"""
Transform lat, lon data to coordinates at a polar steregraphic system

julia> x, y = latlon2polarstereographic(ϕ, λ, ϕc, λ₀)

Where:
* ϕ : latitude in degree North
* λ : longitude in degree East
* x : map coordiante along X-axis [m]
* y : map coordinate along Y-axis [m]
Oprional input parameters:
* ϕ_c: latitude of true scale in degree North, standard parallel (default is -70)
* λ₀: meridian in degree East along positive Y-axis of map (default is 0)
* 𝑎 : Earth radius defined in the projection (default 6378137.0 m, WGS84)
* 𝑒 : eccentricity (default 0.08181919)

Polar stereographic projection used for satellite polar sea ice studies.
Equations adapted from "Map Projections - A Working Manual! by J.P. Snyder (1987)
"""
function latlon2xy(ϕ::T, λ::T; ϕc::T=70.0, λ₀::T=0.0, 𝑎=6378137.0, 𝑒=0.08181919) where T<:Real
    ϕ *= π/180
    ϕc*= π/180
    λ *= π/180
    λ₀*= π/180
    
    t(α) = tan(π/4-α/2)/√((1-𝑒*sin(α))/(1+𝑒*sin(α)))^𝑒  
    m(α) = cos(α)/√(1-(𝑒*sin(α))^2)
    ρ = 𝑎*m(ϕc)*t(ϕ)/t(ϕc)
    x = ρ*sin(λ-λ₀)
    y = -ρ*cos(λ-λ₀)

    return x, y
end

function latlon2xy(ϕ::Vector{T}, λ::Vector{T}; ϕc::T=70.0, λ₀::T=0.0, 𝑎=6378137.0, 𝑒=0.08181919) where T<:Real

    x = T[]
    y = T[]
    for (lat, lon) ∈ zip(ϕ, λ)
        xi, yi = latlon2xy(lat, lon; ϕc=ϕc, λ₀=λ₀, 𝑎=𝑎, 𝑒=𝑒)
        push!(x, xi)
        push!(y, yi)
    end
    return x, y
end
# ----/

"""
Transform x, y data from polar steregraphic system to latitude and longitude

julia> x, y = xy2latlon(x, y, ϕc, λ₀)

Where:
* x : X-axis coordinate of map, scalar or vector [m]
* y : Y-axis coordinate of map, scalar or vector [m]

Optional input parameters:
* ϕ_c: latitude of true scale in degree North, standard parallel (default is -70)
* λ₀: meridian in degree East along positive Y-axis of map (default is 0)
* 𝑎 : Earth radius defined in the projection (default 6378137.0 m, WGS84)
* 𝑒 : eccentricity (default 0.08181919)
OUTPUT:
* ϕ : latitude in degree North
* λ : longitude in degree East

Polar stereographic projection used for satellite polar sea ice studies.
Equations adapted from "Map Projections - A Working Manual! by J.P. Snyder (1987)

"""
function xy2latlon(x::T, y::T; ϕc::T=70.0, λ₀::T=0.0, 𝑎=6378137.0, 𝑒=0.08181919) where T<:Real

    ϕc *= π/180
    λ₀ *= π/180
    
    t(α) = tan(π/4-α/2)/√((1-𝑒*sin(α))/(1+𝑒*sin(α)))^𝑒          # Eq. (15-19)  
    m(α) = cos(α)/√(1-(𝑒*sin(α))^2)    # Eq. (14-15)
    tₖ = t(ϕc)
    mₖ = m(ϕc)
    
    ρ = √(x^2 + y^2)       # Eq. (20-18)
    
    t = ρ*tₖ/mₖ/𝑎          # Eq. (21-40)
    χ = π/2 - 2atan(t)     # Eq. (7-13)

    # Eq. (3-5)
    ϕ = χ + (𝑒^2/2 + 5𝑒^4/24 + 𝑒^6/12 + 13𝑒^8/360)*sin(2χ) +
        (7𝑒^4/48 + 29𝑒^6/240 + 811𝑒^8/11520)*sin(4χ) +
        (7𝑒^6/120 + 81𝑒^8/1120)*sin(6χ) +
        (4279𝑒^8/161280)*sin(8χ)
    
    λ = λ₀ + atan(x, -y)        # Eq. (20-16)
    
    # adapting longitude ∈ [-π, π]
    λ = mod(λ+π, 2π) - π
    
    return rad2deg(ϕ), rad2deg(λ)
end

function xy2latlon(x::Array{T}, y::Array{T}; ϕc::T=70.0, λ₀::T=0.0, 𝑎=6378137.0, 𝑒=0.08181919) where T<:Real

    
    ϕ = similar(x)
    λ = similar(y)

    tmp = @. xy2latlon(x, y, ϕc=ϕc, λ₀=λ₀, 𝑎=𝑎, 𝑒=𝑒)
    foreach(pairs(tmp)) do (i, V)
        ϕ[i] = V[1]
        λ[i] = V[2]
    end
    
    #for (xi, yi) ∈ zip(x, y)
    #    lat, lon = xy2latlon(xi, yi, ϕc=ϕc, λ₀=λ₀, 𝑎=𝑎, 𝑒=𝑒)
    #    push!(ϕ, lat)
    #    push!(λ, lon)
    #end
    
    return ϕ, λ
end
# ----/

# end of file
