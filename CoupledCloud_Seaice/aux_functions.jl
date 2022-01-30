# ************* AUX FUNCTIONS *******************************
function nearest_profile(Xi::Vector{DateTime}, Vt::Vector{DateTime}, Vn::Matrix)

    idx_min = [argmin(abs.(Ti .- Vt)) for Ti ∈ Xi]

    return Vn[:, idx_min]
    
end


function get_variable_at_time(Xi::Vector{DateTime}, Vt::Vector{DateTime}, Vn::AbstractVector)
	
	# base value to convert to milliseconds:
	X0 = floor(Vt[1], Day)
	
	# check whether input Dates are from the same day:
	X0 != floor(Xi[1], Day) && error("Xi and Vt must be same day!")
	
	# Convert Input Dates to Milliseconds and extract the values as Int
	Xn = (Vt .- X0) .|> x->x.value
	X_vn = (Xi .- X0) .|> x->x.value
	
	# Create interpolation function:
	ipt = LinearInterpolation(Xn, Vn, extrapolation_bc=Line())
	
	# return interpolated value for every input DateTime:
	return ipt.(X_vn);
end

function classify_sic(ice)
    all(isnan.(ice)) && return 3
    let seaice = filter(!isnan, ice)
        all(seaice .≤ 15) && return 0
        any(15 .< seaice .< 95 ) && return 1
        all(seaice .≥ 95) && return 2
    end
end


function CleanData(V; minV = [1f-3])
    filter(!isnan, V) |> x->isempty(x) ? minV : x
end

function FunctionCleanData(V; minV = [1f-3], myfun=mean)
    filter(!isnan, V) |> x->isempty(x) ? minV : myfun(x)
end

# obtain the % of lwc and iwc within the profile in a layer between [Hmin:Hmax]
function ratio_IWC_LWC(IWC, LWC; Hmin=2f2, Hmax=8f3, T_in=nothing)
    pin = findall(Hmin .≤ IWC[:height] .≤ Hmax)

    tmpLWC = isnothing(T_in) ? LWC[:lwc] : nearest_profile(T_in, LWC[:time], LWC[:lwc])
    tmpIWC = isnothing(T_in) ? IWC[:iwc] : nearest_profile(T_in, IWC[:time], IWC[:iwc])
    
    Nlwc = [filter(!isnan, V[pin]) |> length for V ∈ eachcol(tmpLWC)]
    Nlwc *= 100/length(pin)

    Niwc = [filter(!isnan, V[pin]) |> length for V ∈ eachcol(tmpIWC)]
    Niwc *= 100/length(pin)

    # indexes where LWC is present within ice cloud:
    i_sclc = @. (Nlwc > 0) & (Niwc > 0)
    
    return Niwc, Nlwc, i_sclc

end
