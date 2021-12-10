#!/opt/julia-1.6.0/bin/julia
###!/home/psgarfias/.local/julia-1.6.1/bin/julia

# script to add the wind direction and SIC to the database based on Lapse-rate

# loading required packages:
using NCDatasets
using Navigation
using Dates
using Printf
using ARMtools
using JLD2
using Statistics
##using GMT
using ATMOStools
const ATM=ATMOStools

# including SEAICE functions:
include(joinpath(homedir(), "LIM/repos/SEAICEtools.jl/src/SEAICEtools.jl"))

# ****
# defining PATHS:
const CAMPAIGN = "arctic-mosaic";
# variables definition for SeaIce:
const BASE_PATH = joinpath(homedir(), "LIM/data/B07/")
const DATA_PATH = joinpath(homedir(), "LIM/remsens", CAMPAIGN);
###const DATA_PATH = "/media/psgarfias/LaCie SSD/LIM/data/";
const SENSOR = "modis_amsr2"; #"amsr2";
const PRODUCT = "SeaIce";
const DATA_SOURCE = "NAV";
const PROC_PATH = joinpath(BASE_PATH, PRODUCT);
const LATLON_FILE = joinpath(PROC_PATH, SENSOR, "lonlat_nsidc_1km.nc");
# Lapse-rate files:
const LPR_PATH = joinpath(BASE_PATH, "LP")


#const LATLON_FILE = joinpath(PROC_PATH, SENSOR, "LongitudeLatitudeGrid-n3125-ChukchiBeaufort.h5");

# variables definition for Radiosondes:
#const BASE_PATH = joinpath("/projekt2/remsens/data_new/site-campaign/", CAMPAIGN);


# *********************************************
### Defining central point e.g. NSA coordinates

nsa_coor = Point(71.323e0, -156.609e0)

# Defining study region parameter:
R_lim = 50f0
θ₀ = 1f0
θ₁ = 360f0

years = (2019)
months =(10, 11,12)
days =(1:31)

!isempty(ARGS) && foreach(ARGS) do argin
	ex = Meta.parse(argin)
	eval(ex)
end

# Reading Sea Ice product lat, lon grid (it doesn't change with day)
sic_coor = SEAICE.read_LatLon_Bremen_product(LATLON_FILE);

# ********** starting loop over time *********************
for yy ∈ years
    # * Reading RV Polarstern NAV data:
    RV = load("data/RVpolarstern_track_$(yy).jld2", "RV");
    ### rvplt = plot(RV[:lon], RV[:lat], lw=1, lc=:blue, label="drift");
    
    for mm ∈ months
        for dd ∈ days

            # selecting RV track from the day:
            tt = try
                DateTime(yy,mm,dd)
            catch
                continue
            end
            println("Working on $dd.$mm.$yy...")

            iday = findall(tt .≤ RV[:time] .< (tt + Day(1))
            Nnav = length(iday)

            # * Loading Thermodinamic information:
            lpr_file = ARMtools.getFilePattern(LPR_PATH, CAMPAIGN, yy, mm, dd)
            isnothing(lpr_file) && (println("no $(lpr_file)"); continue)

            flpr = jldopen(lpr_file, "a+")
            INV = flpr["INVvar"]
            CBH = flpr["CBH"]
            decop_hgt = flpr["decop_hgt"]

    
            ### Obtain wind direction profile based on water vapour flux:
            rs_file = ARMtools.getFilePattern(DATA_PATH, "INTERPOLATEDSONDE", yy ,mm, dd)
            rs = ARMtools.getSondeData(rs_file);
            H_wvt, i_wvt, i_wv50  = ATM.estimate_WVT_peak_altitude(rs) #
            #PBLH = ATM.estimate_Ri_PBLH(rs)

            # * Finding out Wind_dir from porfile and thermodynamic variables:
            wind_dir, wind_spd, wind_range = ATM.Collect_WindDir_using_WVT(INV, rs, CBH,
                                                                           decop_hgt, H_wvt);

            ## READING SIC data for AMSR2 data for $(dd). $(mm). $(yy)
            Pstern = Point(RV[:lat][iday[2]], RV[:lon][iday[2]])
            R_lim = 50e3;
            LonLim, LatLim = SEAICE.estimate_box(Pstern, R_lim, δR=5e3)
            #LonLim = extrema(RV[:lon][iday]) .+ (-3, 3) #(Pstern.λ-3, Pstern.λ+3)
            #LatLim = extrema(RV[:lat][iday]) .+ (-1, 1) #(Pstern.ϕ-1, Pstern.ϕ+1)
	
            idx_box = SEAICE.extract_LonLat_Box(LonLim, LatLim, sic_coor);
            sic_filen = ARMtools.getFilePattern(PROC_PATH, SENSOR, yy, mm, dd);
            isnothing(sic_filen) && (println("No $(sic_filen)..."); continue)
    
            # Reading corresponding SIC data for the selected sector with idx indexes:
            SIC = SEAICE.read_SIC_Bremen_product(sic_filen, idx_box, SICPROD="mersic");

            lon_box, lat_box = SEAICE.Get_LonLat_From_Point(sic_coor[idx_box]);

            # **** PLOTS:

            ### dayplt = scatter(lon_box, lat_box, marker_z=SIC, color=:ice, markerstrokewidth=0,
            ##                 markersize=2, marker=:square, size=(500,500),
	    ##                 clim=(60, 100), label=false, xlim=LonLim, ylim=LatLim,
            # #                title="$(RV[:time][iday[1]])", colorbartitle="SIC %");
            ### plot!(dayplt, rvplt)            
                
            # defining temporal Dictionaries to storage part of the day
            tmpSIC = Dict(:SIC2d => Matrix{Float64}(undef, 10, length(H_wvt)),
                          :errSIC2d => Matrix{Float64}(undef, 10, length(H_wvt)),
                          :range => [])

            tmpWIND = Dict(:WD => wind_dir, :WS => wind_spd,
                           :Idx_W => Vector{Any}(undef, length(H_wvt)) );
            for itime=2:Nnav
                # new position for RV polarstern:
                Pstern = Point(RV[:lat][iday[itime]], RV[:lon][iday[itime]])
                println(Pstern)
                
                # Plotting
                lat = map(p->p.ϕ, sic_coor[idx_box]);
	        lon = map(p->p.λ, sic_coor[idx_box]);
                P_circ = SEAICE.Create_Semi_Circle(Pstern, 1f0, 360f0, R_lim=50e3)
                lon_circ, lat_circ = SEAICE.Get_LonLat_From_Point(P_circ)
                
                ##plot!(dayplt, [RV[:lon][iday[itime]]], [RV[:lat][iday[itime]]], linestyle=:auto,
                ##      marker=:star, markersize=10, mc=RGB(Nnav/itime, .1 ,.1), label="RV Polarstern");

                ##plot!(dayplt, lon_circ, lat_circ, lw=1, lc=RGB(Nnav/itime, .1 ,.1), alpha=0.2, label="radii")
                

                # polar coordinates for the new RV polarstern position:
                θ_all, ρ_all = SEAICE.LonLat_To_CenteredPolar(Pstern, sic_coor);

                idx_nav = let RVtime = (RV[:time][iday[itime-1]], RV[:time][iday[itime]])
                    findall(RVtime[1] .< rs[:time] .≤ RVtime[2])
                end
                # calculating the SIC as 2D matrix from wind_dir and SIC box:
                idx_radial, ρ2d, meSIC2d, sdSIC2d = SEAICE.wind_idx_radial(wind_dir[idx_nav],
                                                                           wind_range[idx_nav],
                                                                           θ_all[idx_box],
                                                                           ρ_all[idx_box],
                                                                           SIC);

                tmpSIC[:SIC2d][:, idx_nav] = meSIC2d;
                tmpSIC[:errSIC2d][:, idx_nav] = sdSIC2d;
                tmpSIC[:range] = ρ2d;
                
                tmpWIND[:Idx_W][idx_nav] = idx_radial;

            end  # over loop of Nnav (number of RV positions along the day)
            # * Saving results to LPR files:
            !haskey(flpr, "WVT") && (flpr["WVT"] = Dict(:HWVT => H_wvt,
                                                        :IdxWVT => i_wvt,
                                                        :Idx50qv => i_wv50));
            
            !haskey(flpr, "SIC") && (flpr["SIC"] = tmpSIC);
            
            !haskey(flpr, "WIND") && (flpr["WIND"] = tmpWIND);
            
            close(flpr)

            ##savefig(dayplt, "plots/rv_sic_$(Dates.format(RV[:time][iday][1], "yyyy-mm-dd")).png")

        end  # end over days
    end  # end overs months
end  # end overs years
# ---/




# end of script
