#!/opt/julia-1.8.0/bin/julia

using NCDatasets
using Navigation
using Dates
using Distributions
using JLD2
using GMT
using Printf
using ARMtools
using CSV, DataFrames
using StatsBase

include("aux_math_functions.jl");

SEAICE = include(joinpath(homedir(), "LIM/repos/SEAICEtools.jl/src/SEAICEtools.jl"));

const PROD_PATH = joinpath(homedir(), "LIM/data/B07/SeaIce");
const RVNAV_PATH = joinpath(homedir(), "LIM/data/B07/arctic-mosaic");
const LATLON_FILE = joinpath(PROD_PATH, "modis_amsr2", "lonlat_nsidc_1km.nc");
const R_lim = 50e3;   # radius around RV polarstern
PRODUCTS = (:DIV, :LF, :SIC)

dist_wdir = Dict(data=>Dict() for data ∈ PRODUCTS)

for data ∈ PRODUCTS
#let data = :SIC
    # relevant constants specific for data product:
    cbfaktor, cblims, cbbins, cbcolor, cbtrunc, pxsize = if data==:DIV
        1f5, (-1.13, 1.13, 0.15), 0.3, :cork, (-1, 1), 0.07 
    elseif data==:LF
        1, (0, 0.4, 0.025), 0.1, :davos, (0, 1), 0.07
    elseif data==:SIC
        1, (75, 100, 2.5), 5, :vik, (-1, 0),  0.15
    else
        @error "given data key $(data) not supported"
    end

    data == :SIC && (DFsic=CSV.read("data/modis_amsr2_osisaf_statsII.csv", DataFrame))
    
    iplt = 0
    
#for (mm, yy) ∈ [(10,2019), (11,2019),(12,2019), (1,2020), (2,2020), (3,2020), (4,2020), (5,2020)]
mm = 11; yy=2019;

    #for dd = (1:31)
dd = 18
        heute = try
            Date(yy, mm, dd);
        catch
            @error "continue"
        end
    # reading RV Polarstern data:
    # @sprintf("../RV_polarstern/data/RVpolarstern_track_%04d.jld2", yy);
    RV_PATH = ARMtools.getFilePattern(RVNAV_PATH, "RVpolarstern", yy, mm, dd);
    isnothing(RV_PATH) && @error "no RV NAV file found for $(dd).$(mm).$(yy)"
    RV = load(RV_PATH, "RV")

## Reading wind direction of the day
winddir = CSV.read("/tmp/winddir_18112019.csv", DataFrame);


        # finding out RV positions on working date:
        iday = findall(==(heute), Date.(RV[:time]))
        isempty(iday) && (@warn "No RV nav for $(heute)";) #continue)
        Nnav = length(iday)

        # loading daily data
    sar = Dict()
    dist_wdir[data] = let tmp=range(cblims[1], stop=cblims[2], step=cblims[3])
        ndat = length(winddir.wvtdir)
        Dict(:bins=>cbfaktor*tmp,
             :freq=>fill(NaN32, length(tmp)-1, ndat),
             :μ => fill(NaN32, ndat),
             :σ => fill(NaN32, ndat),
             :qq=> fill(NaN32, 5, ndat)
             )
    end
    
        dat_coor = if data == :SIC
            lr_filen = ARMtools.getFilePattern(PROD_PATH, "modis_amsr2", yy, mm, dd)
            isnothing(lr_filen) && (@warn "No data for $(heute)") #; continue)
            SEAICE.read_LatLon_Bremen_product(LATLON_FILE)
        else
            
            lr_filen = SEAICE.load.FilePattern(PROD_PATH, "sentinel-1A_luisa/drift_corrected", yy, mm, dd);
            isnothing(lr_filen) && (@warn "No data for $(heute)";)# continue)

            sar_x, sar_y = SEAICE.load.LatLon_AWI_product(lr_filen)
            XX, YY = couple2grid(sar_x, sar_y);
            
            ##sar = SEAICE.load.Data_Divergence_LeadFraction(lr_filen);
            ##XX, YY = couple2grid(sar[:x], sar[:y]);
            lat, lon = xy2latlon(XX, YY, ϕₖ=70.0, λ₀=-45.0);
            Point.(lat, lon) #dat_coor = 
        end
        
        lon, lat = SEAICE.Get_LonLat_From_Point(dat_coor)
        
        println("Working on $(heute) with file $(lr_filen)")

    for i_rv ∈ iday

        
        rv_lat, rv_lon = RV[:lat][i_rv], RV[:lon][i_rv]
        Pstern = Point(rv_lat, rv_lon)
            
        # finding limits of Box:
        LonLim, LatLim = SEAICE.estimate_box(Pstern, R_lim, δR=10e3)

        idx_box = SEAICE.extract_LonLat_Box(LonLim, LatLim, dat_coor);

        if data==:SIC
            fφ = let datum=RV[:time][i_rv]
                df = getSICfixfactor(datum, DFsic)
                df[1, :ratio]
            end
            sar[:SIC] = fix.(SEAICE.read_SIC_Bremen_product(lr_filen,
                                                            idx_box,
                                                            SICPROD="mersic"), fφ)
        else
            sar = SEAICE.load.Data_Divergence_LeadFraction(lr_filen, idx_box)
        end
        # converting to polar for data within the box
        θ, ρ = SEAICE.LonLat_To_CenteredPolar(Pstern, dat_coor);
        
        # Tim = findall(RV[:time][i_rv] .≤ winddir[!, :date] .< RV[:time][i_rv+1])
        # Tim = [argmin(abs.(T.-RV[:time][i_rv])) for T ∈ winddir[!,:date]]
        #for iti ∈ Tim
        let iti=argmin(abs.(winddir[!,:date].-RV[:time][i_rv]))

            # getting wind direction where to extract:
            θᵢ = winddir[iti, :wvtdir]-3; θₑ = winddir[iti, :wvtdir]+3;

            
            idx_wd = findall((θᵢ .≤ θ[idx_box] .≤ θₑ) .& (ρ[idx_box] .≤ 50));

            # calculating statistics of box:
            LF = cbfaktor*sar[data]

            μLF, σLF, qqLF = if data==:LF
                stats_𝑁ₗᵤ(filter(≥(0), LF[idx_wd]))
            elseif data==:DIV
                stats_𝑁ₗᵤ(filter(!isnan, LF[idx_wd]), L=nothing, U=nothing)
            elseif data==:SIC
                
                stats_𝑁ₗᵤ(filter(!isnan, LF[idx_wd]), L=0, U=100)
                # make corrections for SIC here!
            else
                @error "data $(data) not supported!"
            end

            # Creating Histogrmas from the wind sector:
            dist_wdir[data][:freq][:, iti] = let Hdat=fit(Histogram,
                                                          LF[idx_wd],
                                                          dist_wdir[data][:bins],
                                                          closed=:right)
                Hdat.weights
            end
            dist_wdir[data][:μ][iti] = μLF
            dist_wdir[data][:σ][iti] = σLF
            dist_wdir[data][:qq][:, iti] = qqLF
            
            # creating lines with angle limits to plot
            ##Plon, Plat = SEAICE.create_pair_lines(Point(rv_lat, rv_lon), 50e3, [θᵢ, θₑ])

            # creating circle to plot
            P_circ = SEAICE.Create_Semi_Circle(Pstern, 1f0, 360f0, R_lim=R_lim)
            lon_circ, lat_circ = SEAICE.Get_LonLat_From_Point(P_circ)

    # creating wind lines to plot
            lon_wind, lat_wind = SEAICE.create_pair_lines(Pstern, R_lim, winddir[iti, :wvtdir].+[-3, 3])

            # creating title for plot
            plot_name = if data==:SIC
                "/tmp/quicklooks/$(data)/"*replace(basename(lr_filen),
                                              "sic"=>@sprintf("mosaic%04d", iplt),
                                              ".nc"=>"$(data).png");
                                              
            else
                                              
                "/tmp/quicklooks/$(data)/"*replace(basename(lr_filen),
                                              "pairs"=>@sprintf("mosaic%04d", iplt),
                                              "deformation_moved.nc"=>"$(data).png");
            end
            
            utc_title_str = @sprintf("%s %3.2f +/- %3.2f on %s", data, μLF, σLF,
                                     Dates.format(winddir[iti, :date], "dd-uuu HH:MM"))
            
            wind_title_str = "@:9:@%14%"*utc_title_str*"@%%@::"

            cbar = if data==:SIC
                GMT.makecpt(cmap=cbcolor, truncate=cbtrunc, range=cblims)
            else
                GMT.makecpt(cmap=cbcolor, truncate=cbtrunc, range=cblims, inverse=:c)
            end
            
            GMT.scatter(lon[idx_box], lat[idx_box], zcolor=LF,
	                proj=(name=:Stereographic, center=[mean(LonLim), 90]),
                        figsize=(6,6),
                        xaxis=(axes=:Sn, annot=2, ticks="60m", grid=false),
                        yaxis=(axes=:We, annot=0.5, ticks="10m", grid=false),
                        region=(LonLim..., LatLim...), cmap=cbar, markersize=pxsize,
                        marker=:square, title=wind_title_str,
	                par=(FONT_ANNOT_PRIMARY=7, FONT_LABEL=10))
            GMT.basemap!(region=(LonLim..., LatLim...),
                         proj=(name=:Stereographic, center=[mean(LonLim),90]),
                         figsize=(6,6),
	                 #xaxis=(axes=:Sn, annot=5, ticks="1deg", grid=false),
	                 rose = (outside=true, anchor=:TR, width=1., frame=false, offset=(-0.5, -1.0), label=true)
	                 )
            day0 = max(1, i_rv-100)
            GMT.plot!(RV[:lon][day0:i_rv], RV[:lat][day0:i_rv], show=0)
            GMT.plot!(lon_wind, lat_wind, lc=:grey, linestyle=:line, show=0)
            GMT.colorbar!(pos=(anchor=:CR, length=(4,0.15), offset=(0.2 ,0)), cmap=cbar, frame=(annot="$(cbbins)"), par=(FONT_ANNOT_PRIMARY=8, FONT_LABEL=10))#, xlabel="Lead Fraction", font="9p")
            GMT.plot!(lon_circ, lat_circ, linestyle=:dash, lc=:lightred)
            ##GMT.plot!([Plon...], [Plat...], linestyle=:dash, lc=:lightblue)

            GMT.plot!(rv_lon, rv_lat, marker=:star, markersize=0.2, mc=:red, show=0, fmt=:png, savefig=plot_name)

            iplt += 1
        end
    end  # over iday
    
##    end # over daysinmonth
##end # over months and year

end # over variables

save_object("data/seaice_winddir_18112019.jld2", dist_wdir)

# end of script
