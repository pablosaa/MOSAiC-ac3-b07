#!/opt/julia-1.7.2/bin/julia
# script to create SIC timeseries for MOSAiC using data from AMSR2-MODIS and OSISAF.

# loading needed packages:
using ARMtools
using Plots
using Dates
using Navigation, NCDatasets
using StatsBase
using DataFrames, CSV
#using GMT

SEAICE = include(joinpath(homedir(), "LIM/repos/SEAICEtools.jl/src/SEAICEtools.jl"));
include("aux_gmt_functions.jl")
# --

# defining consants/paths
const CAMPAIGN = "arctic-mosaic";
const PRODUCT = "SeaIce";
const BASE_PATH = joinpath(homedir(), "LIM/data/B07");
const SENSOR = ["modis_amsr2", "osisaf"];
const SIC_VAR = ["mersic", "ice_conc"];
const PROC_PATH = joinpath(BASE_PATH, PRODUCT);
# --

## physical constants:
const R_lim = 50e3;   # radius around RV polarstern

# loading time invariant variables:
LATLON_FILE = Vector{String}(undef, 2);
const LATLON_FILE[1] = joinpath(PROC_PATH, SENSOR[1], "lonlat_nsidc_1km.nc");
const LATLON_FILE[2] = ARMtools.getFilePattern(PROC_PATH, SENSOR[2], 2020, 4, 15);

# reading the coordinates for both products (those are fixed on time):
sic_coor = [SEAICE.read_LatLon_Bremen_product(LATLON_FILE[i]) for i ∈ 1:2];

# defining DataFrame to collect data
MOS = DataFrame();

for (mm, yy) ∈ [(10,2019), (11,2019),(12,2019), (1,2020), (2,2020), (3,2020), (4,2020), (5,2020)]

    ## loading RV polarstern NAV data:
    RV = load("../RV_polarstern/data/RVpolarstern_track_$(yy).jld2", "RV");

    for dd=1:daysinmonth(yy, mm)

        # getting the index for the RV position of the day
        iday = let thisday = DateTime(yy, mm, dd)
            findall(thisday .≤ RV[:time] .< (thisday + Day(1) + Minute(1)))
        end
        Nnav = length(iday);

        # finding products file names if available for given date:
        sic_file = [ARMtools.getFilePattern(PROC_PATH, SENSOR[i], yy, mm, dd) for i ∈ 1:2];
        (isnothing.(sic_file) |> any) && continue
        
        for i_rv ∈ iday
            # RV Polarstern coordinates:
            Pstern = Point(RV[:lat][i_rv], RV[:lon][i_rv])

            # R_lim wide box centered at RV location:
            LonLim, LatLim = SEAICE.estimate_box(Pstern, R_lim, δR=5e3)

            # idx_box contains the indexes of SIC coordinates inside LonLim, LatLim box:
            idx_box = [SEAICE.extract_LonLat_Box(LonLim, LatLim, xy) for xy ∈ sic_coor]

            # Reading SIC data from selectec box:
            SIC = map(enumerate(sic_file)) do (i, fsic)
                if i==1
                    SEAICE.read_SIC_Bremen_product(fsic, idx_box[i], SICPROD=SIC_VAR[i])
                else
                    SEAICE.read_SIC_OSISAF_product(fsic, idx_box[i], SICPROD=SIC_VAR[i])
                end
            end

            # estimating the SIC statistics for the box:
            μ_sic = mean.(SIC)
            σ_sic = std.(SIC)
            ϕ_sic = quantile.(SIC)
            df = DataFrame(time= RV[:time][i_rv],
                           asi_μ=μ_sic[1], asi_σ=σ_sic[1],
                           osi_μ=μ_sic[2], osi_σ=σ_sic[2],
                           med_asi=ϕ_sic[1][3], med_osi=ϕ_sic[2][3])
            
            global MOS = vcat(MOS, df)
            #println(DateTime(yy,mm,dd), "->", μ_sic[1],"±",σ_sic[1], ", ", μ_sic[2],"±", σ_sic[2])

            #println(ϕ_sic)
        end
    end
end
CSV.write("data/modis_amsr2_osisaf_stats.csv", MOS);
# --

# reding Sea Ice data for the given date:


# end of script
