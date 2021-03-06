#!/home/psgarfias/.local/julia-1.6.1/bin/julia
###!/opt/julia-1.6.0/bin/julia

# script to extract the RV Polarstern trajectory during the
# MOSAiC expedition and stored as single database

using ARMtools
using Dates
using JLD2
using Printf

const CAMPAIGN = "arctic-mosaic";
const DATA_PATH = "/projekt2/remsens/data_new/site-campaign/"; #"/media/psgarfias/LaCie SSD/LIM/data/";
const DATA_SOURCE = "NAV";
const OUT_PATH = joinpath(pwd(), "data");
#const LATLON_FILE = joinpath(PROC_PATH, "lonlat_nsidc_1km.nc");

jahre = (2020)
monaten = (1:6)
tage = (1:31)

key2keep = (:time, :lon, :lat, :alt, :heave, :roll, :pitch, :yaw);
RV = Dict(x=>[] for x=key2keep)

for yy ∈ jahre
    for mm ∈ monaten
        for dd ∈ tage
            tmp = ARMtools.getFilePattern(DATA_PATH*CAMPAIGN, DATA_SOURCE, yy, mm, dd)
            isnothing(tmp) && (println("No file on $(dd).$(mm).$(yy)"); continue)

            flag = typeof(tmp) <: Vector
            nav_file = flag ? tmp : [tmp]

            foreach(nav_file) do infile
                let nav = ARMtools.getNAVData(infile, time_res=Hour(6));

                    foreach(key2keep) do kk
                        append!(RV[kk], nav[kk])
                    end
                end
            end
            # closing date loops:
        end
    end
end
# Storing the track:
jldsave(joinpath(OUT_PATH, "RVpolarstern_track_2020.jld2"); RV);


# end of script
