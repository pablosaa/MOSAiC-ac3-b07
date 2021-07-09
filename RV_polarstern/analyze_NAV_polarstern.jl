# script to extract the coordinates from RV polarstern during MOSAiC

const BASE_PATH = "/home/psgarfias/LIM/remsens"
const CAMPAIGN="arctic-mosaic"


using NCDatasets
using Dates
using Printf

const PROD_PATH = "NAV"

#include("/home/psgarfias/LIM/repos/ARMtools/src/ARMtools.jl")


years = (2019);
months = (11);
days = (1:10)

name_col = (:time, :lat, :lon, :alt, :pitch, :heave, :yaw, :roll);
RVPS = Dict();
map(x->RVPS[x]=[], name_col)

for yy ∈ years
    for mm ∈ months
        for dd ∈ days
        # defining file name pattern to filter:
        file_pattern = @sprintf("%04d%02d%02d", yy, mm,dd)
        DATA_PATH = joinpath(BASE_PATH, CAMPAIGN, PROD_PATH, "$yy")
        
        # testing whether the data directory exist:
        @assert isdir(DATA_PATH) "$DATA_PATH cannot be found!"
        
        # Reading the files in data directory:
        file_list = readdir(DATA_PATH, join=true)
        file_month_list = filter(x->all(occursin.(file_pattern, x)), file_list)

        isempty(file_month_list) ? continue : println("Working on $mm")

        NCDataset(file_month_list[1]) do nc
            map(name_col) do isym
                append!(RVPS[isym], nc[Symbol(isym)][1:100:end])
            end
        end
        end
    end
end

# end of script
