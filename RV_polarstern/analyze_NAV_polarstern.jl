# script to extract the coordinates from RV polarstern during MOSAiC

const BASE_PATH = "/home/psgarfias/LIM/remsens"
const CAMPAIGN="arctic-mosaic"


using NCDatasets
using Dates
using Printf

const PROD_PATH = "NAV"

include("/home/psgarfias/LIM/repos/ARMtools/src/ARMtools.jl")


years = (2019);
months = (11);
dd = 9

for yy ∈ years
    for mm ∈ months
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
            global lat = nc["lat"][1:50:end]
            global lon = nc["lon"][1:50:end]
            global alt = nc["alt"][1:50:end]
        end
    end
end

# end of script
