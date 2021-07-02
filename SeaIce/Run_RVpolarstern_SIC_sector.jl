# Script to extract the SeaIce Concentration from AMSR2 data for ChukchiBeaufort Sea.
# 
using NCDatasets
using Dates
#using MAT
using JLD2

# Reading the Latitude and Longitude for AMSR2 ChukchiBeaufort region:
const DATA_PATH = "/home/psgarfias/LIM/remsens/utqiagvik-nsa/";
const SENS_PATH = "SeaIce/amsr2/";
const PROC_PATH = joinpath(SENS_PATH, "HDF")
#const WIND_PATH = joinpath(DATA_PATH, "");
const latlon_file = joinpath(DATA_PATH,
                             PROC_PATH,
                             "LongitudeLatitudeGrid-n3125-ChukchiBeaufort.h5");
ncxy = NCDataset(latlon_file, "r");
lat = Array{Float64}(ncxy["Latitudes"][:,:]);
lon = Array{Float64}(ncxy["Longitudes"][:,:]);
close(ncxy)

include("./SEAICE_tools.jl")

# Define coordinates for the North Slope Alaska site:
nsa_lat = 71.323e0;
nsa_lon = -156.609e0;

# Parameters for the semi-circle from west-to-east and R_lim km radius:
R_lim = 50f0;     # Radius for area under interst [km]
θ₀ = -125f0;      # initial azimuth
θ₁ = 110f0;       # end azimuth


θ_all, ρ_all = LonLat_To_CenteredPolar(nsa_lon, nsa_lat, lon, lat);

idx = Get_Sector_Indexes(θ₀, θ₁, R_lim=50f0);

# mat_file = joinpath(WIND_PATH, PROC_PATH, "nsa_sector_lonlat.mat");
mat_file = joinpath(DATA_PATH, SENS_PATH, "JLD", "nsa_sector_lonlat.jld2");

if !isfile(mat_file)
    tmp = Tuple.(idx);
    Sector_Coordinate = Dict(
        "Center" => [nsa_lon, nsa_lat],
        "Radius" => ρ_all[idx],
        "Azimuth"=> θ_all[idx],
        "Indexes" => [first.(tmp) last.(tmp)],
        "Longitudes" => lon[idx],
        "Latitudes"  => lat[idx],
        "Source_file"=> latlon_file
    );
    
    # matwrite(mat_file, Sector_Coordinate)
    @save mat_file Sector_Coordinate
end

# Read radiosonde file and extract winddir,
# Then extract the winddir angle and store as range dependent, for 50 and 100km.
#

for yy ∈ (2018)
    for mm ∈ (11,12)
        for dd ∈ (1:31)

            # Getting SeaIce Concentration file:
            sic_filen = getFilePattern(DATA_PATH, PROC_PATH, yy, mm, dd);
            
            if !isfile(sic_filen)
                println("$sic_filen does not exist!")
                continue
            else
                println("Working on $sic_filen :)")
            end

            # Reading wind direction from Radiosonde datafiles:
            rs_filen = getFilePattern(DATA_PATH, "INTERPOLATEDSONDE", yy, mm, dd);

            rs_data = getSondeData(rs_filen);

            IVT = getIVT(rs_data)
            IVTmax, WDmax, WSmax = getMaxIVT_θ(IVT, rs_data.WD, rs_data.WS)

            idx_wivt = map(x->Get_Azimuthal_Indexes(x, θ_all[idx], ρ_all[idx]),
                          WDmax);

            idx_hgt = 5
            idx_wfix = map(x->Get_Azimuthal_Indexes(x, θ_all[idx], ρ_all[idx]),
                          rs_data.WD[idx_hgt,:]);

            # Reading corresponding SIC data for the selected sector with idx indexes:
            SIC = NCDataset(sic_filen, "r") do ncin
                ncin["ASI Ice Concentration"][idx];
            end

            # creating time dependent radial SIC(IVTmax) and related variables:
            ρ_w, θ_w, SIC_w, lon_w, lat_w, SIC2D_w, range = Get_MATVAR_From_Index(
                idx_wivt,
                lon[idx],
                lat[idx],
                SIC);

            # creating time dependent radial SIC(IVTmax) and related variables:
            ρ_h, θ_h, SIC_h, lon_h, lat_h, SIC2D_h, range = Get_MATVAR_From_Index(
                idx_wfix,
                lon[idx],
                lat[idx],
                SIC);
            
            # Storing to daily data files as MAT files:
            local TMP_PATH = joinpath(DATA_PATH, SENS_PATH, "JLD", @sprintf("%04d",yy))
            # local TMP_PATH = joinpath("/tmp/", @sprintf("%04d",yy))
            if !isdir(TMP_PATH)
                mkdir(TMP_PATH)
            end

            # Sotre the data as DataFrames:
            local mat_file = joinpath(TMP_PATH,
                                      @sprintf("nsasic.c1.%04d%02d%02d.jld2",yy,mm,dd) );
            jldopen(mat_file, "w") do file
                file["time"] = rs_data.time
                file["range"] = range
                file["idx_wind"] = idx_wivt
                file["SIC"] = SIC
                file["SIC_WD"] = SIC_w
                file["SIC_WD_h"] = SIC_h
                file["SIC2D_WD"] = SIC2D_w
                file["SIC2D_WD_h"] = SIC2D_h
                file["rho_w"] = ρ_w
                file["theta_w"] = θ_w
                file["lon_w"] = lon_w
                file["lat_w"] = lat_w
                file["IVT"] = IVT[1,:]
                file["IVTmax"] = IVTmax
                file["WD"] = WDmax #rs_data.WD[idx_hgt,:]
                file["WS"] = WSmax #rs_data.WS[idx_hgt,:]
                file["WD_h"] = rs_data.WD[idx_hgt,:]
                file["WS_h"] = rs_data.WS[idx_hgt,:]
            end
            
            # Store the data as MatLab MAT
            ##local mat_file = joinpath(TMP_PATH,
            ##                          @sprintf("nsasic.c1.%04d%02d%02d.mat",yy,mm,dd) );
            ##
            ##
            ##fmat = matopen(mat_file, "w");
            ##write(fmat, "time", rs_data.time);
            ##write(fmat, "range", range);
            ##write(fmat, "idx_wind", idx_w50);
            ##write(fmat, "SIC", SIC);
            ##write(fmat, "SIC_WD", SIC_w);
            ##write(fmat, "SIC2D_WD", SIC2D_w);
            ##write(fmat, "rho_w", ρ_w);
            ##write(fmat, "theta_w", θ_w);
            ##write(fmat, "lot_w", lon_w);
            ##write(fmat, "lat_w", lat_w);
            ##write(fmat, "WD", rs_data.WD[idx_hgt,:]);
            ##close(fmat);

        end
    end
end

# end of scrip

## extra tools
# const MATLAB_EPOCH = Dates.DateTime(-1,12,31)
# 
# date2num(d::Dates.DateTime) = Dates.value(d-MATLAB_EPOCH)/(1000*60*60*24)
# num2date(n::Number) =  MATLAB_EPOCH + Dates.Millisecond(round(Int64, n*1000*60*60*24))
