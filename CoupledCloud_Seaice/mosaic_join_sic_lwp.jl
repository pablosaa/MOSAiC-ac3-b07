#!/home/psgarfias/.local/julia-1.6.1/bin/julia
#
# script to join dataset of SIC and LWP for NSA site
using JLD2
using ARMtools
using CloudnetTools
using Printf
using Interpolations
using Dates
using Statistics

# defining constant PATHS
const CAMPAIGN = "arctic-mosaic"
const BASE_PATH = "/projekt2/ac3data/B07-data";
const DATA_PATH= joinpath(BASE_PATH, CAMPAIGN); #homedir(), "LIM/data", CAMPAIGN)
const LPR_PATH = joinpath(BASE_PATH, "LP", CAMPAIGN)
const CLT_PATH = joinpath(DATA_PATH, "CloudNet/output")

include("./aux_functions.jl");

# defining period of interest:

years = (2019, 2020)
months_years = ((11,years[1]),(12,years[1]),(1,years[2]), (2, years[2]), (3, years[2]), (4, years[2]))
days = (1:31)

# Creating Dataset dict:
NSAdat = Dict(:datum => [], :lwp =>[], :iwp =>[], :CBH =>[], :DCO=>[],
              :Nlwc=>[], :Niwc =>[], :HWVT =>[],
              :SIC => [], :aveSIC=>[], :sclcflag => [], :sicflag =>[],
              :copflag =>[], :cbhflag =>[]);


for (mm,yy) ∈ months_years

# find all days available for given mm and yy:
for dd ∈ days

    lwp_file = ARMtools.getFilePattern(CLT_PATH, "TROPOS", yy, mm, dd, fileext="_lwc-scaled-adiabatic.nc")
    isnothing(lwp_file) && continue
    
    #sic_file = @sprintf("SeaIce/amsr2/JLD/%04d/nsasic.c1.%04d%02d%02d.jld2",
    #                    yy, yy, mm, dd)
    #sic_file = joinpath(DATA_PATH, sic_file)

    # @sprintf("%04d/nsalp.c0.%04d%02d%02d.jld2", yy, yy, mm, dd)
    lpr_file = ARMtools.getFilePattern(LPR_PATH, "", yy, mm, dd)
    isnothing(lpr_file) && continue

    !(isfile(lwp_file) && isfile(lpr_file)) && continue
    
    ##mwr_file = try
    ##    ARMtools.getFilePattern(DATA_PATH, "MWR/RET", yy, mm, dd)
    ##catch
    ##    continue
    ##end;

    # reading data files:
    sic = try jldopen(lpr_file, "r") do fn
	Dict(
            :time => fn["time"],
            :SIC => fn["SIC"][:SIC2d],
            :errSIC => fn["SIC"][:errSIC2d],
            :range => fn["SIC"][:range],
        )
    end;
    catch e
	
        println("$yy.$mm.$dd has no SIC var");
	continue
    end

    lpr = jldopen(lpr_file, "r") do fn
        Dict(
            :time => fn["time"],
            :CBH  => fn["CBH"],
            :deco => fn["decop_hgt"],
            :idxivt=>fn["WVT"][:IdxWVT],
            :HWVT => fn["WVT"][:HWVT],
        )
    end;

    LWC = CloudnetTools.readLWCFile(lwp_file)

    IWC = CloudnetTools.readIWCFile(replace(lwp_file, "_lwc-scaled-adiabatic." => "_iwc-Z-T-method."))

    Niwc, Nlwc, flag_sclc = ratio_IWC_LWC(IWC, LWC, T_in=sic[:time])
    
    lwp = let tmp = get_variable_at_time(sic[:time], LWC[:time], LWC[:LWP])
	[ismissing(x) ? NaN : 1f3x for x in tmp]
    end;
    lwp[lwp .< 1] .= NaN;
    lwp[lwp .> 1000] .= NaN;

    iwp = let tmp = get_variable_at_time(sic[:time], IWC[:time], IWC[:IWP])
	[ismissing(x) ? NaN : x for x in tmp]
    end;
    iwp[iwp .< 1] .= NaN
    #iwp[iwp .> 1000] .= NaN
    
    # defining SIC_flag for every time stamp:
    sicflag = [classify_sic(ice) for ice ∈ eachcol(sic[:SIC])];
    meSIC = [FunctionCleanData(V, minV=NaN) for V in eachcol(sic[:SIC])]
    stdSIC = [FunctionCleanData(V, minV=NaN, myfun=std) for V in eachcol(sic[:SIC])]

    ##errSIC = let uncer = sic[:errSIC]
    ##    errV2 = [FunctionCleanData(V.^2, minV=NaN, myfun=sum) for V in eachcol(uncer)]
    ##    sqrt.(errV2 .+ stdSIC.^2)
    ##end


    
    cbhflag = @. !isnan(lpr[:CBH]);   #0.2 < lpr[:CBH] ≤ 8 ;  # cloud must be present
    #cbhflag = cbhflag[flag_sclc]

    copflag = lpr[:deco] .≤ lpr[:HWVT]; # .≤ 0.2;  # decopling height must be below 20m

    # appeinding into NSA dataset 172800
    append!(NSAdat[:datum], sic[:time][cbhflag])
    append!(NSAdat[:lwp], lwp[cbhflag])
    append!(NSAdat[:iwp], iwp[cbhflag])
    append!(NSAdat[:CBH], lpr[:CBH][cbhflag])
    append!(NSAdat[:HWVT], lpr[:HWVT][cbhflag])
    append!(NSAdat[:DCO], lpr[:deco][cbhflag])
    append!(NSAdat[:sicflag], sicflag[cbhflag])
    append!(NSAdat[:copflag], copflag[cbhflag])
    append!(NSAdat[:cbhflag], cbhflag[cbhflag])
    append!(NSAdat[:sclcflag], flag_sclc[cbhflag])
    append!(NSAdat[:Nlwc], Nlwc[cbhflag])
    append!(NSAdat[:Niwc], Niwc[cbhflag])
    push!(NSAdat[:SIC], sic[:SIC][:, cbhflag])  #meSIC[cbhflag])
    append!(NSAdat[:aveSIC], meSIC[cbhflag])
end
end

SIC2d = hcat(NSAdat[:SIC]...);

## Plotting stats
using StatsPlots
default(titlefont = (20, "times"), legendfontsize = 18, guidefont = (18, :black), tickfont=(12, :gray3))

# ** SIC FLAGs: 0=> Open water, 1=> Open Seaice, 2=> Cover Ice, 3=> Land 
# For LWP
idx0 = @. NSAdat[:copflag] & (NSAdat[:sicflag]==0);
idx1 = @. NSAdat[:copflag] & (NSAdat[:sicflag]==1);
idx2 = @. NSAdat[:copflag] & (NSAdat[:sicflag]==2); # immer 0 falls coupled
idx3 = @. NSAdat[:copflag] & (NSAdat[:sicflag]==3); # immer 0 falls coupled

dcidx0 = @. (!NSAdat[:copflag]) & (NSAdat[:sicflag]==0);
dcidx1 = @. (!NSAdat[:copflag]) & (NSAdat[:sicflag]==1);
dcidx2 = @. (!NSAdat[:copflag]) & (NSAdat[:sicflag]==2);  # immer 0 falls de-coupled
dcidx3 = @. (!NSAdat[:copflag]) & (NSAdat[:sicflag]==3); # immer 0 falls de-coupled

posco0 = fill(0.99, length(idx0))
posdc0 = fill(1.01, length(dcidx0))
posco1 = fill(1.99, length(idx1))
posdc1 = fill(2.01, length(dcidx1))
posco2 = fill(2.99, length(idx2))
posdc2 = fill(3.01, length(dcidx2))
#posco3 = fill(3.99, length(idx2))
#posdc3 = fill(4.01, length(dcidx2))

str_period = NSAdat[:datum] |> extrema .|> x->Dates.format(x, "yyyy-uu-dd")

### FOR LIQUID WATER:
# open water
pltlwp = violin(posco0, CleanData(NSAdat[:lwp][idx0]), side=:left, label="coupled");
violin!(posdc0, CleanData(NSAdat[:lwp][dcidx0]), side=:right, label="de-coupled", color=RGB(.4, .4, .7));

# ice-open
violin!(posco1, CleanData(NSAdat[:lwp][idx1]), side=:left, label="coupled", color=:lightblue1);
violin!(posdc1, CleanData(NSAdat[:lwp][dcidx1]), side=:right, label="de-doupled", color=:lightblue3);

# ice-covered
violin!(posco2, CleanData(NSAdat[:lwp][idx2]),
        side=:left, label="coupled", color=RGB(.8, .8, .8));
violin!(posdc2, CleanData(NSAdat[:lwp][dcidx2]),
        side=:right, label="de-coupled", color=RGB(.6, .6, .6));

# over land
#violin!(posco3, CleanData(NSAdat[:lwp][idx3]),
#        side=:left, color=:lightgreen, label="coupled");
#violin!(posdc3, CleanData(NSAdat[:lwp][dcidx3]),
#        side=:right, label="de-coupled", color=:green);

str_title = @sprintf(", from %s to %s", str_period[1], str_period[2])
plot!(pltlwp, xticks=((1:3),("SIC ≤ 15%", "15% < SIC < 95%","SIC ≥ 95%")),legend=:outertopright, size=(1200, 700), ylim=(0.9, 200), xlim=(0, 4),
      ylabel="LWP / g m⁻²", ylabelfont=font(20, "Helvetica"), tickdir=:out,
      xlabel="Sea Ice Classification",
      xtickfont = font(14, "Courier"), ytickfont = font(18, "Courier"),
      xrot = 30, title="Liquid cloud water path "*str_title,
      left_margin = 10Plots.mm, bottom_margin=20Plots.mm);

savefig(pltlwp, @sprintf("./plots/%4d-%4d_LWP_surface_class.png",months_years[1][2], months_years[end][2]))

#### FOR ICE WATER CONTENT:
pltiwp = violin(posco0, CleanData(NSAdat[:iwp][idx0]), side=:left, label="coupled");
violin!(posdc0, CleanData(NSAdat[:iwp][dcidx0]), side=:right, label="de-coupled", color=RGB(.4, .4, .7));

# ice-open
violin!(posco1, CleanData(NSAdat[:iwp][idx1]), side=:left, label="coupled", color=:lightblue1);
violin!(posdc1, CleanData(NSAdat[:iwp][dcidx1]), side=:right, label="de-doupled", color=:lightblue3);

# ice-covered
violin!(posco2, CleanData(NSAdat[:iwp][idx2]),
        side=:left, label="coupled", color=RGB(.8, .8, .8));
violin!(posdc2, CleanData(NSAdat[:iwp][dcidx2]),
        side=:right, label="de-coupled", color=RGB(.6, .6, .6));

# over land
#violin!(posco3, CleanData(NSAdat[:iwp][idx3]),
#        side=:left, color=:lightgreen, label="coupled");
#violin!(posdc3, CleanData(NSAdat[:iwp][dcidx3]),
#        side=:right, label="de-coupled", color=:green);

plot!(pltiwp, xticks=((1:3),("SIC ≤ 15%", "15% < SIC < 95%","SIC ≥ 95%")),legend=:outertopright, size=(1200, 700), ylim=(0.9, 250), xlim=(0, 4),
      ylabel="IWP / g m⁻²", ylabelfont=font(20, "Helvetica"), tickdir=:out,
      xlabel="Sea Ice Classification",
      xtickfont = font(14, "Courier"), ytickfont = font(18, "Courier"),
      xrot = 30, title="Ice water path"*str_title,
      left_margin = 10Plots.mm, bottom_margin=20Plots.mm);

savefig(pltiwp, @sprintf("./plots/%4d-%4d_IWP_surface_class.png",months_years[1][2], months_years[end][2]))

## **********************
## FOR CBH .....
pltcbh = violin(posco0, CleanData(NSAdat[:CBH][idx0]), side=:left, label="coupled");
violin!(posdc0, CleanData(NSAdat[:CBH][dcidx0]), side=:right, label="de-coupled", color=RGB(.4, .4, .7));

# ice-open
violin!(posco1, CleanData(NSAdat[:CBH][idx1]), side=:left, label="coupled", color=:lightblue1);
violin!(posdc1, CleanData(NSAdat[:CBH][dcidx1]), side=:right, label="de-doupled", color=:lightblue3);

# ice-covered
violin!(posco2, CleanData(NSAdat[:CBH][idx2]),
        side=:left, label="coupled", color=RGB(.8, .8, .8));
violin!(posdc2, CleanData(NSAdat[:CBH][dcidx2]),
        side=:right, label="de-coupled", color=RGB(.6, .6, .6));

# over land
#violin!(posco3, CleanData(NSAdat[:CBH][idx3]),
#        side=:left, color=:lightgreen, label="coupled");
#violin!(posdc3, CleanData(NSAdat[:CBH][dcidx3]),
#        side=:right, label="de-coupled", color=:green);

plot!(pltcbh, xticks=((1:3),("SIC ≤ 15%", "15% < SIC < 95%","SIC ≥ 95%")),legend=:outertopright, size=(1200, 700), ylim=(0, 6), xlim=(0, 4),
      ylabel="CBH / km", ylabelfont=font(20, "Helvetica"), tickdir=:out,
      xlabel="Sea Ice Classification",
      xtickfont = font(14, "Courier"), ytickfont = font(18, "Courier"),
      xrot = 30, title="Cloud Base Height"*str_title,
      left_margin = 10Plots.mm, bottom_margin=20Plots.mm);

savefig(pltcbh, @sprintf("./plots/%4d-%4d_CBH_surface_class.png",months_years[1][2], months_years[end][2]))
# ----/


## **********************
## FOR peak WVT .....
pltwvt = violin(posco0, CleanData(NSAdat[:HWVT][idx0]), side=:left, label="coupled");
violin!(posdc0, CleanData(NSAdat[:HWVT][dcidx0]), side=:right, label="de-coupled", color=RGB(.4, .4, .7));

# ice-open
violin!(posco1, CleanData(NSAdat[:HWVT][idx1]), side=:left, label="coupled", color=:lightblue1);
violin!(posdc1, CleanData(NSAdat[:HWVT][dcidx1]), side=:right, label="de-doupled", color=:lightblue3);

# ice-covered
violin!(posco2, CleanData(NSAdat[:HWVT][idx2]),
        side=:left, label="coupled", color=RGB(.8, .8, .8));
violin!(posdc2, CleanData(NSAdat[:HWVT][dcidx2]),
        side=:right, label="de-coupled", color=RGB(.6, .6, .6));

plot!(pltwvt, xticks=((1:3),("SIC ≤ 15%", "15% < SIC < 95%","SIC ≥ 95%")),legend=:outertopright, size=(1200, 700), ylim=(0, 3), xlim=(0, 4),
      ylabel="Height / km", ylabelfont=font(20, "Helvetica"), tickdir=:out,
      xlabel="Sea Ice Classification",
      xtickfont = font(14, "Courier"), ytickfont = font(18, "Courier"),
      xrot = 30, title="Height of WVT peak"*str_title,
      left_margin = 10Plots.mm, bottom_margin=20Plots.mm);

savefig(pltwvt, @sprintf("./plots/%4d-%4d_HWVT_surface_class.png",months_years[1][2], months_years[end][2]))
# ----/

#savefig(pwv, "./plots/iwv_$(yy)_decoup.png")

# *****************************************************
# mounthly
pltsic0 = plot();
pltsic1 = plot();
pltsic2 = plot();
pltcbh0 = plot();
pltcbh1 = plot();
pltcbh2 = plot();
pltice0 = plot();
pltice1 = plot();
pltice2 = plot();

t_mm = 0;
t_month = []
for (mm, yy) ∈ months_years
    Nmonths = length(months_years)
   # break
    idxmm = NSAdat[:datum] .|> Month .|> x->x.value==mm;
    idxmmco1 = findall(idxmm .& idx1)
    idxmmdc1 = findall(idxmm .& dcidx1)

    idxmmco2 = idxmm .& idx2
    idxmmdc2 = idxmm .& dcidx2

    append!(t_month, mm);
    str_month = Dates.format.(DateTime.(Month.(t_month)), "uuu")
    global t_mm += 1;

    str_xticks = Dates.format(DateTime(Year(yy), Month(mm)), "uuu/yyyy")

    str_colabel = t_mm==1 ? "coupled" : false
    str_decolabel = t_mm==1 ? "de-coupled" : false
    
    ### FOR LWP:
    # ice-opened
    violin!(pltsic1, fill(t_mm, Nmonths), CleanData(NSAdat[:lwp][idxmmco1]), side=:left, color=:lightblue1, label=str_colabel, title="15% < SIC < 95%");
    violin!(pltsic1, fill(t_mm, Nmonths), CleanData(NSAdat[:lwp][idxmmdc1]), side=:right, color=:lightblue3, label=str_decolabel, tick_dir=:out,
            xticks = ((1:t_mm), str_month), ylim=(0, 200), xlim = (0, Nmonths+1),
            ylabel="LWP / g m⁻²");

    # ice-covered
    violin!(pltsic2, fill(t_mm, Nmonths), CleanData(NSAdat[:lwp][idxmmco2]),
            side=:left, label=str_colabel, color=RGB(.8, .8, .8), title="SIC ≥ 95%" );
    violin!(pltsic2, fill(t_mm, Nmonths), CleanData(NSAdat[:lwp][idxmmdc2]),
            side=:right, label=str_decolabel, color=RGB(.6, .6, .6), tick_dir=:out,
            xticks = ((1:t_mm), str_month), xlim = (0, Nmonths+1), ylim=(0, 200),
            ylabel="LWP / g m⁻²");

    ### FOR ICE
    # ice-opened
    violin!(pltice1, fill(t_mm, Nmonths), CleanData(NSAdat[:iwp][idxmmco1]), side=:left, color=:lightblue1, label=str_colabel, title="15% < SIC < 95%");
    violin!(pltice1, fill(t_mm, Nmonths), CleanData(NSAdat[:iwp][idxmmdc1]), side=:right, color=:lightblue3, label=str_decolabel, tick_dir=:out,
            xticks = ((1:t_mm), str_month), ylim=(0, 250),
            xlim = (0, Nmonths+1), ylabel="IWP / g m⁻²");

    # ice-covered
    violin!(pltice2, fill(t_mm, Nmonths), CleanData(NSAdat[:iwp][idxmmco2]),
            side=:left, label=str_colabel, color=RGB(.8, .8, .8), title="SIC ≥ 95%" );
    violin!(pltice2, fill(t_mm, Nmonths), CleanData(NSAdat[:iwp][idxmmdc2]),
            side=:right, label=str_decolabel, color=RGB(.6, .6, .6), tick_dir=:out,
            xticks = ((1:t_mm), str_month), xlim = (0, Nmonths+1), ylim=(0, 250),
            ylabel="IWP / g m⁻²");

    ### FOR CBH:
    # ice-opened
    violin!(pltcbh1, fill(t_mm, Nmonths), CleanData(NSAdat[:CBH][idxmmco1]), side=:left, color=:lightblue1, label=str_colabel, title="15% < SIC < 95%");
    violin!(pltcbh1, fill(t_mm, Nmonths), CleanData(NSAdat[:CBH][idxmmdc1]), side=:right, color=:lightblue3, label=str_decolabel, tick_dir=:out,
            xticks = ((1:t_mm), str_month), ylim=(0, 6),
            xlim = (0, Nmonths+1), ylabel="Cloud Base / km");

    # ice-covered
    violin!(pltcbh2, fill(t_mm, Nmonths), CleanData(NSAdat[:CBH][idxmmco2]),
            side=:left, label=str_colabel, color=RGB(.8, .8, .8), title="SIC ≥ 95%" );
    violin!(pltcbh2, fill(t_mm, Nmonths), CleanData(NSAdat[:CBH][idxmmdc2]),
            side=:right, label=str_decolabel, color=RGB(.6, .6, .6), tick_dir=:out,
            xticks = ((1:t_mm), str_month), xlim = (0, Nmonths+1), ylim=(0, 6),
            ylabel="Cloud Base / km");
    
end

range_sic = collect(5:5:50)

month_ticks = (floor(NSAdat[:datum][1],Month):Month(1):ceil(NSAdat[:datum][end],Month));
month_str = Dates.format.(month_ticks, "u/yyyy");

# for LWP: SIC2d (NSAdat[:SIC][1:10:end])
plot!(pltsic0, NSAdat[:datum][1:2:end], range_sic, SIC2d[:,1:2:end], st=:heatmap, color=palette(:ice, 5), colorbar_title="SIC(θ_wind) %", ylabel = "radius/km", clim=(75, 100), xticks=(month_ticks, month_str), xlabel="Date", lw=1.5, tick_dir=:out, legend=false);
#plot!(twinx(pltsic0), NSAdat[:datum][1:end], NSAdat[:Nlwc][1:end], ylabel="MPC %", ylim=(0, 50), lc=:red, legend=false);

ll = @layout [a{0.42h}; b{0.42h}; c{0.16h}]
tmp = plot(pltsic1, pltsic2, pltsic0, layout=ll, size=(1100,900), legendfontsize=4, legend=:outertopright, left_margin = 10Plots.mm, bottom_margin=10Plots.mm);

savefig(tmp, @sprintf("./plots/%4d-%4d_LWP_monthy_series.png",months_years[1][2], months_years[end][2]));

# for ICE:
plot!(pltice0, NSAdat[:datum][1:2:end], range_sic, SIC2d[:,1:2:end], st=:heatmap, color=palette(:ice, 5), colorbar_title="SIC(θ_wind) %", ylabel = "radius/km", clim=(75, 100), xticks=(month_ticks, month_str), xlabel="Date", lw=1.5, tick_dir=:out, legend=false);
#plot!(twinx(pltice0), NSAdat[:datum][1:end], NSAdat[:Nlwc][1:end], ylabel="MPC %", ylim=(0, 50), lc=:red, legend=false);

ll = @layout [a{0.42h}; b{0.42h}; c{0.16h}]
tmp = plot(pltice1, pltice2, pltice0, layout=ll, size=(1100,900), legendfontsize=4, legend=:outertopright, left_margin = 10Plots.mm, bottom_margin=10Plots.mm);

savefig(tmp, @sprintf("./plots/%4d-%4d_ICE_monthy_series.png",months_years[1][2], months_years[end][2]));

# for CBH:
plot!(pltcbh0, NSAdat[:datum][1:2:end], range_sic, SIC2d[:,1:2:end], st=:heatmap, color=palette(:ice, 5), colorbar_title="SIC(θ_wind) %", ylabel = "radius/km", clim=(75, 100), xticks=(month_ticks, month_str), xlabel="Date", lw=1.5, tick_dir=:out, legend=false);
#plot!(twinx(pltcbh0), NSAdat[:datum][1:end], NSAdat[:Nlwc][1:end], ylabel="MPC %", ylim=(0, 50), lc=:red, legend=false);

ll = @layout [a{0.42h}; b{0.42h}; c{0.16h}]
tmp = plot(pltcbh1, pltcbh2, pltcbh0, layout=ll, size=(1100,900), legendfontsize=4, legend=:outertopright, left_margin = 10Plots.mm, bottom_margin=10Plots.mm);

savefig(tmp, @sprintf("./plots/%4d-%4d_CBH_monthy_series.png",months_years[1][2], months_years[end][2]));


##    # reading data files:
##    sic = jldopen(sic_file, "r") do fn
##        Dict(
##            :time => fn["time"],
##            :SIC => fn["SIC"],
##            :SIC2d_WD => fn["SIC2D_WD"],
##            :SIC2d_WD_h => fn["SIC2D_WD_h"],
##        )
##    end;
##
##    lpr = jldopen(lpr_file, "r") do fn
##        Dict(
##            :time => fn["time"],
##            :CBH  => fn["CBH"],
##            :deco => fn["decop_hgt"],
##            :idxivt=>fn["idx_ivt"],
##            :Ri => fn["Ri"]
##        )
##    end;
##
##    mwr = ARMtools.getMWRData(mwr_file, onlyvars=["be_lwp", "be_pwv"]);
##

##plt= boxplot(CleanData([:lwp][idx0]), notch=true, label="ice free coupled", outliers=false);
##boxplot!(CleanData([:lwp][dcidx0]), notch=true, label="ice free de-coupled", color=RGB(.6, .6, .9), outliers=false);
##boxplot!(CleanData([:lwp][idx1]), notch=true, label="open ice / coupled", color=:lightblue, outliers=false);
##boxplot!(CleanData([:lwp][dcidx1]), notch=true, label="open ice /de-doupled", color=:lightblue1, outliers=false);
##
##boxplot!(sum(idx2)>2 ? CleanData([:lwp][idx2]) : [1],
##         notch=true, label="ice cover", color=RGB(.5, .5, .5), outliers=false);
##boxplot!(sum(idx3)>2 ? CleanData([:lwp][idx3]) : [1],
##         notch=true, label="land",
##         color=:green, outliers=false);
## 
##plot!(plt, xticks=((1:6),("co","de-co","co","de-co","all", "all")),legend=:outertopright, size=(1000, 600), ylabel="log10( LWP / g m⁻²)", yscale=:log10,title="NSA total cloud water content $(mm).$(yy)")
##

# end of script
