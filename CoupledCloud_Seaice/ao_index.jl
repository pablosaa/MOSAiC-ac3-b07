# Script to read and plot the AO index for a given time period:

using CSV
using DataFrames
using Dates
using StatsPlots

const DATA_PATH = joinpath(homedir(), "LIM/data/synoptic");
const ao_file = joinpath(DATA_PATH, "norm.daily.ao.cdas.z1000.19500101_current.csv");

# loading the data as DF:
df = CSV.read(ao_file, DataFrame, skipto=50*12*30); #, limit=500);

# defining period to plot (in Winter bases):

for years ∈ (2012:2019)
yy = (years, years+1)
mm = (11,12,1,2,3,4)

winter_df = let selectdata(y,m) = @. (y ==yy[1] && 11≤m≤12) || (y==yy[2] && 1≤m≤4)
    local tmpdf = subset(df, [:year, :month]=> selectdata)
    transform(tmpdf, [:year, :month, :day]=> ((y,m,d)->DateTime.(y,m,d)) => :dates)
end

# Starting plotting:
tm_tick = winter_df.dates[1]:Day(14):winter_df.dates[end];
str_tick = Dates.format.(tm_tick, "dd/uuu")
range_years = Dates.format(tm_tick[1], "yyyy-")*Dates.format(tm_tick[end], "yyyy");
@df winter_df plot(:dates, :ao_index_cdas, ylim=(-4, 6.5), ylabel="AO index", xlabel="Winter $(range_years)", legend=false, xticks=(tm_tick, str_tick), xrot=30, size=(900,400), left_margin=5Plots.mm, bottom_margin=9Plots.mm, tickfontsize=11, guidefontsize=13)

    hline!([0], l=:dash, c=:gray, label=false)
    fig_file = joinpath(pwd(), "plots/ao_index", "winter_$(range_years).png")
    savefig(fig_file)
end
# end of script
