function draw_MSAiC_sic(lon, lat, SIC, P0::Dict{Symbol, Any}; pixel_size=0.2, title_str="",
	RV=Dict{Symbol, Any}(), Circ=Dict{Symbol, Any}(), CBar::GMT.GMTcpt=[], SaveFigure::String="")

	# defining default values:
	!haskey(P0, :LonLim) && (P0[:LonLim]=extrema(lon))
	!haskey(P0, :LatLim) && (P0[:LatLim]=extrema(lat))
	!haskey(P0, :symbol) && (P0[:symbol]=:asterisk)
	!haskey(Circ, :label) && (Circ[:label]="")
	!haskey(RV, :label) && (RV[:label] = "Drift")
	
# plotting the SIC pixels:
    GMT.scatter(lon, lat, color=sic_bar, zcolor=SIC, markersize=pixel_size,
		marker=:square, region = (P0[:LonLim]..., P0[:LatLim]...), frame=:g, figsize=(1.5, 7),
		proj = (name=:stereographic, center = [P0[:lon], 90], paralles=30),
		xlabel="lon", ylabel="lat",
		title=title_str, show=0)

    # showing coast (if nearby land)
    GMT.coast!(
	proj =(name=:stereographic, center = [P0[:lon], 90], paralles=30),
	region = (P0[:LonLim]..., P0[:LatLim]...), #frame=:g, #axes=:WS, annot="30m", ticks="15m"),
	xaxis=(axes=:S, annot="120m", ticks="60m"), yaxis=(axes=:W, annot="10m", ticks="5m"),
	res=:high, area=950, land=:lightbrown, figsize=(1.5, 7), #transparency=9, #figscale="1:1900000",
	rose = (outside=true, anchor=:TR, width=1.2, frame=false, offset=(-0.5, -0.9), label=true), 
	shore=:brown, show=0
    )
    #GMT.colorbar!(pos=(anchor=:CR, length=(6,0.2),offset=(1 ,0)), color=sic_bar, frame=(ylabel="%",), show=0)
	GMT.colorbar!(pos=(outside=:CR, length=(5,0.15), offset=(.3 ,-.5)), color=CBar, frame=(ylabel="%",), show=0)
    GMT.plot!(P0[:lon], P0[:lat], symbol=(symb=P0[:symbol], size=.3), mc=:red, show=0)
	GMT.plot!(Circ[:lon], Circ[:lat], lc=:gray1, lw=0.5, ls=:dash, alpha=50, show=0, label=Circ[:label])
    return GMT.plot!(RV[:lon], RV[:lat], linestyle=:line, lw=0.4, lc=:lightred, show=1, legend=(label=RV[:label], pos=:TL), fmt=:png)
    
end



# end of file
