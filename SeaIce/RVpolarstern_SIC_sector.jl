### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 8031e3e4-90a0-4dd6-bc1d-4a7d78eaeef1
begin 
	using NCDatasets
	using Navigation
	using Dates
	using Printf
	using GMT 
	using ARMtools
	using PlutoUI
end

# ╔═╡ 510534b0-cda8-4f91-a8a6-c5ad2d2e888c
include(joinpath(homedir(), "LIM/repos/SEAICEtools.jl/src/SEAICEtools.jl"))

# ╔═╡ 04ad01ac-3e33-11ec-232b-cf4eb346c75d
html"""<style>
main {
    max-width: 1070px;
	font-size: 1.20em;
}
"""

# ╔═╡ 1c5f6ac3-0723-48d4-97c6-c57922d5861b
md"""
# Test-bed for MODIS-AQUA and AMSR2 Sea Ice Concentration extraction
## Visualization for RV Polarstern centered sector SIC wind dependent
"""

# ╔═╡ 72fcd677-cca7-444c-8d7f-501391fc4453
begin
	const CAMPAIGN = "arctic-mosaic";
	const DATA_PATH = "/media/psgarfias/LaCie SSD/LIM/data/";
	const SENSOR = "modis_amsr2";
	const PROC_PATH = joinpath(DATA_PATH, CAMPAIGN, "SeaIce", SENSOR);
	const NSA_PATH = joinpath(DATA_PATH, SENSOR, "JLD");
	const LATLON_FILE = joinpath(PROC_PATH, "lonlat_nsidc_1km.nc");
end;

# ╔═╡ de8ea731-46ac-43bf-9ddf-0e1f9c337fe1
md"""
### Defining central point e.g. RV Polarstern coordinates
"""

# ╔═╡ b84066e0-34d0-45c8-9ceb-c2ce7a832058
rv_coor = Point(81.323e0, 6.609e0)

# ╔═╡ c44f04b5-97be-4d47-a855-5b1a1d50baf2
R_lim = 50f0

# ╔═╡ 65309319-8c55-4fcb-937d-57c8cec95896
function framebox_centered_at(coor::Point{T}, grid_coor::Matrix{Point{T}}; R_lim=50f0, Δr=1f0) where T<:Real
	ρ = map(gr-> abs(coor.ϕ - gr.ϕ) + abs(coor.λ - gr.λ), grid_coor) |> argmin
	return (center=ρ,
		xrange=range(ρ[1]-round(Int, R_lim/Δr), stop=ρ[1]+round(Int, R_lim/Δr)),
		yrange=range(ρ[2]-round(Int, R_lim/Δr), stop=ρ[2]+round(Int, R_lim/Δr)))
end

# ╔═╡ a4508bc3-e802-416d-83e9-cfda5da44791
sic_coor = SEAICE.read_LatLon_Bremen_product(LATLON_FILE);
	#NCDataset(LATLON_FILE, "r") do ds
		#Point.(Float64.(ds["lat"][:,:]),
		#	Float64.(ds["lon"][:,:])
		#)
	#end;

# ╔═╡ 884944db-da75-461c-b48b-eefd1376ef3f
r = framebox_centered_at(rv_coor, sic_coor)

# ╔═╡ 4291f4eb-7646-472c-ba06-62968645e938
xlon,ylat = SEAICE.Get_LonLat_From_Point(sic_coor[r.xrange, r.yrange][:])

# ╔═╡ 799ea447-d917-4d6d-a6c5-5173f539bfc8
begin
	dd = 31
	mm = 12
	yy = 2019
end;

# ╔═╡ ed7da353-838f-4a70-8981-3797398f3c76
md"""
MODIS-AMSR2 data for $(dd). $(mm). $(yy)
"""

# ╔═╡ 6b2475cb-3e13-46b5-aa29-17798e09750d
θ_all, ρ_all = SEAICE.LonLat_To_CenteredPolar(rv_coor, sic_coor);

# ╔═╡ 53a02d50-6161-4eb1-90c4-d951b39b893e
wind_ui = @bind wind_dir Slider(1:360; default=300f0, show_value=true)

# ╔═╡ f9044bad-8ebe-4988-8580-1f0f9b6a44f4
begin 
	θ₀ = wind_dir - 45f0
	θ₁ = wind_dir + 45f0
end;

# ╔═╡ 160489f8-825a-482a-a122-f7ad698a5098
idx_sector = SEAICE.Get_Sector_Indexes(θ₀, θ₁, θ_all, ρ_all, R_lim=50f0);

# ╔═╡ cbb771d6-c4b7-4c66-8f3a-09bf96510d75
typeof(idx_sector)

# ╔═╡ aa9fbbd8-8ce5-4b7e-a250-4815e3947a00
sic_filen = ARMtools.getFilePattern(DATA_PATH, "$(CAMPAIGN)/SeaIce/$(SENSOR)", yy, mm, dd)[1];

# ╔═╡ d733e3c8-295f-4fe6-8484-707072b5c3b9
# Reading corresponding SIC data for the selected sector with idx indexes:
#function read_SIC_product(sic_filen::String, idx_sector::Vector; SICPROD="ASI Ice Concentration")
#	!in(SICPROD, ("ASI Ice Concentration", "mersic", "asic", "msic")) && error("$(SICPROD) not supported!")
#	ncin = NCDataset(sic_filen, "r")
#	if haskey(ncin, SICPROD)
#         SIC = ifelse(isempty(idx_sector), ncin[SICPROD][:,:], ncin[SICPROD][idx_sector])
#	else
#		SIC = nothing	
#	end
#	
#	return SIC
#end;

# ╔═╡ e7bff3f0-f0be-46e4-b24f-36ac7e25fc0e
SIC = SEAICE.read_SIC_Bremen_product(sic_filen, idx_sector, SICPROD="mersic");

# ╔═╡ 12122654-7877-4f2b-8b22-06ec6eac8bc9
idx_radial = SEAICE.Get_Azimuthal_Indexes(wind_dir, θ_all[idx_sector], ρ_all[idx_sector]);

# ╔═╡ 93122f34-c31d-4605-8794-f4c0b1ad620c
begin
	lat = map(p->p.ϕ, sic_coor[idx_sector]);
	lon = map(p->p.λ, sic_coor[idx_sector]);
	P_circ = SEAICE.Create_Semi_Circle(rv_coor, θ₀+360, θ₁, R_lim=51e3)
	lon_circ, lat_circ = SEAICE.Get_LonLat_From_Point(P_circ)
end;

# ╔═╡ 76157aa2-6736-4c2f-8d90-0512ee8a96da
wind_ui

# ╔═╡ 1908443c-235f-4a6f-a95e-777c5966ab6e
begin
	sic_bar = GMT.makecpt(cmap=:abyss, range=(0, 100, 10), continues=true) #color=:abyss
	GMT.coast(region = (extrema(xlon)..., extrema(ylat)...), yaxis=(annot=0.5, ticks=2),
		proj = (name=:stereographic, center = [rv_coor.λ, 90], paralles=5), figsize=8,
		frame=:a10g, res=:high, area=550, land=:darkolivegreen4, water=:gray, #figscale="1:1900000",
		map_scale = "jTL+c40+w60+f+o0.7/0.5",
		rose = (inside=true, anchor=:TR, width=1.2, offset=(0.7, 0.9), label=true),
		shore=:brown, show=0)# 
	
	GMT.scatter!(lon, lat, color=sic_bar, zcolor=SIC, markersize=0.05, marker=:circle,
		xlabel="Longitude / deg", ylabel="Latitude / deg", title="SIC $(SENSOR) WD=$(wind_dir)@+o@+", show=0)
	GMT.colorbar!(pos=(anchor=:CR, length=(6,0.2)), color=sic_bar, frame=(ylabel="%",), show=0)
	GMT.plot!(lon_circ, lat_circ, lc=:red, lw=1.5, show=0)
	GMT.plot!(rv_coor.λ, rv_coor.ϕ, symbol=(symb=:arrow, size=.2), mc=:yellow, show=0)
	GMT.plot!(lon[idx_radial], lat[idx_radial], marker=:square, markeredgecolor=:lightgreen, show=1, fmt=:png)
	#GMT.histogram!(SIC[idx_radial], bin=5, region=(10, 100, 0, 150), xlabel="SIC %",
	#figsize=(3.3,2), x_offset=-1.05, y_offset=5.9, bg=:white, show=0)
end

# ╔═╡ e313fdd8-e014-4958-baca-4d04ef79050a
GMT.histogram(SIC[idx_radial], bin=5, region=(10, 100, 0, 150), xlabel="SIC %",
	figsize=(3.0, 1.5), show=1, fmt=:png)

# ╔═╡ 85968990-b147-40d0-88d5-12cd8e1bb552
#begin
#	#thr_θρ = ((24, 0, 15),(12, 15, 25),(6, 25, 35),(3, 35, 140)) #, (3, 40, 100))
#	θ_π = θ_all[idx_sector]
#	thr_θρ = ((0.30, 0, 15), (0.15, 15, 25), (0.1, 25, 35),(0.05, 35, 140))
#	U = sind.(θ_π)
#	V = cosd.(θ_π)
#	δW = (U .- sind(wind_dir)).^2 .+ (V .- cosd(wind_dir)).^2 .|> sqrt
#
#	#idx_radial = findall(<(.1), δW)
#	idx_radial = []
#	foreach(thr_θρ) do th
#		ii = findall(th[2] .< ρ_all[idx_sector] .≤ th[3])
#		jj=findall(≤(th[1]), δW[ii] ) #abs.(wind_dir .- θ_π[ii]))
#		push!(idx_radial, ii[jj])
#	end
#	idx_radial = vcat(idx_radial...)
#	
#end;

# ╔═╡ 0a71c8e4-c80f-4389-91a7-b01eef2ec1ea
md"""
----
### © 2021, [Pablo Saavedra Garfias](mailto:pablo.saavedra@uni-leipzig.de)
### University of Leipzig
### Faculty of Physics and Geosciences
### LIM
### See LICENSE
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ARMtools = "04fa4220-f7a9-42e2-a909-1083f698c312"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
GMT = "5752ebe1-31b9-557e-87aa-f909b540aa54"
NCDatasets = "85f8d34a-cbdd-5861-8df4-14fed0d494ab"
Navigation = "cc87a3b2-3706-4f35-b5b4-b2c85061916d"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
GMT = "~0.38.2"
NCDatasets = "~0.11.7"
Navigation = "~0.4.1"
PlutoUI = "~0.7.19"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ARMtools]]
deps = ["NCDatasets", "Printf", "Statistics", "Test", "Wavelets"]
git-tree-sha1 = "9430f4c6be7610dc4962ee848aacfa1a8d9e3a38"
repo-rev = "main"
repo-url = "git@github.com:pablosaa/ARMtools.jl.git"
uuid = "04fa4220-f7a9-42e2-a909-1083f698c312"
version = "0.1.0"

[[AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0bc60e3006ad95b4bb7497698dd7c6d649b9bc06"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.1"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[CFTime]]
deps = ["Dates", "Printf"]
git-tree-sha1 = "bca6cb6ee746e6485ca4535f6cc29cf3579a0f20"
uuid = "179af706-886a-5703-950a-314cd64e0468"
version = "0.1.1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f885e7e7c124f8c92650d61b9477b9ac2ee607dd"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.1"

[[ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "9a1d594397670492219635b35a3d830b04730d62"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.1"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Conda]]
deps = ["JSON", "VersionParsing"]
git-tree-sha1 = "299304989a5e6473d985212c28928899c74e9421"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.5.2"

[[DSP]]
deps = ["FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "2a63cb5fc0e8c1f0f139475ef94228c7441dc7d0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.6.10"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Documenter]]
deps = ["Base64", "Dates", "DocStringExtensions", "InteractiveUtils", "JSON", "LibGit2", "Logging", "Markdown", "REPL", "Test", "Unicode"]
git-tree-sha1 = "fb1ff838470573adc15c71ba79f8d31328f035da"
uuid = "e30172f5-a6a5-5a46-863b-614d45cd2de4"
version = "0.25.2"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "463cb335fa22c4ebacfd1faba5fde14edb80d96c"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.5"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[GMT]]
deps = ["Conda", "Dates", "Pkg", "Printf", "Statistics"]
git-tree-sha1 = "91ea96b1ae21a8f8878b7e042f76ad97ff152df1"
uuid = "5752ebe1-31b9-557e-87aa-f909b540aa54"
version = "0.38.2"

[[HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fd83fa0bde42e01952757f01149dd968c06c4dba"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.0+1"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "19cb49649f8c41de7fea32d089d37de917b553da"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.0.1"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "be9eef9f9d78cecb6f262f3c10da151a6c5ab827"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NCDatasets]]
deps = ["CFTime", "DataStructures", "Dates", "NetCDF_jll", "Printf"]
git-tree-sha1 = "5da406d9624f25909a6f556bd8d5c1deaa189ee6"
uuid = "85f8d34a-cbdd-5861-8df4-14fed0d494ab"
version = "0.11.7"

[[Navigation]]
deps = ["Documenter"]
git-tree-sha1 = "93b7f4699f2babb09c9807987827aafbfc7141e4"
uuid = "cc87a3b2-3706-4f35-b5b4-b2c85061916d"
version = "0.4.1"

[[NetCDF_jll]]
deps = ["Artifacts", "HDF5_jll", "JLLWrappers", "LibCURL_jll", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Pkg", "Zlib_jll", "nghttp2_jll"]
git-tree-sha1 = "0cf4d1bf2ef45156aed85c9ac5f8c7e697d9288c"
uuid = "7243133f-43d8-5620-bbf4-c2c921802cf3"
version = "400.702.400+0"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "e071adf21e165ea0d904b595544a8e514c8bb42c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.19"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "OffsetArrays", "RecipesBase"]
git-tree-sha1 = "0b15f3597b01eb76764dd03c3c23d6679a3c32c8"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "1.2.1"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "8de32288505b7db196f36d27d7236464ef50dba1"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.6.2"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VersionParsing]]
git-tree-sha1 = "e575cf85535c7c3292b4d89d89cc29e8c3098e47"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.1"

[[Wavelets]]
deps = ["DSP", "FFTW", "LinearAlgebra", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "e5903fb2bf93697a79d01383618ea0855256a337"
uuid = "29a6e085-ba6d-5f35-a997-948ac2efa89a"
version = "0.9.3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─04ad01ac-3e33-11ec-232b-cf4eb346c75d
# ╟─1c5f6ac3-0723-48d4-97c6-c57922d5861b
# ╠═72fcd677-cca7-444c-8d7f-501391fc4453
# ╠═8031e3e4-90a0-4dd6-bc1d-4a7d78eaeef1
# ╟─de8ea731-46ac-43bf-9ddf-0e1f9c337fe1
# ╠═b84066e0-34d0-45c8-9ceb-c2ce7a832058
# ╠═c44f04b5-97be-4d47-a855-5b1a1d50baf2
# ╠═65309319-8c55-4fcb-937d-57c8cec95896
# ╠═a4508bc3-e802-416d-83e9-cfda5da44791
# ╠═884944db-da75-461c-b48b-eefd1376ef3f
# ╠═4291f4eb-7646-472c-ba06-62968645e938
# ╠═799ea447-d917-4d6d-a6c5-5173f539bfc8
# ╟─ed7da353-838f-4a70-8981-3797398f3c76
# ╠═510534b0-cda8-4f91-a8a6-c5ad2d2e888c
# ╠═6b2475cb-3e13-46b5-aa29-17798e09750d
# ╟─53a02d50-6161-4eb1-90c4-d951b39b893e
# ╠═f9044bad-8ebe-4988-8580-1f0f9b6a44f4
# ╠═160489f8-825a-482a-a122-f7ad698a5098
# ╠═cbb771d6-c4b7-4c66-8f3a-09bf96510d75
# ╠═aa9fbbd8-8ce5-4b7e-a250-4815e3947a00
# ╠═d733e3c8-295f-4fe6-8484-707072b5c3b9
# ╠═e7bff3f0-f0be-46e4-b24f-36ac7e25fc0e
# ╠═12122654-7877-4f2b-8b22-06ec6eac8bc9
# ╠═93122f34-c31d-4605-8794-f4c0b1ad620c
# ╟─76157aa2-6736-4c2f-8d90-0512ee8a96da
# ╠═1908443c-235f-4a6f-a95e-777c5966ab6e
# ╠═e313fdd8-e014-4958-baca-4d04ef79050a
# ╟─85968990-b147-40d0-88d5-12cd8e1bb552
# ╟─0a71c8e4-c80f-4389-91a7-b01eef2ec1ea
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
