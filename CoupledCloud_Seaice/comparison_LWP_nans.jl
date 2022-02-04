### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 928f261d-2037-492d-be45-076923013642
begin
	using Plots
	using ARMtools
	import Pkg
	Pkg.add(path=joinpath(homedir(), "LIM/repos/CloudnetTools.jl"))
	using CloudnetTools
end

# ╔═╡ 7005743d-8a60-48a1-84d5-7675279f30c4
using Dates

# ╔═╡ 76918848-83ce-11ec-2311-1980fb124d5a
html"""<style>
main {
    max-width: 1070px;
	font-size: 1.60em;
}
"""

# ╔═╡ 09e750dc-00a7-44e6-bf38-5bc1cfc3e9b5
begin
	BASE_PATH = joinpath(homedir(), "LIM/data/CloudNet/arctic-mosaic/TROPOS/PS122_1/cloudnet/preliminary_shipversion/polarstern/")
	ARM_PATH = joinpath(homedir(), "LIM/data/arctic-mosaic/MWR/")
end

# ╔═╡ 15a268a6-9c7c-45c5-b245-d7cb83f12ad7
begin
	yy = 2019
	mm = 11
	dd = 14
end

# ╔═╡ 6ddc3731-4aa2-4667-889f-d2b39cbe43c0
Filecate = ARMtools.getFilePattern(BASE_PATH, "processed/categorize/", yy, mm, dd, fileext=".nc")

# ╔═╡ 4a781f1c-0e79-4d2f-8ed8-836524988ebb
Fileadia = ARMtools.getFilePattern(BASE_PATH, "products", yy, mm, dd, fileext="adiabatic.nc")

# ╔═╡ 144d4f92-d28a-479d-9425-f28209dc7842
begin
	cateLWP = CloudnetTools.readCLNFile(Filecate)
end;

# ╔═╡ f48cd22d-4928-4f70-8a88-d656d9df1b3d
adiaLWP = CloudnetTools.readLWCFile(Fileadia);

# ╔═╡ b3027243-6257-4750-9f89-8f2df28ce2fc
pltlwp=plot(adiaLWP[:time], [cateLWP[:LWP] 1f3adiaLWP[:LWP]], label=["Categorize File" "Product File"], xlabel="UTC", ylabel="LWP / g m⁻²", legendposition=:left, title="MOSAiC LWP from $(Date(adiaLWP[:time][1]))")

# ╔═╡ 6ff30684-9091-458f-aef2-1a9eea15e8e7
#savefig(pltlwp, "mosaic_LWP_example_14-11-2019.png")

# ╔═╡ Cell order:
# ╟─76918848-83ce-11ec-2311-1980fb124d5a
# ╠═928f261d-2037-492d-be45-076923013642
# ╠═09e750dc-00a7-44e6-bf38-5bc1cfc3e9b5
# ╠═15a268a6-9c7c-45c5-b245-d7cb83f12ad7
# ╠═6ddc3731-4aa2-4667-889f-d2b39cbe43c0
# ╠═4a781f1c-0e79-4d2f-8ed8-836524988ebb
# ╠═144d4f92-d28a-479d-9425-f28209dc7842
# ╠═f48cd22d-4928-4f70-8a88-d656d9df1b3d
# ╠═7005743d-8a60-48a1-84d5-7675279f30c4
# ╠═b3027243-6257-4750-9f89-8f2df28ce2fc
# ╠═6ff30684-9091-458f-aef2-1a9eea15e8e7
