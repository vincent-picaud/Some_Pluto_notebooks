### A Pluto.jl notebook ###
# v0.18.1

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

# ╔═╡ dedb9336-134d-4913-bff3-869cf3f416e5
md"""
A collection of Pluto tips... 
"""

# ╔═╡ 3ba39661-7654-4f88-a43a-c748de5c57f2
md"""
# Presentation mode

From: [https://github.com/JuliaPluto/PlutoUI.jl/issues/11]()

Click on "Present" to enter mode, use arrows to navigate, click on "Present" to exit mode.
"""

# ╔═╡ 1c0cbd99-c2b2-4102-a098-df67103d4e09
html"<button onclick=present()>Present</button>"

# ╔═╡ 878fd76f-f478-47cf-9279-0d254c60440e
md"""
# Slider with printed value

From: [https://github.com/JuliaPluto/PlutoUI.jl/issues/11]()
"""

# ╔═╡ 1a7764c0-9ee9-11ec-2326-99a31961d2b7
# From: https://github.com/JuliaPluto/PlutoUI.jl/issues/11
begin
	Base.@kwdef struct MySlider 
	    range::AbstractRange
	    default::Number
	end
	function Base.show(io::IO, ::MIME"text/html", slider::MySlider)
	    print(io, """
		    <input type="range" 
		    min="$(first(slider.range))" 
			step="$(step(slider.range))"
			max="$(last(slider.range))" 
			value="$(slider.default)"
			oninput="this.nextElementSibling.value=this.value">
			<output>$(slider.default)</output>
			""")
	end
end

# ╔═╡ 2b022346-1044-4f5d-95b3-1b49edbb15c0
md"""
My example: $(@bind n MySlider(range=1:10,default=5))
"""

# ╔═╡ 493dfa1a-c87e-4a8c-bf95-4f6d9a66a5f2
n

# ╔═╡ 1e86f229-bc91-4871-a883-b23a884478a3
md"""
# Pretty box

!!! example "An example"
    blabla ....

    blabla ....
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╟─dedb9336-134d-4913-bff3-869cf3f416e5
# ╠═3ba39661-7654-4f88-a43a-c748de5c57f2
# ╟─1c0cbd99-c2b2-4102-a098-df67103d4e09
# ╟─878fd76f-f478-47cf-9279-0d254c60440e
# ╟─2b022346-1044-4f5d-95b3-1b49edbb15c0
# ╠═493dfa1a-c87e-4a8c-bf95-4f6d9a66a5f2
# ╠═1a7764c0-9ee9-11ec-2326-99a31961d2b7
# ╠═1e86f229-bc91-4871-a883-b23a884478a3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
