### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 45ee7653-8bdf-472d-aa39-7a5be5f03dbe
begin 
	using PlutoUI

	PlutoUI.TableOfContents(aside=true)
end

# ╔═╡ 15f16d2d-3295-48e2-99a7-8a21e8c28ca0
using StaticArrays

# ╔═╡ 841ab7fb-aefb-434d-8d60-db585132dd19
md"""
!!! example ""
	**Introduction à la Différentiatin Automatique (suite)**

	vincent.picaud@cea.fr
"""

# ╔═╡ 64776627-30c6-4ee1-8df3-6be2829bd9c2
md"""
# Approche directe
"""

# ╔═╡ a3c76871-a056-40b9-8d01-151546b5a7e0
md"""
## Rappel sur les $dx$, $dy$ etc...

On rappelle ici ce que représente $dx$, $dy$ etc... 

Soient $X=\left(\begin{array}{c}a\\b\end{array}\right)\in\mathbb{R}^2$ et les 		applications coordonnées $\mathrm{p}^i$ : 

```math
\mathrm{p}^1(X)=(1\  0)\cdot X
```
```math
\mathrm{p}^2(X)=(0\ 1)\cdot X
```
**Remarque :** si l'on regroupe les $\mathrm{p}^i$, l'application $\mathrm{p}$ est simplement l'identité $\mathrm{p}=\mathrm{Id}$.

Ce sont des applications linéaires, donc :
```math
\mathrm{dp}^1=(1,\  0),\ \ \mathrm{dp}^2=(0,\ 1),\ \mathrm{dp}=\mathrm{Id}
```

Si l'on considère $f:\mathbb{R}^2\rightarrow R$, alors :
```math
f(a,b)=f(\mathrm{p}^1(X),\mathrm{p}^2(X))=f\circ \mathrm{p}(X)
```
la règle de la chaine donne 
```math
\begin{align}
\mathrm{d}(f\circ \mathrm{p})\cdot \delta X &= \mathrm{d}f\circ \mathrm{dp} \cdot \delta X \\
&= 
\left(\begin{array}{cc}\partial_1 f & \partial_2 f \end{array}\right)\cdot \left(\begin{array}{cc}1 & 0 \\ 0 & 1 \end{array}\right) \cdot \delta X \\
    &= \partial_1 f\ \left(\begin{array}{cc}1 & 0 \end{array}\right) \cdot \delta X + \partial_2 f\ \left(\begin{array}{cc}0 & 1 \end{array}\right) \cdot \delta X \\
    &= \partial_1 f\ \underbrace{\mathrm{dp}^1 \cdot \delta X}_{\mathrm{d}x} + \partial_2 f\ \underbrace{\mathrm{dp}^2 \cdot \delta X }_{\mathrm{d}y}
\end{align}
```
ce qui est traditionnellement noté 
```math
\mathrm{d}f = \partial_x f\ \mathrm{d}x + 	\partial_y f\ \mathrm{d}y
```
ou encore [1]
```math
\mathrm{d}f = \partial_1 f\ \mathrm{d}^1 + \partial_2 f\ \mathrm{d}^2
```
Dans tous les cas **les $\mathrm{d}^i$ sont simplement les différentielles des applications coordonnées $\mathrm{p}^i$**.

---

1. [Bossavit, Alain. (1991). Differential Geometry for the student of numerical methods in Electromagnetism](https://www.researchgate.net/publication/200018385_Differential_Geometry_for_the_student_of_numerical_methods_in_Electromagnetism)
"""

# ╔═╡ 2540688c-81be-4fe7-89b5-9976f44152e8
md"""
## Principe de l'approche

On se concentre sur le mode direct et l'on cherche un moyen commode de réaliser les calculs de $d\Phi_k\circ \dots \circ d\Phi_2 \circ d\Phi_1$ en même temps que $\Phi_k\circ \dots \circ \Phi_2 \circ \Phi_1$

Il suffit de regarder ce qu'il se passe pour $\Phi\circ\varphi$ où $\Phi:\mathbb{R}^n\rightarrow\mathbb{R}$ et $\varphi:\mathbb{R}^m\rightarrow\mathbb{R}^n$

La règle de la chaine conduit à :
```math
\begin{align}
\mathrm{d}(\Phi\circ\varphi) &= \sum_{i=1}^m \partial_i (\Phi\circ\phi) \mathrm{d}^ix \\
                 &= \sum_{i=1}^m \left(\sum_{j=1}^n \partial_j \Phi \partial_i \varphi^{\,j}\right) \mathrm{d}^ix \\
&= \sum_{j=1}^n \partial_j \Phi \left(\sum_{i=1}^m  \partial_i \varphi^{\,j} \mathrm{d}^i x\right) \\
&= \sum_{j=1}^n \partial_j \Phi\  \mathrm{d}\varphi^{\,j}
\end{align}
```

La formule $\mathrm{d}(\Phi\circ\varphi) = \sum_{j=1}^n \partial_j \Phi\ \mathrm{d}\varphi^{\,j}$ est le cas général où les applications coordonnées $\mathrm{p}^i$ ont remplacées par des applications quelconques $\varphi^i$. 

On rappelle que pour une application $\Phi: \mathbb{R}^m\rightarrow\mathbb{R}$, alors $\mathrm{d}\Phi$ est une forme linéaire de $\mathbb{R}^m$ sur $\mathbb{R}$. La base (duale) de cet espace est constituée des **covecteurs** $\mathrm{d}^i$, projections de $x\in\mathbb{R}^m$ sur sa i-èm composante, soit $x^i$. Cette remarque est le fondement de l'approche directe proposée ici. Cela va s'éclaircir grâce à l'exemple qui suit.
"""

# ╔═╡ ed6b55b4-882a-475b-8829-a788fa3d6108
md"""
## Exemple

On va appliquer la formule $\mathrm{d}(\Phi\circ\varphi) = \sum_{j=1}^n \partial_j \Phi\  \mathrm{d}\varphi^{\,j}$ à l'exemple $z=x+\sin(xy)$

Les covecteurs $\mathrm{d}x, \mathrm{d}y$ ont pour composantes $(1, 0)$ et $(0, 1)$ dans le dual engendré par $\{\mathrm{d}x, \mathrm{d}y\}$.

La premiere étape est le calcul de 
```math
\mathrm{d}((x,y)\mapsto xy)=y\mathrm{d}x+x\mathrm{d}y
```
En utilisant la formule précédente et en terme de composantes on a :
```math
\mathrm{d}((x,y)\mapsto xy) = y\ (1, 0) + x\ (0, 1) = (y, x)
```

On poursuit de la même façon...

```math
\begin{align}
\mathrm{d}x &\leadsto (1, 0) \\
\mathrm{d}y &\leadsto (0, 1)\\
\mathrm{d}((x,y)\mapsto xy) &\leadsto y\ (1, 0)+x\ (0, 1) = (y, x) \\
\mathrm{d}(xy\mapsto\sin(xy)) &\leadsto \cos(xy)\ (y, x) \\
&= (y\cos(xy), x\cos(xy))\\
\mathrm{d}((x,\sin(xy))\mapsto x+\sin(xy)) &\leadsto 1\ (1\ 0) + 1\ (y\cos(xy)\ \ x\cos(xy)) \\
&= (1+y\cos(xy), x\cos(xy))\\
\end{align}
```
On obtient le résultat attendu :
```math
d(x+\sin(xy)) = (1+y\cos(xy))\mathrm{d}x + x\cos(xy)\mathrm{d}y
```

!!! note "gradient avec le mode forward ?!?"
    On vient le calculer le gradient avec l'approche forward. Cependant bien avoir à l'esprit que le calcul précédent cache 2 calculs simultanés : un par composante, ici 
	```math
	\begin{align}
	\mathrm{d}x &\leadsto (1\ 0) \\
	\mathrm{d}y &\leadsto (0\ 1)\\
	\end{align}
	```
	ce qui posserait un réel problème d'efficacité si l'on avait un grand nombre de variables indépendantes...

!!! remark "Et les \"Dual numbers\"..."
    L'approche précédente est simple et "directe". On peut recourir aux "dual 		numbers", pour faire exactement le même genre de calculs, mais la justification 	emmène beaucoup plus loin [Smooth infinitesimal analysis](https://en.wikipedia.org/wiki/Smooth_infinitesimal_analysis). L'idée est d'introduire un nombre $\epsilon$ tel que $\epsilon^2=0$, on a alors $f(a+\epsilon b)=f(a)+ϵf'(a)b$. C'est une réminiscence des infiniments petits de Leibniz où les termes d'ordre supérieur sont négligés $(dx)^2=0$. 
"""

# ╔═╡ 1989ee2c-3dfb-4fcd-862d-71c0f3a868d9
md"""
# Implémentation jouet

L'implémentation suit exactement le cheminement précédent. À tout scalaire on associe un covecteur. 

En termes d'optimisation, afin d'éviter des allocations dynamiques, on utilise un vecteur statique. Si le nombre de variables est trop important, on découpe ce vecteur en plusieurs morceaux et on relance les calculs plusieurs fois (ou en parallèle, car ce sont des calculs indépendants).
"""

# ╔═╡ 803d756b-2e35-42f3-ae9f-6798e4d82d27
md"""
## Nouveau type pour les scalaires
"""

# ╔═╡ 44263ffa-80c0-4657-a625-28b2c19b23a4
struct ADVar{N} <: Real
	v::Real
	d::SVector{N,Real} # Vecteur statique (dans la stack)

	function ADVar(N::Int, i::Int, v::Real) 
	  new{N}(v,SVector{N,Real}((k==i ? 1 : 0 for k in 1:N)))
	end

	function ADVar(v::Real,d::SVector{N,<:Real}) where N
	  new{N}(v,d)
	end
end

# ╔═╡ 4bd51561-d8eb-4ae8-9e23-e3b5dcd4aa8b
x, y = ADVar(2,1,2.0), ADVar(2,2,3.0)

# ╔═╡ b1c123a1-acb6-4420-9cbc-a147eceb6a5a
md"""
## Surcharge des opérateurs
"""

# ╔═╡ 84616b6d-f84a-49ef-8428-5c58ac04229a
md"""
### Mise en oeuvre pour $\sin()$
"""

# ╔═╡ f902c9a5-33c6-4eb1-87df-db12cd51a6fc
md"""
### Mise en oeuvre pour $+, \times, \dots$
"""

# ╔═╡ 464c0425-5db3-4264-9dd8-a88fbbb6d8c2
begin
	import Base: (+)
	
	function (+)(a::ADVar{N}, b::ADVar{N}) where N
    	v = a.v + b.v
		d = a.d + b.d
	
		ADVar(v, d)
	end

	# For nested demo
	#
	function (+)(a::Real, b::ADVar{N}) where N
    	v = a + b.v
		d = b.d
	
		ADVar(v, d)
	end
end

# ╔═╡ 91cdf562-71bd-421d-bd84-215d74b42074
begin
	import Base: (*)
	
	function (*)(a::ADVar{N}, b::ADVar{N}) where N
    	v = a.v * b.v
		d = b.v * a.d + a.v * b.d
	
		ADVar(v, d)
	end

	# For nested demo
	#
	function (*)(a::ADVar{N}, b::Real) where N
    	v = a.v * b
		d = a.d * b
	
		ADVar(v, d)
	end
end

# ╔═╡ 53ad6cbd-700f-43dd-8bd7-5a641c75f164
begin
	import Base: sin
	
	function sin(a::ADVar{N}) where N
    	v = sin(a.v)
		d = cos(a.v) * a.d
		
    	ADVar(v, d)      
	end

	# For nested demo
	#
	import Base: cos
	
	function cos(a::ADVar{N}) where N
    	v = cos(a.v)
		d = -sin(a.v) * a.d
		
    	ADVar(v, d)      
	end
end

# ╔═╡ 876461be-5460-4d5e-bf99-8428a0334aa9
md"""
## "Syntactic sugar"

Il s'agit d'encapsuler les procédures précédentes afin de rendre leur usage plus pratique. En particulier, il faut connaitre toutes les variables indépendantes dès le début pour être en mesure de les numéroter de $1$ à $N$.
"""

# ╔═╡ d773e0db-d442-42aa-83ad-5e3865e2dbf3
function gradient(f::Function, X::AbstractVector{<:Real})
	N=length(X) 
	X = map(((i,Xᵢ),)->ADVar(N,i,Xᵢ),enumerate(X))
	z=f(X)
	z.v, z.d
end

# ╔═╡ 59a17301-b59f-4602-b03d-7172ca123836
md"""
# Démonstration
"""

# ╔═╡ 0cb422c6-0338-48df-bd3b-b78956e6bc3e
md"""
## Exemple simple
"""

# ╔═╡ 654963db-6487-49b3-b405-95a31b8d8668
foo(X::AbstractVector{<:Real}) = X[1]+sin(X[1]*X[2])

# ╔═╡ 93d0ef14-2b6e-4810-b520-a7008f8a9696
gradient(foo,[2.0, 3.0])

# ╔═╡ bc7a2576-1c57-41d8-aae3-0ca13a1a786a
md"""
Vérification : si $z=x+\sin(xy)$, alors $\nabla z = (1+y\cos(xy),x\cos(xy))$. Si $x=2, y=3$ on obtient $z=1.72058$ et $\nabla z=(3.88051,1.92034)$
"""

# ╔═╡ 6628ad4b-8b16-4aa3-9dfe-50b2a510aecf
md"""
## Derivée directionnelle

Si l'on souhaite calculer une dérivée directionnelle, il est avantageux d'utiliser la formule :
```math
\nabla f\cdot d = \frac{d}{dt}(t\rightarrow f(X+t d))_{|t=0}
```
En effet, la seule variable à definir de type `ADVar` est $t$.

**Note :** dans l'approche utilisant une "cassette" d'enregistrement, cette astuce n'a pas d'intérêt, car il suffit de prendre $v=d$ à la place de $v=e_i$.
"""

# ╔═╡ ec41a684-2075-4a89-b9b1-aa361640f180
function directional_derivative(f::Function,X::AbstractVector,d::AbstractVector)
	t=ADVar(1,1,0)
	X+=t*d
	z=f(X)
end

# ╔═╡ 3a0b7b03-fd62-4818-844b-392b0a33bb11
directional_derivative(foo,[2.0,3.0],[4.0,5.0])

# ╔═╡ 6775c574-75ba-4f2b-a7f9-d5200133a1a7
md"""
## Dérivées d'ordre supérieur
"""

# ╔═╡ 73a4f2de-f4c9-4c5f-ac32-a54da0cb997e
function Hessian(f::Function, X::AbstractVector{<:Real})
    N=length(X) 
	X = map(((i,Xᵢ),)->ADVar(N,i,Xᵢ),enumerate(X))
    X = map(((i,Xᵢ),)->ADVar(N,i,Xᵢ),enumerate(X)) # nested 
	z=f(X)
	z.v.v, map(i->z.d[i].v,1:N), hcat(map(i->z.d[i].d,1:N)...)
end

# ╔═╡ 90182800-66f0-4270-a5af-29da9ba88711
H_z,H_∇,H_H = Hessian(foo,[2.0,3.0])

# ╔═╡ f9a906df-ddb0-4bc4-a39c-43d05fa4f499
H_z

# ╔═╡ a7d99d50-c5ae-4139-8db7-67ff14f25cd1
H_∇

# ╔═╡ b33ee2c2-01e0-4fad-8985-d96e5488c12c
H_H

# ╔═╡ 7ba4ae91-1a00-4b27-b43a-46f4cd4c4739
md"""
Vérification:
```math
H = 
\left(
\begin{array}{cc}
-y^2 \sin (x y) & \cos (x y)-x y \sin (x y) \\
\cos (x y)-x y \sin (x y) & -x^2 \sin (x y) \\
\end{array}
\right)
```
En $x=2, y=3$, on trouve :
```math
H = 
\left(
\begin{array}{cc}
2.51474 & 2.63666 \\
2.63666 & 1.11766 \\
\end{array}
\right)
```
"""

# ╔═╡ fa1e4f38-a3c4-4178-9769-ca7d5491f8dd
md"""
Cela marche, mais **beaucoup de calculs dont certains sont inutiles**. 

En effet, le théorème de Schwarz, nous dit que (si $f$ suffisament lisse) :
```math
\partial^2_{ij}f = \partial^2_{ji}f
```
Donc $1/2N(N-1)$ calculs sont redondants...

Il existe des approches, plus compliquées, qui évitent ces calculs inutiles :
1. R.M. Gower, M.P. Mello, A new framework for the computation of hessians, Optimization Methods and Software. 27 (2012) 251-273. https://doi.org/10.1080/10556788.2011.580098.

2. M. Wang, A. Gebremedhin, A. Pothen, Capitalizing on live variables : New algorithms for efficient hessian computation via automatic differentiation, Mathematical Programming Computation. 8 (2016) 393–433.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[compat]
PlutoUI = "~0.7.37"
StaticArrays = "~1.4.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "bf0a1121af131d9974241ba53f601211e9303a9e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.37"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74fb527333e72ada2dd9ef77d98e4991fb185f04"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─45ee7653-8bdf-472d-aa39-7a5be5f03dbe
# ╟─841ab7fb-aefb-434d-8d60-db585132dd19
# ╟─64776627-30c6-4ee1-8df3-6be2829bd9c2
# ╟─a3c76871-a056-40b9-8d01-151546b5a7e0
# ╟─2540688c-81be-4fe7-89b5-9976f44152e8
# ╟─ed6b55b4-882a-475b-8829-a788fa3d6108
# ╟─1989ee2c-3dfb-4fcd-862d-71c0f3a868d9
# ╟─803d756b-2e35-42f3-ae9f-6798e4d82d27
# ╠═15f16d2d-3295-48e2-99a7-8a21e8c28ca0
# ╠═44263ffa-80c0-4657-a625-28b2c19b23a4
# ╠═4bd51561-d8eb-4ae8-9e23-e3b5dcd4aa8b
# ╟─b1c123a1-acb6-4420-9cbc-a147eceb6a5a
# ╟─84616b6d-f84a-49ef-8428-5c58ac04229a
# ╠═53ad6cbd-700f-43dd-8bd7-5a641c75f164
# ╠═f902c9a5-33c6-4eb1-87df-db12cd51a6fc
# ╠═464c0425-5db3-4264-9dd8-a88fbbb6d8c2
# ╠═91cdf562-71bd-421d-bd84-215d74b42074
# ╟─876461be-5460-4d5e-bf99-8428a0334aa9
# ╠═d773e0db-d442-42aa-83ad-5e3865e2dbf3
# ╟─59a17301-b59f-4602-b03d-7172ca123836
# ╠═0cb422c6-0338-48df-bd3b-b78956e6bc3e
# ╠═654963db-6487-49b3-b405-95a31b8d8668
# ╠═93d0ef14-2b6e-4810-b520-a7008f8a9696
# ╟─bc7a2576-1c57-41d8-aae3-0ca13a1a786a
# ╟─6628ad4b-8b16-4aa3-9dfe-50b2a510aecf
# ╠═ec41a684-2075-4a89-b9b1-aa361640f180
# ╠═3a0b7b03-fd62-4818-844b-392b0a33bb11
# ╟─6775c574-75ba-4f2b-a7f9-d5200133a1a7
# ╠═73a4f2de-f4c9-4c5f-ac32-a54da0cb997e
# ╠═90182800-66f0-4270-a5af-29da9ba88711
# ╠═f9a906df-ddb0-4bc4-a39c-43d05fa4f499
# ╠═a7d99d50-c5ae-4139-8db7-67ff14f25cd1
# ╠═b33ee2c2-01e0-4fad-8985-d96e5488c12c
# ╟─7ba4ae91-1a00-4b27-b43a-46f4cd4c4739
# ╟─fa1e4f38-a3c4-4178-9769-ca7d5491f8dd
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
