### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ ce235e96-6bf8-4601-847b-808cd0d06261
begin
 	using CairoMakie
	using PlutoUI
	PlutoUI.TableOfContents(aside=true)
end

# ╔═╡ 7c079152-43a9-4290-81b4-04788a4eb3bc
using LinearAlgebra,BenchmarkTools

# ╔═╡ c827ef5e-c184-4e06-9d0f-fb4a22a71266
md"""
!!! exemple ""
    **Introduction à la Différentiation Automatique**

    vincent.picaud@cea.fr
"""

# ╔═╡ b8371073-56c9-4075-b16e-8f17b0b93d9b
md"""
# Introduction

La "différentiation automatique" est un **ensemble** de techniques ayant pour but le calcul efficace des dérivées.

La différentiation automatique **n'est pas** :
- la differentiation symbolique
- l'usage des différences finies
"""

# ╔═╡ cffc3edc-1806-444e-aecc-373e869dfddb
md"""
## Différentiation symbolique

Non recommandée en général:
- La coomplexité des expressions croit rapidement
  ```math
  \begin{align*}
  \left(\frac{f}{g}\right)' &= \frac{g(x) f'(x)-f(x) g'(x)}{g(x)^2} \\
  \left(\frac{f}{g}\right)'' &= \frac{g(x)^2 f''(x)-g(x) \left(2 f'(x) g'(x)+f(x)   g''(x)\right)+2 f(x)
   g'(x)^2}{g(x)^3}
  \end{align*}
  ```
  il faut ajouter la lenteur du calcul formel (manipulation d'arbres, heuristiques...)
- Problématique pour différentier un code général: boucles, branchements etc...
```julia
function f(x)
	if x<0
		return 0
	else
		return x*x
	end
end
```
"""

# ╔═╡ 89843eb3-9bac-4924-a98b-63b0c3d42697
md"""
## Différences finies

Deux problèmes principaux:
- stabilité numérique :
  ```math
  \lim_{h\rightarrow 0}\frac{f(x_0+h)-f(x_0)}{h}
  ```

- inefficace dans le cas multidimensionel $f(x_1,\dots,x_n)$ :
  ```math
  i=1\dots n,\ \lim_{h\rightarrow 0}\frac{f(x_0+h e_i)-f(x_0)}{h}
  ```

Par contre c'est un moyen commode lors de phase se test/debuggage
"""

# ╔═╡ 3c9ff7f2-6b4c-4d0a-9bc3-e719504b7f1e
md"""
### Stabilité numérique, illustration 

Il s'agit de calculer $\sin'(1)$ en utilisant
```math
\sin'(1) \approx \frac{\sin(1+h)-\sin(1)}{h}
```



!!! remark
    en pratique si l'on doit vraiment utiliser les différences finies, on prend      
    $h$ de l'ordre de $h\approx x_0 \sqrt{\epsilon}$. Pour les `Float64`, $\epsilon \approx       10^{-16}$, donc $h\approx 10^{-8}$.
"""

# ╔═╡ 579249f9-e7d8-4097-b358-7d410a1b0d87
md"""
## Differentiation automatique

Il y a principalement deux techniques pour mettre en ouvre la différentiation automatique:
- **surcharge des opérateurs**
- transformation du code source

Il y a également deux modes de fonctionnement:
- mode **forward** 
- mode **backward**

!!! remark 
    pour les dérivations **d'ordre supérieur**, il est possible d'alterner les
	modes **forward** et **backward**
"""

# ╔═╡ adfafc42-1dc5-4ba1-8922-bd300423540b
md"""
# Rappels de calcul différentiel
"""

# ╔═╡ a5046c2e-1d37-4a89-9c34-f0c961592f94
md"""
## Définition

!!! definition "Différentielle"
    Une fonction $f$ est differentiable (dérivée de Frechet) en $x$ s'il existe une application linéaire bornée, notée $df_{|x}$, telle que:
    ```math
	f(x+\delta x) = f(x) + df_{|x}\cdot \delta x + o(\|\delta x\|)
    ```
    (rappel: $f\in o(g) \Leftrightarrow \lim_{x\rightarrow 0}\frac{\|f(x)\|}{g(x)}=0$)

    Une fonction $f$ est differentiable si elle est différentiable en tout point $x$.

**Retour sur les notations:**

si $f:\mathbb{R}^m\to\mathbb{R}^n$, alors:
```math
\begin{eqnarray}
 df: \mathbb{R}^m & \to & \mathcal{L}(\mathbb{R}^m,\mathbb{R}^n)\\
     x & \mapsto & df_{|x} \\
\end{eqnarray}
```
et
```math
\begin{eqnarray}
 df_{|x}: \mathbb{R}^m & \to & \mathbb{R}^n\\
     \delta x & \mapsto & df_{|x}\cdot \delta x \\
\end{eqnarray}
```

### Exemple 1
On considère $g:x\mapsto \|x\|_2^2$. 

Comme
```math
\Large{
\underbrace{\|x+\delta x\|_2^2}_{g(x+\delta x)} =  \underbrace{\|x\|_2^2}_{g(x)} + \underbrace{2 \langle  x,\delta x \rangle}_{dg_{|x}\cdot \delta x}  + \underbrace{\|\delta x\|_2^2}_{o(\|\delta x\|)}}
```

La différentielle est donc:
```math
dg =  x\mapsto (\delta x \mapsto 2 \langle  x, \delta x \rangle) \in \left(\mathbb{R}^m \to \mathcal{L}(\mathbb{R}^m,\mathbb{R}) \right)
```
La différentielle en $x_0$ est:
```math
dg_{|x_0} =  \delta x \mapsto 2 \langle  x_0, \delta x \rangle \in \mathcal{L}(\mathbb{R}^m,\mathbb{R})
```
sont action sur le vecteur $\delta x$ est:
```math
dg_{|x_0}\cdot \delta x =  2 \langle  x_0, \delta x \rangle \in \mathbb{R}
```

### Exemple 2
On considère la fonction affine: $f:x\mapsto M\cdot x + y$. 

Comme
```math
\Large{
\underbrace{M\cdot (x+\delta x)-y}_{f(x+\delta x)} =  \underbrace{M\cdot x -y}_{f(x)} + \underbrace{M\cdot \delta x}_{df_{|x}\cdot \delta x}  + \underbrace{0}_{o(\|\delta x\|)}}
```
la differentielle est l'application constante: 
```math
x\mapsto (\delta x \mapsto M\cdot \delta x) 
```
"""

# ╔═╡ 88fc447e-007d-4eca-811f-badf024578cd
md"""
## Règle de la chaine

Il s'agit de trouver la différentielle des fonctions composées

!!! exemple "Règle de la chaine"
    ```math
    d(g \circ f)_{|x} \cdot \delta x = dg_{|f(x)}\circ df_{|x} \cdot \delta x 
    ```

### Exemple

On considère $x\mapsto h(x)=\| M\cdot x - y\|^2_2$. 

Comme $h=g\circ f$ avec $g=\|.\|_2^2$ et $f=M\cdot . -y$, alors:

```math
\Large{
 \begin{align*}
   d(x\mapsto h(x))_{|x_0} &= d(x \mapsto g(x))_{|Mx_0-y}\circ d(x \mapsto f(x))_{|x_0} \\
                                  &= (\delta f \mapsto  2 \langle  M\cdot x_0-y, \delta f \rangle)  \circ (\delta x \mapsto M \delta x) \\
			          &= \delta x \mapsto 2 \langle  M\cdot x_0-y, M \cdot \delta x \rangle \\
\end{align*}
}
```
La différentielle de $h$ est donc:
```math
dh = (x\mapsto (\delta x\mapsto 2 \langle  M\cdot x-y, M\cdot \delta x \rangle))
```
"""

# ╔═╡ f3309a59-d96e-42d3-be7a-3bb124619741
md"""
## Gradient

Contrairement à la différentielle, pour définir un gradient il faut un **produit scalaire** (espace de Hilbert). L'existence du gradient
découle du théorème de représentation de Ritz, qui démontre que toute toute forme linéaire $\mathcal{l}$ alors il existe un vecteur $v$ tel que $l(x)=\langle v,x \rangle$.

!!! definition "Vecteur gradient"
    Pour une fonction differentiable $f:\mathbb{R}^m\to\mathbb{R}$, la diffentielle  	en $x$, $df_{|x}$ est un élément de $\mathcal{L}(\mathbb{R}^m,\mathbb{R})$. C'est 	donc une forme linéaire. En utilisant Ritz, il existe donc un vector $v_x$ tel 		que:
	```math
	df_{|x}\cdot \delta x  =\langle v_x , \delta x \rangle
	```
    Le vecteur $v_x$ est le **gradient** de $f$ en $x$, il est souvent noté $\nabla 		f_{|x}$

### Exemple

Si l'on revient à l'exemple: $x\mapsto h(x)=\| M\cdot x - y\|^2_2$, qui est une application de $\mathbb{R}^m$ dans $\mathbb{R}$, nous avions montré que:

```math
dh_{|x} = (\delta x\mapsto 2 \langle  M\cdot x-y, M\cdot \delta x \rangle) \in \mathcal{L}(\mathbb{R}^m,\mathbb{R})
```
Cette forme linéaire peut s'écrire:
```math
dh_{|x} = \delta x\mapsto \langle 2 M^t(M\cdot x-y), \delta x \rangle
```
le gradient est donc:
```math
\nabla h_{|x} = 2 M^t(M\cdot x-y)
```
"""

# ╔═╡ c08c6755-342b-4c0e-a00e-c15f84c1bf83
md"""
## Matrice Jacobienne

Soit $f:\mathbb{R}^m\to\mathbb{R}^n$, la matrice Jacobienne est la représentation matricielle de l'application linéaire $\delta x\mapsto df_{|x}\cdot \delta x$

```math
J_{|x}=
\left(
\begin{array}{ccc}
\partial_1 f^1_{|x} & \dots & \partial_m f^1_{|x} \\
\vdots & \ddots & \vdots \\
\partial_1 f^n_{|x} & \dots & \partial_m f^n_{|x} \\
\end{array}
\right)
```

!!! example "A retenir..."
    - Pour obtenir la colonne $j$ de $J$ il suffit de calculer $J\cdot e_j$. Ceci 		permet de calculer en une seul produit la dérivée de toutes les quantités par 		rapport à une seule variable $x_j$:
	```math
	\left(
	\begin{array}{c}
	\partial_j f^1_{|x} \\
	\vdots \\
	\partial_j f^n_{|x}  \\
	\end{array}
	\right)
	```
	- pour obtenir la ligne $i$ de $J$ il suffit de calculer $e_i^t J^t = J^t\cdot e_i$. Ceci permet de calculer la dérivée d'une seule grandeur par rapport à toutes les variables (cad une information de "type" gradient)
	```math
	\left(
	\begin{array}{ccc}
	\partial_1 f^i_{|x} & \dots & \partial_m f^i_{|x} \\
	\end{array}
	\right)
    ```

Cette remarque toute simple permettra de cerner l'interet des modes **forward** et **backward**.

### Exemple

Toujours en revenant à l'exemple: $x\mapsto h(x)=\| M\cdot x - y\|^2_2$, les composantes de la différentielle (matrice Jacobienne) et le gradient peuvent s'écrire:

```math
df_{|x}(\delta x) = \underbrace{\left( \partial_1 f_{|x}, \partial_2 f_{|x}, \dots, \partial_m f_{|x} \right)}_J \cdot
\left( \begin{array}{c} \delta x_1 \\ \delta x_2 \\ \vdots \\ \delta x_m \end{array}\right) =
\langle \left( \begin{array}{c} \partial_1 f_{|x} \\ \partial_2 f_{|x} \\ \vdots \\ \partial_m f_{|x} \end{array}\right), \left( \begin{array}{c} \delta x_1 \\ \delta x_2 \\ \vdots \\ \delta x_m \end{array}\right) \rangle =
\langle \nabla f_{|x}, \delta x \rangle 
```
En particulier pour calculer le gradient, il suffit d'évaluer $J^t\cdot \mathbf{1}$

!!! note "Au risque d'insister..."
    Le gradient est un **vecteur**, mais c'est une **ligne** de la Jacobienne.

"""

# ╔═╡ 286f2fe1-443e-430b-8549-fb65c2bc224d
md"""
# Différentiation automatique
"""

# ╔═╡ f2f50bb5-a42d-4701-ae29-d9e2539e5136
md"""
## Idée de base

On souhaite calculer le gradient de $f:(x_1,x_2)\mapsto x_1+\sin(x_1x_2)$.

L'idée est de:
- décomposer le cacul en étapes élémentaires
- utiliser la règle de la chaine

Ici cela consiste a introduire les fonctions:
```math
\Phi_1(x_1,x_2)=
\left(
\begin{array}{c}
x_1 \\
x_2 \\
x_1 x_2
\end{array}
\right)
,\ 
\Phi_2(x_1,x_2,x_3)=
\left(
\begin{array}{c}
x_1 \\
x_2 \\
x_3 \\
\sin(x_3)
\end{array}
\right)
,\ 
\Phi_3(x_1,x_2,x_3,x_4)=
\left(
\begin{array}{c}
x_1 \\
x_2 \\
x_3 \\
x_4 \\
x_3 + x_4
\end{array}
\right)
```
de telle sorte que:
```math
\left(
\begin{array}{c}
x_1 \\
x_2 \\
x_1 x_2 \\
\sin(x_1 x_2) \\
\underbrace{x1 + \sin(x_1 x_2)}_{f(x_1,x_2)} \\
\end{array}
\right) = F(x_1,x_2) = \Phi_3 \circ \Phi_2 \circ \Phi_1(x_1,x_2)
```

La différentielle $dF$ est donnée par:
```math
dF_{|(x_1,x_2)} = (d\Phi_3)_{|\Phi_2(x_1,x_2,x_3)} \circ (d\Phi_2)_{|\Phi_1(x_1,x_2)} \circ (d\Phi_1)_{|(x_1,x_2)}
```
soit plus explicitement:
```math
dF_{|(x_1,x_2)} = 
\left(
\begin{array}{cccc}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1 \\
\hline
1 & 0 & 0 & 1 
\end{array}
\right)
\left(
\begin{array}{ccc}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1 \\
\hline
0 & 0 & \cos(x_3) 
\end{array}
\right)
\left(
\begin{array}{cc}
1 & 0 \\
0 & 1 \\
\hline
x_2 & x_1 
\end{array}
\right)
```
Tous calculs faits nous trouvons:
```math
dF_{|(x_1,x_2)} = 
\left(
\begin{array}{cc}
 1 & 0 \\
 0 & 1 \\
 x_1 & x_2 \\
 \cos \left(x_1 x_2\right) x_1 & \cos \left(x_1 x_2\right) x_2 \\
 \cos \left(x_1 x_2\right) x_1+1 & \cos \left(x_1 x_2\right) x_2 \\
\end{array}
\right)
```
Le gradient $\nabla f_{|(x_1,x_2)}$, se situe à la derniere ligne. Soit
```math
\nabla f_{|(x_1,x_2)} = 
\left(
\begin{array}{c}
\cos \left(x_1 x_2\right) x_1+1 \\
\cos \left(x_1 x_2\right) x_2
\end{array}
\right)
```
!!! note "Attention"
    Il est fondamental de comprendre pourquoi...

Pour calculer ce vecteur nous avons vu qu'il faut calculer $dF^t_{|(x_1,x_2)}\cdot e_5$. Ceci se traduit par:
```math
\begin{align}
dF^t_{|(x_1,x_2)}\cdot e_5 &= ((d\Phi_3)_{|\Phi_2(x_1,x_2,x_3)} \circ (d\Phi_2)_{|\Phi_1(x_1,x_2)} \circ (d\Phi_1)_{|(x_1,x_2)})^t \cdot e_5 \\
&= (d\Phi_1)^t_{|(x_1,x_2)} \circ (d\Phi_2)^t_{|\Phi_2(x_1,x_2)} \circ(d\Phi_3)^t_{|\Phi_2(x_1,x_2,x_3)}  \cdot e_5
\end{align}
```
En explicitant les calculs:
```math
\begin{align}
dF^t_{|(x_1,x_2)}\cdot e_5 &= 
\left(
\begin{array}{cc|c}
 1 & 0 & x_1 \\
 0 & 1 & x_2 \\
\end{array}
\right)
\left(
\begin{array}{ccc|c}
 1 & 0 & 0 & 0 \\
 0 & 1 & 0 & 0 \\
 0 & 0 & 1 & \cos \left(x_1 x_2\right) \\
\end{array}
\right)
\left(
\begin{array}{cccc|c}
 1 & 0 & 0 & 0 & 1 \\
 0 & 1 & 0 & 0 & 0 \\
 0 & 0 & 1 & 0 & 0 \\
 0 & 0 & 0 & 1 & 1 \\
\end{array}
\right)
\left(
\begin{array}{c}
 0 \\
 0 \\
 0 \\
 0 \\
 1 \\
\end{array}
\right) \\
&= 
\left(
\begin{array}{cc|c}
 1 & 0 & x_1 \\
 0 & 1 & x_2 \\
\end{array}
\right)
\left(
\begin{array}{ccc|c}
 1 & 0 & 0 & 0 \\
 0 & 1 & 0 & 0 \\
 0 & 0 & 1 & \cos \left(x_1 x_2\right) \\
\end{array}
\right)
\left(
\begin{array}{c}
 1 \\
 0 \\
 0 \\
 1 \\
\end{array}
\right) \\
&= 
\left(
\begin{array}{cc|c}
 1 & 0 & x_1 \\
 0 & 1 & x_2 \\
\end{array}
\right)
\left(
\begin{array}{c}
 1 \\
 0 \\
 \cos \left(x_1 x_2\right) \\
\end{array}
\right) \\
&= 
\left(
\begin{array}{c}
 \cos \left(x_1 x_2\right) x_1+1 \\
 \cos \left(x_1 x_2\right) x_2 \\
\end{array}
\right) \\
&= \nabla f_{|(x_1,x_2)}
\end{align}
```

Le calcul précédent, de type $J^t\cdot e_i$,  qui calcule le gradient est le mode **backward** de la différentiation automatique. Le mode **forward** consiste à utiliser l'approche direct $J\cdot e_j$.
"""

# ╔═╡ c97b2d9c-fd4f-4502-8482-aa6db5b05325
md"""
## Les modes forward et backward

Si le calcul de $F$ se décompose en une suite d'étapes élémentaires:
```math
F = \Phi_N \circ \dots \circ \Phi_2 \circ \Phi_1 
```
alors la règle de la chaine s'éxprime par:
```math
dF = d\Phi_N \circ \dots \circ d\Phi_2 \circ d\Phi_1 
```

### Mode forward
Le mode **forward** permet de calculer la matrice Jacobienne colonne par colonne en suivant le sens direct des calculs de $F$:
```math
\Large{
\left(
\begin{array}{c}
\partial_j F^1_{|x} \\
\vdots \\
\partial_j F^n_{|x}  \\
\end{array}
\right) = d\Phi_N \circ \dots \circ \underbrace{ d\Phi_3 \circ \underbrace{d\Phi_2 \circ \underbrace{d\Phi_1.e_j}_{v=d\Phi_1.e_j}}_{v=d\Phi_2.v}}_{v=d\Phi_3.v}
}
```

Avec le mode forward, il est également possible de calculer une **dérivée directionnelle** en replacent $e_j$ par $d$:
```math
d(t\rightarrow F(x+t\,d))_{|t=0}=dF_{|x}\cdot d
```

### Mode backward
Le mode **backward** permet de calculer la matrice Jacobienne ligne and ligne en suivant le sens inverse des calculs de $F$. Ceci nécessite donc un **stockage**:
```math
\Large{
\left(
\begin{array}{c}
\partial_1 F^i_{|x} \\
\vdots \\
\partial_m F^i_{|x}  \\
\end{array}
\right) =  d\Phi^t_1 \circ \dots \circ \underbrace{ d\Phi^t_{N-2} \circ \underbrace{d\Phi^t_{N-1} \circ \underbrace{d\Phi^t_N.e_i}_{v=d\Phi^t_N.e_i}}_{v=d\Phi^t_{N-1}.v}}_{v=d\Phi^t_{N-2}.v}
}
```

## Cas d'usages

L'idée principale est que l'on peut calculer une matrice Jacobienne de dimension $n\times m$, ligne par ligne ou colonne par colonne.

- Le mode **forward** est utilisé pour construire $J$ colonne par colonne, ce qui est avantageux quand $m<n$. 
- Le mode **backward** est utilisé pour construire $J$ ligne par ligne, ce qui est avantageux quand $n<m$. Un exemple classique est le cas du grradient.
"""

# ╔═╡ f72ec959-a3aa-420a-bf12-bc61f63a6d65
md"""# Premiere implémentation

Pour l'instant toutes les implémentations présentées sont des implementations "jouet" dont le but est de comprendre la démarche sans chercher ni la perfomance, ni la généricité.
"""


# ╔═╡ 176e0e85-dacc-4c06-ad1a-30adf486ce3d
md"""
## Enregistrement dans une "cassette" (tape)

Si l'on suppose que les transformations élémentaires $\Phi_k$ sont de la forme:
```math
\Phi_k(x_1,\dots,x_{k-1})=
\left(
\begin{array}{c}
x_1 \\
\vdots \\
x_{k-1} \\
x_k=\phi(x_1,\dots,x_{k-1})
\end{array}
\right)
```
alors $d\Phi$ est une matrice creuse de la forme:
```math
d\Phi_k=
\left(
\begin{array}{cccc}
1 & 0 & \dots & 0 \\
0 & 1 & \ddots & \vdots \\
\vdots & \ddots & 1 & 0 \\
0 & \dots & 0 & 1  \\
\hline
\partial_1 \phi_k & \dots & \dots & \partial_{k-1} \phi_{k-1}
\end{array}
\right)
```

L'ensemble des matrices $d\Phi_k, k=1,\dots,N$ est stocké dans un vecteur de vecteurs: la partie supérieure (matrice identité) est implicite et seul le vecteur creux $(\partial_1 \phi_k, \dots, \partial_k \phi_k)$ est stocké à la ligne $k$. 

Cette structure définie la "cassette" où sont enregistrées les informations.
"""

# ╔═╡ 18b26595-be3e-4f29-bcc0-96bdf3fdb6db
# Stocke ∂ᵢϕ
struct ∂iϕ
    i::Int
    ∂::Float64
end

# ╔═╡ d7a484d1-f54c-4ef2-8cde-d48c1d1e3be2
# Global tape
# ( remarque: il est possible d'être plus efficace en utilisant un stockage de type 
#   compressed row stoage (CRS). On ne le fait pas ici dans un but de simplicité )
const tape = Vector{Vector{∂iϕ}}()

# ╔═╡ a48c6a8f-e0c5-4679-9c85-727861ee9f24
md"""
Pour implémenter le mode **direct**, il faut ensuite être en mesure de calculer
```math
\begin{align}
\delta &= d\Phi_1 \cdot \delta \\
\vdots &= \ \ \ \ \vdots \\
\delta &= d\Phi_N \cdot \delta \\
\end{align}
```

La structure des matrices $d\Phi_k$ permet de réaliser cette opération "sur-place" :
"""

# ╔═╡ 65177381-acd6-4f11-ae47-bd942903e608
md"""
Pour implémenter le mode **backward**, il faut ensuite être en mesure de calculer
```math
\begin{align}
\delta &= d\Phi^t_N \cdot \delta \\
\vdots &= \ \ \ \ \vdots \\
\delta &= d\Phi^t_1 \cdot \delta \\
\end{align}
```

La structure des matrices $d\Phi^t_k$,

```math
d\Phi^t_k=
\left(
\begin{array}{cccc|c}
1 & 0 & \dots & 0 & \partial_1 \phi_k \\
0 & 1 & \ddots & \vdots  & \vdots \\
\vdots & \ddots & 1 & 0  & \vdots \\
0 & \dots & 0 & 1  & \partial_{k-1} \phi_k  \\
\end{array}
\right)
```

permet de réaliser cette opération "sur-place" :
"""

# ╔═╡ 8132cc8c-2ca1-49af-8322-b53660d8cbb5
md"""
## Surcharge des opérateurs

Pour remplir la "cassette" d'enregistrement on définit un nouveau type `AFloat64`.

!!! exemple "Remarque"
    Lorsque l'on ajoute une variable, cela se traduit par une ligne vide. Il faut prendre un peu de temps pour réaliser que cela se passe bien. Ceci tient au fait, que la partie matrice identité est implicite.

Le code Julia, ainsi que le contructeur associé est ci-dessous:
"""

# ╔═╡ 3d469d0a-723d-449a-bce2-80833f26b8c8
# A new type of scalar 
struct AFloat64 <: Real
    i::Int
    v::Float64

	function AFloat64(v::Real)
	    # When a variable is created one must add an empty row 
    	# (the "identity" part).
    	push!(tape, [])
    	i = length(tape)

    	new(i, Float64(v))
	end

	AFloat64(i::Int, v::Float64) = new(i,v)
end

# ╔═╡ 98053b26-9545-4f75-a49e-84be153b9af1
md"""
**Illustration:** on déclare deux variables `x=2` et `y=3` et l'on affiche la "cassette", qui ne contient pour l'instant que deux lignes vides.
"""

# ╔═╡ 6770e01c-ca6d-45e3-a3a8-48797e338679
let
	resize!(tape,0)
	x=AFloat64(2.0)
	y=AFloat64(3.0)
	
	x, y, tape
end

# ╔═╡ 76728ec4-b159-4aed-a9d2-c597143546cf
md"""
### Mise en oeuvre pour $\sin()$

Comme
```math
d(\sin(x_i))=\cos(x_i)dx_i
```
la surcharge de la fonction sinus s'écrit:
"""

# ╔═╡ 78a9374b-5543-4bb1-8e0a-ff5d58667b71
begin
	import Base: sin
	
	function sin(a::AFloat64)
    	d = [∂iϕ(a.i, cos(a.v))] # <- dϕ = [0,…,0, ∂ₖsin(xₖ),0,…0]
    	push!(tape, d) # <- enregistrement de la ligne dans la "cassette"

    	i = length(tape) # création de la nouvelle variable xₖ₊₁ = sin(xₖ) 
    	v = sin(a.v)
    	AFloat64(i, v)      
	end
end

# ╔═╡ 8efe8f97-f17a-4163-a89b-934c423999d8
md"""
### Mise en oeuvre pour $+, \times, \dots$
"""

# ╔═╡ a405c213-f0ac-4df6-818a-abd1b52000fb
md"""
Comme $d(x_i + x_j) = dx_i + dx_j$, alors
"""

# ╔═╡ 5b219ef8-0cbe-40a6-beed-335735a07d1b
begin
	import Base: (+)
	
	function (+)(a::AFloat64, b::AFloat64)
    	d = [∂iϕ(a.i, 1.0), ∂iϕ(b.i, 1.0)]
    	push!(tape, d)

    	i = length(tape)
    	v = a.v + b.v
	    AFloat64(i, v)
	end
end

# ╔═╡ a76bc54a-7ed3-49fe-b8b5-82c8ac33b1e7
let
	f(x)=sin(x)
	err(f′,x0)=abs(f′-cos(x0))
	h = 10 .^ range(-17, -1, length=500)
	x0 = 1
	y = map(h->err( (f(x0+h)-f(x0))/h, x0), h)
	
	fig = lines(h, y,
		axis=(xscale=log10, yscale = log10,
		    xlabel="h",ylabel="ϵ",
		    title=L"\epsilon=\mid\frac{\sin(1+h)-\sin(1)}{h}-\cos(1)\mid"))
	
    fig
end

# ╔═╡ 63f891da-12ed-4813-be9e-b08cd6c274cb
md"""
Comme $d(x_i \times x_j)= x_jdx_i + x_idx_j$, alors:
"""

# ╔═╡ 89f2fecf-d07f-44e0-b39d-fbb6c509fbf6
begin
	import Base: (*)
	
	function (*)(a::AFloat64, b::AFloat64)
    	d = [∂iϕ(b.i, a.v), ∂iϕ(a.i, b.v)]
	    push!(tape, d)

	    i = length(tape)
	    v = a.v * b.v
	    AFloat64(i, v)
	end
end

# ╔═╡ e9d11e1f-501c-49f3-bfca-2a9d5fc236b1
# In-place computation of ...∘dΦ₂∘dΦ₁⋅δ
function forward!(δ::Vector)
    n = length(tape)

    for i = 1:n
        for ∂iϕ in tape[i]
            δ[i] += ∂iϕ.∂ * δ[∂iϕ.i]
        end
    end

    δ
end

# ╔═╡ 55a5e351-e1eb-4600-b032-bc144c75b037
# In-place computation of ...Φn-1^t∘Φn^t.δ
function backward!(δ::Vector)
    n = length(tape)

    for i = n:-1:1
        s = δ[i]
        for ∂iϕ in tape[i]
            δ[∂iϕ.i] += ∂iϕ.∂ * s
        end
    end

    δ
end

# ╔═╡ 295942c8-0ef5-4c1a-bf99-561e4b4597c1
md"""
!!! exemple "Remarque"
    Pour être exaustif, il faudrait également définir les opérations "mixtes":
    - `Float64 + AFloat64`, `AFloat64 + Float64`
	- `Float64 × AFloat64`, `AFloat64 × Float64`
    - ...
    Cela ne réprésente pas de difficultés particulières.
"""

# ╔═╡ d2011efa-d26c-4ed1-9688-bca52334c3bf
md"""
## Illustration 
"""

# ╔═╡ 423af97b-79ff-4fc8-8b67-26901b33ac50
begin
	resize!(tape,0)
	
	x=AFloat64(2.0)
	y=AFloat64(3.0)

	z = x + sin(x*y)	
end

# ╔═╡ a7c2573a-19a2-4d4b-91a1-116e6d1e7c6a
md"""
La "cassette" a enregistré, 
1. declaration de `x`
2. declaration de `y`
3. declaration de `x*y`
4. declaration de `sin(x*y)`
5. declaration de `x+sin(x*y)`
"""

# ╔═╡ 99173920-cb95-4a02-9d95-41cf04459183
tape

# ╔═╡ b94d5517-2b3a-4014-87be-8c258c104e02
md"""
Pour calculer le gradient de $z(x,y)$ on crée un vecteur $e_i$ où $i$ est l'indice de $z$ et on utilise le mode backward.
"""

# ╔═╡ c0f4fe11-c149-43e8-afa8-962a035b8b15
let
	v=zeros(length(tape))
	v[z.i]=1

	z, backward!(v)
end

# ╔═╡ d6d76209-a5ba-4f0b-ab30-7632ed091ce5
md"""
Vérification: si $z=x+\sin(xy)$, alors $\nabla z = (1+y\cos(xy),x\cos(xy))$. Si $x=2, y=3$ on obtient $z=1.72058$ et $\nabla z=(3.88051,1.92034)$
"""

# ╔═╡ 5eccbecc-519f-4a3c-80d6-cd77e118d024
md"""
## Encapsulation des routines & "Syntactic sugar"

Il s'agit d'encapsuler les procédures précédentes afin de rendre leur usage plus pratique.

On procède en deux temps:
1. une première encapsulation qui cache la création du vecteur $e_i$ et l'extraction de la valeur des dérivées
2. une seconde encapsulation qui cache toute référence au type `AFloat64`
"""

# ╔═╡ 5691d49f-b2ec-4c24-ba5e-a1eee8f59476
md"""
### Premier temps
"""

# ╔═╡ a2e5fa68-1657-448f-95f1-60d4118f4def
# Calcule ∇z
function gradient(z::AFloat64,wrt::AbstractVector{AFloat64})
	v=zeros(length(tape)) # Création de eᵢ
	v[z.i]=1
	
	backward!(v) # Mode backward

	map(wrt->v[wrt.i],wrt) # Extraction des dérivées
end

# ╔═╡ 84159453-6daa-4304-a9f2-b370391059c7
md"""
Calcul du gradient $\nabla z$ par rapport à $x, y$
"""

# ╔═╡ c9c4f960-2d1b-4131-a3c3-bc722c365aa9
# Calcule ∂ⱼz
function derivatives_wrt(z::AbstractVector{AFloat64},wrt::AFloat64)
   v=zeros(length(tape)) # Création de eᵢ
   v[wrt.i]=1

   forward!(v)	# Mode forward
	
   map(z->v[z.i],z) # Extraction des dérivées
end   

# ╔═╡ 43a0026a-77d7-437d-b154-9baa1bcf0f80
md"""
Calcul des derivées de toutes les grandeurs par rapport à x
- $\partial_x x = 1$
- $\partial_x y = 0$
- $\partial_x z = \partial_x (x+\sin(xy)) = 1+y\cos(xy) = 3.88051$
"""

# ╔═╡ 83745652-35b3-42c8-bc77-61485fb2250f
derivatives_wrt([x,y,z],x)

# ╔═╡ 0d335e76-d3c6-4591-ac43-dae5188e0400
md"""
### Second temps
"""

# ╔═╡ 2952fb74-a351-493f-841f-466392382d89
md"""
Une fonction utilisée pour la démonstation
"""

# ╔═╡ dccc8320-f60f-4849-8ef1-ec3a1990e629
function foo(X::AbstractVector{<:Real})
  X[1]+sin(X[1]*X[2])
end

# ╔═╡ 366f8aa3-696c-4068-804f-28ebbd8b1bbb
md"""
Encapsulation "user friendly" du calcul du gradient. 

Le type `AFloat64` n'apparait plus.
"""

# ╔═╡ 3b072230-fdc4-4153-bf55-30e6dc6e63f5
function gradient(f::Function,X::AbstractVector{<:AbstractFloat})
	resize!(tape,0) # Efface l'ancienne "cassette"
	X = map(Xᵢ->AFloat64(Xᵢ),X) # Transforme les Float64 en AFloat64
	z = f(X) # Execute f en utilisant des scalaires de type AFloat64
	gradient(z,X) # Retourne le gradient
end

# ╔═╡ 6ef4ce4f-8ca5-4c76-9c13-a56fc6099773
gradient(z,[x,y])

# ╔═╡ 8d789020-4c1e-4802-a220-2551af74b373
md"""
**Illustration:** 
l'avantage ici est que toute référence au type AFloat64 est cachée. L'utilisateur ne sait "même pas" que l'on calcule le gradient en utilisant la différentiation automatique.
"""

# ╔═╡ 62d2e70b-c99b-43be-9f39-f67377d1855c
gradient(foo,[2.0,3.0])

# ╔═╡ 1ae1300f-e194-4f13-9dc5-6c94d15fe9e0
md"""
# Avantages/inconvenients

## Avantages

- Implementation facile & unifiée des modes forward et backward
- Il est possible de faire des "cassettes" imbriquées pour calculer les dérivations d'ordre supérieur.

Exemples de librairies:
- [Mission Impossible C++](https://github.com/vincent-picaud/MissionImpossible)
- [Stan](https://mc-stan.org/users/interfaces/math.html)
- [Adept C++](https://github.com/rjhogan/Adept)

## Inconvénients
- pour le mode forward on peut se passer d'enregistrer les opérations dans une "cassette", ce qui peut être plus efficace.
  - [autodiff C++](https://github.com/autodiff/autodiff)
  - [ForwardDiff.jl Julia](https://github.com/JuliaDiff/ForwardDiff.jl)
- travail au niveau des scalaires, ce qui bloque en particulier l'usage de librairies telles que BLAS ou LAPACK 

### 'dot' exemple
Pour l'algèbre linéaire dense, les librairies BLAS, LAPACK sont utilisées.

Elles définissent des fonctions optimisée pour les opérations courantes, telle que le calcul de $\langle x,y \rangle$:

```C
double cblas_ddot (const int n, 
                   const double *x, const int incx, 
                   const double *y, const int incy);
```

!!! note "Attention" 
    La source du problème est que si l'on travaille au niveaux "scalaire", les operations matricielles ne bénificient plus de BLAS & LAPACK. Ceci conduit à une fort perte de performance.
"""

# ╔═╡ 2628ace2-b983-4c97-8598-e118675df304
bench_blas_dot = let 
	x = Float64[2.0,3.0]

	@benchmark dot($x,$x)
end

# ╔═╡ 62dc62c3-ee4b-4148-aa1a-6d048fbec40d
bench_autodiff_dot = let 
	x = AFloat64[2.0,3.0]

	@benchmark dot($x,$x)
end

# ╔═╡ 2bd03e6a-1038-43bd-9065-17938b853374
md"""
On trouve ici un facteur $(minimum(bench_autodiff_dot.times)/minimum(bench_blas_dot.times)) sur les temps d'exécution.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
BenchmarkTools = "~1.3.1"
CairoMakie = "~0.7.4"
PlutoUI = "~0.7.35"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "745233d77146ad221629590b6d82fe7f1ddb478f"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "4.0.3"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA", "StaticArrays"]
git-tree-sha1 = "aedc7c910713eb616391cf95218277b714a7913f"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.7.4"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "3f1f500312161f1ae067abe07d13b40f78f32e07"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.8"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "9d3c0c762d4666db9187f363a76b47f7346e673b"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.49"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "84f04fe68a3176a583b864e492578b9466d87f1e"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "d7ab55febfd0907b285fbf8dc0c73c0825d9d6aa"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.3.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "463cb335fa22c4ebacfd1faba5fde14edb80d96c"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.5"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "80ced645013a5dbdc52cf70329399c35ce007fae"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.13.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "4c7d3757f3ecbcb9055870351078552b7d1dbd2d"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "StaticArrays"]
git-tree-sha1 = "770050893e7bc8a34915b4b9298604a3236de834"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.9.5"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "1c5a84319923bea76fa145d49e93aa4394c73fc2"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "169c3dc5acae08835a573a8a3e25c62f689f8b5c"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.6.5"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "65e4589030ef3c44d3b90bdc5aac462b4bb05567"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.8"

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

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "9a5c62f231e5bba35695a20988fc7cd6de7eeb5a"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.3"

[[deps.ImageIO]]
deps = ["FileIO", "JpegTurbo", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "464bdef044df52e6436f8c018bea2d48c40bb27b"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.1"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b15fc0a95c564ca2e0a7ae12c1f095ca848ceb31"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.5"

[[deps.IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "591e8dc09ad18386189610acafb970032c519707"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.3"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

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

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Serialization", "Showoff", "SignedDistanceFields", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "UnicodeFun"]
git-tree-sha1 = "cd0fd02ab0d129f03515b7b68ca77fb670ef2e61"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.16.5"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "c5fb1bfac781db766f9e4aef96adc19a729bc9b2"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.2.1"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test"]
git-tree-sha1 = "70e733037bbf02d691e78f95171a1fa08cdc6332"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.2.1"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "7e2166042d1698b6072352c74cfd1fca2a968253"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.6"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "eb4dbb8139f6125471aa3da98fb70f02dc58e49c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.14"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "1155f6f937fa2b94104162f01fa400e192e4272f"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.4.2"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a121dfbba67c94a5bec9dde613c3d0cbcf3a12b"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.3+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "13468f237353112a01b2d6b32f3d0f80219944aa"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "85bf3e4bd279e405f91489ce518dedb1e32119cb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.35"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "de893592a221142f3db370f48290e3a2ef39998f"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "afadeba63d90ff223a6a48d2009434ecee2ec9e8"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.1"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "01d341f502250e81f6fec0afe662aa861392a3aa"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SIMD]]
git-tree-sha1 = "39e3df417a0dd0c4e1f89891a281f82f5373ea3b"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.0"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "9cc2955f2a254b18be655a4ee70bc4031b2b189e"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "00b725fffc9a7e9aac8850e4ed75b4c1acbe8cd2"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.5.5"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74fb527333e72ada2dd9ef77d98e4991fb185f04"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "25405d7016a47cf2bd6cd91e66f4de437fd54a07"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.16"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "991d34bbff0d9125d93ba15887d6594e8e84b305"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.5.3"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78736dab31ae7a53540a6b752efc61f77b304c5b"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.8.6+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╟─ce235e96-6bf8-4601-847b-808cd0d06261
# ╟─c827ef5e-c184-4e06-9d0f-fb4a22a71266
# ╟─b8371073-56c9-4075-b16e-8f17b0b93d9b
# ╟─cffc3edc-1806-444e-aecc-373e869dfddb
# ╟─89843eb3-9bac-4924-a98b-63b0c3d42697
# ╟─3c9ff7f2-6b4c-4d0a-9bc3-e719504b7f1e
# ╟─a76bc54a-7ed3-49fe-b8b5-82c8ac33b1e7
# ╟─579249f9-e7d8-4097-b358-7d410a1b0d87
# ╟─adfafc42-1dc5-4ba1-8922-bd300423540b
# ╟─a5046c2e-1d37-4a89-9c34-f0c961592f94
# ╟─88fc447e-007d-4eca-811f-badf024578cd
# ╟─f3309a59-d96e-42d3-be7a-3bb124619741
# ╟─c08c6755-342b-4c0e-a00e-c15f84c1bf83
# ╟─286f2fe1-443e-430b-8549-fb65c2bc224d
# ╟─f2f50bb5-a42d-4701-ae29-d9e2539e5136
# ╟─c97b2d9c-fd4f-4502-8482-aa6db5b05325
# ╟─f72ec959-a3aa-420a-bf12-bc61f63a6d65
# ╟─176e0e85-dacc-4c06-ad1a-30adf486ce3d
# ╠═18b26595-be3e-4f29-bcc0-96bdf3fdb6db
# ╠═d7a484d1-f54c-4ef2-8cde-d48c1d1e3be2
# ╟─a48c6a8f-e0c5-4679-9c85-727861ee9f24
# ╠═e9d11e1f-501c-49f3-bfca-2a9d5fc236b1
# ╟─65177381-acd6-4f11-ae47-bd942903e608
# ╠═55a5e351-e1eb-4600-b032-bc144c75b037
# ╟─8132cc8c-2ca1-49af-8322-b53660d8cbb5
# ╠═3d469d0a-723d-449a-bce2-80833f26b8c8
# ╟─98053b26-9545-4f75-a49e-84be153b9af1
# ╠═6770e01c-ca6d-45e3-a3a8-48797e338679
# ╟─76728ec4-b159-4aed-a9d2-c597143546cf
# ╠═78a9374b-5543-4bb1-8e0a-ff5d58667b71
# ╟─8efe8f97-f17a-4163-a89b-934c423999d8
# ╟─a405c213-f0ac-4df6-818a-abd1b52000fb
# ╠═5b219ef8-0cbe-40a6-beed-335735a07d1b
# ╟─63f891da-12ed-4813-be9e-b08cd6c274cb
# ╠═89f2fecf-d07f-44e0-b39d-fbb6c509fbf6
# ╟─295942c8-0ef5-4c1a-bf99-561e4b4597c1
# ╠═d2011efa-d26c-4ed1-9688-bca52334c3bf
# ╠═423af97b-79ff-4fc8-8b67-26901b33ac50
# ╟─a7c2573a-19a2-4d4b-91a1-116e6d1e7c6a
# ╠═99173920-cb95-4a02-9d95-41cf04459183
# ╟─b94d5517-2b3a-4014-87be-8c258c104e02
# ╠═c0f4fe11-c149-43e8-afa8-962a035b8b15
# ╟─d6d76209-a5ba-4f0b-ab30-7632ed091ce5
# ╟─5eccbecc-519f-4a3c-80d6-cd77e118d024
# ╟─5691d49f-b2ec-4c24-ba5e-a1eee8f59476
# ╠═a2e5fa68-1657-448f-95f1-60d4118f4def
# ╟─84159453-6daa-4304-a9f2-b370391059c7
# ╠═6ef4ce4f-8ca5-4c76-9c13-a56fc6099773
# ╠═c9c4f960-2d1b-4131-a3c3-bc722c365aa9
# ╟─43a0026a-77d7-437d-b154-9baa1bcf0f80
# ╠═83745652-35b3-42c8-bc77-61485fb2250f
# ╟─0d335e76-d3c6-4591-ac43-dae5188e0400
# ╟─2952fb74-a351-493f-841f-466392382d89
# ╠═dccc8320-f60f-4849-8ef1-ec3a1990e629
# ╟─366f8aa3-696c-4068-804f-28ebbd8b1bbb
# ╠═3b072230-fdc4-4153-bf55-30e6dc6e63f5
# ╟─8d789020-4c1e-4802-a220-2551af74b373
# ╠═62d2e70b-c99b-43be-9f39-f67377d1855c
# ╟─1ae1300f-e194-4f13-9dc5-6c94d15fe9e0
# ╠═7c079152-43a9-4290-81b4-04788a4eb3bc
# ╠═2628ace2-b983-4c97-8598-e118675df304
# ╠═62dc62c3-ee4b-4148-aa1a-6d048fbec40d
# ╟─2bd03e6a-1038-43bd-9065-17938b853374
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
