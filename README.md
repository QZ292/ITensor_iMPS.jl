# README

## Introduction

A light weighted Julia package based on [ITensors.jl](https://github.com/ITensor/ITensors.jl), for working with infinite matrix product states (iMPS).

<center> Still woking in progress!</center>

## Installation

The package is **not** registered. In case you want to use it, please add by Pkg command.

```julia
julia> using Pkg; Pkg.add(url="https://github.com/QZ292/ITensor_simple_iMPS.jl.git")
julia> using ITensor_iMPS
```

## Function, Usage and Example

### Mutable Struct: iMPS_canonical

Defined to store translational-invariant iMPS in canonical form with 3 ITensor object

```julia
a = iMPS_canonical()    # initialize
typeof(a.left),typeof(a.center),typeof(a.right)
# output:
(ITensor, ITensor, ITensor)
```

### Function: size (redefined)

For a $(\chi,D,\chi)-$dimensional MPS, the function returns a Tuple{Int64, Int64}$(\chi,D)$​ while asserting correct shape for all three part of a *iMPS_canonical* variable.

##### Usage:

```julia
size(a::iMPS_canonical)
```

##### Examples:

```julia
a,_ = canonical_initialize(3,4)
@show size(a), typeof(size(a))
# output:
(size(a), typeof(size(a))) = ((3, 4), Tuple{Int64, Int64})
```

### Function: deepcopy (redefined)

Just do the same work as the original one.

##### Usage:

```julia
deepcopy(a::iMPS_canonical)
```

### Function: canonical_initialize

Generate random translational-invariant iMPS in canonical form. Currently the initial MPS are drew from a uniform distribution on the $\mathbb{R}^{\chi \times D \times \chi}$​ sphere with a positive-bias for better convergence.

##### Usage:

```julia
canonical_initialize(χ::Int, D::Int, name::String = "a"; tol::Float64 = 1e-6, itr::Int = 30, mode::String = "default", indent::Int64 = 0)
```



##### Examples:

```julia

```



### Function: canonicalize



##### Usage:

```julia
canonicalize(MPS_A::ITensor; tol::Float64 = 1e-6, itr::Int = 30, mode::String = "default", indent::Int64 = 0)
```



##### Examples:

```julia

```



### Function: canonical_gauge



##### Usage:

```julia
canonical_gauge(a::iMPS_canonical; tol::Float64 = 1e-6, mode::String = "default", indent::Int64 = 0)
```



##### Examples:

```julia

```



### Function: biorthonormal_initialize



##### Usage:

```julia
biorthonormal_initialize(χ::Int64, D::Int64; tol::Float64 = 1e-6, itr::Int = 1, mode::String = "default", indent::Int64 = 0)
```



##### Examples:

```julia

```



### Function: biorthonormalize



##### Usage:

```julia
biorthonormalize(a::iMPS_canonical, b::iMPS_canonical; tol::Float64 = 1e-6, itr::Int = 1, mode::String = "default", indent::Int64 = 0)
```



##### Examples:

```julia

```



### Function: checkimps



##### Usage:

```julia
checkimps(a::iMPS_canonical, b::Union{iMPS_canonical, Nothing}=nothing; list::Vector{Int64} = [1,2], printlog::Bool = true, tol::Float64 = 1e-3, indent::Int64 = 0)
```



##### Examples:

```julia

```





