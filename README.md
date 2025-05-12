# README

## Introduction

A light weighted Julia package based on [ITensors.jl](https://github.com/ITensor/ITensors.jl), for working with infinite matrix product states (iMPS).

<center> Still woking in progress!</center>

## Installation

The package is **not** registered. In case you want to use it, please add by Pkg command.

```julia
julia> using Pkg; Pkg.add(url="https://github.com/QZ292/ITensor_simple_iMPS.git")
julia> using ITensor_iMPS
```

## Function, Usage and Example
[Mutable Struct: iMPS_canonical](https://github.com/QZ292/ITensor_simple_iMPS#mutable-struct-imps_canonical)

[Function: size](https://github.com/QZ292/ITensor_simple_iMPS#function-size-redefined)

[Function: deepcopy](https://github.com/QZ292/ITensor_simple_iMPS#function-deepcopy-redefined)

[Function: canonical_initialize](https://github.com/QZ292/ITensor_simple_iMPS#function-canonical_initialize)

[Function: canonicalize](https://github.com/QZ292/ITensor_simple_iMPS?tab=readme-ov-file#function-canonicalize)

[Function: canonical_gauge](https://github.com/QZ292/ITensor_simple_iMPS?tab=readme-ov-file#function-canonical_gauge)

[Function: biorthonormal_initialize](https://github.com/QZ292/ITensor_simple_iMPS?tab=readme-ov-file#function-biorthonormal_initialize)

[Function: biorthonormalize](https://github.com/QZ292/ITensor_simple_iMPS?tab=readme-ov-file#function-biorthonormalize)

[Function: checkimps](https://github.com/QZ292/ITensor_simple_iMPS?tab=readme-ov-file#function-checkimps)

### Mutable Struct: iMPS_canonical

Defined to store translational-invariant iMPS in canonical form with 3 ITensor object

```julia
a = iMPS_canonical()    # initialize
typeof(a.left),typeof(a.center),typeof(a.right)
# output:
(ITensor, ITensor, ITensor)
```

### Function: size (redefined)

For a $(\chi,D,\chi)\-$dimensional MPS, the function returns a Tuple{Int64, Int64}$(\chi,D)$​ while asserting correct shape for all three part of a *iMPS_canonical* variable.

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

Generate random translational-invariant iMPS in canonical form. Currently the initial MPS are drew from a uniform distribution on the $\mathbb{R}^{\chi\times D\times\chi}$​ sphere with a positive-bias for better convergence.

##### Usage:

```julia
canonical_initialize(χ::Int, D::Int, name::String = "a"; mode::String = "default", tol::Float64 = 1e-6, itr::Tuple{Int64, Float64} = (30,1e-8))
```

|  Input   | Parameter |         Type          |                                                              |
| :------: | :-------: | :-------------------: | ------------------------------------------------------------ |
| Argument |     χ     |         Int64         | virtual bond dimension                                       |
| Argument |     D     |         Int64         | physical bond dimension                                      |
| Optional |   name    |        String         | name of output iMPS, naming the ITensor indices              |
| Keyword  |   mode    |        String         | accept values: <br />"default", "debug" with normal or detailed logging and output verification<br />"min" skips output verification to reduce computation time |
| Keyword  |    tol    |        Float64        | max error tolerance for output                               |
| Keyword  |    itr    | Tuple{Int64, Float64} | max iteration steps and max error tolerance for canonicalizing the input |

| Output type    | Tuple{iMPS_canonical, Bool}                                  |
| -------------- | ------------------------------------------------------------ |
| iMPS_canonical | random canonicalized iMPS                                    |
| Bool           | **true** if error of the output is larger than **tol**; <br />if **mode="min"** then output error will not be calculated and this value is not returned, i.e. return **iMPS_canonical** instead of **Tuple{iMPS_canonical, Bool}** |

##### Examples:

```julia
a = canonical_initialize(4,4;mode="min")
# --console--
Random canonical iMPS "a" w/ bond dim = 1 and physical dim = 2 generated.
#--output--
iMPS_canonical
 left  -> ITensor ord=3 (dim=4|id=704|"a1") (dim=4|id=782|"a2") (dim=4|id=174|"a3")
          NDTensors.Dense{Float64, Vector{Float64}}
 center-> ITensor ord=2 (dim=4|id=157|"left") (dim=4|id=504|"right")
          NDTensors.Dense{Float64, Vector{Float64}}
 right -> ITensor ord=3 (dim=4|id=704|"a1") (dim=4|id=782|"a2") (dim=4|id=174|"a3")
          NDTensors.Dense{Float64, Vector{Float64}}
```

```julia
a,_ = canonical_initialize(1,2;mode="default.4")
# --console--
a,_ = canonical_initialize(1,2;mode="default.4")
#--output--
(iMPS_canonical{...}, false)
```

```julia
a,_ = canonical_initialize(1,1;mode="debug.0")
# --console--
Generating random canonical iMPS "a" w/ bond dim = 1 and physical dim = 1
  Check canonicality of given iMPS.
  • norm center=1.0
  • orthonormal left=0.0, orthonormal right=0.0
  • center left, right diff=0.0
#--output--
(iMPS_canonical{...}, false)
```

### Function: canonicalize

Return the canonical form of the given translational-invariant iMPS $\dots-A-A-A-\dots$

##### Usage:

```julia
canonicalize(MPS_A::ITensor; mode::String = "default", tol::Float64 = 1e-6, itr::Tuple{Int64, Float64} = (30,1e-8))
```

|  Input   | Parameter |         Type          |                                                              |
| :------: | :-------: | :-------------------: | ------------------------------------------------------------ |
| Argument |   MPS_A   |        ITensor        | input ITensor, which is the side tensor of the iMPS          |
| Keyword  |   mode    |        String         | accept values: <br />"default", "debug" with normal or detailed logging and output verification<br />"min" skips output verification to reduce computation time |
| Keyword  |    tol    |        Float64        | max error tolerance for output                               |
| Keyword  |    itr    | Tuple{Int64, Float64} | max iteration steps and max error tolerance for canonicalizing the input |

| Output type    | Tuple{iMPS_canonical, Bool}                                  |
| -------------- | ------------------------------------------------------------ |
| iMPS_canonical | canonical form of the given iMPS                             |
| Bool           | **true** if error of the output is larger than **tol**; <br />if **mode="min"** then output error will not be calculated and this value is not returned, i.e. return **iMPS_canonical** instead of **Tuple{iMPS_canonical, Bool}** |

##### Examples:

```julia
a0 = ITensor_iMPS.initialize_x(2,3)
a,_ = canonicalize(a0;mode="debug.0",tol=1e-3,itr=(10,1e-6))
#--console--
Canonicalizing given iMPS w/ bond dim = 2 and physical dim = 3
  Check canonicality of given iMPS.
  • norm center=1.0
  • orthonormal left=4.965068306494546e-16, orthonormal right=4.577566798522237e-16
  • center left, right diff=1.466766694599318e-7
#--output--
(iMPS_canonical{...}, false)
```

### Function: canonical_gauge

Remove the gauge redundancy by inserting unitary matrices $U,V$ that diagonalize the center $C$, on virtual bonds.

##### Usage:

```julia
canonical_gauge(a::iMPS_canonical; mode::String = "default", tol::Float64 = 1e-6)
```

|  Input   | Parameter |      Type      |                                                              |
| :------: | :-------: | :------------: | ------------------------------------------------------------ |
| Argument |     a     | iMPS_canonical | input canonical form infinite MPS                            |
| Keyword  |   mode    |     String     | accept values: <br />"default", "debug" with normal or detailed logging and output verification<br />"min" skips output verification to reduce computation time |
| Keyword  |    tol    |    Float64     | max error tolerance for output                               |

| Output type    | Tuple{iMPS_canonical, Bool}                                  |
| -------------- | ------------------------------------------------------------ |
| iMPS_canonical | output iMPS                                                  |
| Bool           | **true** if error of the output is larger than **tol**; <br />if **mode="min"** then output error will not be calculated and this value is not returned, i.e. return **iMPS_canonical** instead of **Tuple{iMPS_canonical, Bool}** |

##### Examples:

```julia
a1 = canonical_gauge(a;mode="debug")
#--console--
Removing gauge ambiguity of given iMPS.
  Check canonicality of given iMPS.
  • norm center=1.0000000000000002
  • diag center=1.0e50
  • orthonormal left=2.5133742693021536e-16, orthonormal right=2.830524433501838e-16
  • center left, right diff=1.466766694731434e-7
#--output--
(iMPS_canonical{...}, false)
```

### Function: biorthonormal_initialize

Return a pair of random biorthonormal infinite MPS in canonical form.

##### Usage:

```julia
biorthonormal_initialize(χ::Int64, D::Int64; mode::String = "default", tol::Float64 = 1e-6, itr::Int = 1, citr::Tuple{Int64, Float64} = (30,1e-8))
```

|  Input   | Parameter |         Type          |                                                              |
| :------: | :-------: | :-------------------: | ------------------------------------------------------------ |
| Argument |     χ     |         Int64         | virtual bond dimension                                       |
| Argument |     D     |         Int64         | physical bond dimension                                      |
| Keyword  |   mode    |        String         | accept values: <br />"default", "debug" with normal or detailed logging and output verification<br />"min" skips output verification to reduce computation time |
| Keyword  |    tol    |        Float64        | max error tolerance for output                               |
| Keyword  |    itr    |         Int64         | max iteration steps during biorthonormalizing generated iMPS pair, 1 is good enough in most cases |
| Keyword  |   citr    | Tuple{Int64, Float64} | max iteration steps and max error tolerance for canonicalizing the iMPS pair during random iMPS generation |

| Output type    | Tuple{iMPS_canonical, iMPS_canonical, Bool}                  |
| -------------- | ------------------------------------------------------------ |
| iMPS_canonical | random biorthonormal infinite MPS "a" in canonical form      |
| iMPS_canonical | random biorthonormal infinite MPS "b" in canonical form      |
| Bool           | **true** if error of the output is larger than **tol**; <br />if **mode="min"** then output error will not be calculated and this value is not returned, i.e. return **Tuple{iMPS_canonical, iMPS_canonical}** instead of **Tuple{iMPS_canonical, iMPS_canonical, Bool}** |

##### Examples:

```julia
biorthonormal_initialize(2,2)
#--console--
Random canonical iMPS "a" w/ bond dim = 2 and physical dim = 2 generated.
 
Random canonical iMPS "b" w/ bond dim = 2 and physical dim = 2 generated.
 
Random biorthonormal iMPS pair w/ bond dim = 2 and physical dim = 2 generated in 1 iteration.
  Check biorthonormality of given iMPS pair.
  • left biorthonormal = 6.131407058289218e-16
  • right biorthonormal = 1.066334393690002e-15
  • center left, right diff = [2.506090114876823e-9, 1.650317478538157e-9]
  • diag left, right = [5.577606463939291e11, 5.577606463939291e11]
  • RDM left, right diff = 9.385611898852907e-17
#--output--
(iMPS_canonical{...}, false)
```

### Function: biorthonormalize

Biorthonormalize the given iMPS pair.

##### Usage:

```julia
biorthonormalize(a::iMPS_canonical, b::iMPS_canonical; mode::String = "default", tol::Float64 = 1e-6, itr::Int = 1)
```

|  Input   | Parameter |      Type      |                                                              |
| :------: | :-------: | :------------: | ------------------------------------------------------------ |
| Argument |     a     | iMPS_canonical | input canonical form infinite MPS "a"                        |
| Argument |     b     | iMPS_canonical | input canonical form infinite MPS "b"                        |
| Keyword  |   mode    |     String     | accept values: <br />"default", "debug" with normal or detailed logging and output verification<br />"min" skips output verification to reduce computation time |
| Keyword  |    tol    |    Float64     | max error tolerance for output                               |
| Keyword  |    itr    |     Int64      | max iteration steps during biorthonormalizing generated iMPS pair, 1 is good enough in most cases |

| Output type    | Tuple{iMPS_canonical, iMPS_canonical, Bool}                  |
| -------------- | ------------------------------------------------------------ |
| iMPS_canonical | biorthonormalized infinite MPS "a" in canonical form         |
| iMPS_canonical | biorthonormalized infinite MPS "b" in canonical form         |
| Bool           | **true** if error of the output is larger than **tol**; <br />if **mode="min"** then output error will not be calculated and this value is not returned, i.e. return **Tuple{iMPS_canonical, iMPS_canonical}** instead of **Tuple{iMPS_canonical, iMPS_canonical, Bool}** |

##### Examples:

```julia
a = canonical_initialize(8,8;mode="min.2")
b = canonical_initialize(8,8;mode="min.2")
biorthonormalize(a, b; mode = "debug.2", itr = 2)
#--console--
  Random canonical iMPS "a" w/ bond dim = 8 and physical dim = 8 generated.
 
  Random canonical iMPS "a" w/ bond dim = 8 and physical dim = 8 generated.
 
  Biorthonormalizing given iMPS w/ bond dim = 8 and physical dim = 8
 
    Biortho itr = 0
      Check biorthonormality of given iMPS pair.
      • left biorthonormal = 3.1643902826377253
      • right biorthonormal = 3.161461300481157
      • center left, right diff = [5.876180273663705e-10, 8.62994931238692e-10]
      • diag left, right = [0.0001501676952522035, 0.0001501676952522035]
      • RDM left, right diff = 0.05955716138673125
 
    Biortho itr = 1
    > after biortho
      Check biorthonormality of given iMPS pair.
      • left biorthonormal = 3.508564750841953e-15
      • right biorthonormal = 2.6320208773439526e-15
      • center left, right diff = [2.2168596955197628e-10, 3.241000611588641e-10]
      • diag left, right = [1.6030959480656694e-5, 1.6030959480656694e-5]
      • RDM left, right diff = 0.3475979759910587
    > after gauge
      • unitarity of GL, GR = 5.811665519990758e-15, 4.421619646969187e-15
      Check biorthonormality of given iMPS pair.
      • left biorthonormal = 4.182813166126303e-14
      • right biorthonormal = 1.7470494651503157e-14
      • center left, right diff = [7.353428898599791e-10, 1.1309286696382305e-9]
      • diag left, right = [4.897009216352439e8, 4.897009216352439e8]
      • RDM left, right diff = 2.857315103756729e-16
 
    Biortho itr = 2
    > after biortho
      Check biorthonormality of given iMPS pair.
      • left biorthonormal = 2.0176143453228117e-14
      • right biorthonormal = 1.7271125508266392e-14
      • center left, right diff = [2.463266027289115e-10, 3.788407078497747e-10]
      • diag left, right = [0.0012833790245745142, 0.0012833790245745142]
      • RDM left, right diff = 1.4332117117319505
    > after gauge
      • unitarity of GL, GR = 5.12369307251331e-16, 5.741725067399372e-16
      Check biorthonormality of given iMPS pair.
      • left biorthonormal = 2.386163283742619e-14
      • right biorthonormal = 1.5617057558918153e-14
      • center left, right diff = [2.4632659927419515e-10, 3.7884070718100597e-10]
      • diag left, right = [1.860885821094701e8, 1.860885821094701e8]
      • RDM left, right diff = 5.338614275121067e-16
 
  Done
 #--output--
(iMPS_canonical{...}, iMPS_canonical{...}, false)
```

### Function: checkimps

Check properties of input iMPS

##### Usage:

```julia
checkimps(a::iMPS_canonical, b::Union{iMPS_canonical, Nothing}=nothing; list::Vector{Int64} = [1,2], tol::Float64 = 1e-3, printlog::Bool = true, indent::Int64 = 0)
```

|  Input   | Parameter |                 Type                 |                                                              |
| :------: | :-------: | :----------------------------------: | ------------------------------------------------------------ |
| Argument |     a     |            iMPS_canonical            | input canonical form infinite MPS "a"                        |
| Argument |     b     | Union{<br />iMPS_canonical, Nothing} | input canonical form infinite MPS "b", <br />keep empty(nothing) when checking the properties of one single iMPS |
| Keyword  |   list    |            Vector{Int64}             | list of properties to be examined, both 1 and 2 iMPS cases accepts values [1,2,3] |
| Keyword  |    tol    |               Float64                | max error allowed for the input iMPS                         |
| Keyword  | printlog  |                 Bool                 | whether print the detailed results, or simply return the final true/false result. |
| Keyword  |  indent   |                Int64                 | the number of indent spaces, for the print output            |

| Output type | Bool                                                   |
| ----------- | ------------------------------------------------------ |
| Bool        | **true** if error of the output is larger than **tol** |

**Case 1 of 1 iMPS:** 
		check the left/right canonicalization of  the left/right local tensor $A_L/A_R$, return values are the Frobenius norm of $[A_L^\dagger A_L -I,A_R^\dagger A_R -I]$; 
		check the difference between $A_L C$ and $C A_R$, which is expected to be the same $A_C$, return value is the Frobenius norm of $A_L C - C A_R$.

**Case 2 of 1 iMPS:** check normalization of the center $C$, return value is $\text{norm}(C)$.

**Case 3 of 1 iMPS:** check diagonalization of the center $C$, return values is $\frac{\text{min}({\text{diagonal emement}})}{\text{max}(\text{off diag element})+10^{-50}}$.

**Case 1 of 2 iMPS:** check left/right biorthonormalization of $A_L,B_L$ and $A_R,B_R$,  return values are the Frobenius norm of $[A_L B_L^\dagger -I,A_R B_R^\dagger -I]$;

**Case 2 of 2 iMPS:** check the difference between $A_L C, C A_R$  and $B_L D, D B_R$, return values are the Frobenius norm of $[A_L C - C A_R,B_L D - D B_R]$.

**Case 3 of 2 iMPS:** check the diagonalization and difference of $\rho_L = CD^\dagger$ and $\rho_R = D^\dagger C$, when the gauge ambiguity of the biorthonormal pair is removed, it is expected that $\rho_L=\rho_R$  are the same diagonal matrix.

##### Examples:

```julia
checkimps(a;list=[3,2,1],tol=1e-3,printlog=true)
#--console--
Check canonicality of given iMPS.
• diag center=2.116927050643497e10
• norm center=0.38146487350496094
• orthonormal left=10.42079212727393, orthonormal right=4.2455449590347385
• center left, right diff=2.3877612699297965e-9
#--output--
true
```

```julia
checkimps(a,b;list=[1,2,3],tol=1e-3,printlog=true)
#--console--
Check biorthonormality of given iMPS pair.
• left biorthonormal = 1.0456705022016911e-14
• right biorthonormal = 6.083223048914707e-15
• center left, right diff = [2.387761269929797e-9, 3.9275077336323146e-10]
• diag left, right = [6.845856112346311e8, 6.845856112346311e8]
• RDM left, right diff = 2.4944522588680643e-16
#--output--
false
```
