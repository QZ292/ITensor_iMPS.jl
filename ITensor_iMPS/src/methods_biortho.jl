# biorthonormalize iMPS
function biortho_normal(A::ITensor, B::ITensor, flag::Int64)
    if size(A)==size(B) && size(A)[1] == size(A)[3]
        χ = size(A)[1]    # MPS bond dimension
        D = size(A)[2]    # physical dimension
        iM = [
            Index(χ, "left")
            Index(χ, "right")
        ]
    else
        println("ERROR@ tensor size mismatch: ",size(A),"≠",size(B))
        return nothing
    end
    # rehsape
    c1 = combiner(inds(B)[1], inds(A)[1], tags="left")
    c2 = combiner(inds(B)[3], inds(A)[3], tags="right")
    transfer = matrix( dag(B) * A * delta(inds(A)[2], inds(B)[2]) * c1 * c2 )
    # left eigen 
    if flag == 1
        egval, egvec = Arpack.eigs(transfer', nev=1, which=:LR)
        egval = real(egval)[1]
        egvec = real(egvec')
        egvec = ITensor(reshape(egvec, χ, χ), (inds(B)[1], inds(A)[1]))
        # factorize X into new AL, BL
        ux,sx,vx = svd(egvec, inds(egvec)[1])
        sx = sqrt.(sx)
        px = sx * vx
        qx = dag(ux * sx)
        px = replaceinds( px, (inds(px)[1]=>iM[1],inds(px)[2]=>inds(A)[1]) )
        qx = replaceinds( qx, (inds(qx)[2]=>iM[1],inds(qx)[1]=>inds(B)[1]) )
        # compute new A, B
        px_inv = ITensor( inv(matrix(px)), (inds(A)[3], iM[2]) )
        qx_inv = ITensor( inv(matrix(qx)), (iM[2], inds(B)[3]) )
        An = (egval^(-1)) * px * A * px_inv
        Bn = qx * B * qx_inv
        An = replaceinds( An, (inds(An)[1]=>inds(A)[1], inds(An)[2]=>inds(A)[2], inds(An)[3]=>inds(A)[3]))
        Bn = replaceinds( Bn, (inds(Bn)[1]=>inds(B)[1], inds(Bn)[2]=>inds(B)[2], inds(Bn)[3]=>inds(B)[3]))
    # right eigen
    elseif flag == 2
        egval, egvec = Arpack.eigs(transfer, nev=1, which=:LR)
        egval = real(egval)[1]
        egvec = real(egvec)
        egvec = ITensor(reshape(egvec, χ, χ), (inds(B)[3], inds(A)[3]))
        # factorize X into new AR, BR
        ux,sx,vx = svd(egvec, inds(egvec)[1])
        sx = sqrt.(sx)
        px = sx * vx
        qx = dag(ux * sx)
        px = replaceinds( px, (inds(px)[1]=>iM[2], inds(px)[2]=>inds(A)[3]) )
        qx = replaceinds( qx, (inds(qx)[2]=>iM[2], inds(qx)[1]=>inds(B)[3]) )
        # compute new A, B
        px_inv = ITensor( inv(matrix(px)), (inds(A)[1], iM[1]) )
        qx_inv = ITensor( inv(matrix(qx)), (iM[1], inds(B)[1]) )
        An = (egval^(-1)) * px_inv * A * px
        Bn = qx_inv * B * qx
        An = replaceinds( An, (inds(An)[1]=>inds(A)[1], inds(An)[2]=>inds(A)[2], inds(An)[3]=>inds(A)[3]))
        Bn = replaceinds( Bn, (inds(Bn)[1]=>inds(B)[1], inds(Bn)[2]=>inds(B)[2], inds(Bn)[3]=>inds(B)[3]))
    end
    return An, Bn, px, qx
end
function biortho_iteration(a::iMPS_canonical, b::iMPS_canonical)
    χ = size(a.left)[1]
    a.left, b.left, pl, ql = biortho_normal(a.left, b.left,1)
    a.right, b.right, pr, qr = biortho_normal(a.right, b.right,2)
    # rename indices
    iC = [
        Index(χ, "cl"),  #1
        Index(χ, "cr"),  #2
        Index(χ, "mid left"),  #3
        Index(χ, "mid right"), #4
        Index(χ, "dl"),  #5
        Index(χ, "dr")   #6  
    ]
    a.center = replaceinds( a.center, (inds(a.center)[1]=>iC[3], inds(a.center)[2]=>iC[4]) )
    b.center = replaceinds( b.center, (inds(b.center)[1]=>iC[3], inds(b.center)[2]=>iC[4]) )
    pl = replaceinds( pl, (inds(pl)[1]=>iC[1], inds(pl)[2]=>iC[3]) )
    pr = replaceinds( pr, (inds(pr)[1]=>iC[2], inds(pr)[2]=>iC[4]) )
    ql = replaceinds( ql, (inds(ql)[1]=>iC[3], inds(ql)[2]=>iC[5]) )
    qr = replaceinds( qr, (inds(qr)[1]=>iC[4], inds(qr)[2]=>iC[6]) )
    # contract new Pl Pr Ql Qr
    a.center = pl*a.center*pr
    b.center = ql*b.center*qr
    return a, b
end
# remove gauge ambiguity
function biortho_gauge(a::iMPS_canonical, b::iMPS_canonical; printlog::Bool = true, indent::Int64 = 2)
    χ = size(a.left)[1]
    iG = [
        Index(χ, "gl-l"),    #1
        Index(χ, "gl-r"),    #2
        Index(χ, "gli-l"),   #3
        Index(χ, "gli-r"),   #4
        Index(χ, "gr-l"),    #5
        Index(χ, "gr-r"),    #6
        Index(χ, "gri-l"),   #7
        Index(χ, "gri-r")    #8
    ]
    # diagonalize RDM: L=CD†, R=D†C
    rhoL = matrix(a.center * dag(b.center) * delta(inds(a.center)[2],inds(b.center)[2]))
    rhoR = matrix(dag(b.center) * a.center * delta(inds(a.center)[1],inds(b.center)[1]))
    # diagonalize to find G
    FG = eigen(matrix(rhoL))
    GL = ITensor(FG.vectors,(iG[1],iG[2]))
    GLinv = ITensor(inv(FG.vectors),(iG[3],iG[4]))
    FR = eigen(matrix(rhoR))
    GR = ITensor(FR.vectors,(iG[5],iG[6]))
    GRinv = ITensor(inv(FR.vectors),(iG[7],iG[8]))
    # check unitary
    if printlog
        uni_l = matrix(GL * GLinv * delta(iG[2],iG[3]))-I(χ)
        uni_r = matrix(GR * GRinv * delta(iG[6],iG[7]))-I(χ)
        println(repeat(" ", indent), "• unitarity of GL, GR = ", norm(uni_l),", ", norm(uni_r))
    end
    # new tensors
    an = iMPS_canonical()
    bn = iMPS_canonical()
    # update A_L
    an.left = GLinv * a.left * delta(iG[4],inds(a.left)[1])
    an.left = an.left * GL *delta(inds(an.left)[3],iG[1])
    # update C
    an.center = GLinv * a.center * delta(iG[4],inds(a.center)[1])
    an.center = an.center * GR * delta(iG[5],inds(an.center)[2])
    # update A_R
    an.right = GRinv * a.right * delta(iG[8],inds(a.right)[1])
    an.right = an.right * GR * delta(inds(an.right)[3],iG[5])
    # update B_L
    bn.left = dag(GL) * b.left * delta(iG[1],inds(b.left)[1])
    bn.left = bn.left * dag(GLinv) *delta(inds(bn.left)[3],iG[4])
    # update D
    bn.center = dag(GL) * b.center * delta(iG[1],inds(b.center)[1])
    bn.center = bn.center * dag(GRinv) * delta(iG[8],inds(bn.center)[2])
    # update B_R
    bn.right = dag(GR) * b.right * delta(iG[5],inds(b.right)[1])
    bn.right = bn.right * dag(GRinv) * delta(inds(bn.right)[3],iG[8])
    # rename indices
    an.left = replaceinds(an.left, (inds(an.left)=>inds(a.left)))
    an.center = replaceinds(an.center, (inds(an.center)=>inds(a.center)))
    an.right = replaceinds(an.right, (inds(an.right)=>inds(a.right)))
    bn.left = replaceinds(bn.left, (inds(bn.left)=>inds(b.left)))
    bn.right = replaceinds(bn.right, (inds(bn.right)=>inds(b.right)))
    bn.center = replaceinds(bn.center, (inds(bn.center)=>inds(b.center)))
    return an, bn
end
# Initialize two biorthonormal iMPS
function biorthonormal_initialize(χ::Int64, D::Int64; mode::String = "default", tol::Float64 = 1e-6, itr::Int = 1, citr::Tuple{Int64, Float64} = (30,1e-8))
    parts = split(mode, '.'; limit=2)
    mode = parts[1]
    if length(parts) == 1
        indent = 0
    else
        indent = tryparse(Int, parts[2])
        if indent === nothing
            indent = 0
        else
            indent = max(indent,0)
        end
    end
    # output
    if mode=="min" || mode=="default"
        a,erra = canonical_initialize(χ,D,"a";tol=tol,itr=citr,mode=string("default.",string(indent)))
        b,errb = canonical_initialize(χ,D,"b";tol=tol,itr=citr,mode=string("default.",string(indent)))
        if erra
            println("ERROR@ generate iMPS \"a\"")
            return iMPS_canonical(),iMPS_canonical(),true
        elseif errb
            println("ERROR@ generate iMPS \"b\"")
            return iMPS_canonical(),iMPS_canonical(),true
        else
            nothing
        end
    elseif mode == "debug"
        println(repeat(" ", indent), "Generating random biorthonormal iMPS pair w/ bond dim = $χ and physical dim = $D.")
        println(" ")
        a,erra = canonical_initialize(χ,D,"a";tol=tol,itr=citr,mode=string("debug.",string(indent+2)))
        b,errb = canonical_initialize(χ,D,"b";tol=tol,itr=citr,mode=string("debug.",string(indent+2)))
        if erra
            println("ERROR@ generate iMPS \"a\"")
            return iMPS_canonical(),iMPS_canonical(),true
        elseif errb
            println("ERROR@ generate iMPS \"b\"")
            return iMPS_canonical(),iMPS_canonical(),true
        else
            nothing
        end
    else
        println("ERROR@ mode: \"$mode\"")
        println(" ")
        return iMPS_canonical(),iMPS_canonical(),true
    end
    # biorthonormailization
    bitr = 0
    if mode == "debug"
        println(repeat(" ", indent+2), "Biortho itr = ", bitr)
        biortho_test(a, b, [1,2,3]; tol=tol, printlog = true, indent=indent+4)
        println(" ")
    end
    while bitr < itr
        bitr += 1
        a, b = biortho_iteration(a, b)           
        if mode == "debug"
            println(repeat(" ", indent+2), "Biortho itr = ", bitr)
            println(repeat(" ", indent+2), "> after biortho")
            biortho_test(a, b, [1,2,3];tol=tol, printlog = true, indent=indent+4)
            println(repeat(" ", indent+2), "> after gauge")
            a, b = biortho_gauge(a, b; printlog = true, indent=indent+4)
            errorflag = biortho_test(a, b, [1,2,3];tol=tol, printlog = true, indent=indent+4)
            println(" ")
        else
            a, b = biortho_gauge(a, b; printlog = false)
        end
    end
    # output
    if mode == "min"
        println(repeat(" ", indent), "Random biorthonormal iMPS pair w/ bond dim = $χ and physical dim = $D generated in $bitr iteration.")
        println(" ")
        return a,b
    elseif mode == "default"
        println(repeat(" ", indent), "Random biorthonormal iMPS pair w/ bond dim = $χ and physical dim = $D generated in $bitr iteration.")
        errorflag = biortho_test(a, b, [1,2,3]; tol=tol, printlog = true, indent=indent+2)
        println(" ")
    elseif mode == "debug"
        println(repeat(" ", indent), "Done")
        println(" ")
    end
    return a, b, errorflag
end
# Biorthonormalize two canonical iMPS
function biorthonormalize(a::iMPS_canonical, b::iMPS_canonical; mode::String = "default", tol::Float64 = 1e-6, itr::Int = 1)
    @assert size(a) == size(b)
    χ = size(a)[1]
    D = size(a)[2]
    errorflag = false
    parts = split(mode, '.'; limit=2)
    mode = parts[1]
    if length(parts) == 1
        indent = 0
    else
        indent = tryparse(Int, parts[2])
        if indent === nothing
            indent = 0
        else
            indent = max(indent,0)
        end
    end
    # output
    if mode=="min" || mode=="default"
        nothing
    elseif mode == "debug"
        println(repeat(" ", indent), "Biorthonormalizing given iMPS w/ bond dim = ",χ, " and physical dim = ",D)
        println(" ")
    else
        println("ERROR@ mode: \"$mode\"")
        println(" ")
        return a,b,true
    end
    # biorthonormailization
    bitr = 0
    if mode == "debug"
        println(repeat(" ", indent+2), "Biortho itr = ", bitr)
        biortho_test(a, b, [1,2,3]; printlog = true,tol=tol,indent=indent+4)
        println(" ")
    end
    while bitr < itr
        bitr += 1
        a, b = biortho_iteration(a, b)
        if mode == "debug"
            println(repeat(" ", indent+2), "Biortho itr = ", bitr)
            println(repeat(" ", indent+2), "> after biortho")
            biortho_test(a, b, [1,2,3]; printlog = true,tol=tol,indent=indent+4)
            println(repeat(" ", indent+2), "> after gauge")
            a, b = biortho_gauge(a, b; printlog = true ,indent=indent+4)
            errorflag = biortho_test(a, b, [1,2,3]; printlog = true,tol=tol,indent=indent+4)
            println(" ")
        else
            a, b = biortho_gauge(a, b; printlog = false)
        end
    end
    # output
    if mode == "min"
        println(repeat(" ", indent), "Given iMPS pair w/ bond dim = ",χ, " and physical dim = ",D," biorthonormalized in $bitr iterations.")
        println(" ")
        return a,b
    elseif mode == "default"
        println(repeat(" ", indent), "Given iMPS pair w/ bond dim = ",χ, " and physical dim = ",D," biorthonormalized in $bitr iterations.")
        errorflag = biortho_test(a, b, [1,2,3]; tol=tol, printlog = true, indent=indent+2)
        println(" ")
    elseif mode == "debug"
        println(repeat(" ", indent), "Done")
        println(" ")
    end
    return a, b, errorflag
end
