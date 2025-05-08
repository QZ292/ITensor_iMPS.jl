# iMPS initialization
function initialize_x(χ::Int, D::Int, name::String = "A")
    iX = [
        Index(χ, "$name 1"),
        Index(D, "$name 2"), 
        Index(χ, "$name 3"),
    ]
    tensorX = randn(χ,D,χ).+3
    tensorX = ITensor(tensorX,iX)
    factor = sqrt( (tensorX * tensorX)[1] )
    tensorX /= factor
    return tensorX
end
# canonicalize iMPS
function canon_iteration(A::ITensor,R::ITensor)
    M = R * A * delta(inds(A)[1], inds(R)[2])
    M = replaceinds(M, (inds(M)[1]=>inds(A)[1], inds(M)[2]=>inds(A)[2], inds(M)[3]=>inds(A)[3]))
    AL, Rn = qr(M, (inds(M)[1],inds(M)[2]), positive=true)
    AL = replaceinds(AL, (inds(AL)[3]=>inds(A)[3],))
    Rn = replaceinds(Rn, (inds(Rn)[1]=>inds(R)[1], inds(Rn)[2]=>inds(R)[2]))
    factor = (Rn * delta(inds(Rn)[1], inds(Rn)[2]))[1] #/dim(inds(Rn)[1])
    Rn /= factor
    AL = factor * AL
    frobenius = sqrt(sum((Rn - R) .^ 2))/sqrt(sum((R) .^ 2))
    return AL, Rn, frobenius
end
# Generate canonical form random iMPS
function canonical_initialize(χ::Int, D::Int, name::String = "a"; tol::Float64 = 1e-6, itr::Int = 30, printlog::Bool = true, indent::Int64 = 0)
    # initialization
    MPS_AL = nothing
    MPS_AR_rev = nothing
    MPS_error = 1.0
    MPS_iteration = 0
    # random initial R
    R_initial = rand(χ,χ)
    _, R_initial = qr(R_initial)
    # random initial MPS
    MPS_A = initialize_x(χ, D, name)
    # left canonicalization
    indR = [Index(χ, "left"), Index(χ, "r2")]
    MPS_R = ITensor(R_initial, indR)
    while MPS_error > tol && MPS_iteration < itr
        MPS_iteration += 1
        MPS_AL, MPS_R, MPS_error = canon_iteration(MPS_A, MPS_R)
    end
    # right canonicalization
    MPS_iteration = 0
    MPS_error = 1.0
    indL = [
        Index(χ, "right"),
        Index(χ, "l1"),
    ]
    MPS_L_rev = ITensor(R_initial, indL)
    MPS_A_rev = permute(MPS_A, [inds(MPS_A)[3], inds(MPS_A)[2], inds(MPS_A)[1]])
    while MPS_error > tol && MPS_iteration < itr
        MPS_iteration += 1
        MPS_AR_rev, MPS_L_rev, MPS_error = canon_iteration(MPS_A_rev, MPS_L_rev)
    end
    MPS_AR = permute(MPS_AR_rev, [inds(MPS_AR_rev)[3], inds(MPS_AR_rev)[2], inds(MPS_AR_rev)[1]])
    # center
    MPS_C = MPS_R * MPS_L_rev * delta(inds(MPS_L_rev)[2], inds(MPS_R)[2])
    # normailzation
    norm_c = norm(matrix(MPS_C))
    trace_al = sqrt(tensor_three(MPS_AL,1))
    trace_ar = sqrt(tensor_three(MPS_AR,2))
    MPS_C /= norm_c
    MPS_AL /= trace_al
    MPS_AR /= trace_ar
    # results
    mpsA = iMPS_canonical()
    mpsA.left = MPS_AL
    mpsA.center = MPS_C
    mpsA.right = MPS_AR
    #test normailzation
    println(repeat(" ", indent), "Generating random canonical iMPS \"",name,"\" w/ bond dim = ",χ, " and physical dim = ",D)
    errordetect = canon_test(mpsA, [2,1]; printlog = printlog, indent = indent+2, tol = tol)
    if printlog
        println(" ")
    end
    return mpsA, errordetect
end
# Generate canonical form given translational invariant iMPS
function canonicalize(MPS_A::ITensor; tol::Float64 = 1e-6, itr::Int = 30, printlog::Bool = true, indent::Int64 = 0)
    # size
    χ = size(MPS_A)[1]
    @assert size(MPS_A)[1] == size(MPS_A)[3]
    # initialization
    MPS_AL = nothing
    MPS_AR_rev = nothing
    MPS_error = 1.0
    MPS_iteration = 0
    # random initial R
    R_initial = rand(χ,χ)
    _, R_initial = qr(R_initial)
    # left canonicalization
    indR = [Index(χ, "left"), Index(χ, "r2")]
    MPS_R = ITensor(R_initial, indR)
    while MPS_error > tol && MPS_iteration < itr
        MPS_iteration += 1
        MPS_AL, MPS_R, MPS_error = canon_iteration(MPS_A, MPS_R)
    end
    # right canonicalization
    MPS_iteration = 0
    MPS_error = 1.0
    indL = [
        Index(χ, "right"),
        Index(χ, "l1"),
    ]
    MPS_L_rev = ITensor(R_initial, indL)
    MPS_A_rev = permute(MPS_A, [inds(MPS_A)[3], inds(MPS_A)[2], inds(MPS_A)[1]])
    while MPS_error > tol && MPS_iteration < itr
        MPS_iteration += 1
        MPS_AR_rev, MPS_L_rev, MPS_error = canon_iteration(MPS_A_rev, MPS_L_rev)
    end
    MPS_AR = permute(MPS_AR_rev, [inds(MPS_AR_rev)[3], inds(MPS_AR_rev)[2], inds(MPS_AR_rev)[1]])
    # center
    MPS_C = MPS_R * MPS_L_rev * delta(inds(MPS_L_rev)[2], inds(MPS_R)[2])
    # normailzation
    norm_c = norm(matrix(MPS_C))
    trace_al = sqrt(tensor_three(MPS_AL,1))
    trace_ar = sqrt(tensor_three(MPS_AR,2))
    MPS_C /= norm_c
    MPS_AL /= trace_al
    MPS_AR /= trace_ar
    # results
    mpsA = iMPS_canonical()
    mpsA.left = MPS_AL
    mpsA.center = MPS_C
    mpsA.right = MPS_AR
    #test normailzation
    println(repeat(" ", indent), "Canonicalizing given iMPS w/ bond dim = ",χ, " and physical dim = ",size(MPS_A)[2])
    errordetect = canon_test(mpsA, [2,1]; printlog = printlog, indent = indent+2, tol = tol)
    if printlog
        println(" ")
    end
    return mpsA, errordetect
end