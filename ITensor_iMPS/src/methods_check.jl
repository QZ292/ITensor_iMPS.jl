# operations for three leg tensors
function tensor_three(A::ITensor, flag::Int)
    Ap = dag(deepcopy(A))
    indAp = [
        Index(dim(inds(A)[1]); tags="dual1"),
        Index(dim(inds(A)[2]); tags="dual2"),
        Index(dim(inds(A)[3]); tags="dual3"),
    ]
    Ap = replaceinds(Ap, [inds(A)[1]=>indAp[1], inds(A)[2]=>indAp[2], inds(A)[3]=>indAp[3]])
    if flag == 1    # left_normalization
        factor = Ap * A * delta(inds(A)[1], inds(Ap)[1]) * delta(inds(A)[2], inds(Ap)[2])
        factor = factor[1][1]
    elseif flag == 2    # right normalization
        factor = Ap * A * delta(inds(A)[3], inds(Ap)[3]) * delta(inds(A)[2], inds(Ap)[2])
        factor = factor[1][1]
    elseif flag == 0    # return full trace
        factor = Ap * A * delta(inds(A)[1], inds(Ap)[1]) * delta(inds(A)[2], inds(Ap)[2]) * delta(inds(A)[3], inds(Ap)[3])
        factor = factor[1]
    elseif flag == -1    # return matrix left
        factor = Ap * A * delta(inds(A)[1], inds(Ap)[1]) * delta(inds(A)[2], inds(Ap)[2])
    elseif flag == -10    # return left error
        factor = Ap * A * delta(inds(A)[1], inds(Ap)[1]) * delta(inds(A)[2], inds(Ap)[2])
        factor = matrix(factor) - I(dim(inds(A)[3]))
        err = norm(sqrt(sum(factor*factor')))
        factor = err
    elseif flag == -2    # return matrix right
        factor = Ap * A * delta(inds(A)[3], inds(Ap)[3]) * delta(inds(A)[2], inds(Ap)[2])
    elseif flag == -20    # return right error
        factor = Ap * A * delta(inds(A)[3], inds(Ap)[3]) * delta(inds(A)[2], inds(Ap)[2])
        factor = matrix(factor) - I(dim(inds(A)[1]))
        err = norm(sqrt(sum(factor*factor')))
        factor = err
    elseif flag == 3    # return transfer matrix
        c1 = combiner(inds(A)[1], inds(Ap)[1]; tags="aa1")
        c2 = combiner(inds(A)[3], inds(Ap)[3]; tags="aa3")
        factor = Ap * A * delta(inds(A)[2], inds(Ap)[2]) * c1 * c2
    else    # type not defined
        println("ERROR@ normalization type = ",flag)
        factor = 0
    end
    return factor
end
# trace of four leg tensors
function tensor_trace(A::ITensor)
    trace = 0
    delta_mat_1 = delta(inds(A)[1],inds(A)[3])
    delta_mat_2 = delta(inds(A)[2],inds(A)[4])
    trace = A*delta_mat_1*delta_mat_2
    return diag(trace)[1]
end
# Check biorthonormality of two iMPS
function max_offdiag_min_diag(A::AbstractMatrix)
    @assert size(A, 1) == size(A, 2)
    n = size(A, 1)
    diag_elements = diag(A)
    diag_norms = abs.(diag_elements)
    diag_min = diag_elements[argmin(diag_norms)]
    offdiag_mask = .!I(n)
    offdiag_elements = A[offdiag_mask]
    offdiag_norms = abs.(offdiag_elements)
    max_offdiag = offdiag_elements[argmax(offdiag_norms)]
    ratio = norm(max_offdiag) / norm(diag_min)
    return 1/(ratio+1e-50)
end
function biortho_test_flag(a::iMPS_canonical, b::iMPS_canonical, flag::Int; printlog::Bool = true, indent::Int64 = 2, tol::Float64 = 1e-3)
    χ = size(a.left)[1]
    _ = size(a.left)[2]
    errordetect = false
    spaces = repeat(" ", indent)
    if flag == 1    # check biorthonormality
        bioleft = matrix( a.left * dag(b.left) * delta(inds(a.left)[1],inds(b.left)[1]) * delta(inds(a.left)[2],inds(b.left)[2]) ) - I(χ)
        bioright = matrix( a.right * dag(b.right) * delta(inds(a.right)[3],inds(b.right)[3]) * delta(inds(a.right)[2],inds(b.right)[2]) ) - I(χ)
        frobeniusleft = norm(sqrt(sum(bioleft*bioleft')))
        frobeniusright = norm(sqrt(sum(bioright*bioright')))
        if printlog
            println(spaces, "• left biorthonormal = ", frobeniusleft)
            println(spaces, "• right biorthonormal = ", frobeniusright)
        end
        if frobeniusleft > tol || frobeniusright > tol
            errordetect = true
        end
    elseif flag == 2    # check center difference
        test_l = a.left * a.center * delta( inds(a.left)[3], inds(a.center)[1] )
        test_r = a.center * a.right * delta( inds(a.center)[2], inds(a.right)[1] )
        array_l = Array(test_l, inds(test_l))
        array_r = Array(test_r, inds(test_r))
        err_lr = array_l - array_r
        biocenterA = sqrt(sum(abs.(err_lr).^2))
        test_l = b.left * b.center * delta( inds(b.left)[3], inds(b.center)[1] )
        test_r = b.center * b.right * delta( inds(b.center)[2], inds(b.right)[1] )
        array_l = Array(test_l, inds(test_l))
        array_r = Array(test_r, inds(test_r))
        err_lr = array_l - array_r
        biocenterB = sqrt(sum(abs.(err_lr).^2))
        if printlog
            println(spaces, "• center left, right diff = [", biocenterA, ", ", biocenterB, "]")
        end
        if biocenterA > tol || biocenterB > tol
            errordetect = true
        end
    elseif flag == 3    # check gauge
        rholeft = a.center * dag(b.center) * delta(inds(a.center)[2],inds(b.center)[2])
        rhoright = dag(b.center) * a.center * delta(inds(a.center)[1],inds(b.center)[1])
        ratio_left = max_offdiag_min_diag(matrix(rholeft))
        ratio_right = max_offdiag_min_diag(matrix(rholeft))
        err = matrix(rholeft) - matrix(rhoright)
        mean = (matrix(rholeft) + matrix(rhoright))/2
        relative = sqrt(real(sum(err * err'))) / sqrt(real(sum(mean * mean')))
        if printlog
            println(spaces, "• diag left, right = [", ratio_left, ", ", ratio_right,"]")
            println(spaces, "• RDM left, right diff = ", relative)
        end
        if ratio_left < 1/tol || ratio_right < 1/tol || relative> tol
            errordetect = true
        end
    else
        println("ERROR@ check type = ",flag)
        errordetect = true
    end
    return errordetect
end
function biortho_test(a::iMPS_canonical, b::iMPS_canonical, flags::Vector{Int64} = [1,2]; printlog::Bool = true, indent::Int64 = 2, tol::Float64 = 1e-3)
    errorflag = false
    errordetect = false
    if printlog
        println(repeat(" ", indent), "Check biorthonormality of given iMPS pair.")
    end
    for flag in flags
        errordetect = biortho_test_flag(a, b, flag; printlog = printlog, indent = indent, tol = tol)
        errorflag = errorflag || errordetect
    end
    return errorflag
end
# canonical test
function canon_test_flag(a::iMPS_canonical, flag::Int64; printlog::Bool = true, indent::Int64 = 2, tol::Float64 = 1e-3)
    errordetect = false
    spaces = repeat(" ", indent)
    if flag ==1
        acleft = a.left*a.center * delta(inds(a.left)[3],inds(a.center)[1])
        acright = a.center*a.right * delta(inds(a.right)[1],inds(a.center)[2])
        acleft = array(acleft,inds(acleft))
        acright = array(acright,inds(acright))
        err = norm(acleft-acright)
        ortholeft = tensor_three(a.left,-10)
        orthoright = tensor_three(a.right,-20)
        if printlog
            println(spaces, "• orthonormal left=", ortholeft, ", orthonormal right=", orthoright)
            println(spaces, "• center left, right diff=",err)
        end
        if err>tol || ortholeft>tol || orthoright>tol
            errordetect = true
        end
    elseif flag == 2
        norm_c = norm(matrix(a.center))
        if printlog
            println(spaces, "• norm center=", norm_c)
        end
        if norm(norm_c-1)>tol
            errordetect = true
        end
    elseif flag == 3
        ratio = max_offdiag_min_diag(matrix(a.center))
        if printlog
            println(spaces, "• diag center=", ratio)
        end
        if norm(ratio)<1/tol
            errordetect = true
        end
    else
        println("ERROR@ check type = ",flag)
        errordetect = true
    end
    return errordetect
end
function canon_test(a::iMPS_canonical, flags::Vector{Int64} = [1,2]; printlog::Bool = true, indent::Int64 = 2, tol::Float64 = 1e-3)
    errorflag= false
    if printlog
        println(repeat(" ", indent), "Check canonicality of given iMPS.")
    end
    for flag in flags
        errordetect = canon_test_flag(a, flag; printlog = printlog, indent = indent, tol = tol)
        errorflag = errorflag || errordetect
    end
    return errorflag
end
# general check
function checkimps(a::iMPS_canonical, b::Union{iMPS_canonical, Nothing}=nothing; list::Vector{Int64} = [1,2], printlog::Bool = true, tol::Float64 = 1e-3, indent::Int64 = 0)
    errordetect = false
    if b ≠ nothing    # checking pair
        flags = [1,2,3]
        if size(a) == size(b)    #all(x -> x in flags, list)
            errordetect = biortho_test(a, b, list; printlog = printlog, indent = indent, tol = tol)
        elseif size(a) ≠ size(b)
            println("ERROR@ mps shape ",size(a),"≠",size(b))
            errordetect = true
        end
    else    # check single
        errordetect = canon_test(a, list; printlog = printlog, indent = indent, tol = tol)
    end
    return errordetect
end
