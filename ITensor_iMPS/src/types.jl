mutable struct iMPS_canonical
    left::ITensor
    center::ITensor
    right::ITensor
    iMPS_canonical() = new(ITensor(0),ITensor(0),ITensor(0))
end
function sizemps(F::iMPS_canonical)
    @assert size(F.left) == size(F.right)
    @assert size(F.center)[1] == size(F.center)[2]
    @assert size(F.left)[1] == size(F.left)[3]
    @assert size(F.center)[1] == size(F.left)[1]
    return (size(F.center)[1],size(F.left)[2])
end
Base.size(F::iMPS_canonical) = sizemps(F::iMPS_canonical)
function deepcopymps(F::iMPS_canonical)
    G = iMPS_canonical()
    G.left = deepcopy(F.left)
    G.center = deepcopy(F.center)
    G.right = deepcopy(F.right)
    return G
end
Base.deepcopy(F::iMPS_canonical) = deepcopymps(F::iMPS_canonical)
