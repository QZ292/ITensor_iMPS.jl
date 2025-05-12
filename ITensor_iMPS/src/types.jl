mutable struct iMPS_canonical
    left::ITensor
    center::ITensor
    right::ITensor
    iMPS_canonical() = new(ITensor(0),ITensor(0),ITensor(0))
end
# redefine Base.size
function sizemps(F::iMPS_canonical)
    @assert size(F.left) == size(F.right)
    @assert size(F.center)[1] == size(F.center)[2]
    @assert size(F.left)[1] == size(F.left)[3]
    @assert size(F.center)[1] == size(F.left)[1]
    return (size(F.center)[1],size(F.left)[2])
end
Base.size(F::iMPS_canonical) = sizemps(F::iMPS_canonical)
# redefine Base.deepcopy
function deepcopymps(F::iMPS_canonical)
    G = iMPS_canonical()
    G.left = deepcopy(F.left)
    G.center = deepcopy(F.center)
    G.right = deepcopy(F.right)
    return G
end
Base.deepcopy(F::iMPS_canonical) = deepcopymps(F::iMPS_canonical)
# redefine Base.show
function capture_display_output(obj)
    buffer = IOBuffer()
    custom_display = TextDisplay(buffer)
    pushdisplay(custom_display)
    local captured_string::String = ""
    try
        display(obj)
        captured_string = String(take!(buffer))
    catch e
        rethrow(e)
    finally
        popdisplay(custom_display)
        close(buffer)
    end
    return captured_string
end
function Base.show(io::IO, mime::MIME"text/plain", a::iMPS_canonical)
    oleft = capture_display_output(a.left)
    oleft = split(oleft, '\n'; limit=2)
    ocenter = capture_display_output(a.center)
    ocenter = split(ocenter, '\n'; limit=2)
    oright = capture_display_output(a.right)
    oright = split(oright, '\n'; limit=2)
    omps = string("iMPS_canonical\n left  -> ",oleft[1],"\n          ",oleft[2],
        " center-> ",ocenter[1],"\n          ",ocenter[2],
        " right -> ",oright[1],"\n          ",oright[2]
    )
    println(io, omps)
    return nothing
end
function Base.show(io::IO, a::iMPS_canonical)
    print(io, "iMPS_canonical{\n")
    print(io, "> Left ", a.left)
    print(io, "\n> Center ", a.center)
    print(io, "\n> Right ", a.right,"}")
end
