# Compute the longest streak of `nuc` in vector `seq`
function longeststreak(seq::Vector{UInt8},nuc::UInt8)
    len=0
    current=0
    last=0
    count=0
    for i in linearindices(seq)
        if seq[i]==nuc
            count+=1
            current+=1
            if current > len
                len = current
                last = i
            end
        else
            current=0
        end
    end
    return len, last, count
end

# Compute a Hamming distance between `id1` and `id2`, assuming that
# both are of length `len`
function hamming(id1::UInt,id2::UInt,len)
    id1-=1
    id2-=1
    dist::UInt8 = 0
    for n in 1:len
        if (id1&0b11) $ (id2&0b11) != 0
            dist+=1
        end
        id1>>=2
        id2>>=2
    end
    return dist
end


function getminhammingdistance(codes::AbstractVector{UInt},N)
    nearest = zeros(UInt8,length(codes))
    @inbounds for i1 in 1:length(codes)
        nearest[i1] = N       # start with a maximal distance
        for i2 in 1:(i1-1)
            dist = hamming(codes[i1],codes[i2],N)
            nearest[i1] = min(nearest[i1],dist)
            # the minimum has been reached
            if dist == 1
                continue
            end
        end
    end
    return nearest
end
