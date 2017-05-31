module FASTQDemultiplexerQC

include("utils.jl")

import FASTQDemultiplexer: Output, Protocol, Barcode, InterpretedRecord,
    mergeoutput
using DataFrames
using Weave


type OutputQC{N,C,U} <: Output
    results::DataFrame
    poly::NTuple{N,DataFrame}
    nreads::Float64
    countpoly::Bool
end


function OutputQC{N,C,U}(protocol::Protocol{N,C,U};
                         outputdir::String = ".",
                         maxreads = Inf,
                         countpoly::Bool = true,
                         hamming::Bool = true,
                         kwargs...)

    results = DataFrame(cellid = UInt[],
                        umiid = UInt[],
                        groupname = PooledDataArray(String,UInt8,0))
    poly = [
        DataFrame(A=UInt8[],C=UInt8[],T=UInt8[],G=UInt8[])
        for i in 1:N ]
    return OutputQC{N,C,U}(results,(poly...),Float64(maxreads),countpoly)
end


function Base.write{N}(o::OutputQC{N},ir::InterpretedRecord{N})

    if o.nreads <= 0
        return
    else
        o.nreads-=1
    end

    push!(o.results,(ir.cellid.val,ir.umiid.val,ir.groupname))
    if o.countpoly
        for i = 1:N
            seq = ir.records[i].data[ir.records[i].sequence]
            push!(o.poly[i],
                  map(UInt8['A','C','T','G']) do nuc
                  longeststreak(seq,nuc)
                  end)
        end
    end
end


function mergeoutput{N,C,U}(outputs::Vector{OutputQC{N,C,U}};
                            outputdir::String = ".",
                            hamming::Bool = false,
                            kwargs...)
    @time if length(outputs) > 1
        results = vcat((o.results for o in outputs)...)
        poly = map(1:N) do i
            vcat((o.poly[i] for o in outputs)...)
        end
    else
        results = first(outputs).results
        poly = first(outputs).poly
    end

    mkpath(outputdir)

    report(results,poly,outputdir,hamming,C)
end

function report(results,poly,outputdir,hamming,celllen)

    args = Dict(:results=>results,
                :poly=>poly)

    @time if hamming
        cellids = unique(results[:cellid])
        distances = getminhammingdistance(cellids,celllen)
        args[:distances] =
            DataFrame(cellid=cellids,
                      mindistance=distances)
    end

    weave(Pkg.dir("FASTQDemultiplexerQC","src","weave","QC.jmd"),
          args = args, out_path = outputdir)
end

end # module
