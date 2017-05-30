module FASTQDemultiplexerQC


import FASTQDemultiplexer: Output, Protocol, Barcode, InterpretedRecord,
    mergeoutput
using DataFrames
using Weave


type OutputQC{N} <: Output
    results::DataFrame
    poly::NTuple{N,DataFrame}
    nreads::Float64
    countpoly::Bool
end


function OutputQC{N}(protocol::Protocol{N};
                     outputdir::String = ".",
                     maxreads = Inf,
                     countpoly::Bool = true,
                     kwargs...)

    results = DataFrame(cellid = UInt[],
                        umiid = UInt[],
                        groupname = PooledDataArray(String,UInt8,0))
    poly = [
        DataFrame(A=UInt8[],C=UInt8[],T=UInt8[],G=UInt8[])
        for i in 1:N ]
    return OutputQC(results,(poly...),Float64(maxreads),countpoly)
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


function mergeoutput{N}(outputs::Vector{OutputQC{N}};
                        outputdir::String = ".",
                        kwargs...)
    if length(outputs) > 1
        results = vcat((o.results for o in outputs)...)
        poly = map(1:N) do i
            vcat((o.poly[i] for o in outputs)...)
        end
    else
        results = outputs[1].results
        poly = outputs[1].poly
    end

    mkpath(outputdir)

    report(results,poly,outputdir)
end

function report(results, poly, outputdir)

    args = Dict(:results=>results,
                :poly=>poly)

    weave(Pkg.dir("FASTQDemultiplexerQC","src","weave","QC.jmd"),
          args = args, out_path = outputdir)
end


function longeststreak(seq::Vector{UInt8},nuc::UInt8)
    len=0
    current=0
    for a in seq
        if a==nuc
            len=max(len,current)
            current+=1
        else
            current=0
        end
    end
    return len
end

end # module
