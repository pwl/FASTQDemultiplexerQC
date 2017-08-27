module FASTQDemultiplexerQC

import FASTQDemultiplexer:
    Output, Protocol, Barcode, Barcodes, InterpretedRecord,
    mergeoutput, cellidtoname, idtoname, umiid, cellid
using DataFrames
using Weave
using DataStructures

include("utils.jl")
include("polycounter.jl")

# TODO: why does it have to be mutable?
mutable struct OutputQC{N,C,U} <: Output
    results::DataFrame
    poly::Vector{PolyCounter}
    nreads::Float64
    countpoly::Bool
    protocol::Protocol{N,C,U}
end


function OutputQC{N,C,U}(protocol::Protocol{N,C,U},
                         barcodes::Barcodes;
                         outputdir::String = ".",
                         maxreads = Inf,
                         countpoly::Bool = true,
                         hamming::Bool = true,
                         polyspecs=[('A',1),('T',1)],
                         kwargs...)

    results = DataFrame(cellid = UInt[],
                        umiid = UInt[],
                        groupid = Int8[])
    poly = map(polyspecs) do specs
        PolyCounter(specs...,C)
    end
    return OutputQC{N,C,U}(results,poly,Float64(maxreads),countpoly,protocol)
end


function Base.write{N}(o::OutputQC{N},ir::InterpretedRecord{N})

    if o.nreads <= 0
        return
    else
        o.nreads-=1
    end

    push!(o.results,(cellid(ir),umiid(ir),ir.groupid))

    if o.countpoly && ir.groupid > -1
        for p in o.poly
            push!(p,ir)
        end
    end
end


function mergeoutput{N,C,U}(outputs::Vector{OutputQC{N,C,U}};
                            outputdir::String = ".",
                            hamming::Bool = false,
                            writestats::Bool = true,
                            writereport::Bool = true,
                            kwargs...)
    protocol = outputs[1].protocol
    if length(outputs) > 1
        results = vcat((o.results for o in outputs)...)

        poly = map(zip((o.poly for o in outputs)...)) do p
            merge(p...)
        end
    else
        results = first(outputs).results
        poly = first(outputs).poly
    end

    mkpath(outputdir)

    if writestats
        for p in poly
            write(outputdir,p)
        end
    end

    if writereport
        report(results,poly,outputdir,hamming,C,protocol.readnames)
    end
end

function report(results,poly,outputdir,hamming,celllen,readnames)

    args = Dict(:results=>results,
                :poly=>poly,
                :readnames=>readnames)

    if hamming
        cellids = unique(results[:cellid])
        distances = getminhammingdistance(cellids,celllen)
        args[:distances] =
            DataFrame(cellid=cellids,
                      mindistance=distances)
    end

    weave(joinpath(@__DIR__,"weave","QC.jmd"),
          args = args, out_path = outputdir)
end

end # module
