# cellid, streak
const CountIdentifier = Tuple{UInt,UInt8}

struct PolyCounter{NUC,R,C}
    streakcounter::Accumulator{CountIdentifier,Int}
    lastcounter::Accumulator{CountIdentifier,Int}
    nuccounter::Accumulator{CountIdentifier,Int}
end


function PolyCounter(nuc::Char,read::Int,bclen::Int)
    PolyCounter{nuc,read,bclen}(
        counter(CountIdentifier),
        counter(CountIdentifier),
        counter(CountIdentifier))
end


function Base.push!(pc::PolyCounter{NUC,R},ir::InterpretedRecord) where {NUC,R}
    seq = ir.records[R].data[ir.records[R].sequence]
    streak, last, count = longeststreak(seq,UInt8(NUC))
    push!(pc.streakcounter,(cellid(ir),UInt8(streak)))
    push!(pc.lastcounter,(cellid(ir),UInt8(last)))
    push!(pc.nuccounter,(cellid(ir),UInt8(count)))
end


function Base.merge!(pc::PolyCounter{NUC,R},
                     counters::PolyCounter{NUC,R}...) where {NUC,R}
    merge!(pc.streakcounter,(c.streakcounter for c in counters)...)
    merge!(pc.lastcounter,(c.lastcounter for c in counters)...)
    merge!(pc.nuccounter,(c.nuccounter for c in counters)...)
    return pc
end


function Base.merge(counters::PolyCounter{NUC,R,C}...) where {NUC,R,C}
    pc = PolyCounter(NUC,R,C)
    return merge!(pc,counters...)
end


function Base.write(name::String,pc::PolyCounter{NUC,R,C}) where {NUC,R,C}
    if isdir(name)
        filename = joinpath(name,"poly$NUC-read$R")
    else
        filename = name
    end
    bc(id) = idtoname(id,C)
    streaks = DataFrame(
        cellid = [bc(id) for ((id,streak),count) in pc.streakcounter],
        streak = [streak for ((id,streak),count) in pc.streakcounter],
        count  = [count  for ((id,streak),count) in pc.streakcounter])
    writetable(filename*"-streaks.csv",unstack(streaks,:cellid,:streak,:count))
    last = DataFrame(
        cellid = [bc(id) for ((id,last),count)  in pc.lastcounter],
        last   = [last   for ((id,last),count)  in pc.lastcounter],
        count  = [count  for ((id,last),count)  in pc.lastcounter])
    writetable(filename*"-last.csv",unstack(last,:cellid,:last,:count))
    nuc = DataFrame(
        cellid = [bc(id) for ((id,nuc),count)   in pc.nuccounter],
        nuc    = [nuc    for ((id,nuc),count)   in pc.nuccounter],
        count  = [count  for ((id,nuc),count)   in pc.nuccounter])
    writetable(filename*"-nuc.csv",unstack(nuc,:cellid,:nuc,:count),nastring="0")
end
