# this is a small library to make it easier to checkpoint the
# internal state of long subroutines like in the legacy SLAB
# code, making it easier to build unit tests and subdivide
# large subroutines into smaller ones

# usage:
# at the top of the file you want to checkpoint in add
# include("../test/testutils/tracel.jl")
# 
# inside the function you want to checkpoint
# trace_log = ACTIVE_TRACE[]
# checkpoint!(trace_log, :label, object_you_want_to_save)
#
# when running the code to generate some data run
# trace_log = TraceLog("/path/to/checkpoint/data.jls")
# with_trace(trace_log) do
#  my_function()
# end
# close(trace_log)
#
# can then read the data to test against using
# my_data = readtrace("/path/to/checkpoint/data.jls")
# my_data will be a vector of named tuples (label, state)

using Serialization

mutable struct TraceLog
    io::IO
end

const ACTIVE_TRACE = Ref{Union{Nothing,TraceLog}}(nothing)

TraceLog(path::AbstractString) = TraceLog(open(path, "w"))
Base.close(log::TraceLog) = close(log.io)

function record!(log::TraceLog, label::Symbol, state)
    serialize(log.io, (label=label, state=state))
    flush(log.io)
    return nothing
end

function readtrace(path::AbstractString)
    open(path, "r") do io
        events = Any[]
        while !eof(io)
            push!(events, deserialize(io))
        end
        events
    end
end

function with_trace(f::Function, log::TraceLog)
    previous = ACTIVE_TRACE[]
    ACTIVE_TRACE[] = log
    try
        return f()
    finally
        ACTIVE_TRACE[] = previous
    end
end

function checkpoint!(log, label, snapshot)
    log === nothing && return nothing
    record!(log, label, snapshot)
end

