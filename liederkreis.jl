#!/usr/bin/env julia
using MFCC: mfcc
using LinearAlgebra: norm
using DataStructures
using Clustering
using ArgParse
using WAV
using Crayons

struct SeamPoint
    startidx::Int
    endidx::Int
    loss::Float32
end

import Base.isless
isless(a::SeamPoint, b::SeamPoint) = a.loss < b.loss

function restride(input::AbstractArray{T, 2}) where T
    out = Vector{T}(undef, length(input))
    m, n = size(input)

    for idx = 1:m, jdx = 1:n
        @inbounds out[(idx-1) * n + jdx] = input[idx, jdx]
    end

    return out
end

# the number of packs per one second
const s = 32

# the number of MFC coefficients
const numcep = 30

# length of the pack (measured in samples)
const τ = 1

function extract(cepstrum::Vector{Float32}, lowerbound::Int, upperbound::Int, numpacks::Int)::Vector{SeamPoint}
    capacity = 128
    heap = BinaryMaxHeap{SeamPoint}([SeamPoint(1, 1, Inf)])

    for a = 1:numpacks - lowerbound - τ, b = a + lowerbound : min(a + upperbound - τ, numpacks - τ)
        pack₁ = @view cepstrum[(a - 1) * numcep + 1 : (a + τ) * numcep]
        pack₂ = @view cepstrum[(b - 1) * numcep + 1 : (b + τ) * numcep]

        loss = norm(pack₁ - pack₂, 1)

        if loss < top(heap).loss
            length(heap) > capacity && pop!(heap)

            push!(heap, SeamPoint(a, b, loss))
        end
    end

    return extract_all!(heap)
end

function main(args)
    rules = ArgParseSettings("Extract precise loops from the audio track")

    @add_arg_table! rules begin
        "audio-path"
        arg_type = String
        required = true
        help = "path to the source audio file"

        "--min-seconds", "--min"
        nargs = '?'
        arg_type = Int
        default = 10
        help = "minimum length of the loop"

        "--max-seconds", "--max"
        nargs = '?'
        arg_type = Int
        default = 300
        help = "minimum length of the loop"

        "--quantity", "-n"
        nargs = '?'
        arg_type = Int
        default = 1
        help = "the number of loops"

        "--quiet", "-q"
        action = :store_true
        help = "silence the output"
    end

    settings = parse_args(args, rules)

    !settings["quiet"] && println(crayon"red", "reading the file...")
    source, sr = wavread(settings["audio-path"])
    source = source[:, 1]

    !settings["quiet"] && println(crayon"yellow", "computing mfcc...")
    cepstrum = Float32.(first(mfcc(source, sr; wintime=1/s, steptime=1/s, numcep=numcep)))

    # since coefficients are strided, restride gives ~2x speedup
    cepstrum = restride(cepstrum)

    # number of samples in a single pack
    numsamples = Int(sr) ÷ s

    # length of the track (measured in the number of packs)
    numpacks = length(source) ÷ numsamples

    !settings["quiet"] && println(crayon"green", "extracting seaming points...")
    seams = extract(cepstrum, settings["min-seconds"] * s, settings["max-seconds"] * s, numpacks)

    matrix = Array{Float32}(undef, 2, length(seams))

    for idx = 1:length(seams)
        matrix[1, idx] = seams[idx].startidx
        matrix[2, idx] = seams[idx].endidx
    end

    !settings["quiet"] && println(crayon"blue", "trimming results...")
    indices = kmeans(matrix, settings["quantity"]) |> assignments

    sourcename = split(settings["audio-path"], '/') |> last |> s -> split(s, '.') |> first

    for idx = 1:settings["quantity"]
        seam = minimum(seams[findall(isequal(idx), indices)])

        loop = source[seam.startidx * numsamples : seam.endidx * numsamples]

        loopname = if settings["quantity"] == 1
            "$sourcename-loop.wav"
        else
            "$sourcename-loop-$idx.wav"
        end

        wavwrite(loop, loopname, Fs=sr, nbits=16, compression=WAVE_FORMAT_PCM)
        !settings["quiet"] && println(crayon"magenta", "congratulations, now run\nmpv $loopname --loop")
    end
end

main(ARGS)
