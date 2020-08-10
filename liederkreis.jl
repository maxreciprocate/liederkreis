#!/usr/bin/env julia
using MFCC: mfcc
using LinearAlgebra: norm
using DataStructures
using Clustering
using WAV

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

function main(sourcefile::String, lowerbound=10, upperbound=100, quantity=1)
    println("reading the file")
    source, sr = wavread(sourcefile)
    source = source[:, 1]

    println("computing mfcc")
    cepstrum = Float32.(first(mfcc(source, sr; wintime=1/s, steptime=1/s, numcep=numcep)))

    # since coefficients are strided, restride gives ~2x speedup
    cepstrum = restride(cepstrum)

    # number of samples in a single pack
    numsamples = Int(sr) ÷ s

    # length of the track (measured in the number of packs)
    numpacks = length(source) ÷ numsamples

    println("extracting seaming points")
    seams = extract(cepstrum, lowerbound * s, upperbound * s, numpacks)

    matrix = Array{Float32}(undef, 2, length(seams))

    for idx = 1:length(seams)
        matrix[1, idx] = seams[idx].startidx
        matrix[2, idx] = seams[idx].endidx
    end

    println("trimming results")
    indices = kmeans(matrix, quantity) |> assignments

    sourcename = split(sourcefile, '/') |> last |> s -> split(s, '.') |> first

    for idx = 1:quantity
        seam = minimum(seams[findall(isequal(idx), indices)])

        loop = source[seam.startidx * numsamples : seam.endidx * numsamples]

        loopname = if quantity == 1
            "$sourcename-loop.wav"
        else
            "$sourcename-loop-$idx.wav"
        end

        wavwrite(loop, loopname, Fs=sr, nbits=16, compression=WAVE_FORMAT_PCM)
    end
end

if length(ARGS) < 1 || length(ARGS) > 4
    println("usage: julia liederkreis.jl <track.wav> &optional: <min-seconds> <max-seconds> <the number of loops>")
    exit(1)
end

if length(ARGS) == 1
    main(ARGS[1])
else
    main(ARGS[1], map(arg -> parse(Int, arg), ARGS[2:end])...)
end
