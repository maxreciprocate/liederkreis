using FFTW
using WAV
using Statistics
using Plots
# ■

import Base.match
using DataStructures

mutable struct Vertex
    children::Dict{UInt16, Vertex}
    key::Union{Int, Nothing}
    suffix::Union{Vertex, Nothing}
    isroot::Bool

    Vertex() = new(Dict{UInt16, Vertex}(), nothing, nothing, false)
end

function add!(vx::Vertex, string, key::Int)
    for char in string
        if !haskey(vx.children, char)
            vx.children[char] = Vertex()
        end

        vx = vx.children[char]
    end

    vx.key = key
end

function build!(root::Vertex)
    queue = Queue{Vertex}()
    enqueue!(queue, root)

    for vx in queue, (char, child) in vx.children
        child.suffix = if vx.isroot
            vx
        else
            search(vx.suffix, char)
        end

        enqueue!(queue, child)
    end
end

function search(vx::Vertex, char::UInt16)
    if haskey(vx.children, char)
        vx.children[char]
    elseif vx.isroot
        vx
    else
        search(vx.suffix, char)
    end
end
# ■

function match(source::Vector{UInt16}, markers)
    fsm = Vertex()
    fsm.isroot = true

    duplicates = Dict{UInt64, Vector{UInt16}}()
    marked_mapping = Dict{Vector{UInt16}, UInt64}()

    for (idx, marker) in enumerate(markers)
        if haskey(marked_mapping, marker)
            push!(duplicates[marked_mapping[marker]], idx)
        else
            add!(fsm, marker, idx)
            marked_mapping[marker] = idx
            duplicates[idx] = []
        end
    end

    build!(fsm)
    output = zeros(UInt64, length(markers))

    for (idx, char) in enumerate(source)
        fsm = search(fsm, char)
        vx = fsm

        while true
            if vx.key != nothing
                output[vx.key] = idx
            end

            vx.isroot && break
            vx = vx.suffix
        end
    end

    for (idx, markersids) in duplicates
        output[idx] == 0 && continue
        for jdx in markersids
            output[jdx] = output[idx]
        end
    end

    output
end
# ■

sourcefn = "roll.wav"
source, samplerate = wavread(sourcefn)

# one channel is sufficient
source = source[:, 1]

# ■

# packs in on second
θ = 16
# samples in one pack
δ = Int(samplerate) ÷ θ
# length of packs
ℓ = length(source) ÷ δ

maxfreqs = zeros(UInt16, ℓ)

frequencies = rfftfreq(δ, samplerate)
higherbound = findfirst(f -> f > 1000, frequencies) - 1

maxenergy = 1

# plot()

for idx = 1:ℓ
    transform = rfft(source[δ*(idx-1)+1:δ*idx+1])
    freqs = real.(transform)[1:higherbound]

    # plot!(frequencies[1:length(freqs)], freqs)

    # how does this behave with a zero signal?
    maxfreqs[idx] = argmax(freqs)
    maxenergy = max(maxenergy, maximum(freqs))
end

for jdx = 1:length(maxfreqs)
    if maxfreqs[jdx] < maxenergy / 16
        maxfreqs[jdx] = 0
    end
end

# maxfreqs[100:104]
# plot!()
# maxfreqs
# ■

markers = []

τ = θ
for from = 1:length(maxfreqs) - τ
    window = maxfreqs[from:from+τ]
    push!(markers, window)
end

markers
matches = match(maxfreqs, markers)

for (idx, m) in enumerate(matches)
    if (idx + first(matches)) < m + 1
        println(maxfreqs[idx:idx+τ-1], maxfreqs[m-τ:m-1])
    end
end

# ■

maxcorr = 0.0

# how many pack to take for a 4 seconds window
τ = 4 * θ
bestfrom = 0
bestto = 0

for fromidx = 1:τ:ℓ-τ, toidx = fromidx+τ:ℓ-τ
    X = @view maxfreqs[fromidx:fromidx+τ]
    Y = @view maxfreqs[toidx:toidx+τ]

    corr = cov(X, Y) / std(X) / std(Y)

    if corr > maxcorr
        maxcorr = corr
        bestfrom = fromidx
        bestto = toidx
    end
end

println(maxcorr)
maxfreqs[bestfrom:bestfrom+τ] |> println
maxfreqs[bestto:bestto+τ] |> println

looped = @view source[bestfrom * δ:bestto * δ]
length(looped) / samplerate

loop = vcat(looped, looped, looped, looped, looped)

wavwrite(loop, "loop.wav", Fs=samplerate, nbits=16, compression=WAVE_FORMAT_PCM)
