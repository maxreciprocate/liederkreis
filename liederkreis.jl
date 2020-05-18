using WAV
using MFCC: mfcc
using LinearAlgebra
import Base.MathConstants: φ

# packs in one second
const θ = 16

# size of windows
const τ = 4θ

# spacing
const δ = θ ÷ 4

const Timings = Dict{Tuple{Int, Int}, Float64}

function removefiles(dir::String)
    rm.(map(x -> "stash/" * x, split(read(`ls stash`, String))))
end

function gettimings(sourcefile::String; disk=true)
    source, samplerate = wavread(sourcefile)

    channel = ceil(Int, last(size(source)) * rand())
    source = source[Int(samplerate) >> 2 : end, channel]

    # samples in one pack
    β = Int(samplerate) ÷ θ

    # length of packs
    ℓ = length(source) ÷ β

    xs = first(mfcc(source, samplerate; wintime=1/θ, steptime=1/θ))
    xs = transpose(xs)

    timings = Timings()

    for fromidx = 8θ:4:ℓ-8θ, toidx = fromidx+θ:4:ℓ-8θ
        X = @view xs[:, fromidx:fromidx+τ]
        Y = @view xs[:, toidx:toidx+τ]

        if fromidx + 8θ < toidx
            timings[fromidx * β, toidx * β] = norm(X - Y)
        end
    end

    segments = sort(collect(zip(values(timings), keys(timings))))

    threshold = φ * first(first(segments))
    filter!(x -> first(first(x)) < threshold, segments)
    sort!(segments, by=x -> first(last(x)))

    quantity = 24
    lastfrom = 1
    for (distance, (from, to)) in segments
        # moderate spacing
        abs(lastfrom - from) < 2θ * β && continue
        lastfrom = from

        timing = "$from-$to-$(round(from / samplerate, digits=2))-$(round(to / samplerate, digits=2))"
        println(timing)

        disk && wavwrite(
            source[from : to], "stash/$sourcefile-$timing.wav",
                Fs=samplerate, nbits=16, compression=WAVE_FORMAT_PCM)

        quantity -= 1
        quantity <= 0 && break
    end
end

if !isfile(ARGS[1])
    println("$(ARGS[1]) is a non-existent file")
    exit(1)
end

gettimings(ARGS[1], disk=false)
