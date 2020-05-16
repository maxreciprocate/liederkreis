using WAV
using MFCC
using LinearAlgebra
# ■

source, samplerate = wavread("ravel-jeux-deau.wav")

channel = ceil(Int, last(size(source)) * rand())

source = source[Int(samplerate) >> 2 : end, channel]
# ■

# packs in one second
θ = 16

# spacing
τ = 1θ

# samples in one pack
δ = Int(samplerate) ÷ θ

# length of packs
ℓ = length(source) ÷ δ

xs = first(mfcc(source, samplerate; wintime=1/θ, steptime=1/θ))
xs = transpose(xs)

# ■

short = Dict{Tuple{Int, Int}, Float64}()
large = Dict{Tuple{Int, Int}, Float64}()

for fromidx = 2τ:ℓ-τ, toidx = fromidx+τ:ℓ-2τ
    X = @view xs[:, fromidx:fromidx+τ]
    Y = @view xs[:, toidx:toidx+τ]

    if fromidx + 8τ > toidx
        short[fromidx * δ, toidx * δ] = norm(X - Y)
    else
        large[fromidx * δ, toidx * δ] = norm(X - Y)
    end
end

rm.(map(x -> "loops/" * x, split(read(`ls loops`, String))))

function flesh(matches::Dict{Tuple{Int, Int}, Float64}, tag::String, quantity::Int)
    println("[fleshing $tag]:")
    segments = sort(collect(zip(values(matches), keys(matches))))

    lastto = 1
    for (distance, (from, to)) in segments
        # moderate spacing
        lastto + τ * δ > to && continue
        lastto = to

        timing = "$from-$to-$(from / samplerate)-$(to / samplerate)"
        println(timing)

        wavwrite(
            source[from : to], "loops/$tag-$timing.wav",
                Fs=samplerate, nbits=16, compression=WAVE_FORMAT_PCM)

        quantity -= 1
        quantity <= 0 && break
    end
end

flesh(short, "short", 12)
flesh(large, "large", 4)
