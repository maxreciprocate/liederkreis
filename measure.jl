using LinearAlgebra
using MFCC
using WAV
using Plots
# ■

songname = "newworld.wav"
source, samplerate = wavread(songname)
source = source[:, 1]
# ■

# packs in on second
θ = 16
# samples in one pack
δ = Int(samplerate) ÷ θ
# length of packs
ℓ = length(source) ÷ δ

xs, _, _  = mfcc(source, samplerate; wintime=1/θ, steptime=1/θ)
xs = xs'

# ■
τ = 4θ
distances = Dict{Tuple{Int, Int}, Float64}()

for fromidx = 2τ:τ:ℓ-τ, toidx = fromidx+τ:τ:ℓ-2τ
    X = @view xs[:, fromidx:fromidx+τ]
    Y = @view xs[:, toidx:toidx+τ]

    distances[fromidx * δ, toidx * δ] = norm(X - Y)
end

matches = sort(collect(zip(values(distances), keys(distances))))

rm.(map(x -> "loops/" * x, split(read(`ls loops`, String))))

for idx = 1:10
    distance, (from, to) = matches[idx]

    loop = vcat(source[to - 2τ * δ : to], source[from : from + 2τ * δ])

    wavwrite(
        loop, "loops/$(first(split(name, '.')))-$from-$to-$(from ÷ samplerate)-$(to ÷ samplerate).wav",
             Fs=samplerate, nbits=16, compression=WAVE_FORMAT_PCM)
end

_, (bestfrom, bestto) = first(matches)

looped = source[bestfrom : bestto]
loop = vcat(looped, looped, looped, looped, looped)

wavwrite(loop, "loops/$name-loop.wav", Fs=samplerate, nbits=16, compression=WAVE_FORMAT_PCM)
