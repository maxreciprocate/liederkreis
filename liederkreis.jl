using FFTW
using WAV
using Statistics

function concat(sourcefn::String)
    xs, samplerate = wavread(sourcefn)
    xs = xs[:, 1]
    println("done reading...")
    maxfreqs = []

    ℓ = Int(samplerate) >> 4
    from = ℓ
    to = from + ℓ

    frequencies = rfftfreq(ℓ, samplerate)
    higherbound = findfirst(f -> f > 1e3, frequencies) - 1

    while to < length(xs)
        fs = real.(rfft(xs[from:to]))[1:higherbound]

        for idx = 1:length(fs)
            if fs[idx] < maximum(fs) / 4
                fs[idx] = 0
            end
        end

        push!(maxfreqs, argmax(fs))

        from = to
        to += ℓ
    end

    println("done analysis...")
    maxcorr = 0.0
    bfrom = 0
    bto = 0

    for from = 200:100:length(maxfreqs)-500, to = from+500:length(maxfreqs)-500
        X = maxfreqs[from:from+500]
        Y = maxfreqs[to:to+500]

        corr = cov(X, Y) / std(X) / std(Y)

        if corr > maxcorr
            maxcorr = corr
            bfrom = from
            bto = to
        end
    end

    println("done searching... $bfrom $bto $maxcorr")
    timed = xs[bfrom * ℓ:bto * ℓ]
    ys = vcat(timed, timed, timed)

    wavwrite(ys, "see.wav", Fs=samplerate, nbits=16, compression=WAVE_FORMAT_PCM)
end

concat(ARGS[1])
