## Try ^ↀᴥↀ^

```bash
julia --project -e 'using Pkg; Pkg.activate(); Pkg.instantiate()' # install dependencies

julia --project liederkreis.jl mequetrefe.wav

mpv mequetrefe-loop.wav --loop
```

to specify minimal and maximal length of the loop (in seconds)

```bash
julia --project liederkreis.jl mequetrefe.wav --min 8 --max 16
```

to find more than one loop (after a certain point they may start to repeat)

```bash
julia --project liederkreis.jl mequetrefe.wav --min 8 --max 16 -n 4
```
