# try me ^ↀᴥↀ^

```bash
julia liederkreis.jl mequetrefe.wav

mpv mequetrefe-loop.wav --loop
```

to specify minimal and maximal length of the loop (in seconds)

```bash
julia liederkreis.jl mequetrefe.wav 8 16
```

to find more than one loop (after a certain point they may start to repeat)

```bash
julia liederkreis.jl mequetrefe.wav 8 16 4
```
