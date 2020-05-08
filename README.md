# Corona.jl 
Corona.jl performs the assimilation of epidemiological data using an augmented SIR model.

## Getting started
You can get started with Corona.jl in a few steps:

1. [Install Corona.jl](#installing-coronajl)
2. [Install the dependecies](#installing-dependencies)
3. [Define data directories](#defining-data-directories)
3. [Run the data assimilation](#running-the-data-assimilation)
4. [Analyze the results](#analyzing-the-results)

### Installing Corona.jl
1. Clone Corona.jl (replace `/path/to/julia/modules` with your desired path):

  ```bash
  git clone https://github.com/sesterhenn/Corona.jl /path/to/julia/modules/Corona
  ```

2. Add to `~/.julia/config/startup.jl` the following line (unless `/path/to/julia/modules` is already added to your Julia `LOAD_PATH` variable):

  ```julia
  push!(LOAD_PATH,"/path/to/julia/modules/Corona")
  ```

### Installing dependencies
Execute
  
  ```julia
    ] add Dates TimeSeries DataStructures DifferentialEquations LinearAlgebra Flux Interpolations FFTW FileIO JLD2 DataFrames Unicode CSV Formatting Plots LaTeXStrings
  ```

from the Julia REPL.

### Defining data directories
For a quick start, use `config/example-paths.jl` in your Corona.jl installation folder.

### Running the data assimilation
Move to the installation folder and execute

  ```julia
    include("scripts/run_DA.jl")
  ```

from the Julia REPL.

### Analyzing the results 
Move to the installation folder (if you are not already there) and execute

  ```julia
    include("scripts/run_DA.jl")
  ```

from the Julia REPL.

## License
This package is free to use for noncommercial purposes and for commercial purposes during a trial period under the terms of the [Prosperity Public License](LICENSE.md).

