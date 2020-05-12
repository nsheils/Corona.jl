# Corona.jl 
Corona.jl performs the assimilation of epidemiological data using an augmented SIR model.

## Getting started
You can get started with Corona.jl in a few steps:

1. [Install Corona.jl](#installing-coronajl)
2. [Install the dependecies](#installing-dependencies)
3. [Prepare data](#preparing-data)
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

### Preparing data

1. Provide your preferred data paths in `config/paths.jl` as shown in `config/example-paths.jl`.
2. Provide information on the data sources in `config/data.jl` following the structure of `config/example-data.jl`.

For a quick start, you can use the fake data in the files `raw/native/example-Laputa.csv` and `raw/example-population.csv`, by copying them to `raw/native/Laputa.csv` and `raw/population.csv`, respectively.

### Running the data assimilation
Move to the installation folder, modify the Julia script `script/run_DA.jl` as desired and execute

  ```julia
    include("script/run_DA.jl")
  ```

from the Julia REPL.

### Analyzing the results 
Move to the installation folder (if you are not already there), modify the Julia script `script/analyze_DA.jl` as desired and execute

  ```julia
    include("script/analyze_DA.jl")
  ```

from the Julia REPL.

## License
This package is free to use for noncommercial purposes and for commercial purposes during a trial period under the terms of the [Prosperity Public License](LICENSE.md).

