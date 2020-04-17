using Corona
using Plots
using FileIO, JLD2
using TimeSeries
using LaTeXStrings
using Formatting

### Load configuration
dataconfig = Corona.DataConfig();
include("../config/data.jl");
if isfile(joinpath(Base.@__DIR__,"../config/paths.jl"))
    include("../config/paths.jl")
end;

## Decide region and title
region = "New York";
title = "New York City";

## Select region and paths
filename = Corona.build_filename(dataconfig,region);
figdir = joinpath("figs",filename);
mkpath(figdir);
datapath = joinpath(dataconfig.results_path,filename,"da.jld2");


## Load data
da = load(datapath,"da");
lastday = timestamp(da.data)[meta(da.data)["last_day_idxs"]];

## Compute residual
J = Corona.residual(da, relative=true);

## Plot results
plot(da,region)
title!(format(title*", J = {:.3f} %",J*100))

##
savefig(joinpath(figdir,"da.pdf"))

## Plot results in log scale
plot(da,region,yaxis=:log10,legend=:bottomright)
title!(format(title*", J = {:.3f} %",J*100))

##
savefig(joinpath(figdir,"da_log.pdf"))

## Plot model parameters
ptime = Date(2020,3,24):Day(1):lastday;
plot(da.p.β,color=:red, lw=3, legend=false)
scatter!(da.p.β,color=:red, lw=3)
ylabel!("beta")
title!(format(title*", J = {:.3f} %",J*100))

##
savefig(joinpath(figdir,"beta.pdf"))

## Plot model parameters
σ=da.p.β./(da.p.γ .+ da.p.δ);
plot(σ,color=:orange, lw=3, legend=false)
scatter!(σ,color=:orange, lw=3)
ylabel!("sigma")
title!(format(title*", J = {:.3f} %",J*100))

##
savefig(joinpath(figdir,"sigma.pdf"))
