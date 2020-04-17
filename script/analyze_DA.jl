using Corona
using Plots
using FileIO, JLD2
using TimeSeries
using LaTeXStrings
using Formatting

### Load configuration
dataconfig = Corona.DataConfig();
include("../config/data.jl");

## Select region and paths
region = "New York";
filename = Corona.build_filename(dataconfig,region);
figdir = joinpath("figs",filename);
mkpath(figdir);
datapath = joinpath("results",filename,"da.jld2");

## Load data
da = load(datapath,"da");
lastday = timestamp(da.data)[meta(da.data)["last_day_idxs"]];

## Compute residual
J = Corona.residual(da, relative=true);

## Plot results
plot(da,region)
title!(format("King County, Washington, J = {:.3f} %",J*100))
savefig(joinpath(figdir,"da.pdf"))

## Plot results in log scale
plot(da,region,yaxis=:log10,legend=:bottomright)
title!(format("King County, Washington, J = {:.3f} %",J*100))
savefig(joinpath(figdir,"da_log.pdf"))

## Plot model parameters
ptime = Date(2020,3,24):Day(1):lastday;
plot(da.p.β[ptime],color=:red, lw=3, legend=false)
scatter!(da.p.β[ptime],color=:red, lw=3)
ylabel!("beta")
title!(format("King County, Washington, J = {:.3f} %",J*100))
savefig(joinpath(figdir,"beta.pdf"))

## Plot model parameters
σ=da.p.β./(da.p.γ .+ da.p.δ);
plot(σ[ptime],color=:orange, lw=3, legend=false)
scatter!(σ[ptime],color=:orange, lw=3)
ylabel!("sigma")
title!(format("King County, Washington, J = {:.3f} %",J*100))
savefig(joinpath(figdir,"sigma.pdf"))
