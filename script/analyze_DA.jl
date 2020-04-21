using Corona
using LinearAlgebra
using Plots
using FileIO, JLD2
using TimeSeries
using LaTeXStrings
using Formatting
using CSV
using FFTW

## Decide region and title
region = "New York";
title = "New York City";

## What to do
display_fig = true;
save_fig = true;
show_window = false;
show_actions = true;

### Load configuration
dataconfig = Corona.DataConfig();
include("../config/data.jl");
if isfile(joinpath(Base.@__DIR__,"../config/paths.jl"))
    include("../config/paths.jl")
end;

## Select paths
fname = Corona.build_filename(dataconfig,region);
figdir = joinpath(dataconfig.paths["figs"],fname);
mkpath(figdir);
fig_ext = ["pdf", "png"];

## Load actions
act = CSV.read(joinpath(dataconfig.paths["raw"],"GOV",fname,"actions.csv"));

## Load data
data, da = load(joinpath(dataconfig.paths["results"],fname,"da.jld2"), "data", "da");

## Compute residual
J_min = Corona.residual(Corona.baseline(da), relative=true);

## Set title
title_str = format(title*", J = {:.3f} %",J_min*100);

## Plot results
let plt = Plots.Plot()
        datarange = data.outbreakdate:Day(1):data.lastdate;
        scatter!(plt, DateTime.(timestamp(data.cases[datarange])), values(data.cases.Confirmed[datarange]),
                    color=:orange, label="C (data)", Atickfontsize=12, legend=:topleft,
                    minorticks=true);
        scatter!(plt, DateTime.(timestamp(data.cases[datarange])), values(data.cases.Deaths[datarange]),
                    color=:black, label="D (data)");

        plot!(plt, da.u.I .+ da.u.R .+ da.u.D, color=:orange, lw=3, label="C");
        plot!(plt, da.u.I, color=:red, lw=3);
        plot!(plt, da.u.R, color=:green, lw=3);
        plot!(plt, da.u.D, color=:black, lw=3);
        title!(plt, title_str);

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action] , lw=3)
            end
        end

        show_window && plot!(plt, ylims(plt)[2] .* da.σ.S, fillrange=ylims(plt)[1],
                                α=0.2, color=:grey, label="control");

        display_fig && display(plt)
        if save_fig
            for ext in fig_ext
                savefig(plt, joinpath(figdir,"da."*ext))
            end
        end

        plt
end

## Plot results in log scale
let plt = Plots.Plot()
        datarange = data.outbreakdate:Day(1):data.lastdate;
        scatter!(plt, DateTime.(timestamp(data.cases[datarange])), values(data.cases.Confirmed[datarange]) .+ 1,
                    color=:orange, label="C (data)", Atickfontsize=12, legend=:none, yaxis=:log10);
        scatter!(plt, DateTime.(timestamp(data.cases[datarange])), values(data.cases.Deaths[datarange]) .+ 1,
                    color=:black, label="D (data)");

        plot!(plt, da.u.I .+ da.u.R .+ da.u.D .+ 1.0, color=:orange, lw=3, label="C");
        plot!(plt, da.u.I .+ 1.0, color=:red, lw=3);
        plot!(plt, da.u.R .+ 1.0, color=:green, lw=3);
        plot!(plt, da.u.D .+ 1.0, color=:black, lw=3);
        title!(plt, title_str);

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action] , lw=3)
            end
        end

        show_window && plot!(plt, ylims(plt)[2] .* da.σ.S .+ ylims(plt)[1], fillrange=ylims(plt)[1],
                                α=0.2, color=:grey, label="control");

        display_fig && display(plt)
        if save_fig
            for ext in fig_ext
                savefig(plt, joinpath(figdir,"da_log."*ext))
            end
        end

        plt
end

## Plot β
let plt = Plots.Plot()
        plot!(plt, da.p.β, color=:orange, lw=3, label="\\beta");

        show_window && plot!(plt, ylims(plt)[2] .* da.σ.S, fillrange=ylims(plt)[1],
                                α=0.2, color=:grey, label="control", legend=:topleft);

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action] , lw=3)
            end
        end

        title!(plt, title_str);

        display_fig && display(plt)
        if save_fig
            for ext in fig_ext
                savefig(plt, joinpath(figdir,"da_beta."*ext))
            end
        end

        plt
end

## Plot δ and γ
let plt = Plots.Plot()
        plot!(plt, da.p.γ, color=:blue, lw=3, label="\\gamma", legend=:right);
        plot!(plt, da.p.δ, color=:black, lw=3, label="\\delta");

        show_window && plot!(plt, ylims(plt)[2] .* da.σ.S, fillrange=ylims(plt)[1],
                                α=0.2, color=:grey, label="control");

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action] , lw=3)
            end
        end

        title!(plt, title_str);

        display_fig && display(plt)
        if save_fig
            for ext in fig_ext
                savefig(plt, joinpath(figdir,"da_delta_gamma."*ext))
            end
        end

        plt
end

## Plot sigma
let plt = Plots.Plot()
        plot!(plt, da.p.β./(da.p.γ .+ da.p.δ), color=:red, lw=3, label="\\sigma" , legend=:topright)
        ylims!(0,3)

        show_window && plot!(plt, ylims(plt)[2] .* da.σ.S, fillrange=ylims(plt)[1],
                                α=0.2, color=:grey, label="control");

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action] , lw=3)
            end
        end

        title!(plt, title_str);

        display_fig && display(plt)
        if save_fig
            for ext in fig_ext
                savefig(plt, joinpath(figdir,"da_recovery_rate."*ext))
            end
        end

        plt
end;

## Compute eigenvalues of A and extract A22
Λ = TimeArray(timestamp(da.u), complex(zeros(size(da.u))), [:λ₁, :λ₂, :λ₃, :λ₄]);
A22 = TimeArray(timestamp(da.u), zeros(size(da.u,1)), [:A22]);

for i=1:size(da.u,1)
        A = Corona.dfdu(values(da.u)[i,:], values(da.p)[i,:])
        values(A22)[i] = A[2,2]
        values(Λ)[i,:] = eigvals(A)
end;

## Plot eigenvalues of A
let plt = Plots.Plot()
        plot!(plt, maximum(real.(Λ), dims=2), label="\\lambda max", lw=3, color=:blue)
        ylims!(plt, -.01, 0.4)

        show_window && plot!(plt, ylims(plt)[2] .* da.σ.S, fillrange=ylims(plt)[1],
                                α=0.2, color=:grey, label="control");

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action] , lw=3)
            end
        end

        title!(plt, title_str);

        display_fig && display(plt)
        if save_fig
            for ext in fig_ext
                savefig(plt, joinpath(figdir,"da_lambda_max."*ext))
            end
        end

        plt
end;

## Plot eigenvalues of A22
let plt = Plots.Plot()
        plot!(plt, A22, label="A22", lw=3, color=:blue)
        # ylims!(plt, -.01, 0.4)

        show_window && plot!(plt, ylims(plt)[2] .* da.σ.S, fillrange=ylims(plt)[1],
                                α=0.2, color=:grey, label="control");

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action] , lw=3)
            end
        end

        title!(plt, title_str);

        display_fig && display(plt)
        if save_fig
            for ext in fig_ext
                savefig(plt, joinpath(figdir,"da_A22."*ext))
            end
        end

        plt
end;

## Compute FFT of A22
freq = 24*(0:(size(A22,1)/2))/size(A22,1);
P_A22 = abs.(fft(values(A22)));

let plt = Plots.Plot()
        sticks!(plt, freq[7:60].^(-1), P_A22[7:60], lw=3, label="FFT(A22)",
                xaxis="T (days)", legend=:topleft, color=:blue)
        # ylims!(plt, -.01, 0.4)

        title!(plt, title_str);

        display_fig && display(plt)
        if save_fig
            for ext in fig_ext
                savefig(plt, joinpath(figdir,"da_FFT_A22."*ext))
            end
        end

        plt
end;
