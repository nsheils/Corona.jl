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
region = "Laputa";
plt_title = "Laputa";

## What to do
display_fig = true;
save_fig = true;
show_window = false;
show_actions = true;

### Load configuration
dataconfig = Corona.DataConfig();
include(joinpath(pwd(),"config/data.jl"));
if isfile(joinpath(pwd(),"config/paths.jl"))
    include(joinpath(pwd(),"config/paths.jl"))
end;

## Select paths
fname = Corona.build_filename(dataconfig,region);
figdir = joinpath(dataconfig.paths["figs"],fname);
resdir = joinpath(dataconfig.paths["results"],fname);
mkpath(figdir);
fig_ext = ["pdf", "png"];

## Load actions
act = CSV.read(joinpath(dataconfig.paths["raw"],"GOV",fname,"actions.csv"));
act_colors =  [:darkolivegreen, :palegreen2, :khaki, :darkslategray4, :sienna4];

## Load data
data, da = load(joinpath(resdir,"da.jld2"), "data", "da");

## Compute residual
J_min = Corona.residual(Corona.baseline(da), relative=true);

## Set title
# title_str = format(plt_title*", J = {:.3f} %",J_min*100);
title_str = plt_title*", $(data.lastdate)";
# title_str = plt_title;

## Determine datarange
datarange = data.outbreakdate:Day(1):data.lastdate;
datarange_dt = DateTime(data.outbreakdate):step(da.time):DateTime(data.lastdate);

## Determine result range
resrange = data.outbreakdate:Day(1):data.lastdate+Day(14);
resrange_dt = DateTime(resrange.start):step(da.time):DateTime(resrange.stop);

## Plot results
let plt = Plots.Plot()

        scatter!(plt, DateTime.(timestamp(data.cases[datarange])), values(data.cases.Confirmed[datarange]),
                    color=:orange, label="C (Confirmed) - data", Atickfontsize=12, legend=:topleft,
                    minorticks=true);
        scatter!(plt, DateTime.(timestamp(data.cases[datarange])), values(data.cases.Deaths[datarange]),
                    color=:black, label="D (Deaths) - data");

        plot!(plt, da.u.I[resrange_dt] .+ da.u.R[resrange_dt] .+ da.u.D[resrange_dt], color=:orange, lw=3, label="C (Confirmed)");
        plot!(plt, da.u.I[resrange_dt], color=:red, lw=3, label="I (Infectives)");
        plot!(plt, da.u.R[resrange_dt], color=:green, lw=3, label="R (Recovered)");
        plot!(plt, da.u.D[resrange_dt], color=:black, lw=3, label="D (Deaths)");
        title!(plt, title_str);

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action] , lw=2, ls=:dash, color = act_colors[i])
            end
        end

        show_window && plot!(plt, ylims(plt)[2] .* da.σ.Confirmed, fillrange=ylims(plt)[1],
                                α=0.2, color=:grey, label="control");

        display_fig && display(plt)
        if save_fig
            for ext in fig_ext
                savefig(plt, joinpath(figdir,"da."*ext))
            end
        end

        plt
end;

## Plot results in log scale
let plt = Plots.Plot()

        scatter!(plt, DateTime.(timestamp(data.cases[datarange])), values(data.cases.Confirmed[datarange]) .+ 1,
                    color=:orange, label="C (data)", Atickfontsize=12, legend=:none, yaxis=:log10);
        scatter!(plt, DateTime.(timestamp(data.cases[datarange])), values(data.cases.Deaths[datarange]) .+ 1,
                    color=:black, label="D (data)");

        plot!(plt, da.u.I[datarange_dt] .+ da.u.R[datarange_dt] .+ da.u.D[datarange_dt] .+ 1.0, color=:orange, lw=3, label="C");
        plot!(plt, da.u.I[datarange_dt] .+ 1.0, color=:red, lw=3);
        plot!(plt, da.u.R[datarange_dt] .+ 1.0, color=:green, lw=3);
        plot!(plt, da.u.D[datarange_dt] .+ 1.0, color=:black, lw=3);
        title!(plt, title_str);

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action] , lw=2, ls=:dash, color = act_colors[i])
            end
        end

        show_window && plot!(plt, ylims(plt)[2] .* da.σ.Confirmed[datarange_dt] .+ ylims(plt)[1], fillrange=ylims(plt)[1],
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

        plot!(plt, da.p.β[datarange_dt], color=:orange, lw=3, label="\\beta (Contact rate)", legend=:topright);

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action] , lw=2, ls=:dash, color = act_colors[i])
            end
        end

        show_window && plot!(plt, ylims(plt)[2] .* da.σ.Confirmed, fillrange=ylims(plt)[1],
                                α=0.2, color=:grey, label="control");

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

        plot!(plt, da.p.γ[datarange_dt], color=:blue, lw=3, label="\\gamma", legend=:topleft);
        plot!(plt, da.p.δ[datarange_dt], color=:black, lw=3, label="\\delta");

        show_window && plot!(plt, ylims(plt)[2] .* da.σ.Confirmed[datarange_dt], fillrange=ylims(plt)[1],
                                α=0.2, color=:grey, label="control");

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action], lw=2, ls=:dash, color = act_colors[i])
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

        plot!(plt, da.p.β[datarange_dt]./(da.p.γ[datarange_dt] .+ da.p.δ[datarange_dt]), color=:red, lw=3, label="\\sigma (Contact number)" , legend=:topright)
        ylims!(0,3)

        plot!(plt, timestamp(da.p.β[datarange_dt]), ones(size(da.p.β[datarange_dt],1)) , fillrange=ylims(plt)[1],
                                α=0.2, color=:green, label=:none)


        show_window && plot!(plt, ylims(plt)[2] .* da.σ.Confirmed, fillrange=ylims(plt)[1],
                                α=0.2, color=:grey, label="control");

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action], lw=2, ls=:dash, color = act_colors[i])
            end
        end

        title!(plt, title_str);

        display_fig && display(plt)
        if save_fig
            for ext in fig_ext
                savefig(plt, joinpath(figdir,"da_contact_number."*ext))
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
        plot!(plt, real.(Λ.λ₁)[datarange_dt], label="\\lambda 1", lw=3, color=:orange)
        plot!(plt, real.(Λ.λ₂)[datarange_dt], label="\\lambda 2", lw=3, color=:blue)
        plot!(plt, real.(Λ.λ₃)[datarange_dt], label="\\lambda 3", lw=3, color=:green)
        plot!(plt, real.(Λ.λ₄)[datarange_dt], label="\\lambda 4", lw=3, color=:red)

        ylims!(plt, -.1, 0.4)

        show_window && plot!(plt, ylims(plt)[2] .* da.σ.Confirmed, fillrange=ylims(plt)[1],
                                α=0.2, color=:grey, label="control");

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action] , lw=2, ls=:dash, color = act_colors[i])
            end
        end

        title!(plt, title_str);

        display_fig && display(plt)
        if save_fig
            for ext in fig_ext
                savefig(plt, joinpath(figdir,"da_lambda."*ext))
            end
        end

        plt
end;

## Plot A22
let plt = Plots.Plot()
        plot!(plt, A22[datarange_dt], label="A22", lw=3, color=:blue)
        # ylims!(plt, -.01, 0.4)

        plot!(plt, timestamp(A22[datarange_dt]), ylims(plt)[1]*ones(size(A22[datarange_dt],1)) , fillrange=0.0,
                                α=0.2, color=:green, label=:none)

        show_window && plot!(plt, ylims(plt)[2] .* da.σ.Confirmed, fillrange=ylims(plt)[1],
                                α=0.2, color=:grey, label="control");

        if show_actions
            for i=1:size(act,1)
                    vline!(plt, [DateTime(act[i,:Date])], label=act[i,:Action] ,  lw=2, ls=:dash, color = act_colors[i])
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
        sticks!(plt, freq.^(-1), P_A22, lw=3, label="FFT(A22)",
                xaxis="T (days)", legend=:top, color=:blue)
        # ylims!(plt, -.01, 0.4)

        xlims!(plt,2,9)
        ylims!(plt,0,10)

        title!(plt, title_str);

        display_fig && display(plt)
        if save_fig
            for ext in fig_ext
                savefig(plt, joinpath(figdir,"da_FFT_A22."*ext))
            end
        end

        plt
end;

## Export simulation result to CSV
let tah = when(da.u, hour, 0)
    tad = TimeArray(Date.(timestamp(tah)), values(tah), colnames(tah))
    writetimearray(tad, joinpath(resdir,"da_SIRD.csv"))
end;

## Export parameters to CSV
let tah = when(da.p, hour, 0)
    tad = TimeArray(Date.(timestamp(tah)), values(tah), colnames(tah))
    writetimearray(tad, joinpath(resdir,"da_param.csv"))
end;
