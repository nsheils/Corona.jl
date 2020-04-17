using Plots

import Plots.plot
function plot(da::DA,region::String;yaxis=:none::Symbol,legend=:topleft::Symbol)
    if yaxis == :log10
        add = 1.0
    else
        add = 0.0
    end
    datatime = timestamp(da.data)[1:meta(da.data)["last_day_idxs"]]
    scatter(da.data.Confirmed[datatime],color=:orange,
                Atickfontsize=12,title=region,legend=legend,yaxis=yaxis)
    scatter!(da.data.Deaths[datatime].+add,color=:black)
    plot!(da.u.I .+ da.u.R .+ da.u.D,color=:orange,lw=3,label="C")
    plot!(da.u.I,color=:red, lw=3)
    plot!(da.u.R .+ add,color=:green, lw=3)
    plot!(da.u.D .+ add,color=:black, lw=3)
end
