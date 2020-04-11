using Corona
using Plots
using FileIO, JLD2
using TimeSeries
using LaTeXStrings

##
region="New York City"
regiondir="US-New_York-New_York"
figpath=joinpath("figs",regiondir)
filename=joinpath("results",regiondir,"da.jld2")

##
da = load(filename)

##
daC=rename(da["result"].I .+ da["result"].R .+ da["result"].D,:I_R_D => :C)

globtime = timestamp(da["result"])
datatime = globtime[1:meta(da["data"])["last_day_idxs"]]

##
scatter(da["data"][datatime].Confirmed,color=:orange,Atickfontsize=12,title=region,legend=:topleft)
scatter!(da["data"][datatime].Deaths,color=:black)

plot!(daC[datatime],legend=:topleft,color=:orange,lw=3)
plot!(da["result"][datatime].I,color=:red, lw=3)
plot!(da["result"][datatime].R,color=:green, lw=3)
plot!(da["result"][datatime].D,color=:black, lw=3)

##
savefig(joinpath(figpath,"da.pdf"))


##
scatter(da["data"][datatime].Confirmed,color=:orange,
 tickfontsize=12,title=region,legend=:topleft,yaxis=:log10)
scatter!(da["data"][datatime].Deaths .+ 1,color=:black,yaxis=:log10)
plot!(daC[datatime],legend=:topleft,color=:orange,lw=3)
plot!(da["result"][datatime].I,color=:red, lw=3)
plot!(da["result"][datatime].R .+ 1 ,color=:green, lw=3)
plot!(da["result"][datatime].D .+ 1,color=:black, lw=3)

##
# savefig("/home/jls/data/2020-Corona/figs/Italy/da_log.pdf")
savefig(joinpath(figpath,"da_log.pdf"))

##
scatter(da["model_params"][datatime].β,lw=3,label=L"\beta",title=region,color=:blue)
scatter!(da["model_params"][datatime].γ,lw=3,label=L"\gamma",title=region,color=:red)
scatter!(da["model_params"][datatime].δ,lw=3,label=L"\delta",title=region,color=:orange)


##
t=Date(2020,03,21):Day(1):Date(2020,04,08)
σ = da["model_params"][t].β./(da["model_params"][t].γ.+da["model_params"][t].δ)
scatter(σ,lw=3,label=L"\sigma",title=region,color=:blue)

##
savefig(joinpath(figpath,"da_sigma.pdf"))


##
scatter(da["data"][datatime].Confirmed,color=:orange,tickfontsize=12,title=region,legend=:topleft)
scatter!(da["data"][datatime].Deaths,color=:black)

plot!(daC,legend=:topleft,color=:orange,lw=3)
plot!(da["result"].I,color=:red, lw=3)
plot!(da["result"].R,color=:green, lw=3)
plot!(da["result"].D,color=:black, lw=3)

##
savefig(joinpath(figpath,"da_fc.pdf"))

""
