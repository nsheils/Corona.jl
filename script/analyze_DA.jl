using Corona
using Plots
using FileIO, JLD2
using TimeSeries
using LaTeXStrings
region="Italy"

##
it = load("/home/jls/data/2020-Corona/results/Italy/da.jld2")
#
it = load("results/Italy/da.jld2")

##
itC=rename(it["result"].I .+ it["result"].R .+ it["result"].D,:I_R_D => :C)

globtime = timestamp(it["result"])
datatime = globtime[1:meta(it["data"])["last_day_idxs"]]

##
scatter(it["data"][datatime].Confirmed,color=:orange,Atickfontsize=12,title=region,legend=:topleft)
scatter!(it["data"][datatime].Deaths,color=:black)

plot!(itC[datatime],legend=:topleft,color=:orange,lw=3)
plot!(it["result"][datatime].I,color=:red, lw=3)
plot!(it["result"][datatime].R,color=:green, lw=3)
plot!(it["result"][datatime].D,color=:black, lw=3)

##
savefig("/home/jls/data/2020-Corona/figs/Italy/da.pdf")


##
scatter(it["data"][datatime].Confirmed,color=:orange,
 tickfontsize=12,title=region,legend=:topleft,yaxis=:log10)
scatter!(it["data"][datatime].Deaths .+ 1,color=:black,yaxis=:log10)

plot!(itC[datatime],legend=:topleft,color=:orange,lw=3)
plot!(it["result"][datatime].I,color=:red, lw=3)
plot!(it["result"][datatime].R .+ 1 ,color=:green, lw=3)
plot!(it["result"][datatime].D .+ 1,color=:black, lw=3)

savefig("/home/jls/data/2020-Corona/figs/Italy/da_log.pdf")



plot(it["model_params"][datatime].β,lw=3,label=L"\beta",color=:red)
savefig("/home/jls/data/2020-Corona/figs/Italy/beta.pdf")






##
scatter(it["data"][datatime].Confirmed,color=:orange,tickfontsize=12,title=region,legend=:topleft)
scatter!(it["data"][datatime].Deaths,color=:black)

plot!(itC,legend=:topleft,color=:orange,lw=3)
plot!(it["result"].I,color=:red, lw=3)
plot!(it["result"].R,color=:green, lw=3)
plot!(it["result"].D,color=:black, lw=3)

##
savefig("/home/jls/data/2020-Corona/figs/Italy/da_pre.pdf")

""
