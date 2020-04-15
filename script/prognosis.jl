using corona
using Plots 
using JLD2
using DifferentialEquations
using CSV
using Dates
using TimeSeries
using LinearAlgebra
using DataFrames
using june
using LaTeXStrings
############################################################
rootdir="/home/jls/prog/corona/"
cd(rootdir)
Population_0=CSV.read(rootdir*"data/population.csv",header=false);
Countries=String.(Population_0[:Column1])
Confirmed=Float64.(corona.merge_datasets(corona.read("Confirmed"),Countries))
Deaths=   Float64.(corona.merge_datasets(corona.read("Deaths"),Countries));
Trend=corona.trend(Confirmed,14,40)
Growths=Trend[end]
###########################################
S₀=Float64(Population_0[Population_0[!,:Column1].=="Germany",:][!,:Column2][1])

Today=Dates.today()
@load "data/outbreak.jld2" Outbreak
FirstDate=Outbreak[:Germany]
LastDate=timestamp(Deaths)[end]
DataTime=FirstDate:Day(1):LastDate

tspan=(0.0,Float64(length(DataTime))-1)
@load "data/Germany_final.jld" 


u₀=[values(Germany.S[end])[1],
    values(Germany.I[end])[1],
    values(Germany.R[end])[1],
    values(Germany.D[end])[1]]
NDays=180
β=values(Germany.β[end])[1]*ones(NDays+10)
γ=values(Germany.γ[end])[1]*ones(NDays+10)
δ=values(Germany.δ[end])[1]*ones(NDays+10)
uₓ=values(merge(Confirmed[:Germany],Deaths[:Germany])[DataTime]);
##

SimulRange=DataTime[end]:Day(1):DataTime[end]+Day(NDays)
tsimul=(0.0,Float64(NDays))
u=corona.forward(β,γ,δ,u₀,tsimul)
begin
    scatter(DataTime, uₓ,label=["Confirmed" "Deaths"],title="Prognosis Germany",legend=:left,color=[:orange :black])
    plot!(SimulRange, u,label=["S" "I" "R" "D"],lw=3,color=[:blue :red :green :black])
#    lens!([DataTime[1], DataTime[end]+Day(30)], [0.0, 50000.0], inset=(1, bbox(0.1, 0.5, 0.25, 0.25)))
end
    savefig("figs/Germany_Prognosis.pdf")


maxInfect=maximum(u[:,2])
casulties=u[end,4]

open("data/Germany_Prognosis.csv","w") do prognos
    write(prognos, "Germany, $maxInfect, $casulties")
end

    println("Prognosis for Germany on $LastDate: $maxInfect, $casulties")







