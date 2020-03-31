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
Today=Dates.today()
## Populatuion size and Country names
Population_0=CSV.read(rootdir*"data/population.csv",header=false);
S₀=Float64(Population_0[Population_0[!,:Column1].=="United Kingdom",:][!,:Column2][1])
Countries=String.(Population_0[:Column1])
################### Data and Data Intervall ####################################
Confirmed=Float64.(corona.merge_datasets(corona.read("Confirmed"),Countries))
Deaths=   Float64.(corona.merge_datasets(corona.read("Deaths"),Countries));
Data=merge(
    TimeSeries.rename(Confirmed[Symbol("United Kingdom")],:Confirmed),
    TimeSeries.rename(Deaths[Symbol("United Kingdom")],:Deaths ))
Trend=corona.trend(Confirmed,14,40)
Growths=Trend[end]
#
FirstDate=timestamp(Data)[1]
LastDate=timestamp(Data)[end]
DataTime=FirstDate:Day(1):LastDate
nData=length(DataTime)
DataTimeSpan=(0.0,Float64(nData)-1)
DataTimeRange=1:nData

################################################################################
@load "data/outbreak.jld2" Outbreak
OutbreakDate=Outbreak[Symbol("United Kingdom")]
#@save "data/outbreak.jld2" Outbreak
#
nDays=180
SimulLast=FirstDate+Day(nDays)
GlobalTime=FirstDate:Day(1):SimulLast
nGlobal=length(GlobalTime)
SimulTime=OutbreakDate:Day(1):SimulLast
nSimul=length(SimulTime)
SimulTimeSpan=(0.0,Float64(nSimul-1))
SimulTimeRange=1:nSimul
##########################################
A=zeros(nGlobal,2);A[DataTimeRange,:]=values(Data)
uₓ=TimeArray(GlobalTime,A,TimeSeries.colnames(Data))
u₀=[S₀,Float64(values(Confirmed[Symbol("United Kingdom")][OutbreakDate])[1]),0.0,0.0]
β₀=1/4*ones(nSimul+1)
γ₀=1/7*ones(nSimul+1)
δ₀=0.001*1/7*ones(nSimul+1);
##
u₈=corona.forward(β₀,γ₀,δ₀,u₀,SimulTimeSpan)
Simul=plot(uₓ[DataTime],lw=3)
Simul=plot!(DataTime, u₈[DataTimeRange,2])
plot(Confirmed,legend=:topleft,lw=3)
savefig("figs/Confirmed.pdf")
plot(Trend,lw=3)
savefig("figs/Trend_avg.pdf")
########################################
β,γ,δ=β₀,γ₀,δ₀
@load "data/United Kingdom_model_parameters.jld"

AssimTime=OutbreakDate:Day(1):LastDate
#AssimTime=OutbreakDate:Day(1):Date(2020,03,14)
nAssim=length(AssimTime)
AssimTimeSpan=(0.0,Float64(nAssim-1))
AssimTimeRange=1:nAssim

u=corona.forward(β,γ,δ,u₀,AssimTimeSpan)
Cₓ=values(uₓ[AssimTime][:Confirmed])
Dₓ=values(uₓ[AssimTime][:Deaths])
J=[norm([Cₓ;Dₓ]-[sum(u[AssimTimeRange,:],dims=2);u[AssimTimeRange,4]])/norm([Cₓ;Dₓ])]
J₀=norm([Cₓ;Dₓ])
α=1.0e-7/J₀
println("α: $α")
for i=1:20000
    global u
    a=corona.backward(u,values(uₓ[AssimTime]),β,δ,γ,reverse(AssimTimeSpan));
    global β,γ,δ
    β,γ,δ=corona.linesearch(β,γ,δ,a,values(uₓ[AssimTime]),u₀,AssimTimeSpan,α,100)
    β[AssimTimeRange]=∂⁰(β)[AssimTimeRange]
    γ[AssimTimeRange]=∂⁰(γ)[AssimTimeRange]
    δ[AssimTimeRange]=∂⁰(δ)[AssimTimeRange]
    u=corona.forward(β,γ,δ,u₀,AssimTimeSpan)
    if mod(i,1000) == 0
        @save "data/United Kingdom_model_parameters.jld" FirstDay β γ δ J
        I=u[:,2]
        R=u[:,3]
        D=u[:,4]
        global J=[J; norm([Cₓ;Dₓ]-[u[:,2]+u[:,3]+u[:,4];u[:,4]])/norm([Cₓ;Dₓ])]
        println(i," ",J[end])
        global P=scatter(uₓ[AssimTime][:Confirmed],legend=:topleft,color=:orange,title="United Kingdom")
        P=scatter!(uₓ[AssimTime][:Deaths],label="Deaths",color=:black,lw=3)
        P=plot!(AssimTime,sum(u[:,2:4],dims=2),label="C",color=:orange,lw=3)
        P=plot!(AssimTime,u[:,2],label="I",color=:red,lw=3)
        P=plot!(AssimTime,u[:,3],label="R",color=:green,lw=3)
        P=plot!(AssimTime,u[:,4],label="D",color=:black,lw=3)
        display(P)
        savefig("figs/United Kingdom.pdf")
        global lP=scatter(uₓ[AssimTime][:Confirmed].+1,legend=:topleft,
            color=:orange,title="United Kingdom",yaxis=:log10)
        lP=scatter!(uₓ[AssimTime][:Deaths] .+1,label="Deaths",color=:black,lw=3)
        lP=plot!(AssimTime,sum(u[AssimTimeRange,2:4],dims=2),label="C",color=:orange,lw=3)
        lP=plot!(AssimTime,u[AssimTimeRange,2] .+1,label="I",color=:red,lw=3)
        lP=plot!(AssimTime,u[AssimTimeRange,3] .+1.0,label="R",color=:green,lw=3)
#        savefig(lP,"figs/United Kingdom_CIR_log.pdf")
        lP=plot!(AssimTime,u[AssimTimeRange,4] .+1.0,label="D",color=:black,lw=3)
        savefig(lP,"figs/United Kingdom_log.pdf")
        pP=plot(AssimTime,β[AssimTimeRange],label=L"\beta",legend=:left,
                lw=3,title="United Kingdom",color=:red)
        savefig(pP,"figs/United Kingdom_beta.pdf")
        pP=plot!(AssimTime,γ[AssimTimeRange],label=L"\gamma",lw=3,color=:green)
        pP=plot!(AssimTime,δ[AssimTimeRange],label=L"\delta",lw=3,color=:black)
        display(pP)
        savefig(pP,"figs/United Kingdom_parameters3.pdf")
   end
end
## prolong parameters
β[nAssim+1:end].=β[nAssim]
γ[nAssim+1:end].=γ[nAssim]
δ[nAssim+1:end].=δ[nAssim]
#Save data
################### Prognosis ################
u=corona.forward(β,γ,δ,u₀,SimulTimeSpan)
Prognosis=TimeArray(SimulTime,u,["S","I","R","D"])
@save "data/United Kingdom_final.jld" Prognosis uₓ β γ δ
PrognosisPlot=scatter(uₓ[DataTime],label=["Confirmed" "Deaths"],
            title="Prognosis United Kingdom",legend=:left,color=[:orange :black])
PrognosisPlot=plot!(SimulTime, u,label=["S" "I" "R" "D"],lw=3,color=[:blue :red :green :black])

savefig(PrognosisPlot,"figs/United Kingdom_Prognosis.pdf")
Parameter=TimeArray(SimulTime,[β[SimulTimeRange] γ[SimulTimeRange] δ[SimulTimeRange]],
    ["β", "γ", "δ"])

writetimearray(Prognosis, "data/United Kingdom_Prognosis.csv")
writetimearray(Parameter, "data/United Kingdom_Parameter.csv")
