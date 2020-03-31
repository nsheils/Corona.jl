using corona
using Plots
using JLD2
using CSV
using Dates
using TimeSeries
using LinearAlgebra
using DataFrames
using LaTeXStrings
############################################################
rootdir="/home/jls/prog/corona/"
cd(rootdir)
Today=Dates.today()
################################################################################
province="Germany"
@load "data/outbreak.jld2" Outbreak
OutbreakDate=Outbreak[Symbol(province)]
@load "data/$province"*"Germany_final.jld" Prognosis Data
S=values(Prognosis.S)
I=values(Prognosis.I)
R=values(Prognosis.R)
D=values(Prognosis.D)
β=values(Prognosis.β)
γ=values(Prognosis.γ)
δ=values(Prognosis.δ)
N=S+I+R+D
C=TimeSeries.rename(Prognosis.I .+Prognosis.R .+Prognosis.D,:I_R_D=>:C)
σ=(β.*S./N).-(γ + δ )
SimulTime=timestamp(Prognosis)
DataTime=timestamp(Data)


WeekTime=SimulTime[1]:Day(1):DataTime[end]+Day(8)
nWeek=length(WeekTime)
WeekTimeSpan=(0.0,Float64(nWeek)-1)
WeekTimeRange=1:nWeek

MonthTime=SimulTime[1]:Day(1):DataTime[end]+Day(30)
nMonth=length(MonthTime)
MonthTimeSpan=(0.0,Float64(nMonth)-1)
MonthTimeRange=1:nMonth


#####################  Infection Rate  ###############################
plot(WeekTime,σ[WeekTimeRange],tickfontsize=12,lw=3,label=L"\lambda",title="Growth Rate of Infections in Germany ",color=:red)
#yaxis!((-0.01,0.3))
savefig("Germany_WeekInfectionrate.pdf")
plot(MonthTime,σ[MonthTimeRange],tickfontsize=12,lw=3,label=L"\lambda",title="Growth Rate of Infections in Germany ",color=:red)
#yaxis!((-0.01,0.3))
savefig("Germany_MonthInfectionrate.pdf")
#####################  Prognosis  ###############################
WeekPrognosisPlot=scatter(Data,label=["Confirmed","Deaths"],
                          title="One Week Prognosis Germany",legend=:left,color=[:orange :black],tickfontsize=12)
WeekPrognosisPlot=plot!(C[WeekTime] ,lw=3,color=:orange)
WeekPrognosisPlot=plot!(Prognosis.I[WeekTime] ,label="I",lw=3,color=:red)
WeekPrognosisPlot=plot!(Prognosis.R[WeekTime] ,label="I",lw=3,color=:green)
WeekPrognosisPlot=plot!(Prognosis.D[WeekTime] ,label="D",lw=3,color=:black)
savefig(WeekPrognosisPlot,"figs/Germany_WeekPrognosis.pdf")


WeekPrognosisPlot=scatter(Data,label=["Confirmed","Deaths"],
                          title="One Week Prognosis Germany",legend=:left,color=[:orange :black],tickfontsize=12)
WeekPrognosisPlot=plot!(C[WeekTime] ,lw=3,color=:orange)
WeekPrognosisPlot=plot!(Prognosis.I[WeekTime] ,label="I",lw=3,color=:red)
WeekPrognosisPlot=plot!(Prognosis.R[WeekTime] ,label="I",lw=3,color=:green)
WeekPrognosisPlot=plot!(Prognosis.D[WeekTime] ,label="D",lw=3,color=:black)
savefig(WeekPrognosisPlot,"figs/Germany_WeekPrognosis.pdf")
WeekPrognosisPlot=scatter(Data,label=["Confirmed","Deaths"],
                          title="One Week Prognosis Germany",legend=:left,color=[:orange :black],tickfontsize=12)
WeekPrognosisPlot=plot!(C[WeekTime] ,lw=3,color=:orange)
WeekPrognosisPlot=plot!(Prognosis.I[WeekTime] ,label="I",lw=3,color=:red)
WeekPrognosisPlot=plot!(Prognosis.R[WeekTime] ,label="I",lw=3,color=:green)
WeekPrognosisPlot=plot!(Prognosis.D[WeekTime] ,label="D",lw=3,color=:black)
savefig(WeekPrognosisPlot,"figs/Germany_WeekPrognosis.pdf")


MonthPrognosisPlot=scatter(Data,label=["Confirmed","Deaths"],
                          title="One Month Prognosis Germany",legend=:left,color=[:orange :black],tickfontsize=12)
MonthPrognosisPlot=plot!(C[MonthTime] ,lw=3,color=:orange)
MonthPrognosisPlot=plot!(Prognosis.I[MonthTime] ,label="I",lw=3,color=:red)
MonthPrognosisPlot=plot!(Prognosis.R[MonthTime] ,label="I",lw=3,color=:green)
MonthPrognosisPlot=plot!(Prognosis.D[MonthTime] ,label="D",lw=3,color=:black)
savefig(MonthPrognosisPlot,"figs/Germany_WeekPrognosis.pdf")

