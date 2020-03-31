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
Cases=readtimearray("data/bayern.csv",header=false)
Confirmed=Cases[:B]
Deaths=Cases[:C]
DataTime=timestamp(Cases)
###########################################
S₀=13076721.0
u₀=[S₀,values(Confirmed[1])[1],0.0,values(Deaths[1])[1]]

################### Data and Data Intervall ####################################
Data=merge(
    TimeSeries.rename(Confirmed,:Confirmed),
    TimeSeries.rename(Deaths,:Deaths ))

#
FirstDate=timestamp(Data)[1]
LastDate=timestamp(Data)[end]
DataTime=FirstDate:Day(1):LastDate
nData=length(DataTime)
DataTimeSpan=(0.0,Float64(nData)-1)
DataTimeRange=1:nData

################################################################################
@load "data/outbreak.jld2" Outbreak
OutbreakDate=Outbreak[Symbol("Bavaria")]
#Outbreak[Symbol("Bavaria")]=Date(2020,02,24)
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

########################################


β₀=1/4*ones(nSimul+1)
γ₀=1/7*ones(nSimul+1)
δ₀=0.001*1/7*ones(nSimul+1);
β=β₀;γ=γ₀;δ=δ₀
@load "data/Bavaria_model_parameters.jld"

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
for i=1:30000
    global u
    a=corona.backward(u,values(uₓ[AssimTime]),β,δ,γ,reverse(AssimTimeSpan));
    global β,γ,δ
    β,γ,δ=corona.linesearch(β,γ,δ,a,values(uₓ[AssimTime]),u₀,AssimTimeSpan,α,100)
#    β[nAssim+1:end].=β[nAssim]
#    γ[nAssim+1:end].=γ[nAssim]
#    δ[nAssim+1:end].=δ[nAssim]
#    β[AssimTimeRange]=∂⁰(β)[AssimTimeRange]#
#    γ[AssimTimeRange]=∂⁰(γ)[AssimTimeRange]
#    δ[AssimTimeRange]=∂⁰(δ)[AssimTimeRange]
    β[AssimTimeRange]=smooth(β)[AssimTimeRange]
    γ[AssimTimeRange]=smooth(γ)[AssimTimeRange]
    δ[AssimTimeRange]=smooth(δ)[AssimTimeRange]

    u=corona.forward(β,γ,δ,u₀,AssimTimeSpan)

    if mod(i,1000) == 0
#        @save "data/Bavaria_model_parameters.jld"  β γ δ J
        I=u[:,2]
        R=u[:,3]
        D=u[:,4]
        global J=[J; norm([Cₓ;Dₓ]-[u[:,2]+u[:,3]+u[:,4];u[:,4]])/norm([Cₓ;Dₓ])]
        println(i," ",J[end])
        global P=scatter(uₓ[AssimTime][:Confirmed],legend=:topleft,color=:orange,title="Bavaria")
        P=scatter!(uₓ[AssimTime][:Deaths],label="Deaths",color=:black,lw=3,tickfontsize=12)
        P=plot!(AssimTime,sum(u[:,2:4],dims=2),label="C",color=:orange,lw=3)
        P=plot!(AssimTime,u[:,2],label="I",color=:red,lw=3)
        P=plot!(AssimTime,u[:,3],label="R",color=:green,lw=3)
        P=plot!(AssimTime,u[:,4],label="D",color=:black,lw=3)
        display(P)
        savefig("figs/Bavaria.pdf")
        global lP=scatter(uₓ[AssimTime][:Confirmed].+1,legend=:topleft,
            color=:orange,title="Bavaria",yaxis=:log10)
        lP=scatter!(uₓ[AssimTime][:Deaths] .+1,label="Deaths",color=:black,lw=3,tickfontsize=12)
        lP=plot!(AssimTime,sum(u[AssimTimeRange,2:4],dims=2),label="C",color=:orange,lw=3)
        lP=plot!(AssimTime,u[AssimTimeRange,2] .+1,label="I",color=:red,lw=3)
        lP=plot!(AssimTime,u[AssimTimeRange,3] .+1.0,label="R",color=:green,lw=3)
#        savefig(lP,"figs/Bavaria_CIR_log.pdf")
        lP=plot!(AssimTime,u[AssimTimeRange,4] .+1.0,label="D",color=:black,lw=3,tickfontsize=12)
        savefig(lP,"figs/Bavaria_log.pdf")
        pP=plot(AssimTime,β[AssimTimeRange],label=L"\beta",legend=:left,
                lw=3,title="Bavaria",color=:red,tickfontsize=12)
        savefig(pP,"figs/Bavaria_beta.pdf")
        pP=plot!(AssimTime,γ[AssimTimeRange],label=L"\gamma",lw=3,color=:green)
        pP=plot!(AssimTime,δ[AssimTimeRange],label=L"\delta",lw=3,color=:black)
        display(pP)
        savefig(pP,"figs/Bavaria_parameters3.pdf")
   end
end
β[nAssim+1:end].=extrapolate(β[nAssim-3:nAssim])
γ[nAssim+1:end].=extrapolate(γ[nAssim-3:nAssim])
δ[nAssim+1:end].=extrapolate(δ[nAssim-3:nAssim])

#Save data
################### Prognosis ################
u=corona.forward(β,γ,δ,u₀,SimulTimeSpan)$
Prognosis=TimeArray(SimulTime,u,["S","I","R","D"])
@save "data/Bavaria_final.jld" Prognosis uₓ β γ δ u₀ SimulTime SimulTimeSpan SimulTimeRange
PrognosisPlot=scatter(uₓ[DataTime],label=["Confirmed" "Deaths"],
            title="Prognosis Bavaria",legend=:right,color=[:orange :black])
PrognosisPlot=plot!(SimulTime, u,label=["S" "I" "R" "D"],lw=3,color=[:blue :red :green :black])
savefig(PrognosisPlot,"figs/Bavaria_Prognosis.pdf")
Parameter=TimeArray(SimulTime,[β[SimulTimeRange] γ[SimulTimeRange] δ[SimulTimeRange]],
    ["β", "γ", "δ"])
writetimearray(Prognosis, "data/Bavaria_Prognosis.csv")
writetimearray(Parameter, "data/Bavaria_Parameter.csv")

WeekTime=SimulTime[1]:Day(1):DataTime[end]+Day(8)
nWeek=length(WeekTime)
WeekTimeSpan=(0.0,Float64(nWeek)-1)
WeekTimeRange=1:nWeek
WeekPrognosisPlot=scatter(uₓ[DataTime],label=["Confirmed","Deaths"],
            title="One Week Prognosis Bavaria",legend=:left,color=[:orange :black])
WeekPrognosisPlot=plot!(WeekTime, u[WeekTimeRange,2:4] ,label=["I" "R" "D"],lw=3,color=[:red :green :black])
WeekPrognosisPlot=plot!(WeekTime, sum(u[WeekTimeRange,2:4],dims=2) ,label="C",lw=3,color=[:orange])
savefig(WeekPrognosisPlot,"figs/Bavaria_WeekPrognosis.pdf")

MonthTime=SimulTime[1]:Day(1):DataTime[end]+Day(30)
nMonth=length(MonthTime)
MonthTimeSpan=(0.0,Float64(nMonth)-1)
MonthTimeRange=1:nMonth
MonthPrognosisPlot=scatter(uₓ[DataTime],label=["Confirmed","Deaths"],
            title="One Month Prognosis Bavaria",legend=:left,color=[:orange :black])
MonthPrognosisPlot=plot!(MonthTime, u[MonthTimeRange,2:4] ,label=["I" "R" "D"],lw=3,color=[:red :green :black])
MonthPrognosisPlot=plot!(MonthTime, sum(u[MonthTimeRange,2:4],dims=2) ,label="C",lw=3,color=[:orange])
savefig(MonthPrognosisPlot,"figs/Bavaria_MonthPrognosis.pdf")


writetimearray(Prognosis[MonthTime], "data/Bavaria_MonthPrognosis.csv")
writetimearray(Parameter[MonthTime], "data/Bavaria_MonthParameter.csv")
writetimearray(Prognosis[WeekTime], "data/Bavaria_WeekPrognosis.csv")
writetimearray(Parameter[WeekTime], "data/Bavaria_WeekParameter.csv")
