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
println("C O R O N A   ",String(:Germany))
#rootdir="/home/tms-archiv/Daten/2020-Corona/"
rootdir="/home/jls/data/2020-Corona/"
startdir=pwd()
cd(rootdir)
Today=Dates.today()

println(Today)
## Population size and Country names

println("Data read")
D=corona.data(:Germany)
S₀=D.population
Data=D.cases
OutbreakDate=D.outbreak
DataTime=D.TimeRange
DataIndexRange=D.TimeIndexRange
##
FirstDate=timestamp(D)[1]
LastDate=timestamp(D)[end]

nData=length(D)
DataTimeSpan=D.TimeIntervall
DataTimeRange=D.TimeRange
println("Outbreak : $OutbreakDate")
println("Last Date: $LastDate ",  Int(values(D.cases[:Confirmed])[end]),
    "  ",Int(values(D.cases[:Deaths])[end]))


##
nDays=400
SimulLast=FirstDate+Day(nDays)
GlobalTime=FirstDate:Day(1):SimulLast
nGlobal=length(GlobalTime)
SimulTime=OutbreakDate:Day(1):SimulLast
nSimul=length(SimulTime)
SimulTimeSpan=(0.0,Float64(nSimul-1))
SimulTimeRange=1:nSimul
##########################################
A=zeros(nGlobal,2);A[DataIndexRange,:]=values(D)
uₓ=TimeArray(GlobalTime,A,TimeSeries.colnames(Data))
########################################



β₀=1/4*ones(nSimul+1)
γ₀=1/10*ones(nSimul+1)
δ₀=zeros(nSimul+1);

β=β₀;γ=γ₀;δ=δ₀
#@load "data/Germany_model_parameters.jld"

sponge=1
AssimTime=OutbreakDate:Day(1):LastDate+Day(sponge)
nAssim=length(AssimTime)
AssimTimeSpan=(0.0,Float64(nAssim-1))
AssimTimeRange=1:nAssim
nDataWindow=length(OutbreakDate:Day(1):LastDate)
W=[ones(nDataWindow); zeros(sponge)]
#DataWindow=TimeArray(collect(AssimTime),W,["Data Window"])
##
u₀=[S₀,values(D.cases[:Confirmed][OutbreakDate])[1],0.0,values(D.cases[:Deaths][OutbreakDate])[1]]
u=corona.forward(β,γ,δ,u₀,AssimTimeSpan)
Cₓ=values(uₓ[AssimTime][:Confirmed])
Dₓ=values(uₓ[AssimTime][:Deaths])
J₀=norm([Cₓ.*W;Dₓ.*W])
J=[norm([Cₓ.*W;Dₓ.*W]-[sum(u[AssimTimeRange,:],dims=2).*W;u[AssimTimeRange,4].*W])/J₀]
α=1.0e-1/J₀
println("α: $α")
a=corona.backward(u,values(uₓ[AssimTime]),β,δ,γ,reverse(AssimTimeSpan),W);
##
#β,γ,δ=corona.linesearch(β,γ,δ,a,values(uₓ[AssimTime]),u₀,AssimTimeSpan,α,1024,W)
##
for i=1:10000
    global u
    a=corona.backward(u,values(uₓ[AssimTime]),β,δ,γ,reverse(AssimTimeSpan),W);
    global β,γ,δ
    β,γ,δ=corona.linesearch(β,γ,δ,a,values(uₓ[AssimTime]),u₀,AssimTimeSpan,α,2^13,W)
    β[nAssim+1:end].=β[nAssim]
    γ[nAssim+1:end].=γ[nAssim]
    δ[nAssim+1:end].=δ[nAssim]
    β[DataIndexRange]=∂⁰(β)[DataIndexRange]
    γ[DataIndexRange]=∂⁰(γ)[DataIndexRange]
    δ[DataIndexRange]=∂⁰(δ)[DataIndexRange]
#    @save "data/Germany_model_parameters.jld"  β γ δ
    u=corona.forward(β,γ,δ,u₀,AssimTimeSpan)

    if mod(i,100) == 0
        global J=[J; norm([Cₓ.*W;Dₓ.*W]-[sum(u[:,2:4],dims=2).*W;u[:,4].*W])/J₀]
        println(i," ",J[end])
    end
    if mod(i,100) == 0
#        @save "data/Germany_model_parameters.jld"  β γ δ J
        I=u[:,2]
        R=u[:,3]
        D=u[:,4]
        global P=scatter(uₓ[OutbreakDate:Day(1):LastDate][:Confirmed],legend=:topleft,color=:orange,title="Germany")
        P=scatter!(uₓ[OutbreakDate:Day(1):LastDate][:Deaths],label="Deaths",color=:black,lw=3,tickfontsize=12)
        P=plot!(AssimTime,sum(u[:,2:4],dims=2),label="C",color=:orange,lw=3)
        P=plot!(AssimTime,u[:,2],label="I",color=:red,lw=3)
        P=plot!(AssimTime,u[:,3],label="R",color=:green,lw=3)

        P=plot!(AssimTime,u[:,4],label="D",color=:black,lw=3)
        display(P)
        savefig("figs/Germany.pdf")
        global lP=scatter(uₓ[OutbreakDate:Day(1):LastDate][:Confirmed],legend=:topleft,
            color=:orange,title="Germany",yaxis=:log10)
        lP=scatter!(uₓ[OutbreakDate:Day(1):LastDate][:Deaths] .+1,label="Deaths",color=:black,lw=3,tickfontsize=12)
        lP=plot!(AssimTime,sum(u[AssimTimeRange,2:4],dims=2),label="C",color=:orange,lw=3)
        lP=plot!(AssimTime,u[AssimTimeRange,2],label="I",color=:red,lw=3)
        lP=plot!(AssimTime,u[AssimTimeRange,3] .+1.0,label="R",color=:green,lw=3)
#        savefig(lP,"figs/Germany_CIR_log.pdf")
        lP=plot!(AssimTime,u[AssimTimeRange,4] .+1.0,label="D",color=:black,lw=3,tickfontsize=12)
        savefig(lP,"figs/Germany_log.pdf")
        βₘ=maximum(β)
#        pP=plot(DataWindow.*βₘ,lw=0,color=:whitesmoke,fill=(0,:whitesmoke),α=0.9,legend=:topleft)
        pP=plot(DataTime,β[DataIndexRange],label=L"\beta",legend=:left,
                lw=3,title="Germany",color=:red,tickfontsize=12)
        savefig(pP,"figs/Germany_beta.pdf")
        pP=plot!(DataTime,γ[DataIndexRange],label=L"\gamma",lw=3,color=:green)
        pP=plot!(DataTime,δ[DataIndexRange],label=L"\delta",lw=3,color=:black)
        display(pP)
        savefig(pP,"figs/Germany_parameters3.pdf")
        debuggP=plot(β[AssimTimeRange[1]:AssimTimeRange[1]+12],label=L"\beta",legend=:left,
                     lw=3,title="Germany",color=:red,tickfontsize=12)
        savefig(debuggP,"figs/Germany_dbug.pdf")

   end
end
##



β[nAssim+1:end].=β[nAssim]
γ[nAssim+1:end].=γ[nAssim]
δ[nAssim+1:end].=δ[nAssim]

β[nAssim+1:end].=corona.extrapolate(β[nAssim-3:nAssim])
γ[nAssim+1:end].=corona.extrapolate(γ[nAssim-3:nAssim])
δ[nAssim+1:end].=corona.extrapolate(δ[nAssim-3:nAssim])
@save "data/Germany_model_parameters.jld"  β γ δ
Parameter=TimeArray(SimulTime,[β[SimulTimeRange] γ[SimulTimeRange] δ[SimulTimeRange]],["β", "γ", "δ"])
writetimearray(Parameter, "data/Germany_Parameter.csv")#Save data
################### Prognosis ################
u=corona.forward(β,γ,δ,u₀,SimulTimeSpan)
Prognosis=TimeArray(SimulTime,u,["S","I","R","D"])
@save "data/Germany_final.jld" Prognosis uₓ β γ δ u₀ SimulTime SimulTimeSpan SimulTimeRange
writetimearray(Prognosis, "data/Germany_Prognosis.csv")
cd(startdir)
