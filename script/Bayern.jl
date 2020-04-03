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
rootdir="/home/tms-archiv/Daten/2020-Corona/"
cd(rootdir)
Today=Dates.today()
## Population size and Country names
Cases=readtimearray("raw/RKI/bayern.csv",header=false)
Confirmed=Cases[:B]
FirstConfirmed=values(Confirmed[1])[1]
LastConfirmed=values(Confirmed[end])[1]
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
nDays=400
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

Δt=DataTimeSpan[2]-DataTimeSpan[1]
β₀=((log(LastConfirmed)-log(FirstConfirmed))/Δt)*ones(nSimul+1)
γ₀=1/7*ones(nSimul+1)
δ₀=0.001*1/7*ones(nSimul+1);

β₀=1/4*ones(nSimul+1)
γ₀=1/7*ones(nSimul+1)
δ₀=zeros(nSimul+1);

β=β₀;γ=γ₀;δ=δ₀
#@load "data/Bavaria/model_parameters.jld"

sponge=1
AssimTime=OutbreakDate:Day(1):LastDate+Day(sponge)
nAssim=length(AssimTime)
AssimTimeSpan=(0.0,Float64(nAssim-1))
AssimTimeRange=1:nAssim
nDataWindow=length(OutbreakDate:Day(1):LastDate)
W=[ones(nDataWindow); zeros(sponge)]
#DataWindow=TimeArray(collect(AssimTime),W,["Data Window"])
u=corona.forward(β,γ,δ,u₀,AssimTimeSpan)
Cₓ=values(uₓ[AssimTime][:Confirmed])
Dₓ=values(uₓ[AssimTime][:Deaths])
J₀=norm([Cₓ.*W;Dₓ.*W])
J=[norm([Cₓ.*W;Dₓ.*W]-[sum(u[AssimTimeRange,:],dims=2).*W;u[AssimTimeRange,4].*W])/J₀]
α=1.0e-7/J₀
Iterations=1000
println("Start with α= $α  and $Iterations Iterations ")
for i=1:Iterations
    global u
    a=corona.backward(u,values(uₓ[AssimTime]),β,δ,γ,reverse(AssimTimeSpan),W);
    global β,γ,δ
    β,γ,δ=corona.linesearch(β,γ,δ,a,values(uₓ[AssimTime]),u₀,AssimTimeSpan,α,2^10,W)
    β[nAssim+1:end].=β[nAssim]
    γ[nAssim+1:end].=γ[nAssim]
    δ[nAssim+1:end].=δ[nAssim]
    β[AssimTimeRange]=∂⁰(β[AssimTimeRange])
    γ[AssimTimeRange]=∂⁰(γ[AssimTimeRange])
    δ[AssimTimeRange]=∂⁰(δ[AssimTimeRange])
    @save "data/Bavaria/model_parameters.jld"  β γ δ
    u=corona.forward(β,γ,δ,u₀,AssimTimeSpan)

    if mod(i,100) == 0
        global J=[J; norm([Cₓ.*W;Dₓ.*W]-[sum(u[:,2:4],dims=2).*W;u[:,4].*W])/J₀]
        println(i," ",J[end])
    end
    if mod(i,1000) == 0
        @save "data/Bavaria/model_parameters.jld"  β γ δ J
        I=u[:,2]
        R=u[:,3]
        D=u[:,4]
        global P=scatter(uₓ[AssimTime][:Confirmed],legend=:topleft,color=:orange,title="Bavaria")
        P=scatter!(uₓ[AssimTime][:Deaths],label="Deaths",color=:black,lw=3,tickfontsize=12)
        P=plot!(AssimTime,sum(u[:,2:4],dims=2),label="C",color=:orange,lw=3)
        P=plot!(AssimTime,u[:,2],label="I",color=:red,lw=3)
        P=plot!(AssimTime,u[:,3],label="R",color=:green,lw=3)
        P=plot!(AssimTime,u[:,4],label="D",color=:black,lw=3)
        display(P)
        savefig("figs/Bavaria/Development.pdf")
        global lP=scatter(uₓ[DataTime][:Confirmed].+1,legend=:topleft,
            color=:orange,title="Bavaria",yaxis=:log10)
        lP=scatter!(uₓ[DataTime][:Deaths] .+1,label="Deaths",color=:black,lw=3,tickfontsize=12)
        lP=plot!(DataTime,sum(u[DataTimeRange,2:4],dims=2),label="C",color=:orange,lw=3)
        lP=plot!(DataTime,u[DataTimeRange,2] .+1,label="I",color=:red,lw=3)
        lP=plot!(DataTime,u[DataTimeRange,3] .+1.0,label="R",color=:green,lw=3)
#        savefig(lP,"figs/Bavaria/CIR_log.pdf")
        lP=plot!(DataTime,u[DataTimeRange,4] .+1.0,label="D",color=:black,lw=3,tickfontsize=12)
        savefig(lP,"figs/Bavaria/Developmentlog.pdf")
        βₘ=maximum(β)
#        pP=plot(DataWindow.*βₘ,lw=0,color=:whitesmoke,fill=(0,:whitesmoke),α=0.9,legend=:topleft)
        pP=plot(DataTime,β[DataTimeRange],label=L"\beta",legend=:left,
                lw=3,title="Bavaria",color=:red,tickfontsize=12)
        savefig(pP,"figs/Bavaria/beta.pdf")
        pP=plot!(DataTime,γ[DataTimeRange],label=L"\gamma",lw=3,color=:green)
        pP=plot!(DataTime,δ[DataTimeRange],label=L"\delta",lw=3,color=:black)
        display(pP)
        savefig(pP,"figs/Bavaria/parameters3.pdf")
        debuggP=plot(β[AssimTimeRange[1]:AssimTimeRange[1]+12],label=L"\beta",legend=:left,
                     lw=3,title="Bavaria",color=:red,tickfontsize=12)
        savefig(debuggP,"figs/Bavaria/dbug.pdf")

   end
end
β[nAssim+1:end].=β[nAssim]
γ[nAssim+1:end].=γ[nAssim]
δ[nAssim+1:end].=δ[nAssim]

β[nAssim+1:end].=corona.extrapolate(β[nAssim-3:nAssim])
γ[nAssim+1:end].=corona.extrapolate(γ[nAssim-3:nAssim])
δ[nAssim+1:end].=corona.extrapolate(δ[nAssim-3:nAssim])
@save "data/Bavaria/model_parameters.jld"  β γ δ
Parameter=TimeArray(SimulTime,[β[SimulTimeRange] γ[SimulTimeRange] δ[SimulTimeRange]],["β", "γ", "δ"])
writetimearray(Parameter, "data/Bavaria/Parameter.csv")#Save data
################### Prognosis ################
u=corona.forward(β,γ,δ,u₀,SimulTimeSpan)
Prognosis=TimeArray(SimulTime,u,["S","I","R","D"])
@save "data/Bavaria/final.jld" Prognosis uₓ β γ δ u₀ SimulTime SimulTimeSpan SimulTimeRange
writetimearray(Prognosis, "data/Bavaria/Prognosis.csv")
