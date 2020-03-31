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
## Population size and Country names
Cases=readtimearray("data/bayern.csv",header=false)
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
@load "data/Bavaria_model_parameters.jld"

sponge=180
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
α=1.0e-8/J₀
println("α: $α")
for i=1:100000
    global u
    a=corona.backward(u,values(uₓ[AssimTime]),β,δ,γ,reverse(AssimTimeSpan),W);
    global β,γ,δ
    β,γ,δ=corona.linesearch(β,γ,δ,a,values(uₓ[AssimTime]),u₀,AssimTimeSpan,α,100,W)
    β[nAssim+1:end].=β[nAssim]
    γ[nAssim+1:end].=γ[nAssim]
    δ[nAssim+1:end].=δ[nAssim]
    β[AssimTimeRange]=∂⁰(β[AssimTimeRange])
    γ[AssimTimeRange]=∂⁰(γ[AssimTimeRange])
    δ[AssimTimeRange]=∂⁰(δ[AssimTimeRange])
    @save "data/Bavaria_model_parameters.jld"  β γ δ
    u=corona.forward(β,γ,δ,u₀,AssimTimeSpan)

    if mod(i,100) == 0
        global J=[J; norm([Cₓ.*W;Dₓ.*W]-[sum(u[:,2:4],dims=2).*W;u[:,4].*W])/J₀]
        println(i," ",J[end])
    end
    if mod(i,100) == 0
        @save "data/Bavaria_model_parameters.jld"  β γ δ J
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
        βₘ=maximum(β)
#        pP=plot(DataWindow.*βₘ,lw=0,color=:whitesmoke,fill=(0,:whitesmoke),α=0.9,legend=:topleft)
        pP=plot(AssimTime,β[AssimTimeRange],label=L"\beta",legend=:left,
                lw=3,title="Bavaria",color=:red,tickfontsize=12)
        savefig(pP,"figs/Bavaria_beta.pdf")
        pP=plot!(AssimTime,γ[AssimTimeRange],label=L"\gamma",lw=3,color=:green)
        pP=plot!(AssimTime,δ[AssimTimeRange],label=L"\delta",lw=3,color=:black)
        display(pP)
        savefig(pP,"figs/Bavaria_parameters3.pdf")
        debuggP=plot(β[AssimTimeRange[1]:AssimTimeRange[1]+12],label=L"\beta",legend=:left,
                     lw=3,title="Bavaria",color=:red,tickfontsize=12)
        savefig(debuggP,"figs/Bavaria_dbug.pdf")

   end
end
β[nAssim+1:end].=β[nAssim]
γ[nAssim+1:end].=γ[nAssim]
δ[nAssim+1:end].=δ[nAssim]

β[nAssim+1:end].=corona.extrapolate(β[nAssim-3:nAssim])
γ[nAssim+1:end].=corona.extrapolate(γ[nAssim-3:nAssim])
δ[nAssim+1:end].=corona.extrapolate(δ[nAssim-3:nAssim])
@save "data/Bavaria_model_parameters.jld"  β γ δ
Parameter=TimeArray(SimulTime,[β[SimulTimeRange] γ[SimulTimeRange] δ[SimulTimeRange]],["β", "γ", "δ"])
writetimearray(Parameter, "data/Bavaria_Parameter.csv")#Save data
################### Prognosis ################
u=corona.forward(β,γ,δ,u₀,SimulTimeSpan)
Prognosis=TimeArray(SimulTime,u,["S","I","R","D"])
@save "data/Bavaria_final.jld" Prognosis uₓ β γ δ u₀ SimulTime SimulTimeSpan SimulTimeRange
writetimearray(Prognosis, "data/Bavaria_Prognosis.csv")
