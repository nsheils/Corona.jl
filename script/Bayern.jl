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
rootdir="/home/jls/data/2020-Corona/"

Iterations=1000000
FilterFreq=100
ScreenFreq=100
PlotFreq=10000
sponge=2
############################################################
ColdStart=true
β₀=1/4*ones(nSimul+1)
γ₀=1/7*ones(nSimul+1)
δ₀=zeros(nSimul+1);
############################################################
parameter=DataFrame(
parameter=["rootdir", "ColdStart", "Iterations", "FilterFreq", "ScreenFreq","PlotFreq","sponge"],
value=[rootdir, ColdStart, Iterations, FilterFreq, ScreenFreq,PlotFreq,sponge])
println(parameter)



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
OutbreakDate=Date(2020,2,24)
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




if ColdStart == true
    printstyled("Cold Start";color=:red)
    β=β₀;γ=γ₀;δ=δ₀
else
    printstyled("Warm Start";color=:green)
    @load "data/Bavaria/model_parameters.jld"
end
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
α=1.0/J₀
J=[norm([Cₓ.*W;Dₓ.*W]-[sum(u[AssimTimeRange,:],dims=2).*W;u[AssimTimeRange,4].*W])/J₀]
println(" with α= $α  and $Iterations Iterations ")
for i=1:Iterations
    global u
    v=corona.backward(u,values(uₓ[AssimTime]),β,δ,γ,reverse(AssimTimeSpan),W);
    global β,γ,δ
    β,γ,δ=corona.linesearch(β,γ,δ,v,values(uₓ[AssimTime]),
                            u₀,AssimTimeSpan,α,2^12,W)
    u=corona.forward(β,γ,δ,u₀,AssimTimeSpan)
    if mod(i,FilterFreq) == 0
        β=∂⁰(β)
        γ=∂⁰(γ)
        δ=∂⁰(δ)
        β[nAssim-sponge+1:end].=β[nAssim-sponge]
        γ[nAssim-sponge+1:end].=γ[nAssim-sponge]
        δ[nAssim-sponge+1:end].=δ[nAssim-sponge]
        U=TimeArray(collect(AssimTime),u,["S", "I" ,"R" ,"D"])
        V=TimeArray(collect(AssimTime),v,["S", "I" ,"R" ,"D"])
        P=TimeArray(collect(AssimTime),[ β[AssimTimeRange] γ[AssimTimeRange] δ[AssimTimeRange] ],["β", "γ" ,"δ"])
        @save "data/Bavaria/solution.jld" AssimTime U V P
        @save "data/Bavaria/model_parameters.jld"  β γ δ
    end
    if mod(i,ScreenFreq) == 0
        global J=[J; norm([Cₓ.*W;Dₓ.*W]-[sum(u[:,2:4],dims=2).*W;u[:,4].*W])/J₀]
        if J[end]==minimum(J)
            color=:green 
        elseif J[end]>=J[end-1]
            color=:red
        else
            color=:orange
        end
        print(i," ")
        printstyled(J[end],"\n";color=color)
    end
    if mod(i,PlotFreq) == 0
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
        debuggP=plot(β[AssimTimeRange[1]:AssimTimeRange[end]+30],label=L"\beta",legend=:left,
                     lw=3,title="Bavaria",color=:red,tickfontsize=12)
        savefig(debuggP,"figs/Bavaria/dbug.pdf")

   end
end
#prolong data
β[nAssim-sponge+1:end].=β[nAssim-sponge]
γ[nAssim-sponge+1:end].=γ[nAssim-sponge]
δ[nAssim-sponge+1:end].=δ[nAssim-sponge]
