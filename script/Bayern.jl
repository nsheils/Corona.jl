vusing corona
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
rootdir      = "/home/jls/data/2020-Corona/"
Iterations   = 1000000
FilterFreq   = 10000000
ScreenFreq   = 100
PlotFreq     = 100
sponge       = 1
LineIterMax  = 2^12
LineRangeMax = 1e-6
ColdStart    = false
############################################################
Actions=DataFrame(Date=Date(2020,03,16),Action="Schulschließung")
push!(Actions,[Date(2020,03,20),"Vorläufige Ausgangsbeschränkung"])
parameter=DataFrame(
    parameter=["rootdir", "ColdStart", "Iterations", "FilterFreq", "ScreenFreq","PlotFreq","sponge","LineIterMax","LineRangeMax"],
    value=[rootdir, ColdStart, Iterations, FilterFreq, ScreenFreq,PlotFreq,sponge,LineIterMax,LineRangeMax])
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
    β₀=1/4*ones(nSimul+1)
    γ₀=1/7*ones(nSimul+1)
    δ₀=zeros(nSimul+1);
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
α=LineRangeMax/J₀
J=[norm([Cₓ.*W;Dₓ.*W]-[sum(u[AssimTimeRange,:],dims=2).*W;u[AssimTimeRange,4].*W])/J₀]
println(" with α= $α  and $Iterations Iterations ")
for i=1:Iterations
    global u=corona.forward(β,γ,δ,u₀,AssimTimeSpan)
    global v=corona.backward(u,values(uₓ[AssimTime]),β,δ,γ,reverse(AssimTimeSpan),W);
    global β,γ,δ
    β,γ,δ,success,Jₑ,αₑ=corona.linesearch(β,γ,δ,v,values(uₓ[AssimTime]),
                                          u₀,AssimTimeSpan,α,LineIterMax,W,
                                          "linesearch")
    if success
        β[nAssim-sponge+1:end].=β[nAssim-sponge]
        γ[nAssim-sponge+1:end].=γ[nAssim-sponge]
        δ[nAssim-sponge+1:end].=δ[nAssim-sponge]
        @save "data/Bavaria/model_parameters.jld"  β γ δ
        corona.save(:Bavaria,AssimTime,u,v,Data,β,γ,δ,"data/Bavaria/solution.jld")
    else
        β,γ,δ,success,Jₑ,αₑ=corona.linesearch(
            β,γ,δ,v,values(uₓ[AssimTime]),
            u₀,AssimTimeSpan,α,LineIterMax,W,
            "brute_force")
    end


    u=corona.forward(β,γ,δ,u₀,AssimTimeSpan)
    if mod(i,ScreenFreq) == 0
        global J=[J; Jₑ/J₀]
        if J[end]==minimum(J)
            color=:green 
        elseif J[end]>minimum(J)
            color=:red
        else
            color=:blue
        end
        print(i," ")
        printstyled(J[end],"\n";color=color)
    end
    if mod(i,FilterFreq) == 0
        println("Filtering")
        β=∂⁰(β)
        γ=∂⁰(γ)
        δ=∂⁰(δ)
        corona.save(:Bavaria,AssimTime,u,v,Data,β,γ,δ,"data/Bavaria/solution.jld")
        @save "data/Bavaria/model_parameters.jld"  β γ δ
    end

    U=TimeArray(collect(AssimTime),u,["S", "I" ,"R" ,"D"])
    if mod(i,PlotFreq) == 0
        P=corona.plot_solution(U,Data,"Bavaria")
        savefig(P,"figs/Bavaria/Development.pdf")
        P=corona.plot_solution(U,Data,"Bavaria",:log10)
        savefig(P,"figs/Bavaria/Developmentlog.pdf")

        pP=scatter(DataTime,β[DataTimeRange],label=L"\beta",legend=:left,
                lw=3,title="Bavaria",color=:red,tickfontsize=12)
        for action in eachrow(Actions)
            d=action[1]
            plot!([d,d],[minimum(β[DataTimeRange]),maximum(β[DataTimeRange])] ,label=action[2],lw=3)
        end
        savefig(pP,"figs/Bavaria/beta.pdf")
        pP=plot!(DataTime,γ[DataTimeRange],label=L"\gamma",lw=3,color=:green)
        pP=plot!(DataTime,δ[DataTimeRange],label=L"\delta",lw=3,color=:black)
        display(pP)
        savefig(pP,"figs/Bavaria/parameters3.pdf")
        debuggP=plot(β[AssimTimeRange[1]:AssimTimeRange[end]+30],label=L"\beta",legend=:left,
                     lw=3,title="Bavaria",color=:red,tickfontsize=12)


   end







end


