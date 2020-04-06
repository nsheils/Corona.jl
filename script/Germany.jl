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
############################################################
rootdir      = "/home/jls/data/2020-Corona/"
Iterations   = 10000
FilterFreq   = 10
ScreenFreq   = 100
PlotFreq     = 100
sponge       = 60
LineIterMax  = 2^12
LineRangeMax = 1e-6
ColdStart    = false
Region       = :Germany
############################################################
cd(rootdir)
govactfile="raw/GOV/Germany/actions.csv"
if isfile(govactfile) 
    Actions=CSV.read(govactfile)
    println(Actions)
else
    wd=pwd()
    println("No GOV Actions $wd*$govactfile")
end
parameter=DataFrame(
    parameter=["rootdir", "ColdStart", "Iterations", "FilterFreq", "ScreenFreq","PlotFreq","sponge","LineIterMax","LineRangeMax","Region"],
    value=[rootdir, ColdStart, Iterations, FilterFreq, ScreenFreq,PlotFreq,sponge,LineIterMax,LineRangeMax,Region])
println(parameter)
Today=Dates.today()
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
##################################################################################
nDays=400
SimulLast=FirstDate+Day(nDays)
GlobalTime=FirstDate:Day(1):SimulLast
nGlobal=length(GlobalTime)
SimulTime=OutbreakDate:Day(1):SimulLast
nSimul=length(SimulTime)
SimulTimeSpan=(0.0,Float64(nSimul-1))
SimulTimeRange=1:nSimul
##########################################
A=zeros(nGlobal,2);A[DataIndexRange,:]=values(Data)
uₓ=TimeArray(GlobalTime,A,TimeSeries.colnames(Data))
u₀=[S₀,values(D.cases[:Confirmed][OutbreakDate])[1],0.0,values(D.cases[:Deaths][OutbreakDate])[1]]
########################################
if ColdStart == true
    β₀=1/4*ones(nSimul+1)
    γ₀=1/7*ones(nSimul+1)
    δ₀=zeros(nSimul+1);
    printstyled("Cold Start";color=:red)
    β=β₀;γ=γ₀;δ=δ₀
else
    printstyled("Warm Start";color=:green)
    @load "data/Germany/model_parameters.jld"
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
diagnostic_message=true
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
        @save "data/Germany/model_parameters.jld"  β γ δ
#        #corona.save(:Germany,AssimTime,u,v,Data,β,γ,δ,"data/Germany/solution.jld")
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
    if mod(i,FilterFreq) ==  0
        if diagnostic_message == true
            print("Filtering γ,δ\n")
            global diagnostic_message= false
        end
#        β=∂⁰(β)
        γ=∂⁰(γ)
        δ=∂⁰(δ)
        #corona.save(:Germany,AssimTime,u,v,Data,β,γ,δ,"data/Germany/solution.jld")
        @save "data/Germany/model_parameters.jld"  β γ δ
    end

    U=TimeArray(collect(AssimTime),u,["S", "I" ,"R" ,"D"])
    if mod(i,PlotFreq) == 0
        P=corona.plot_solution(U,Data,"Germany")
        savefig(P,"figs/Germany/Development.pdf")
        P=corona.plot_solution(U,Data,"Germany",:log10)
#        for action in eachrow(Actions)
#            d=[action[1],action[1]]
#            v=[1,10000]
#            P=plot!(d,v,lw=3,label=action[2])
#        end

        savefig(P,"figs/Germany/Development_log.pdf")

        pP=scatter(DataTime,β[DataIndexRange],label=L"\beta",legend=:left,
                lw=3,title="Germany",color=:red,tickfontsize=12)
#        for action in eachrow(Actions)
#            d=action[1]
#            plot!([d,d],[minimum(β[DataIndexRange]),maximum(β[DataIndexRange])] ,label=action[2],lw=3)
#        end
        savefig(pP,"figs/Germany/beta.pdf")
        pP=plot!(DataTime,γ[DataIndexRange],label=L"\gamma",lw=3,color=:green)
        pP=plot!(DataTime,δ[DataIndexRange],label=L"\delta",lw=3,color=:black)
        display(pP)
        savefig(pP,"figs/Germany/parameters3.pdf")

   end
end
