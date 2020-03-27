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
S₀=Float64(Population_0[Population_0[!,:Column1].=="United Kingdom",:][!,:Column2][1])
Today=Dates.today()
FirstDate=timestamp(Deaths)[1]
FirstDate=Date(2020,1,31)
LastDate=timestamp(Deaths)[end]
DataTime=FirstDate:Day(1):LastDate
u₀=[S₀,Float64(values(Confirmed[Symbol("United Kingdom")][FirstDate])[1]),0,1e-3]
tspan=(0.0,Float64(length(DataTime))-1)
NDays=400
β₀=1/4*ones(NDays+10)
γ₀=1/7*ones(NDays+10)
δ₀=0.0001*1/7*ones(NDays+10);
uₓ=values(merge(Confirmed[Symbol("United Kingdom")],Deaths[Symbol("United Kingdom")])[DataTime]);
##
SimulRange=FirstDate:Day(1):FirstDate+Day(NDays)
tsimul=(0.0,Float64(NDays))
u₈=corona.forward(β₀,γ₀,δ₀,u₀,tsimul)
Simul=plot(DataTime, uₓ,lw=3)
Simul=plot(SimulRange, u₈)
plot(Confirmed,legend=:topleft,lw=3)
savefig("figs/Confirmed.pdf")
plot(Trend,lw=3)
savefig("figs/Trend_avg.pdf")
##
β,γ,δ=β₀,γ₀,δ₀
γ=γ₀/2
u=corona.forward(β,γ,δ,u₀,tspan)
Cₓ=uₓ[:,1]
Dₓ=uₓ[:,2]
J=[norm([Cₓ;Dₓ]-[u[:,2]+u[:,3]+u[:,4];u[:,4]])/norm([Cₓ;Dₓ])]
α=2.0e-13


#@load "United Kingdom_model_parameters.jld" 
for i=1:9000
    global u
    a=corona.backward(u,uₓ,β,δ,γ,reverse(tspan))
    global β,γ,δ
    β,γ,δ=corona.linesearch(β,γ,δ,a,uₓ,u₀,tspan,α,100)
    β=∂⁰(β)
    u=corona.forward(β,γ,δ,u₀,tspan)
    if mod(i,100) == 0
        I=u[:,2]
        R=u[:,3]
        D=u[:,4]
        Cₓ=uₓ[:,1]
        Dₓ=uₓ[:,2]
        Iₓ=Cₓ-R-Dₓ
        global J=[J; norm([Cₓ;Dₓ]-[u[:,2]+u[:,3]+u[:,4];u[:,4]])/norm([Cₓ;Dₓ])]
        println(i," ",J[end])
        global P=scatter(DataTime,uₓ[:,1],legend=:topleft,
                         label="Confirmed",color=1,title="United Kingdom")
        P=scatter!(DataTime,uₓ[:,2],label="Deaths",color=2,lw=3)
        P=plot!(DataTime,sum(u[:,2:4],dims=2),label="C",color=1,lw=3)
        P=plot!(DataTime,u[:,2],label="I",color=3,lw=3)
        P=plot!(DataTime,u[:,3],label="R",color=4,lw=3)
        P=plot!(DataTime,u[:,4],label="D",color=2,lw=3)
        display(P)
        savefig("figs/United_Kingdom.pdf")
        global lP=scatter(DataTime,uₓ[:,1],legend=:topleft,
                         label="Confirmed",color=1,title="United Kingdom",yaxis=:log10)
        lP=scatter!(DataTime,uₓ[:,2] .+1,label="Deaths",color=2,lw=3)
        lP=plot!(DataTime,sum(u[:,2:4],dims=2),label="C",color=1,lw=3)
        lP=plot!(DataTime,u[:,2] .+1,label="I",color=3,lw=3)
        lP=plot!(DataTime,u[:,3] .+1.0,label="R",color=4,lw=3)
        lP=plot!(DataTime,u[:,4] .+1.0,label="D",color=2,lw=3)
        display(lP)
        savefig(lP,"figs/United_Kingdom_log.pdf")
        @save "data/United Kingdom_model_parameters.jld" DataTime β γ δ J
   end
end
savefig(P,"figs/United_Kingdom.pdf")
savefig(lP,"figs/United_Kingdom_log.pdf")
n=(1:length(J))*100
plot(n,J,xlabel="Iterations",ylabel="J",yaxis=:log10,leg=false,lw=3,title="United Kingdom")
savefig("figs/United_Kingdom_convergence.pdf")
pP=plot(DataTime,β[1:length(DataTime)],label=L"\beta",legend=:left,
        lw=3,title="United Kingdom")
pP=plot!(DataTime,γ[1:length(DataTime)],label=L"\gamma",lw=3)
pP=plot!(DataTime,δ[1:length(DataTime)],label=L"\delta",lw=3)
savefig(pP,"figs/United_Kingdom_parameters3.pdf")
pP=plot(DataTime,β[1:length(DataTime)],label=L"\beta",legend=:left,
        lw=3,title="United Kingdom")
pP=plot!(DataTime,γ[1:length(DataTime)],label=L"\gamma",lw=3)
savefig(pP,"figs/United_Kingdom_parameters.pdf")
σ=β./γ;
pS=plot(DataTime,σ[1:length(DataTime)],label=L"\sigma",legend=:left,
        lw=3,title="United Kingdom")
savefig(pS,"figs/United_Kingdom_Sigma.pdf")








