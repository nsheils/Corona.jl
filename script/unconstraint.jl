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
Trend=corona.trend(Confirmed,14,40)
Growths=Trend[end]
#
FirstDate=timestamp(Deaths)[1]
LastDate=timestamp(Deaths)[end]
################################################################################
@load "data/outbreak.jld2" Outbreak
OutbreakDate=Outbreak[Symbol("United Kingdom")]
DataTime=OutbreakDate:Day(1):LastDate
nData=length(DataTime)
tspan=(0.0,Float64(nData)-1)
NDays=400
SimulRange=OutbreakDate:Day(1):OutbreakDate+Day(NDays)
tsimul=(0.0,Float64(NDays))
###########################################

u₀=[S₀,Float64(values(Confirmed[Symbol("United Kingdom")][OutbreakDate])[1]),0.0,0.0]

β₀=1/4*ones(NDays+10)
γ₀=1/7*ones(NDays+10)
δ₀=0.001*1/7*ones(NDays+10);
uₓ=zeros(NDays+10)
uₓ[DataTime]=values(merge(Confirmed[Symbol("United Kingdom")],Deaths[Symbol("United Kingdom")])[DataTime]);
##
u₈=zeros(NDays+10)
u₈[DataTime]=corona.forward(β₀,γ₀,δ₀,u₀,tsimul)
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


α=1.0e-12
@load "data/United_Kingdom_model_parameters.jld" 
for i=1:10000
    global u
    a=corona.backward(u,uₓ,β,δ,γ,reverse(tspan))
    global β,γ,δ
    β,γ,δ=corona.linesearch(β,γ,δ,a,uₓ,u₀,tspan,α,100)
    β=∂⁰(β)
    γ=∂⁰(γ)
    δ=∂⁰(δ)
    u=corona.forward(β,γ,δ,u₀,tspan)
    if mod(i,1000) == 0
        I=u[:,2]
        R=u[:,3]
        D=u[:,4]
        Cₓ=uₓ[:,1]
        Dₓ=uₓ[:,2]
        Iₓ=Cₓ-R-Dₓ
        global J=[J; norm([Cₓ;Dₓ]-[u[:,2]+u[:,3]+u[:,4];u[:,4]])/norm([Cₓ;Dₓ])]
        println(i," ",J[end])
        global P=scatter(DataTime,uₓ[:,1],legend=:topleft,
                         label="Confirmed",color=:orange,title="United Kingdom")
        P=scatter!(DataTime,uₓ[:,2],label="Deaths",color=:black,lw=3)
        P=plot!(DataTime,sum(u[:,2:4],dims=2),label="C",color=:orange,lw=3)
        P=plot!(DataTime,u[:,2],label="I",color=:red,lw=3)
        P=plot!(DataTime,u[:,3],label="R",color=:green,lw=3)
        P=plot!(DataTime,u[:,4],label="D",color=:black,lw=3)
        display(P)
        savefig("figs/United_Kingdom.pdf")
        global lP=scatter(DataTime,uₓ[:,1],legend=:topleft,
                         label="Confirmed",color=:orange,title="United Kingdom",yaxis=:log10)
        lP=scatter!(DataTime,uₓ[:,2] .+1,label="Deaths",color=:black,lw=3)
        lP=plot!(DataTime,sum(u[:,2:4],dims=2),label="C",color=:orange,lw=3)
        lP=plot!(DataTime,u[:,2] .+1,label="I",color=:red,lw=3)
        lP=plot!(DataTime,u[:,3] .+1.0,label="R",color=:green,lw=3)
        lP=plot!(DataTime,u[:,4] .+1.0,label="D",color=:black,lw=3)
        savefig(lP,"figs/United_Kingdom_log.pdf")
        FirstDay=OutbreakDate
        @save "data/United_Kingdom_model_parameters.jld" FirstDay β γ δ J
        pP=plot(DataTime,β[1:nData],label=L"\beta",legend=:left,
                lw=3,title="United Kingdom",color=:red)
        pP=plot!(DataTime,γ[1:nData],label=L"\gamma",lw=3,color=:green)
        pP=plot!(DataTime,δ[1:nData],label=L"\delta",lw=3,color=:black)
        display(pP)
        savefig(pP,"figs/United_Kingdom_parameters3.pdf")
   end
end
SIM=DataFrame()
SIM.dates=collect(DataTime)
SIM.S=u[:,1]
SIM.I=u[:,2]
SIM.R=u[:,3]
SIM.D=u[:,4]
SIM.β=β[1:nData]
SIM.γ=γ[1:nData]
SIM.δ=δ[1:nData]
United Kingdom=TimeArray(SIM,timestamp=:dates) 
@save "data/United_Kingdom_final.jld" United Kingdom
@save "data/United_Kingdom_model_parameters.jld" β γ δ
savefig(P,"figs/United_Kingdom.pdf")
savefig(lP,"figs/United_Kingdom_log.pdf")
n=(1:length(J))*100
plot(n,J,xlabel="Iterations",ylabel="J",yaxis=:log10,leg=false,lw=3,title="United Kingdom")
savefig("figs/United_Kingdom_convergence.pdf")
pP=plot(DataTime,β[1:nData],label=L"\beta",legend=:left,
        lw=3,title="United Kingdom",color=:red)
pP=plot!(DataTime,γ[1:nData],label=L"\gamma",lw=3,color=:green)
pP=plot!(DataTime,δ[1:nData],label=L"\delta",lw=3,color=:black)
display(pP)
savefig(pP,"figs/United_Kingdom_parameters3.pdf")
pP=plot(DataTime,β[1:nData],label=L"\beta",legend=:left,
        lw=3,title="United Kingdom",color=:red)
pP=plot!(DataTime,γ[1:nData],label=L"\gamma",lw=3,color=:green)
savefig(pP,"figs/United_Kingdom_parameters.pdf")
σ=β./γ;
pS=plot(DataTime,σ[1:nData],label=L"\sigma",legend=:left,
        lw=3,title="United Kingdom",color=:red)
savefig(pS,"figs/United_Kingdom_Sigma.pdf")
μ=δ./γ
pM=plot(DataTime,μ[1:nData],label=L"\mu",legend=:left,
        lw=3,title="United Kingdom",color=:black)
savefig(pM,"figs/United_Kingdom_Mortality.pdf")







