busing corona
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

Cases=readtimearray("data/bayern.csv",header=false)
Confirmed=Cases[:B]
Deaths=Cases[:C]
DataTime=timestamp(Cases)
###########################################
S₀=13076721.0
u₀=[S₀,values(Confirmed[1])[1],0.0,values(Deaths[1])[1]]
tspan=(0.0,Float64(length(DataTime))-1)
NDays=400
β₀=1/4*ones(NDays+10)
γ₀=1/7*ones(NDays+10)
δ₀=0.001*1/7*ones(NDays+10);
uₓ=[values(Confirmed) values(Deaths)]
##
β,γ,δ=β₀,γ₀,δ₀
u=corona.forward(β,γ,δ,u₀,tspan)
Cₓ=uₓ[:,1]
Dₓ=uₓ[:,2]
J=[norm([Cₓ;Dₓ]-[u[:,2]+u[:,3]+u[:,4];u[:,4]])/norm([Cₓ;Dₓ])]
plot(DataTime,u[:,2],lw=3,color=:red)
plot!(DataTime,u[:,3],lw=3,color=:green)
plot!(DataTime,u[:,4],lw=3,color=:black)
scatter!(DataTime,uₓ[:,1],color=:orange)
scatter!(DataTime,u[:,4],lw=3,color=:black)


α=1.0e-12
#@load "data/Bavaria_model_parameters.jld" 
for i=1:50000
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
                         label="Confirmed",color=:orange,title="Bavaria")
        P=scatter!(DataTime,uₓ[:,2],label="Deaths",color=:black,lw=3)
        P=plot!(DataTime,sum(u[:,2:4],dims=2),label="C",color=:orange,lw=3)
        P=plot!(DataTime,u[:,2],label="I",color=:red,lw=3)
        P=plot!(DataTime,u[:,3],label="R",color=:green,lw=3)
        P=plot!(DataTime,u[:,4],label="D",color=:black,lw=3)
        display(P)
        savefig("figs/Bavaria.pdf")
        global lP=scatter(DataTime,uₓ[:,1],legend=:topleft,
                         label="Confirmed",color=:orange,title="Bavaria",yaxis=:log10)
        lP=scatter!(DataTime,uₓ[:,2] .+1,label="Deaths",color=:black,lw=3)
        lP=plot!(DataTime,sum(u[:,2:4],dims=2),label="C",color=:orange,lw=3)
        lP=plot!(DataTime,u[:,2] .+1,label="I",color=:red,lw=3)
        lP=plot!(DataTime,u[:,3] .+1.0,label="R",color=:green,lw=3)
        lP=plot!(DataTime,u[:,4] .+1.0,label="D",color=:black,lw=3)
        savefig(lP,"figs/Bavaria_log.pdf")
        @save "data/Bavaria_model_parameters.jld"  β γ δ J
        pP=plot(DataTime,β[1:length(DataTime)],label=L"\beta",legend=:left,
                lw=3,title="Bavaria",color=:red)
        pP=plot!(DataTime,γ[1:length(DataTime)],label=L"\gamma",lw=3,color=:green)
        pP=plot!(DataTime,δ[1:length(DataTime)],label=L"\delta",lw=3,color=:black)
        display(pP)
        savefig(pP,"figs/Bavaria_parameters3.pdf")
   end
end
SIM=DataFrame()
SIM.dates=collect(DataTime)
SIM.S=u[:,1]
SIM.I=u[:,2]
SIM.R=u[:,3]
SIM.D=u[:,4]
SIM.β=β[1:length(DataTime)]
SIM.γ=γ[1:length(DataTime)]
SIM.δ=δ[1:length(DataTime)]
Bavaria=TimeArray(SIM,timestamp=:dates) 
@save "data/Bavaria_final.jld" Bavaria
@save "data/Bavaria_model_parameters.jld" β γ δ
savefig(P,"figs/Bavaria.pdf")
savefig(lP,"figs/Bavaria_log.pdf")
n=(1:length(J))*100
plot(n,J,xlabel="Iterations",ylabel="J",yaxis=:log10,leg=false,lw=3,title="Bavaria")
savefig("figs/Bavaria_convergence.pdf")
pP=plot(DataTime,β[1:length(DataTime)],label=L"\beta",legend=:left,
        lw=3,title="Bavaria",color=:red)
pP=plot!(DataTime,γ[1:length(DataTime)],label=L"\gamma",lw=3,color=:green)
pP=plot!(DataTime,δ[1:length(DataTime)],label=L"\delta",lw=3,color=:black)
display(pP)
savefig(pP,"figs/Bavaria_parameters3.pdf")
pP=plot(DataTime,β[1:length(DataTime)],label=L"\beta",legend=:left,
        lw=3,title="Bavaria",color=:red)
pP=plot!(DataTime,γ[1:length(DataTime)],label=L"\gamma",lw=3,color=:green)
savefig(pP,"figs/Bavaria_parameters.pdf")
σ=β./γ;
pS=plot(DataTime,σ[1:length(DataTime)],label=L"\sigma",legend=:left,
        lw=3,title="Bavaria",color=:red)
savefig(pS,"figs/Bavaria_Sigma.pdf")
μ=δ./γ
pM=plot(DataTime,μ[1:length(DataTime)],label=L"\mu",legend=:left,
        lw=3,title="Bavaria",color=:black)
savefig(pM,"figs/Bavaria_Mortality.pdf")







