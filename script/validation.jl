using corona
using Plots
using DataFrames
using CSV
using Dates
using LinearAlgebra
rootdir="/home/jls/prog/corona/"
Population_0=CSV.read(rootdir*"data/population.csv",header=false);
S₀=Float64.(Population_0[!,2][2])
FirstDate=Date(2020,1,28)
NDays=400
SimulRange=FirstDate:Day(1):FirstDate+Day(NDays)
Today=Dates.today()
u₀=[S₀,2.0,0.0,0.0]
β₀=1/4*ones(NDays+10)
γ₀=1/7*ones(NDays+10)
δ₀=0.04*1/7*ones(NDays+10);
tsimul=(0.0,Float64(NDays))
u₈=corona.forward(β₀,γ₀,δ₀,u₀,tsimul)
plot(u₈)
###########################
du=zeros(4)
u=u₈[20,:]
v=rand(4)

t=20.2
ϵ=1e-12
fₑ=copy(corona.sir!(du,u+ϵ*v,[β₀ γ₀ δ₀ ],t));
f₀=copy(corona.sir!(du,u,[β₀ γ₀ δ₀ ],t));
Δ=(fₑ.-f₀)/ϵ

S=u₈[:,1];
I=u₈[:,2];
R=u₈[:,3];
D=u₈[:,4];
A=corona.dfdu([S I R D β₀[1:NDays+1] γ₀[1:NDays+1] δ₀[1:NDays+1]] ,t)
Δ-A*v
corona.dfda(u₈)


### TEST
J=[]
α=1e-12
β=β₀ .+ 0.02*(rand(length(β₀)) .-0.5)
γ=γ₀
δ=δ₀
tsimul=(0.0,30.0)
u₈=corona.forward(β₀,γ₀,δ₀,u₀,tsimul)
u=corona.forward(β,γ,δ,u₀,tsimul)
plot(u₈[:,2:4],legend=:topleft,lw=3)
plot!(u[:,2:4])
uₓ=u₈[:,3:4]
uₓ[:,1]=u₈[:,2]+u₈[:,3]+u₈[:,4]
plot!(uₓ,ls=:dash,lw=3)
v=corona.backward(u,uₓ,β,δ,γ,reverse(tsimul))

J=[0.0]
for i=1:10000
    println("$i ",J[end])
    α=1e-9
    global     β,γ,δ
    u=corona.forward(β,γ,δ,u₀,tsimul)
    v=corona.backward(u,uₓ,β,δ,γ,reverse(tsimul))
    β,γ,δ=corona.linesearch(β,γ,δ,v,uₓ,u₀,tsimul,α,100)
    u=corona.forward(β,γ,δ,u₀,tsimul)
    v=corona.backward(u,uₓ,β,δ,γ,reverse(tsimul))
    global    J=[J;[norm(β-β₀)/norm(β₀)]]
end
plot(J[2:end])
plot(u₈[:,3:4],legend=:topleft,lw=3)
plot(u[:,4])
u[3,4]

plot(u[:,2:4])
plot((δ))
plot(β)
plot!(γ)

