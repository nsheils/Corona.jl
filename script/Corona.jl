using Revise
using corona
using Plots
using TimeSeries
using CSV
using Statistics
using LinearAlgebra
using JLD2
Countries=["Italy","Germany","United Kingdom","US"]
Population_0=CSV.read("/home/jls/prog/corona/data/population.csv",header=false)
FirstDate=Date(2020,2,1)
LastDate=Date(2020,3,22)
LastDate=Dates.today()-Day(2)
DataTime=FirstDate:Day(1):LastDate
Confirmed=Float64.(corona.merge_datasets(corona.read("Confirmed"),Countries)[DataTime])
Deaths=   Float64.(corona.merge_datasets(corona.read("Deaths"),Countries)[DataTime])
Recovered=Float64.(corona.merge_datasets(corona.read("Recovered"),Countries)[DataTime])


Infected=TimeArray(timestamp(Confirmed),values(Confirmed) .- values(Deaths) .- values(Recovered),Countries)
Population=TimeArray(timestamp(Confirmed),Population_0[!,2]' .- values(Deaths),Countries)
Susceptible=TimeArray(timestamp(Confirmed),values(Population) .- values(Deaths) .- values(Recovered),Countries)
gamma=TimeArray(timestamp(Confirmed)[2:end],values(diff(Recovered) ./ (Infected .+ 1)),Countries)
beta=TimeArray(timestamp(Confirmed)[2:end],values(((diff(Infected) .+ diff(Recovered)).*Population)./(Infected .+1) ./Susceptible),Countries)
γ₀=1/14
β₀=mean(beta[end-14])
#########################################
Ndays=360

TimeSpan=FirstDate:Day(1):FirstDate+Day(Ndays)
β=ones(Ndays+1)*values(β₀[:Germany])[1]
γ=ones(Ndays+1)*γ₀

Param=TimeArray(TimeSpan,[β γ] ,[:β  ,:γ])




tspan=(0.0,Float64((LastDate-FirstDate).value))
S=Susceptible[:Germany][FirstDate:Day(1):LastDate]
I=Infected[:Germany][FirstDate:Day(1):LastDate]
R=Recovered[:Germany][FirstDate:Day(1):LastDate]
uₓ=values(merge(S,I,R))
u₀=uₓ[1,:]

u=corona.forward(β,γ,u₀,tspan)
J=[norm(uₓ-u)]
P=plot(uₓ[:,2],legend=false,yaxis=:log10)
P= plot!(u[:,2])
for i=1:10000
    global u
    a=corona.backward(u,uₓ,β,γ,reverse(tspan))
    global β,γ
    β,γ₀=corona.linesearch(β,γ,a,uₓ,u₀,tspan)
    u=corona.forward(β,γ,u₀)
    global J=[J; norm(uₓ-u)]
    if mod(i,1000) == 0
        println(i," ",J[end])
#        P=plot(uₓ[:,2],legend=:topleft,yaxis=:log10)
        P=plot!(u[:,2])
        display(P)
        @save "Param_Germany.jld" β γ
   end
end
@save "Param_Germany.jld" β γ
Param=TimeArray(TimeSpan,[β γ] ,[:β  ,:γ])
CSV.write("Param_Germany.csv", Param)
conv=plot(J,xaxis="iter",yaxis="J")
P=plot(uₓ[:,2],legend=:topleft,yaxis=:log10)
P=plot!(u[:,2])

