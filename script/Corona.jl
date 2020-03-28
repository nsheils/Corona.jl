using corona
using Plots 
using TimeSeries
using CSV
using Statistics
using LinearAlgebra
using JLD2
Countries=["Italy","Germany","United Kingdom","US"];
Population_0=CSV.read("/home/jls/prog/corona/data/population.csv",header=false);
FirstDate=Date(2020,2,1);
LastDate=Date(2020,3,22);
LastDate=Dates.today()#-Day(2);
DataTime=FirstDate:Day(1):LastDate;
Confirmed=Float64.(corona.merge_datasets(corona.read("Confirmed"),Countries)[DataTime]);
Deaths=   Float64.(corona.merge_datasets(corona.read("Deaths"),Countries)[DataTime]);
Recovered=Float64.(corona.merge_datasets(corona.read("Recovered"),Countries)[DataTime[end-3]]);
Infected=TimeArray(timestamp(Confirmed),values(Confirmed) .- values(Deaths) .- values(Recovered),Countries);
Population=TimeArray(timestamp(Confirmed),Population_0[!,2]' .- values(Deaths),Countries);
Susceptible=TimeArray(timestamp(Confirmed),values(Population) .- values(Deaths) .- values(Recovered),Countries);
#gamma=TimeArray(timestamp(Confirmed)[2:end],values(diff(Recovered) ./ (Infected .+ 1)),Countries);
#beta=TimeArray(timestamp(Confirmed)[2:end],values(((diff(Infected) .+ diff(Recovered)).*Population)./(Infected .+1) ./Susceptible),Countries);
γ₀=1/14
#β₀=mean(beta[end-14]);
β₀=1/4
#########################################
Ndays=400

TimeSpan=FirstDate:Day(1):FirstDate+Day(Ndays);
β=ones(Ndays+1)*β₀;
γ=ones(Ndays+1)*γ₀;
δₒ=1/100

tspan=(0.0,Float64((LastDate-FirstDate).value));
S=Susceptible[:Germany][FirstDate:Day(1):LastDate];
I=Infected[:Germany][FirstDate:Day(1):LastDate];
R=Recovered[:Germany][FirstDate:Day(1):LastDate];
uₓ=values(merge(S,I,R));
u₀=uₓ[1,:]
Parameters="Param_Germany.jld" 
println("load paremeters $Parameters")
@load "Param_Germany.jld" β γ
u=corona.forward(β,γ,u₀,tspan)
J=[norm(uₓ[:,2]-u[:,2])]
P=scatter(uₓ[:,2],legend=false,yaxis=:log10)
P= plot!(u[:,2])
for i=1:50000
    global u
    a=corona.backward(u,uₓ,β,γ,reverse(tspan))
    global β,γ
    β,γ=corona.linesearch(β,γ,a,uₓ,u₀,tspan,1.0e-12)
    u=corona.forward(β,γ,u₀)
    if mod(i,1000) == 0
        global J=[J; norm(uₓ[:,2]-u[:,2])]
        println(i," ",J[end])
        P=scatter(DataTime,uₓ[:,2],legend=false,yaxis=:log10,
        xlabel="Date",ylabel="I (Germany)")
        P=plot!(DataTime,u[:,2])
        display(P)
        @save "Param_Germany.jld" β γ
        @save "u_Germany.jld" u uₓ
    end
end
@save "Param_Germany.jld" β γ
Param=TimeArray(TimeSpan,[β γ] ,[:β  ,:γ])
CSV.write("Param_Germany.csv", Param)
conv=plot(J,xaxis="iter",yaxis="J")

tspan=(0.0,Float64(Ndays+1))
u=corona.forward(β,γ,u₀,tspan)
U=TimeArray(TimeSpan,u)
plot!(U,legend=:topleft)

