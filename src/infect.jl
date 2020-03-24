using DataFrames
using CSV
using Plots
using LinearAlgebra
using Dates
using TimeSeries
using corona
# %%

Confirmed=corona.read("Confirmed")
Recovered=corona.read("Recovered")
Deaths=corona.read("Deaths")
DateRange=datarange(Confirmed)
UK_confirmed=extract(Confirmed,"United Kingdom")
Germany_confirmed=extract(Confirmed,"Germany")
Germany_recovered=extract(Recovered,"Germany" )
Germany_deaths=extract(Deaths,"Germany" )
Spain_confirmed=reshape(convert(Matrix,filter(g -> g[2] =="Spain", Confirmed)[5:end]),:,1)
Italy_confirmed=reshape(convert(Matrix,filter(g -> g[2] =="Italy", Confirmed)[5:end]),:,1)
Italy_recovered=reshape(convert(Matrix,filter(g -> g[2] =="Italy", Recovered)[5:end]),:,1)
Italy_deaths=reshape(convert(Matrix,filter(g -> g[2] =="Italy", Deaths)[5:end]),:,1)
Ch_confirmed=reshape(convert(Matrix,filter(g -> g[2] =="Switzerland", Confirmed)[5:end]),:,1)
France_confirmed=convert(Matrix,filter(g -> g[2] =="France", Confirmed)[5:end])
France_confirmed=reshape(sum(France_confirmed,dims=1),:,1)
China_confirmed=convert(Matrix,filter(g -> g[2] == "China", Confirmed)[5:end])
China_confirmed=reshape(sum(China_confirmed,dims=1),:,1)
US_confirmed=convert(Matrix,filter(g -> g[2] =="US", Confirmed)[5:end])
US_confirmed=reshape(sum(US_confirmed,dims=1),:,1)
Iran_confirmed=reshape(convert(Matrix,filter(g -> g[2] =="Iran", Confirmed)[5:end]),:,1)
Cases=TimeArray(DateRange,[Germany_confirmed UK_confirmed Spain_confirmed Italy_confirmed France_confirmed US_confirmed Iran_confirmed Ch_confirmed China_confirmed])
Deceased=TimeArray(DateRange,[Germany_deaths Italy_deaths])
Country=["D","UK","E","I","F","US","Iran","Ch","C"];
Cases=TimeSeries.rename(Cases,Country)
RecoveredCases=TimeArray(DateRange,Italy_recovered)
D= values(Cases.D[end])[1]
UK=values(Cases.UK[end])[1]
E=values(Cases.E[end])[1]
It=values(Cases.I[end])[1]
F=values(Cases.F[end])[1]
US=values(Cases.US[end])[1]
Iran=values(Cases.Iran[end])[1]
C=values(Cases.C[end])[1]
Ch=values(Cases.Ch[end])[1]



# %%
function growth(u,Δ=10,Ω=0)
#    (
#    log10(u[end])-log10(u[end-Δ])
#    )/Δ
    t=0:Δ
    P=[ones(length(t)) t]
    c=P\log10.(u[end - Ω .- t])
    -c[2]
end
# %%
function trend(C,N=30)
    g=zeros(N+1)
    ΔD=timestamp(C[end-N])[1]:Day(1):timestamp(C[end])[1]
    for    Ω=0:N
        g[end-Ω] = growth(values(C),7,Ω)
    end
    return TimeArray(ΔD,log10(2)./g,colnames(C))
end

pv=plot(trend(Cases[:I],29),label="I",lw=3,yaxis=(0,5),legend=:topleft,ylabel="Verdopplung")
pv=plot!(trend(Cases[:D],23),label="D",lw=3,yaxis=(0,5))
pv=plot!(trend(Cases[:F],23),label="F",lw=3,yaxis=(0,5))
pv=plot!(trend(Cases[:E],25),label="E",lw=3,yaxis=(0,5))
pv=plot!(trend(Cases[:US],18),label="US",lw=3,yaxis=(0,5))
pv=plot!(trend(Cases[:Ch],26),label="CH",lw=3,yaxis=(0,5))
savefig(pv,"Verdopplungsverlauf.pdf")

pD=plot(1.0 .+ Deceased,yaxis=:log10,legend=:bottomleft,lw=3)
savefig(pD,"Tote.pdf")
pE=plot(1.0 .+ Cases,yaxis=:log10,legend=:bottomleft,lw=3)

savefig(pE,"Fälle.pdf")

summa=plot(pE,pv);
savefig(summa,"Summa.pdf")

plot(Cases.I)
RS=TimeArray(DateRange-Day(12),values(RecoveredCases.A))
plot!(RS)
plot!(DateRange,1000*ones(length(DateRange)))
savefig("Italy_Recovery_delay.pdf")
cd(opw)

