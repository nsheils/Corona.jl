using corona
using Plots
using JLD2
using CSV
using Dates
using TimeSeries
using LinearAlgebra
using DataFrames
using june
using LaTeXStrings
############################################################
rootdir      = "/home/jls/data/2020-Corona/"
cd(rootdir)
Today=Dates.today()
@load "data/Bavaria/solution.jld" 
Actions=DataFrame(Date=Date(2020,03,16),Action="Schulschließung")
push!(Actions,[Date(2020,03,20),"Vorläufige Ausgangsbeschränkung"])
λ=corona.growth(u,p)
P=plot(log(2)./λ,lw=3,label=L"\lambda",tickfontsize=12)
for action in eachrow(Actions)
    d=[action[1],action[1]]
    v=[minimum(values(log(2)./λ))[1],maximum(values(log(2)./λ))[1]]
    global P=plot!(d,v,lw=3)
end
display(P)
P=plot(λ,lw=3,label=L"\lambda",tickfontsize=12,legend=:bottomleft)
d=[Today,Today]
v=[minimum(values(λ))[1],maximum(values(λ))[1]]
P=plot!(d,v,lw=3,label="today")
for action in eachrow(Actions)
    d=[action[1],action[1]]
    v=[minimum(values(λ))[1],maximum(values(λ))[1]]
    global P=plot!(d,v,lw=3,label=action[2])
end



P=corona.plot_solution(u,uₓ,"Bavaria",:log10)
for action in eachrow(Actions)
    d=[action[1],action[1]]
    v=[1,10000]
    global P=plot!(d,v,lw=3,label=action[2],legend=false)
end
display(P)




