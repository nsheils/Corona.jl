__precompile__()
module corona

using CSV
using DataFrames
using Dates
using TimeSeries
using DifferentialEquations
using LinearAlgebra
using JLD2
using Plots
using LaTeXStrings
using FileIO
function read(dataset::String,data="raw/JHU/COVID-19/csse_covid_19_data/csse_covid_19_time_series/")
    if dataset=="Confirmed"
        #        Confirmed=CSV.read(data*"time_series_19-covid-Confirmed.csv");
        Confirmed=CSV.read(data*"time_series_covid19_confirmed_global.csv");
    elseif dataset=="Recovered"
        Recovered=CSV.read(data*"time_series_19-covid-Recovered.csv")
    elseif dataset=="Deaths"
        Deaths=CSV.read(data*"time_series_covid19_deaths_global.csv")
    else
        Unknown=CSV.read(data*dataset)
    end
end
function datarange(ds::DataFrame)
    FirstDate=Date(String(names(ds)[5]),"m/d/y")+Dates.Year(2000)
    LastDate=Date(String(names(ds)[end]),"m/d/y")+Dates.Year(2000)
    DataRange=FirstDate:Day(1):LastDate
end

function extract(df::DataFrame,State::String,Province="sum")
    DS=convert(Matrix,filter(g -> g[2] ==State, df)[5:end])
    nprovinces,days=size(DS)
    if nprovinces == 1
        return reshape(DS,:,1)
    elseif nprovinces > 1 || Province != "sum"
        return reshape(sum(DS,dims=1),:,1)
    else
        return DS
    end
end

function merge_datasets(df::DataFrame,names::Array{String,1})
    table=extract(df,names[1])
    for i=2:length(names)
        table=[table extract(df,names[i])]
    end
    return TimeArray(datarange(df),table,names)
end

function population(Country=:Italy::Symbol,rootdir="raw/")
    Population=CSV.read(rootdir*"population.csv",header=false);
    Float64(Population[Population[!,:Column1].==String(Country),:][!,:Column2][1])
end
function confirmed(Country=:Italy::Symbol,rootdir="raw/")
    Population=CSV.read(rootdir*"population.csv",header=false);
    Countries=String.(Population[:Column1])
    Confirmed=TimeArray(Float64.(corona.merge_datasets(corona.read("Confirmed"),Countries))[Country])
    TimeSeries.rename(Confirmed,Country => :Confirmed)
end
function deaths(Country=:Italy::Symbol,rootdir="raw/")
    Population=CSV.read(rootdir*"population.csv",header=false);
    Countries=String.(Population[:Column1])
    Deaths=TimeArray(Float64.(corona.merge_datasets(corona.read("Deaths"),Countries))[Country])
    TimeSeries.rename(Deaths,Country => :Deaths)
end
function outbreak(Country=:Italy::Symbol,rootdir="raw/")
    D=confirmed(Country,rootdir)
    timestamp(D[D.>0])[1]
end

struct data
    name::Symbol
    outbreak::Date
    population::Float64
    TimeRange::StepRange{Date,Day}
    TimeIntervall::Tuple{Float64,Float64}
    TimeIndexRange::UnitRange{Int64}
    cases::TimeArray
end
function data(S=:Italy::Symbol)

    C=corona.confirmed(S)
    D=corona.deaths(S)
    FirstDate=timestamp(C)[1]
    LastDate=timestamp(C)[end]
    Range=FirstDate:Day(1):LastDate
    nData=length(Range)
    Intervall=(0.0,Float64(nData-1))
    IndexRange=1:nData
    
    data(S,
         corona.outbreak(S),
         corona.population(S),
         Range,
         Intervall,
         IndexRange,
         merge(C,D))
end
import Base.+
function +(D::data,x::Float64)
    data(D.name,
         D.outbreak,
         D.population,
         D.TimeRange,
         D.TimeIntervall,
         D.TimeIndexRange,
         merge(D.cases[:Confirmed] .+ x,D.cases[:Deaths] .+ x))
end
import Base.length
length(D::data)=length(D.cases)
import Base.values
import TimeSeries.timestamp
values(D::data)=values(D.cases)
export timestamp
timestamp(D::data)=timestamp(D.cases)
@recipe function f(D::corona.data)
    linecolor --> [:orange :black]
    linewidth --> 3
    D.cases
end
function growth(u,p,t=timestamp(u))
    S=values(u[:S])
    I=values(u[:I])
    R=values(u[:R])
    D=values(u[:D])
    N=S+R+I+D
    β=p[:β]
    γ=p[:γ]
    δ=p[:δ]
    λ=values(β) .*S ./N  - (values(γ)+values(δ))
    λ=TimeArray(t,λ[1:length(t)],["λ"])
end
"Classic Epidemology Model  Hethcote (2000), added deaths"
function sir!(du::Array{Float64,1},u::Array{Float64,1},p::Array{Float64,2},t::Float64)
    β,γ,δ=p[Int64(floor(t)+1),1:3]
#     β,γ,δ=lagrange(p[:,1:3],t)
##
    S,I,R,D=u[1:4]
    N=S+I+R+D
    du[1]=-β*I*S/N
    du[2]= β*I*S/N - (γ+δ)*I
    du[3]= γ*I
    du[4]= δ*I
    du
end

function dfdu(p::Array{Float64,2},t::Float64)
    S,I,R,D,β,γ,δ=p[Int64(floor(t)+1),1:7]
#     S,I,R,D,β,γ,δ=lagrange(p[:,1:7],t,0)
    N=S+I+R+D
    A=zeros(4,4)
    A[1,1]=-β*I*(N+S)/N^2
    A[1,2]=-β*S*(N+I)/N^2
    A[1,3]=+β*I*S/N^2
    A[1,4]=+β*I*S/N^2

    A[2,1]= β*I*(N+S)/N^2
    A[2,2]= β*S*(N+I)/N^2 - (γ+δ)
    A[2,3]=-β*I*S/N^2
    A[2,4]=-β*I*S/N^2

    A[3,2]=γ
    A[4,2]=δ
    return A
end


function sir_adj(v::Array{Float64,1},p::Array{Float64,2},t::Float64)
    S,I,R,D,β,γ,δ,Cₓ,Dₓ,window=p[Int64(floor(t)+1),1:10]
#    S,I,R,D,β,γ,δ,Cₓ,Dₓ,window=lagrange(p[:,1:10],t,0)
    A=dfdu(p,t)
    # Observations
    OBS=[0 0;1 0;1 0;1 1]
    # weights
    g=((OBS'*[S,I,R,D] - [Cₓ,Dₓ])' * OBS')'

    dv= -A'*v -g.*window

    return dv
end

function dfda(u::Array{Float64,2})
    S=u[:,1]
    I=u[:,2]
    R=u[:,3]
    D=u[:,4]
    N=S+I+R+D
    f=zeros(length(S),4,3)
    f[:,1,1]=-I.*S./N
    f[:,1,2].= 0.0
    f[:,1,3].= 0.0

    f[:,2,1]= I.*S./N
    f[:,2,2]=-I
    f[:,2,3]=-I

    f[:,3,1].=0.0
    f[:,3,2]= I
    f[:,3,3].=0.0

    f[:,4,1].=0.0
    f[:,4,2].=0.0
    f[:,4,3]= I

    return f
end


function forward(β,γ,δ,u₀,tspan)
    p=[β γ δ]
    problem=ODEProblem(sir!,u₀,tspan,p,maxiters=1e6)
    solution=solve(problem)
    u=collect(solution(0:tspan[2])')
end
function backward(u,uₓ,β,γ,δ,tspan=(50.0,0.0),window=ones(size(u)[1]))
    nt,ne=size(u)
    p=[u β[1:nt] γ[1:nt] δ[1:nt] uₓ window]
    v₀=[0.0 ,0.0, 0.0, 0.0]
    adj=ODEProblem(corona.sir_adj,v₀,tspan,p,maxiters=1e6)
    adj_sol=solve(adj)
    v=collect(adj_sol(0:tspan[1])')
end
function update(u,v,α,β,γ,δ)
    if α==0.0
        return  β,γ,δ
    else
        nt,ne=size(u)
        f=dfda(u)
        δβ=zeros(nt)
        δγ=zeros(nt)
        δδ=zeros(nt)
        for i=1:4
            δβ=δβ+ v[:,i] .* f[:,i,1]
            δγ=δγ+ v[:,i] .* f[:,i,2]
            δδ=δδ+ v[:,i] .* f[:,i,3]
        end
        βₙ=copy(β)
        γₙ=copy(γ)
        δₙ=copy(δ)
        βₙ[1:nt]=max.(β[1:nt]-α*δβ,0.0)
        γₙ[1:nt]=max.(γ[1:nt]-α*δγ,0.0)
        δₙ[1:nt]=max.(δ[1:nt]-α*δδ,0.0)
        return  βₙ,γₙ,δₙ
    end
end

function linesearch(β,γ,δ,v,uₓ,u₀,tspan,
                    α=1.0e-12,itMax=2048,window=ones(size(u)[1]),
                    method="bisection")
    u=corona.forward(β,γ,δ,u₀,tspan)
    function probe(α)
        β₁,γ₁,δ₁=update(u,v,α,β,γ,δ)
        u₁=forward(β₁,γ₁,δ₁,u₀,tspan)
        I=u₁[:,2]
        R=u₁[:,3]
        D=u₁[:,4]
        Cₓ=uₓ[:,1]
        Dₓ=uₓ[:,2]
        C=I+R+D
        norm([C .* window;D .* window]-[Cₓ .* window;Dₓ .* window])
    end
    function brute_force(αₐ,αₑ,nΔ=itMax)
        J₀=probe(0.0)
        Δα=(αₑ-αₐ)
        Δ=αₐ .+ Δα*sort(rand(nΔ+1))
        J=zeros(nΔ+1,2)
        for  i=1:nΔ+1
            J[i,:]=[Δ[i] probe(Δ[i])]
        end
        αₘ=J[argmin(J[:,2]),1]
        J₁=minimum(J[:,2])
        β₁,γ₁,δ₁=update(u,v,αₘ,β,γ,δ)
        return β₁,γ₁,δ₁,true,J₁,αₘ
    end
    function bisection(αₐ= 0.0,αₑ= α,itMax=itMax)
        i=1
        ϵ = 1/3
        J₀=probe(0.0)
        Jₑ=probe(αₑ)
        if Jₑ<J₀
            β₁,γ₁,δ₁=update(u,v,αₑ,β,γ,δ)
            return β₁,γ₁,δ₁,true,Jₑ,αₑ
        end
        Jₐ=J₀
        while i < itMax
            #        αₘ=(αₐ+αₑ)/2
            Δα=(αₑ-αₐ)
            αₘ=αₐ+ rand()*Δα
            Jₘ=probe(αₘ)
            if (Jₘ<J₀) & (Δα/α<ϵ)
                β₁,γ₁,δ₁=update(u,v,αₘ,β,γ,δ)
                return β₁,γ₁,δ₁,true,Jₘ,αₘ
            end
            if Jₘ < Jₐ
                αₐ=αₘ
                Jₐ=Jₘ
            else
                αₑ=αₘ
                Jₑ=Jₘ
            end
            i=i+1
        end
        #    println("$i failed")
        return β,γ,δ,false,J₀,0.0
    end
    ###########################
    # calculation start
    if method=="bisection"
        β₁,γ₁,δ₁,success,Jfinal,α₁ = bisection(0.0,α,itMax)
    else
        β₁,γ₁,δ₁,success,Jfinal,α₁ = brute_force(0.0,α,itMax)
    end
end
function extrapolate(y,Δx=1.0)
    y⁰=y[end]
    y¹= [1/2,-2,3/2]'*y[end-2:end]
    y²= [-1, 4,-5,2]'*y[end-3:end]
    y⁰ + y¹*Δx + + y²*Δx^2/2
end

function smooth(y,Δx=1.0)
    f⁰=(circshift(y,-1) +   2*y +circshift(y,1))/4
    f⁰[1]=y[1]
    f⁰[end]=y[end]
    return f⁰
end


function save(Region,AssimTime,U,V,Uₓ,β,γ,δ,filename)
    nAssim=length(AssimTime)
    AssimTimeRange=1:nAssim
    u=TimeArray(collect(AssimTime),U,["S", "I" ,"R" ,"D"])
    v=TimeArray(collect(AssimTime),V,["S", "I" ,"R" ,"D"])
    uₓ=Uₓ[AssimTime]
    p=TimeArray(collect(AssimTime),[ β[AssimTimeRange] γ[AssimTimeRange] δ[AssimTimeRange] ],["β", "γ" ,"δ"])
    FileIO.save(filename,
                "Region",Region,
                "AssimTime",AssimTime,
                "u",u,
                "v",v,
                "p",p,
                "uₓ",uₓ)
end




function plot_solution(u,uₓ,title="Corona growth",yscale=:identity)
    I=u[:I]
    R=u[:R]
    D=u[:D]
    C=I .+R .+D
    if yscale==:identity
        P=scatter(uₓ[:Confirmed],legend=:topleft,color=:orange,title=title)
        P=scatter!(uₓ[:Deaths],label="Deaths",color=:black,lw=3,tickfontsize=12)
        P=plot!(C,label="C",color=:orange,lw=3)
        P=plot!(I,label="I",color=:red,lw=3)
        P=plot!(R,label="R",color=:green,lw=3)
        P=plot!(D,label="D",color=:black,lw=3)
    elseif yscale==:log10
        P=scatter(uₓ[:Confirmed].+1,legend=:topleft,
                  color=:orange,title=title,yaxis=:log10)
        P=scatter!(uₓ[:Deaths] .+1,label="Deaths",color=:black,lw=3,tickfontsize=12)
        P=plot!(C,label="C",color=:orange,lw=3)
        P=plot!(I,label="I",color=:red,lw=3)
        P=plot!(R.+1,label="R",color=:green,lw=3)
        P=plot!(D.+1,label="D",color=:black,lw=3)
    end
    return P
end


"Lagrange Interpolation, t=0:n"
function lagrange(y::Array{Float64,1},t::Float64,order=2)
    tᵢ = floor(t)  ;nᵢ =Int(tᵢ )+1
    if order == 0
        return y[nᵢ]
    end
    tᵢ₊= floor(t+1);nᵢ₊=Int(tᵢ₊)+1
    tᵢ₋= floor(t-1);nᵢ₋=Int(tᵢ₋)+1
    if nᵢ₋ <  1
        yᵢ =y[nᵢ]
        yᵢ₊=y[nᵢ₊]
        return yᵢ + (t-tᵢ) *  (yᵢ₊-yᵢ)/(tᵢ₊-tᵢ)
    elseif nᵢ >= length(y)
        yᵢ =y[end]
        yᵢ₋=y[end-1]
        tᵢ = Float64(length(y)-1)
        tᵢ₋= Float64(length(y)-2)
        return yᵢ + (t-tᵢ) *  (yᵢ-yᵢ₋)/(tᵢ-tᵢ₋)
    else
        yᵢ =y[nᵢ]
        yᵢ₋=y[nᵢ₋]
        yᵢ₊=y[nᵢ₊]
        l =
            yᵢ₋*(t-tᵢ )/(tᵢ₋-tᵢ ) *(t-tᵢ₊)/(tᵢ₋-tᵢ₊) + 
            yᵢ *(t-tᵢ₋)/(tᵢ -tᵢ₋) *(t-tᵢ₊)/(tᵢ -tᵢ₊) + 
            yᵢ₊*(t-tᵢ₋)/(tᵢ₊-tᵢ₋) *(t-tᵢ )/(tᵢ₊-tᵢ )
    end
end
function lagrange(y::Array{Float64,2},t::Float64,order=2)
    tᵢ = floor(t)  ;nᵢ =Int(tᵢ )+1
    if order == 0
        return y[nᵢ,:]
    end
    tᵢ₊= floor(t+1);nᵢ₊=Int(tᵢ₊)+1
    tᵢ₋= floor(t-1);nᵢ₋=Int(tᵢ₋)+1
    if nᵢ₋ <  1
        yᵢ =y[nᵢ,:]
        yᵢ₊=y[nᵢ₊,:]
        return yᵢ + (t-tᵢ) *  (yᵢ₊-yᵢ)./(tᵢ₊-tᵢ)
    elseif nᵢ >= length(y)
        yᵢ =y[:,end]
        yᵢ₋=y[:,end-1]
        tᵢ = Float64(length(y)-1)
        tᵢ₋= Float64(length(y)-2)
        return yᵢ + (t-tᵢ) *  (yᵢ-yᵢ₋)./(tᵢ-tᵢ₋)
    else
        yᵢ =y[:,nᵢ]
        yᵢ₋=y[:,nᵢ₋]
        yᵢ₊=y[:,nᵢ₊]
        l =
            yᵢ₋*(t-tᵢ )/(tᵢ₋-tᵢ ) *(t-tᵢ₊)/(tᵢ₋-tᵢ₊) + 
            yᵢ *(t-tᵢ₋)/(tᵢ -tᵢ₋) *(t-tᵢ₊)/(tᵢ -tᵢ₊) + 
            yᵢ₊*(t-tᵢ₋)/(tᵢ₊-tᵢ₋) *(t-tᵢ )/(tᵢ₊-tᵢ )
    end
end




end
