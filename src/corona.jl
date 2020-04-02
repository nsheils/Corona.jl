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
    DataTime::StepRange{Date,Day}
    cases::TimeArray
end
function data(S=:Italy::Symbol)
    C=corona.confirmed(S)
    D=corona.deaths(S)
    FirstDate=timestamp(C)[1]
    LastDate=timestamp(C)[end]
    DataTime=FirstDate:Day(1):LastDate
    data(S,corona.outbreak(S),
         corona.population(S),
         DataTime,
         merge(C,D))
end
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
import Base.+
function +(D::data,x::Float64)
    data(D.name,D.outbreak,
         D.population,
         D.DataTime,
         merge(D.cases[:Confirmed] .+ x,D.cases[:Deaths] .+ x))
end

    


# %%
function growth(u,Δ=10,Ω=0)
    t=0:Δ
    P=[ones(length(t)) t]
    if ndims(u) == 1
        c=P\log10.(u[end - Ω .- t])
        -c[2]
    elseif ndims(u) == 2
        c=P\log10.(u[end - Ω .- t,:])
        -c[2,:]
    end
end
# %%
function trend(C,Δ=10,N=30)
    nt,nc=size(C)
    ΔD=timestamp(C[end-N])[1]:Day(1):timestamp(C[end])[1]
    Countries=String.(colnames(C))
    if ndims(C) == 1
        g=zeros(N+1)
        for    Ω=0:N
            g[end-Ω] = growth(values(C),Δ,Ω)
        end
    elseif ndims(C) == 2
        g=zeros(N+1,nc)
        for    Ω=0:N
            g[end-Ω,:] = growth(values(C),Δ,Ω)'
        end
    end
    return  TimeArray(ΔD,g,Countries)
end



"Classic Epedemic Model  Hethcote (2000), added deaths"
function sir!(du::Array{Float64,1},u::Array{Float64,1},p::Array{Float64,2},t::Float64)
    it=Int(floor(t)+1.0)
    Δt=t-floor(t)
    nt,ne=size(p)
    if it == 1
        β,γ,δ=p[it,:] + Δt*(p[it+1,:] -p[it,:])
    elseif it == nt
        β,γ,δ=p[it,:] + Δt*(p[it,:] -p[it-1,:])
    else
        β,γ,δ=p[it,:] + Δt*(p[it+1,:] -p[it-1,:])/2
    end
    S,I,R,D=u[1:4]
    N=S+I+R+D
    du[1]=-β*I*S/N
    du[2]= β*I*S/N - (γ+δ)*I
    du[3]= γ*I
    du[4]= δ*I
    du
end

function dfdu(p::Array{Float64,2},t::Float64)
    it=Int(floor(t)+1)
    Δt=t-floor(t)
    nt,ne=size(p)
    if it == nt
        S,I,R,D,β,γ,δ=p[it,:] + Δt*(p[it,:] -p[it-1,:])
    elseif it == 1
        S,I,R,D,β,γ,δ=p[it,:] + Δt*(p[it+1,:] -p[it,:])
    else
        S,I,R,D,β,γ,δ=p[it,:] + Δt*(p[it+1,:] -p[it-1,:])/2
    end
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
    it=Int(floor(t)+1)
    Δt=t-floor(t)
    nt,ne=size(p)
    if it == nt
        S,I,R,D,β,γ,δ,Cₓ,Dₓ,window=p[it,:] + Δt*(p[it,:] -p[it-1,:])
    elseif it == 1
        S,I,R,D,β,γ,δ,Cₓ,Dₓ,window=p[it,:] + Δt*(p[it+1,:] -p[it,:])
    else
        S,I,R,D,β,γ,δ,Cₓ,Dₓ,window=p[it,:] + Δt*(p[it+1,:] -p[it-1,:])/2
    end
    S,I,R,D,β,γ,δ,Cₓ,Dₓ,window=p[it,:]
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
    problem=ODEProblem(sir!,u₀,tspan,p)
    solution=solve(problem)
    u=collect(solution(0:tspan[2])')
end
function backward(u,uₓ,β,γ,δ,tspan=(50.0,0.0),window=ones(size(u)[1]))
    nt,ne=size(u)
    p=[u β[1:nt] γ[1:nt] δ[1:nt] uₓ window]
    v₀=[0.0 ,0.0, 0.0, 0.0]
    adj=ODEProblem(corona.sir_adj,v₀,tspan,p)
    adj_sol=solve(adj)
    v=collect(adj_sol(0:tspan[1])')
end
function update(u,v,α,β,γ,δ)
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
    βₙ[1:nt]=β[1:nt]-α*δβ
    γₙ[1:nt]=γ[1:nt]-α*δγ
    δₙ[1:nt]=δ[1:nt]-α*δδ
#    βₙ,γₙ,δₙ
    abs.(βₙ),abs.(γₙ),abs.(δₙ)
end

function linesearch(β,γ,δ,v,uₓ,u₀,tspan,
                    α=1.0e-12,itMax=1024,window=ones(size(u)[1]))
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
    function probe(αₐ,αₑ,nΔ=10)
        Δα=(αₑ-αₐ)
        if (αₐ==αₑ) & (nΔ == 1)
            αₘ=(αₐ+αₑ)/2
            return αₘ, probe(αₘ)
        else
            Δ=αₐ .+ Δα*sort(rand(nΔ+1))
            J=zeros(nΔ+1,2)
            for  i=1:nΔ+1
                J[i,:]=[Δ[i] probe(Δ[i])]
            end
            αₘ=J[argmin(J[:,2]),1]
            Jₘ=minimum(J[:,2])
            return αₘ,Jₘ
        end
    end
    ## Bisection
    i=1
    αₐ= 0.0
    αₑ= α
    ϵ = α/10.0
    J₀=probe(0.0)
    Jₐ=probe(αₐ)
    Jₑ=probe(αₑ)
    while i < itMax
        αₘ=(αₐ+αₑ)/2
#        αₘ=αₐ+ rand()*(αₑ-αₐ)
#        print(αₘ," => ")
        αₘ,Jₘ=probe(αₘ-ϵ/50,αₘ+ϵ/50)
#        println(αₘ)
#        readline()
        if (Jₘ<J₀) & ((αₐ-αₑ)/2<ϵ)
#            println("found $αₘ, $Jₘ")
            return β₁,γ₁,γ₁=update(u,v,αₘ,β,γ,δ)
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
    println("$i failed")
    β,γ,δ
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


end
