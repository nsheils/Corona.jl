module corona

export extract
using CSV
using DataFrames
using Dates
using TimeSeries
using DifferentialEquations
using LinearAlgebra
function read(dataset::String,data="/home/jls/data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/")
    if dataset=="Confirmed"
        Confirmed=CSV.read(data*"time_series_19-covid-Confirmed.csv");
    elseif dataset=="Recovered"
        Recovered=CSV.read(data*"time_series_19-covid-Recovered.csv")
    elseif dataset=="Deaths"
        Deaths=CSV.read(data*"time_series_19-covid-Deaths.csv")
    else
        throw(ArgumentError("no datset: $dataset"))
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



"Classic Epedemic Model  Hethcote (2000)"
function sir!(du,u,p,t)
    N=sum(u[1:3])
    it=Int(round(t)+1)
    β,γ=p[it,1:2]
    # contact rate, average infectous period
    S,I,R=u[1:3]
    du[1]=-β*I*S/N
    du[2]= β*I*S/N - γ*I
    du[3]= γ*I
end
function sir_adj(u,p,t)
    it=Int(round(t)+1)
    N=sum(p[it,1:3])
    S,I,R=p[it,1:3]
    g=p[it,4:6]
    β,γ=p[it,7:8]
    A=zeros(3,3)
    A[1,1]=-β*I/N
    A[1,2]=-β*S/N
    A[1,3]=+β*I*S/N^2

    A[2,1]= β*I/N
    A[2,2]= β*S/N - γ
    A[2,3]=-β*I*S/N^2

    A[3,2]=γ
    du=-A'*u - g
end

function force(u)
    N=sum(u,dims=2)
    S=u[:,1]
    I=u[:,2]
    R=u[:,3]
    f=zeros(length(S),3,2)
    f[:,1,1]=-I.*S./N
    f[:,2,1]= I.*S./N
    f[:,2,2]=-I
    f[:,3,2]= I
    return f
end


function forward(β,γ,u₀,tspan=(0.0,50.0))
    p=[β γ]
    problem=ODEProblem(sir!,u₀,tspan,p)
    solution=solve(problem)
    u=collect(solution(0:tspan[2])')
end
function backward(u,uₓ,β,γ,tspan=(50.0,0.0))
    nt,ne=size(u)
    g=u .- uₓ[1:nt,:]
    p=[u g β[1:nt] γ[1:nt]]
    v₀=[0.0 ,0.0, 0.0]
    adj=ODEProblem(corona.sir_adj,v₀,tspan,p)
    adj_sol=solve(adj)
    as=collect(adj_sol(0:tspan[1])')
end
function update(u,as,α,β,γ)
    nt,ne=size(u)
    f=force(u)
    δβ=f[:,1,1].*as[:,1]+f[:,2,1].*as[:,2]+f[:,3,1].*as[:,3]
    δγ=f[:,1,2].*as[:,1]+f[:,2,2].*as[:,2]+f[:,3,2].*as[:,3]
    βₙ=copy(β)
    γₙ=copy(γ)
    βₙ[1:nt]=β[1:nt]-α*δβ
    γₙ[1:nt]=γ[1:nt]+α*δγ
    βₙ[nt+1:end] .=βₙ[nt]
    γₙ[nt+1:end] .=γₙ[nt]
    βₙ,γₙ
end

function linesearch(β,γ,a,uₓ,u₀,tspan,α=1.0e-12)
    J=zeros(500,2)
    u=corona.forward(β,γ,u₀)
    actual=norm(uₓ-u)
    for i=1:500
        β₁,γ₁=corona.update(u,a,i*α,β,γ)
        u₁=corona.forward(β₁,γ₁,u₀,tspan)
        global J[i,:]=[i*α norm(uₓ-u₁)]
    end
    α₁=argmin(J[:,2])*α
    β₁,γ₁=corona.update(u,a,α₁,β,γ)
end

end
