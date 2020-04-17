using DifferentialEquations
using LinearAlgebra
using Flux
using Formatting

import Flux.Optimise: apply!

###### global constants ###########
const varnames = [:S,:I,:R,:D]
const parnames = [:β,:γ,:δ]

###### structures #################

mutable struct Baseline{T<:Real,D<:TimeType}
    data::TimeArray{T,2,D}
    C::Matrix{T}
    σ::TimeArray{T,2,D}
    u::TimeArray{T,2,D}
    p::TimeArray{T,2,D}
    g::TimeArray{T,2,D}

    function Baseline(data::TimeArray{T,2,D},
                      C::Matrix{T},
                      σ::TimeArray{T,2,D},
                      u::TimeArray{T,2,D},
                      p::TimeArray{T,2,D},
                      g::TimeArray{T,2,D},
                      ) where {T,D}

       if timestamp(data) != timestamp(σ) ||
          timestamp(data) != timestamp(u) ||
          timestamp(data) != timestamp(p)
            error("input arguments timestamps are difform")
       end
       size_C = (size(data,2),length(varnames))
       if size(C) != size_C
           error("C must have size $size_C")
       end
       if size(σ,2) != length(varnames)
           error("σ must have $(length(varnames)) columns")
       end
       if colnames(σ) != varnames
           error("columns of σ must have names $varnames")
       end
       if size(u,2) != length(varnames)
           error("u must have $(length(varnames)) columns")
       end
       if colnames(u) != varnames
           error("columns of u must have names $varnames")
       end
       if size(p,2) != length(parnames)
           error("u must have $(length(parnames)) columns")
       end
       if colnames(p) != parnames
           error("columns of p must have names $parnames")
       end
       if size(g,2) != length(varnames)
           error("u must have $(length(parnames)) columns")
       end
       if colnames(g) != varnames
           error("columns of g must have names $varnames")
       end
       new{T,D}(data,C,σ,u,p,g)
    end

    function Baseline(
                data::TimeArray{<:Real,2,D},
                u₀::Dict{Symbol,<:Real},
                p::Union{Dict{Symbol,<:Real},TimeArray{<:Real,2,D}};
                C=missing::Union{Missing,Matrix{<:Real}},
                start=timestamp(data)[1]::D,
                stop=start+Day(30)::Period,
                step=Day(1)::Period,
                ) where {D<:TimeType}

        if length(u₀) != length(varnames)
            error("u₀ must have length $(length(varnames))")
        end
        if ismissing(C)
            if size(data,2) == length(varnames)
                _C = Diagonal(ones(length(varnames)))
            else
                error("C must be always specified when the number of observed and model parameters differ")
            end
        else
            if size(C) != (size(data,2),length(varnames))
                error("datamap must have size $(size(data,2))×$(length(varnames))")
            else
                _C = float(C)
            end
        end
        if start < timestamp(data)[1]
            error("start date must be after the last date available in data")
        end
        if start > timestamp(data)[end]
            error("start date must be before the last date available in data")
        end
        if stop < start
            error("stop date must be after the start date")
        end
        if step <= Day(0)
            error("step must be a positive time period")
        end

        # time
        t = start:step:stop

        # Copy data to working data array ux
        _data = TimeArray(t,zeros(length(t),size(data,2)),colnames(data),Dict{String,Any}())
        overlap = intersect(t,timestamp(data))
        lrows = indexin(overlap,t)
        rrows = indexin(overlap,timestamp(data))
        values(_data)[lrows,:] = values(data)[rrows,:]
        meta(_data)["last_day_idxs"] = lrows[end]

        # Initialize window
        _σ = TimeArray(t,zeros(length(t),length(varnames)),varnames)
        values(_σ)[lrows,:] .= 1.0

        # Initialize model parameters
        _p = TimeArray(t,zeros(length(t),length(parnames)),parnames)
        if p isa Dict
            for (k,v) in enumerate(parnames)
                haskey(p,v) || error("p does not have key $v")
                values(_p)[:,k] .= p[v]
            end
        else
            overlap = intersect(t,timestamp(p))
            lrows = indexin(overlap,t)
            rrows = indexin(overlap,timestamp(p))
            for (k,v) in enumerate(parnames)
                in(v,colnames(p)) || error("p does not have column $v")
                values(_p)[lrows,k] = values(p)[rrows,indexin([v],parnames)]
            end
        end

        # Initialize baseline data array
        _u₀ = OrderedDict{Symbol,Float64}()
        for v in varnames
            haskey(u₀,v) || error("u₀ does not have key $v")
            _u₀[v] = u₀[v]
        end
        _u = forward(_p,_u₀)

        new{Float64,D}(_data,_C,_σ,_u,_p,_g)
    end
end

mutable struct DA{T<:Real,D<:TimeType}
    data::TimeArray{T,2,D}
    C::Matrix{T}
    σ::TimeArray{T,2,D}
    u::TimeArray{T,2,D}
    v::TimeArray{T,2,D}
    p::TimeArray{T,2,D}
    δp::TimeArray{T,2,D}
    g::TimeArray{T,2,D}

    function DA(base::Baseline{T,D}) where {T,D}
        v = backward(base)
        δp = gradient(base,v)
        g = forcing(base)
        new{T,D}(base.data,base.C,base.σ,base.u,v,base.p,δp,g)
    end

    function DA(da::DA, α::Real)
        DA(da,α.*da.δp)
    end

    function DA(da::DA, δp::TimeArray{<:Real,2})
        DA(da,values(δp))
    end

    function DA(da::DA, δp::Array{<:Real,2})
        p = max.(da.p .- δp, 0.0)
        u₀ = OrderedDict(colnames(da.u) .=> values(da.u)[1,:])
        u = forward(p,u₀)
        base = Baseline(da.data,da.C,da.σ,u,p)
        DA(base)
    end
end

function baseline(da::DA)
    Baseline(da.data,da.C,da.σ,da.u,da.p)
end

function residual(da::DA; relative=false::Bool, norm=LinearAlgebra.norm::Function)
    J = norm(values(da.g))
    if relative
        J = J/datanorm(da, norm=norm)
    end
    J
end

function datanorm(da::DA; norm=LinearAlgebra.norm::Function)
    norm(values(da.data).*values(da.σ))
end

function extend_solution!(da::DA)
    idxs = meta(da.data)["last_day_idxs"]-1
    values(da.p)[idxs+1:end,:] .= values(da.p)[idxs:idxs,:]
end

function _interpolation(p::AbstractArray{Float64},it::Int,Δt::Float64)
    # if it == size(p,1)
    #     p[it,:] + Δt*(p[it,:] -p[it-1,:])
    # elseif it == 1
    #     p[it,:] + Δt*(p[it+1,:] -p[it,:])
    # else
    #     p[it,:] + Δt*(p[it+1,:] -p[it-1,:])/2
    # end
    if it < size(p,1)
        p=(p[it,:] + p[it+1,:])/2
    else
        p[it,:]
    end
end

"Classic Epidemic Model  Hethcote (2000), added deaths"
function sir!(du::Vector{Float64}, u::Vector{Float64},
              p::TimeArray{Float64,2}, t::Float64)
    it = floor(Int,t) + 1
    Δt = t - floor(t)

    β,γ,δ = _interpolation(values(p),it,Δt)

    S,I,R,D=u[:]
    N=sum(u)
    du[1]=-β*I*S/N
    du[2]= β*I*S/N - (γ+δ)*I
    du[3]= γ*I
    du[4]= δ*I
    du
end

function dfdu(u::Vector{Float64}, p::Vector{Float64}, t::Float64)
    S,I,R,D = u[:]
    β,γ,δ = p[:]

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

function sir_adj(v::Vector{Float64}, base::Baseline, t::Float64)
    it = floor(Int,t) + 1
    Δt = t - floor(t)

    data = values(base.data)[it,:]
    u = values(base.u)[it,:]
    σ = values(base.σ)[it,:]
    g = values(base.g)[it,:]

    p = _interpolation(values(base.p),it,Δt)

    A = dfdu(u,p,t)
    dv = -A'*v - g.*σ

    return dv
end

function dfda(u::Matrix{Float64})
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

function forward(p::TimeArray{Float64,2}, u₀::OrderedDict{Symbol,Float64};
                 maxiters=1000000::Int)
    tspan = (0.0,size(p,1)-1.0)
    trange = range(tspan...,length=size(p,1))
    _u₀ = collect(values(u₀))
    problem = ODEProblem(sir!,_u₀,tspan,p,maxiters=maxiters)
    solution = solve(problem)
    u = collect(solution(trange)')
    TimeArray(timestamp(p),u,collect(keys(u₀)))
end

function backward(base::Baseline; maxiters=1000000::Int)
    tspan = (size(base.u,1)-1.0,0.0)
    trange = range(tspan...,length=size(base.u,1))
    v₀ = zeros(size(base.u,2))
    adj = ODEProblem(sir_adj,v₀,tspan,base,maxiters=maxiters)
    adj_sol = solve(adj)
    v = collect(adj_sol(reverse(trange))')
    TimeArray(timestamp(base.u),v,colnames(base.u))
end

function gradient(base::Baseline, v::TimeArray)
    u = values(base.u)
    p = values(base.p)
    f = dfda(u)
    δp = zeros(size(p))
    for i=1:size(u,2)
        δp = δp + values(v)[:,i] .* f[:,i,:]
    end
    TimeArray(timestamp(base.p),δp,colnames(base.p))
end

function forcing(data::TimeArray{Float64,2}, u::TimeArray{Float64,2},
                 C::Matrix{Float64})
    data = values(base.data)
    u = values(base.u)
    C = base.C
    g = (values(u) * C' - values(data)) * C
    TimeArray(timestamp(base.u),g,colnames(base.u))
end

function apply!(opt, da::DA)
    Δp = apply!(opt,values(da.p),values(da.δp))
    DA(da,Δp)
end
