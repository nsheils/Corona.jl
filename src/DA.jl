using DifferentialEquations
using LinearAlgebra

###### structures #################

mutable struct Baseline{T<:Real,D<:TimeType}
    data::TimeArray{T,2,D}
    C::Matrix{T}
    σ::TimeArray{T,1,D}
    u::TimeArray{T,2,D}
    p::TimeArray{T,2,D}

    function Baseline(da::DA{T,D}) where {T,D}
        new{T,D}(da.data,da.C,da.σ,da.u,da.p)
    end

    function Baseline(
                data::TimeArray{<:Real,2,D},
                u₀::Vector{<:Real},
                p::Union{Dict{Symbol,<:Real},TimeArray{<:Real,2,D}};
                C=missing::Union{Missing,Matrix{<:Real}},
                start=timestamp(data)[1]::D,
                stop=start+Day(30)::Period,
                step=Day(1)::Period
                ) where {D<:TimeType}

        var_names = [:S,:I,:R,:D]
        param_names = [:β,:γ,:δ]

        if length(init_vals) != length(var_names)
            error("init_vals must have length $(length(var_names))")
        end
        if ismissing(datamap)
            if size(data,2)==length(vars)
                Diagonal(ones(length(vars)))
            else
                error("datamap must be always specified when the number of observed and model parameters differ")
            end
        end
        if size(datamap) != (size(data,2),length(var_names))
            error("datamap must have size $(size(data,2))×$(length(var_names))")
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
        uₓ = TimeArray(t,zeros(length(t),size(data,2)),colnames(data),Dict{String,Any}())
        overlap = intersect(t,timestamp(data))
        lrows = indexin(overlap,t)
        rrows = indexin(overlap,timestamp(data))
        values(uₓ)[lrows,:] = values(data)[rrows,:]
        meta(uₓ)["last_day_idxs"] = lrows[end]

        # Initialize window
        σ = TimeArray(t,zeros(length(t)),[:W])
        values(σ)[lrows] .= 1.0

        # Initialize model parameters
        p = TimeArray(t,zeros(length(t),length(param_names)),param_names)
        if model_params isa Dict
            for (k,v) in enumerate(param_names)
                haskey(model_params,v) || error("model_params does not have key $v")
                values(p)[:,k] .= model_params[v]
            end
        else
            overlap = intersect(t,timestamp(model_params))
            lrows = indexin(overlap,t)
            rrows = indexin(overlap,timestamp(model_params))
            for (k,v) in enumerate(param_names)
                in(v,colnames(model_params)) || error("model_params does not have column $v")
                values(p)[lrows,k] = values(model_params)[rrows,indexin([v],param_names)]
            end
        end

        # Initialize baseline data array
        u = forward(p,u₀)

        new{Float64,D}(uₓ,C,σ,p,u)
    end
end



mutable struct DA{T<:Real,D<:TimeType}
    data::TimeArray{T,2,D}
    C::Matrix{T}
    σ::TimeArray{T,1,D}
    u::TimeArray{T,2,D}
    v::TimeArray{T,2,D}
    p::TimeArray{T,2,D}
    δp::TimeArray{T,2,D}
    g::TimeArray{T,2,D}
    J::T

    function DA(
                data::TimeArray{<:Real,2,D},
                u₀::Vector{<:Real},
                p::Union{Dict{Symbol,<:Real},TimeArray{<:Real,2,D}};
                C=missing::Union{Missing,Matrix{<:Real}},
                start=timestamp(data)[1]::D,
                stop=start+Day(30)::Period,
                step=Day(1)::Period
                ) where {D<:TimeType}

        new{Float64,D}(data,datamap,W,u,p)
    end
end



function DA(base::Baseline)
    v = backward(base)
    δp= ∂J∂α(base,v)
    J = residual(base)
    new{Float64,D}(base.data,base.C,base.σ,base.u,v,base.p,δp,J)
end

function DA(da::DA,α::Float64)
    DA(da,-α*da.δp)
end

function DA(da::DA,δp::Array{Float64,2})
    p=da.p + δp
    u₀=da.u[1,:]
    u=forward(p,u₀)
    base=Baseline(da.data,da.C,da.σ,u,p)
    return da′=DA(base)
end






reconstruct_data(da::DA) =
    TimeArray(timestamp(da.data), values(da.result)*da.datamap', colnames(da.data))

function propagate_solution!(da::DA)
    idxs = meta(da.data)["last_day_idxs"]
    values(da.model_params)[idxs+1:end,:] .= values(da.model_params)[idxs:idxs,:]
end

function add_diff(p::AbstractArray{Float64},it::Int,Δt::Float64)
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
function sir!(du::AbstractArray{Float64,1},
              u::AbstractArray{Float64,1},
              params::AbstractArray{Float64,2},
              t::Float64)
    it = floor(Int,t) + 1
    Δt = t - floor(t)

    β,γ,δ = add_diff(values(params),it,Δt)

    S,I,R,D=u[:]
    N=S+I+R+D
    du[1]=-β*I*S/N
    du[2]= β*I*S/N - (γ+δ)*I
    du[3]= γ*I
    du[4]= δ*I
    du
end

function dfdu(da::DA,t::Float64)
    it=floor(Int,t) + 1
    Δt=t-floor(t)

    #S,I,R,D = add_diff(values(da.result),it,Δt)
    S,I,R,D = values(da.result)[it,:]
    β,γ,δ = add_diff(values(da.model_params),it,Δt)

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

function sir_adj!(v::AbstractArray{Float64,1},
                  da::DA,
                  t::Float64)
    it = floor(Int,t) + 1
    Δt = t - floor(t)

    data = values(da.data)[it,:]
    result = values(da.result)[it,:]
    window = values(da.window)[it,:]

    params = add_diff(values(da.model_params),it,Δt)

    A=dfdu(da,t)
    g = da.datamap' * (da.datamap * result - data)
    dv = -A'*v - g.*window

    return dv
end

function dfda(u::AbstractArray{Float64,2})
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

# function forward(da::DA,
#                  params::AbstractArray{Float64,2};
#                  maxiters=1000000::Int)
#     tspan = (0.0,size(da.data,1)-1.0)
#     trange = range(tspan...,length=size(da.data,1))
#     problem = ODEProblem(sir!,values(da.result)[1,:],tspan,params,maxiters=maxiters)
#     solution = solve(problem)
#     collect(solution(trange)')
# end

function forward(da::DA;maxiters=1000000::Int)
    tspan = (0.0,size(da.data,1)-1.0)
    trange = range(tspan...,length=size(da.data,1))
    u0 = values(da.result)[1,:]
    problem = ODEProblem(sir!,u0,tspan,da,maxiters=maxiters)
    solution = solve(problem)
    collect(solution(trange)')
end

function forward(p,u₀)
    tspan = (0.0,size(p,1)-1.0)
    trange = range(tspan...,length=size(p,1))
    problem = ODEProblem(sir!,u₀,tspan,p,maxiters=maxiters)
    solution = solve(problem)
    collect(solution(trange)')
end


function forward!(da::DA;maxiters=1000000::Int)
    values(da.result) = forward(da,maxiters=maxiters)
end

# function forward!(da::DA;maxiters=1000000::Int)
#     values(da.result)[:,:] = forward(da,values(da.model_params),maxiters=maxiters)
# end

function backward(da::DA;maxiters=1000000::Int)
    tspan = (size(da.result,1)-1.0,0.0)
    trange = range(tspan...,length=size(da.result,1))
    v0 = values(da.adjoint)[1,:]
    adj = ODEProblem(sir_adj!,v0,tspan,da,maxiters=maxiters)
    adj_sol = solve(adj)
    collect(adj_sol(reverse(trange))')
end

function backward(;maxiters=1000000::Int)
    tspan = (size(u,1)-1.0,0.0)
    trange = range(tspan...,length=size(da.result,1))
    v₀ = zeros(size(u,2))
    adj = ODEProblem(sir_adj!,v₀,tspan,base,maxiters=maxiters)
    adj_sol = solve(adj)
    collect(adj_sol(reverse(trange))')
end

function backward!(da::DA;maxiters=1000000::Int)
    values(da.adjoint) = backward(da,maxiters=maxiters)
end

function ∂J∂α(base,V)
    u=values(base.u)
    v=values(V)
    p=values(base.p)

    f=dfda(u)
    δp=zeros(size(p))
    for i=1:size(u,2)
        δp = δp + v[:,i] .* f[:,i,:]
    end
    TimeArray(timestamp(base.p),δp,colnames(base.p))
end

function residual(base)
    data = values(base.data)
    u = values(base.u)
    σ = values(base.σ)
    C = base.C
    J = norm( (u*C' - data) .* σ )
end

import Base.+
function +(da₀::DA,δp)
    da = deepcopy(da₀)
    da.p = da.p .+ δp
    forward!(da)
    backward!(da)
    ∂J∂α!(da)
    residual!(da)
    da
end

function update(da₀::DA,α::Float64)
    δp = values(da.δp)
    da = da₀ + (-α*δp)
end

function update(u,v,α,p)
    if α > 0
        f=dfda(u)
        δp=zeros(size(p))
        for i=1:4
            δp = δp + v[:,i] .* f[:,i,:]
        end
        max.(p-α*δp,0.0)
    else
        p
    end
end

function linesearch!(da::DA,
                     v::AbstractArray{Float64};
                     α=1.0e-12::Float64,
                     maxiters=2048::Int,
                     method="bisection",
                     norm=LinearAlgebra.norm::Function)
    p = values(da.model_params)
    u = forward(da,p)
    function probe(α)
        p₁=update(u,v,α,p)
        u₁=forward(da,p₁)
        q₁=u₁*da.datamap'
        norm((q₁-values(da.data)).*values(da.window))
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
        p₁=update(u,v,αₘ,p)
        return p₁,true,J₁,αₘ
    end
    function bisection(αₐ=0.0,αₑ=α,itMax=itMax)
        i = 1
        ϵ = 1/3
        J₀=probe(0.0)
        Jₑ=probe(αₑ)
        if Jₑ<J₀
            p₁=update(u,v,αₑ,p)
            return p₁,true,Jₑ,αₑ
        end
        Jₐ=J₀
        while i < itMax
            #        αₘ=(αₐ+αₑ)/2
            Δα=(αₑ-αₐ)
            αₘ=αₐ+ rand()*Δα
            Jₘ=probe(αₘ)
            if (Jₘ<J₀) && (Δα/α<ϵ)
                p₁=update(u,v,αₘ,p)
                return p₁,true,Jₘ,αₘ
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
        return p,false,J₀,0.0
    end
    ###########################
    # calculation start
    if method=="bisection"
        p₁,success,J₁,α₁ = bisection(0.0,α,maxiters)
    elseif method=="brute_force"
        p₁,success,J₁,α₁ = brute_force(0.0,α,maxiters)
    elseif method=="plot"
        J=zeros(maxiters,2)
        αrange = range(0,α,length=maxiters)
        for (i,αᵢ) in enumerate(αrange)
            J[i,:] = [αᵢ probe(αᵢ)]
        end
        return J
    else
        throw(ArgumentError("method `$method' unknown"))
    end
    values(da.model_params)[:,:] = p₁
    success,J₁,α₁
end



function diis(da::DA,
               α=1.0e-12::Float64,
               maxiters=10::Int,
               method="bisection",
               norm=LinearAlgebra.norm::Function)
    p₀ = values(da.model_params)
    ntParam,nParam=size(p)
    ntEqn,nEqn=size(values(da.data))
    δp=zeros(ntParam,nParam,maxiters)
    da=DA(da,δp[:,:,i])
    q=u₁*da.datamap'
    e=zeros(ntEqn,  nEqn,  maxiters)
    e[:,:,1]=(q-values(da.data)).*values(da.window)
    for i =2:maxiters
        u = forward(da,p₀)
        v = backward(da);
        δp[:,:,i]=gradient(u,v,δp[:,:,i-1])
        da.model_params=p₀-α*δp[:,:,i]
        u₁=forward(da,p₀-α*δp[:,:,i])
        q₁=u₁*da.datamap'
        e[:,:,i]=da.g
    end
    B=zeros(maxiters+1,maxiters+1)
    B[end,1:maxiters]=ones(maxiters)
    B[1:maxiters,end]=ones(maxiters)
    x=[zeros(maxiters); 1]
    for i=1:maxiters
        for j=1:maxiters
            B[i,j]=dot(e[:,:,j],e[:,:,i])
        end
    end
    c=B\x
    P=reshape(δp,ntParam*nParam,maxiters)
    δp′=reshape(P*c,ntParam,nParam)
    da.model_params=p₀-α*δp′
    return da
end
