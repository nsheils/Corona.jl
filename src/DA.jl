using DifferentialEquations
using LinearAlgebra
using Formatting

###### global constants ###########
const varnames = [:S,:I,:R,:D]
const parnames = [:β,:γ,:δ]

###### structures #################

mutable struct Baseline{T<:Real,D<:TimeType}
    data::TimeArray{T,2,D}
    C::Matrix{T}
    σ::TimeArray{T,1,D}
    u::TimeArray{T,2,D}
    p::TimeArray{T,2,D}

    function Baseline(data::TimeArray{T,2,D},
                      C::Matrix{T},
                      σ::TimeArray{T,1,D},
                      u::TimeArray{T,2,D},
                      p::TimeArray{T,2,D},
                      ) where {T,D}

       if timestamp(data) != timestamp(σ) ||
          timestamp(data) != timestamp(u) ||
          timestamp(data) != timestamp(p)
            error("Input arguments timestamps are difform")
       end
       size_C = (size(data,2),length(varnames))
       if size(C) != size_C
           error("C must have size $size_C")
       end
       if colnames(u) != varnames
           error("Columns of u must have names $varnames")
       end
       if colnames(p) != parnames
           error("Columns of p must have names $parnames")
       end
       new{T,D}(data,C,σ,u,p)
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
        _σ = TimeArray(t,zeros(length(t)),[:σ])
        values(_σ)[lrows] .= 1.0

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

        new{Float64,D}(_data,_C,_σ,_u,_p)
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

    function DA(base::Baseline{T,D}) where {T,D}
        v = backward(base)
        δp = gradient(base,v)
        g = forcing(base)
        J = norm(values(g))
        new{T,D}(base.data,base.C,base.σ,base.u,v,base.p,δp,g,J)
    end

    function DA(da::DA, α::Real)
        DA(da,-α.*da.δp)
    end

    function DA(da::DA, δp::TimeArray{<:Real,2})
        DA(da,values(δp))
    end

    function DA(da::DA, δp::Array{<:Real,2})
        p = da.p .+ δp
        u₀ = OrderedDict(colnames(da.u) .=> values(da.u)[1,:])
        u = forward(p,u₀)
        base = Baseline(da.data,da.C,da.σ,u,p)
        DA(base)
    end
end

function baseline(da::DA)
    Baseline(da.data,da.C,da.σ,da.u,da.p)
end

function extend_solution!(da::DA)
    idxs = meta(da.data)["last_day_idxs"]
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
    σ = values(base.σ)[it]

    p = _interpolation(values(base.p),it,Δt)

    A = dfdu(u,p,t)
    g = base.C' * (base.C * u - data)
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

function forcing(base::Baseline)
    data = values(base.data)
    u = values(base.u)
    σ = values(base.σ)
    C = base.C
    g = (u*C' - data) .* σ
    TimeArray(timestamp(base.data),g,colnames(base.data))
end

function linesearch(da::DA; α=1.0e-12::Float64,
                    maxiters=2048::Int, method="bisection",
                    norm=LinearAlgebra.norm::Function)
    p = values(da.p)
    u = values(da.u)
    v = values(da.v)
    function brute_force(αₐ,αₑ,nΔ=itMax)
        J₀=da.J
        Δα=(αₑ-αₐ)
        Δ=αₐ .+ Δα*sort(rand(nΔ+1))
        J=zeros(nΔ+1,2)
        for  i=1:nΔ+1
            J[i,:]=[Δ[i] DA(da,Δ[i]).J]
        end
        αₘ=J[argmin(J[:,2]),1]
        J₁=minimum(J[:,2])
        return αₘ,true
    end
    function bisection(αₐ=0.0,αₑ=α,itMax=itMax)
        i = 1
        ϵ = 1/3
        J₀= da.J
        # Jₑ=DA(da,αₑ).J
        # if Jₑ<J₀
        #     return αₑ,true
        # end
        Jₐ=J₀
        while i < itMax
            Δα=(αₑ-αₐ)
            αₘ=αₐ+ rand()*Δα
            Jₘ=DA(da,αₘ).J
            if (Jₘ<J₀) && (Δα/α<ϵ)
                return αₘ,true
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
        return 0.0,false
    end
    ###########################
    # calculation start
    if method=="bisection"
        return bisection(0.0,α,maxiters)
    elseif method=="brute_force"
        return brute_force(0.0,α,maxiters)
    elseif method=="plot"
        J=zeros(maxiters,2)
        αrange = range(0,α,length=maxiters)
        for (i,αᵢ) in enumerate(αrange)
            J[i,:] = [αᵢ DA(da,αᵢ).J]
        end
        return J
    else
        throw(ArgumentError("method `$method' unknown"))
    end
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

# function linesearch!(da::DA,
#                      v::AbstractArray{Float64};
#                      α=1.0e-12::Float64,
#                      maxiters=2048::Int,
#                      method="bisection",
#                      norm=LinearAlgebra.norm::Function)
#     p = values(da.p)
#     u = forward(da,p)
#     function probe(α)
#         p₁=update(u,v,α,p)
#         u₁=forward(da,p₁)
#         q₁=u₁*da.datamap'
#         norm((q₁-values(da.data)).*values(da.window))
#     end
#     function brute_force(αₐ,αₑ,nΔ=itMax)
#         J₀=probe(0.0)
#         Δα=(αₑ-αₐ)
#         Δ=αₐ .+ Δα*sort(rand(nΔ+1))
#         J=zeros(nΔ+1,2)
#         for  i=1:nΔ+1
#             J[i,:]=[Δ[i] probe(Δ[i])]
#         end
#         αₘ=J[argmin(J[:,2]),1]
#         J₁=minimum(J[:,2])
#         p₁=update(u,v,αₘ,p)
#         return p₁,true,J₁,αₘ
#     end
#     function bisection(αₐ=0.0,αₑ=α,itMax=itMax)
#         i = 1
#         ϵ = 1/3
#         J₀=probe(0.0)
#         Jₑ=probe(αₑ)
#         if Jₑ<J₀
#             p₁=update(u,v,αₑ,p)
#             return p₁,true,Jₑ,αₑ
#         end
#         Jₐ=J₀
#         while i < itMax
#             #        αₘ=(αₐ+αₑ)/2
#             Δα=(αₑ-αₐ)
#             αₘ=αₐ+ rand()*Δα
#             Jₘ=probe(αₘ)
#             if (Jₘ<J₀) && (Δα/α<ϵ)
#                 p₁=update(u,v,αₘ,p)
#                 return p₁,true,Jₘ,αₘ
#             end
#             if Jₘ < Jₐ
#                 αₐ=αₘ
#                 Jₐ=Jₘ
#             else
#                 αₑ=αₘ
#                 Jₑ=Jₘ
#             end
#             i=i+1
#         end
#         #    println("$i failed")
#         return p,false,J₀,0.0
#     end
#     ###########################
#     # calculation start
#     if method=="bisection"
#         p₁,success,J₁,α₁ = bisection(0.0,α,maxiters)
#     elseif method=="brute_force"
#         p₁,success,J₁,α₁ = brute_force(0.0,α,maxiters)
#     elseif method=="plot"
#         J=zeros(maxiters,2)
#         αrange = range(0,α,length=maxiters)
#         for (i,αᵢ) in enumerate(αrange)
#             J[i,:] = [αᵢ probe(αᵢ)]
#         end
#         return J
#     else
#         throw(ArgumentError("method `$method' unknown"))
#     end
#     values(da.p)[:,:] = p₁
#     success,J₁,α₁
# end
#
#

function heavy_ball(da::DA, α::Float64, β::Float64;
                    ϵ=1e-3::Float64,maxiters=1000::Int,
                    screenfreq=10::Int)
    daᵢ = deepcopy(da)
    J₀ = norm(values(da.data) .* values(da.σ))
    J = Vector{Float64}()
    δp = values(daᵢ.δp)
    if norm(δp) < ϵ
        return da,true
    end
    #αₘ,_ = Corona.linesearch(daᵢ; α=α, maxiters=100, method="brute_force")
    αₘ = α
    Δp = -αₘ.*δp
    for i=1:maxiters
        daᵢ = DA(daᵢ,Δp)
        δp = values(daᵢ.δp)
        if norm(δp)/norm(values(da.p)) < ϵ
            return daᵢ,true
        end
        if daᵢ.J/J₀ < 0.008
            return daᵢ,true
        end
        if mod(i,screenfreq) == 0
            push!(J,daᵢ.J/J₀)
            if length(J) == argmin(J)
                color=:green
            else
                color=:red
            end
            print(i," ")
            printstyled(format("{:.8f}",J[end]),"\n";color=color)
        end
        #αₘ,_ = Corona.linesearch(daᵢ; α=α, maxiters=100, method="brute_force")
        Δp = -αₘ*δp + β*Δp
    end
    return daᵢ,false
end

function diis(da::DA,
               α=1.0e-12::Float64,
               nspace=10::Int,
               method="bisection",
               norm=LinearAlgebra.norm::Function)
    ntParam,nParam=size(da.p)
    ntEqn,nEqn=size(da.data)
    Δp=zeros(ntParam,nParam,nspace)
    e=zeros(ntEqn,nEqn,nspace)
    p₀=da.p
    da′=deepcopy(da)
    for i=1:nspace
        αᵢ,_ = linesearch(da′)
        da′ = DA(da′,αᵢ)
        Δp[:,:,i] = values(da′.p .- p₀)
        e[:,:,i] = values(da′.g)
    end
    B=zeros(nspace+1,nspace+1)
    B[end,1:nspace]=ones(nspace)
    B[1:nspace,end]=ones(nspace)
    x=[zeros(nspace); 1]
    for i=1:nspace
        for j=1:nspace
            B[i,j]=dot(e[:,:,j],e[:,:,i])
        end
    end
    c=B\x
    P=reshape(Δp,ntParam*nParam,nspace)
    δp=reshape(P*c[1:end-1],ntParam,nParam)
    DA(da,TimeArray(timestamp(da.p),δp,colnames(da.p)))
    B,c,δp
end
