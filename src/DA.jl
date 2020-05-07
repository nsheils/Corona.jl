using DifferentialEquations
using LinearAlgebra
using Flux
using Interpolations
# using SparseArrays

import Flux.Optimise: apply!

###### global constants ###########
const varnames = [:S,:I,:R,:D]
const parnames = [:β,:γ,:δ]

###### structures #################

mutable struct Baseline{T<:Real,D<:TimeType}
    time::StepRange{D}
    data::TimeArray{T,2,D}
    C::Matrix{T}
    σ::TimeArray{T,2,D}
    μ::TimeArray{T,2,D}
    u::TimeArray{T,2,D}
    g::TimeArray{T,2,D}
    p::TimeArray{T,2,D}

    function Baseline(time::StepRange{D},
                      data::TimeArray{T,2,D},
                      C::Matrix{T},
                      σ::TimeArray{T,2,D},
                      μ::TimeArray{T,2,D},
                      u::TimeArray{T,2,D},
                      g::TimeArray{T,2,D},
                      p::TimeArray{T,2,D},
                      ) where {T,D}

       @assert timestamp(data) == timestamp(σ) == timestamp(μ) ==
            timestamp(u) == timestamp(g) == timestamp(p) == time
            "input arguments timestamps are difform"

       size_C = (size(data,2),length(varnames))
       @assert size(C) == size_C "C must have size $size_C"

       @assert size(σ,2) == size(data,2) "σ must have $(size(data,2)) columns"
       @assert colnames(σ) == colnames(data) "columns of σ must have names $(colnames(data))"

       @assert size(μ,2) == length(parnames) "μ must have $(length(parnames)) columns"
       @assert colnames(μ) == parnames "columns of μ must have names $parnames"

       @assert size(u,2) == length(varnames) "u must have $(length(varnames)) columns"
       @assert colnames(u) == varnames "columns of u must have names $varnames"

       @assert size(p,2) == length(parnames) "u must have $(length(parnames)) columns"
       @assert colnames(p) == parnames "columns of p must have names $parnames"

       new{T,D}(time, data, C, σ, μ, u, g, p)
    end

    function Baseline(
                data::TimeArray{<:Real,2,<:TimeType},
                u₀::Dict{Symbol,<:Real},
                p::Union{Dict{Symbol,<:Real},TimeArray{<:Real,2,<:TimeType}};
                C=missing::Union{Missing,Matrix{<:Real}},
                start=timestamp(data)[1]::TimeType,
                stop=start+Day(30)::Period,
                step=Day(1)::Period,
                σ=missing::Union{Missing,TimeArray{<:Real,2,<:TimeType}},
                μ=missing::Union{Missing,TimeArray{<:Real,2,<:TimeType}},
                # dataw=step::Period,
                interptype=BSpline(Linear())::Union{Interpolations.InterpolationType,Tuple{Vararg{Interpolations.InterpolationType}}}
                )

        @assert length(u₀) == length(varnames) "u₀ must have length $(length(varnames))"

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

        @assert start >= timestamp(data)[1] "start date must be after the last date available in data"
        @assert start < timestamp(data)[end] "start date must be before the last date available in data"
        @assert stop > start "stop date must be after the start date"
        @assert step > Day(0) "step must be a positive time period"

        @assert size(σ,2) == size(data,2) "σ must have $(size(data,2)) columns"
        @assert colnames(σ) == colnames(data) "columns of σ must have names $(colnames(data))"

        @assert size(μ,2) == length(parnames) "μ must have $(length(parnames)) columns"
        @assert colnames(μ) == parnames "columns of μ must have names $parnames"

        if interptype isa Tuple
            @assert length(interptype) == size(data,2) "interptype must be a tuple with $(size(data,2)) elements"
        else
            interptype = ntuple(i -> interptype, size(data,2))
        end

        # time
        if step < Day(1)
            t = DateTime(start):step:DateTime(stop)
        else
            t = start:step:stop
        end

        # Copy data to working data array ux
        _data = TimeArray(t,zeros(length(t),size(data,2)),colnames(data),Dict{String,Any}())
        # for j=1:size(data,2)
        #     values(_data)[:,j] = _interpolation(timerange(timestamp(data)),
        #                             values(data)[:,j], t)
        # end
        for j=1:size(data,2)
            values(_data)[:,j] = _interpolation(timerange(timestamp(data)),
                                    values(data)[:,j], t + Millisecond(step)/2,
                                    #BSpline(Constant()), Flat())
                                    interptype[j], Flat())
        end
        data_t = convert.(eltype(t), timestamp(data))
        # _data_idxs, data_idxs = Corona.overlap(t, data_t)
        # values(_data)[_data_idxs,:] = values(data)[data_idxs,:]
        meta(_data)["lastdate"] = convert(eltype(t), timestamp(data)[end])

        # expfun = exp(-32.0*(Dates.toms(t .- data_t[i])./Dates.toms(datawidth)).^2)
        # values(_data) .+=  values(data)[i,:].*expfun
        # values(_σ)

        # Initialize forcing window
        _σ = TimeArray(t,zeros(size(_data)),colnames(data))
        if ismissing(σ)
            σ = TimeArray(timestamp(data),ones(size(data,2)),varnames)
        end
        for j=1:size(σ,2)
            values(_σ)[:,j] = _interpolation(timerange(timestamp(σ)),
                                values(σ)[:,j], t, BSpline(Linear()), 0.0)
        end
        #σ_t = convert.(eltype(t), timestamp(σ))
        # _σ_idxs, σ_idxs, _ = Corona.overlap(t, σ_t, data_t)
        # values(_σ)[_σ_idxs,:] = values(σ)[σ_idxs,:]
        # σ_data = similar(values(data))
        # for j=1:size(data,2)
        #     σ_data[:,j] = _interpolation(timerange(timestamp(σ)),
        #                             values(σ)[:,j], data_t, BSpline(Linear()), 0.0)
        # end

        # not_σ = ones(size(_data))
        # for i=1:size(data,1)
        #     not_σ .*= 1.0 .- exp.(-32.0*(Dates.toms.(t .- data_t[i])./Dates.toms.(dataw)).^2)
        # end
        # values(_σ) .*= 1.0 .- not_σ

        # Initialize control window
        _μ = TimeArray(t,zeros(length(t),length(parnames)),parnames)
        if ismissing(μ)
            μ = TimeArray(timestamp(data),ones(size(data,1),length(parnames)),parnames)
        end
        # for j=1:size(μ,2)
        #     values(_μ)[:,j] = _interpolation(timerange(timestamp(μ)),
        #                         values(μ)[:,j], t, smooth = 1.0, extrapscheme = 0.0)
        # end
        for j=1:size(μ,2)
            values(_μ)[:,j] = _interpolation(timerange(timestamp(μ)),
                                values(μ)[:,j], t, BSpline(Linear()), 0.0)
        end

        # Initialize model parameters
        _p = TimeArray(t,zeros(length(t),length(parnames)),parnames)
        if p isa Dict
            for (k,v) in enumerate(parnames)
                haskey(p,v) || error("p does not have key $v")
                values(_p)[:,k] .= p[v]
            end
        else
            # for j=1:size(p,2)
            #     values(_p)[:,j] = _interpolation(timerange(timestamp(p)),
            #                         values(p)[:,j], t, smooth = 1.0)
            # end
            for j=1:size(p,2)
                values(_p)[:,j] = _interpolation(timerange(timestamp(p)),
                                    values(p)[:,j], t, BSpline(Linear()), Flat())
            end
        end

        # Initialize baseline data array
        _u₀ = OrderedDict{Symbol,Float64}()
        for v in varnames
            haskey(u₀,v) || error("u₀ does not have key $v")
            _u₀[v] = u₀[v]
        end
        _u = forward(t, values(_p), _u₀)

        # Initialize forcing
        _g = forcing(_data,_C,_σ,_u)

        new{Float64,eltype(t)}(t, _data, _C, _σ, _μ, _u, _g, _p)
    end
end

struct BaselineInterpolation{T<:Real}
    data::AbstractInterpolation{T,2}
    u::AbstractInterpolation{T,2}
    g::AbstractInterpolation{T,2}
    p::AbstractInterpolation{T,2}

    function BaselineInterpolation(base::Baseline{T}, trange::AbstractRange{<:Real}) where {T}
        interptype = (BSpline(Linear()), NoInterp())
        data = scale(interpolate(values(base.data), interptype), trange, 1:size(base.data,2))
        u = scale(interpolate(values(base.u), interptype), trange, 1:size(base.u,2))
        g = scale(interpolate(values(base.g), interptype), trange, 1:size(base.g,2))
        p = scale(interpolate(values(base.p), interptype), trange, 1:size(base.p,2))
        new{T}(data, u, g, p)
    end
end

mutable struct DA{T<:Real,D<:TimeType}
    time::StepRange{D}
    data::TimeArray{T,2,D}
    C::Matrix{T}
    σ::TimeArray{T,2,D}
    μ::TimeArray{T,2,D}
    u::TimeArray{T,2,D}
    g::TimeArray{T,2,D}
    p::TimeArray{T,2,D}
    v::TimeArray{T,2,D}
    δp::TimeArray{T,2,D}

    function DA(base::Baseline{T,D}) where {T,D}
        v = backward(base)
        δp = gradient(base,v)
        new{T,D}(base.time, base.data, base.C, base.σ, base.μ,
                    base.u, base.g,base.p, v, δp)
    end

    function DA(da::DA, p::TimeArray{<:Real,2})
        u₀ = OrderedDict(colnames(da.u) .=> values(da.u)[1,:])
        u = forward(da.time, values(p), u₀)
        g = forcing(da.data, da.C ,da.σ ,u)
        base = Baseline(da.time, da.data, da.C, da.σ, da.μ, u, g, p)
        DA(base)
    end
end

function baseline(da::DA)
    Baseline(da.time, da.data, da.C, da.σ, da.μ, da.u, da.g, da.p)
end

function update(da::DA, α::Real)
    update(da, α.*da.δp)
end

function update(da::DA, δp::TimeArray{<:Real,2})
    update(da, values(δp))
end

function update(da::DA, δp::Array{<:Real,2})
    p = max.(da.p .- δp, 0.0)
    DA(da, p)
end

function residual(base::Baseline; relative=false::Bool)
    J = 0.5 * norm(values(base.g), 2)
    # d = max.(values(base.data), 1)
    # J = 0.5 * norm((values(base.u) * base.C' - values(base.data))
    #                 .* values(base.σ) ./ d * base.C, 2)
    if relative
        J = J/datanorm(base)
    end
    J
end

function datanorm(base::Baseline)
    0.5 * norm(values(base.data) .*values(base.σ) * base.C, 2)
end

function extend_solution!(da::DA)
    lastdate = meta(da.data)["lastdate"] - Day(1)
    idxs = findfirst(lastdate .== timestamp(da.data))
    values(da.p)[idxs+1:end,:] .= values(da.p)[idxs:idxs,:]
end

function timerange(d::AbstractArray{<:TimeType})
    drange = range(d[1], d[end], step = d[2] - d[1])
    @assert d == drange "d must be equispaced in time"
    drange
end

function timestep(trange::StepRange{<:TimeType,<:Period})
    Dates.toms(step(trange))/Dates.toms(Day(1))
end

import Interpolations.interpolate
interpolate(range::AbstractRange, vs::AbstractVector, interptype::BSpline) = scale(interpolate(vs, interptype), range)

function _interpolation(t::AbstractRange{<:TimeType}, u::AbstractVector{<:Real},
                    tq::AbstractVector{<:TimeType}, interptype, extrapscheme)
# function _interpolation(t::AbstractRange{<:TimeType}, u::AbstractVector{<:Real},
#                     tq::AbstractVector{<:TimeType}; smooth=nothing::Union{Nothing,<:Real},
#                     extrapscheme=:fill::Union{Symbol,<:Real})
        @assert step(t) > Millisecond(0) "step of t must be positive"
        x = datetime2unix(DateTime(t.start)):Second(t.step).value:datetime2unix(DateTime(t.stop))
        xq = datetime2unix.(DateTime.(tq))
        #itp = scale(interpolate(u, interptype), x)
        itp = interpolate(x, u, interptype)
        etp = extrapolate(itp, extrapscheme)
        etp.(xq)
        # spl = fit(SmoothingSpline, float(x), float(u), 2000.0)
        # SmoothingSplines.predict(spl, float(xq))
        # spl = csaps.CubicSmoothingSpline(x, u, smooth = smooth)
        # lrng = xq .< x[1]
        # urng = xq .> x[end]
        # irng = .!lrng .& .!urng
        # uq = similar(xq, Float64)
        # uq[irng] = spl(xq[irng])
        # if extrapscheme isa Symbol
        #     if extrapscheme == :fill
        #         uq[lrng] .= spl([x[1]])
        #         uq[urng] .= spl([x[end]])
        #     else
        #         error("unknown extrapscheme $(extrapscheme)")
        #     end
        # else
        #     uq[lrng] .= extrapscheme
        #     uq[urng] .= extrapscheme
        # end
        # uq
end

"Classic Epidemic Model  Hethcote (2000), added deaths"
function sir!(du::Vector{Float64}, u::Vector{Float64},
              p::AbstractInterpolation{Float64,2}, t::Float64)
    β,γ,δ = p(t,1:3)
    S,I,R,D = u[:]
    N = sum(u)
    du[1]=-β*I*S/N
    du[2]= β*I*S/N - (γ+δ)*I
    du[3]= γ*I
    du[4]= δ*I
    return du
end

function dfdu(u::Vector{Float64}, p::Vector{Float64})
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

function sir_adj(v::Vector{Float64}, base::BaselineInterpolation, t::Float64)
    u = base.u(t,1:4)
    p = base.p(t,1:3)
    g = base.g(t,1:4)
    A = dfdu(u,p)
    dv = -A'*v - g
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

function forward(t::StepRange{<:TimeType}, p::Array{Float64,2}, u₀::OrderedDict{Symbol,Float64})
    Δt = timestep(t)
    tspan = (0.0, (size(p,1) - 1)*Δt)
    trange = range(tspan..., length=size(p,1))
    _u₀ = collect(values(u₀))
    interptype = (BSpline(Linear()), NoInterp())
    itp = scale(interpolate(p, interptype), trange, 1:size(p,2))
    problem = ODEProblem(sir!, _u₀, tspan, itp)
    solution = solve(problem, RK4(), dt=1.0, adaptive=false)
    u = collect(solution(trange)')
    TimeArray(t, u, collect(keys(u₀)))
end

function backward(base::Baseline)
    Δt = timestep(base.time)
    tspan = ((size(base.u,1) - 1)*Δt, 0.0)
    trange = reverse(range(tspan...,length=size(base.u,1)))
    v₀ = zeros(size(base.u,2))
    base_itp = BaselineInterpolation(base, trange)
    adj = ODEProblem(sir_adj, v₀, tspan, base_itp)
    adj_sol = solve(adj, tstops=trange)
    adj_sol = solve(adj, RK4(), dt=1.0, adaptive=false)
    v = collect(adj_sol(trange)')
    TimeArray(base.time, v, colnames(base.u))
end

function gradient(base::Baseline, v::TimeArray)
# function gradient(base::Baseline, v::TimeArray; λ=0.001::Real)
    u = values(base.u)
    p = values(base.p)
    f = dfda(u)
    δp = zeros(size(p))
    for i=1:size(u,2)
        δp = δp + values(v)[:,i] .* f[:,i,:]
    end
    δp ./= maximum(abs.(δp))
    # N = size(p,1)
    # D = spdiagm(N, N, 0 => -1*ones(N - 1), 1 => ones(N - 1))
    # δp = δp/maximum(abs.(δp)) - D'*sign.(D*p)*λ
    TimeArray(timestamp(base.p), δp, colnames(base.p))
end

function forcing(data::TimeArray{Float64,2}, C::Matrix{Float64},
                 σ::TimeArray{Float64,2}, u::TimeArray{Float64,2})
    g = (values(u) * C' - values(data)) .* values(σ) * C
    # λ = 100000.0
    # DD(N) = spdiagm(N, N, -1 => -ones(N - 2), 0 => 2*ones(N - 1), 1 => -ones(N - 2))
    # g = (values(u) * C' - values(data)) .* values(σ) * C
            # + λ * DD(size(u,1)) * values(u)
    # d = max.(values(data), 1)
    # g = (values(u) * C' - values(data)) .* (values(σ) ./ d).^2 * C
    TimeArray(timestamp(u), g, colnames(u))
end

function apply!(opt, da::DA)
    δp = values(da.δp) .* values(da.μ)
    Δp = apply!(opt, values(da.p), δp)
    update(da, Δp)
end
