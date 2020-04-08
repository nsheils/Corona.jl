using Base: @propagate_inbounds

import Base: getindex, setindex!

###### getindex #################
@propagate_inbounds getindex(ta::TimeArray{T,N,D},d::D,cols::AbstractVector{Symbol}) where {T,N,D<:TimeType} =
    getindex(ta,[d],cols)

@propagate_inbounds getindex(ta::TimeArray{T,N,D},dates::AbstractVector{D},c::Symbol) where {T,N,D<:TimeType} =
    getindex(ta,dates,[c])

@propagate_inbounds getindex(ta::TimeArray{T,N,D},dates::StepRange{D},cols::AbstractVector{Symbol}) where {T,N,D<:TimeType} =
    getindex(ta,collect(dates),cols)

@propagate_inbounds function getindex(ta::TimeArray{T,N,D},dates::Vector{D},cols::AbstractVector{Symbol}) where {T,N,D<:TimeType}
    rows, _ = TimeSeries.overlap(timestamp(ta), dates)
    getindex(ta,rows,cols)
end

###### setindex #################

setindex!(ta::TimeArray) = throw(BoundsError(typeof(ta), []))

@propagate_inbounds function setindex!(ta::TimeArray{T,1,D},v::T,d::D) where {T,D<:TimeType}
    idxs = searchsorted(timestamp(ta), d)
    length(idxs) == 0 && error("date $d not found")
    values(ta)[idxs[1]] = v
    ta
end

@propagate_inbounds function setindex!(ta::TimeArray{T,1,D},vals::AbstractVector{T},dates::Vector{D}) where {T,D<:TimeType}
    idxs, _ = TimeSeries.overlap(timestamp(ta), dates)
    values(ta)[idxs] = vals
    ta
end

@propagate_inbounds function setindex!(ta::TimeArray{T,2,D},vals::AbstractArray{T,2},d::D) where {T,D<:TimeType}
    idxs = searchsorted(timestamp(ta), d)
    values(ta)[idxs.start,:] = vals
    ta
end

@propagate_inbounds function setindex!(ta::TimeArray{T,2,D},vals::AbstractArray{T,2},dates::Vector{D}) where {T,D<:TimeType}
    idxs, _ = TimeSeries.overlap(timestamp(ta), dates)
    values(ta)[idxs,:] = vals
    ta
end

@propagate_inbounds setindex!(ta::TimeArray{T,N,D},vals::AbstractArray{T},dates::StepRange{D}) where {T,N,D<:TimeType} =
    setindex!(ta,vals,collect(dates))

@propagate_inbounds function setindex!(ta::TimeArray{T,N,D},v::T,d::D,c::Symbol) where {T,N,D<:TimeType}
    row = searchsorted(timestamp(ta), d)
    length(row) == 0 && error("date $d not found")
    col = TimeSeries.findcol(ta, c)
    values(ta)[row[1],col] = v
    ta
end

@propagate_inbounds setindex!(ta::TimeArray{T,N,D},vals::AbstractArray{T},d::D,cols::AbstractVector{Symbol}) where {T,N,D<:TimeType} =
    setindex!(ta,vals,[d],cols)

@propagate_inbounds setindex!(ta::TimeArray{T,N,D},vals::AbstractArray{T},dates::AbstractVector{D},c::Symbol) where {T,N,D<:TimeType} =
    setindex!(ta,vals,dates,[c])

@propagate_inbounds function setindex!(ta::TimeArray{T,N,D},vals::AbstractArray{T},dates::Colon,cols::AbstractVector{Symbol}) where {T,N,D<:TimeType}
    values(ta)[:, map(s -> TimeSeries.findcol(ta, s), cols)] = vals
    ta
end

setindex!(ta::TimeArray{T,N,D},vals::AbstractArray{T},dates::StepRange{D},cols::AbstractVector{Symbol}) where {T,N,D<:TimeType} =
    setindex!(ta,vals,collect(dates),cols)

@propagate_inbounds function setindex!(ta::TimeArray{T,N,D},vals::AbstractArray{T},dates::Vector{D},cols::AbstractVector{Symbol}) where {T,N,D<:TimeType}
    rows, _ = TimeSeries.overlap(timestamp(ta), dates)
    values(ta)[rows, map(s -> TimeSeries.findcol(ta, s), cols)] = vals
    ta
end
