using JLD2, FileIO
using DataFrames
using DataStructures
using Unicode

###### structures #################

mutable struct DataConfig
  rawdata_path::String
  results_path::String
  sources::DefaultDict{AbstractString,String}
  loaders::DefaultDict{AbstractString,String}
  data_dirs::DefaultDict{AbstractString,String}
  args_filename_builder::DefaultDict{AbstractString,Any}
  args_data_loader::DefaultDict{AbstractString,Any}
  args_population_loader::DefaultDict{AbstractString,Any}
  data_maps::DefaultDict{AbstractString,Matrix{<:Real}}
  outbreak_conditions::DefaultDict{AbstractString,Function}
    function DataConfig()
      new(
        "raw", # rawdata_path
        "results", # results_path
        DefaultDict{AbstractString,String}("native"), # sources
        DefaultDict{AbstractString,String}("native"), # loaders
        DefaultDict{AbstractString,String}(""), # data_dirs
        DefaultDict{AbstractString,Any}(k -> k; passkey=true), # args_filename_builder
        DefaultDict{AbstractString,Any}(k -> k; passkey=true), # args_data_loader
        DefaultDict{AbstractString,Any}(k -> k; passkey=true), # args_population_loader
        DefaultDict{AbstractString,Matrix{<:Real}}(Diagonal(ones(Int,4))), # data_maps
        DefaultDict{AbstractString,Function}(c -> findfirst(values(c[:,1].>0))), # outbreak_conditions
      )
    end
end

struct DataStruct{T<:Integer,D<:TimeType}
    cases::TimeArray{T,2,D}
    outbreakdate::Union{D,Missing}
    lastdate::D
    population::T
    map::Matrix{<:Real}
    source::String
    function DataStruct(cases::TimeArray{T,2,D},
                        outbreakdate::Union{D,Missing},
                        lastdate::D,
                        population::T,
                        datamap::Matrix{<:Real},
                        source::AbstractString) where {T<:Integer,D<:TimeType}
        new{T,D}(cases,outbreakdate,lastdate,population,datamap,source)
    end
end

###### filename builder #################
function build_filename(config::DataConfig,region::AbstractString)
    _build_filename(config.args_filename_builder[region])
end

function _build_filename(region::NTuple{N,AbstractString}) where {N}
    fn = _build_filename(region[1])
    for i=2:N
        fn = string(_build_filename(region[i]),"-",fn)
    end
    fn
end

_build_filename(region::NTuple{2,AbstractString}) =
        _build_filename(region[2])*"-"*_build_filename(region[1])

function _build_filename(region::AbstractString)
    char_map = Dict('ä' => "ae", 'ö' => "oe", 'ü' => "ue", 'ß' => "ss", ' ' => "_")
    fn = region
    fn = Unicode.normalize(fn,compat=true,stripcc=true,rejectna=true)
    fn = join(haskey(char_map,c) ? char_map[c] : c for c in fn)
    fn = titlecase(fn,strict=false)
    fn = Unicode.normalize(fn,stripmark=true)
    fn = replace(fn, r"[\W]" => "")
    return fn
end

###### data loaders #################

function load_data(config::DataConfig,region::AbstractString)
    source = config.sources[region]
    path = config.rawdata_path
    loader = config.loaders[source]
    cases = _load_data(Val(Symbol(loader)),config.args_data_loader[region],
                        joinpath(path,config.data_dirs[loader]))
    outbreakdate = get_outbreak_date(config,cases,loader)
    lastdate = timestamp(cases)[end]
    population = _load_population(config.args_population_loader[region],path)
    datamap = config.data_maps[loader]
    DataStruct(cases,outbreakdate,lastdate,population,datamap,source)
end

include("data_loaders/native.jl")
include("data_loaders/population.jl")
include("data_loaders/JHU.jl")

###### outbreak date #################
function get_outbreak_date(config::DataConfig,
                           cases::TimeArray,
                           loader::AbstractString)
    outbreakindex = config.outbreak_conditions[loader](cases)
    if isnothing(outbreakindex)
        outbreakdate = missing
    else
        outbreakdate = timestamp(cases)[outbreakindex]
    end
    return outbreakdate
end

###### save/load results #################

function save(config::DataConfig,region::AbstractString,da::DA;
              interactive=false::Bool)
    path = joinpath(config.results_path,build_filename(config,region))
    mkpath(path)
    filename = joinpath(path,"da.jld2")
    if interactive && isfile(filename)
        if input("overwrite `$filename'? ") != "y"
            return
        end
    end
    _save(filename,da)
end

function _save(filename::AbstractString,da::DA)
    FileIO.save( filename,
                   "data", da.data,
                "datamap", da.datamap,
                 "window", da.window,
                 "result", da.result,
           "model_params", da.model_params )
end

function load_model_params(config::DataConfig,region::AbstractString)
    filename = joinpath(config.results_path,build_filename(config,region),"da.jld2")
    FileIO.load(filename,"model_params")
end
