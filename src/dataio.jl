using JLD2, FileIO
using DataFrames
using DataStructures
using Unicode

###### structures #################

mutable struct DataConfig
  sources::DefaultDict{AbstractString,String}
  loaders::DefaultDict{AbstractString,String}
  args_filename_builder::DefaultDict{AbstractString,Tuple}
  args_data_loader::DefaultDict{AbstractString,Tuple}
  args_population_loader::DefaultDict{AbstractString,Tuple}
  data_paths::DefaultDict{AbstractString,String}
  data_maps::DefaultDict{AbstractString,Matrix{<:Real}}
  outbreak_conditions::DefaultDict{AbstractString,Function}
  save_path::String
    function DataConfig()
      new(
        DefaultDict{AbstractString,String}("native"), # sources
        DefaultDict{AbstractString,String}("native"), # loaders
        DefaultDict{AbstractString,Tuple}(k -> (k,); passkey=true), # args_filename_builder
        DefaultDict{AbstractString,Tuple}(k -> (k,); passkey=true), # args_data_loader
        DefaultDict{AbstractString,Tuple}(k -> (k,); passkey=true), # args_population_loader
        DefaultDict{AbstractString,String}(""), # data_paths
        DefaultDict{AbstractString,Matrix{<:Real}}(Diagonal(ones(Int,4))), # data_maps
        DefaultDict{AbstractString,Function}(c -> findfirst(values(c[:,1].>0))), # outbreak_conditions
        "data/results", # save_path
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
function build_filename(config::DataConfig,
                        region::AbstractString;
                        path=""::AbstractString,
                        ext=""::AbstractString)
    fn = _build_filename(config.args_filename_builder[region]...)
    if !isempty(path)
        fn = rstrip(path,'/')*"/"*fn
    end
    if !isempty(ext)
        fn = fn*"."*lstrip(ext)
    end
    return fn
end

_build_filename(region::AbstractString,subregion::AbstractString) =
        _build_filename(region)*"-"*_build_filename(subregion)

function _build_filename(region::AbstractString)
    char_map = Dict('ä' => "ae", 'ö' => "oe", 'ü' => "ue", 'ß' => "ss", ' ' => "_")
    fn = region
    fn = Unicode.normalize(fn,compat=true,stripcc=true,rejectna=true)
    fn = join(haskey(char_map,c) ? char_map[c] : c for c in fn)
    fn = titlecase(fn)
    fn = Unicode.normalize(fn,stripmark=true)
    fn = replace(fn, r"[\W]" => "")
    return fn
end

###### data loaders #################

function load_data(config::DataConfig,region::AbstractString)
    source = config.sources[region]
    path = config.data_paths[source]
    loader = config.loaders[source]
    cases = _load_data(Val(Symbol(loader)),config.args_data_loader[region]...;
            (isempty(path) ? () : (path = path,))...)
    outbreakdate = get_outbreak_date(config,cases,loader)
    lastdate = timestamp(cases)[end]
    population = _load_population(config.args_population_loader[region]...)
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
    fn = build_filename(config,region)
    fp = joinpath(config.save_path,fn)*".jld2"
    if interactive && isfile(fp)
        if input("overwrite `$fp'? ") != "y"
            return
        end
    end
    _save(fp,da)
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
    fn = build_filename(config,region)
    fp = joinpath(config.save_path,fn)*".jld2"
    load(fp,"model_params")
end
