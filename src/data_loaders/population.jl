using CSV

###### population loader #################

_load_population(region,args...;path="data/raw"::AbstractString) = __load_population(region,args...,path)

__load_population(region::Union{AbstractString,Regex},path::AbstractString) = __load_population(region,r"^$",path)

function __load_population(region::Union{AbstractString,Regex},subregion::Union{AbstractString,Regex},path::AbstractString)
    header = ["Subregion","Region", "Population"];
    types = Dict(1=>Union{Missing,String}, 2=>Union{Missing,String}, 3=>Int)
    fp = rstrip(path,'/')*"/population.csv"
    df = CSV.read(fp,header=header,copycols=true,ignoreemptylines=true,types=types)
    for col in [:Subregion,:Region]
        df[ismissing.(df[:,col]),col] .= ""
        df[:,col] = strip.(df[:,col])
    end
    sub = occursin.(region,df.Region) .& occursin.(subregion,df.Subregion)
    try
        df[sub,:Population][1]
    catch e
        isa(e, BoundsError) || rethrow(e)
        error("no entry in `$fp' for $region")
    end
end
