using CSV
using TimeSeries
using DataFrames

###### John Hopkins University data loader #################

_load_data(source::Val{:JHU}, args...; path="data/raw/JHU") = __load_data(source, args..., path)
__load_data(source::Val{:JHU}, region::Union{AbstractString,Regex}, path::AbstractString) = __load_data(source, region, r".*", path)

function __load_data(source::Val{:JHU}, region::Union{AbstractString,Regex}, subregion::Union{AbstractString,Regex}, path::AbstractString)

    basepath = path*"/csse_covid_19_data/csse_covid_19_time_series"
    fns = Dict("Confirmed" => "time_series_covid19_confirmed_global.csv",
                  "Deaths" => "time_series_covid19_deaths_global.csv")

    ta = Dict{String,TimeArray}()
    for (key,fn) in fns
        fp = basepath*"/"*fn
        df = CSV.read(fp,copycols=true)
        df[ismissing.(df[:,Symbol("Province/State")]),Symbol("Province/State")] .= ""

        dates = Date.(String.(names(df)[5:end]),"m/d/y")+Dates.Year(2000)
        sub = occursin.(region,df[:,Symbol("Country/Region")]) .& occursin.(subregion,df[:,Symbol("Province/State")])
        if any(sub)
            ta[key] = TimeArray(dates,sum(Array{Int64}(df[sub,5:end]),dims=1)[1,:],[Symbol(key),])
        else
            error("no entry in `$fp' for $subregion, $region")
        end
    end

    merge(values(ta)...)
end
