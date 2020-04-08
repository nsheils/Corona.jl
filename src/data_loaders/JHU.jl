using CSV
using TimeSeries
using DataFrames

###### John Hopkins University data loader #################

_load_data(source::Val{:JHU}, region::Union{AbstractString,Regex}, path::AbstractString) = _load_data(source, (r".*", region), path)

function _load_data(source::Val{:JHU}, region::NTuple{2,Union{AbstractString,Regex}}, path::AbstractString)

    basepath = path*"/COVID-19/csse_covid_19_data/csse_covid_19_time_series"
    fns = Dict("Confirmed" => "time_series_covid19_confirmed_global.csv",
                  "Deaths" => "time_series_covid19_deaths_global.csv")

    ta = Dict{String,TimeArray}()
    for (key,fn) in fns
        fp = basepath*"/"*fn
        df = CSV.read(fp,copycols=true)
        df[ismissing.(df[:,Symbol("Province/State")]),Symbol("Province/State")] .= ""

        dates = Date.(String.(names(df)[5:end]),"m/d/y")+Dates.Year(2000)
        sub = occursin.(region[2],df[:,Symbol("Country/Region")]) .& occursin.(region[1],df[:,Symbol("Province/State")])
        if any(sub)
            ta[key] = TimeArray(dates,sum(Array{Int64}(df[sub,5:end]),dims=1)[1,:],[Symbol(key),])
        else
            error("no entry in `$fp' for $region")
        end
    end

    merge(values(ta)...)
end
