using CSV
using TimeSeries
using DataFrames

###### John Hopkins University data loader #################

_load_data(source::Val{:JHU}, region::Union{AbstractString,Regex}, path::AbstractString) = _load_data(source, (r".*", region, "global"), path)

_load_data(source::Val{:JHU}, region::NTuple{2,Union{AbstractString,Regex}}, path::AbstractString) = _load_data(source, (region..., "global"), path)

function _load_data(source::Val{:JHU}, region::Tuple{Union{AbstractString,Regex},Union{AbstractString,Regex},AbstractString}, path::AbstractString)

    basepath = joinpath(path,"COVID-19/csse_covid_19_data/csse_covid_19_time_series")
    if region[3] == "global"
        fns = Dict(:Confirmed => "time_series_covid19_confirmed_global.csv",
                      :Deaths => "time_series_covid19_deaths_global.csv")
        region_cols = map(Symbol,("Province/State","Country/Region"))
        date_cols_s = 5
    elseif region[3] == "US"
        fns = Dict(:Confirmed => "time_series_covid19_confirmed_US.csv",
                      :Deaths => "time_series_covid19_deaths_US.csv")
        region_cols = map(Symbol,("Admin2","Province_State"))
        date_cols_s = 13
    else
        error("region $region unknown")
    end

    ta = Dict{Symbol,TimeArray}()
    for (key,fn) in fns
        fp = joinpath(basepath,fn)
        df = CSV.read(fp,copycols=true)
        let col = region_cols[1]
            rows = ismissing.(df[!,col]) .| (df[!,col] .== "Unassigned")
            df[rows,col] .= ""
        end

        dates = Date.(String.(names(df)[date_cols_s:end]),"m/d/y") + Dates.Year(2000)
        sub = trues(size(df,1))
        for (i,col) in enumerate(region_cols)
            sub = sub .& occursin.(region[i],df[!,col])
        end
        if any(sub)
            ta[key] = TimeArray(dates,sum(Array{Int64}(df[sub,date_cols_s:end]),dims=1)[1,:],[key,])
        else
            error("no entry in `$fp' for $region")
        end
    end

    merge(values(ta)...)
end
