using CSV
using TimeSeries
using DataFrames

###### native data loader #################

_load_data(source::Val{:native}, args...; path="data/raw/native"::AbstractString) = __load_data(source, args..., path)

__load_data(source::Val{:native}, region::AbstractString, path::AbstractString) = __load_data(source, region, r".*", path)

function __load_data(source::Val{:native}, region::AbstractString, subregion::Union{AbstractString,Regex}, path::AbstractString)

    fn = _build_filename(region)
    fp = rstrip(path,'/')*"/"*fn*".csv"
    if fn != region
        @warn("loading data from `"*fp*"'")
    end

    header = ["Date", "Region", "Confirmed", "Deaths"];
    df = CSV.read(fp,header=header,copycols=true,dateformat="yyyy-mm-dd",ignoreemptylines=true)
    df[ismissing.(df.Region),:Region] .= ""
    df.Region = strip.(df.Region)
    sub = occursin.(subregion,df.Region)
    if any(sub)
        cases = TimeArray(df.Date[sub],Array{Int64}(df[sub,3:4]),Symbol.(header[3:4]))
    else
        error("no entry in `$fp' for $region")
    end

    cases
end
