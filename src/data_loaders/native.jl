using CSV
using TimeSeries
using DataFrames

###### native data loader #################


_load_data(source::Val{:native}, region::AbstractString, path::AbstractString) = _load_data(source, (r".*", region), path)

function _load_data(source::Val{:native}, region::Tuple{Union{AbstractString,Regex},AbstractString}, path::AbstractString)

    fn = _build_filename(region[2])
    fp = rstrip(path,'/')*"/"*fn*".csv"
    if fn != region[2]
        @warn("loading data from `"*fp*"'")
    end

    header = ["Date", "Region", "Confirmed", "Deaths"];
    df = CSV.read(fp,header=header,copycols=true,dateformat="yyyy-mm-dd",ignoreemptylines=true)
    df[ismissing.(df.Region),:Region] .= ""
    df.Region = strip.(df.Region)
    sub = occursin.(region[1],df.Region)
    if any(sub)
        cases = TimeArray(df.Date[sub],Array{Int64}(df[sub,3:4]),Symbol.(header[3:4]))
    else
        error("no entry in `$fp' for $region")
    end

    cases
end
