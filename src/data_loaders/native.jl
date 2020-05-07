using CSV
using TimeSeries
using DataFrames

###### native data loader #################

function _load_data(source::Val{:native}, region::AbstractString, path::AbstractString)
    # fn = _build_filename(region)
    fn = region
    fp = rstrip(path,'/')*"/"*fn*".csv"
    if fn != region
        @warn("loading data from `"*fp*"'")
    end
    header = ["Date", "Confirmed", "Deaths"];
    df = CSV.read(fp,header=header,datarow=2,dateformat="yyyy-mm-dd",ignoreemptylines=true)
    TimeArray(df.Date,Array{Int64}(df[!,2:3]),Symbol.(header[2:3]))
end
