#! /usr/bin/julia
using TimeSeries
using HTTP
root="/home/jls/data/2020-Corona/raw/RKI/"

function download_rki(source)
    r = HTTP.request("GET", source)
    String(r.body)
end

function parse_rki(pattern,file)
    m=match(pattern,file).match
    all=collect(eachmatch(r">([\d\.]+)<",m))
    confirmed=parse(Float64,replace(all[1][1],"." => ""))
    deaths=parse(Float64,replace(all[3][1],"." => ""))
    return confirmed,deaths
end

function date_rki(file)
    datum=replace(split(match(r"(Stand:.*),",file).match)[2],"," =>"")
    dd=parse(Int,match(r"^(\d+)",datum)[1])
    mm=parse(Int,match(r"\.(\d+)\.",datum)[1])
    yyyy=parse(Int,match(r"(\d+)$",datum)[1])
    d=Date(yyyy,mm,dd)
end




RKI="https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Fallzahlen.html"
rki_file=download_rki(RKI);
d=date_rki(rki_file)



laender=Dict([
    ("Baden-Würtemberg"    ,r"Baden.*Bayern")
    ("Bavaria"             ,r"Bayern.*Berlin")
    ("Berlin"              ,r"Berlin.*Brandenburg")
    ("Rheinland-Pfalz"     ,r"Rhein.*Saar")
    ("Nordrhein-Westphalen",r"Nordrhein-West.*Rhein")
])

for land in laender
    datafile=land[1]*".csv"
#    println(datafile)
    if isfile(root*datafile)
        dataarray=readtimearray(root*datafile,delim=',')
        LastDate=timestamp(dataarray)[end]
    else
        dataarray=TimeArray()
        LastDate=Date(2000)
    end
#    println(LastDate)
    confirmed,deaths=parse_rki(land[2],rki_file)
    #    println(datafile," ",LastDate," $confirmed, $deaths ")
    if d > LastDate
        dataarray=update(dataarray,d ,[confirmed deaths])
        writetimearray(dataarray,root*datafile)
    end
end



