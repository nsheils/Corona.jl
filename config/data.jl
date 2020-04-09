##### Configuration file #####

# results path
dataconfig.results_path = "results"

# data path
dataconfig.rawdata_path = "raw"

# source names by region
push!(dataconfig.sources,
    "Germany" => "JHU",
    "Bayern" => "native",
    "Italy" => "JHU",
    "France" => "JHU",
    "China" => "JHU",
    "US" => "JHU",
    "King" => "JHU",
)

# data loaders by source
push!(dataconfig.loaders,
    "JHU" => "JHU",
    "RKI" => "native",
)

# data dirs by loader
push!(dataconfig.data_dirs,
    "JHU" => "JHU",
    "RKI" => "RKI",
)

# filename builder configuration
push!(dataconfig.args_filename_builder,
    "Bayern" => ("Bayern","Germany"),
    "King" => ("King","Washington","US"),
)

# data loader arguments by region
push!(dataconfig.args_data_loader,
    "Bayern" => "Bayern",
    "France" => (r"^$","France"),
    "King" => ("King","Washington","US"),
)

# population loader arguments by region
push!(dataconfig.args_population_loader,
    "Bayern" => ("Bayern","Germany"),
    "France" => (r"^$","France"),
    "King" => ("King","Washington"),
)

# data maps by loader
push!(dataconfig.data_maps,
    "JHU" => [0 1 1 1; 0 0 0 1],
    "native" => [0 1 1 1; 0 0 0 1],
)

# outbreak conditions by loader
push!(dataconfig.outbreak_conditions,
    "JHU" => cases -> findfirst(values(cases)[:,1].>0),
    "native" => cases -> findfirst(values(cases)[:,1].>0)
)
