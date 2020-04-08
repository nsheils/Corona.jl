##### Configuration file #####

# dump path
dataconfig.dump_path = "dump"

# data path
dataconfig.data_path = "data"

# source names by region
push!(dataconfig.sources,
    "Germany" => "JHU",
    "Bayern" => "native",
    "Italy" => "JHU",
    "France" => "JHU",
    "China" => "JHU",
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
    "Bayern" => ("Germany","Bayern"),
)

# data loader arguments by region
push!(dataconfig.args_data_loader,
    "Bayern" => "Bayern",
    "France" => ("France",r"^$"),
)

# population loader arguments by region
push!(dataconfig.args_population_loader,
    "Bayern" => ("Germany","Bayern"),
    "France" => ("France",r"^$"),
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
