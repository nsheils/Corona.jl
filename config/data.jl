##### Data configuration file #####

# source names by region
push!(dataconfig.sources,
    "Germany" => "JHU",
    "Bayern" => "RKI",
    "NRW" => "RKI",
    "Italy" => "JHU",
    "France" => "JHU",
    "China" => "JHU",
    "US" => "JHU",
    "King" => "JHU",
    "New York" => "JHU",
    "United Kingdom" => "JHU",
    "Miami" => "JHU"
)

# data loaders by source
push!(dataconfig.loaders,
    "JHU" => "JHU",
    "RKI" => "native",
)

# data dirs by source
push!(dataconfig.data_dirs,
    "native" => "native",
    "JHU" => "JHU",
    "RKI" => "RKI",
)

# filename builder configuration
push!(dataconfig.args_filename_builder,
    "Bayern" => ("Bayern","Germany"),
    "NRW" => ("Nordrhein-Westphalen","Germany"),
    "King" => ("King","Washington","US"),
    "New York" => ("New York","New York","US"),
    "Miami" => ("Miami-Dade","Florida","US"),
    "Laputa" => ("Laputa","Balnibarbi"),
)

# data loader arguments by region
push!(dataconfig.args_data_loader,
    "Bayern" => "Bayern",
    "NRW" => "Nordrhein-Westphalen",
    "France" => (r"^$","France"),
    "King" => ("King","Washington","US"),
    "New York" => ("New York","New York","US"),
    "Miami" => ("Miami-Dade","Florida","US"),
)

# population loader arguments by region
push!(dataconfig.args_population_loader,
    "Bayern" => ("Bayern","Germany"),
    "NRW" => ("Nordrhein-Westphalen","Germany"),
    "France" => (r"^$","France"),
    "King" => ("King","Washington"),
    "New York" => ("New York","New York"),
    "Miami" => ("Miami-Dade","Florida"),
    "Laputa" => ("Laputa","Balnibarbi"),
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
