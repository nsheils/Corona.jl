using Corona
using Flux
using Dates
using LinearAlgebra
using TimeSeries
using Formatting

### Parameters
maxiters     = 1000000;
tolerance    = 0.001;
screenfreq   = 100;
sponge       = 30;
opt          = Momentum(1e-13, 0.99);
coldstart    = true;
region       = "Bayern";

############################################################
printstyled("C O R O N A",bold=true,color=:blue)
print("   $region\n\n")

### Load configuration
dataconfig = Corona.DataConfig();
include("../config/data.jl");
if isfile(joinpath(Base.@__DIR__,"../config/paths.jl"))
    include("../config/paths.jl")
end;

### Load data
data = Corona.load_data(dataconfig,region);

# Print some information
print("Data source : $(data.source)\n")
print("Outbreakdate : $(data.outbreakdate)\n")
print("Last date : $(data.lastdate)\n\n")

###  Set initial state vector
S₀ = data.population;
C₀,D₀ = values(data.cases[data.outbreakdate]);
init_u = Dict(:S => S₀, :I => C₀, :R => 0, :D => D₀);

### Set or load initial model parameters
if coldstart
    init_p = Dict(:β => 1//4, :γ => 1//7, :δ => 0)
else
    init_p = Corona.load_model_params(dataconfig,region)
end;

### Initialize data assimilation
base = Corona.Baseline(data.cases, init_u, init_p,
                C = data.map, start = data.outbreakdate,
                stop = timestamp(data.cases)[end] + Day(sponge));

### Print some other information
if coldstart
    printstyled("Cold Start";color=:red)
else
    printstyled("Warm Start";color=:green)
end
printfmt(" with tolerance {}\n",tolerance)
printfmt("| maximum number of iterations {}\n",maxiters)
print("└ and optimiser $opt\n\n")

### Exectue data assimilation
da_run = Corona.DA(base);
da_opt = da_run;

### Initialize residual
J = Vector{Union{Missing,Float64}}(missing,maxiters);
J_ini = norm(values(da_run.data).*values(da_run.σ));
J_min = Inf;
J_argmin = 1;

### Optimise data assimilation
try
    for i=1:maxiters
        global da_run, da_opt, J_min, J_argmin, convegence

        da_run = Corona.apply!(opt,da_run)
        J[i] = Corona.residual(da_run)/J_ini
        Corona.extend_solution!(da_run)

        if J[i] < J_min
            convegence = (log10(J_min) - log10(J[i])) / (log10(i) - log10(J_argmin))
            J_min = J[i]
            J_argmin = i
            da_opt = da_run
        end
        if mod(i,screenfreq) == 0
            iscreen = i-screenfreq+1:i
            if argmin(skipmissing(J)) in iscreen
                color=:green
                n_old=i;
            else
                color=:red
            end
            print(i," ")
            printstyled(format("{:.8f}",minimum(J[iscreen]));color=color)
            printstyled(format(" {:.8f}\n",convegence);color=:blue)
        end
        if J[i] < tolerance
            printstyled("\nSTOP: ";bold=true,color=:green)
            printstyled("Residual within prescribed tolerance\n";color=:green)
            break
        elseif i == maxiters
            printstyled("\nSTOP: ";bold=true,color=:yellow)
            printstyled("Maximum number of iterations reached\n";color=:yellow)
        end
    end
catch e
    isa(e, InterruptException) || rethrow(e);
    printstyled("\nSTOP: ";bold=true)
    print("Keyboard interrupt\n")
end

### Save data assimilation result
Corona.save(dataconfig,region,da_opt,interactive=true)
