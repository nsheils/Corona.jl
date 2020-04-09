using Corona
using Dates
using LinearAlgebra
using TimeSeries
using Formatting
using june

### Parameters
maxiters     = 10000000;
tolerance    = 0.001;
filterfreq   = 10;
screenfreq   = 100;
sponge       = 5;
maxlineiters = 4096;
maxlinerange = 1e-6;
coldstart    = true;
region       = "King";

############################################################
printstyled("C O R O N A",bold=true,color=:blue)
print("   $region\n\n")

### Load configuration
dataconfig = Corona.DataConfig();
include("../config/data.jl");

### Load data
data = Corona.load_data(dataconfig,region);

# Print some information
print("Data source : $(data.source)\n")
print("Outbreakdate : $(data.outbreakdate)\n")
print("Last date : $(data.lastdate)\n\n")

###  Set initial state vector
S0 = data.population;
C0,D0 = values(data.cases[data.outbreakdate]);
init_vals = [S0, C0, 0, D0];

### Set or load initial model parameters
if coldstart
    init_mp = Dict(:β => 1//4, :γ => 1//7, :δ => 0);
else
    init_mp = Corona.load_model_params(dataconfig,region);
end

### Initialize data assimilation
da = Corona.DA(data.cases, init_vals,
               init_mp, datamap = data.map,
               start = data.outbreakdate,
               stop = timestamp(data.cases)[end] + Day(sponge)
              );


### Initialize residual
J = Array{Float64,1}()
J0 = norm(values(da.data).*values(da.window))

α=maxlinerange/J0;

### Print some other information
if coldstart
    printstyled("Cold Start";color=:red)
else
    printstyled("Warm Start";color=:green)
end
printfmt(" with α= {:.4e} and tolerance {}\n",α,tolerance)
printfmt("└ maximum number of iterations is {}\n\n",maxiters)

## Test
# Corona.forward!(da);
# v = Corona.backward(da);
# J = Corona.linesearch!(da,v,α=α,method="plot",maxiters=maxlineiters)
# error()

### Exectue data assimilation
Corona.forward!(da);
try
    for i=1:maxiters
        v = Corona.backward(da);
        success,Je,_ = Corona.linesearch!(da,v,α=α,method="bisection",maxiters=maxlineiters);
        if !success
            success,Je,_ = Corona.linesearch!(da,v,α=α,method="brute_force",maxiters=10*maxlineiters);
        end
        Corona.propagate_solution!(da);
        Corona.forward!(da);

        if mod(i,screenfreq) == 0
            global J = [J; Je/J0]
            if length(J) == argmin(J)
                color=:green
            else
                color=:red
            end
            print(i," ")
            printstyled(format("{:.8f}",J[end]),"\n";color=color)
        end
        if mod(i,filterfreq) == 0
            values(da.model_params)[:,:] = ∂⁰(values(da.model_params))
        end
        if Je < J0*tolerance
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
Corona.save(dataconfig,region,da,interactive=true)
