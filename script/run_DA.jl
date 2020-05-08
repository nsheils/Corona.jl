using Corona
using Flux
using Dates
using LinearAlgebra
using Interpolations
using TimeSeries
using Formatting

### Parameters
maxiters     = 100000;
tolerance    = 0.001;
screenfreq   = 10;
sponge       = 90;
opt          = Momentum(1e-3, 0.99);
# opt          = NADAM(1e-12, (0.89, 0.995))
coldstart    = true;
c₀           = Day(0);
#c₀           = missing; -> lastdate
Δc           = Day(1000);
region       = "Laputa";

############################################################
printstyled("C O R O N A",bold=true,color=:blue)
print("   $region\n\n")

### Load configuration
dataconfig = Corona.DataConfig();
include(joinpath(pwd(),"config/data.jl"));
if isfile(joinpath(pwd(),"config/paths.jl"))
    include(joinpath(pwd(),"config/paths.jl"))
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

### Create windows
if ismissing(c₀)
    d₀ = data.lastdate
else
    d₀ = data.outbreakdate + c₀
end;

assimtime = DateTime(d₀)-13*Millisecond(Δc)/10:Hour(1):DateTime(d₀)+13*Millisecond(Δc)/10;
w = exp.(-((assimtime.-DateTime(d₀))./Millisecond(Δc)).^10);
w[assimtime .> DateTime(data.lastdate)] .= 0.0;

σ = TimeArray(assimtime,repeat(w,1,2),[:Confirmed,:Deaths]);
μ = TimeArray(assimtime,repeat(w,1,3),[:β,:γ,:δ]);

### Initialize data assimilation
base = Corona.Baseline(data.cases, init_u, init_p;
                C = data.map, start = data.outbreakdate,
                stop = timestamp(data.cases)[end] + Day(sponge),
                step = Hour(1), σ = σ, μ = μ,
                interptype = SteffenMonotonicInterpolation());

### Print some other information
if coldstart
    printstyled("Cold Start";color=:red)
else
    printstyled("Warm Start";color=:green)
end
printfmt(" with tolerance {}\n",tolerance)
printfmt("| maximum number of iterations {}\n",maxiters)
print("└ and optimiser $opt\n\n")

println("Assimilation time: $assimtime\n")

### Initialize some variables
J = Vector{Union{Missing,Float64}}(missing,maxiters);
J_ini = Corona.datanorm(base);
J_min = Inf;
J_argmin = 1;
i_rof = 0;

### Exectue data assimilation
da_run = Corona.DA(base);
da_opt = da_run;

### Optimise data assimilation
try
    for i=1:maxiters
        global da_run, da_opt, J_min, J_argmin, convegence, i_rof

        da_run = Corona.apply!(opt,da_run)
        J[i] = Corona.residual(Corona.baseline(da_run))/J_ini
        Corona.extend_solution!(da_run)

        if J[i] < J_min
            J_min = J[i]
            J_argmin = i
            da_opt = da_run
        end
        if mod(i,screenfreq) == 0
            iscreen = i-screenfreq+1:i
            if argmin(skipmissing(J)) in iscreen
                color=:green
            else
                color=:red
            end
            print(i," ")
            printstyled(format("{:.8f}",minimum(J[iscreen]));color=color)
            printstyled(format(" {:.8f}\n",J_min);color=:blue)
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

### Useful for plotting
da = da_opt;

### Save data assimilation result
Corona.save(dataconfig, region, da, data = data, interactive = true);
