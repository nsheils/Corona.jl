"Lagrange Interpolation, t=0:n"
function lagrange(y,t)
    tᵢ = floor(t)  ;nᵢ =Int(tᵢ )+1
    tᵢ₊= floor(t+1);nᵢ₊=Int(tᵢ₊)+1
    tᵢ₋= floor(t-1);nᵢ₋=Int(tᵢ₋)+1
    if nᵢ₋ <  1
        yᵢ =y[nᵢ]
        yᵢ₊=y[nᵢ₊]
        return yᵢ + (t-tᵢ) *  (yᵢ₊-yᵢ)/(tᵢ₊-tᵢ)
    elseif nᵢ >= length(y)
        yᵢ =y[end]
        yᵢ₋=y[end-1]
        tᵢ = Float64(length(y)-1)
        tᵢ₋= Float64(length(y)-2)
        return yᵢ + (t-tᵢ) *  (yᵢ-yᵢ₋)/(tᵢ-tᵢ₋)
    else
        yᵢ =y[nᵢ]
        yᵢ₋=y[nᵢ₋]
        yᵢ₊=y[nᵢ₊]
        l =
            yᵢ₋*(t-tᵢ )/(tᵢ₋-tᵢ ) *(t-tᵢ₊)/(tᵢ₋-tᵢ₊) + 
            yᵢ *(t-tᵢ₋)/(tᵢ -tᵢ₋) *(t-tᵢ₊)/(tᵢ -tᵢ₊) + 
            yᵢ₊*(t-tᵢ₋)/(tᵢ₊-tᵢ₋) *(t-tᵢ )/(tᵢ₊-tᵢ )
    end
end
