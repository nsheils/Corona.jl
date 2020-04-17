"""
   input(prompt::AbstractString="")::String

Read a string from STDIN. The trailing newline is stripped.

The prompt string, if given, is printed to standard output without a
trailing newline before reading input.
"""
function input(prompt::AbstractString="")
   print(prompt)
   chomp(readline())
end

function overlap(a::AbstractArray, b::AbstractArray)
   ab = intersect(a,b)
   loca = indexin(ab,a)
   locb = indexin(ab,b)
   loca,locb
end
