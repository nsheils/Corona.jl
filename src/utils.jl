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


"""
   overlap(x1, x2,..., xn) -> i1, i2,..., in

Return the (first) indices of the common elements in the vectors
`x1, x2 ,..., xn` so that `x1[i1] = x2[i2] = x3[i3] = ... = xn[in]`.
"""
function overlap(X::Vararg{AbstractVector})
   I = intersect(X...)
   (indexin(I,x) for x in X)
end
