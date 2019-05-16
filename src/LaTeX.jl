
#===========================================
Code to export parameters table `t`
and potentially optimal parameters `popt`
to a LaTeX format
===========================================#
module LaTeX

using Printf



function print_header(str_lengths)
l_L, l_D, l_U = str_lengths
headtop =
"""
\\begin{table*}[t]
    \\caption{
    }
    \\begin{tabular}{lp{10cm}lll}
        \\tophline
"""
str_L1 = rpad("", l_L+2)
str_L2 = rpad("Symbol", l_L+2)
str_D1 = rpad("", l_D)
str_D2 = rpad("Description", l_D)
str_U1 = rpad("", l_U)
str_U2 = rpad("Unit", l_U)
str_IV1 = rpad("Initial", 8)
str_IV2 = rpad("value", 8)
str_FV1 = rpad("Optimal", 8)
str_FV2 = rpad("value", 8)
str_1 = string("        ", str_L1, " & ", str_D1, " & ", str_IV1, " & ", str_FV1, " & ", str_U1, " \\\\\n")
str_2 = string("        ", str_L2, " & ", str_D2, " & ", str_IV2, " & ", str_FV2, " & ", str_U2, " \\\\\n")
headbot =
"""
        \\middlehline
"""
return string(headtop, str_1, str_2, headbot)
end

function print_footer()
"""
        \\bottomhline
    \\end{tabular}
\\end{table*}
"""
end

function print_row(t, p, i, str_lengths)
    l_L, l_D, l_U = str_lengths
    str_L = rpad(string("\$", t[i, :LaTeX], "\$"), l_L+2)
    str_D = rpad(t[i, :description], l_D)
    str_U = rpad(latexify(t[i, :printunit]), l_U)
    str_IV = @sprintf "%8.2f" ustrip(t[i, :value] * t[i, :unit] |> t[i, :printunit])
    str_FV = @sprintf "%8.2f" ustrip(getfield(p, fieldnames(typeof(p))[i]) * t[i, :unit] |> t[i, :printunit])
    str_FV = t[i, :optimizable] ? str_FV : "        "
    return string("        ", str_L, " & ",  str_D, " & ", str_IV, " & ", str_FV, " & ", str_U, " \\\\\n")
end

function string_lengths(t)
    l_LaTeX = maximum(length.(t[:LaTeX]))
    l_Definition = maximum(length.(t[:description]))
    l_unit = maximum(length.([latexify(u) for u in t[:printunit]]))
    return l_LaTeX, l_Definition, l_unit
end



function print_table(t, p)
    body = ""
    str_lengths = string_lengths(t)
    for irow in 1:size(t, 1)
        body = string(body, print_row(t, p, irow, str_lengths))
    end
    return string(print_header(str_lengths), body, print_footer())
end

function latexify(U)
    str = string(U)
    str = replace(str, r"\^-(?<exp>\d+)" => s"^{-\g<exp>}") # add brackets around exponents
    str = replace(str, r"\s" => s"\\,") # replace spaces with small spaces
    str = replace(str, r"yr" => s"a") # replace spaces with small spaces
    str = replace(str, r"μ" => s"\\mu{}") # replace spaces with small spaces
    return string("\\unit{", str, "}")
end



end # module
