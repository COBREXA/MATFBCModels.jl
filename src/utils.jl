
guesskeys(id, model) = first(intersect(keys(model.mat), getfield(key_names, id)))

function parse_formula(x::Maybe{String})
    isnothing(x) && return nothing
    x == "" && return nothing

    res = Dict{String,Int}()
    pattern = @r_str "([A-Z][a-z]*)([1-9][0-9]*)?"
    for m in eachmatch(pattern, x)
        res[m.captures[1]] = isnothing(m.captures[2]) ? 1 : parse(Int, m.captures[2])
    end
    return res
end

function parse_charge(x)::Maybe{Int}
    if isa(x, Int)
        x
    elseif isa(x, Float64)
        Int(x)
    elseif isa(x, String)
        Int(parse(Float64, x))
    elseif isnothing(x)
        nothing
    else
        throw(DomainError(x, "cannot parse charge"))
    end
end

function unparse_formula(x::Maybe{A.MetaboliteFormula})
    isnothing(x) && return nothing
    ks = sort(collect(keys(x)))
    join(k * string(x[k]) for k in ks)
end
