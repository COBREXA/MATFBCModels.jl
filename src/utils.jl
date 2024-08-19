
guesskeys_maybe(id, model) =
    let keys = [k for k in getfield(key_names, id) if k in keys(model.mat)]
        isempty(keys) ? nothing : first(keys)
    end

guesskeys(id, model) =
    let key = guesskeys_maybe(id, model)
        isnothing(key) ? throw(KeyError, "no applicable model key found for $id") : key
    end

function parse_compartment(x::String)
    isempty(x) && return nothing
    return x
end

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
        return x
    elseif isa(x, Float64)
        return isnan(x) ? nothing : Int(x)
    elseif isa(x, String)
        return Int(parse(Float64, x))
    elseif isnothing(x)
        return nothing
    else
        throw(DomainError(x, "cannot parse charge"))
    end
end

unparse_compartment(::Nothing) = ""
unparse_compartment(x::String) = x

function unparse_formula(x::Maybe{A.MetaboliteFormula})
    isnothing(x) && return nothing
    ks = sort(collect(keys(x)))
    join(k * string(x[k]) for k in ks)
end

unparse_charge(::Nothing) = NaN
unparse_charge(x::Int) = float(x)
