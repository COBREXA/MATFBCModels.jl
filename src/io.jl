
function A.load(::Type{MATFBCModel}, path::String)
    model_pair = first(matread(path))
    @info "Loading MAT: taking a model with ID $(model_pair.first)"
    return MATFBCModel(model_pair.second)
end

function A.save(model::MATFBCModel, path::String; model_name = "model")
    m =
        typeof(model) == MATFBCModel ? model :
        begin
            @warn "Automatically converting $(typeof(model)) to MATFBCModel for saving, information may be lost."
            convert(MATFBCModel, model)
        end
    matwrite(path, Dict(model_name => m.mat))
end

A.filename_extensions(type::Type{MATFBCModel}) = ["mat", "MAT"]
