
function A.load(::Type{MATFBCModel}, path::String)
    id, model = first(matread(path))
    @info "Loading MAT: taking a model with ID $(id)"
    return MATFBCModel(id, model)
end

A.save(model::MATFBCModel, path::String) = matwrite(path, Dict(model.id => model.mat))

A.filename_extensions(type::Type{MATFBCModel}) = ["mat", "MAT"]
