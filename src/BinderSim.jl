module BinderSim
export compute_binder
compute_binder(; L::Int, lambda_x::Float64, lambda_zz::Float64) =
    2/3 + 0.3*tanh((lambda_x - lambda_zz)*L/50)
end
