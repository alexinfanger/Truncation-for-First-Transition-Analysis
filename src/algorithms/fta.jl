function fta_lb(
    G::AbstractArray,
    Aidx::Vector{Int64}, 
    f::Vector{Float64})
    u_l = (I-G[Aidx, Aidx])\f[Aidx]
    return u_l
end

function fta_ub_complete(
    G::AbstractArray, 
    Aidx::Vector{Int64}, 
    Acidx::Vector{Int64}, 
    f::Vector{Float64}, 
    h::Vector{Float64})
    
    u_l = fta_lb(G, Aidx, f)
    k = G[Aidx, Acidx]*h[Acidx]    
    u_u = u_l + (I-G[Aidx, Aidx])\k
    return u_u
end


function fta_ub_incomplete(
    G::AbstractArray, 
    Aidx::Vector{Int64}, 
    Acidx::Vector{Int64}, 
    f::Vector{Float64}, 
    w::Vector{Float64})

    u_l = fta_lb(G, Aidx, f)
    k = G[Aidx, Acidx]*w[Acidx]
    e_S = ones(length(w))
    gamma = norm((I-G[Aidx, Aidx])\(G[Aidx, Acidx]*e_S[Acidx]), Inf)
    @assert gamma < 1
    beta = (1/(1-gamma))*norm((I-G[Aidx, Aidx])\(f[Aidx]+k), Inf)
    u_u = u_l + (I-G[Aidx, Aidx])\(k+beta*G[Aidx, Acidx]*e_S[Acidx])
    return u_u
end
