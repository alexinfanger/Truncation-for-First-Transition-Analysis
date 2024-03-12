using LinearAlgebra



struct OpenJacksonNetwork <:InfiniteMarkovChain
    num_queues::Int64
    num_states::Array{Int64,1}
    num_servers::Array{Int64,1}
    mus::AbstractArray
    Arrivals::AbstractArray
    RoutingMatrix::AbstractArray
end

function get_lambda(jn::OpenJacksonNetwork)
    ld = (I-jn.RoutingMatrix')\jn.Arrivals
    return ld
end

function get_stat_dist(jn::OpenJacksonNetwork)
    ld = get_lambda(jn)
    if any(ld.>=jn.mus)
        ArgumentError("At least one lambda >= mu.")
    end
    stat_dists = zeros(jn.num_queues, maximum(jn.num_states)+1)
    for i = 1:jn.num_queues
        stat_dists[i,:]=get_stat_dist_Q(MMs_Q(ld[i], jn.mus[i], jn.num_states[i], jn.num_servers[i]))
    end
    return stat_dists
end

struct TruncatedOpenJacksonNetwork <: FiniteMarkovChain
    nqueues::Int64
    array_of_sizes::Array{Int64,1}
    nservers::Array{Int64,1}
    mus::Vector{Float64}
    arrivals::Vector{Float64}
    RoutingMatrix::Matrix{Float64}
    vss::VectorStateSpace
    Q::AbstractArray
end


function product_form_to_vec(vss::VectorStateSpace, product_array::AbstractArray)
    piv = zeros(vss.size)
    for i = 1:vss.size
        vec = ind_to_vec(vss, i)
        prod = 1.0
        for jp = 1:size(product_array, 1)
            prod *= product_array[jp, vec[jp]]
        end
        piv[i] = prod
    end
    return piv
end


function TruncatedOpenJacksonNetwork(array_of_sizes::Vector{Int64},
                                    nservers::Vector{Int64},
                                    mus::Vector{Float64},
                                    arrivals::Vector{Float64},
                                    RoutingMatrix::Matrix{Float64})

    vss = VectorStateSpace(array_of_sizes)
    rows = []
    cols = []
    vals = Array{Float64, 1}(UndefInitializer(), 0)

    rsums = vec(sum(RoutingMatrix,dims=2))

    for i=1:vss.size

        s1 = ind_to_vec(vss,i)
        sm = 0.0

        for kp = 1:vss.ndims
            alpha = min(s1[kp]-1.0, nservers[kp])*mus[kp]
            for k2p=1:vss.ndims
                if k2p == kp || s1[kp]-1.0 == 0 ||
                             s1[k2p] == array_of_sizes[k2p]
                    continue
                end
                s2 = s1+e_i(k2p, vss.ndims)-e_i(kp, vss.ndims)
                j = vec_to_ind(vss, s2)
                append!(rows, i)
                append!(cols, j)
                append!(vals, alpha*RoutingMatrix[kp,k2p])
                sm += alpha*RoutingMatrix[kp,k2p]
            end

            if s1[kp]-1.0 != 0 && (1.0-rsums[kp])>0
                s2 = s1 - e_i(kp, vss.ndims)
                j = vec_to_ind(vss, s2)
                append!(rows, i)
                append!(cols, j)
                append!(vals, alpha*(1.0-rsums[kp]))
                sm += alpha*(1.0-rsums[kp])
            end

            if s1[kp] <= array_of_sizes[kp]-1 && arrivals[kp]>0
                s2 = s1 + e_i(kp, vss.ndims)
                j = vec_to_ind(vss, s2)
                append!(rows, i)
                append!(cols, j)
                append!(vals, arrivals[kp])
                sm += arrivals[kp]
            end

        end


        append!(rows, i)
        append!(cols, i)
        append!(vals, -sm)
    end

    Q = sparse(rows,cols, vals, vss.size, vss.size)
    # P,ld = uniformize(Q)

    return TruncatedOpenJacksonNetwork(size(array_of_sizes,1), array_of_sizes, nservers, mus,
                                arrivals, RoutingMatrix, vss, Q)
end










# array_of_sizes = [25,25,25]
# nservers = [1,1,1]
# mus = [2.0,2.0,2.0]
# arrivals = [0.5,0.5,0.0]
# R = [0.0 0.0 1.0; 0.0 0.0 1.0; 0.0 0.0 0]
# jn = TruncatedOpenJacksonNetwork(array_of_sizes, nservers, mus, arrivals, R)
# pit = get_stat_dist(jn)
#
# jni = OpenJacksonNetwork(jn.nqueues, array_of_sizes, nservers, mus, arrivals, R)
# prod_array = get_stat_dist(jni)
# piv = product_form_to_vec(jn.vss, prod_array)
#
# tv(piv,pit)
#
# "
