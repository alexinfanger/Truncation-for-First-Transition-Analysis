using DataStructures

function lambda_func(t, p)
    if p["lambda_func_type"] == "sinusoidal"
        base = p["sin_base"]
        rel_amp = p["sin_rel_amp"]
        rate = p["sin_rate"]
        phase_shift = p["phase_shift"]
        return base .+ base.*rel_amp.* cos.(rate*t.+phase_shift)
    elseif p["lambda_func_type"] == "linear"
        return p["linear_slope"]*t + p["linear_intercept"]
    elseif p["lambda_func_type"] == "linear_vec"
        return vec_linear_interpolation(t, p["lambda_tvec"],p["lambda_yvec"])
    elseif p["lambda_func_type"] == "constant"
        return p["lambda_constant"]
    elseif p["lambda_func_type"] == "constant_vec"
        return vec_constant_interpolation(t, p["lambda_tvec"],p["lambda_yvec"])
    else
        print("Error: Not known lambda function type.")
    end
end

function mu_func(t,p)
    if p["mu_func_type"] == "linear_vec"
        return vec_linear_interpolation(t, p["mu_tvec"], p["mu_yvec"])
    elseif p["mu_func_type"] == "constant"
        return p["mu_constant"]
    elseif p["mu_func_type"] == "constant_vec"
        return vec_constant_interpolation(t, p["mu_tvec"], p["mu_yvec"])
    else
        print("Error: Not known mu function type.")
    end
end

function serv_func(t,p)
    if p["serv_func_type"] == "constant_vec"
        return vec_constant_interpolation(t, p["serv_tvec"], p["serv_yvec"])
    elseif p["serv_func_type"] == "constant"
        return p["s"]
    else
        print("Error: Not known server function type.")
    end

end


function Finite_Server_Q(t, p)
    S = Int(p["S"])
    s = serv_func(t,p)
    mu = mu_func(t,p)
    lambd = lambda_func(t,p)
    v = collect(1:S)
    br = ones(S)*lambd
    dr = min.(v, s)*mu .+ p["gamma"]*max.(0,v.-s)
    dg = zeros(S+1)
    dg[1] = -br[1]
    dg[2:end-1] = -br[2:end]-dr[1:end-1]
    dg[end] = -dr[end]
    dr = float(dr)
    Q = Tridiagonal(dr,dg,br)
    return Q
end

function Finite_Server_dudt(mu_t, p, t)
    if p["forwards"] == true
        return Finite_Server_Q(t, p)'*mu_t
    else
        return Finite_Server_Q(1440.0-t, p)*mu_t
    end
end

function Finite_Server_dudt!(du, mu_t, p, t)
    if p["forwards"] == true
        mul!(du, Finite_Server_Q(t, p)',mu_t)
    else
        mul!(du, Finite_Server_Q(1440.0-t, p),mu_t)
    end
end


function Finite_Server_Number_Waiting_r(p)
    S = Int(p["S"])
    s = Int(p["s"])
    a = zeros(s)
    b = range(0, stop=S-s)
    return vcat(a,b)
end


function Infinite_Server_closed_form(S, t, lambd, mu)
    alpha = lambd/mu * (1-exp.(-mu*t))
    a = collect(1:S)
    a = pushfirst!(a,1)
    b = cumsum(log.(a))
    r = collect(0:1.0:S)
    mu_t = exp.(-alpha .+ log.(alpha)*r .- b)
    return mu_t
end

function Infinite_Server_forward(S, t, lambd, mu, mu_0)
    prob = ODEProblem(Infinite_Server_dudt, mu_0, (0.0,t),[mu, lambd, S])
    solve(prob, reltol=1e-8,abstol=1e-8,alg_hints = "stiff")    # alg_hints = "stiff",
end

function test_function(a,b)
    if a < 1 && b==3
        c = 5
    else
        c=0
    end
    return c
end


# US Callcenter Model





# function S_Selector(rho_max, eps=1e-8)
#     S=50
#     while((1-rho_max)*(rho_max^S)>eps)
#         S = S+1
#     end
#     return S
# end


function generate_primatives(chain_dict)
    if p["type"] == "sinusoidal"
        rho_max = chain_dict["rho_bar"] *  (1 + chain_dict["sin_rel_amp"])
        chain_dict["S"] = S_Selector(rho_max)
        chain_dict["sin_base"] = chain_dict["s"]*chain_dict["mu"]*chain_dict["rho_bar"]
        if chain_dict["sin_type"] == "rev-up"
            chain_dict["phase_shift"] = pi
            chain_dict["sin_rate"] = (1*pi)/1440.0
            print("dude")
        elseif chain_dict["sin_type"] == "Linda_Green"
            chain_dict["phase_shift"] = 0
            chain_dict["sin_rate"] = (2*pi)/1440.0
        elseif chain_dict["sin_type"] == "rev-up-four-times"
            chain_dict["phase_shift"] = pi
            chain_dict["sin_rate"] = (8*pi)/1440.0
        end
    end
    return chain_dict
end


# reward functions

function prob_of_wait_r(S, num_servers)
    a = zeros(num_servers+1)
    b = ones(S-num_servers)
    return vcat(a,b)
end

function number_waiting_r(S, num_servers)
    a = zeros(num_servers)
    b = range(0, stop=S-num_servers)
    return vcat(a,b)
end


function MM1_Q(ld::AbstractFloat, mu::AbstractFloat, S::Integer)
    p = OrderedDict(); p["lambda_constant"] = ld; p["mu_constant"] = mu;
    p["s"] = 1; p["S"] = S
    p["lambda_func_type"] = "constant"; p["mu_func_type"] = "constant";
    p["serv_func_type"] = "constant";
    p["gamma"] = 0; n = p["S"] +1;

    return Finite_Server_Q(0,p)
end

function MM1_P(ld::AbstractFloat, mu::AbstractFloat, S::Integer)
    P, ld = uniformize(MM1_Q(ld, mu, S))
    return P;
end

function MMs_Q(ld::AbstractFloat, mu::AbstractFloat, S::Integer, s::Integer)
    p = OrderedDict(); p["lambda_constant"] = ld; p["mu_constant"] = mu;
    p["s"] = s; p["S"] = S
    p["lambda_func_type"] = "constant"; p["mu_func_type"] = "constant";
    p["serv_func_type"] = "constant";
    p["gamma"] = 0; n = p["S"] +1;

    return Finite_Server_Q(0,p)
end

function MMs_P(ld::AbstractFloat, mu::AbstractFloat, S::Integer, s::Integer)
    P, ld = uniformize(MMs_Q(ld, mu, S-1, s))
    return P;
end

struct MM1K <: FiniteMarkovChain
    ld::Number
    mu::Number
    S::Number
    P::AbstractArray
    vss::VectorStateSpace
end


function MM1K(ld::Number, mu::Number, S::Number)
    return MM1K(ld, mu, S, MM1_P(ld, mu, S), VectorStateSpace([S+1]))
end

function to_string(mm1k::MM1K)
    return @sprintf("M(ld=%0.2f)/M(mu=%0.2f)/1/%d",mm1k.ld, mm1k.mu, mm1k.S-2)
end

struct MMsK<:FiniteMarkovChain
    ld::Number
    mu::Number
    S::Number
    s::Number
    P::AbstractArray
end

function MMsK(ld::Number, mu::Number, S::Number, s::Number)
    return MMsK(ld, mu, S, s, MMs_P(ld, mu, S, s))
end

function to_string(mmsk::MMsK)
    return @sprintf("M(ld=%0.2f)/M(mu=%0.2f)/%d/%d",mmsk.ld, mmsk.mu, mmsk.s, mm1k.S-2)
end


struct MM1 <: InfiniteMarkovChain
    ld::Number
    mu::Number
    rv::Sampleable
end

function MM1(ld::Number, mu::Number)
    rv = Geometric(1-ld/mu)
    return MM1(ld, mu, rv)
end

function mean(mm1::MM1)
    return Statistics.mean(mm1.stat)
end

function to_string(mm1::MM1)
    return @sprintf("M(ld=%0.2f)/M(mu=%0.2f)/1",mm1k.ld, mm1k.mu)
end

function get_stat_dist(mm1::MM1, eps=1e-16)
    q = Integer(quantile(mm1.rv,1-eps))
    stat_dist = zeros(q)
    for i = 1:q
        stat_dist[i] = pdf(mm1.rv, i-1)
    end
    return stat_dist
end


struct MMinfty <: InfiniteMarkovChain
    ld::Number
    mu::Number
    stat::Sampleable
end

function MMinfty(ld::Number, mu::Number)
    return MMinfty(ld,mu, Poisson(ld/mu))
end

function get_stat_dist(mmi::MMinfty; eps=1e-16)
    q = Integer(quantile(mmi.stat,1-eps))
    stat_dist = zeros(q)
    for i = 1:q
        stat_dist[i] = pdf(mmi.stat, i-1)
    end
    return stat_dist
end


struct MMinftyK <: FiniteMarkovChain
    ld::Number
    mu::Number
    K::Number
    P::AbstractArray
end

function MMinftyK(ld::Number, mu::Number, K::Number)
    return MMsK(ld, mu, K, K)
end
