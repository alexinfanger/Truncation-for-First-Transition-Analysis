# Model
include("../src/models/queues/JacksonNetwork/JacksonNetwork.jl")
ld = 5.0 
mu = 1.0 

n_universe = 300
jnt = TruncatedOpenJacksonNetwork([n_universe, n_universe],[n_universe, n_universe],[mu,mu],
                                  [ld,0.0],[0 1.0; 0 0 ])
Q = jnt.Q
R, d_S = get_R_matrix(Q)
vss = jnt.vss
Sidx = collect(1:vss.size)
zero_idx = vec_to_ind(vss, [1,1])
e_S = ones(vss.size)

# Define FTA Linear system
G = R
f = e_S./d_S

# Lyapunov function
c = 1
g_func(x) = 2c*(x[1]-1.0)+c*(x[2]-1.0)
g_S = vss_func_to_vec(vss, g_func)

closed_form_Qg_func(x) = c*(2ld-mu*((x[1]-1.0)+(x[2]-1.0)))
closed_form_Qg_S = vss_func_to_vec(vss, closed_form_Qg_func)

diff_S = Q*g_S - closed_form_Qg_S
println("Confirm closed form Qg is in agreement Qg: ",
        all(diff_S[get_linear_sublevel_set(vss, n_universe-2)].== 0))


# Compute Reference Solution
t_theory = Integer((2ld+1)/mu)
println("t_theory: ", t_theory)
t_reference = 60
println("t_reference: ", t_reference)
@assert t_reference > t_theory

Aidx_reference = get_linear_sublevel_set(vss, t_reference)
Acidx_reference = setdiff(Sidx, Aidx_reference)
Aidx_reference, Acidx_reference = drop_states(Aidx_reference, Acidx_reference, [zero_idx])

@assert is_Lyapunov_function_rate_full_multiply(Q, Acidx_reference, e_S, g_S) 

u_l_reference = fta_lb(G, Aidx_reference, f)
u_u_reference = fta_ub_incomplete(G, Aidx_reference, Acidx_reference, f, g_S)

u_l_reference_S = Aidx_to_Sidx(vss, u_l_reference, Aidx_reference)
u_u_reference_S = Aidx_to_Sidx(vss, u_u_reference, Aidx_reference)

# End Compute Reference Solution

# Check correctness by comparing to M/G/\infty queue interpretation
busy_period = (exp(2*ld/mu)-1)/ld
println(busy_period)
println(u_l_reference_S[vec_to_ind(vss,[2,1])])
println("|True value of u(1,0)-Reference lb of u(1,0)|:     ",
         abs(busy_period-u_l_reference_S[vec_to_ind(vss,[2,1])]))
println("Relative Error:     ", 
        abs(busy_period-u_l_reference_S[vec_to_ind(vss,[2,1])])/abs(busy_period))

println("|True value of u(1,0)-Reference ub of u(1,0)|:     ", 
        abs(busy_period-u_u_reference_S[vec_to_ind(vss,[2,1])]))
println("Relative Error:     ", 
        abs(busy_period-u_u_reference_S[vec_to_ind(vss,[2,1])])/abs(busy_period))

# Compute approximation and compute error
t_approx = 40

Aidx_approx = get_linear_sublevel_set(vss, t_approx)
Acidx_approx = setdiff(Sidx, Aidx_approx)
Aidx_approx, Acidx_approx = drop_states(Aidx_approx, Acidx_approx, [zero_idx])

println("Error Guarantee of Reference Solution on A_approx:   ", 
        norm(u_u_reference_S[Aidx_approx]-u_l_reference_S[Aidx_approx], Inf))
println("Relative Error Guarantee of Reference Solution on A_approx:   ", 
        norm((u_u_reference_S[Aidx_approx]-u_l_reference_S[Aidx_approx])./
                u_l_reference_S[Aidx_approx], Inf))

@assert is_Lyapunov_function_rate_full_multiply(Q, Acidx_approx, e_S, g_S) 

u_l_approx = fta_lb(G, Aidx_approx, f)
u_u_approx = fta_ub_incomplete(G, Aidx_approx, Acidx_approx, f, g_S)

rel_err_gap_both = Aidx_to_Sidx(vss,(u_u_approx-u_l_approx)./
                                    u_u_reference_S[Aidx_approx], Aidx_approx)
rel_err_gap_upper = Aidx_to_Sidx(vss, (u_u_approx-u_u_reference_S[Aidx_approx])./
                                        u_u_reference_S[Aidx_approx], Aidx_approx)
rel_err_gap_lower = Aidx_to_Sidx(vss, (u_u_reference_S[Aidx_approx]-u_l_approx)./
                                        u_u_reference_S[Aidx_approx], Aidx_approx)

# End compute approximation and compute error

int_t_approx = Integer(t_approx)
p_both = heatmap(0:int_t_approx-1, 0:int_t_approx-1, log10.(vec_to_twodmatrix(vss,rel_err_gap_both)[1:int_t_approx,1:int_t_approx]'),
            xlabel=L"x_1",
            ylabel=L"x_2",
            title="Relative Error Gap",
            xticks = 10:10:30,
            ytick = 10:10:30,
            topmargin=6mm
        )
p_upper = heatmap(0:int_t_approx-1, 0:int_t_approx-1,log10.(vec_to_twodmatrix(vss,rel_err_gap_upper)[1:int_t_approx,1:int_t_approx]'),
            xlabel=L"x_1",
            ylabel=L"x_2",
            title=L"\log_{10}(\Delta)",
            xticks = 10:10:30,
            ytick = 10:10:30,
            topmargin=6mm
        )
p_lower = heatmap(0:int_t_approx-1, 0:int_t_approx-1, log10.(vec_to_twodmatrix(vss,rel_err_gap_lower)[1:int_t_approx,1:int_t_approx]'),
            xlabel=L"x_1",
            ylabel=L"x_2",
            title=L"\log_{10}(\delta)",
            xticks = 10:10:30,
            ytick = 10:10:30,
            topmargin=6mm
        )

# Paper
# dir = string("../../../../../../../../Apps/Overleaf/Truncation-FTA/Figures/")
# savefig(p_both, string(dir, "Tandem-Both",".png"))  
# savefig(p_upper, string(dir, "Tandem-Upper",".png")) 
# savefig(p_lower, string(dir, "Tandem-Lower",".png")) 
# dir = string(
#  "../../../../../../../../Apps/Overleaf/Truncation-FTA/ORLetters Version/reviewmode/Figures/")
# savefig(p_both, string(dir, "Tandem-Both",".png"))  
# savefig(p_upper, string(dir, "Tandem-Upper",".png")) 
# savefig(p_lower, string(dir, "Tandem-Lower",".png")) 

# Repo
dir = string("assets/")
savefig(p_both, string(dir, "Tandem-Both",".png"))
savefig(p_upper, string(dir, "Tandem-Upper",".png"))
savefig(p_lower, string(dir, "Tandem-Lower",".png"))


# Convergence
minval = 25
tvals = [minval+5*i for i in 0:3]
rel_err_gaps_upper = zeros(length(tvals))
rel_err_gaps_lower = zeros(length(tvals))
A_0 = get_linear_sublevel_set(vss, tvals[1])

for (i,t_i) in enumerate(tvals)

  Aidx_i = get_linear_sublevel_set(vss, t_i)
  Acidx_i = setdiff(Sidx, Aidx_i)
  Aidx_i, Acidx_i = drop_states(Aidx_i, Acidx_i, [zero_idx])
  
  println("----------------------------------------")
  println("t_i:   ", t_i, "  ")

  println("Error Guarantee of Reference Solution on A:   ",
            norm(u_u_reference_S[Aidx_i]-u_l_reference_S[Aidx_i], Inf))
  println("Relative Error Guarantee of Reference Solution on A:   ", 
            norm((u_u_reference_S[Aidx_i]-u_l_reference_S[Aidx_i])./
                  u_u_reference_S[Aidx_i], Inf))
  @assert is_Lyapunov_function_rate_full_multiply(Q, Acidx_i, e_S, g_S)

  u_l_i = fta_lb(G, Aidx_i, f)
  u_u_i = fta_ub_incomplete(G, Aidx_i, Acidx_i, f, g_S)

  rel_err_gap_lower = Aidx_to_Sidx(vss, (u_u_reference_S[Aidx_i]-u_l_i)./
                                          u_u_reference_S[Aidx_i], Aidx_i)
  rel_err_gap_upper = Aidx_to_Sidx(vss, (u_u_i-u_u_reference_S[Aidx_i])./
                                          u_u_reference_S[Aidx_i], Aidx_i)

  rel_err_gaps_lower[i] = norm(rel_err_gap_lower[A_0], Inf)
  rel_err_gaps_upper[i] = norm(rel_err_gap_upper[A_0], Inf)
end

p_tandem_convergence = plot(tvals, log10.(rel_err_gaps_lower), 
                            label=L"\log_{10}\left(\max_{x\in A_0}\ \delta^i(x)\right)",
                            xlabel=L"$t_i$ for $A_i=\{x: x_1+x_2 \leq t_i\}$",
                            legend=:bottomleft, marker = :circle)
xticks!(tvals)
plot!(tvals, log10.(rel_err_gaps_upper), 
      label=L"\log_{10}\left(\max_{x\in A_0}\ \Delta^i(x)\right)", marker = :circle)

# Paper
# dir = string("../../../../../../../../Apps/Overleaf/Truncation-FTA/Figures/")
# savefig(p_tandem_convergence, string(dir, "TandemQueue-relative-err-convergence",".png"))
# dir = string(
#  "../../../../../../../../Apps/Overleaf/Truncation-FTA/ORLetters Version/reviewmode/Figures/")
# savefig(p_tandem_convergence, string(dir, "TandemQueue-relative-err-convergence",".png"))

# Repo
dir = string("assets/")
savefig(p_tandem_convergence, string(dir, "TandemQueue-relative-err-convergence",".png"))


