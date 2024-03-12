
# Model
ld = 0.9
mu = 1.0
discount = 0.1
K = 1000

include("../src/models/queues/mms.jl")
Q = MM1_Q(ld, mu, K-1)
R, d_S = get_R_matrix(Q)
vss = VectorStateSpace([K])
Sidx = collect(1:vss.size)

e_S = ones(K)
lin_S = vss_abspow_vec(vss,1)
sq_S = vss_abspow_vec(vss,2)

# Define FTA Linear system
G = (d_S./(discount.+d_S)) .* R
f = sq_S./(discount.+d_S)

# Lyapunov Function
c = 1/discount
d = max(2/discount^2 *(ld-mu),0)
e = max((ld/discount)*(c+d), (ld-mu)*d/discount + c*(ld+mu)/discount)
g_S = c * sq_S .+ d * lin_S .+ e

@assert all(G*g_S .< g_S - f)

# Reference Solution
function MM1_discount_closed_form_solution(x)
  z1 = (discount + ld + mu - sqrt(-4*ld*mu+(discount+ld+mu)^2))/(2ld)
  a2 = 1/discount
  a1 = 2(ld-mu)/discount^2
  a0 = -(-discount*ld-2*ld^2-discount*mu+4*ld*mu-2mu^2)/discount^3
  c1 = (-discount*mu +  2ld*mu-2mu^2)/(discount^2*(discount+ld+0.5*
                (-discount-ld-mu+sqrt(-4*ld*mu+(discount+ld+mu)^2))))
  return c1*z1.^x.+a2*x.^2+a1.*x.+a0
end

# Reference Solution
u_reference = (I-G)\f
u_ref2 = (discount*I-Q)\sq_S # Should be an equivalent linear system
u_ref_closed_form = MM1_discount_closed_form_solution(0:vss.size-1)


# Reference solutions based on Q use a truncated and augmented system, so we
# expect them to be in disagreement for larger states. So we only compare them
# on the indices we will be approximating (Aidx_approx).

Aidx_approx = collect(1:81)
Acidx_approx = setdiff(Sidx, Aidx_approx)

println("Confirm various reference solutions are in agreement:")
println(norm(u_reference[Aidx_approx]-u_ref2[Aidx_approx],Inf))
println(norm(u_reference[Aidx_approx]-u_ref_closed_form[Aidx_approx],Inf))
println(norm(u_ref_closed_form[Aidx_approx]-u_ref2[Aidx_approx],Inf))

# Approximation


u_l = fta_lb(G, Aidx_approx, f)
u_u = fta_ub_complete(G, Aidx_approx, Acidx_approx, f, g_S)


p = plot(Aidx_approx.-1, log10.(u_reference[Aidx_approx]), 
  title="Bounds and Reference Solution", 
  label=L"\log_{10}(u)", 
  xlabel="State",
  legend=:topleft)
plot!(Aidx_approx.-1, log10.(u_u), label=L"\log_{10}(\tilde u)")
plot!(Aidx_approx.-1, log10.(u_l), label=L"$\log_{10}$(utilde($u$))")

# Paper
# dir = string("../../../../../../../../Apps/Overleaf/Truncation-FTA/Figures/")
# savefig(p, string(dir, "MM1",".png"))  
# dir = string("../../../../../../../../Apps/Overleaf/Truncation-FTA/ORLetters Version/reviewmode/Figures/")
# savefig(p, string(dir, "MM1",".png"))  

## Repo
dir = string("assets/")
savefig(p, string(dir, "MM1-reference-and-bounds",".png"))  


rel_err_gap_both = (u_u-u_l)./(u_ref_closed_form[Aidx_approx])
rel_err_gap_upper = (u_u-u_ref_closed_form[Aidx_approx])./(u_ref_closed_form[Aidx_approx])
rel_err_gap_lower = (u_ref_closed_form[Aidx_approx]-u_l)./(u_ref_closed_form[Aidx_approx])


p = plot(Aidx_approx.-1,log10.(rel_err_gap_lower), 
  title="Relative Error Gap", 
  label=L"\log_{10}(\delta)", 
  xlabel="State",
  legend=:topleft)

plot!(Aidx_approx.-1, log10.(rel_err_gap_upper), label=L"\log_{10}(\Delta)")

# # Paper
# dir = string("../../../../../../../../Apps/Overleaf/Truncation-FTA/Figures/")
# savefig(p, string(dir, "MM1-rel-err-gap",".png"))  
# dir = string(
#   "../../../../../../../../Apps/Overleaf/Truncation-FTA/ORLetters Version/reviewmode/Figures/")
# savefig(p, string(dir, "MM1-rel-err-gap",".png"))  

## Repo
dir = string("assets/")
savefig(p, string(dir, "MM1-rel-err-gap",".png"))  


# Convergence
tvals = [80, 100, 120, 140].+1

rel_err_gaps_lower = zeros(length(tvals))
rel_err_gaps_upper = zeros(length(tvals))

A_0 = get_linear_sublevel_set(vss, tvals[1])

for (i,t_i) in enumerate(tvals)
    
  Aidx_i = collect(1:t_i)
  Acidx_i = setdiff(Sidx, Aidx_i)

  u_l_i = fta_lb(G, Aidx_i, f)
  u_u_i = fta_ub_complete(G, Aidx_i, Acidx_i, f, g_S)

  println("Test: Confirming reference solutions are in agreement with simple truncation:")
  println(norm(u_reference[Aidx_i]-u_ref2[Aidx_i],Inf))
  println(norm(u_reference[Aidx_i]-u_ref_closed_form[Aidx_i],Inf))
  println(norm(u_ref_closed_form[Aidx_i]-u_ref2[Aidx_i],Inf))

  rel_err_gap_lower = Aidx_to_Sidx(vss, (u_ref_closed_form[Aidx_i]-u_l_i)./
                                   u_ref_closed_form[Aidx_i], Aidx_i)
  rel_err_gap_upper = Aidx_to_Sidx(vss, (u_u_i-u_ref_closed_form[Aidx_i])./
                                   u_ref_closed_form[Aidx_i], Aidx_i)

  rel_err_gaps_lower[i] = norm(rel_err_gap_lower[A_0], Inf)
  rel_err_gaps_upper[i] = norm(rel_err_gap_upper[A_0], Inf)
end

p_convergence = plot(tvals.-1, log10.(rel_err_gaps_lower), 
                     label=L"\log_{10}\left(\max_{x\in A_0}\ \delta^i(x)\right)",
                     marker = :circle, xlabel=L"$t_i$ for $A_i=\{x: x \leq t_i\}$")
plot!(tvals.-1, log10.(rel_err_gaps_upper),
                    label=L"\log_{10}\left(\max_{x\in A_0}\ \Delta^i(x)\right)",
                    marker = :circle)

# ## Paper
# dir = string("../../../../../../../../Apps/Overleaf/Truncation-FTA/Figures/")
# savefig(p2, string(dir, "MM1-relative-err-convergence",".png"))  
# dir = string(
#  "../../../../../../../../Apps/Overleaf/Truncation-FTA/ORLetters Version/reviewmode/Figures/")
# savefig(p2, string(dir, "MM1-relative-err-convergence",".png"))  

## Repo
dir = string("assets/")
savefig(p_convergence, string(dir, "MM1-relative-err-convergence",".png"))


