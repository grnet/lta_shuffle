"""
References are to the version included with this folder:
./as-implemented-Sep-2019/sub_shuffle.pdf
"""

from petlib.bn import Bn

from crs import CRS, POLYS_RHO, CRS_SM, CRS_PM, CRS_CON, CRS_UV, CRS_PKV

def mk_crs(gk, n, Chi):
    """
    Monolithic computation of CRS. Used to test multi-party computations
    """
    # Direct computation of involved parameters
    q = gk.q
    g1 = gk.g1
    g2 = gk.g2
    gt = gk.gt
    chi = Chi.chi                                                               # χ
    theta = Chi.theta                                                           # θ
    beta = Chi.beta                                                             # β
    betahat = Chi.betahat                                                       # β_hat
    rho = Chi.rho                                                               # ρ
    rhohat = Chi.rhohat                                                         # ρ_hat
    chipows = [chi ** i for i in range(0, 2 * n + 1)]                           # χ ^ i, 0<=i<=2n
    polyhats = [theta ** (2 * i) % q for i in range(1, n + 1)]                  # P_hat_i = θ ^ (2 * i), 1 <=i<= n
    els = generate_els(chipows, n, q)                                           # l_i, 1<=i<=n+1
    el_nplus1 = els[-1]                                                         # l_(n+1)
    poly_zero = el_nplus1 - Bn(1)                                               # P_0 = l_(n+1) - 1
    polys = [2 * els[i - 1] + el_nplus1 for i in range(1, n + 1)]               # P_i = 2 * l_i + l_(n+1), 1<=i<=n
    ques = [(poly + poly_zero) ** 2 - Bn(1) for poly in polys]                  # Q_i = (P_i + P_0) ^ 2 - 1
    beta_polys = [beta * poly + betahat * polyhat
            for (poly, polyhat) in zip(polys, polyhats)]                        # β * P_i + β_hat * P_hat_i, 1<=i<=n
    b_beta_polys = [beta * beta_poly for beta_poly in beta_polys]               # β ^ 2 * P_i  + β * β_hat * P_hat_i, 1<=i<=n

    # Convert to g_ζ elements, ζ = 1, 2
    g1_chi = chi * g1                                                           # [χ]_1
    g2_chi = chi * g2                                                           # [χ]_2
    g1_theta = theta * g1                                                       # [θ]_1
    g2_theta = theta * g2                                                       # [θ]_2
    g1_beta = beta * g1                                                         # [β]_1
    g2_beta = beta * g2                                                         # [β]_2
    g2_betasquare = (beta ** 2) * g2                                            # [β ^ 2]_2
    g1_betahat = betahat * g1                                                   # [β_hat]_1
    g2_betahat = betahat * g2                                                   # [β_hat]_2
    g1_beta_betahat = (beta * betahat) * g1                                     # [β * β_hat]_1
    g2_beta_betahat = (beta * betahat) * g2                                     # [β * β_hat]_2
    g1_rho = rho * g1                                                           # [ρ]_1
    g2_rho = rho * g2                                                           # [ρ]_2
    g1_rhohat = rhohat * g1                                                     # [ρ_hat]_1
    g1_beta_rho = g1_rho.mul(beta)                                              # [β ρ]_1
    g1_betasquare_rho = (beta ** 2 * rho) * g1                                  # [β ^ 2 ρ]_1
    g1_betahat_rhohat = g1_rhohat.mul(betahat)                                  # [β_hat ρ]_1
    g1_thetapows = [theta ** i * g1 for i in range(0, 2 * n + 1)]               # [θ ^ i]_1, 0<=i<=2n
    g1_poly_zero = poly_zero * g1                                               # [P_0]_1
    g2_poly_zero = poly_zero * g2                                               # [P_0]_2
    g1_polys = [_ * g1 for _ in polys]                                          # [P_i]_1 = 2 * [l_i]_1 + [l_(n+1)]_1, 1<=i<=n
    g2_polys = [_ * g2 for _ in polys]                                          # [P_i]_2 = 2 * [l_i]_2 + [l_(n+1)]_2, 1<=i<=n
    g1_poly_sum = sum(g1_polys, 0 * g1)                                         # [P_1 + ... + P_n]_1
    g2_poly_sum = sum(g2_polys, 0 * g2)                                         # [P_1 + ... + P_n]_2
    g1_polyhats = [_ * g1 for _ in polyhats]                                    # [P_hat_i]_1 = [θ ^ (2 * i)]_1, 1<=i<= n
    g1_polyhats_sum = sum(g1_polyhats, 0 * g1)                                  # [P_hat_1 + ... + P_hat_n]_1
    g1_beta_polys = [_ * g1 for _ in beta_polys]                                # [β * P_i + β_hat * P_hat_i]_1, 1<=i<=n
    g1_ques_rho = [(_ * g1).mul(rho.mod_inverse(q)) for _ in ques]              # [Q_i/ρ]_1 = [Q_i]/ρ, 1<=i<=n
    g1_b_beta_polys = [_ * g1 for _ in b_beta_polys]                            # [β ^ 2 * P_i  + β * β_hat * P_hat_i]_1, 1<=i<=n

    # Generate CRS fields
    g1_polys_rho = POLYS_RHO(g1_polys, g1_rho)
    g2_polys_rho = POLYS_RHO(g2_polys, g2_rho)
    crs_sm = CRS_SM(g1_beta_polys, g1_beta_rho,
        g1_betahat_rhohat, g2_beta, g2_betahat)
    crs_pm = CRS_PM(g1_poly_zero, g1_ques_rho, g1_poly_sum, g1_polyhats_sum,
        g2_poly_zero, g2_poly_sum)
    crs_con = CRS_CON(g1_polyhats, g1_rhohat)
    crs_uv = CRS_UV(g1_betasquare_rho, g1_beta_betahat, g1_b_beta_polys,
        g2_betasquare, g2_beta_betahat)
    crs_pkv = CRS_PKV(g1_beta, g1_betahat, g1_thetapows, g1_chi, g2_chi,
        g1_theta, g2_theta, g2_beta, g2_betahat)

    return CRS(gk, g1_polys_rho, g2_polys_rho, crs_sm, crs_pm, crs_con, crs_uv, crs_pkv)


def generate_els(chipows, n, q):
    """
    l_i, 1<=i<=n+1 (Theorem 7, C.2, pg. 26)
    """
    omega = (q - 1) // (n + 1)
    div_nplus1 = lambda _: Bn(n + 1).mod_inverse(q) * _
    els = []
    for i in range(1, n + 2):
        el_i = 0
        for j in range(0, n + 1):
            omega_ij = pow(omega, i * j, q)
            el_i += omega_ij.mod_inverse(q) * chipows[j]
        els.append(div_nplus1(el_i))
    return els
