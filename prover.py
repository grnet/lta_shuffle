"""
References are to the version included with this folder:
./as-implemented-Sep-2019/sub_shuffle.pdf
"""

from functools import reduce
from collections import namedtuple
from bplib.bp import Bn

from crs import extract_gk, extract_con, extract_pm, extract_polys_rhos, extract_uv

PROOF = namedtuple('PROOF', ['cprimes', 'pi_sh'])

def extract_proof(proof):
    cprimes = proof.cprimes
    g1_alpha_hats, g1_s, g2_N, pi_per = proof.pi_sh
    return cprimes, g1_alpha_hats, g1_s, g2_N, pi_per

def prove(n, pk_e, crs, ciphers, sigma, t_randoms):
    """
    Fig. 3, pg. 14
    """
    r_hats, g1_alpha_hats, sigma_inv = step_1(n, crs, sigma)
    rn_hat = step_2(n, crs, r_hats)
    pi_per = step_3(n, crs, g1_alpha_hats, sigma_inv, r_hats)
    r_hat, g1_s = step_4(n, crs, t_randoms)
    g2_tprimes = step_5(n, crs, pk_e, t_randoms)
    g2_N = step_6(n, crs, pk_e, ciphers, r_hats, r_hat)
    cprimes = step_7(n, crs, ciphers, sigma, g2_tprimes)
    pi_con = step_8(n, crs, g1_s, g2_N)

    return PROOF(cprimes, pi_sh=(g1_alpha_hats[:-1], g1_s, g2_N, pi_per))

def step_1(n, crs, sigma):
    _, q, g1, _, _, _ = extract_gk(crs)
    g1_polyhats, _ = extract_con(crs)

    r_hats = [q.random() for _ in range(1, n)]
    g1_alpha_hats = []
    sigma_inv = inverse_perm(sigma)
    for (perm_i, ri_hat) in zip(sigma_inv, r_hats):
        g1_alpha_hats.append(g1_polyhats[perm_i] + ri_hat * g1)
    return r_hats, g1_alpha_hats, sigma_inv

def step_2(n, crs, r_hats):
    _, q, _, _, _, _ = extract_gk(crs)

    rn_hat = (- sum(r_hats)).mod(q)
    r_hats.append(rn_hat)
    return rn_hat

def step_3(n, crs, g1_alpha_hats, sigma_inv, r_hats):
    return P_per(n, crs, g1_alpha_hats, sigma_inv, r_hats)

def step_4(n, crs, t_randoms):
    _, q, g1, _, _, _ = extract_gk(crs)
    g1_polyhats, _ = extract_con(crs)

    r_hat = q.random()
    g1_s = sum([t * p for (t, p) in zip(t_randoms, g1_polyhats)], r_hat * g1)
    return r_hat, g1_s

def step_5(n, crs, pk_e, t_randoms):
    g2_tprimes = [(_ * pk_e[0], _ * pk_e[1]) for _ in t_randoms]
    return g2_tprimes

def step_6(n, crs, pk_e, ciphers, r_hats, r_hat):
    aux = tuple((r * c[0], r * c[1]) for (r, c) in zip(r_hats, ciphers))
    aux += ((r_hat * pk_e[0], r_hat * pk_e[1]),)
    g2_N = reduce(lambda x, y: (x[0] + y[0], x[1] + y[1]), aux)
    return g2_N

def step_7(n, crs, ciphers, sigma, g2_tprimes):
    aux = lambda x, y: (x[0] + y[0], x[1] + y[1])
    cprimes = [aux(ciphers[index], g2_tprime)
        for index, g2_tprime in zip(sigma, g2_tprimes)]
    return cprimes

def step_8(n, crs, g1_s, g2_N):
    pi_con = (g1_s, g2_N)
    return pi_con

def P_per(n, crs, g1_alpha_hats, sigma_inv, r_hats):
    """
    Permutation argument, Fig. 6, pg. 17
    """
    _, _, g1, _, _, _ = extract_gk(crs)
    _, _, _, g1_polyhats_sum, _, _ = extract_pm(crs)

    g1_alpha_n = g1_polyhats_sum - sum(g1_alpha_hats, 0 * g1)
    g1_alpha_hats.append(g1_alpha_n)
    pi_per = []
    for i in range(0, n):
        pi_uvi = P_uv(n, crs, g1_alpha_hats[i], sigma_inv[i], r_hats[i])
        pi_per.append(pi_uvi)
    return pi_per

def P_uv(n, crs, g1_alpha_hat, index, r_hat):
    """
    New unit verctor argument, Fig. 5, pg 16
    """
    _, q, g1, _, _, _ = extract_gk(crs)
    g1_polys, g2_polys, g1_rho, g2_rho = extract_polys_rhos(crs)
    g1_poly_zero, g1_ques_rho, _, _, _, _ = extract_pm(crs)
    g1_betasquare_rho, g1_beta_betahat, g1_b_beta_polys, _, _ = extract_uv(crs)

    # step 1.
    r = q.random()
    g1_rprime = r * g1_rho
    g1_d = g1_b_beta_polys[index] + r * g1_betasquare_rho + r_hat * g1_beta_betahat

    # step 2.
    g1_a = g1_polys[index] + g1_rprime
    g2_b = g2_polys[index] + r * g2_rho

    # step 3.
    g1_e = r * (2 * (g1_a + g1_poly_zero) - g1_rprime) + g1_ques_rho[index]

    # step 4.
    return (g1_d, g1_a, g2_b, g1_e)

def inverse_perm(sigma):
    sigma_inv = [None] * len(sigma)
    for index, value in enumerate(sigma):
        sigma_inv[value] = index
    return sigma_inv
