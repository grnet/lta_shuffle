"""
References are to the version included with this folder:
./as-implemented-Sep-2019/sub_shuffle.pdf
"""

from functools import reduce
from operator import mul
from bplib.bp import GTElem
from prover import extract_proof

from crs import (extract_gk, extract_pm, extract_con, extract_polys_rhos,
    extract_uv, extract_pkv)

def verify(n, crs, pk_e, ciphers, proof):
    """
    Fig. 3, pg. 14
    """
    pk_ok = step_0(n, crs, proof)
    assert pk_ok

    cprimes, g1_alpha_hats, g1_s, g2_N, pi_per = step_1(n, crs, proof)
    g1_alpha_hats = step_2(n, crs, g1_alpha_hats)
    perm_ok = step_3(n, crs, g1_alpha_hats, pi_per)
    cons_ok = step_4(n, crs, pk_e, ciphers, cprimes, g1_alpha_hats, g1_s, g2_N)

    return perm_ok, cons_ok

def step_0(n, crs, proof):
    return V_pk(n, crs)

def step_1(n, crs, proof):
    cprimes, g1_alpha_hats, g1_s, g2_N, pi_per = extract_proof(proof)
    return cprimes, g1_alpha_hats, g1_s, g2_N, pi_per

def step_2(n, crs, g1_alpha_hats):
    _, _, g1, _, _, _ = extract_gk(crs)
    _, _, _, g1_polyhats_sum, _, _ = extract_pm(crs)

    alpha_hat_n = g1_polyhats_sum - sum(g1_alpha_hats, 0 * g1)
    g1_alpha_hats.append(alpha_hat_n)
    return g1_alpha_hats

def step_3(n, crs, g1_alpha_hats, pi_per):
    return V_per(n, crs, g1_alpha_hats, pi_per)

def step_4(n, crs, pk_e, ciphers, cprimes, g1_alpha_hats, g1_s, g2_N):
    _, _, g1, _, _, pair = extract_gk(crs)
    _, _, _, g1_polyhats_sum, _, _ = extract_pm(crs)
    g1_polyhats, _ = extract_con(crs)

    left_1 = multi_product(pair, g1_polyhats, cprimes)
    left_2 = multi_product(pair, g1_alpha_hats, ciphers)
    left_side = (left_1[0] * left_2[0].inv(), left_1[1] * left_2[1].inv())

    right_1 = (pair(g1_s, pk_e[0]), pair(g1_s, pk_e[1]))
    right_2 = (pair(g1, g2_N[0]), pair(g1, g2_N[1]))
    right_side = (right_1[0] * right_2[0].inv(), right_1[1] * right_2[1].inv())

    return left_side == right_side


def V_pk(n, crs):
    """
    Fig. 8, pg. 20
    """
    _, _, g1, g2, gt, pair = extract_gk(crs)
    g1_polys, g2_polys, g1_rho, g2_rho = extract_polys_rhos(crs)
    g1_betasquare_rho, g1_beta_betahat, g1_b_beta_polys, \
        g2_betasquare, g2_beta_betahat = extract_uv(crs)
    g1_polyhats, _ = extract_con(crs)
    g1_poly_zero, g1_ques_rho, _, _, g2_poly_zero, _ = extract_pm(crs)

    # step 1
    g1_beta, g1_betahat, g1_thetapows, g1_chi, g2_chi, \
        g1_theta, g2_theta, g2_beta, g2_betahat = extract_pkv(crs)

    # steps 2-11
    return all((
        g1_rho != 0 * g1,
        all((pair(g1__, g2) == pair(g1, g2__)
                for (g1__, g2__) in ((g1_chi, g2_chi),
                                     (g1_theta, g2_theta),
                                     (g1_beta, g2_beta),
                                     (g1_betahat, g2_betahat),
                                     (g1_rho, g2_rho)))),
        gt == pair(g1, g2),
        all((pair(g1_thetapows[i], g2) == pair(g1_thetapows[i - 1], g2_theta)
                for i in range(2, 2 * n + 1))),
        pair(g1, g2_betasquare) == pair(g1_beta, g2_beta),
        pair(g1_betasquare_rho, g2) == pair(g1_rho, g2_betasquare),
        pair(g1_beta_betahat, g2) == pair(g1_beta, g2_betahat),
        pair(g1, g2_beta_betahat) == pair(g1_beta_betahat, g2),
        pair(g1, g2_poly_zero) == pair(g1_poly_zero, g2),
        all(((
            all((pair(g1, g2_polys[i - 1]) == pair(g1_polys[i - 1], g2),
                pair(g1_b_beta_polys[i - 1], g2) == \
                    pair(g1_polys[i - 1], g2_betasquare) * \
                    pair(g1_polyhats[i - 1], g2_beta_betahat),
                pair(g1_ques_rho[i - 1], g2_rho) == \
                    pair(g1_polys[i - 1] + g1_poly_zero, g2_polys[i - 1] + \
                    g2_poly_zero).mul(gt.inv()))))
                for i in range(1, n + 1)))
        ))

def V_per(n, crs, g1_alpha_hats, pi_per):
    """
    Permutation argument, Fig. 6, pg. 17
    """
    for pi_uvi, g1_alpha_hat in zip(pi_per, g1_alpha_hats):
        if not V_uv(n, crs, g1_alpha_hat, pi_uvi):
            return False
    return True

def V_uv(n, crs, g1_alpha_hat, pi):
    """
    New unit verctor argument, Fig. 5, pg 16
    """
    _, q, g1, g2, gt, pair = extract_gk(crs)
    _, _, _, g2_rho = extract_polys_rhos(crs)
    _, _, _, g2_betasquare, g2_beta_betahat = extract_uv(crs)
    g1_poly_zero, _, _, _, g2_poly_zero, _ = extract_pm(crs)

    g1_d, g1_a, g2_b, g1_e = pi
    A = q.random()
    if pair(g1_d, g2) != vector_product(pair, (g1_a, g1_alpha_hat), (g2_betasquare, g2_beta_betahat)):
        return False
    left_side = pair(g1_a + A * g1 + g1_poly_zero, g2_b - A * g2 + g2_poly_zero)
    right_side = pair(g1_e, g2_rho) * gt ** ((1 - A ** 2) % q)
    return left_side == right_side

def vector_product(pair, g1_vector, g2_vector):
    """
    Applies componentwise bilinear pairing and multiplies the results:

    (a_1, ..., b_1), (c_2, ..., d_2) |---> pair(a_1, c_2) * ... * pair(b_1, d_2)
    """
    # # Equivalent to:
    #
    # result = one
    # for (g1_, g2_) in zip(g1_vector, g2_vector):
    #     result *= pair(g1_, g2_)
    # return result
    #
    return reduce(mul, map(lambda _: pair(_[0], _[1]) , zip(g1_vector, g2_vector)))

def multi_product(pair, g1_vector, g2_vectors):
    """
    Returns in respective order the vector products of the provided g1_vector with
    each one of the g2_vectors comprising the provided collection
    """
    return tuple(vector_product(pair, g1_vector, g2_vector) for g2_vector in zip(*g2_vectors))
