"""
References are to the version included with this folder:
./as-implemented-Sep-2019/sub_shuffle.pdf
"""

import datetime
from math import log
from collections import namedtuple
from petlib.bn import Bn

from md import MPCMD
from utils import MD_Elems


Chi_share = namedtuple('Chi_share', ['chi',
                                     'theta',
                                     'beta',
                                     'betahat',
                                     'rho',
                                     'rhohat'])

def mk_Chi_share(q):
    randstar = lambda: 1 + (q - 1).random()

    chi = randstar()
    theta = randstar()
    beta = randstar()
    betahat = randstar()
    rho = randstar()
    rhohat = randstar()

    return Chi_share(chi, theta, beta, betahat, rho, rhohat)


Chi = namedtuple('Chi', ['g1_chi', 'g2_chi',
                         'g1_theta',
                         'g1_beta', 'g2_beta',
                         'g1_betahat', 'g2_betahat',
                         'g1_rho', 'g2_rho',
                         'g1_rhohat'])

CRS = namedtuple('CRS', ['gk',
                         'g1_polys_rho',
                         'g2_polys_rho',
                         'crs_sm',
                         'crs_pm',
                         'crs_con',
                         'crs_uv',
                         'crs_pkv'])

POLYS_RHO = namedtuple('POLYS_RHO', ['polys',
                                     'rho'])
CRS_SM = namedtuple('CRS_SM', ['g1_beta_polys',
                               'g1_beta_rho',
                               'g1_betahat_rhohat',
                               'g2_beta',
                               'g2_betahat'])
CRS_PM = namedtuple('CRS_PM', ['g1_poly_zero',
                               'g1_ques_rho',
                               'g1_poly_sum',
                               'g1_polyhats_sum',
                               'g2_poly_zero',
                               'g2_poly_sum'])
CRS_CON = namedtuple('CRS_CON', ['g1_polyhats',
                                 'g1_rhohat'])
CRS_UV = namedtuple('CRS_UV', ['g1_betasquare_rho',
                               'g1_beta_betahat',
                               'g1_b_beta_polys',
                               'g2_betasquare',
                               'g2_beta_betahat'])
CRS_PKV = namedtuple('CRS_PKV', ['g1_beta',
                                 'g1_betahat',
                                 'g1_thetapows',
                                 'g1_chi',
                                 'g2_chi',
                                 'g1_theta',
                                 'g2_theta',
                                 'g2_beta',
                                 'g2_betahat'])
# CRS extraction

def extract_gk(crs):
    """
    Extracts the underlying bilinear group G, its order q, the generators g1, g2
    of its first resp. second component, their pairing gt and the operation of
    bilinear pairing itself
    """
    gk = crs.gk
    G, q, g1, g2, gt, pair = gk.G, gk.q, gk.g1, gk.g2, gk.gt, gk.pair
    return G, q, g1, g2, gt, pair

def extract_polys_rhos(crs):
    """
    Extract parameters

    [P_(1<=i<=n)]_1, [P_(1<=i<=n)]_2, [ρ]_1, [ρ]_2
    """
    g1_polys = crs.g1_polys_rho.polys
    g2_polys = crs.g2_polys_rho.polys
    g1_rho = crs.g1_polys_rho.rho
    g2_rho = crs.g2_polys_rho.rho

    return g1_polys, g2_polys, g1_rho, g2_rho

def extract_sm(crs):
    """
    Extract crs_sm parameters

    [(β P_i + β_hat P_hat_i)_(1<=i<=n)]_1, [β ρ]_1, [β_hat ρ_hat]_1,
    [β]_2, [β_hat]_2
    """
    crs_sm = crs.crs_sm

    g1_beta_polys = crs_sm.g1_beta_polys
    g1_beta_rho = crs_sm.g1_beta_rho
    g1_betahat_rhohat = crs_sm.g1_betahat_rhohat
    g2_beta = crs_sm.g2_beta
    g2_betahat = crs_sm.g2_betahat

    return (g1_beta_polys, g1_beta_rho, g1_betahat_rhohat, g2_beta, g2_betahat)

def extract_pm(crs):
    """
    Extract crs_pm parameters

    [P_0]_1, [Q_(1<=i<=n)/ρ]_1, [P_1 + ... + P_n]_1,
    [P_hat_1 + ... + P_hat_n]_1, [P_0]_2, [P_1 + ... + P_n]_2
    """
    crs_pm = crs.crs_pm

    g1_poly_zero = crs_pm.g1_poly_zero
    g1_ques_rho = crs_pm.g1_ques_rho
    g1_poly_sum = crs_pm.g1_poly_sum
    g1_polyhats_sum = crs_pm.g1_polyhats_sum
    g2_poly_zero = crs_pm.g2_poly_zero
    g2_poly_sum = crs_pm.g2_poly_sum

    return (g1_poly_zero, g1_ques_rho, g1_poly_sum, g1_polyhats_sum,
        g2_poly_zero, g2_poly_sum)

def extract_con(crs):
    """
    Extract crs_con parameters

    [P_hat_(1<=i<n)]_1, [ρ_hat]_1
    """
    crs_con = crs.crs_con

    g1_polyhats = crs_con.g1_polyhats
    g1_rhohat = crs_con.g1_rhohat

    return (g1_polyhats, g1_rhohat)

def extract_uv(crs):
    """
    Extract selected crs_uv parameters

    [β^2 ρ]_1, [β β_hat]_1, [β^2 P_i + β β_hat P_i]_1, 1<=i<=n, [β^2]_2, [β β_hat]_2

    (cf. Fig 3, pg 14)
    """
    crs_uv = crs.crs_uv

    g1_betasquare_rho = crs_uv.g1_betasquare_rho
    g1_beta_betahat = crs_uv.g1_beta_betahat
    g1_b_beta_polys = crs_uv.g1_b_beta_polys
    g2_betasquare = crs_uv.g2_betasquare
    g2_beta_betahat = crs_uv.g2_beta_betahat

    return (g1_betasquare_rho, g1_beta_betahat, g1_b_beta_polys,
        g2_betasquare, g2_beta_betahat)

def extract_pkv(crs):
    """
    Extract selected crs_pkv parameters

    ...

    (cf. Fig 3, pg. 14)
    """
    crs_pkv = crs.crs_pkv

    g1_beta = crs_pkv.g1_beta
    g1_betahat = crs_pkv.g1_betahat
    g1_thetapows = crs_pkv.g1_thetapows
    g1_chi = crs_pkv.g1_chi
    g2_chi = crs_pkv.g2_chi
    g1_theta = crs_pkv.g1_theta
    g2_theta = crs_pkv.g2_theta
    g2_beta = crs_pkv.g2_beta
    g2_betahat = crs_pkv.g2_betahat

    return (g1_beta, g1_betahat, g1_thetapows, g1_chi, g2_chi,
        g1_theta, g2_theta, g2_beta, g2_betahat)

# CRS generation

def mk_crs(gk, n, Chi_share, config):
    """
    CRS computation as a circuit, Fig. 9, App. A, pg. 25
    """
    peers = config['peers']
    name = config['name']
    spec = config['spec']
    mpc = MPCMD(peers, name, spec)

    # Run layers

    Chi = layer_C1_1(gk, Chi_share, mpc)

    g1_chipows, g1_thetapows = layer_C1_2(gk, Chi_share, Chi, n, mpc)

    g2_theta, g1_polyhats = layer_C1_3(gk, Chi_share, g1_thetapows, n, mpc)

    g2_chipows, g1_beta_chipows, g1_betahat_polyhats = \
        layer_C1_4(gk, Chi_share, Chi, g1_chipows[:n + 1], g1_polyhats, n, mpc)

    (g2_betasquare, g1_beta_rho, g1_betasquare_rho, g1_beta_betahat,
        g2_beta_betahat, g1_betahat_rhohat) = layer_C1_5(gk, Chi_share, Chi, mpc)

    g1_els, g2_els, g1_el_squares, g1_beta_els, g1_elinplus1_prods = \
        layer_L1_1(gk, g1_chipows, g1_beta_chipows, g2_chipows, n)

    g1_el_nplus1 = g1_els[-1]
    g2_el_nplus1 = g2_els[-1]
    g1_poly_zero, g2_poly_zero = layer_L1_2(gk, g1_el_nplus1, g2_el_nplus1)

    g1_el_nplus1_sq = g1_el_squares[-1]
    g1_beta_el_nplus1 = g1_beta_els[-1]
    g1_polys, g2_polys, g1_ques, g1_beta_polys = \
    layer_L1_3(gk, n, Chi_share, g1_els, g1_el_nplus1, g2_els, g1_el_squares,
        g1_beta_el_nplus1, g1_el_nplus1_sq, g1_elinplus1_prods, g1_betahat_polyhats, g1_beta_els)

    g1_poly_sum, g2_poly_sum, g1_polyhats_sum = layer_L1_4(gk, g1_polys, g2_polys, g1_polyhats)

    g1_ques_rho, g1_b_beta_polys = layer_C2(gk, Chi_share, g1_ques, g1_beta_polys, mpc)

    mpc.quit()
    end = datetime.datetime.now()
    print()
    print(' * MPC successfully completed - elapsed: %s' % (end - start))

    # Set CRS fields
    g1_polys_rho = POLYS_RHO(g1_polys, Chi.g1_rho)
    g2_polys_rho = POLYS_RHO(g2_polys, Chi.g2_rho)
    crs_sm = CRS_SM(g1_beta_polys, g1_beta_rho, g1_betahat_rhohat, Chi.g2_beta, Chi.g2_betahat)
    crs_pm = CRS_PM(g1_poly_zero, g1_ques_rho, g1_poly_sum, g1_polyhats_sum,
        g2_poly_zero, g2_poly_sum)
    crs_con = CRS_CON(g1_polyhats, Chi.g1_rhohat)
    crs_uv = CRS_UV(g1_betasquare_rho, g1_beta_betahat, g1_b_beta_polys,
            g2_betasquare, g2_beta_betahat)
    crs_pkv = CRS_PKV(Chi.g1_beta, Chi.g1_betahat, g1_thetapows, Chi.g1_chi,
        Chi.g2_chi, Chi.g1_theta, g2_theta, Chi.g2_beta, Chi.g2_betahat)

    crs = CRS(gk, g1_polys_rho, g2_polys_rho, crs_sm, crs_pm, crs_con, crs_uv, crs_pkv)
    return crs


# Layers

def layer_C1_1(gk, Chi_share, mpc):
    print(' * Start Multi-Party Computation')
    global start
    start = datetime.datetime.now()

    print("\n Layer C1.1")
    md_list = []

    md_list.append(MD_Elems(gk, Chi_share.chi, 1, gk.g1))
    md_list.append(MD_Elems(gk, Chi_share.chi, 1, gk.g2))
    md_list.append(MD_Elems(gk, Chi_share.theta, 1, gk.g1))

    md_list.append(MD_Elems(gk, Chi_share.beta, 1, gk.g1))
    md_list.append(MD_Elems(gk, Chi_share.beta, 1, gk.g2))

    md_list.append(MD_Elems(gk, Chi_share.betahat, 1, gk.g1))
    md_list.append(MD_Elems(gk, Chi_share.betahat, 1, gk.g2))

    md_list.append(MD_Elems(gk, Chi_share.rho, 1, gk.g1))
    md_list.append(MD_Elems(gk, Chi_share.rho, 1, gk.g2))

    md_list.append(MD_Elems(gk, Chi_share.rhohat, 1, gk.g1))

    (g1_chi, g2_chi, g1_theta, g1_beta, g2_beta, g1_betahat,
        g2_betahat, g1_rho, g2_rho, g1_rhohat) = mpc.parallel_md(md_list)

    return Chi(g1_chi, g2_chi, g1_theta, g1_beta, g2_beta,
        g1_betahat, g2_betahat, g1_rho, g2_rho, g1_rhohat)


def layer_C1_2(gk, Chi_share, Chi, n, mpc):
    print("\n Layer C1.2")
    g1_chipows, g1_thetapows = [gk.g1], [gk.g1]
    share_chi, share_theta = Chi_share.chi, Chi_share.theta
    g1_chi, g1_theta = Chi.g1_chi, Chi.g1_theta
    share_chi_pow, share_theta_pow = Bn(1), Bn(1)
    q = gk.q
    md_list_x = []
    md_list_u = []
    for i in range(1, 2 * n + 1):
        md_list_x.append(MD_Elems(gk, share_chi_pow, 1, g1_chi))
        share_chi_pow = share_chi * share_chi_pow % gk.q
        md_list_u.append(MD_Elems(gk, share_theta_pow, 1, g1_theta))
        share_theta_pow = share_theta * share_theta_pow % gk.q
    g1_chipows.extend(mpc.parallel_md(md_list_x))
    g1_thetapows.extend(mpc.parallel_md(md_list_u))

    return g1_chipows, g1_thetapows


def layer_C1_3(gk, Chi_share, g1_thetapows, n, mpc):
    print("\n Layer C1.3")
    (g2_theta,) = mpc.parallel_md([MD_Elems(gk, Chi_share.theta, 1, gk.g2)])
    g1_polyhats = [g1_thetapows[2 * i] for i in range(1, n + 1)]

    return g2_theta, g1_polyhats


def layer_C1_4(gk, Chi_share, Chi, g1_chipows, g1_polyhats, n, mpc):
    print("\n Layer C1.4")
    g2_chipows = [gk.g2]
    g1_beta_chipows = []
    g1_betahat_polyhats = []

    share_chi = Chi_share.chi
    share_beta = Chi_share.beta
    share_betahat = Chi_share.betahat

    g2_chi = Chi.g2_chi

    share_chi_pow = Bn(1)
    md_list_x  = []
    md_list_bx = [MD_Elems(gk, share_beta, 1, g1_chipows[0])]
    md_list_bp = []
    q = gk.q
    share_chi = Chi_share.chi
    for i in range(1, n + 1):
        md_list_x.append(MD_Elems(gk, share_chi_pow, 1, g2_chi))
        share_chi_pow = share_chi * share_chi_pow % q

        md_list_bx.append(MD_Elems(gk, share_beta, 1, g1_chipows[i]))
        md_list_bp.append(MD_Elems(gk, share_betahat, 1, g1_polyhats[i - 1]))

    g2_chipows.extend(mpc.parallel_md(md_list_x))
    g1_beta_chipows.extend(mpc.parallel_md(md_list_bx))
    g1_betahat_polyhats.extend(mpc.parallel_md(md_list_bp))

    return g2_chipows, g1_beta_chipows, g1_betahat_polyhats


def layer_C1_5(gk, Chi_share, Chi, mpc):
    print("\n Layer C1.5")
    md_list = []

    md_list.append(MD_Elems(gk, Chi_share.beta, 1, Chi.g2_beta))
    md_list.append(MD_Elems(gk, Chi_share.beta, 1, Chi.g1_rho))
    md_list.append(MD_Elems(gk, Chi_share.beta ** 2, 1, Chi.g1_rho))
    md_list.append(MD_Elems(gk, Chi_share.beta, 1, Chi.g1_betahat))
    md_list.append(MD_Elems(gk, Chi_share.beta, 1, Chi.g2_betahat))
    md_list.append(MD_Elems(gk, Chi_share.betahat, 1, Chi.g1_rhohat))

    (g2_betasquare, g1_beta_rho, g1_betasquare_rho,
        g1_beta_betahat, g2_beta_betahat, g1_betahat_rhohat) = mpc.parallel_md(md_list)

    return (g2_betasquare, g1_beta_rho, g1_betasquare_rho,
                g1_beta_betahat, g2_beta_betahat, g1_betahat_rhohat)


def layer_L1_1(gk, g1_chipows, g1_beta_chipows, g2_chipows, n):
    print("\n Layer L1.1")
    (g1_els, g1_elinplus1_prods, g1_el_squares, g2_els, g1_beta_els) = \
        generate_locally(gk, g1_chipows, g2_chipows, g1_beta_chipows, n)

    return g1_els, g2_els, g1_el_squares, g1_beta_els, g1_elinplus1_prods


def layer_L1_2(gk, g1_el_nplus1, g2_el_nplus1):
    print(" Layer L1.2")
    g1_poly_zero = lc((1, -1), (g1_el_nplus1, gk.g1))
    g2_poly_zero = lc((1, -1), (g2_el_nplus1, gk.g2))
    return g1_poly_zero, g2_poly_zero


def layer_L1_3(gk, n, Chi_share, g1_els, g1_el_nplus1, g2_els, g1_el_squares,
        g1_beta_el_nplus1, g1_el_nplus1_sq, g1_elinplus1_prods,
        g1_betahat_polyhats, g1_beta_els):
    print(" Layer L1.3")

    g1_polys = [lc((2, 1), (g1_el, g1_els[-1])) for g1_el in g1_els[:-1]]
    g2_polys = [lc((2, 1), (g2_el, g2_els[-1])) for g2_el in g2_els[:-1]]

    g1_ques = [lc((4, 4, 8, -4, -4), (
        g1_el_square,
        g1_el_nplus1_sq,
        g1_elinplus1_prod,
        g1_el,
        g1_el_nplus1)) for (g1_el_square, g1_elinplus1_prod, g1_el) in
            zip(g1_el_squares[:-1], g1_elinplus1_prods, g1_els[:-1])]

    g1_beta_polys = [lc((2, 1, 1),(
        g1_beta_el,
        g1_beta_el_nplus1,
        g1_betahat_polyhat))
    for (g1_beta_el, g1_betahat_polyhat) in
        zip(g1_beta_els[:-1], g1_betahat_polyhats)]

    return g1_polys, g2_polys, g1_ques, g1_beta_polys


def layer_L1_4(gk, g1_polys, g2_polys, g1_polyhats):
    print(" Layer L1.4")
    g1_poly_sum = sum(g1_polys, 0 * gk.g1)
    g2_poly_sum = sum(g2_polys, 0 * gk.g2)
    g1_polyhats_sum = sum(g1_polyhats, 0 * gk.g1)
    return g1_poly_sum, g2_poly_sum, g1_polyhats_sum


def layer_C2(gk, Chi_share, g1_ques, g1_beta_polys, mpc):
    print("\n Layer C2")
    md_list = []
    md_list.extend([MD_Elems(gk, 1, Chi_share.rho, _) for _ in g1_ques])
    md_list.extend([MD_Elems(gk, Chi_share.beta, 1, _) for _ in g1_beta_polys])
    result = mpc.parallel_md(md_list)

    g1_ques_rho = result[:len(g1_ques)]
    g1_b_beta_polys = result[len(g1_ques):]

    return g1_ques_rho, g1_b_beta_polys


# Helpers

def lc(cs, elems):
    """
    Implementation of lincomb, 3.1, pg. 11
    """
    assert elems
    g_zero = 0 * elems[0]
    return sum((c * elem for (c, elem) in zip(cs, elems)), g_zero)


def generate_locally(gk, g1_chipows, g2_chipows, g1_beta_chipows, n):
    """
    Computes (Θ(n ^ 2)) the following collections:

    [l_i]_1,           1<=i<=n+1
    [l_i * l_(n+1)]_1, 1<=i<=n+1
    [l_i ^ 2]_1,       1<=i<=n+1
    [l_i]_2,           1<=i<=n+1
    [β * l_i]_1,       1<=i<=n+1

    Computations follow directly from Theorem 7, C.2, pg. 26
    """
    q = gk.q
    omega = (q - 1) // (n + 1)
    div_nplus1 = lambda _: _.mul(Bn(n + 1).mod_inverse(q))
    div_nplus1_sq = lambda _: _.mul(Bn((n + 1) ** 2).mod_inverse(q))

    g1_zero = 0 * gk.g1
    g2_zero = 0 * gk.g2

    def compute_g1_elij_prod(i, j):
        g1_el_ij = g1_zero
        for s in range(0, n + 1):
            for t in range(0, n + 1):
                omega_pow_inv = pow(omega, i * s + j * t, q).mod_inverse(q)
                g1_el_ij += g1_chipows[s + t].mul(omega_pow_inv)
        return g1_el_ij

    g1_elinplus1_prods = []
    g1_el_squares = []
    g1_els = []
    g2_els = []
    g1_beta_els = []
    for i in range(1, n + 2):
        g1_el_i = g1_zero
        g2_el_i = g2_zero
        g1_beta_el_i = g1_zero

        for j in range(0, n + 1):
            omega_ij_inv = pow(omega, i * j, q).mod_inverse(q)

            g1_el_i += g1_chipows[j].mul(omega_ij_inv)
            g2_el_i += g2_chipows[j].mul(omega_ij_inv)
            g1_beta_el_i += g1_beta_chipows[j].mul(omega_ij_inv)

        g1_el_square = compute_g1_elij_prod(i, i)
        g1_elinplus1_prod = compute_g1_elij_prod(i, n + 1)

        g1_els.append(div_nplus1(g1_el_i))
        g1_el_squares.append(div_nplus1_sq(g1_el_square))
        g1_elinplus1_prods.append(div_nplus1_sq(g1_elinplus1_prod))
        g2_els.append(div_nplus1(g2_el_i))
        g1_beta_els.append(div_nplus1(g1_beta_el_i))

    return g1_els, g1_elinplus1_prods, g1_el_squares, g2_els, g1_beta_els
