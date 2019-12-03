from hashlib import sha256
from bplib import bp
from petlib import bn, pack
from collections import namedtuple


KeyPair = namedtuple('KeyPair', ['public', 'secret'])
Commitment = namedtuple('Commitment', ['cert', 'alpha_open', 'beta_open'])
MD_Elems = namedtuple('MD_Elems', ['gk', 'alpha', 'beta', 'elem'])


def to_Bn(num):
    return bn.Bn.from_decimal(str(num))


gk_T = namedtuple('gk_T', ['G', 'q', 'g1', 'g2', 'gt', 'pair'])

def mk_gk():
    G = bp.BpGroup()
    q = G.order()
    g1 = G.gen1()
    g2 = G.gen2()
    gt = G.pair(g1, g2)
    return gk_T(G, q, g1, g2, gt, G.pair)
