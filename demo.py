import random
import datetime
import argparse
import ast
import operator
import sys
from functools import reduce

import crs
import fixedcrs
import utils
import encdec
import prover
import verifier


system_random = random.SystemRandom()

def secure_shuffle(lst):
    random.shuffle(lst, random=system_random.random)

def random_permutation(n):
    s = list(range(n))
    secure_shuffle(s)
    return s

def mk_t_randoms(n, q):
    return [q.random() for _ in range(n)]

def make_key(gk):
    secret = gk.q.random()
    public = gk.g2, secret * gk.g2
    return utils.KeyPair(public, secret)

def encrypt_messages(order, pk, messages):
    return [encdec.encrypt(order, pk, message) for message in messages]

def decrypt_messages(secret, table, ciphertexts):
    return [encdec.decrypt(cs, secret, table) for cs in ciphertexts]

def check_namedtuples(t1, t2):
    d1 = t1._asdict()
    d2 = t2._asdict()
    previous_is_namedtuple = False
    for key, value in d1.items():
        other = d2[key]
        if hasattr(value, '_asdict'):
            print('\n Tuple: %s %s\n' % (key, (28 - len(key)) * '-'))
            check_namedtuples(value, other)
        else:
            sign, message = '+', 'ok'
            if value != other:
                global failed
                failed = True
                sign, message = '-', 'FAILED'
            print(' [%s] Checking %s: %s' % (sign, key, message))
            previous_is_namedtuple = False

def equal_lists(lst_1, lst_2):
    if len(lst_1) != len(lst_2):
        return False
    for elem_1, elem_2 in zip(lst_1, lst_2):
        if elem_1 != elem_2:
            return False
    return True

def spawn_mpc(gk, n, local_share, config):
    print('\n * Waiting for all peers to connect')
    CRS = crs.mk_crs(gk, n, local_share, config)
    return CRS

def run_fixed(gk, n, global_share):
    fixedCRS = fixedcrs.mk_crs(gk, n, Chi=global_share)
    return fixedCRS

def initialize(n, config):
    gk = utils.mk_gk()
    local_share, global_share = get_shares(config, gk.q)
    return gk, local_share, global_share

def run_computations(gk, n, local_share, global_share, config, fixed):
    CRS = spawn_mpc(gk, n, local_share, config)
    print(' * Verify computations:')
    fixedCRS = run_fixed(gk, n, global_share)
    global failed
    failed = False
    check_namedtuples(CRS, fixedCRS)
    print('\n * %s' % ('All checks PASSED' if not failed else 'Some checks FAILED'))
    return CRS if not fixed else fixedCRS

def get_shares(config, q):
    local_name = config['name']
    peers = config['peers']

    shares = []
    local_share = None
    for peer in peers:
        if 'chi_share' in peer:
            cur_chi_share = list(map(utils.to_Bn, peer['chi_share']))
            shares.append(cur_chi_share)
            if local_name == peer['name']:
                local_share = crs.Chi_share(*cur_chi_share)
    if len(peers) != len(shares) or local_share is None:
        share = crs.mk_Chi_share(q)
        return (share, share)

    print('Using shares found in config file: %s' % CONFIG)
    return (local_share, combine(shares))

def combine(shares):
    element_collections = zip(*shares)
    combined_share = []
    for elements in element_collections:
        combined_share.append(reduce(operator.mul, elements))
    return crs.Chi_share(*combined_share)

def demo(n, messages, config, batch=True, fixed=False):
    print('Maximum number of messages: n + 1 = %d' % (n + 1))

    gk, local_share, global_share = initialize(n, config)
    CRS = run_computations(gk, n, local_share, global_share, config, fixed)

    sigma = random_permutation(n)
    print()
    print(' * Proof and verification:')
    print('\n SIGMA = %s\n' % sigma)
    key = make_key(gk)
    pk_e = key.public
    ciphers = encrypt_messages(gk.q, key.public, messages)
    t_randoms = mk_t_randoms(n, gk.q)
    start = datetime.datetime.now()
    proof = prover.prove(n, pk_e, CRS, ciphers, sigma, t_randoms)
    perm_ok, cons_ok = verifier.verify(n, CRS, pk_e, ciphers, proof)
    end = datetime.datetime.now()
    failed = False
    for (verified, label) in ((perm_ok, 'permutation'), (cons_ok, 'consistency')):
        sign, message = '+', 'ok'
        if not verified:
            failed = True
            sign, message = '-', 'FAILED'
        print(' [%s] Check %s: %s' % (sign, label, message))
    print()
    print(' * Proof %sVERIFIED - elapsed: %s' % ('NOT ' if failed else '', end - start))

    print(' * Decryption:')
    start = datetime.datetime.now()
    shuffled_ciphertexts = proof.cprimes
    TABLE = encdec.make_table(pk_e, n)
    shuffled_ms = decrypt_messages(key.secret, TABLE, shuffled_ciphertexts)
    end = datetime.datetime.now()
    print('\n ', shuffled_ms, '\n')
    print(' * Decryption completed - elapsed: % s' % (end - start))


DEFAULT_N = 2 ** 3 - 1
CONFIG = './bbparty.config.local'
parser = argparse.ArgumentParser(description='Sub shuffler')
parser.add_argument('n', metavar='N', type=int, nargs='?', default=DEFAULT_N,
                    help='number of messages, default %s' % DEFAULT_N)
parser.add_argument('--fixed', action='store_true',
                    help='use fixed CRS version')
parser.add_argument('--nb', action='store_true',
                    default=False,
                    help="Don't batch verification (default is batching)")
parser.add_argument('--name', metavar='name', type=str,
                    help='Peer name', required=True)
parser.add_argument('--config', metavar='config', type=str, default=CONFIG,
                    help='Peers config, default %s' % CONFIG)


if __name__ == '__main__':
    args = parser.parse_args()
    CONFIG = args.config
    with open(CONFIG) as f:
        config = ast.literal_eval(f.read())
    config['name'] = args.name
    messages = list(range(args.n))
    demo(len(messages), messages, config, batch=not args.nb, fixed=args.fixed)
    sys.exit(0)
