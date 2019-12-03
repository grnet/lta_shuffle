from .elgamal import enc, mk_table, dec

def encrypt(q, pk, m):
    s = q.random()
    return enc(pk, s, m)

def make_table(pk, n):
    return mk_table(pk[0], n)

def decrypt(ctext, secret, table):
    return dec(secret, ctext, table)

__all__ = ('encrypt', 'make_table', 'decrypt')
