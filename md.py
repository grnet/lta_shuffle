import ast
from utils import MD_Elems, to_Bn, mk_gk
from bbparty import BBParty, Quit
from petlib import bn, pack
from bplib import bp



def mpcMD(gk, alpha, beta, elem):
    if isinstance(beta, int):
        beta = bn.Bn(beta)
    inv_beta = beta.mod_inverse(gk.q)
    alpha_open, beta_open = create_openings(gk, alpha, beta, elem)
    return (alpha * inv_beta * elem, alpha_open, beta_open)


def vmpcMD(orig_elem, certs, alpha_openings, beta_openings, gk):
    prev = orig_elem
    for cert, alpha, beta in zip(certs, alpha_openings, beta_openings):
        pair1 = gk.pair(cert, beta) if isinstance(cert, bp.G1Elem)\
                                 else gk.pair(beta, cert)
        pair2 = gk.pair(prev, alpha) if isinstance(cert, bp.G1Elem)\
                                  else gk.pair(alpha, prev)
        if pair1 != pair2:
            return False
        prev = cert
    return True


def create_openings(gk, alpha, beta, elem):
    if isinstance(elem, bp.G1Elem):
        gx = gk.g2
    else:
        gx = gk.g1
    return (alpha * gx, beta * gx)


def simulate_mpc(mpc_list, transcript=None):
    if transcript is None:
        transcript = []

    stepno = len(transcript)
    mpc_iter = iter(mpc_list)
    for mpc in mpc_iter:
        owner = mpc.get_step_owner(stepno)
        break

    for mpc in mpc_iter:
        assert owner == mpc.get_step_owner(stepno)

    mpc_owner = mpc_list[owner]
    step_output = mpc_owner.produce_step_output()
    transcript.append(step_output)
    for i, mpc in enumerate(mpc_list):
        if i == owner:
            continue
        mpc.consume_step_output(step_output)

    return transcript


class MPCMD(BBParty):
    def parallel_md(self, md_list):
        self.md_list = md_list
        self.elem_state = [x.elem for x in self.md_list]
        self.transc_len_before = len(self.transcript)
        BBParty.run(self)
        return self.get_result()

    def get_step_owner(self, step_index):
        if self.md_is_finished():
            if self.verify_mpcmd():
                print(" [+] Successful MD")
                raise Quit
            else:
                print("Multi party computation can't be verified. Aborting.")
                raise SystemExit(1)
        return BBParty.get_step_owner(self, step_index)

    def md_is_finished(self):
        return len(self.transcript) - self.transc_len_before == len(self.peers)

    def produce_step_output(self):
        commitments = []
        for (md_index, md) in enumerate(self.md_list):
            commitment = mpcMD(md.gk, md.alpha, md.beta, self.elem_state[md_index])
            commitments.append(commitment)
            self.elem_state[md_index] = commitment[0]
        return pack.encode(commitments)

    def consume_step_output(self, step_output):
        commitments = pack.decode(step_output)
        for (md_index, c) in enumerate(commitments):
            self.elem_state[md_index] = c[0]
        return step_output

    def verify_mpcmd(self):
        cur_md_transcript = self.transcript[self.transc_len_before:]
        if len(cur_md_transcript) != len(self.peers):
            return False
        commitments = map(pack.decode, cur_md_transcript)
        for md in range(len(self.md_list)):
            certs = [c[md][0] for c in commitments]
            alpha_openings = [c[md][1] for c in commitments]
            beta_openings = [c[md][2] for c in commitments]
            if not vmpcMD(self.md_list[md].elem, certs, alpha_openings,
                    beta_openings, self.md_list[md].gk):
                return False
        return True

    def get_result(self):
        return self.elem_state


def print_help():
    import sys
    print("Usage: {0} <config_file> <peer_name>".format(sys.argv[0]))
    raise SystemExit(1)



if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print_help()

    config_file = sys.argv[1]
    with open(config_file) as f:
        config = ast.literal_eval(f.read())

    peers = config['peers']
    name = sys.argv[2]
    spec = {}
    starting_elem = to_Bn(11466060538852421779614017077833160237)

    gk = mk_gk()
    q = gk.q
    elem = starting_elem * gk.g1
    alpha = 1 + (q - 1).random()
    beta = 1 + (q - 1).random()
    md_1 = MD_Elems(gk, alpha, beta, elem)
    alpha = 1 + (q - 1).random()
    md_2 = MD_Elems(gk, alpha, beta, elem)

    mpc = MPCMD(peers, name, spec)
    result = mpc.parallel_md([md_1, md_2])
    print("Successful Multi Party Computation")
