import mdtraj as mdt
import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from argparse import ArgumentParser

# trajfile = sys.argv[1]
# topfile = sys.argv[2]
# outfile = sys.argv[3]
# water_resname = sys.argv[4]
# n_cycles = int(sys.argv[5])


def run():
    parser = ArgumentParser()
    parser.add_argument('--trajin', help='input trajectory')
    parser.add_argument('--topology', help='topology file')
    parser.add_argument('--trajout', help='output trajectory')
    parser.add_argument('--water_name', help='water residue name')
    parser.add_argument('--n_cycles', type=int, help='number of iterations')

    args = parser.parse_args()

    t = mdt.load(args.trajin, top=args.topology)
    w = t.topology.select(f'resname {args.water_name} and mass > 2')
    n_w = len(w)
    if n_w == 0:
        print(f'Error: no water molecules with residue name {args.water_name}')
        exit(1)
    dw = w[1:] - w[:-1]
    if not np.all(dw == dw[0]):
        print('Error - waters must be contiguous in the file')
        exit(1)
    n_aw = dw[0]
    if n_aw > 1:
        ro = t.topology.atom(w[0]).residue.index
        for j in range(1, n_aw):
            rh = t.topology.atom(w[0]+j).residue.index
            if rh != ro:
                print('Error - waters must have their O atoms listed first')
                exit(1)
    w_start = w[0]
    w_end = w[-1] + n_aw
    tw = t.atom_slice(w)
    x_opt = tw.xyz.copy()
    x_avg = x_opt[0]
    w_ind = np.zeros(n_w * n_aw, dtype=int)
    for i in range(args.n_cycles):
        cost = 0.0
        for j, x in enumerate(x_opt):
            dm = cdist(x_avg, x_opt[j])
            cost += np.diag(dm).mean()
            row_ind, col_ind = linear_sum_assignment(dm)
            x_opt[j] = x_opt[j][col_ind]
            w_ind[::n_aw] = col_ind * n_aw + w_start
            for k in range(1, n_aw):
                w_ind[k::n_aw] = w_ind[::n_aw] + k
            t.xyz[j, w_start:w_end] = t.xyz[j, w_ind]
        print(f'cycle {i}: mean cost: {cost/len(x_opt):8.2f}')
        x_avg = x_opt.mean(axis=0)
    t.save(args.trajout)


if __name__ == '__main__':
    run()
