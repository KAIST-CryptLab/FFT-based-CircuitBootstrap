import argparse
from prettytable import PrettyTable

load("var.sage")

def optimize_log_B_ks(N, k, n, log_q, Var_LWE, l_ks):
    q = 2^log_q
    cur = log_q // 2

    # Optimize B_ks
    while cur < log_q:
        Var_cur_gadget = get_var_lwe_ks_gadget(N, k, q, 2^cur, l_ks)
        Var_cur_key = get_var_lwe_ks_key(N, k, q, Var_LWE, 2^cur, l_ks)
        Var_cur = Var_cur_gadget + Var_cur_key

        next = cur + 1 if Var_cur_gadget > Var_cur_key else cur - 1
        Var_next = get_var_lwe_ks(N, k, q, Var_LWE, 2^next, l_ks)

        if Var_next > Var_cur:
            break

        cur = next

    return cur

def optimize_log_B_pbs(N, k, n, log_q, Var_GLWE, l_pbs):
    q = 2^log_q
    cur = log_q // 2

    # Optimize B_pbs by Var_pbs
    while cur < log_q:
        Var_cur_gadget = get_var_pbs_gadget(N, k, n, q, 2^cur, l_pbs)
        Var_cur_key = get_var_pbs_key(N, k, n, q, Var_GLWE, 2^cur, l_pbs)
        Var_cur = Var_cur_gadget + Var_cur_key

        next = cur + 1 if Var_cur_gadget > Var_cur_key else cur - 1
        Var_next = get_var_pbs(N, k, n, q, Var_GLWE, 2^next, l_pbs)

        if Var_next > Var_cur:
            break

        cur = next

    # Optimize B_pbs by Var_pbs + Var_fft_pbs
    Var_fft_cur = get_var_fft_pbs(N, k, n, 2^cur, l_pbs)
    if Var_cur < Var_fft_cur:
        Var_cur_tot = 0
        Var_next_tot = -1
        next = cur

        while Var_cur_tot > Var_next_tot:
            cur = next
            Var_cur_tot = get_var_pbs(N, k, n, q, Var_GLWE, 2^cur, l_pbs)
            Var_cur_tot += get_var_fft_pbs(N, k, n, 2^cur, l_pbs)

            next = cur - 1
            Var_next_tot = get_var_pbs(N, k, n, q, Var_GLWE, 2^next, l_pbs)
            Var_next_tot += get_var_fft_pbs(N, k, n, 2^next, l_pbs)

    return cur

def optimize_log_B_tr(N, k, log_q, Var_GLWE, l_tr, split_limit):
    q = 2^log_q
    cur = log_q // 2

    # Optimize B_tr
    while cur < log_q:
        Var_cur_gadget = get_var_glwe_ks_gadget(N, k, q, 2^cur, l_tr)
        Var_cur_key = get_var_glwe_ks_key(N, k, q, Var_GLWE, 2^cur, l_tr)
        Var_cur = get_var_tr(N, k, q, Var_GLWE, 2^cur, l_tr)

        next = cur + 1 if Var_cur_gadget > Var_cur_key else cur - 1
        Var_next = get_var_tr(N, k, q, Var_GLWE, 2^next, l_tr)

        if Var_next > Var_cur:
            break

        cur = next
    log_B_tr = cur
    B_tr = 2^log_B_tr

    # Optimize b_tr
    cur = log_q // 2
    _, fp_split_fft = get_fp_split_fft_glwe_ks(N, k, q, B_tr, l_tr, 2^cur)
    log_fp_split_fft = log(fp_split_fft, 2).n(10000)

    if log_fp_split_fft < split_limit:
        next = cur
        while log_fp_split_fft < split_limit and cur > 0:
            cur = next
            next = cur - 1
            _, fp_split_fft = get_fp_split_fft_glwe_ks(N, k, q, B_tr, l_tr, 2^next)
            log_fp_split_fft = log(fp_split_fft, 2).n(10000)
    else:
        while log_fp_split_fft > split_limit and cur < log_q:
            next = cur + 1
            cur = next
            _, fp_split_fft = get_fp_split_fft_glwe_ks(N, k, q, B_tr, l_tr, 2^cur)
            log_fp_split_fft = log(fp_split_fft, 2).n(10000)

    log_b_tr = cur
    b_tr = 2^log_b_tr
    _, fp_split_fft = get_fp_split_fft_glwe_ks(N, k, q, B_tr, l_tr, b_tr)
    log_fp_split_fft = log(fp_split_fft, 2)

    if log_fp_split_fft > split_limit:
        log_b_tr = None

    return log_B_tr, log_b_tr

def optimize_log_B_ss(N, k, log_q, Var_GLWE, l_ss):
    q = 2^log_q
    cur = log_q // 2

    # Optimize B_ss by Var_ss
    while cur < log_q:
        Var_cur_gadget = get_var_ss_gadget(N, k, q, 2^cur, l_ss)
        Var_cur_key = get_var_ss_inc(N, k, q^2 * Var_GLWE, 2^cur, l_ss)
        Var_cur = Var_cur_gadget + Var_cur_key

        next = cur + 1 if Var_cur_gadget > Var_cur_key else cur - 1
        Var_next = get_var_ss(N, k, q, q^2 * Var_GLWE, 2^next, l_ss)

        if Var_next > Var_cur:
            break

        cur = next

    # Optimize B_ss by Var_ss + Var_fft_ss
    Var_fft_cur = get_var_fft_ss(N, k, q, 2^cur, l_ss)
    if Var_cur < Var_fft_cur:
        Var_cur_tot = 0
        Var_next_tot = -1
        next = cur

        while Var_cur_tot > Var_next_tot:
            cur = next
            Var_cur_tot = get_var_ss(N, k, q, q^2 * Var_GLWE, 2^cur, l_ss)
            Var_cur_tot += get_var_fft_ss(N, k, q, 2^cur, l_ss)

            next = cur - 1
            Var_next_tot = get_var_ss(N, k, q, q^2 * Var_GLWE, 2^next, l_ss)
            Var_next_tot += get_var_fft_ss(N, k, q, 2^next, l_ss)

    return cur

def optimize_log_B_cbs(N, k, q, Var_cbs, l_cbs):
    q = 2^log_q
    cur = log_q // 2

    # Optimize B_cbs by Var_add
    while cur < log_q:
        Var_add_gadget = get_var_ext_prod_gadget(N, k, q, 2^cur ,l_cbs)
        Var_add_inc = get_var_ext_prod_inc(N, k, Var_cbs, 2^cur, l_cbs)
        Var_cur = Var_add_gadget + Var_add_inc

        next = cur + 1 if Var_add_gadget > Var_add_inc else cur - 1
        Var_next = get_var_ext_prod(N, k, q, Var_cbs, 2^next, l_cbs)

        if Var_next > Var_cur:
            break

        cur = next

    # Optimize B_cbs by Var_add + Var_fft_add
    Var_fft_cur = get_var_fft_ext_prod(N, k, q, 2^cur, l_cbs)
    if Var_cur < Var_fft_cur:
        Var_cur_tot = 0
        Var_next_tot = -1
        next = cur

        while Var_cur_tot > Var_next_tot:
            cur = next
            Var_cur_tot = get_var_ext_prod(N, k, q, Var_cbs, 2^cur, l_cbs)
            Var_cur_tot += get_var_fft_ext_prod(N, k, q, 2^cur, l_cbs)

            next = cur - 1
            Var_next_tot = get_var_ext_prod(N, k, q, Var_cbs, 2^next ,l_cbs)
            Var_next_tot += get_var_fft_ext_prod(N, k, q, 2^next ,l_cbs)
    return cur

parser = argparse.ArgumentParser(prog='sage finding_param.sage', description='finding appropriate decomposition base logs under the given TFHE parameters and decomposition levels')
parser.add_argument('-n', nargs=1, type=int, default=[571], help='n, default to 571')
parser.add_argument('-lwe', nargs=1, type=float, default=[1.953125e-4], help='LWE std dev, default to 1.953125e-4')
parser.add_argument('-k', nargs=1, type=int, default=[1], help='k, default to 1')
parser.add_argument('-N', nargs=1, type=int, default=[2048], help='N, default to 2048')
parser.add_argument('-glwe', nargs=1, type=float, default=[1.7763568e-16], help='GLWE std dev, default to 1.7763568e-16')
parser.add_argument('-log_q', nargs=1, type=int, default=[64], help='log q, default to 64')
parser.add_argument('-t', nargs=1, type=int, default=[2], help='vartheta, default to 2')
parser.add_argument('-l', nargs=5, type=int, help='decomposition levels: [pbs, ks, tr, ss, cbs]')
parser.add_argument('-B', nargs=6, type=int, help='log of decomposition bases: [pbs, ks, tr, split, ss, cbs]')
parser.add_argument('-sp', nargs=1, type=int, default=[-2000], help='log of f.p. of split fft, default to 2000')
parser.add_argument('-thrs', nargs='?', type=int, default=[-32, -80, -128], help='log of threshold f.p., default to [-32, -80, -128]')

if __name__ == '__main__':
    args = parser.parse_args()
    n = args.n[0]
    sigma_lwe = args.lwe[0]
    k = args.k[0]
    N = args.N[0]
    sigma_glwe = args.glwe[0]
    log_q = args.log_q[0]
    q = 2^log_q
    theta = args.t[0]
    split_limit = args.sp[0]
    is_opt = args.B is None
    log_fp_thrs_list = args.thrs

    print(f"LWE Dimension (-n):\t{n}")
    print(f"LWE Std. Dev (-lwe):\t{sigma_lwe:.5e}")
    print(f"GLWE Dimension (-k):\t{k}")
    print(f"Polynomial Size (-N):\t{N}")
    print(f"GLWE Std. Dev (-glwe):\t{sigma_glwe:.5e}")
    print(f"q (-log_q):\t\t2^{log_q}")
    print(f"theta (-t):\t\t{theta}")
    print()

    if args.l is None or len(args.l) != 5:
        print("Give 5 gadget lengths by the -l flag: -l l_ks l_pbs l_tr l_ss l_cbs")
        exit()
    else:
        [l_ks, l_pbs, l_tr, l_ss, l_cbs] = args.l

    Var_LWE = sigma_lwe^2
    Var_GLWE = sigma_glwe^2
    # log_fp_thrs_list = [-32, -80, -128]

    if is_opt:
        log_B_ks = optimize_log_B_ks(N, k, n, log_q, Var_LWE, l_ks)
        B_ks = 2^log_B_ks

        log_B_pbs = optimize_log_B_pbs(N, k, n, log_q, Var_GLWE, l_pbs)
        B_pbs = 2^log_B_pbs

        log_B_tr, log_b_tr = optimize_log_B_tr(N, k, log_q, Var_GLWE, l_tr, split_limit)
        B_tr = 2^log_B_tr
        b_tr = 2^log_b_tr

        log_B_ss = optimize_log_B_ss(N, k, log_q, Var_GLWE, l_ss)
        B_ss = 2^log_B_ss

        Var_pbs_tot = get_var_pbs(N, k, n, q, Var_GLWE, B_pbs, l_pbs)
        Var_pbs_tot += get_var_fft_pbs(N, k, n, B_pbs, l_pbs)
        Var_tr_tot = get_var_tr(N, k, q, Var_GLWE, B_tr, l_tr)
        Var_tr_tot += get_var_fft_tr(N, k, B_tr, l_tr, b_tr)
        Var_ss_tot = get_var_ss(N, k, q, q^2 * Var_GLWE, B_ss, l_ss)
        Var_ss_tot += get_var_fft_ss(N, k, q, B_ss, l_ss)

        Var_cbs = (Var_pbs_tot + Var_ss_tot) + (N/2) * Var_tr_tot
        log_B_cbs = optimize_log_B_cbs(N, k, q, Var_cbs, l_cbs)
        B_cbs = 2^log_B_cbs
    else:
        [log_B_ks, log_B_pbs, log_B_tr, log_b_tr, log_B_ss, log_B_cbs] = args.B
        [B_ks, B_pbs, B_tr, b_tr, B_ss, B_cbs] = [2^log_B_ks, 2^log_B_pbs, 2^log_B_tr, 2^log_b_tr, 2^log_B_ss, 2^log_B_cbs]

    tab = PrettyTable(['', 'l', 'B', 'b'])
    tab.add_row(['LWE KS', l_ks, f'2^{log_B_ks}', ''])
    tab.add_row(['PBS', l_pbs, f'2^{log_B_pbs}', ''])
    tab.add_row(['Trace', l_tr, f'2^{log_B_tr}', f'2^{log_b_tr}'])
    tab.add_row(['SS', l_ss, f'2^{log_B_ss}', ''])
    tab.add_row(['CBS', l_cbs, f'2^{log_B_cbs}', ''])
    print(tab)
    print()

    Var_lwe_ks = get_var_lwe_ks(N, k, q, Var_LWE, B_ks, l_ks)
    print(f"Var_lwe_ks:\t\t2^{log(Var_lwe_ks, 2).n():.4f}")

    Var_pbs = get_var_pbs(N, k, n, q, Var_GLWE, B_pbs, l_pbs)
    Var_fft_pbs = get_var_fft_pbs(N, k, n, B_pbs, l_pbs)
    Var_pbs_tot = Var_pbs + Var_fft_pbs

    print(f"Var_pbs_tot:\t\t2^{log(Var_pbs_tot, 2).n():.4f}")
    print(f"  - Var_pbs:\t\t2^{log(Var_pbs, 2).n():.4f}")
    print(f"  - Var_fft_pbs:\t2^{log(Var_fft_pbs, 2).n():.4f}")

    Var_tr = get_var_tr(N, k, q, Var_GLWE, B_tr, l_tr)
    Var_fft_tr = get_var_fft_tr(N, k, B_tr, l_tr, b_tr)
    Var_tr_tot = Var_tr + Var_fft_tr
    _, fp_split_fft = get_fp_split_fft_glwe_ks(N, k, q, B_tr, l_tr, b_tr)

    print(f"Var_tr_tot:\t\t2^{log(Var_tr_tot, 2).n():.4f}")
    print(f"  - Var_tr:\t\t2^{log(Var_tr, 2).n():.4f}")
    print(f"  - Var_fft_tr:\t\t2^{log(Var_fft_tr, 2).n():.4f}")
    print(f"  - F.P. of split FFT:\t2^{log(fp_split_fft, 2).n(10000):.4f}")

    Var_ss = get_var_ss(N, k, q, q^2 * Var_GLWE, B_ss, l_ss)
    Var_fft_ss = get_var_fft_ss(N, k, q, B_ss, l_ss)
    Var_ss_tot = Var_ss + Var_fft_ss

    print(f"Var_ss_tot:\t\t2^{log(Var_ss_tot, 2).n():.4f}")
    print(f"  - Var_ss:\t\t2^{log(Var_ss, 2).n():.4f}")
    print(f"  - Var_fft_ss:\t\t2^{log(Var_fft_ss, 2).n():.4f}")

    Var_cbs = (Var_pbs_tot + Var_ss_tot) + (N/2) * Var_tr_tot
    Var_add = get_var_ext_prod(N, k, q, Var_cbs, B_cbs, l_cbs)
    Var_fft_add = get_var_fft_ext_prod(N, k, q, B_cbs, l_cbs)
    Var_add_tot = Var_add + Var_fft_add

    print(f"Var_add_tot:\t\t2^{log(Var_add_tot, 2).n():.4f}")
    print(f"  - Var_cbs:\t\t2^{log(Var_cbs, 2).n():.4f}")
    print(f"  - Var_add:\t\t2^{log(Var_add, 2).n():.4f}")
    print(f"  - Var_fft_add:\t2^{log(Var_fft_add, 2).n():.4f}")

    print(f"max-depth for")
    for log_fp_thrs in log_fp_thrs_list:
        delta_in = q / 2
        log_var_thrs = find_var_thrs(n, q, N, theta, delta_in, log_fp_thrs)
        Var_thrs = 2^log_var_thrs
        max_depth = ((Var_thrs - Var_lwe_ks) / Var_add_tot).n()
        print(f"  - F.P. of 2^{log_fp_thrs}:\t{max_depth:.4f}")
    print()
