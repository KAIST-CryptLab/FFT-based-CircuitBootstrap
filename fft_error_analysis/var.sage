# -------- PBS -------- #
def get_var_pbs(N, k, n, q, Var_GLWE, B_pbs, l_pbs):
    Var_PBS = 0
    Var_PBS += get_var_pbs_gadget(N, k, n, q, B_pbs, l_pbs)
    Var_PBS += get_var_pbs_key(N, k, n, q, Var_GLWE, B_pbs, l_pbs)

    return Var_PBS

def get_var_pbs_gadget(N, k, n, q, B_pbs, l_pbs):
    Bp_2l_pbs = B_pbs^(2*l_pbs)

    Var_gadget = 0
    Var_gadget += n*(q^2-Bp_2l_pbs)/(24*Bp_2l_pbs) * (1+k*N/2)
    Var_gadget += (n*k*N)/32
    Var_gadget += (n/16)*(1-k*N/2)^2

    return Var_gadget

def get_var_pbs_key(N, k, n, q, Var_GLWE, B_pbs, l_pbs):
    B2_12_pbs = (B_pbs^2 + 2)/12
    Var_BSK = Var_GLWE * q^2

    return n*l_pbs*(k+1)*N*B2_12_pbs*Var_BSK

def get_var_fft_pbs(N, k, n, B_pbs, l_pbs):
    return n * 2^(22 - 2.6) * l_pbs * B_pbs^2 * N^2 * (k+1)

def get_fp_pbs(n, q_prime, N, theta, delta_in, Var_in):
    w = 2*N/2^theta
    Gamma = (delta_in / 2) * (Var_in + q_prime^2/(12*w^2) - 1/12 + n*q_prime^2/(24*w^2) + n/48)^(-1/2)
    sq2 = 2^(1/2)
    fp = 1 - erf(Gamma/sq2)

    return Gamma, fp

def get_min_fp_pbs(n, q_prime, N, theta, delta_in):
    w = 2*N/2^theta
    Gamma = (delta_in / 2) * (q_prime^2/(12*w^2) - 1/12 + n*q_prime^2/(24*w^2) + n/48)^(-1/2)
    sq2 = 2^(1/2)
    fp = 1 - erf(Gamma/sq2)

    return Gamma, fp

def find_var_thrs(n, q_prime, N, theta, delta_in, log_fp_thrs, prec=3):
    def get_log_fp(var_in):
        _, fp = get_fp_pbs(n, q_prime, N, theta, delta_in, var_in)
        return log(fp, 2).n(1000)

    _, min_fp = get_min_fp_pbs(n, q_prime, N, theta, delta_in)
    min_log_fp = log(min_fp, 2).n(1000)

    if log_fp_thrs < min_log_fp:
        return None

    int_lower = 0
    int_upper = 128
    pivot = 64
    log_fp = get_log_fp(2^pivot)
    while int_lower != int_upper - 1:
        if log_fp < log_fp_thrs:
            int_lower = pivot
        else:
            int_upper = pivot
        pivot = (int_lower + int_upper) // 2
        log_fp = get_log_fp(2^pivot)

    int_part = pivot
    frac_lower = 0
    frac_upper = (10^prec) - 1
    pivot = (frac_lower + frac_upper) // 2
    log_fp = get_log_fp(2^(int_part + pivot / 10^prec))
    while frac_lower != frac_upper - 1:
        if log_fp < log_fp_thrs:
            frac_lower = pivot
        else:
            frac_upper = pivot
        pivot = (frac_lower + frac_upper) // 2
        log_fp = get_log_fp(2^(int_part + pivot / 10^prec))

    return int_part + pivot / 10^prec

# -------- External Product -------- #
def get_var_ext_prod(N, k, q, Var_in, B_ep, l_ep):
    Var_ep = 0
    Var_ep += get_var_ext_prod_gadget(N, k, q, B_ep, l_ep)
    Var_ep += get_var_ext_prod_inc(N, k, Var_in, B_ep, l_ep)

    return Var_ep

def get_var_ext_prod_gadget(N, k, q, B_ep, l_ep):
    Bp_2l_ep = B_ep^(2*l_ep)
    return (1 + k*N)*((q^2 - Bp_2l_ep)/(24 * Bp_2l_ep) + 1/16)

def get_var_ext_prod_inc(N, k, Var_in, B_ep, l_ep):
    B2_12_ep = (B_ep^2 + 2)/12
    return (k+1)*l_ep*N * B2_12_ep * Var_in

def get_var_fft_ext_prod(N, k, q, B_ep, l_ep):
    return 2^(-2*53-2.6) * (k+1) * l_ep * B_ep^2 * q^2 * N^2

# -------- LWE KS -------- #
def get_var_lwe_ks(N, k, q, Var_LWE, B_ksk, l_ksk):
    Var_KS = 0
    Var_KS += get_var_lwe_ks_gadget(N, k, q, B_ksk, l_ksk)
    Var_KS += get_var_lwe_ks_key(N, k, q, Var_LWE, B_ksk, l_ksk)

    return Var_KS

def get_var_lwe_ks_gadget(N, k, q, B_ksk, l_ksk):
    Bp_2l_ksk = B_ksk^(2*l_ksk)

    return k*N*((q^2-Bp_2l_ksk)/(24*Bp_2l_ksk) + 1/16)

def get_var_lwe_ks_key(N, k, q, Var_LWE, B_ksk, l_ksk):
    B2_12_ksk = B_ksk^2 / 12
    Var_KSK = q^2 * Var_LWE

    return k*N*l_ksk*Var_KSK * (B2_12_ksk + 1/6)

# -------- GLWE KS -------- #
def get_var_glwe_ks(N, k_src, q, Var_dst, B_ksk, l_ksk):
    Var_KS = 0
    Var_KS += get_var_glwe_ks_gadget(N, k_src, q, B_ksk, l_ksk)
    Var_KS += get_var_glwe_ks_key(N, k_src, q, Var_dst, B_ksk, l_ksk)

    return Var_KS

def get_var_glwe_ks_gadget(N, k_src, q, B_ksk, l_ksk):
    return k_src * N * ((q^2 - B_ksk^(2*l_ksk)) / (24 * B_ksk^(2*l_ksk)) + 1/16)

def get_var_glwe_ks_key(N, k_src, q, Var_dst, B_ksk, l_ksk):
    Var_KSK = Var_dst * q^2
    return k_src * N * l_ksk * Var_KSK * (B_ksk^2 + 2) / 12

def get_var_fft_glwe_ks(N, k, B_ksk, l_ksk, b_fft):
    return 2^(-2*53-2.6) * k * l_ksk * B_ksk^2 * b_fft^2 * N^2

# -------- HomTrace -------- #
def get_var_tr(N, k, q, Var_ak, B_auto, l_auto):
    Var_auto = get_var_glwe_ks(N, k, q, Var_ak, B_auto, l_auto)
    Var_tr = (N^2 - 1)/3 *  Var_auto

    return Var_tr

def get_var_fft_tr(N, k, B_ksk, l_ksk, b_fft):
    return (N^2 - 1)/3 * get_var_fft_glwe_ks(N, k, B_ksk, l_ksk, b_fft)

def get_fp_split_fft_glwe_ks(N, k, q, B_ks, l_ks, b_fft):
    Var_fft_upper = get_var_fft_glwe_ks(N, k, B_ks, l_ks, q / b_fft)
    Gamma = 1/(2 * Var_fft_upper^(1/2))
    sq2 = 2^(1/2)
    fp = 1 - erf(Gamma/sq2)

    return Gamma, fp

# -------- SchemeSwitch -------- #
def get_var_ss(N, k, q, Var_in, B_ss, l_ss):
    Var_ss = 0
    Var_ss += get_var_ss_gadget(N, k, q, B_ss, l_ss)
    Var_ss += get_var_ss_inc(N, k, Var_in, B_ss, l_ss)

    return Var_ss

def get_var_ss_gadget(N, k, q, B_ss, l_ss):
    return get_var_ext_prod_gadget(N, k, q, B_ss, l_ss) * N/2

def get_var_ss_inc(N, k, Var_in, B_ss, l_ss):
    return get_var_ext_prod_inc(N, k, Var_in, B_ss, l_ss)

def get_var_fft_ss(N, k, q, B_ss, l_ss):
    return get_var_fft_ext_prod(N, k, q, B_ss, l_ss)
