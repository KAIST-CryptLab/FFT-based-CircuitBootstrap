load("var.sage")

q = 2^64
n = 768
k = 2
logN = 10
N = 2**logN
theta = 3

N_common = 256
k_src = (k * N) // N_common

Var_GLWE_k1 = (0.00000000000000029403601535432533)^2
Var_GLWE_k2 = (0.0000000000000000002168404344971009)^2

Var_GLWE = Var_GLWE_k1
Var_LWE = (0.00000702047462940120)^2
B_hc_tr = 2^6
l_hc_tr = 8
B_hc_ss = 2^19
l_hc_ss = 2
B_hc_ggsw = 2^5
l_hc_ggsw = 5

B_tr = 2^13
l_tr = 3
B_ss = 2^17
l_ss = 2
B_ggsw = 2^2
l_ggsw = 7
B_ds = 2^4
l_ds = 3
B_pbs = 2^23
l_pbs = 1

amp_by_int_msg = 5
split_size = 1

Var_ssk = q^2 * Var_GLWE

print(f"n: {n}, N: {N}, k: {k}, Var_LWE: 2^{log(Var_LWE, 2).n():.4f}, Var_GLWE: 2^{log(Var_GLWE, 2).n():.4f}")
print(f"---- HalfCBS Parameters ----")
print(f"B_tr: 2^{log(B_hc_tr, 2)}, l_tr: {l_hc_tr}, B_ss: 2^{log(B_hc_ss, 2)}, l_ss: {l_hc_ss}, B_ggsw: 2^{log(B_hc_ggsw, 2)}, l_ggsw: {l_hc_ggsw}")

print(f"---- CBS Parameters ----")
print(f"B_tr: 2^{log(B_tr, 2)}, l_tr: {l_tr}, B_ss: 2^{log(B_ss, 2)}, l_ss: {l_ss}, B_ggsw: 2^{log(B_ggsw, 2)}, l_ggsw: {l_ggsw}, B_ds: 2^{log(B_ds, 2)}, l_ds: {l_ds}")
print()

print("==== 1st Round ====")
Var_keyed_lut = q^2 * Var_GLWE
Var_linear = 4 * Var_keyed_lut
print(f"Var_keyed_lut: 2^{log(Var_keyed_lut, 2).n():.4f}")
print(f"Var_linear: 2^{log(Var_linear, 2).n():.4f}")

print(f"\n==== 2nd Round (HalfCBS) ===")
Var_tr = get_var_tr(N, k, q, Var_GLWE, B_hc_tr, l_hc_tr)
Var_auto_gadget = get_var_glwe_ks_gadget(N, k, q, B_hc_tr, l_hc_tr)
Var_auto_key = get_var_glwe_ks_key(N, k, q, Var_GLWE, B_hc_tr, l_hc_tr)

print(f"Var_tr: 2^{log(Var_tr, 2).n():.4f}")
print(f"  - Var_auto_gadget: 2^{log(Var_auto_gadget, 2).n():.4f}")
print(f"  - Var_auto_key   : 2^{log(Var_auto_key, 2).n():.4f}")

Var_ss_gadget = get_var_ss_gadget(N, k, q, B_hc_ss, l_hc_ss)
Var_ss_key = get_var_ss_inc(N, k, Var_ssk, B_hc_ss, l_hc_ss)
Var_ss = get_var_ss(N, k, q, Var_ssk, B_hc_ss, l_hc_ss)

print(f"Var_ss: 2^{log(Var_ss, 2).n():.4f}")
print(f"  - Var_ss_gadget: 2^{log(Var_ss_gadget, 2).n():.4f}")
print(f"  - Var_ss_key   : 2^{log(Var_ss_key, 2).n():.4f}")

Var_cbs = Var_linear + (N/2) * Var_tr + Var_ss
Var_cbs_amp = (N/2) * Var_tr

print(f"Var_cbs: 2^{log(Var_cbs, 2).n():.4f}")
print(f"  - (N/2) * Var_tr: 2^{log(Var_cbs_amp, 2).n():.4f}")
print(f"  -         Var_ss: 2^{log(Var_ss, 2).n():.4f}")

Var_lut_ep_gadget = amp_by_int_msg * get_var_ext_prod_gadget(N, k, q, B_hc_ggsw, l_hc_ggsw)
Var_lut_ep_inc = get_var_ext_prod_inc(N, k, Var_cbs, B_hc_ggsw, l_hc_ggsw)

print("8-bit keyed-LUT")
print(f"  - Integer amplification in average: {amp_by_int_msg:.1f}")
print(f"  - Var_lut_ep_gadget: 2^{log(Var_lut_ep_gadget, 2).n():.4f} (integer message)")
print(f"  - Var_lut_ep_inc   : 2^{log(Var_lut_ep_inc, 2).n():.4f}")

Var_lut_out = 0
for i in range(8):
    Var_lut_out *= amp_by_int_msg
    print(f"  [{i+1}] Amp: 2^{log(Var_lut_out, 2).n():.4f}")
    Var_lut_out += Var_lut_ep_gadget + Var_lut_ep_inc
    print(f"      Inc: 2^{log(Var_lut_out, 2).n():.4f}")

Var_linear = 4 * Var_lut_out + q^2 * Var_GLWE
print(f"Var_linear: 2^{log(Var_linear, 2).n():.4f}")

Var_lwe_ks = get_var_glwe_ks(N_common, k_src, q, Var_LWE, B_ds, l_ds)
Var_lwe_ks_gadget = get_var_glwe_ks_gadget(N_common, k_src, q, B_ds, l_ds)
Var_lwe_ks_key = get_var_glwe_ks_key(N_common, k_src, q, Var_LWE, B_ds, l_ds)
print(f"Var_lwe_ks: 2^{log(Var_lwe_ks, 2).n():.4f}")
print(f"  - Var_lwe_ks_gadget: 2^{log(Var_lwe_ks_gadget, 2).n():.4f}")
print(f"  - Var_lwe_ks_key   : 2^{log(Var_lwe_ks_key, 2).n():.4f}")


Var_pbs_in = Var_linear + Var_lwe_ks
print()
print(f"Var_pbs_in: 2^{log(Var_pbs_in, 2).n():.4f}")

w = 2*N/(2^theta)
q_prime = q
delta_in = 2^63

Gamma, fp = get_fp_pbs(n, q_prime, N, theta, delta_in, Var_linear)
log_fp = log(fp, 2).n(1000)
print(f"F.P. of next PBS with theta = {theta}: Gamma = {Gamma}, f.p. = 2^{log_fp:.4f}")

print(f"\n==== 3rd Round (CBS) ====")
Var_pbs = get_var_pbs(N, k, n, q, Var_GLWE, B_pbs, l_pbs)
Var_pbs_fft = get_var_fft_pbs(N, k, n, B_pbs, l_pbs)
print(f"Var_pbs    : 2^{log(Var_pbs, 2).n():.4f}")
print(f"Var_pbs_fft: 2^{log(Var_pbs_fft, 2).n():.4f}")

Var_tr = get_var_tr(N, k, q, Var_GLWE, B_tr, l_tr)
Var_auto_gadget = get_var_glwe_ks_gadget(N, k, q, B_tr, l_tr)
Var_auto_key = get_var_glwe_ks_key(N, k, q, Var_GLWE, B_tr, l_tr)

print(f"Var_tr: 2^{log(Var_tr, 2).n():.4f}")
print(f"  - Var_auto_gadget: 2^{log(Var_auto_gadget, 2).n():.4f}")
print(f"  - Var_auto_key   : 2^{log(Var_auto_key, 2).n():.4f}")

Var_ss_gadget = get_var_ss_gadget(N, k, q, B_ss, l_ss)
Var_ss_key = get_var_ss_inc(N, k, Var_ssk, B_ss, l_ss)
Var_ss = get_var_ss(N, k, q, Var_ssk, B_ss, l_ss)

print(f"Var_ss: 2^{log(Var_ss, 2).n():.4f}")
print(f"  - Var_ss_gadget: 2^{log(Var_ss_gadget, 2).n():.4f}")
print(f"  - Var_ss_key   : 2^{log(Var_ss_key, 2).n():.4f}")

Var_cbs = Var_pbs + Var_pbs_fft + (N/2) * Var_tr + Var_ss
Var_cbs_amp = (N/2) * Var_tr

print(f"Var_cbs: 2^{log(Var_cbs, 2).n():.4f}")
print(f"  - Var_pbs (+ fft): 2^{log(Var_pbs + Var_pbs_fft, 2).n():.4f}")
print(f"  -  (N/2) * Var_tr: 2^{log(Var_cbs_amp, 2).n():.4f}")
print(f"  -          Var_ss: 2^{log(Var_ss, 2).n():.4f}")

Var_lut_ep = get_var_ext_prod(N, k, q, Var_cbs, B_ggsw, l_ggsw)
Var_lut_ep_gadget = get_var_ext_prod_gadget(N, k, q, B_ggsw, l_ggsw)
Var_lut_ep_inc = get_var_ext_prod_inc(N, k, Var_cbs, B_ggsw, l_ggsw)
print(f"Var_ep: 2^{log(Var_lut_ep, 2).n():.4f}")
print(f"  - Var_ep_gadget: 2^{log(Var_lut_ep_gadget, 2).n():.4f}")
print(f"  - Var_ep_inc   : 2^{log(Var_lut_ep_inc, 2).n():.4f}")

Var_8_lut = 8 * Var_lut_ep
Var_linear = 4 * Var_8_lut
print(f"Var_8_lut : 2^{log(Var_8_lut, 2).n():.4f}")
print(f"Var_linear: 2^{log(Var_linear, 2).n():.4f}")
print(f"Var_lwe_ks: 2^{log(Var_lwe_ks, 2).n():.4f}")

Var_pbs_in = Var_linear + Var_lwe_ks + q^2 * Var_GLWE
print()
print(f"Var_pbs_in: 2^{log(Var_pbs_in, 2).n():.4f}")

Gamma, fp = get_fp_pbs(n, q_prime, N, theta, delta_in, Var_linear)
log_fp = log(fp, 2).n(1000)
print(f"F.P. of next PBS with theta = {theta}: Gamma = {Gamma}, f.p. = 2^{log_fp:.4f}")
