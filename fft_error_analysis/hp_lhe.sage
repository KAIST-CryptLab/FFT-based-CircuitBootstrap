load("var.sage")

q = 2^64
N = 2048
k = 1
Var_GLWE = (0.00000000000000029403601535432533)^2
k_large = 2
Var_large_GLWE = (0.0000000000000000002168404344971009)^2



wopbs_2_2 = (
    "wopbs_param_message_2_carry_2_ks_pbs",
    (769, 0.0000043131554647504185^2), # (LWE dim, LWE var)
    (2^15, 2), # (PBS base, PBS level)
    (2^5, 3), # (CBS base, CBS level) (note: (2^6, 3) is better)
    (2^6, 2), # (KS base, KS level)
    # (0, 2^116.567, 2^117.270, 2^118.670), # (theta, Var_thrs_128, Var_thrs_80, Var_thrs_32)
    0, # theta
    4, # log_modulus
    1, # num_extract
)

wopbs_3_3 = (
    "wopbs_param_message_3_carry_3_ks_pbs",
    (873, 0.0000006428797112843789^2), # (LWE dim, LWE var)
    (2^9, 4), # (PBS base, PBS level)
    (2^6, 3), # (CBS base, CBS level)
    (2^10, 1), # (KS base, KS level)
    # (0, 2^116.566, 2^117.269, 2^118.669), # (theta, Var_thrs_128, Var_thrs_80, Var_thrs_32)
    0, # theta
    6, # log_modulus
    1, # num_extract
)

wopbs_4_4 = (
    "wopbs_param_message_4_carry_4_ks_pbs",
    (953, 0.0000001486733969411098^2), # (LWE dim, LWE var)
    (2^9, 4), # (PBS base, PBS level)
    (2^4, 6), # (CBS base, CBS level)
    (2^11, 1), # (KS base, KS level)
    # (0, 2^116.565, 2^117.269, 2^118.669), # (theta, Var_thrs_128, Var_thrs_80, Var_thrs_32)
    0, # theta
    8, # log_modulus
    1, # num_extract
)


param_list = [
    wopbs_2_2,
    wopbs_3_3,
    wopbs_4_4,
]

for param in param_list:
    name = param[0]
    (n, Var_LWE) = param[1]
    (B_pbs, l_pbs) = param[2]
    (B_cbs, l_cbs) = param[3]
    (B_ksk, l_ksk) = param[4]
    theta = param[5]
    log_modulus = param[6]
    max_num_extract = param[7]

    print(f"========================= {name} =========================")
    print(f"n: {n}, N: {N}, k: {k}, B_pbs: 2^{log(B_pbs, 2)}, l_pbs: {l_pbs}\n")

    Var_PBS = get_var_pbs(N, k, n, q, Var_GLWE, B_pbs, l_pbs)
    Var_FFT = get_var_fft_pbs(N, k, n, B_pbs, l_pbs)
    Var_KS = get_var_lwe_ks(N, k, q, Var_LWE, B_ksk, l_ksk)

    print(f"Var_PBS: 2^{log(Var_PBS, 2).n():.4f}")
    print(f"Var_FFT: 2^{log(Var_FFT, 2).n():.4f}")
    print(f"Var_KS: 2^{log(Var_KS, 2).n():.4f}")
    print()
    print(f"(B_cbs, l_cbs): (2^{log(B_cbs, 2)}, {l_cbs})")
    print(f"(B_ksk, l_ksk): (2^{log(B_ksk, 2)}, {l_ksk})")
    print()

    for num_extract in range(1, max_num_extract+1):
        Var_Add = get_var_ext_prod(N, k, q, Var_PBS + Var_FFT, B_cbs, l_cbs)
        # MV-PBS [CIM19] to extract num_extract bits
        if num_extract == 3:
            Var_Add *= 3
        Var_scaled_in = 2^(2*(log_modulus - num_extract)) * Var_Add

        print("-------------------------------------------------------------------------")
        print(f"# extracting bits: {num_extract} (MV-PBS [CIM19] is used after PBSmanyLUT [CLOT21] for LWEtoLev Conversion)")
        print(f"  - Var_Add: 2^{log(Var_Add, 2).n():.4f}")
        print(f"  - Var_scaled_in: 2^{log(Var_scaled_in, 2).n():.4f}")
        print()

        q_prime = q
        delta_in = 2^(64 - num_extract)
        print(f"Delta_in: 2^{log(delta_in, 2)}")
        print("theta:", theta)
        _, min_fp = get_min_fp_pbs(n, q_prime, N, theta, delta_in)
        log_min_fp = log(min_fp, 2).n(1000)
        print(f"Min f.p.: 2^{log_min_fp:.4f}")

        for log_fp_thrs in [-128, -80, -32]:
            if log_min_fp > log_fp_thrs:
                print(f"  - Var_thrs_{-log_fp_thrs:.0f}: impossible")
                print()
            else:
                log_Var_thrs = find_var_thrs(n, q_prime, N, theta, delta_in, log_fp_thrs)
                Var_thrs = 2^log_Var_thrs
                Gamma, fp = get_fp_pbs(n, q_prime, N, theta, delta_in, Var_thrs)
                print(f"  - Var_thrs_{-log_fp_thrs:.0f}: 2^{log_Var_thrs.n():.3f}, Gamma: {Gamma.n():.3f}, fp: 2^{log(fp, 2).n(1000):.4f}")
                if log(fp, 2).n(1000) > log_fp_thrs:
                    print(f"  - Invalid Var_thrs_{-log_fp_thrs:.0f}")
                    break
                max_depth = ((Var_thrs - Var_KS) / Var_scaled_in).n()
                print(f"  - max-depth: {max_depth:.2f}, (max-depth / log_modulus): {max_depth/log_modulus:.2f}")
                print()
        print()
    print()
    print()



wopbs_2_2_wo_refresh = (
    "wopbs_param_message_2_carry_2_ks_pbs w/o refresh",
    (769, 0.0000043131554647504185^2), # (LWE dim, LWE var)
    (2^15, 2), # (PBS base, PBS level)
    (2^7, 7, 2^36), # (Tr base, Tr level, split base)
    (2^16, 2), # (SS base, SS level)
    (2^4, 4), # (CBS base, CBS level) with increased level
    (2^4, 3), # (KS base, KS level) with increased level
    2, # theta
    4, # log_modulus
    2, # num_extract
)

param_list = [
    wopbs_2_2_wo_refresh,
]

for param in param_list:
    name = param[0]
    (n, Var_LWE) = param[1]
    (B_pbs, l_pbs) = param[2]
    (B_tr, l_tr, b_tr) = param[3]
    (B_ss, l_ss) = param[4]
    (B_cbs, l_cbs) = param[5]
    (B_ksk, l_ksk) = param[6]
    theta = param[7]
    log_modulus = param[8]
    max_num_extract = param[9]

    print(f"========================= {name} =========================")
    print(f"n: {n}, N: {N}, k: {k}, B_pbs: 2^{log(B_pbs, 2)}, l_pbs: {l_pbs}\n")

    Var_pbs = get_var_pbs(N, k, n, q, Var_GLWE, B_pbs, l_pbs)
    Var_fft_pbs = get_var_fft_pbs(N, k, n, B_pbs, l_pbs)
    Var_pbs_tot = Var_pbs + Var_fft_pbs
    Var_ks = get_var_lwe_ks(N, k, q, Var_LWE, B_ksk, l_ksk)

    Var_tr = get_var_tr(N, k, q, Var_GLWE, B_tr, l_tr)
    Var_auto_gadget = get_var_glwe_ks_gadget(N, k, q, B_tr, l_tr)
    Var_auto_key = get_var_glwe_ks_key(N, k, q, Var_GLWE, B_tr, l_tr)
    Var_fft_tr = get_var_fft_tr(N, k, B_tr, l_tr, b_tr)
    Var_tr_tot = Var_tr + Var_fft_tr

    Var_ss = get_var_ss(N, k, q, q^2 * Var_GLWE, B_ss, l_ss)
    Var_fft_ss = get_var_fft_ext_prod(N, k, q, B_ss, l_ss)
    Var_ss_tot = Var_ss + Var_fft_ss

    Var_cbs = Var_pbs_tot + Var_ss_tot + (N/2) * Var_tr_tot
    Var_cbs_additive = Var_pbs_tot + Var_ss_tot
    Var_cbs_amp = (N/2) * Var_tr_tot

    Var_Add = get_var_ext_prod(N, k, q, Var_cbs, B_cbs, l_cbs) + get_var_fft_ext_prod(N, k, q, B_cbs, l_cbs)
    Var_add_gadget = get_var_ext_prod_gadget(N, k, q, B_cbs, l_cbs)
    Var_add_inc = get_var_ext_prod_inc(N, k, Var_cbs, B_cbs, l_cbs)
    Var_fft_add = get_var_fft_ext_prod(N, k, q, B_cbs, l_cbs)

    print(f"Var_ks: 2^{log(Var_ks, 2).n():.4f}")
    print(f"Var_pbs_tot: 2^{log(Var_pbs_tot, 2).n():.4f}")
    print(f"  - Var_pbs: 2^{log(Var_pbs, 2).n():.4f}")
    print(f"  - Var_fft_pbs: 2^{log(Var_fft_pbs, 2).n():.4f}")
    print(f"Var_tr_tot : 2^{log(Var_tr_tot, 2).n():.4f}")
    print(f"  - Var_tr: 2^{log(Var_tr, 2).n():.4f}")
    print(f"     - Var_auto_gadget: 2^{log(Var_auto_gadget, 2).n():.4f}")
    print(f"     - Var_auto_key: 2^{log(Var_auto_key, 2).n():.4f}")
    print(f"  - Var_fft_tr: 2^{log(Var_fft_tr, 2).n():.4f}")
    print(f"Var_ss_tot: 2^{log(Var_ss_tot, 2).n():.4f}")
    print(f"  - Var_ss: 2^{log(Var_ss, 2).n():.4f}")
    print(f"  - Var_fft_ss: 2^{log(Var_fft_ss, 2).n():.4f}")
    print(f"Var_cbs: 2^{log(Var_cbs, 2).n():.4f}")
    print(f"  - Var_cbs_additive: 2^{log(Var_cbs_additive, 2).n():.4f}")
    print(f"  - Var_cbs_amp     : 2^{log(Var_cbs_amp, 2).n():.4f}")
    print(f"Var_Add: 2^{log(Var_Add, 2).n():.4f}")
    print(f"  - Var_Add_gadget: 2^{log(Var_add_gadget, 2).n():.4f}")
    print(f"  - Var_Add_inc   : 2^{log(Var_add_inc, 2).n():.4f}")
    print(f"  - Var_fft_Add   : 2^{log(Var_fft_add, 2).n():.4f}")
    print()
    print(f"(B_cbs, l_cbs): (2^{log(B_cbs, 2)}, {l_cbs})")
    print(f"(B_ksk, l_ksk): (2^{log(B_ksk, 2)}, {l_ksk})")
    print()

    for num_extract in range(1, max_num_extract+1):
        # MV-PBS [CIM19] to extract num_extract bits
        if num_extract == 3:
            Var_Add *= 3
        Var_scaled_in = 2^(2*(log_modulus - num_extract)) * Var_Add

        print("-------------------------------------------------------------------------")
        print(f"# extracting bits: {num_extract} (MV-PBS [CIM19] is used after PBSmanyLUT [CLOT21] for LWEtoLev Conversion)")
        print(f"  - Var_Add: 2^{log(Var_Add, 2).n():.4f}")
        print(f"  - Var_scaled_in: 2^{log(Var_scaled_in, 2).n():.4f}")
        print()

        q_prime = q
        delta_in = 2^(64 - num_extract)
        print(f"Delta_in: 2^{log(delta_in, 2)}")
        print("theta:", theta)
        _, min_fp = get_min_fp_pbs(n, q_prime, N, theta, delta_in)
        log_min_fp = log(min_fp, 2).n(1000)
        print(f"Min f.p.: 2^{log_min_fp:.4f}")

        for log_fp_thrs in [-128, -80, -32]:
            if log_min_fp > log_fp_thrs:
                print(f"  - Var_thrs_{-log_fp_thrs:.0f}: impossible")
                print()
            else:
                log_Var_thrs = find_var_thrs(n, q_prime, N, theta, delta_in, log_fp_thrs)
                Var_thrs = 2^log_Var_thrs
                Gamma, fp = get_fp_pbs(n, q_prime, N, theta, delta_in, Var_thrs)
                print(f"  - Var_thrs_{-log_fp_thrs:.0f}: 2^{log_Var_thrs.n():.3f}, Gamma: {Gamma.n():.3f}, fp: 2^{log(fp, 2).n(1000):.4f}")
                if log(fp, 2).n(1000) > log_fp_thrs:
                    print(f"  - Invalid Var_thrs_{-log_fp_thrs:.0f}")
                    break
                max_depth = ((Var_thrs - Var_ks) / Var_scaled_in).n()
                print(f"  - max-depth: {max_depth:.2f}, (max-depth / log_modulus): {max_depth/log_modulus:.2f}")
                print()
        print()
    print()
    print()



wopbs_3_3_wo_refresh = (
    "wopbs_param_message_3_carry_3_ks_pbs w/o refresh",
    (873, 0.0000006428797112843789^2), # (LWE dim, LWE var)
    (2^11, 3), # (PBS base, PBS level)
    (2^12, 4, 2^41), # (Tr base, Tr level, split fft)
    (2^10, 4), # (SS base, SS level),
    (2^15, 3, 2^44), # (KS_to_large base, KS_to_large level, split fft),
    (2^13, 3, 2^42), # (KS_from_large base, KS_from_large level, split fft),
    (2^5, 4), # (CBS base, CBS level) with increased level
    (2^7, 2), # (KS base, KS level) with increased level
    2, # theta
    6, # log_modulus
    4, # num_extract
)

wopbs_4_4_wo_refresh = (
    "wopbs_param_message_4_carry_4_ks_pbs w/o refresh",
    (953, 0.0000001486733969411098^2), # (LWE dim, LWE var)
    (2^9, 4), # (PBS base, PBS level)
    (2^9, 6, 2^39), # (Tr base, Tr level, split fft)
    (2^10, 4), # (SS base, SS level),
    (2^15, 3, 2^44), # (KS_to_large base, KS_to_large level, split fft),
    (2^10, 4, 2^39), # (KS_to_large base, KS_to_large level, split fft),
    (2^3, 8), # (CBS base, CBS level) with increased level
    (2^7, 2), # (KS base, KS level) with increased level
    3, # theta
    8, # log_modulus
    3, # num_extract
)



param_list = [
    wopbs_3_3_wo_refresh,
    wopbs_4_4_wo_refresh,
]

for param in param_list:
    name = param[0]
    (n, Var_LWE) = param[1]
    (B_pbs, l_pbs) = param[2]
    (B_tr, l_tr, b_tr) = param[3]
    (B_ss, l_ss) = param[4]
    (B_to_large, l_to_large, b_to_large) = param[5]
    (B_from_large, l_from_large, b_from_large) = param[6]
    (B_cbs, l_cbs) = param[7]
    (B_ksk, l_ksk) = param[8]
    theta = param[9]
    log_modulus = param[10]
    max_num_extract = param[11]

    print(f"========================= {name} =========================")
    print(f"n: {n}, N: {N}, k: {k}, B_pbs: 2^{log(B_pbs, 2)}, l_pbs: {l_pbs}\n")

    Var_pbs = get_var_pbs(N, k, n, q, Var_GLWE, B_pbs, l_pbs)
    Var_fft_pbs = get_var_fft_pbs(N, k, n, B_pbs, l_pbs)
    Var_pbs_tot = Var_pbs + Var_fft_pbs
    Var_ks = get_var_lwe_ks(N, k, q, Var_LWE, B_ksk, l_ksk)

    Var_to_large = get_var_glwe_ks(N, k, q, Var_large_GLWE, B_to_large, l_to_large)
    Var_fft_to_large = get_var_fft_glwe_ks(N, k, B_to_large, l_to_large, b_to_large)
    Var_to_large_tot = Var_to_large + Var_fft_to_large

    Var_tr = get_var_tr(N, k_large, q, Var_large_GLWE, B_tr, l_tr)
    Var_auto_gadget = get_var_glwe_ks_gadget(N, k_large, q, B_tr, l_tr)
    Var_auto_key = get_var_glwe_ks_key(N, k_large, q, Var_large_GLWE, B_tr, l_tr)
    Var_fft_tr = get_var_fft_tr(N, k, B_tr, l_tr, b_tr)
    Var_tr_tot = Var_tr + Var_fft_tr

    Var_from_large = get_var_glwe_ks(N, k_large, q, Var_GLWE, B_from_large, l_from_large)
    Var_fft_from_large = get_var_fft_glwe_ks(N, k_large, B_from_large, l_from_large, b_from_large)
    Var_from_large_tot = Var_from_large + Var_fft_from_large

    Var_ss = get_var_ss(N, k, q, q^2 * Var_GLWE, B_ss, l_ss)
    Var_fft_ss = get_var_fft_ext_prod(N, k, q, B_ss, l_ss)
    Var_ss_tot = Var_ss + Var_fft_ss

    Var_cbs = Var_pbs_tot + Var_to_large_tot + Var_ss_tot + (N/2) * (Var_tr_tot + Var_from_large_tot)
    Var_cbs_additive = Var_pbs_tot + Var_to_large_tot + Var_ss_tot
    Var_cbs_amp = (N/2) * (Var_tr_tot + Var_from_large_tot)

    Var_Add = get_var_ext_prod(N, k, q, Var_cbs, B_cbs, l_cbs) + get_var_fft_ext_prod(N, k, q, B_cbs, l_cbs)
    Var_add_gadget = get_var_ext_prod_gadget(N, k, q, B_cbs, l_cbs)
    Var_add_inc = get_var_ext_prod_inc(N, k, Var_cbs, B_cbs, l_cbs)
    Var_fft_add = get_var_fft_ext_prod(N, k, q, B_cbs, l_cbs)

    print(f"Var_ks: 2^{log(Var_ks, 2).n():.4f}")
    print(f"Var_pbs_tot: 2^{log(Var_pbs_tot, 2).n():.4f}")
    print(f"  - Var_pbs: 2^{log(Var_pbs, 2).n():.4f}")
    print(f"  - Var_fft_pbs: 2^{log(Var_fft_pbs, 2).n():.4f}")
    print(f"Var_to_large_tot : 2^{log(Var_to_large_tot, 2).n():.4f}")
    print(f"  - Var_to_large: 2^{log(Var_to_large, 2).n():.4f}")
    print(f"  - Var_fft_to_large: 2^{log(Var_fft_to_large, 2).n():.4f}")
    print(f"Var_tr_tot : 2^{log(Var_tr_tot, 2).n():.4f}")
    print(f"  - Var_tr: 2^{log(Var_tr, 2).n():.4f}")
    print(f"     - Var_auto_gadget: 2^{log(Var_auto_gadget, 2).n():.4f}")
    print(f"     - Var_auto_key: 2^{log(Var_auto_key, 2).n():.4f}")
    print(f"  - Var_fft_tr: 2^{log(Var_fft_tr, 2).n():.4f}")
    print(f"Var_from_large_tot : 2^{log(Var_from_large_tot, 2).n():.4f}")
    print(f"  - Var_from_large: 2^{log(Var_from_large, 2).n():.4f}")
    print(f"  - Var_fft_from_large: 2^{log(Var_fft_from_large, 2).n():.4f}")
    print(f"Var_ss_tot: 2^{log(Var_ss_tot, 2).n():.4f}")
    print(f"  - Var_ss: 2^{log(Var_ss, 2).n():.4f}")
    print(f"  - Var_fft_ss: 2^{log(Var_fft_ss, 2).n():.4f}")
    print(f"Var_cbs: 2^{log(Var_cbs, 2).n():.4f}")
    print(f"  - Var_cbs_additive: 2^{log(Var_cbs_additive, 2).n():.4f}")
    print(f"  - Var_cbs_amp     : 2^{log(Var_cbs_amp, 2).n():.4f}")
    print(f"Var_Add: 2^{log(Var_Add, 2).n():.4f}")
    print(f"  - Var_Add_gadget: 2^{log(Var_add_gadget, 2).n():.4f}")
    print(f"  - Var_Add_inc   : 2^{log(Var_add_inc, 2).n():.4f}")
    print(f"  - Var_fft_Add   : 2^{log(Var_fft_add, 2).n():.4f}")
    print()
    print(f"(B_cbs, l_cbs): (2^{log(B_cbs, 2)}, {l_cbs})")
    print(f"(B_ksk, l_ksk): (2^{log(B_ksk, 2)}, {l_ksk})")
    print()

    for num_extract in range(1, max_num_extract+1):
        # MV-PBS [CIM19] to extract num_extract bits
        if num_extract == 3:
            Var_Add *= 3
        Var_scaled_in = 2^(2*(log_modulus - num_extract)) * Var_Add

        print("-------------------------------------------------------------------------")
        print(f"# extracting bits: {num_extract} (MV-PBS [CIM19] is used after PBSmanyLUT [CLOT21] for LWEtoLev Conversion)")
        print(f"  - Var_Add: 2^{log(Var_Add, 2).n():.4f}")
        print(f"  - Var_scaled_in: 2^{log(Var_scaled_in, 2).n():.4f}")
        print()

        q_prime = q
        delta_in = 2^(64 - num_extract)
        print(f"Delta_in: 2^{log(delta_in, 2)}")
        print("theta:", theta)
        _, min_fp = get_min_fp_pbs(n, q_prime, N, theta, delta_in)
        log_min_fp = log(min_fp, 2).n(1000)
        print(f"Min f.p.: 2^{log_min_fp:.4f}")

        for log_fp_thrs in [-128, -80, -32]:
            if log_min_fp > log_fp_thrs:
                print(f"  - Var_thrs_{-log_fp_thrs:.0f}: impossible")
                print()
            else:
                log_Var_thrs = find_var_thrs(n, q_prime, N, theta, delta_in, log_fp_thrs)
                Var_thrs = 2^log_Var_thrs
                Gamma, fp = get_fp_pbs(n, q_prime, N, theta, delta_in, Var_thrs)
                print(f"  - Var_thrs_{-log_fp_thrs:.0f}: 2^{log_Var_thrs.n():.3f}, Gamma: {Gamma.n():.3f}, fp: 2^{log(fp, 2).n(1000):.4f}")
                if log(fp, 2).n(1000) > log_fp_thrs:
                    print(f"  - Invalid Var_thrs_{-log_fp_thrs:.0f}")
                    break
                max_depth = ((Var_thrs - Var_ks) / Var_scaled_in).n()
                print(f"  - max-depth: {max_depth:.2f}, (max-depth / log_modulus): {max_depth/log_modulus:.2f}")
                print()
        print()
    print()
    print()