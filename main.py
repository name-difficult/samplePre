import numpy as np

# =========================================================
# 0) Utility
# =========================================================

def disc_gauss_round(shape, s: float, rng: np.random.Generator) -> np.ndarray:
    """Teaching version: round(N(0, s^2)) -> Z."""
    return np.rint(rng.normal(0.0, s, size=shape)).astype(np.int64)

# =========================================================
# 1) Gadget matrix G = I_n ⊗ g^T, g=(1,2,4,...,2^{k-1})
# =========================================================

def gen_G_Matrix(n: int, k: int) -> np.ndarray:
    g = np.array([1 << i for i in range(k)], dtype=np.int64)  # (k,)
    I = np.eye(n, dtype=np.int64)
    G = np.kron(I, g)  # (n, n*k)
    return G

# =========================================================
# 2) Sample R ~ P^{bar_m x w}, P: 0(1/2), +1(1/4), -1(1/4)
# =========================================================

def sample_R(bar_m: int, w: int, rng: np.random.Generator) -> np.ndarray:
    u = rng.random((bar_m, w))
    R = np.zeros((bar_m, w), dtype=np.int64)
    R[u < 0.25] = -1
    R[(u >= 0.25) & (u < 0.50)] = 1
    return R

# =========================================================
# 3) Deterministic GenTrap (trapGen): pass rng IN, no internal RNG
#    Output A, R, G, H, q s.t. A [R;I] = H G (mod q)
# =========================================================

def trapGen(n: int, k: int, rng: np.random.Generator):
    q = 2 ** k
    w = n * k
    bar_m = n * k
    m = bar_m + w

    G = gen_G_Matrix(n, k)                # (n, w)
    H = np.eye(n, dtype=np.int64)
    HG = G % q                            # H=I 时可简化；保留写法也行

    bar_A = rng.integers(0, q, size=(n, bar_m), dtype=np.int64)

    R = sample_R(bar_m, w, rng)           # (bar_m, w)

    # A = [bar_A | HG - bar_A R] mod q
    AR = (bar_A @ (R % q)) % q
    A_right = (HG - AR) % q
    A = np.hstack([bar_A % q, A_right]).astype(np.int64)

    return A, R, G, H, q

def Sk_for_power_of_two(k: int) -> np.ndarray:
    S = np.zeros((k, k), dtype=np.int64)
    for i in range(k - 1):
        S[i, i] = 2
        S[i, i + 1] = -1
    S[k - 1, k - 1] = 2
    return S

# =========================================================
# 4) Gadget-side preimage sampler for q=2^k
#    z = bitdecomp(v) + S_k^T t  (IMPORTANT: S.T)
# =========================================================
def bitdecomp_u(u: int, k: int) -> np.ndarray:
    return np.array([(u >> i) & 1 for i in range(k)], dtype=np.int64)

def sample_preimage_gadget_block(S: np.ndarray, v_i: int, k: int, s_t: float, rng: np.random.Generator) -> np.ndarray:
    q = 2 ** k
    z0 = bitdecomp_u(int(v_i) % q, k)          # g^T z0 = v_i (as integer)
    t = disc_gauss_round((k,), s_t, rng)
    z = z0 + (S.T @ t)                         # KEY FIX
    return z.astype(np.int64)

def sample_preimage_G(v: np.ndarray, n: int, k: int, s_t: float, rng: np.random.Generator) -> np.ndarray:
    q = 2 ** k
    S = Sk_for_power_of_two(k)
    v = np.asarray(v, dtype=np.int64) % q
    blocks = [sample_preimage_gadget_block(S, int(v[i]), k, s_t, rng) for i in range(n)]
    return np.concatenate(blocks, axis=0).astype(np.int64)  # length w=nk

# =========================================================
# 5) Section 3.4 SampleD (teaching version, deterministic via rng)
# =========================================================

def sampleD_34(A: np.ndarray, R: np.ndarray, G: np.ndarray, H: np.ndarray, q: int,
               u: np.ndarray, s_p: float, s_t: float,
               rng: np.random.Generator) -> np.ndarray:
    A = np.asarray(A, dtype=np.int64)
    R = np.asarray(R, dtype=np.int64)
    G = np.asarray(G, dtype=np.int64)
    H = np.asarray(H, dtype=np.int64) % q
    u = np.asarray(u, dtype=np.int64) % q

    n, m = A.shape
    bar_m, w = R.shape
    assert m == bar_m + w
    assert G.shape == (n, w)
    assert u.shape == (n,)

    # This deterministic version assumes H = I (same as your trapGen)
    if not np.array_equal(H % q, np.eye(n, dtype=np.int64) % q):
        raise NotImplementedError("This deterministic teaching code supports H=I only.")

    I_w = np.eye(w, dtype=np.int64)
    Rbar = np.vstack([R, I_w])               # (m, w)

    # 1) sample p
    p = disc_gauss_round((m,), s_p, rng)

    # 2) v = u - A p (mod q)
    v = (u - (A @ (p % q)) % q) % q

    # 3) sample z with G z = v (mod q)
    k = w // n
    z = sample_preimage_G(v, n=n, k=k, s_t=s_t, rng=rng)

    # 4) y = [R;I] z
    y = (Rbar @ z).astype(np.int64)

    # 5) x = p + y
    x = (p + y).astype(np.int64)

    return x

# =========================================================
# 6) Demo: all randomness is FIXED and REPRODUCIBLE
# =========================================================
import time

if __name__ == "__main__":
    rng = np.random.default_rng()

    n = 200
    k = 12  # q = 2^k

    # ---- trapGen timing ----
    rng_trap = np.random.default_rng(rng)
    A, R, G, H, q = trapGen(n, k, rng_trap)

    w = n * k
    HG = (H @ G) % q
    I_w = np.eye(w, dtype=np.int64)
    Y = np.vstack([R % q, I_w])          # (m, w)
    ok = np.array_equal((A @ Y) % q, HG % q)
    print("trapGen verify A*[R;I] == H*G (mod q):", ok)

    # ---- sampleD_34 timing ----
    rng_sd = np.random.default_rng(rng)

    u = rng_sd.integers(0, q, size=(n,), dtype=np.int64)
    s_p = 2.0
    s_t = 2.0

    x = sampleD_34(A, R, G, H, q, u, s_p=s_p, s_t=s_t, rng=rng_sd)

    # ---- correctness check ----
    Ax_mod_q = (A @ (x % q)) % q
    print("verify A x ≡ u (mod q):", np.array_equal(Ax_mod_q, u % q))

