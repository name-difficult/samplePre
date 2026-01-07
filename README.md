# G-陷门（G-Trapdoor）生成、原像采样与陷门委托（工程原型）

本仓库实现了论文 **“Trapdoors for Lattices: Simpler, Tighter, Faster, Smaller”** 中若干核心算法的**工程化原型**，用于验证结构关系与端到端流程可跑通，并提供一个可扩展的代码框架，便于后续替换为严格的离散高斯采样器与更严谨的参数选择。

> 重要声明：本实现**不是**生产级密码实现。当前版本在若干关键步骤中**未采用严格的离散高斯采样**，而是用高斯浮点采样 + 四舍五入等方式做了“能跑通/便于调试”的近似替代。详见本文末尾“离散高斯采样位置与简化说明”。

------

## 1. 目标与实现范围

本代码覆盖论文中以下三个核心模块：

1. **Primitive / Gadget 矩阵构造**（论文第 2 节 Primitive Lattices）
    - 构造 gadget 向量 $g=(1,2,4,\ldots,2^{k-1})$ 与 $G := I_n \otimes g^T \in \mathbb{Z}_q^{n\times nk}.$
    - 构造 $\Lambda^\perp(g^T)$ 的一个短基 $S_k$，并用于验证 $G(I_n\otimes S_k)\equiv 0 \pmod q$。
2. **G-陷门生成（GenTrap）**（论文第 3.2 节，算法 1）
    - 实现统计/计算实例化中共同的矩阵形式：$A=[\bar A \mid HG-\bar A R] \pmod q$，当前原型默认取 (H=I)，即：$A=[\bar A \mid G-\bar A R] \pmod q$.
3. **基于陷门的原像采样与陷门委托**（论文第 3.4 与 3.5 节）
    - 用一个工程版 `samplePre` 作为“oracle”，实现类似论文中 **SampleD / coset sampling** 的流程骨架。
    - 用该 oracle 实现 **DelTrap**（算法 3）的 batched 版本：一次对多列综合向量并行求解，输出委托后的 (R')。

------

## 2. 参数与整体流程概述

### 2.1 参数设置（工程示例）

```python
n, k, q, w, bar_m, m = get_secure_param(256, 16)
# q = prevprime(2^k) 选取接近 2^k 的素数（便于实验与 mod 运算）
# w = n*k
# bar_m = 2n
# m = bar_m + w
```

- 论文中常见设定是 $k=\lceil \log_2 q\rceil$，w=nk，$m=\bar m + w$。

### 2.2 端到端流程（main 中的调用关系）

1. **生成 gadget 矩阵与基：**
    - `G = gen_gadget_G(n,k,q)`
    - `S_k = build_Sk_from_G(G,q)`（构造 (S_k)）
2. **生成带陷门的矩阵 (A)：**
    - `A, R, G, S_k = gen_trapdoor_G_trapdoor(...)`
    - 其中 `R` 是 G-陷门（论文定义 1），并满足（当 (H=I)）：$A\begin{bmatrix}R\\I\end{bmatrix} \equiv G \pmod q$.
    - 提供了 `verify_G_trapdoor_H_is_I` 用于检查该等式是否成立。
3. **给定综合向量 (u)，求一个原像 (x)，使 (Ax\equiv u\pmod q)：**
    - `x = samplePre(A,R,u,S_k,q,...)`
    - `verify_preimage(A,x,u,q)` 可用于验证。
4. **陷门委托（DelTrap）：**
    - 随机扩展块 `A1` 后，令 (A'=[A\mid A_1])。
    - 目标是输出 (R')，使得（当 (H'=I)）：$A R' \equiv G - A_1 \pmod q$.
    - `deltrap_HI_batched` 调用 `samplePre_batch` 批量求解多列综合向量，输出 R'。

------

## 3. 代码结构与论文对应关系

### 3.1 Primitive / Gadget（论文第 2 节）

- `gen_gadget_G(n,k,q)`
    - 对应：$G := I_n \otimes g^T$，其中 $g=(1,2,\dots,2^{k-1})$。
- `build_Sk_from_G(G,q)`
    - 对应：给出 $\Lambda^\perp(g^T)$ 的一个显式基 $S_k$。
    - 并配套验证函数：
        - `verify_G_S(G,S_k,q)`：验证 $G\cdot(I_n\otimes S_k)\equiv 0\pmod q$。

> 注：这里 `build_Sk_from_G` 用的是常见的“2 下对角 -1 + 最后一列 q e1”的结构化基写法，用于体现论文中“低维 primitive lattice 可显式给出短基”的思想。

### 3.2 GenTrap：G-陷门生成（论文第 3.2 节，算法 1）

- `gen_trapdoor_G_trapdoor(...)`
    - 对应论文算法 1 的实例（固定 H=I）：
        1. 采样 $\bar A \leftarrow \mathbb{Z}_q^{n\times\bar m}$
        2. 采样 $R \leftarrow D$（论文中 D 应满足次高斯+正则性/或 LWE 形式）
        3. 输出 $A=[\bar A\mid G-\bar A R]\pmod q$
    - 这里的 `R` 采用 `sample_R` 生成：`round(N(0,sigma^2))`，属于工程近似（非严格离散高斯，见末尾说明）。

### 3.3 SampleD / 原像采样骨架（论文第 3.4 节）

- `samplePre(A,R,u,S_k,q,...)`
    - 对应论文 3.4 的“卷积（convolution）思路”的**工程骨架**：
        1. 采样扰动 p（论文里应来自某个离散高斯/随机化取整流程）
        2. 计算 $v=u-Ap\pmod q$
        3. 在 primitive lattice 上构造 z（你这里用 bit-decomposition 得到 $z_0$，再加 $S_k t$）
        4. 计算 $y=\begin{bmatrix}R\\I\end{bmatrix}z$
        5. 输出 x=p+y，从而 $Ax \equiv u\pmod q$

> 注：严格论文版需要确保输出分布接近目标陪集离散高斯分布；当前版本主要保证“代数关系正确、结构可运行”。

### 3.4 DelTrap：陷门委托（论文第 3.5 节，算法 3）

- `deltrap_HI_batched(...)`
    - 对应算法 3（固定 (H'=I) 的特例）：
        - 先计算 $U = G-A_1 \pmod q$
        - 对每一列 $u_j$，调用 oracle 从 $\Lambda^\perp(A)$ 的相应陪集中采样得到列向量 $r'_j$，使 $A r'_j \equiv u_j\pmod q$
    - 这里的 oracle 是 `samplePre_batch`，一次处理 B 列，显著减少 Python 循环开销。

------

## 4. 当前实现的“工程化简化”与安全性声明

本仓库的实现重点是**结构正确与流程可运行**，因此在多处把论文要求的**离散高斯采样**替换为更简单的近似采样：

### 4.1 需要“严格离散高斯采样”的位置（关键点）

以下位置在论文中原则上应使用离散高斯（或等价的随机化取整+拒绝采样）以满足统计距离与安全性证明要求：

1. **GenTrap 中采样陷门 (R)**
    - 论文：$R \leftarrow D_{\mathbb{Z},s}^{\bar m\times w}$（计算实例化时甚至是 LWE 形式，$s=\alpha q$）
    - 代码位置：`sample_R(bar_m,w,sigma)`
    - 当前简化：`R_ij = round(N(0, sigma^2))`
2. **SampleD/卷积采样中的扰动向量 (p)**
    - 论文：(p) 来自某个与目标协方差匹配的离散高斯（或随机化取整生成的离散分布），用于修正倾斜分布并实现“球形化”
    - 代码位置：`sample_p(m,s,...)`
    - 当前简化：`p = round(N(0, (c*s)^2 I))`
3. **构造 (z = z_0 + S_k t) 中的 (t)**
    - 论文：(t)（或等价变量）应来自离散高斯，以保证 primitive lattice 上的采样分布满足 (\eta_\varepsilon) 等条件与统计隐藏性
    - 代码位置：
        - `samplePre` 中的 `t = round(N(0, sigma_t^2))`
        - `samplePre_batch` 中的 `t = round(N(0, sigma_t^2))`
    - 当前简化：同样用 `np.random.normal` + `round`
4. **DelTrap（算法 3）中的 oracle：对陪集采样每一列 (r'_j)**
    - 论文：oracle $\mathcal O$ 必须能对 $\Lambda^\perp(A)$ 的相应陪集进行离散高斯采样（参数 $s'\ge \eta_\varepsilon(\Lambda)$）
    - 代码位置：`deltrap_HI` / `deltrap_HI_batched` 内部调用的 `samplePre` / `samplePre_batch`
    - 当前简化：用工程版 `samplePre` 充当 oracle，仅保证同余方程成立，不保证分布性质

### 4.2 我们当前如何简化实现

当前版本的统一简化策略是：

- 用 `np.random.normal(...)` 生成连续高斯样本
- 用 `np.rint(...)`（四舍五入）映射到整数
- 在需要模 q 的地方做 `mod q`
- 在需要“陪集正确性”的地方，用构造保证 $A x \equiv u\pmod q$

该策略优点是：实现简单、速度快、便于调试验证矩阵关系；缺点是：**不满足论文证明所需的严格分布假设**，不能用于任何安全结论。

------

## 5. 运行方式（示例）

直接运行主程序：

```bash
python samplePre.py
```

你会看到类似输出：

- `[gen_trapdoor_G_trapdoor] ... ms`
- `[samplePre] ... ms`
- `[deltrap_batched] ... ms`

如需额外验证，可打开：

- `verify_G_trapdoor_H_is_I(A, R, G, q)`：验证 $A[R;I]\equiv G$
- `verify_G_S(G, S_k, q)`：验证 $G(I\otimes S_k)\equiv 0$
- `verify_preimage(A, x, u, q)`：验证 $Ax\equiv u$

# Engineering Prototype for G-Trapdoor Generation, Preimage Sampling, and Trapdoor Delegation

This repository provides an **engineering prototype** of several core algorithms from the paper **“Trapdoors for Lattices: Simpler, Tighter, Faster, Smaller.”** It is intended to verify key algebraic relations and ensure the end-to-end pipeline runs correctly, while offering an extensible code framework that can later be upgraded with **rigorous discrete Gaussian samplers** and **more principled parameter selection**.

> Important note: this implementation is **not** production-grade cryptographic software. In several critical steps, it does **not** use strict discrete Gaussian sampling; instead, it adopts a pragmatic approximation (continuous Gaussian sampling in floating point followed by rounding) to keep the prototype simple and easy to debug. See “Where discrete Gaussian sampling is required” near the end.

------

## 1. Goals and Scope

This code covers three main modules from the paper:

1. **Primitive / Gadget matrix construction** (Paper Section 2: *Primitive Lattices*)
    - Construct the gadget vector (g=(1,2,4,\ldots,2^{k-1})) and the gadget matrix
      [
      G := I_n \otimes g^T \in \mathbb{Z}_q^{n\times nk}.
      ]
    - Construct a short basis (S_k) for (\Lambda^\perp(g^T)), and verify
      [
      G,(I_n\otimes S_k)\equiv 0 \pmod q.
      ]
2. **G-trapdoor generation (GenTrap)** (Paper Section 3.2, Algorithm 1)
    - Implement the common matrix form used in both the statistical and computational instantiations:
      [
      A=[\bar A \mid HG-\bar A R] \pmod q.
      ]
      The current prototype fixes (H=I), i.e.,
      [
      A=[\bar A \mid G-\bar A R] \pmod q.
      ]
3. **Trapdoor-based preimage sampling and trapdoor delegation** (Paper Sections 3.4 and 3.5)
    - Use an engineering implementation `samplePre` as an “oracle” to realize a workflow skeleton similar to **SampleD / coset sampling** in the paper.
    - Use this oracle to implement a **batched** version of **DelTrap** (Algorithm 3): solve multiple syndrome columns in parallel and output the delegated trapdoor (R').

------

## 2. Parameters and End-to-End Workflow

### 2.1 Parameter setting (engineering example)

```python
n, k, q, w, bar_m, m = get_secure_param(256, 16)
# q = prevprime(2^k): choose a prime close to 2^k (convenient for experiments and mod arithmetic)
# w = n*k
# bar_m = 2n
# m = bar_m + w
```

- A common setting in the paper is (k=\lceil \log_2 q\rceil), (w=nk), and (m=\bar m + w).

### 2.2 End-to-end workflow (main call graph)

1. **Generate the gadget matrix and basis**
    - `G = gen_gadget_G(n,k,q)`
    - `S_k = build_Sk_from_G(G,q)` (construct (S_k))
2. **Generate a trapdoor-equipped matrix (A)**
    - `A, R, G, S_k = gen_trapdoor_G_trapdoor(...)`
    - Here `R` is a G-trapdoor (Definition 1 in the paper). When (H=I), it satisfies:
      [
      A\begin{bmatrix}R\I\end{bmatrix} \equiv G \pmod q.
      ]
    - `verify_G_trapdoor_H_is_I` is provided to check this relation.
3. **Given a syndrome (u), find a preimage (x) such that (Ax\equiv u\pmod q)**
    - `x = samplePre(A,R,u,S_k,q,...)`
    - `verify_preimage(A,x,u,q)` can be used for validation.
4. **Trapdoor delegation (DelTrap)**
    - After sampling an extension block `A1`, define (A'=[A\mid A_1]).
    - The goal is to output (R') such that (for (H'=I)):
      [
      A R' \equiv G - A_1 \pmod q.
      ]
    - `deltrap_HI_batched` calls `samplePre_batch` to solve multiple syndrome columns in batches and outputs (R').

------

## 3. Code Structure vs. Paper Mapping

### 3.1 Primitive / Gadget (Paper Section 2)

- `gen_gadget_G(n,k,q)`
    - Corresponds to (G := I_n \otimes g^T), where (g=(1,2,\dots,2^{k-1})).
- `build_Sk_from_G(G,q)`
    - Provides an explicit basis (S_k) for (\Lambda^\perp(g^T)).
    - Verification helper:
        - `verify_G_S(G,S_k,q)`: checks (G\cdot(I_n\otimes S_k)\equiv 0\pmod q).

> Note: `build_Sk_from_G` uses a common structured basis form (“2 on the diagonal, -1 on the subdiagonal, and the last column as (q e_1)”), reflecting the paper’s idea that low-dimensional primitive lattices admit explicit short bases.

### 3.2 GenTrap: G-trapdoor generation (Paper Section 3.2, Algorithm 1)

- `gen_trapdoor_G_trapdoor(...)`
    - Implements an instance of Algorithm 1 with (H=I):
        1. Sample (\bar A \leftarrow \mathbb{Z}_q^{n\times\bar m})
        2. Sample (R \leftarrow D) (in the paper, (D) must satisfy subgaussianity + regularity, or have an LWE form)
        3. Output (A=[\bar A\mid G-\bar A R]\pmod q)
    - In this prototype, `R` is generated via `sample_R` as `round(N(0,sigma^2))`, which is an engineering approximation (not a strict discrete Gaussian).

### 3.3 SampleD / preimage-sampling skeleton (Paper Section 3.4)

- `samplePre(A,R,u,S_k,q,...)`
    - Implements an engineering skeleton of the paper’s convolution-based idea:
        1. Sample a perturbation vector (p) (in the paper, this should come from discrete Gaussian / randomized rounding)
        2. Compute (v=u-Ap\pmod q)
        3. Construct (z) over the primitive lattice (here: bit-decomposition gives (z_0), then add (S_k t))
        4. Compute (y=\begin{bmatrix}R\I\end{bmatrix}z)
        5. Output (x=p+y), hence (Ax \equiv u\pmod q)

> Note: the strict paper version requires the output distribution to be statistically close to the target coset discrete Gaussian. The current prototype mainly ensures the **algebraic correctness** and a runnable structure.

### 3.4 DelTrap: trapdoor delegation (Paper Section 3.5, Algorithm 3)

- `deltrap_HI_batched(...)`
    - Implements the special case of Algorithm 3 with (H'=I):
        - Compute (U = G-A_1 \pmod q)
        - For each column (u_j), call an oracle to sample from the corresponding coset of (\Lambda^\perp(A)) to obtain (r'_j) such that (A r'_j \equiv u_j\pmod q)
    - The oracle here is `samplePre_batch`, which processes (B) columns at once and substantially reduces Python-loop overhead.

------

## 4. Engineering Simplifications and Security Disclaimer

This repository prioritizes **structural correctness and end-to-end executability**, so it replaces discrete Gaussian sampling required by the paper with simpler approximations in multiple places.

### 4.1 Where strict discrete Gaussian sampling is required (key points)

In principle, the following points should use **discrete Gaussian sampling** (or equivalent randomized rounding + rejection sampling) to meet the statistical-distance and security-proof requirements in the paper:

1. **Sampling the trapdoor (R) in GenTrap**
    - Paper: (R \leftarrow D_{\mathbb{Z},s}^{\bar m\times w}) (in computational instantiation, it can even be LWE-form with (s=\alpha q))
    - Code: `sample_R(bar_m,w,sigma)`
    - Current simplification: `R_ij = round(N(0, sigma^2))`
2. **Perturbation vector (p) in SampleD / convolution**
    - Paper: (p) should come from a discrete Gaussian (or a distribution derived from randomized rounding) matched to the target covariance, to “sphericalize” the distribution and correct skew
    - Code: `sample_p(m,s,...)`
    - Current simplification: `p = round(N(0, (c*s)^2 I))`
3. **Sampling (t) in (z = z_0 + S_k t)**
    - Paper: (t) (or an equivalent variable) should be discrete Gaussian to ensure the primitive-lattice sampling meets conditions such as (\eta_\varepsilon) and achieves statistical hiding
    - Code:
        - `samplePre`: `t = round(N(0, sigma_t^2))`
        - `samplePre_batch`: `t = round(N(0, sigma_t^2))`
    - Current simplification: again `np.random.normal` + rounding
4. **The DelTrap (Algorithm 3) oracle: sampling each coset column (r'_j)**
    - Paper: oracle (\mathcal O) must sample from the appropriate coset of (\Lambda^\perp(A)) under a discrete Gaussian with parameter (s'\ge \eta_\varepsilon(\Lambda))
    - Code: internal calls to `samplePre` / `samplePre_batch` in `deltrap_HI` / `deltrap_HI_batched`
    - Current simplification: use engineering `samplePre` only to satisfy the congruence, without distribution guarantees

### 4.2 What we do in this prototype

The unified simplification strategy is:

- Use `np.random.normal(...)` to sample continuous Gaussian values
- Use `np.rint(...)` (rounding) to map to integers
- Apply `mod q` where needed
- Use constructions that guarantee (A x \equiv u\pmod q) when coset correctness is required

Advantages: simple implementation, fast execution, easy debugging of matrix relations.
Disadvantages: it **does not satisfy** the distributional assumptions required by the paper’s proofs and therefore supports **no security claims**.

------

## 5. How to run (example)

Run the main script directly:

```bash
python samplePre.py
```

You should see output similar to:

- `[gen_trapdoor_G_trapdoor] ... ms`
- `[samplePre] ... ms`
- `[deltrap_batched] ... ms`

Optional checks you can enable:

- `verify_G_trapdoor_H_is_I(A, R, G, q)`: verify (A[R;I]\equiv G)
- `verify_G_S(G, S_k, q)`: verify (G(I\otimes S_k)\equiv 0)
- `verify_preimage(A, x, u, q)`: verify (Ax\equiv u)