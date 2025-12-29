# README — Simplified Implementation of “Trapdoors for Lattices: Simpler, Tighter, Faster, Smaller” (Sec. 3.1–3.4)

本代码实现了论文 **Trapdoors for Lattices: Simpler, Tighter, Faster, Smaller** 中第 **3.1–3.2** 节的核心构造（G-trapdoor / GenTrap），并提供了一个与第 **3.4** 节思路一致的**简化版** `SampleD`：用于生成满足同余约束 \(Ax\equiv u\pmod q\) 的原像向量。

------

## 1. 代码完成了什么

### 1.1 生成 gadget 矩阵 (G)

- 取 \(q=2^k\)，令 $g=(1,2,4,\dots,2^{k-1}),\qquad G = I_n\otimes g^T \in \mathbb{Z}_q^{n\times w},; w=nk$.
- 代码：`gen_G_Matrix(n, k)`，返回 `G.shape == (n, n*k)`。

这对应论文第 2 节“primitive matrix / gadget”构造，并是第 3 节陷门生成与采样的基础。

------

### 1.2 生成 G-trapdoor（GenTrap 的一个实例）

对应论文第 3.1 定义的 **G-trapdoor** 与第 3.2/算法 1 的构造：

- 采样 \(\bar A \leftarrow \mathbb{Z}_q^{n\times \bar m}\)（均匀）
- 采样 \(R\leftarrow \mathcal{P}^{\bar m\times w}\)（短矩阵分布 (\mathcal P)：0/±1）
- 输出 $A = [\bar A \mid HG-\bar A R] \pmod q$ 使得 $A\begin{bmatrix}R\\I\end{bmatrix} \equiv HG \pmod q$.

代码：`trapGen(n,k,rng)` 输出 `(A,R,G,H,q)`，并在 main 中验证：

```python
Y = np.vstack([R % q, I_w])
ok = np.array_equal((A @ Y) % q, (H @ G) % q)
```

------

### 1.3 生成满足 \(Ax\equiv u\pmod q\) 的原像（简化版 SampleD）

代码的 `sampleD_34` 实现了第 3.4 节“卷积思路”的**框架**，保证最终输出满足同余方程：

1. 采样 \(p\in \mathbb{Z}^m\)（简化：逐坐标独立采样）
2. 计算 \(v = u - A p \pmod q\)
3. 采样 \(z\in\mathbb{Z}^{w}\) 使 \(Gz\equiv v\pmod q\)（按块做 gadget 预像）
4. 令 \(y = \begin{bmatrix}R\\I\end{bmatrix}z\)
5. 输出 \(x=p+y\)

因为（在本实现里）(H=I) 且 \(A\begin{bmatrix}R\\I\end{bmatrix} \equiv G \pmod q\)，所以： $A x = A(p+[R;I]z) \equiv Ap + Gz \equiv Ap + (u-Ap) \equiv u \pmod q$.

main 中验证：

```python
Ax_mod_q = (A @ (x % q)) % q
np.array_equal(Ax_mod_q, u % q)
```

------

## 2. 与论文不完全一致的地方

下面说明哪些地方**没有按论文完整实现**，而是做了简化，以及这样做的原因。

### 2.1 固定标签 (H=I)，未实现一般 (H)

论文允许任意可逆 \(H\in\mathbb{Z}_q^{n\times n}\)，并在后续算法中涉及 \(H^{-1}\)、\(H^{-t}\) 等。

本实现：

- `trapGen()` 固定 `H = I`
- `sampleD_34()` 对非单位 `H` 直接报错

------

### 2.2 离散高斯采样器被简化为 “连续高斯 + 四舍五入”

论文第 3.4 的采样目标是输出分布接近 $D_{\Lambda^\perp_u(A), s}$ 并且需要满足平滑参数、随机化取整（randomized rounding）、统计距离界等条件。

本实现使用：

```python
np.rint(rng.normal(0, s))
```

作为整数采样器。

**原因：**

- 严格离散高斯（可证明统计距离）需要更复杂的采样器与误差控制（例如更严格的 1D sampler、拒绝采样/精确 CDF、随机化取整等）；
- 本代码把重点放在：陷门构造与同余约束正确、矩阵维度与拼接关系正确。

------

### 2.3 “卷积”校正未按论文协方差公式实现

论文 3.4 的关键之一是：为了让最终输出接近球形高斯，需要让 p 的协方差满足 $\mathrm{Cov}(p) = s^2 I - [R;I]\Sigma_G[R;I]^T$, 并要求该差为正定/半正定，从而实现分布校正。

本实现：

- `p` 直接用标量 `s_p` 做逐坐标独立采样；
- 没有显式构造上述协方差，也没有做 PSD 检查与矩阵采样。

------

### 2.4 gadget 侧预像采样实现为 “bit decomposition + 格向量扰动”的简化形式

本实现对每个 \(v_i\)：

- 先构造一个特解 \(z_0\in{0,1}^k) 使 (g^Tz_0=v_i\)
- 再加 \(S^T t\)（其中 (t) 由简化采样器生成）以得到更多随机性

这保证了同余约束，但并没有完整覆盖论文对 oracle/参数选择（如 (\tilde S=2I)、平滑要求、统计距离等）的全部条件与证明要求。

------

## 3. 如何运行

```bash
python main.py
```

预期输出包含两条验证：

1. `trapGen verify A*[R;I] == H*G (mod q): True`
   说明 GenTrap 输出满足 G-trapdoor 关系。
2. `verify A x ≡ u (mod q): True`
   说明 `sampleD_34` 输出的 (x) 满足 (Ax\equiv u\pmod q)。

------

## 4. 参数与维度（当前 demo）

- `n = 200`
- `k = 12` → `q = 2^k = 4096`
- `w = n*k = 2400`
- `bar_m = n*k = 2400`
- `m = bar_m + w = 4800`

采样参数：

- `s_p = 2.0`（用于 `p`）
- `s_t = 2.0`（用于 gadget 扰动 `t`）

这些参数用于功能验证与演示流程，并不直接对应论文中的安全参数选取。

------

## 5. 安全性声明

本实现用于：

- 复现论文第 3.1–3.2 的核心陷门矩阵关系；
- 跑通第 3.4 的流程骨架并验证同余正确性。

本实现不应直接用于生产级密码系统，因为：

- 采样器与卷积校正未严格实现论文的分布要求；
- 参数选择与安全性证明条件未在代码中体现与验证。