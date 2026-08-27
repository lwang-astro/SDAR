# Hierarchical BLogH: Theory Reference

> **用途**: 本文档提供 BTLogH (g_func=4) 的理论基础与开发历史记录。实际代码实现见 [`hierarchical_blogh_impl_notes.md`](./hierarchical_blogh_impl_notes.md)。
> 原始计划中关于梯度无需修正、`multiplyOuterNodePotentials()` 用 `_bin.semi` 等方案在实测中被发现有问题，
> 已被 `processOuterNode()` 替代。本文档仅保留经验证正确的理论部分。
> 2026-08-26：吸收已完成的《双曲层级 ds 规范统一计划》（原独立文件已删除），见 §3；§2 标注为历史设计。
> §3.5 为 q-cap 之后的束缚段跳变归因与对照矩阵判决（分辨率/重构/LogH 三思路全部证伪，事件型 floor ~1e-6 为当前框架上限）。

---

## 0. BLogH Theory Foundation

### 0.1 What is time-transformed symplectic integration?

The standard equations of motion $\dot{Q}=P/m, \dot{P}=-\nabla U$ are stiff for hierarchical
systems (Kepler orbits with vastly different periods). The **time transformation** trick
introduces a fictitious time $s$:

$$dt = \frac{ds}{g(Q)}, \quad \frac{dQ}{ds} = \frac{P}{g(Q)}, \quad \frac{dP}{ds} = -\frac{\nabla U(Q)}{g(Q)}$$

where $g(Q) > 0$ is the **inverse time transformation function**. If $g \propto 1/r$ (the
pair potential), then $ds \propto dt/r$ — the fictitious time step is uniform along a
Kepler orbit, removing the stiffness. This is the **LogH method** (Mikkola & Aarseth 2002).

The symplectic DKD (Drift-Kick-Drift) integrator splits the Hamiltonian:
- **Drift** (kinetic only): $dQ/ds = P/g$
- **Kick** (potential only): $dP/ds = -\nabla U/g$

The DKD error scales as $\Delta t^2$ and depends on Poisson brackets of the generating
functions. For LogH ($g = \sum U_{ij}$), the error is:

$$H_{\text{DKD,err}} \propto \left(\frac{U}{T}\right)^2$$

This becomes large when $|U| \gg |T|$ (deep potential well), which is exactly the case
for hierarchical systems.

### 0.2 Why product instead of sum? (BLogH = MulPot)

Standard LogH uses $g_i = \sum_{j\neq i} U_{ij}$ (sum of pair potentials). 
In a hierarchical triple/quadruple, $U_{in} \gg U_{out}$, so the sum is dominated by
the inner binary — the outer structure is lost in the time transformation.

BLogH (Binary LogH) uses the **product** instead:

$$g_i^{\text{Mul}} = \prod_{j\neq i} U_{ij}$$

The gradient that drives the drift becomes:

$$\frac{\nabla g_i^{\text{Mul}}}{g_i^{\text{Mul}}} = \sum_{j\neq i} \frac{\nabla U_{ij}}{U_{ij}} = \sum_{j\neq i} \nabla(\ln U_{ij})$$

**Key insight**: Each $\nabla(\ln U_{ij}) = -\hat{r}_{ij}/r_{ij}$ scales as $1/r_{ij}$
REGARDLESS of the magnitude of $U_{ij}$. So inner and outer pairs contribute with
equal weight in the logarithmic gradient — the hierarchical scale disparity is compressed.

### 0.3 Dimensional analysis and ds

| Method | $g$ dimension | $ds$ needed for $dt=ds/g$ to be time |
|--------|--------------|--------------------------------------|
| Standard LogH | energy | $ds \sim$ energy·time |
| BLogH binary-only (K inner pairs) | energy^K | $ds \sim$ energy^K·time |
| BLogH all-pairs (C(N,2) pairs) | energy^{C(N,2)} | $ds \sim$ energy^{C(N,2)}·time |

For a B-B quadruple with binary-only (K=2): $ds = ds_1 \cdot ds_2 / \sqrt{P_1 P_2}$.
For "all" mode (K=6): ds would need ~energy^6·time — dimensionally explosive.

> **实现注意**: 最终 ds 公式已由此处的近似版本升级为基于 $P_{\rm eff,min}$ 的公式
> (paper.tex Eq.~ds\_combined)，详见 impl notes。

### 0.4 Existing hybrid_switch modes

The code in `symplectic_integrator.h` supports these modes for `AR_G_FUNC_MUL_POT`:

| switch | CLI flag | $g$ | gtgrad source | `calc_gt_cross` | Use case |
|--------|----------|-----|---------------|-----------------|----------|
| 0 | `off` | $\sum U_{ij}$ (standard LogH) | all pairs, sum of $\nabla U$ | true | Baseline |
| 1 | `binary` | $\prod_{inner} U_{ij}$ | inner pairs only, sum of $\nabla\ln U$ | false | B-B quadruple |
| 2 | `normal-binary` | $(\prod_{inner} U_{ij})^{1/K}$ | same as 1 | false | Like 1 but g~energy |
| 3 | `all` | $\prod_{all} U_{ij}$ | ALL pairs, sum of $\nabla\ln U$ | true | Small-N systems |
| -1 | `auto` | dynamic 0↔1 | depends | varies | Unstable triple |
| **4** | **`hierarchical`** | **$\prod_{\rm n} U_{\rm n}$ (tree-level)** | **inner + outer nodes, $\nabla\ln U$ per mass fraction** | **false** | **B-B quadruple+, BTLogH** |

Key code variables:
- `gt_kick_inv_.value` — scalar $g$ for the kick step ($dt_{kick} = ds / g$)
- `force_[i].gtgrad[3]` — per-particle drift gradient $\nabla\ln U$ (or $\nabla U$ for switch=0)
- `gt_drift_inv_` — accumulated $g$ for drift, updated via predictor-corrector
- `calc_gt_cross` — whether cross-tree pairs contribute to gtgrad/gt_kick_inv

### 0.5 The "all" mode problem

For N=4, switch=3 computes $g = \prod_6 U_{ij}$. Dimension is $[\text{energy}]^6$.
The 4 cross-pair $U_{ij}$ values are typically < 1 (in G=1 units), making $g_{all} < g_{bin}$
and actually REDUCING the time step — counterproductive. For N>4 the dimensional
explosion makes ds estimation impractical.

---

## 1. 从计划到实现：关键变更

以下总结了 plan 阶段的设计与最终实现之间的差异，避免后续开发重踩旧坑。

### 1.1 梯度修正（原计划最大遗漏）

**Plan (§0.6)**: "梯度计算与 switch=1 完全一致，只需 $U_{\rm out}$ 乘入 gt_kick_inv_.value"

**实际**: 缺少外层节点的 $\nabla\ln U_{\rm n}$ 贡献会导致 DKD drift-kick 格式不一致。
最终用 `processOuterNode()` 一次 tree walk 同时完成 $U_{\rm n}$ 乘入和梯度分发（按质量分数 $\pm m_i/M_{\rm m} \cdot \hat{r}_{\rm s}/r_{\rm s}$ 写到各粒子 `force_[i].gtgrad`）。

### 1.2 距离计算

**Plan (§1.3)**: `multiplyOuterNodePotentials()` 用 `_bin.semi`（半长轴）

**实际**: `_bin.semi` 在积分期间不更新。`processOuterNode()` 用瞬时成员位置 `pos0 - pos1`，
正确反映离心率变化。`_bin.r` 同样过时，不可用。

### 1.3 nbin 不递增

**Plan (§1.3)**: `gt_kick_inv_.nbin++`

**实际**: nbin 仅用于 switch=2 几何平均归一化，switch=4 不读取此值，不递增。

### 1.4 ds 公式

**Plan (§2)**: `ds = ds1*ds2/sqrt(P1*P2) * U_o_init / sdiv`

**实际**: 被 paper.tex Eq.~(ds\_combined) 取代，基于 $P_{\rm eff,min} = \min_k(\kappa_k P_k)$ 和
$N_{\rm s}=32$ substeps/orbit。详见 impl notes。

---

## 2. 自动切换 (4↔0) — 开发计划

> **状态（2026-08-26）**: 本节为历史设计记录，与现行实现有实质差异：
> 1. **`g_func_user==4` 时判据被旁路**——`checkGFuncCriterionIter()` 对 BTLogH 直接 `continue`（"always use"），
>    即 auto 模式下 **4→0 的 LogH 回退从不发生**（ustabtri 全程 g_func=4 实证）。这是刻意设计：
>    BTLogH 的外层节点经 `processOuterNode` 进 g，双曲层级由 q-cap 处理（§3），无需退回 LogH。
> 2. **hysterisis / check_interval 已取消**——现行实现在 `syncTreeSlowDownAndDs` Step 6 每调用点评估，
>    无确认计数。
> 3. §2.3 的 `switchHierarchicalMethod()` 已并入 `syncTreeSlowDownAndDs()`（含 gt_drift_inv_ 重置与 ds 重算）。
> 下文保留原始设计供参考；若未来恢复 4↔0 切换（如 KL 循环场景），需重新评估判据与切换代价。

### 2.1 目标与设计原则

当层次结构因强扰动、共振交互等被破坏时，BTLogH 应从 `hybrid_switch=4` 自动 fallback 到标准 LogH（`hybrid_switch=0`）；当层次结构恢复时再切回。新增 CLI sentinel 值 `-2` 对应此 auto 模式。

**设计原则**：

1. **复用现有判据**：`checkHybridMethodCriterionIter()` 的扰动比判据对 4↔0 同样适用，无需重新设计。
2. **内置于 `integrateToTime()`**（方案 A）：对 standalone AR、Hermite integrator、PeTar 统一生效，调用方无需额外逻辑。
3. **不依赖 tree 重建**：normal integration 期间 binary tree 结构不变，每次检查不需重建 tree。
4. **切换后必须重算 ds**：BTLogH 与 LogH 的 ds 公式不同，切换时必须调用 `calcDsAndStepOption()`。
5. **防止振荡**：加 hysterisis 避免边界情况频繁切换。

### 2.2 判据

`checkHybridMethodCriterionIter()`（`symplectic_integrator.h` 1626行）逐层检查每个 inner binary：

$$\text{pert\_ratio} = \frac{M_{\text{out},1}M_{\text{out},2}}{M_{\text{in},1}M_{\text{in},2}} \cdot \left(\frac{a_{\text{in}}(1+e_{\text{in}})}{a_{\text{out}}(1-e_{\text{out}})}\right)^3$$

- 所有 inner binary 的 `pert_ratio < 1` → 用 BTLogH（switch=4）
- 任意 inner binary 的 `pert_ratio >= 1` → 退回 LogH（switch=0）
- 任意 inner binary 为 hyperbolic（`semi <= 0`）→ 退回 LogH

该判据同时被现有 `switchHybridMethod()`（0↔1）使用，对 4↔0 完全适用——因为 switch=4 和 switch=1 同属 product 类方法，层次结构破坏时同样需要退回 sum 类方法。

### 2.3 新函数：`switchHierarchicalMethod()`

在 `symplectic_integrator.h` 中 `switchHybridMethod()` 旁新增：

```cpp
void switchHierarchicalMethod() {
    int hybrid_switch_bk = hybrid_switch;
    auto& bin_root = info.getBinaryTreeRoot();
    if (checkHybridMethodCriterionIter(bin_root)) hybrid_switch = 4;
    else hybrid_switch = 0;

#ifdef AR_TTL
    if (hybrid_switch_bk != hybrid_switch) {
        // 1. Recalc forces and g (with new hybrid_switch)
        Float gt_kick_inv_bk = gt_kick_inv_.value;
        calcAccPotAndGTKickInv();

        // 2. Adjust gt_drift_inv_: use reset when change is large
        Float dg = gt_kick_inv_.value - gt_kick_inv_bk;
        if (fabs(dg) / std::max(fabs(gt_kick_inv_bk), fabs(gt_kick_inv_.value)) > 1e-3)
            gt_drift_inv_ = gt_kick_inv_.value;
        else
            gt_drift_inv_ += dg;

        // 3. Recalc ds (CRITICAL: BTLogH and LogH ds formulas differ)
        calcDsAndStepOption(
            manager->step.getOrder(),
            manager->interaction.gravitational_constant,
            manager->ds_scale,
            hybrid_switch);
    }
#endif
}
```

**与 `switchHybridMethod()` 的差异**：

| 项目 | `switchHybridMethod()` (0↔1) | `switchHierarchicalMethod()` (4↔0) |
|------|------------------------------|-------------------------------------|
| 目标 switch | 1 或 0 | 4 或 0 |
| gt_drift_inv_ 调整 | 仅 `+= diff` | 变化 > 0.1% 时直接重置（参考 interrupt 处理 2711行策略） |
| ds 重算 | ❌ 不重算（已知缺陷） | ✅ 调用 `calcDsAndStepOption()` |

### 2.4 集成到 `integrateToTime()`

在 `integrateToTime()` 主循环中，`updateSlowDownAndCorrectEnergy()` 之后、`integrateOneStep()` 之前插入检查。

**不每步重建 binary tree**：与 ar.cxx 的 0↔1 auto 模式不同，这里不调用 `generateBinaryTree()`。因为：
- 扰动判据依赖的 `semi`/`ecc`/`m1`/`m2` 通过 `updateBinarySemiEccPeriodIter` 定期更新，不需 tree 重建。
- 层次结构破坏是缓慢渐进的过程，切换不会高频发生。
- 每步重建 tree 对 N 较大的系统性能开销不可接受。

**检查频率**：每 `check_interval` 步检查一次（建议默认 100）。首次在 `initialIntegration()` 中初始确定。

伪代码：

```cpp
// in integrateToTime() while loop, before integrateOneStep():
if (hybrid_switch == -2 && step_count % check_interval == 0) {
    switchHierarchicalMethod();
}
```

**关于 ds 切换**：切换发生的步，用旧 ds 完成当前步；`calcDsAndStepOption()` 设置的新 ds 从下一步开始生效。这避免步中突变。

### 2.5 CLI sentinel 值

在 `sample/AR/ar.cxx` 中添加 `hybrid_switch=-2` 的 CLI 映射：

| switch 值 | CLI flag | 说明 |
|-----------|----------|------|
| -2 | `auto-hierarchical` | 自动在 4（hierarchical）与 0（standard）间切换 |

`ar.cxx` 的 `hybrid_option` 参数扩展：字符串 `"auto-hierarchical"` 映射到 `-2`。

Hermite integrator 中对应的配置传递路径也需支持 `-2`。

### 2.6 初始化处理

`initialIntegration()` 依次执行：

```
generateBinaryTree → calcDsAndStepOption → calcAccPotAndGTKickInv → gt_drift_inv_ = gt_kick_inv_.value → calcAccPotAndGTKickInv
```

当 `hybrid_switch=-2` 时，必须在首次 `calcDsAndStepOption` 和 `calcAccPotAndGTKickInv` 之前确定实际 switch 值（4 或 0）：

```cpp
// In initialIntegration(), before first calcDsAndStepOption:
if (hybrid_switch == -2) {
    auto& bin_root = info.getBinaryTreeRoot();
    if (checkHybridMethodCriterionIter(bin_root)) hybrid_switch = 4;
    else hybrid_switch = 0;
}
```

之后 `calcDsAndStepOption` 和 `calcAccPotAndGTKickInv` 使用确定的 switch 值，`gt_drift_inv_` 初始化一致性得到保证。主循环再通过 `switchHierarchicalMethod()` 监测变化。

### 2.7 Hysteresis（防振荡）

当系统处于层次结构边界时，扰动比可能在 1 附近来回穿越，导致频繁切换。策略：

- **连续确认**：要求连续 `confirm_count` 次（建议 3 次）判据一致才触发切换。
- **状态跟踪变量**：`int hybrid_switch_pending_ = 0` 和 `int switch_confirm_count_ = 0`。
  每次检查时：若目标 switch 与当前不同且与 `hybrid_switch_pending_` 一致 → `confirm_count++`；否则重置 `hybrid_switch_pending_`。`confirm_count == 3` 时执行切换。

### 2.8 风险与缓解

| 风险 | 缓解措施 |
|------|----------|
| `gt_drift_inv_` 量级突变 | 变化 > 0.1% 时直接重置（参考 interrupt 处理的 2711行策略） |
| 切换瞬间力瞬态 | `calcAccPotAndGTKickInv()` 重新计算所有力，`calc_gt_cross` 自动随新 switch 变化 |
| ds 步中突变 | 切换步用旧 ds 完成，新 ds 从下一步起生效 |
| PeTar 下的兼容性 | 由于内置于 `integrateToTime()`，PeTar 无需额外修改；`hybrid_switch=-2` 通过 HermiteIntegrator 传递即可 |
| hyperbolic inner binary | `checkHybridMethodCriterionIter()` 已有处理：`semi<=0` 时返回 false |
| tree 在积分期间不变 | 正常积分中 binary tree 拓扑不重建（仅在 merge/destroy 中断时 `generateBinaryTree()`），不影响判据 |

### 2.9 测试设计要点

以下结论来自 2026-08-04 开发讨论，供后续测试设计参考：

**推荐测试场景**：层级三体 KL 循环。

- 三体 tree 固定（1-2 为 inner binary，3 为 outer），但 $e_{\rm in}$ 在 KL 机制下周期性剧烈振荡
- $e_{\rm in}$ 振荡 → pert_ratio 在 1 附近穿越 → 自然触发 4→0→4 切换
- 可复现、非混沌、不瞬间瓦解

**不推荐的场景**：

- **等质量 free-fall 三体**（如黄金三角）：高度混沌、难以复现；初始 pert_ratio>1 即从 LogH 开始，无法验证 BTLogH→LogH 切换
- **依赖 tree 变化的场景**：`integrateToTime()` 中 tree 仅在 merge/destroy 中断时重建，close encounter 不触发重建——这是设计特性而非缺陷

**测试关键指标**：

1. 能量守恒 dE/E：-2 模式应不差于固定 4，接近固定 0
2. hybrid 列时间序列：验证 pert_ratio>1 时 hybrid 从 4 切到 0
3. $e_{\rm in}$ 与 hybrid 列的相关性

**注意事项**：

- `HIERARCHICAL_CHECK_INTERVAL=100`，`CONFIRM_REQUIRED=3`——对于非常短暂的 pert_ratio>1 窗口可能错过切换，可临时降低用于测试
- 测试数据中列格式包含 `hybrid_flag` 列（`hybrid_switch=4` 特有），需用 Python SDARData `hybrid=True` 读入
### 2.10 实现步骤清单

1. [x] 三变量设计：`g_func` / `g_func_user` / `g_func_switch`（2026-08-05 重构）
2. [x] `switchGFuncAuto()` 泛化为通用自动切换（`g_func_user ↔ 0`）
3. [x] 全局重命名：`hybrid_switch` → `g_func`，`checkHybridMethodCriterionIter` → `checkGFuncCriterionIter`
4. [x] `initialIntegration()`：三变量初始化
5. [x] `integrateToTime()`：`g_func_switch == GFUNC_AUTO` 周期检查 + ds 同步
6. [x] `ar.cxx`：`--g-func` + `--g-func-switch fixed|auto` CLI
7. [x] `hermite_integrator.h`：`hybrid_switch` → `g_func`
8. [x] `tools/ar.py`：`hybrid_flag` → `g_func` 列名
9. [x] 编译验证：所有 variants 零 warning
10. [x] 冒烟测试：fixed + auto 模式均运行正常
11. [ ] 深度测试：层级三体 KL 循环（g_func 切换验证）——注：现行实现对 `g_func_user==4` 旁路判据（§2 状态注），此测试只覆盖 mode 1/3 的 auto 路径；若要测 4↔0 需先恢复判据

---

## 3. 双曲层级 ds/g 规范失配修复（2026-08-26，已完成）

> 本节由已完成的独立计划《双曲层级 ds 规范（gauge）统一》迁移合并（原文件已删除）。

### 3.1 问题与机制

Unstable triple（inner: m=(0.1,0.9), a=1e-3, e=0.9; outer: m=(1,1), a=0.01, e=0.9）在 `--g-func 4 --g-func-switch auto -m orbit` 下发生两次双星交换。交换 #2 的密近相遇本身被正确分辨（ds→4.09e-5，dE~1e-6，15c 的 P_eff_min 锚点有效）；**误差在其后的逃逸尾期累积**：

- 交换后新双星 + 逃逸星构成**双曲根节点**：ds 按近心口径 q=|a|(e−1) 计算，而双曲 (a,e) 是 Kepler 常数 → ds 为公式级常数，tree 不变期间**逐位冻结**（0.0127135 恒定至终点）；
- 运行时 g(t) 的外层因子 G·M_L·M_R/r_esc 按 1/r 单调衰减（r_esc: ~2e-3 → 1.31）；
- **规范失配 → dt=ds/g 无界增长** → 内双星分辨率从 ~40 步/轨道跌至个位数（末端真实 dt≈ds/g 已达 ~1 步/轨道；输出间隔 6.1e-5 封顶掩盖了部分恶化），dE 包络从 1e-6 单调爬到 3.8e-2。

机制的一般表述：束缚层级 U_node 有轨道平均（semi 口径，g 周期回归自校，15b）；**双曲层级不存在轨道平均**——ds 的 U 因子是常数而 g 的对应因子是瞬时 1/r。误差通道归因（Phase 0 脚本 `sample/test/analysis_ustabtri_btlogh.py` 固化）：主导通道是 DKD 分裂误差；尾期生效 κ≡1（κ_org~1e-12 被钳制），κ_max=timescale/period 棘轮 1→53 仅为同源症状。

### 3.2 被否决的方案（记录以免反复）

| 方案 | 结果 | 否决原因 |
|---|---|---|
| ds 侧 r_ref=max(q, r_inst) + semi<0 每步重算 ds（原计划 Fix A/B1） | dE 达标 9.9e-7 | **违背"tree 不变则 ds 冻结"设计原则**（破坏扩展相空间结构与 time symmetry）；实测过保守（尾段分辨率超需 ~50×，+12.8% 步数） |
| g 侧 cap = 2\|a\| | dE 跳至 1.67，H_sd 恶化到 6e-3 | 密近三相舞期间瞬态配对的 \|a\| 可任意小 → cap 使 U 虚胀上百倍，树重建时 g 跳变巨大 → 一次性能量破坏。**教训：g 作为状态函数对树重建必须连续有界** |
| 跳过双曲非叶 U_node（v1） | — | dt=ds/g 少一个能量因子，逃逸时发散更严重；对束缚群内 flyby 直接错误 |
| 破碎后切回 LogH（v1） | — | g_func=4 不回退为刻意设计，LogH 并不比 BTLogH 精确（见 §2 状态注） |

### 3.3 采纳方案：g 函数侧 q-cap

`processOuterNode()`（`symplectic_integrator.h`）对双曲非叶节点（semi<0 且 ecc>1）：

```
r_eff = min(r_sep, q),   q = |a|·(e−1)
U_node = G·m1·m2 / r_eff        （cap 生效时梯度分发同步跳过）
```

- **q 与 ds 的 15b 口径严格同规范** → 分辨率锚定设计值、dt 有界；
- **瞬态 flyby 期间 r≤q 恒成立**（q 即近心距）→ 相遇期 cap 根本不触发，无接合跳变（对比 2|a| 的致命缺陷）；
- 接合点 r=q 处值连续（仅梯度折点，可积）；椭圆节点不受影响（闭合轨道 r 有界）；
- **cap 生效时必须同时抑制该节点的梯度分发**（`grad_active=false`），否则 `gt_drift_inv_`（经 gtgrad 演化）与 g 值不一致，破坏 TTL 扩展哈密顿量守恒；
- **ds 逻辑完全不动**——保持"tree 不变则 ds 冻结"原则。

### 3.4 验证（ustabtri；输出 `/home/lwang/localdata/SDAR_BLogH/fixg_test/`）

| 指标 | 基线 | q-cap |
|---|---|---|
| 终点 \|dE\| | 3.8e-2（max 3.77e-2） | **9.9e-7**（H_sd=-3.5e-9，与相遇期同 level） |
| 交换 #2 前 978 行 | — | **位级一致** |
| logh / btlogh 固定模式 | — | **位级一致**（仅计时列 18–20 不同） |
| Nstep_sum | 774265 | 873238（+12.8%） |

- 尾段 +99.5k 步均匀分布于逃逸段（rebuild 后 50 个输出区间内 11k、其后 88.5k），ds 经误差控制器在 ~10 个区间内从 0.0127 自适应至 0.03593 后冻结；有效分辨率 ~1.6e3 步/轨道（设计值 32 的 ~50×，继承 node_scale 与 P_eff_min 锚点在逃逸构型下的保守化，方向正确但过保守——若需省步数为后续独立课题，受"tree 不变 ds 冻结"原则约束）；
- 生产中逃逸尾由 Hermite `checkBreak` 截断，长尾成本不存在；standalone ar 用 `--break-check` 记录事件（默认仅记录不终止；验证 t=0.0598 记录 r=1.06·r_crit，轨迹零影响）；
- 基线注：原 8 月 15 日基线文件已被 2026-08-26 用户以新二进制重跑覆盖（位级复现 fixg_test 结果），基线数字以本节表格与 Phase 0 脚本输出留档为准；
- 分析脚本读回参数：`SDARData(g_func=True, N_particle=3, slowdown=True, N_sd=2, time_measure=True)`；κ 用 `sd` 列（生效值），`sd_max` 列 = timescale/period 上限。

### 3.5 束缚段跳变分析与对照矩阵判决（2026-08-26，q-cap 之后）

q-cap 解决逃逸尾后，束缚段的能量跳变成为剩余误差源。分析（auto s256 日志）：

**事件时间线**（tree 拓扑变化仅 4 次，用 SD 成员对索引检测）：

| t | 事件 | 当行 dE |
|---|---|---|
| 0.002197 / 0.002319 | 交换 #1 + 复配 | ~1e-11 / +2.4e-9 |
| 0.036438 / 0.036499 | 共振飞掠翻转 (0,1)→(0,2)→(0,1)（churn 仅存活 1 个输出间隔） | ~1e-10 |
| **0.03680** | **大跳变 +2.2e-7**：飞掠把内双星 a 从 5.2e-3 硬化到 5.7e-4（9×），跳变与硬化后内双星极近心同时 | — |

**通道排除**：跳变处 ΔdE_SDC_cum=0（slowdown 修正无关）、生效 κ 无变化、跳变时 ds 恒定（不在重建行）、q-cap 未触发（根束缚）。

**三个候选机制被对照实验逐一判决**（s256/s512 = `--ds-scale` 0.125/0.0625；fixed = 不带 `-m orbit`，无 tree 重构与 ds 重算；数据 2026-08-26，auto/fixed/logh × s256/s512 六组）：

| 判决 | 证据 | 结论 |
|---|---|---|
| **分辨率细化无效** | auto s512 同事件跳变不缩小反变差（交换#1 跳变 ×6 差、终点 6.8e-6 vs 9.9e-7、Nstep ×2.05）；若是 6 阶截断应缩小 64× | **事件型误差 floor**：每次强相互作用的贡献取决于穿越相位（不同分辨率=混沌散射的不同抽样），不可收敛 |
| **翻转抑制方向反了** | fixed 同段包络 2.5e-5 vs auto 2.8e-7（事件后窗口差 90–120×，事件前底子差 2.8e3×）；事件瞬间的 ~1e-7 增量两种模式都在（重构不降低事件本身，但①保持事件间干净底子②防止事件后冻结 ds 对 9× 硬化失配的持续累积） | **tree 重构+ds 重算是保护机制，必须保留**；churn 行自身误差 ~1e-10 无害 |
| **切 LogH 代价惨重** | logh 同段 0.16 vs auto 2.8e-7（差 ~6 个量级）；虽 s512 收敛快（~100×改善）但同成本完全劣势 | ΣU 型 g 被内双星主导，共振段分辨率根本不够；4↔0 切换思路正式排除 |

**fixed s512 终点 dE=218 的机制警示**：base 模式（默认 `-m`，无 `-m orbit`）不调用 `updateBinarySemiEccPeriodIter`，tree semi/ecc 停在初始值（外轨道初始束缚 semi>0）；交换 #2 后系统真实变为双曲，stale semi>0 使 **q-cap 的 `semi<0` 条件永不触发** → g 无界衰减 → dt 爆炸。**standalone base 模式不适用于发生交换/逃逸的场景**；生产路径（Hermite/PeTar 经 `integrateToTime`）有根数更新，不受影响。

**结论**：auto（orbit 模式）实际做了一次有利交换——把 fixed 的可收敛但量级大的截断误差（2.5e-5）压成不可收敛但小 90 倍的事件型 floor（~1e-7–1e-6，结构性上限）。束段终点 ~1e-6、逃逸尾零误差、h4 参考同段 ~0.33——已好 5–6 个量级，**接受现状**。若未来确需突破 1e-7，方向不在步长/结构切换，而在事件时刻的相位精确穿越（近心检测+局部细分到穿越点对齐），属另一量级的工程投入。

---

## 4. 其他后续方向

1. **更深层级测试**: 当前仅验证了 B-B 四星（2 层）与 unstable triple（含双曲逃逸，§3）。对 3+1 四星、5 体等更深
   binary tree，需验证 `processOuterNode()` 的递归梯度分发正确性。

2. **与 PeTar 集成**: BTLogH 当前仅在 standalone SDAR 中可用。PeTar 中使用
   SDAR 处理 close encounters 时，需传递 `hybrid_switch=4` 并确保
   binary tree 结构在 PeTar 的 group 检测后仍然有效。头文件变更（§3 的 q-cap）会传播到 PeTar 构建，
   需重编 + functional smoke。

3. **Slowdown + 自动 ds 调优**: 当前 ds 用 Eq.~(ds\_combined) 一次性估计。
   对于 $\kappa$ 随时间变化的情况（如 KL 循环中 inner binary 周期变化），
   可能需要动态调整 ds。**约束（2026-08-26 实测教训）**：任何"ds 随状态重算"的方案都破坏
   "tree 不变则 ds 冻结"原则与 time symmetry（§3.2 第一行），除非先建立 ds 跳变的
   能量修正框架；§3.4 显示逃逸段分辨率过保守 ~50×，优化应优先考虑 g 侧或
   node_scale/P_eff_min 锚点的口径，而非 ds 时变化。另注意 §3.5 判决：对强相互作用
   事件型误差（~1e-7），ds 细化已被 S512 对照证明无效——ds 调优只影响事件之间的
   累积通道，不触碰事件 floor。

