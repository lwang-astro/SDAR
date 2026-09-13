# Hierarchical BLogH: Theory Reference

> **用途**: 本文档提供 BTLogH（历史运行期编号 g_func=4；2026-08-27 宏体系重构后为独立编译宏 `AR_G_FUNC_BTLOGH`、运行期选项 `--g-func 1`，重构记录见 [`hierarchical_blogh_impl_notes.md`](./hierarchical_blogh_impl_notes.md) 的"g-func 重构记录"节）的理论基础与开发历史记录。实际代码实现见 [`hierarchical_blogh_impl_notes.md`](./hierarchical_blogh_impl_notes.md)。
> 原始计划中关于梯度无需修正、`multiplyOuterNodePotentials()` 用 `_bin.semi` 等方案在实测中被发现有问题，
> 已被 `processOuterNode()` 替代。本文档仅保留经验证正确的理论部分。
> 2026-08-26：吸收已完成的《双曲层级 ds 规范统一计划》（原独立文件已删除），见 §3；§2 标注为历史设计。
> §3.5 为 q-cap 之后的束缚段跳变归因与对照矩阵判决（分辨率/重构/LogH 三思路全部证伪，事件型 floor ~1e-6 为当前框架上限）。
> 2026-09-12：§5.6 时间对称性分析裁定 R1 不可对称化，主路线改为 g 侧双侧 gauge clamp（方案 G）；§5.7 为其实施计划。

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


## 5. 逃逸段 ds 过保守：P_eff,min 取双曲交会时标 —— 问题定位与修复计划（2026-09-12）

> **一句话**：ejection 后的 escape 段，BTLogH-Adjust 每剩余双星轨道走 ~448 步（设计 N_s=128，
> 超出 ×3.5）。原因**不是** q-cap 与 ds 失配（两者精确抵消），也**不是**扰动太强
> （κ=1、χ≈1），而是 `calcEffectivePeriod` 给双曲节点的"交会时标"
> T_h = 2π|a|^{3/2}/√(G·M) 抢占了 P_eff,min —— 一个**已经结束的交会**仍在
> 要求分辨率。本节给出已验证的机制代数、核心数据、修复方向与新 session 的资源地图。

### 5.1 机制（2026-09-12 从源码 + 日志定量验证）

ejection 后树结构 = [剩余双星 (elliptic leaf)] + [逃逸星 (hyperbolic root)]。

**ds 侧**（`src/AR/information.h`）：
- `calcBLogHDsIter`: ds = Π(ds_i) · P_eff,min / Π(P_eff) · Π_nodes(U^g·χ) · (ds_scale/32)；
  P_eff 候选来自 `calcEffectivePeriod`（leaves **和** internal nodes 都参与）：
  elliptic → P·κ；hyperbolic → **T_h = 2π|a_h|^{3/2}/√(G(m1+m2))**。
- `multiplyDsByNodePotentials`: elliptic node U^g = G m1 m2/a（semi 轨道平均）；
  hyperbolic node U^g = **G m1 m2/q**（近心点，保守上界）。

**g 侧**（`symplectic_integrator.h::processOuterNode` q-cap）：hyperbolic 节点
r > q 时因子冻结在 G m1 m2/q（梯度同步抑制，TTL 一致）。

**代数**：cap 激活时 dt = ds·|U_b(r)|·(G M_b m_e/q)，U_o 因子与 ds 侧 q-gauge
**精确相消** ⇒ 每 P_eff,min 步数恒等于 32/ds_scale。于是

    每 P_r 步数 = N_s · P_r / min(P_r, T_h)     （T_h < P_r 时超出 N_s）

T_h 保护的是交会**逼近段**（plunging 是最快运动，必须分辨）——设计正确；
问题在 outgoing 段：交会已结束、r 单调增长，但冻结的 gauge 让这个"已完成交会"
的分辨率合同一直付到运行结束（同一 tree → 同一 ds；adjust 实际 gauge-frozen）。

### 5.2 核心数据（ustabtri S128，v2.1 二进制 2026-09-02 12:15）

| 量 | 值 | 说明 |
|---|---|---|
| 剩余双星 (0.9+1.0) | a=2.52e-3, e=0.76, P_r=2.92 P_in0 | 交换后形成（比初始内双星更宽） |
| 逃逸星 (0.1) 双曲根数 | \|a_h\|=1.146e-3, e=1.087, q=1.0e-4 | q = 深交会距离 = 0.1 a_in |
| T_h | 0.297·P_r (1.72e-4) | < P_r ⇒ P_eff,min = T_h 分支 |
| committed ds | 0.8607（冻结） | = 公式 0.909·χ_leaf(≈0.947)；若 P_r 分支应为 3.06 → 排除 |
| 实测步率（轨道平均） | **448 步/P_r** | 预测 128·P_r/T_h·(1/χ) ≈ 453 ✓ |
| 扰动因子 | κ=1, χ_root=1（guard）, χ_leaf≈0.95 | **与扰动强度无关** |
| LogH（固定 ds） | 690 步/P_r | = ds_orbit(remnant)/ds_fixed=699 预测 ✓（remnant m1m2 大 10×） |
| H4 | **AR 步率 = 0**（整个 escape 段） | 全部 3,313 步在 ejection 前；group break 后从未重建 |
| Adjust 总步数拆分 | 28,704 = 24,182（t<3.7 P_out）+ 4,522（escape 段） | escape 段占 16% |

**H4 步数 ≪ Adjust 的主因是结构性的**：ejection 后 H4 的 (0.9,1.0) pair
（近心点 6e-4）不再重建 AR group（`checkNewGroup` 的 κ_org>1e-2 等条件），
由 Hermite block steps 处理（dt≤P_in0/4，每次近心约十步）；Adjust 则全程以
含 T_h 缩放的 ds 积分完整三体链。按原始 P_in0 归一：Adjust 153.6/P_in0
（1.2× floor，温和）、LogH 236.5、H4 29.6 —— ×3.5 只在按**当前**剩余轨道归一时出现。

**一般标度**：超出因子 = P_r/T_h = P_r·√(G·M_tot)/(2π|a_h|^{3/2})。
快速逃逸（e≫1，|a_h|→q）最严重（可达 ≫10×）；勉强逃逸（e≈1+，|a_h|≫q）
则 T_h>P_r → 无超出（每轨道恰好 N_s）。

### 5.3 修复方向（按推荐排序；全部是 ds 侧改动，g 侧勿动 —— variant C 教训）

> **状态（2026-09-12）**：本节方向已被 §5.6 的时间对称性分析修订——R1 的退行判据（dr/dt）或刷新阀均为反转奇/历史依赖，原理上不可对称化；主路线改为 g 侧双侧 gauge clamp（方案 G，只动 `processOuterNode`，ds 机制零改动），R1 降级为备选，R2 的 T_live 思想被吸收进 X 的定义。实施计划见 §5.7。

**R1（推荐）双曲层级"退役"（timescale 候选移除）**
过近心且正在退行时（r>q **且** dr/dt>0，或更保守 r>f·q），把该 hyperbolic level
从 P_eff,min 候选中剔除。U^g gauge **保留**在 ds 里（q-cap 保证 g 一致 → 相消仍成立），
只退役时间尺度 ⇒ P_eff,min 回落到 P_r，ds 增大 P_r/T_h 倍，escape 段步数 → N_s。
- 代码锚点：`calcBLogHDsIter` internal-node 分支（information.h :123 附近），
  在 `if (P_node < _P_eff_min)` 前加退行判据（member 位置/速度现成，参考
  `calcPertScale` hyperbolic 分支的取法）。
- 风险：若逃逸星返回（束缚）或另一成员逼近，tree 重构/判据翻转 —— 需连续性论证；
  dr/dt<0（逼近）时判据天然不触发，交会段保护不受影响。

**R2 双曲层级实时时标（R1 的连续版）**
r>q 时用 T_live = 2π·r^{3/2}/√(G·M)（或 r/v_r）替代轨道拟合的 T_h ⇒ P_eff,min
随退行平滑增长，ds 平滑松弛。**注意**：这使 ds 在 tree epoch 内变为时变 ——
与"tree 不变则 ds 冻结"的时间对称性论述冲突（§3.2）；但 slowdown 的 κ 在
epoch 内同样漂移且已被接受（P_eff=P·κ 本来就让 ds 依赖 κ），可援引同一先例。
必须用 Γ/TTL 监控 + 能量误差对照验证（先例：Fix-2 ceiling 的 freeze-then-decay 病）。

**R3 层级从 tree 退役（coupling 判据）**
B_ij = m_i m_j/r³ 低于内层级某阈值时移出 tree（类比 H4 group release；
tree monitoring 已算 pert 量）。改动最大（拓扑变化走重构机制），且 PeTar 里
host 会接管逃逸星 —— standalone 场景优先级低，cluster 场景再考虑。

**R4（不推荐）放宽 one-way cost valve**
现阀 C_v=30 不触发（3.5× < 30×）。调低 C_v 是创可贴：会误伤**有意的**近心
分辨率集中。仅当 R1/R2 均不可行时作为兜底。

**实现约束（历史教训，务必遵守）**：
1. 只改 ds 侧；g 侧对 elliptic/transient 节点加 cap 会破坏 TTL（2026-08-29
   variant C：de_max 退化到 2.36e-3）。q-cap 能存活正因为它只在弱逃逸尾生效。
2. 静默门（quiescence gate v2）与单向阀（v2.1）在
   `syncTreeSlowDownAndDs` 中与估计器交互 —— 任何估计器改动后需全回归。
3. 二进制先装 scratch 名（如 `ar.btlogh.ttl.sd.cm.r1`），回归通过再考虑覆盖
   ~/bin（`make install` 会覆盖标准名，历史日志将不可复现）。

### 5.4 验证协议（新 session）

1. 重建 btlogh 二进制 → scratch 名安装。
2. 回归四套：`btlogh_tri`（稳定三体，要求 value-identical —— 估计器改动不触界）、
   `quad`、`ustabtri`（escape 段步率 → ~128/P_r；de 不劣化）、`ustabquin`
   （burst 窗口 de_max 不劣化 —— R1/R2 不得在逼近段退役时标）+ H4 对照。
3. 重跑验证脚本（持久副本：`localdata/SDAR_BLogH/check_ustabtri_escape.py`），
   `pred/P_r` 列应收敛到 N_s（当前 448/128）。
4. 论文影响：§4.5（freezing/regeneration 小节）+ ustabtri 表格数字需更新；
   同步更新 NOTES。

### 5.5 新 session 资源地图

| 资源 | 位置 |
|---|---|
| Notebook（分析 cell + reader 模式） | `/home/lwang/src/python/SDAR_blogh_method.ipynb`，节 "Automatic ds estimator after tree changes"（ustabtri） |
| Reader kwargs（AR） | `sdar.SDARData(N_particle=3, slowdown=True, g_func=True, N_sd=2, time_measure=True, float_type=np.float64)` |
| Reader kwargs（H4） | 先过滤 `#` 开头行；`sdar.HermiteData(N_particle=3, N_sd=1, time_measure=True)`；时间用 `time+time_offset` |
| 历史笔记（完整机制链+历次实验） | `/home/lwang/src/python/SDAR_blogh_method_NOTES.md`：§ "TODO-session setup: can BTLogH's ds relax after an ejection?" 与 § "2026-09-12: mechanism re-verified" |
| 测试数据 | `/home/lwang/localdata/SDAR_BLogH/ustabtri.{logh,btlogh,btlogh_auto,h4}.s{32,64,128}.{log,err}` |
| IC 参数 | inner 0.1+0.9, a=1e-3, e=0.9；outer 1.0+1.0, a=0.01, e=0.95（apocenter 起）；P_in0=1.987e-4, P_out=4.443e-3；tend=5 P_out；tout=P_in0/16 |
| 二进制 | `~/bin/ar.btlogh.ttl.sd.cm`（v2.1, 2026-09-02 12:15）、`~/bin/ar.blogh.ttl.sd.cm`（LogH/--g-func 0）、`~/bin/hermite` |
| SDAR 代码 | `/home/lwang/code/SDAR`，branch `experiment` @598abd9 + worktree 未提交 `information.h`（= 当前 v2.1 状态） |
| 关键函数 | `information.h`: `calcEffectivePeriod`(:91)、`calcBLogHDsIter`(:111)、`multiplyDsByNodePotentials`(:169)、`calcPertScale`(:331)、`calcDsAndStepOption`(:435)；`symplectic_integrator.h`: `processOuterNode`(:975, q-cap)、`syncTreeSlowDownAndDs`(gate/valve) |
| 验证脚本（持久副本） | `/home/lwang/localdata/SDAR_BLogH/check_ustabtri_escape.py`（逐输出行诊断表）与 `check_ustabtri_escape2.py`（轨道平均步率 + κ） |
| 事件时标参考 | encounter t≈0.5 P_out（首次外近心）；ejection t≈3.7 P_out；escape 段 3.7–5 P_out |

### 5.6 时间对称性分析与方案 G：双侧 gauge clamp（2026-09-12）

> **缘起**：质疑——BTLogH 在 ds 恒定时 time symmetric，而现行 ds 切换机制（regen 提交冻结 + 静默门 + 单向阀）已破坏该性质；若从 escape 状态反积回 t=0，ds 按 escape 状态定，当前问题即不会出现。本节将该直觉形式化并导出主路线修订。结论：**估计器输入本身是反转偶的，不对称性全部住在提交/触发逻辑里**；R1 无法对称化（原理性），R2 方向对但破坏 epoch 不变量；采纳**方案 G**（g 侧双侧 clamp）。

#### 5.6.1 对称性框架：什么必须反转偶

扩展相空间（TTL/LogH + 恒定 ds）的时间对称性要求三个组件都是动量反转（$v\to-v$）下的偶函数（或常数）：

| 组件 | 现状 | 判定 |
|---|---|---|
| $g(Q)$（含 q-cap） | $r>q$ 钉在 q-gauge，纯位置函数 | ✅ 偶 |
| 估计器输入（fits、$T_h$、$q$、semi/ecc） | Kepler 根数在 $v\to-v$ 下不变 | ✅ 偶 |
| ds 提交值 | regen 时刻由状态算出，epoch 内冻结 | ⚠️ 值偶，"提交时机"是历史 |
| quiescence gate | 比较连续两次候选（跨 epoch 历史） | ❌ 依赖方向 |
| 单向阀 | 只升不降（显式方向性） | ❌ |
| R1 的 $dr/dt>0$ | 速度投影 | ❌ 奇 |

#### 5.6.2 反向积分判据的修正

对"从 escape 状态反积回 t=0"的推演给出一个关键修正：反向初始化会重新做 Kepler 拟合，而双曲根数 $(a_h, e)$ 本身反转不变——$T_h$ 是拟合常数，不随 $r$ 衰减，所以反向运行同样把 $P_{eff,min}$ 定为 $T_h$，escape 段照样 448 步/轨道。问题比"提交时机不对"更深：**$T_h(q)$ 是拟合预言的近心时标，而近心穿越物理上只属于相邻两个 epoch 中的一个**（本例属 dance epoch，不属 post-ejection epoch）。resolution 需求是 epoch 局域且方向无关的（同一 epoch 正走反走包含同一状态集），正确的对称对象是：

> **epoch-realized 契约**——ds 作为 epoch 不变量，其分辨率合同应覆盖 epoch 内实际实现的状态：椭圆层级每轨道实现全部相位 → 轨道平均 $P\cdot\kappa$（现行为正确）；双曲层级 epoch 内只实现单边段 $r\in[r_{epoch\,min},\to]$ → 合同应是 realized 活时标 $T_{live}(r)=2\pi r^{3/2}/\sqrt{\mu}$，而非拟合近心 $T_h(q)$。

#### 5.6.3 R1/R2 审判与方案 G

- **R1 无法对称化**（原理性）：区分"近心已过"必须用 $dr/dt$（奇）或刷新阀（历史）；且 R1 依赖 ejection 时刻 regen 的提交路径，静默门可能 DEFER（候选跳变 ×3.37 恰在 ×3 边界附近）而 escape 段树稳定、无后续 regen。
- **R2 的 $T_{live}$ 是偶函数**，方向正确，但把时变性引入 ds，破坏 epoch 不变量（§3.2/Skeel 变步长：无单一影子哈密顿量，扩展能漂移）。
- **出路**：松弛放进 $g$ 侧——$g$ 本来就被设计为逐步重估的偶位置函数（§3.3 的 q-cap 已开先例，当时只加了下界封住 dt 发散；缺的是上界以释放过度分辨率）。

**方案 G**：对双曲非叶节点，把 q-cap 推广为双侧 clamp：

$$r_{eff} = \mathrm{clamp}(r,\; q,\; X), \qquad X = q\cdot\frac{P_{r,eff}}{T_h}$$

$P_{r,eff}$ = 全树最快椭圆层级的 $P\cdot\kappa$（即剔除双曲候选后的 $P_{eff,min}$），$X$ 全由 epoch 常数构成。代数（沿用 §5.1 相消框架）：

- $r\ge X$（钉在上沿）：$\text{steps}/P_r = N_s\cdot\frac{P_r}{T_h}\cdot\frac{q}{X} = N_s$ 精确成立；
- 带内 $q<r<X$（真 $1/r$ 势）：$\text{steps}/P_r = N_s\cdot\frac{P_r}{T_h}\cdot\frac{q}{r} \ge N_s$；
- $r\le q$：与现行完全相同（plunge 区域零行为变化）。

性质：(1) ds 机制零改动，gate/valve 与方案正交；(2) 点wise 不劣于现行——$r\le q$ 全同，带内只在现行可证过度分辨率处省步，双星每有效轨道步数全域 $\ge N_s$，逃逸星 per-$T_{live}$ 下界与现行相同；(3) 退化自洽：$e\to1^+$（$T_h>P_r$）$\Rightarrow X<q$，clamp 退化为现行 q-cap（无问题处零行为变化）；$e\gg1$ 带最宽（正是最浪费情形）；(4) 完全时间对称：$g(Q)$ 偶、ds 常数、零新增触发器，反向积分遍历同一 clamp 带。

| | R1 | R2 | 方案 G |
|---|---|---|---|
| ds 保持 epoch 不变量 | ✅ | ❌ | ✅ |
| 触发/判据全反转偶 | ❌ | 部分 | ✅（零新增触发） |
| escape 段 → $N_s$ | ✅ | ✅ | ✅ |
| 需动 gate/valve | 是（有 DEFER 死锁风险） | 是 | 否 |

#### 5.6.4 残余不对称与二阶段路线

gate/valve 的不对称性作为遗留保留（新机制下 escape 尾不再依赖阀门触发——cap 激活时 g 冻结使 n_live 恒为设计值，阀门本就不触发）；二阶段可对称化：固定 s-格点上的双向刷新，替换"单向 + 历史 gate"。

#### 5.6.5 对称性验证实验（反向积分协议）

把缘起中的思想实验变成可测量协议（详见 §5.7.2 与 Phase 4）：取正向末态、速度取反、同机制反积回 t=0，测量 (a) escape 段步率双向均 $\approx N_s$；(b) dE 包络镜像；(c) 可逆性——反向末行 vs 正向首行的逐粒子偏差与 $|\Delta E|/|E_0|$。对 baseline 与 G 各跑一对，量化 G 移除多少不对称。已知边界：反向穿过 resonant encounter 段含 Lyapunov 放大，属被量化的对象本身，不作 bit 级期待。

**路线裁定**：方案 G 为主路线；R1 降级为备选；R2 的 $T_{live}$ 思想吸收进 $X$ 的定义。

### 5.7 方案 G 实施计划（2026-09-12）

> 执行 §5.6 裁定的方案 G。回归矩阵已按用户裁定排除 btlogh_tri、quad 等全束缚测试（与本问题无关）；无 clamp 路径的 bit-identity 由 R1（ustabtri 冻结树）与 R3（LogH 宏门）覆盖。

#### 5.7.1 代码改动点（唯一文件：`src/AR/symplectic_integrator.h`）

**1.1 `processOuterNode`（:975 起）— q-cap 块（现 :1008–1017）改为双侧 clamp**

- 签名加参：`void processOuterNode(AR::BinaryTree<Tparticle>& _bin, const Float _P_r_eff_min)`；递归调用（:1053）与 `calcAccPotAndGTKickInv` 调用点（:1109）同步传参。
- 逻辑（hyperbolic 分支 `semi<0 && ecc>1` 内）：
  - `q = (-semi)·(ecc−1)`（不变）；
  - `T_h = 2π·|semi|^{3/2}/√(G·(m1+m2))`（G 已在函数内取得；与 `calcEffectivePeriod` 双曲分支同式）；
  - `X = (_P_r_eff_min < NUMERIC_FLOAT_MAX && T_h>0 && _P_r_eff_min > T_h) ? q·_P_r_eff_min/T_h : q`；
  - `r_eff = min(r_sep, X)`（仅当 r_sep>q 生效；r≤q 仍为真 1/r）；
  - 梯度规则不变式：r_sep ≤ X → 激活（in-band 真 1/r 势，与未 cap 一致）；r_sep > X → 抑制（冻结因子无梯度，TTL 一致，同今天 r>q）。
- 值连续：r=q 与 r=X 处仅梯度折点（同今天 r=q 折点，可积）。代数（已验证）：r=X 钉住时 steps/P_r = N_s·(P_r/T_h)·(q/X) = N_s；in-band (q<r<X) 为 N_s·(P_r/T_h)·(q/r) ≥ N_s；T_h ≥ P_r_eff（勉强逃逸）→ X=q → 与现行为逐点相同。
- 更新函数头注释块（记录 Scheme G 代数、kinks、TTL 规则、回退条件）。

**1.2 P_r_eff 预扫 — `calcAccPotAndGTKickInv` 内、`#ifdef AR_G_FUNC_BTLOGH if (g_func_on)` 块（:1107–1111）中、调用 `processOuterNode` 之前**

- 一次 O(N) 遍历 `info.binarytree` 平表：`m1>0 && m2>0 && semi>0 && period>0` → `slowdown.getEffectivePeriod()`（=P·κ，live κ），取 min（模式同 `syncTreeSlowDownAndDs` Step 5 的 flat-list 遍历，但含全部 elliptic 层级）。
- 无 elliptic 层级（全双曲瞬态树）→ 保持 `NUMERIC_FLOAT_MAX` → X 退化为 q（现行为，硬约束 3 的 fallback）。
- 成本：每次力重算 O(N_particle)（三体 ~3 次比较），无分配、无 I/O 格式影响；相对 O(N²) 相互作用可忽略。`processOuterNode` 本就在每次 kick 的力重算路径上。
- 语义：P_r_eff 即"剔除 hyperbolic 候选后的 P_eff,min"——正是 ds 设计分辨率所对应的量；κ 的 epoch 内漂移沿 ds 侧先例（P_eff=P·κ 本就时变）。

**1.3 不改动清单（硬约束）**

- `information.h`（calcBLogHDsIter / multiplyDsByNodePotentials / calcEffectivePeriod / gate / valve / `calcDsAndStepOption`）：零改动（其中 2026-08-26 注释提到的 "r cap in processOuterNode" 语义仍成立，不必动）。
- `BinarySlowDown` 二进制/ASCII I/O、`slow_down.h`：零改动（P_r_eff 瞬态计算，不入快照）。
- 宏门：改动全部在 `AR_G_FUNC_BTLOGH` + `g_func_on` 内 → `--g-func 0` 与 `ar.blogh`/`ar.logh` 构建 bit 不变。

#### 5.7.2 反向积分（backward）对称性验证协议（新）

机制（已从 `ar.cxx` 核实）：
1. **正向**：scratch 二进制跑 ustabtri S128 Adjust（`-m orbit`），加 `--print-precision 17`（无 `-f`——`-f` 会关闭 stdout ASCII 列输出，`ar.cxx:643–648` 的 else 分支）；stdout 每个输出间隔打印完整列（含粒子 mass/pos/vel/radius/id，`printColumnAscii` → `particles.printColumnAscii`）。
2. **末态提取**：`sdar.SDARData(N_particle=3, slowdown=True, g_func=True, N_sd=2, time_measure=True, float_type=np.float64).loadtxt(log, skiprows=1)`；取末行 `snap.particles[-1]` 的 mass/pos/vel（不用 `.last` 二进制 dump，避免 reader 格式风险；末行时间 `snap.time[-1]`）。
3. **反向 IC**：新文件（副本，绝不改原文）`ustabtri.back`：首行 `3`；3 行 `mass x y z −vx −vy −vz radius`（radius 抄原值；CM-frame 打印无碍——初始化自行平移到 CM 系，反向态总动量 ≈0）。
4. **反向运行**：同二进制、同 flags（`--g-func 1 -m orbit --ds-scale 0.25 --slowdown-ref 1e-20 --slowdown-timescale-max 1e10 -o 6.103515625e-05 --print-precision 17`），`-t` = 正向末行时间；tree/ds 估计器从反向态重新起步（same mechanism）。注意 `.last` 会写到 IC 副本名下，IC 必须在 run dir 内副本上操作。
5. **指标**（新脚本 `check_backward.py`，持久化到 localdata）：(a) escape 段 steps/P_r 双向（正向 [3.7,5] P_out ↔ 反向 [0, 5−3.7] P_out，t_b = t_end−t_f 映射；两者均 ≈ N_s）；(b) de 包络镜像（de_b(t_end−t) vs de_f(t)）；(c) 可逆性：反向末行 vs 正向首行（v 取反）逐粒子 |Δr|/r_scale、|Δv|/v_scale 及 |ΔE|/|E0|。
6. **对 baseline（.base）与 Scheme G（.g）各跑一对**，报告 G 移除了多少不对称。
- 已知边界：反向积分穿过 resonant encounter 段（0.5–1.24 P_out，3 次互穿）→ 末态偏差含 Lyapunov 放大，这是被量化的对象本身，不期待 bit 级还原；ASCII 17 位转写噪声 ~1e-16 rel，可忽略。

#### 5.7.3 回归矩阵

| # | 测试（输入, 模式） | 分辨率 | 二进制 | 期望结果 | 性质 |
|---|---|---|---|---|---|
| R1 | ustabtri `btlogh`（Fixed, `-i 0` 树冻结=初始全椭圆） | S32/64/128 | .g | **byte-identical**（冻结树无双曲节点，clamp 不触发路径的不变性证明） | 硬门 |
| R2 | ustabtri `btlogh_auto`（Adjust, `-m orbit`）— 行为目标 | S32/64/128 | .g vs .base | escape steps/P_r 448→**~128**（S128；X≈3.37e-4）；总步数 ~28.7k→~25.5k（escape 4,522→~1,300）；de_max/endpoint 不劣于 baseline 量级（S128 baseline de_max≈2.1e-9，以 Phase 0 实测为准，容 ≤×3）；ejection 事件 dE 量级不变 | 行为 |
| R3 | ustabtri `logh`（`--g-func 0`） | S32/64/128 | .g 构建 | **byte-identical**（宏门证明） | 硬门 |
| R4 | ustabquin `btlogh`（Adjust; 瞬态双曲拟合吃到 band）— **KEY RISK** | S128/S256/S512 | .g | de_max 不劣化（baseline 8.60e-7 / 1.05e-9 / 1.62e-11 @ 360,713/731,173/1,452,446 步，容 ≤×2）；N_intact 校准不变（128.9/256.4/512.8）；步数允许下降 | 硬门 |
| R5 | ustabquin `logh` | S128/256/512 | .g 构建 | byte-identical | 硬门 |
| R6 | ustabquin short-run 对（`ustabquin_short_new.sh` 设定） | S128 | .g | bin1 近星点浓度保持（de_peak ~5e-12 级 @ 匹配分辨率） | 辅助 |
| R7 | ustabtri H4（hermite 以新头文件重建） | S128 | hermite.g | 汇总量不变（AR group 均束缚，clamp 预计不触发；历史先例：数值不变） | 硬门 |
| R8 | backward 对称性（5.7.2 协议） | S128 | .base 与 .g | 诊断性：双向 escape rate ≈N_s；G 的不对称 ≤ baseline | 诊断 |

**范围排除（2026-09-12 用户裁定）**：btlogh_tri（稳定三体）、quad（B-B 四体）等全束缚测试与本问题（双曲逃逸段过分辨率）无关，不纳入本矩阵——clamp 不触发路径的 bit-identity 由 R1（冻结树）与 R3（LogH 宏门）覆盖。

#### 5.7.4 分阶段执行

**Phase 0 — 基线钉死（不改代码）**

- **目标**：证明当前 worktree（experiment @598abd9 + 未提交 information.h 清理 = v2.1）可确定性复现存档日志，锁定全部 baseline 数字。
- **行动/命令**：
  - `cd /home/lwang/code/SDAR && make -C sample/AR && make -C sample/Hermite`
  - `install -m 755 sample/AR/build/ar.btlogh.ttl.sd.cm ~/bin/ar.btlogh.ttl.sd.cm.base`（logh/blogh 同理 `.base`；hermite → `hermite.base`）；记录 md5 与日期于 commands.log。
  - run dir：`mkdir -p /home/lwang/localdata/SDAR_BLogH/schemeG && cd $_ && cp ../ustabtri ../ustabtri.orbit ../ustabquin ../ustabquin.orbit .`
  - 复跑（示例 S128）：`~/bin/ar.btlogh.ttl.sd.cm.base --g-func 1 -m orbit -t 0.022214414690791832 --ds-scale 0.25 --slowdown-ref 1e-20 --slowdown-timescale-max 1e10 -o 6.103515625e-05 ustabtri > base.ustabtri.btlogh_auto.s128.log 2> base.ustabtri.btlogh_auto.s128.err`（btlogh/logh 去 `-m orbit`；ustabquin: `-t 0.9125502020940626 -o 1.52587890625e-05 -m orbit --ds-scale {0.265,0.1332,0.0666}`）。全部命令记入 `schemeG/commands.log`。
- **通过标准**：R1、R3、R5 对应 base 复跑与存档日志 `diff` 为空（确定性先例成立）；R2/R4 base 数字（步数、de_max、escape 率）与 §5.2/NOTES 一致。
- **风险**：若不可复现 → 未提交清理非纯 cosmetic 或存档来自异构二进制 → 中止并查 diff（Phase 0 是唯一拦截点）。

**Phase 1 — 实现 Scheme G + 无 clamp 路径 bit-identity**

- **目标**：落地 5.7.1 改动；证明 clamp 不触发路径（冻结树 / LogH 宏门）逐位不变。
- **行动**：编辑 `src/AR/symplectic_integrator.h`（5.7.1 的 1.1/1.2）；`make -C sample/AR`；`install -m 755 sample/AR/build/ar.btlogh.ttl.sd.cm ~/bin/ar.btlogh.ttl.sd.cm.g`；确认 `ar.logh/ar.blogh` 产物与 Phase 0 构建 md5 一致（宏无泄漏）；`make -C sample/Hermite` 编译通过。
- **检查**：R1（ustabtri `btlogh`，`-i 0` 冻结树）、R3（ustabtri `logh`）用 `.g` 构建复跑 → 与 base 输出 `diff` 为空。
- **通过标准**：R1/R3 byte-identical；LogH 二进制 md5 不变；H4 编译通过。
- **风险**：低——clamp 只在 `semi<0 && ecc>1 && r_sep>q` 分支内改值。

**Phase 2 — ustabtri 行为验证**

- **目标**：R2 达标（448→~128）且精度不退化；R1/R3 bit-identity。
- **行动**：`schemeG/` 内跑 `.g` × {btlogh, btlogh_auto, logh} × {s32,s64,s128}（s32/s64 至少 btlogh_auto；命令同 Phase 0 换 `.g` 与 `g.` 前缀）。拷贝并改 `path`/`fname` 两行得到 `check_ustabtri_escape_g.py` / `check_ustabtri_escape2_g.py`（原文不动；脚本预测列需增 clamp-aware 公式：steps/P_r = N_s·(P_r/T_h)·(q/max(q, min(r_sep, X)))）。Python：`/home/lwang/.pyenv/versions/general3.12/bin/python`（sdar 经 `make -C tools` 装于 ~/include/sdar）。
- **通过标准**：escape 段 orbit-averaged steps/P_r = 128±15%（S128）；总步数 ≈25.5k（允 ±10%）；de_max 与 endpoint de ≤ baseline×3；ejection 事件（t≈3.7 P_out）dE 跃变量级不变；R1/R3 diff 为空。
- **风险**：若 de 退化 → 检查 X 是否被瞬态拟合放大（band 过宽）→ 转入 D3 备选（X 上限系数 C_x）。

**Phase 3 — ustabquin burst 回归（KEY RISK 硬门）**

- **目标**：瞬态双曲拟合段吃到 band 后 de_max 不劣化。
- **行动**：`.g` × ustabquin {S128,S256,S512}（btlogh+logh，参数见 Phase 0）+ short-run 对（R6）；与 base/存档比 de_max、nsteps、N_intact。
- **通过标准**：de_max ≤ baseline×2（三分辨率）；N_intact 偏差 <1%；R5 diff 为空；R6 de_peak 同量级。
- **风险**：多点同 clamp 的乘积相互作用未逐节点证明——若恶化：备选 (i) `X = q·min(P_r_eff/T_h, C_x)`，C_x≈3–10；(ii) 仅当该节点为树中唯一 hyperbolic 非叶节点时启用 X 扩展。任一备选需重走 Phase 1–3。

**Phase 4 — 反向积分对称性协议（新）**

- **目标**：量化 baseline 的方向不对称并证明 G 不增大（理想：减小）。
- **行动**：按 5.7.2 机制，对 `.base` 与 `.g` 各做 正向(`--print-precision 17`) → 末态提取 → 反向 IC → 反向运行；`check_backward.py` 输出 (a)(b)(c) 三指标对比表。
- **通过标准**（诊断性）：双向 escape rate 均 ≈N_s；G 的 |ΔE|/|E0| 与 Δr 指标 ≤ baseline 同类量。
- **风险**：混沌放大使指标噪声大 → 以 base-vs-G 差分呈现，不设绝对阈值。

**Phase 5 — H4 + 收尾 + 交付决策**

- **目标**：头文件传播面闭合；归档；安装/提交决策。
- **行动**：`make -C sample/Hermite` → `hermite.g`；跑 `ustabtri.h4.s128`（`--r-group 0.002 --r-neighbor-over-group 20 --ds-scale 0.25 --dt-max-power 14 --slowdown-ref 1e-20 -e 1e-4 -o 14 -t 0.022214414690791832`）→ R7；全矩阵终跑归档（schemeG/ 目录 + commands.log）；用户决策：(a) 覆盖 `~/bin/ar.btlogh.ttl.sd.cm`（先备份旧二进制，先例 backup_* 目录），(b) commit 到 experiment 分支，(c) release.note/VERSION bump。
- **通过标准**：R7 汇总不变；矩阵全绿；文档清单（5.7.6）完成。
- **风险**：~/bin 标准名覆盖使历史日志不可复现 → 未获明确确认前保持 scratch 名。

#### 5.7.5 未决设计点与建议

- **D1 P_r_eff 用 live κ 还是 epoch 冻结 κ**：建议 live（`getEffectivePeriod()` 与 `calcEffectivePeriod` 同源；κ epoch 内漂移已有 ds 侧先例）；冻结版需额外存储且与 ds 语义分裂。
- **D2 P_r_eff 候选是否含 elliptic internal nodes**：建议含（flat 全扫；与 `calcBLogHDsIter` 的 P_eff 候选集合同源，且更保守）。
- **D3 X 是否加绝对上限**：建议不加（代数自洽：e→1⁺ 深拟合时 q 与 T_h 同步缩小 → X<q 自动回落 q-cap；e≫1 时 X 有界），以 Phase 3 经验门兜底；若 quin 恶化再引入 C_x。
- **D4 backward 是否多分辨率**：建议仅 S128（诊断协议，成本优先）。
- **D5 安装/提交/版本**：全矩阵通过 + 用户确认后一次性执行；scratch 名 `.base/.g` 期间 ~/bin 标准名不动。
- **D6 是否补跑 ustabtri2**：建议不跑（已被 ustabtri 取代），除非用户要求。

#### 5.7.6 文档更新清单（收尾）

- `NOTES`（`/home/lwang/src/python/SDAR_blogh_method_NOTES.md`）：新 session 节（Scheme G 机制/代数/回归数字/backward 协议结果）；改写 "TODO-session setup: can BTLogH's ds relax after an ejection?" 节结论；§File naming 增补 schemeG/ 与 `.base/.g` 约定。
- 论文 `paper.tex`：§4（ds estimation 的节点 cap 描述 → 双侧 clamp；eq:nk_resolution 逃逸尾讨论）；§5.4.1 ustabtri 表 Adjust 总步数（28,728→~25.5k）与 "measured settled rate 445 = N_s×3.5" 段（→~128）；§5.4.2 quin 数字（若变）；abstract/summary 中 "1.18× LogH" 类步数比句。
- `docs/hierarchical_blogh_impl_notes.md`：`processOuterNode` 实现描述同步（两处：q-cap 段与 2026-08-26 注释所指）。
- localdata：`check_ustabtri_escape{,2}.py` 头部机制注释 + clamp-aware 预测列（改副本或升级原文由用户定）；新增 `check_backward.py` 持久化。
- `release.note` + `VERSION`（若提交，v2.1→v2.2 由用户定）；**跨仓库**：SDAR 头文件变更会传播至 PeTar AR interface（PeTar 侧重建+回归为后续独立任务，非本计划范围）。

#### 5.7.7 总体风险与缓解

1. quin burst 恶化（最大风险）→ Phase 3 硬门 + D3 备选带系数版。
2. baseline 不可复现 → Phase 0 拦截，禁止带病前进。
3. κ 漂移使 X 时变（epoch 内折点缓移）→ 先例成立；监控 de 包络与 backward 指标。
4. 多 hyperbolic 节点乘积相互作用 → 仅经验门（quin）覆盖；如需理论补证留给 §5.6。
5. 覆盖标准二进制/提交过早 → 全程 scratch 名，Phase 5 用户确认制。
