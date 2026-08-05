# Hierarchical BLogH: Theory Reference

> **用途**: 本文档提供 BTLogH (hybrid_switch=4) 的理论基础与后续开发计划（自动切换等）。
> 实际代码实现见 [`hierarchical_blogh_impl_notes.md`](./hierarchical_blogh_impl_notes.md)。
> 原始计划中关于梯度无需修正、`multiplyOuterNodePotentials()` 用 `_bin.semi` 等方案在实测中被发现有问题，
> 已被 `processOuterNode()` 替代。本文档仅保留经验证正确的理论部分。

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
11. [ ] 深度测试：层级三体 KL 循环（g_func 切换验证）

---

## 3. 其他后续方向

1. **更深层级测试**: 当前仅验证了 B-B 四星（2 层）。对 3+1 四星、5 体等更深
   binary tree，需验证 `processOuterNode()` 的递归梯度分发正确性。

2. **与 PeTar 集成**: BTLogH 当前仅在 standalone SDAR 中可用。PeTar 中使用
   SDAR 处理 close encounters 时，需传递 `hybrid_switch=4` 并确保
   binary tree 结构在 PeTar 的 group 检测后仍然有效。

3. **Slowdown + 自动 ds 调优**: 当前 ds 用 Eq.~(ds\_combined) 一次性估计。
   对于 $\kappa$ 随时间变化的情况（如 KL 循环中 inner binary 周期变化），
   可能需要动态调整 ds。

