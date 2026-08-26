# Hierarchical BLogH 实现笔记

> **最后更新**: 2026-08-26
> - 2026-08-26：**双曲层级 ds/g 规范失配修复（g 函数侧 q-cap）**——详见下方专节与 `hierarchical_blogh_plan.md` §3。要点：`processOuterNode` 对双曲非叶节点把 g 的距离因子 cap 在 q=|a|(e−1)（与 ds 的 15b 口径同规范），cap 生效时同时抑制该节点梯度（TTL 一致性）；ds 逻辑不动（保持 tree 不变即冻结的设计原则）。ustabtri 终点 |dE| 3.8e-2 → 9.9e-7，交换前位级不变，logh/btlogh 固定模式位级不变。
> - 2026-08-15c：BTLogH ds 双曲修复——双曲叶 ds 恢复 2π/256（原丢失 8× 分辨率）；`P_eff_min` 锚点覆盖全部层级（非叶节点贡献 `P·κ`/遭遇时标），遭遇瞬态由真实最快层驱动
> - 重构：`updateSlowDownAndCorrectEnergy` + `switchGFuncAuto` → 统一为 `syncTreeSlowDownAndDs`
> - 新增辅助函数：`calcBinaryTreeSlowDown`, `applyStableCheckAndSlowDown`, `correctSlowDownEnergy`
> - 中断处理简化：`_is_interrupt` 标志消除重复 tree/force/energy 代码
> - `stableCheckIter` 条件化：仅 tree 重建或轨道更新后触发
> - pert 判据修正：`pert_out/pert_in > 1`（无量纲比值）
> - g_func-aware ds 统一：`update_flag` 路径改用 `calcDsAndStepOption`
> - `synch_flag` → `integration_mode` 重命名，短选项 `-S` → `-m`（ar.cxx）
> - 2026-08-11：tree 配对度量 `r²` → `r³/(m_i·m_j)`（`calcMinDisList` → `calcBindingList`）
> - 2026-08-14：`calcBinaryTreeSlowDown` 的 tree-stale 判据改用瞬时度量（`calcPertFromMR`），与 `generateBinaryTree` 自洽；排除组外摄动
> - 2026-08-15：pert 度量函数（`calcPertFromMR`/`calcPertFromBinary`/`calcPertFromForcePot`）从 interaction 类迁移到 `COMM::Binary`，单一事实源；`calcBindingList` 改调用同一函数（`pairLess` → `pairGreater`）
> - 2026-08-15b：BTLogH ds 引入 perturbation 影响——节点势改轨道平均（semi 口径）+ 节点级扰动缩放；`calcPertRatio` 统一 pert 比值（双曲用瞬时 MR 度量）
> - 2026-08-15c：BTLogH ds 双曲修复——双曲叶 ds 恢复 2π/256（原丢失 8× 分辨率）；`P_eff_min` 锚点覆盖全部层级（非叶节点贡献 `P·κ`/遭遇时标），遭遇瞬态由真实最快层驱动

## Summary

1. **BTLogH (g_func=4)**: 树级层级 BLogH，`--g-func 4`，支持 B-B 四星及以上。✅ 已实现并验证。
2. **Auto-Switching**: `--g-func 4 --g-func-switch auto`，自动在 4↔0 间切换。✅ 已实现，待深度测试。

## Files Changed

- **`src/AR/symplectic_integrator.h`**:
  - `syncTreeSlowDownAndDs()` — 统一函数：slowdown → tree 重建 → stab → g_func → 力 → energy → ds
  - `calcBinaryTreeSlowDown()` — root + inner pert/slowdown 合并
  - `applyStableCheckAndSlowDown()` — 仅设 slowdown factor（不调 stableCheckIter）
  - `correctSlowDownEnergy()` — 能量修正（全量重算 vs sd_backup 缩放）
  - `_is_interrupt` 参数 — 中断路径跳过 pert 检查，tree 重建内化
  - `stableCheckIter` 条件化 — 仅 tree 重建或 `update_flag` 时触发
  - pert 判据：`pert_out/pert_in > 1.0`（无量纲）
  - g_func 每调用点评估（取消 hysterisis / 100 步间隔）
  - `processOuterNode()` — BTLogH 梯度分发
  - `updateBinarySemiEccPeriodIter` 提升为 public

- **`sample/AR/ar.cxx`**: 
  - `--g-func 0-4` + `--g-func-switch fixed|auto` CLI
  - `synch_flag` → `integration_mode` 重命名，短选项 `-S` → `-m`
  - 每步积分前：`updateBinarySemiEccPeriodIter` + `stableCheckIter` + `syncTreeSlowDownAndDs`

- **`src/AR/information.h`**: `calcDsAndStepOption` 统一 ds 接口（g_func-aware）
- **`src/Hermite/hermite_integrator.h`**: `g_func` 传递给 `calcDsAndStepOption`
- **`tools/ar.py`**: `g_func` 列名

## 编译 Flag 速查

| Flag | 含义 | g_func 范围 |
|------|------|:---:|
| `AR_G_FUNC_MUL_POT` | 乘积型 g 函数 | 0-4 |
| `AR_G_FUNC_MAX_POT` | 最大 pair potential | 0-1 |
| `AR_G_FUNC_ADD_POT` | inner pair 之和 | 0-1 |

三者均通过 umbrella flag `AR_G_FUNC`（自动定义）启用 g-func 相关代码路径。

## 关键设计决策

1. **距离用瞬时成员位置**：`pos0 - pos1`，非 `_bin.semi`，正确反映离心率变化
2. **`_bin.r` 过时**：integration 期间不更新（仅在 tree generation 时），故从成员位置实时计算
3. **`gt_kick_inv_.nbin` 不递增**：nbin 仅用于 switch=2 几何平均，switch=4 不读取
4. **ds substep 系数**：`ds_i` 用 coeff=2π（per-orbit），最终 ds 乘 1/32（非 2π/32，否则多一个 2π）
5. **成员位置是当前的**：`processOuterNode` 在 drift 步骤后被调用
6. **USE_CM_FRAME 无关**：`pos0 - pos1` 在 CM frame 和绝对坐标下都成立

## CLI 用法

```bash
# 固定模式
--g-func 4 --g-func-switch fixed     # 始终 BTLogH
--g-func 0                           # 标准 LogH

# 自动切换
--g-func 4 --g-func-switch auto      # BTLogH ↔ LogH 自动切换
--g-func 1 --g-func-switch auto      # BLogH ↔ LogH 自动切换（通用）
```

## g_func 模式速查

| g_func | $g$ 定义 | 适用场景 |
|--------|----------|----------|
| 0 | $\sum U_{ij}$（标准 LogH） | 无层级结构 |
| 1 | $\prod_{\rm inner} U_{ij}$（BLogH） | S-B 三星 |
| 2 | $(\prod_{\rm inner} U_{ij})^{1/K}$（几何平均） | 类似 1，g~energy 量纲 |
| 3 | $\prod_{\rm all} U_{ij}$ | 小 N 系统 |
| 4 | $\prod_{\rm n} U_{\rm n}$（BTLogH） | B-B 四星及以上 |
| -1 | `--g-func 1 --g-func-switch auto` | 动态 0↔1 | 不稳定三星 |
| **-2** | **`--g-func 4 --g-func-switch auto`** | **动态 4↔0（auto (generic)）** | **层次破坏场景** |

## 梯度修正（processOuterNode）

`processOuterNode()` 的核心逻辑：

1. **递归遍历 binary tree**：每个内部节点 n 计算 $U_{\rm n} = G M_{\rm L} M_{\rm R} / r_{\rm s}$（用瞬时成员位置 `pos0-pos1`）
2. **乘入 gt_kick_inv_**: `gt_kick_inv_.value *= U_n`
3. **梯度分发**: 对节点的左右子树成员，按质量分数 $m_i/M_{\rm m}$ 加上 $\pm \hat{r}_{\rm s}/r_{\rm s}$ 到各自的 `force_[i].gtgrad`，符号左 + 右 −

这一步修正了以下问题：plan 阶段（`hierarchical_blogh_plan.md` §0.6）认为"梯度计算与 switch=1 完全一致"，但实际发现缺少外层节点的 $\nabla\ln U_{\rm n}$ 贡献会导致 DKD 格式不一致。`processOuterNode` 在一次 tree walk 中同时完成 $U_{\rm n}$ 乘入和梯度分发。

## ds 公式（最终版）

plan 中的 `ds = ds1*ds2/sqrt(P1*P2)` 已被替换为基于 $P_{\rm eff,min}$ 的公式（见 `paper.tex` Eq. ds_combined）：

$$ds = \frac{\prod_k ds_k \cdot P_{\rm eff,min}}{\prod_k P_{{\rm eff},k}} \cdot \frac{1}{N_{\rm s}} \cdot \prod_{\rm o} U_{\rm n}$$

其中 $N_{\rm s}=32$，$P_{{\rm eff},k} = \kappa_k P_k$。Python notebook 中对应的手动估计：

```python
ds_mix = ds1 * ds2 / np.sqrt(period_in1 * period_in2)  # 近似
ds_auto = ds_mix * P_eff_min / sqrt(P1*P2) * U_o_init / 32  # Eq.(ds_combined)
```

## Python 数据读回

g-func 输出（`g_func!=0`）比标准 LogH 多一个 `g_func` 列。**必须**用 `g_func=True` 构造 reader。

> **旧代码迁移 (pre-2026-08)**：`hybrid=True` → `g_func=True`，`hybrid_flag` 列 → `g_func` 列。

```python
# g-func 输出（g_func=1~4）
snap = sdar.SDARData(g_func=True, N_particle=4, slowdown=True, N_sd=3, time_measure=True)
snap.loadtxt(path, skiprows=1)

# Hermite 输出（时间轴需加 time_offset）
snap = sdar.HermiteData(N_particle=4, N_sd=2, time_measure=True)
snap.loadtxt(path, skiprows=1)
time = snap.time + snap.time_offset
```

不用 `g_func=True` 会导致列错位，所有后续分析数据错误。

## 已验证的测试场景

### 2026-08-26 双曲层级 ds/g 规范失配（g 函数侧 q-cap 修复）

**问题**（unstable triple，`--g-func 4 --g-func-switch auto -m orbit`）：交换 #2 后新双星 + 逃逸星构成双曲根节点，ds 按近心口径 q 冻结（Kepler 常数，逐位恒定），而运行时 g(t) 的外层因子 Gm₁m₂/r 随逃逸按 1/r 衰减 → dt=ds/g 无界增长 → 内双星分辨率跌至 ~8 步/轨道（设计 32），dE 从 1e-6 爬到 3.8e-2。基线归因（Phase 0 脚本 `sample/test/analysis_ustabtri_btlogh.py` 固化）：主导通道是 DKD 分裂误差（尾期生效 κ≡1，κ_max 棘轮 1→53 仅为同源症状）。

**被否决的两版方案**（记录以免反复）：
1. **ds 侧 r_ref = max(q, r_inst)**（计划 v2 Fix A）+ 每步重算 ds（Fix B1）：dE 达 9.9e-7，但违背 BTLogH 设计原则——tree 不变时 ds 应冻结（扩展相空间结构、time symmetry），且实测过保守（尾段分辨率超需 ~40×，+12.8% 步数全花在过度分辨）。
2. **g 侧 cap = 2|a|**：密近三相舞期间瞬态配对的 |a| 可任意小 → cap 使 U 虚胀上百倍，每次树重建产生巨大 g 跳变 → 一次性能量破坏（dE 跳至 1.67、H_sd 恶化到 6e-3）。教训：g 作为状态函数必须对树重建连续有界。

**采纳方案**：`processOuterNode`（`symplectic_integrator.h`）双曲非叶节点（semi<0 且 ecc>1）取 `r_eff = min(r_sep, q)`，q=|a|(e−1)：
- 与 ds 的 15b 口径**同规范** → 分辨率锚定设计值、dt 有界；
- 瞬态 flyby 期间 r≤q 恒成立（q 即近心距）→ 相遇期间 cap 不触发，无接合跳变；
- 接合点 r=q 处值连续（仅梯度折点，可积）；
- **cap 生效时必须同时抑制该节点的梯度分发**（`grad_active=false`），否则 `gt_drift_inv_`（经 `gtgrad` 演化）与 g 值不一致，破坏 TTL 扩展哈密顿量守恒。

**验证**（`/home/lwang/localdata/SDAR_BLogH/fixg_test/`，基线在上级目录）：
| 指标 | 基线 | q-cap |
|---|---|---|
| 终点 \|dE\| | 3.8e-2 | **9.9e-7**（与相遇期末同 level） |
| 交换 #2 前（978 行） | — | **位级一致** |
| logh / btlogh 固定模式 | — | **位级一致**（排除计时列 18–20） |
| Nstep_sum | 774265 | 873238（+12.8%，误差控制器自适应 ds 0.0127→0.0359 部分摊销；生产中逃逸尾由 break 截断，长尾成本不存在） |

- `--break-check`（`ar.cxx` 新增，Fix C）：镜像 `H4::checkBreak` 双曲逃逸分支，默认仅记录（t=0.0598 记录到 r=1.06·r_crit 事件），轨迹零影响。
- 分析脚本 `sample/test/analysis_ustabtri_btlogh.py`：读回需 `SDARData(g_func=True, N_particle=3, slowdown=True, N_sd=2, time_measure=True)`；κ 用 `sd` 列（生效值）、`sd_max` 列=timescale/period 上限。

- **B-B quadruple, quad_sd2**: 两个 inner binary 都触发 slowdown
  - Inner 1: $m=(0.01,0.09)$, $a=10^{-4}$, $e=0.9$, $90^\circ$ 倾角 → KL 离心率振荡
  - Inner 2: $m=(3,7)$, $a=10^{-3}$, $e=0.9$, $90^\circ$ 倾角
  - Outer: $m=(0.1,10)$, $a=0.5$, $e=0.1$
  - $t_{\rm KZ}\approx4.15\times10^3$, $P_{\rm o}\approx0.70$
  - 8 模型对比，两 group: (Group 1) no-sd S32/S64 vs SD6-S64 vs H4-SD6; (Group 2) SD6-S32/S64/S128 vs SD2-S32
- **核心发现**: slowdown 增大每有效轨道分辨率（$n_{\rm step}\kappa/P$），因此 BTLogH+SD 精度反超无 slowdown；强 slowdown (SD2) 在某些诊断量上选择性改善但全局能量误差变差
- 数据路径: `/home/lwang/localdata/SDAR_BLogH/quad_sd2*.log`
- 分析 notebook: `/home/lwang/src/python/SDAR_blogh_method.ipynb`
- 论文: `/home/lwang/write/Astrophysics/BTLogH/paper.tex`

## 已知局限

- **仅测试 B-B 四星**: 更深层级（如 3+1 四星、5 体等）尚未验证
- **ds 公式对 hyperbolic 内双星**: Eq.(ds_combined) 用 $T_{\rm eq}$ 替代 $P_{\rm eff}$，此路径尚未充分测试
- **Auto-switching 未深度测试**: 编译通过、冒烟通过，但尚未在真实层次破坏场景（如 KL 循环）中验证 4→0→4 切换行为

---

## g_func 三变量设计

| 变量 | 含义 | 取值 |
|------|------|------|
| `g_func` | 当前生效的 g 函数 | 0-4 |
| `g_func_user` | 用户 CLI 选择的方法 | 0-4 |
| `g_func_switch` | 自动切换模式 | `GFUNC_FIXED=0`, `GFUNC_AUTO=1` |

`checkGFuncCriterionIter()` 遍历 tree 检查扰动比——所有内层 binary `pert_ratio < 1` 时可用 g_func，否则退为 0。BTLogH (g_func=4) 跳过 pert_ratio 检查（`processOuterNode` 处理外层节点）。

## 架构 — 2026-08-07 重构

### `syncTreeSlowDownAndDs` 统一函数

将原来的 `updateSlowDownAndCorrectEnergy` + `switchGFuncAuto` 合并为单一入口：

```
Step 1: calcBinaryTreeSlowDown      ← root + inner pert/slowdown
Step 2: pert_ratio_max > 1 → tree 重建
Step 3: stableCheckIter（仅重建时）
Step 4: applyStableCheckAndSlowDown  ← 基于已有 stab 设 slowdown factor
Step 5: BTLogH κ-capping
Step 6: g_func 评估 + 切换
Step 7: 力同步 + correctSlowDownEnergy + calcDsAndStepOption
```

### 辅助函数

| 函数 | 职责 |
|------|------|
| `calcBinaryTreeSlowDown(Float* ratio)` | root pert + inner slowdown，可选返回 max pert_out/pert_in |
| `applyStableCheckAndSlowDown(flag)` | 基于已有 stab 设 root slowdown factor |
| `correctSlowDownEnergy(sd_backup, force_recalc, is_interrupt)` | 能量修正：全量重算 或 sd_backup 缩放 |

### 调用点

| 位置 | 调用 | 特点 |
|------|------|------|
| `initialIntegration` | `calcBinaryTreeSlowDown` + `stableCheckIter` + `applyStableCheckAndSlowDown` | 无 tree/g_func/ds |
| `integrateToTime` 循环前 | `syncTreeSlowDownAndDs(true, true)` | 在 ds 初始化之前 |
| `integrateToTime` orbit 更新 | `updateBinarySemiEccPeriodIter` + `stableCheckIter`(if updated) + `syncTreeSlowDownAndDs` | ds 变化时传播到数组 |
| `integrateToTime` 中断 | `syncTreeSlowDownAndDs(true, true, true)` | 内化 tree 重建，跳过 pert 检查 |
| `ar.cxx` 手动循环 | 同 orbit 更新 | 每步 |

### `_is_interrupt` 标志

当 `true` 时：跳过 pert 检查，强制 tree 重建（因中断引起质量变化），能量修正额外更新 `de_sd_change_interrupt_`。

### 中断处理优化

中断处理中的 `generateBinaryTree`、g_func 评估、`calcAccPotAndGTKickInv`、`calcEKin`、slowdown 能量记账全部删除——由 `syncTreeSlowDownAndDs(true, true, true)` 一次完成。中断处理器仅保留非 slowdown 能量记账和 merger 检测。

### 判据与切换

g_func 在每次 `syncTreeSlowDownAndDs` 调用时评估（`GFUNC_AUTO` 模式），不设间隔限制。

### Tree-stale 判据度量统一（2026-08-14）

**问题**：旧判据中 `pert_in` 用 apo 口径（`calcPertFromBinary`：$m_1m_2/\mathrm{apo}^3$），而 `pert_out` 与 tree 配对都用瞬时距离口径，导致内双星近日点附近 ratio 被高估 $(\mathrm{apo}/r)^3$ 倍（$e=0.9$ 时约 7000 倍），频繁假重建。

**修改**（`calcBinaryTreeSlowDown`）：
- `pert_in` 改用 `calcPertFromMR(r_12, m1, m2)`，$r_{12}$ 取两成员**实时**分离（叶子用粒子 `pos`，子树用 `getMemberAsTree()->pos`，两种坐标系下均为活量）；
- `pert_out` 排除组外部分：`bini.slowdown.pert_out - bin_root.slowdown.pert_out`（组外摄动不能改变组内配对，但 PeTar 中会抬升 ratio 造成假重建）；
- 边界保护：`r²>0`、`m1,m2>0`、`pert_out_internal>0`。

**结果**：判据与 `generateBinaryTree` 的 `r³/(m_i·m_j)` 配对度量完全自洽——`pert_ratio_max>1` 即外部潮汐力超过内束缚力，正是配对被翻转的物理条件。apo 口径 `pert_in` 仅保留给 (保守的) slowdown factor 本身。

### 度量函数迁移到 `COMM::Binary`（2026-08-15）

**动机**：`calcBindingList` 原先内联 `r³/(m_i·m_j)`，与 `calcPertFromMR`（interaction 类内）存在 R4 宏口径分裂隐患。且度量是纯牛顿潮汐力，无用户自定义自由度；真正可定制的累积钩子（`calcSlowDownPert*`、PN 修正等）仍在 interaction 类。

**修改**：
- `calcPertFromMR` / `calcPertFromBinary` / `calcPertFromForcePot(G, ...)` 移入 `COMM::Binary`（`binary_tree.h`），成为 tree 配对、tree-stale 判据、slowdown 三者的**单一事实源**（同一函数 → 口径永久一致，R4 分支保留亦可自洽）；
- `calcBindingList` 改调 `COMM::Binary::calcPertFromMR`，方向从 min(r³/m₁m₂) 翻转为 max(m₁m₂/r³)，排序比较器 `pairLess` → `pairGreater`（降序，最紧束缚优先）；
- SDAR sample 的 `interaction.h` / `ar_interaction.h` 保留 deprecated 转发包装（兼容外部用户）；PeTar `ar_interaction.hpp` 直接删除（内部全部调用点已更新）；Hermite 侧 `calcPertFromForcePot` 新增 `G` 参数；
- 兼容性注意：PeTar 粒子 `pos` 为 `F64vec`，访问需通过基类引用而非裸指针。

### BTLogH ds 引入 perturbation 影响（2026-08-15b）

**问题**：tree 重构自洽后，实测发现外天体近心点附近树多次变动时，某次重构使 ds 从 0.03 跳升到 0.36（12×），能量误差从 10⁻¹⁰ 恶化到 10⁻⁷。机理：新配对使叶子 m₁m₂ 放大 10× + `multiplyDsByNodePotentials` 用瞬时 `r_sep` 在近心点高估 U_node（1/(1-e) 倍）；然后 integrator 的 error-based 增长机制看见"误差低于阈值一半"就放行 ds 增长，永久退化。

**修改**（`information.h`，均在 `calcDsAndStepOption` / ds 计算链内）：
- **a) 节点势改轨道平均**：`multiplyDsByNodePotentials` 的 U_node 改用 semi 口径 $Gm_1m_2/a$（双曲/退化轨道用近心点 $|a|(e-1)$ 保守值，再退化为瞬时 `r_sep`）；ds 是"每轨道分辨率"的量，轨道平均势才是正确量纲来源；
- **b) 节点级扰动缩放**：每个非叶节点乘 $\min(1,\ (c\cdot \mathrm{pert_{in}/pert_{out}})^{1/n_{\mathrm{order}}})$（与 LogH 叶子同配方，数据来自 `slowdown.pert_in/pert_out`，即 apo 口径）——层级破碎时 ds 自动退化回内层主导的保守值；
- **d) 双曲叶子缩放修复**：新 `calcPertRatio(_bin)` 统一计算 pert 比值——椭圆用 apo 口径 `slowdown.pert_in`，双曲（semi≤0，apo 口径为负导致 `pow(NaN)`→scale 恒为 1）改用活距离的 `calcPertFromMR`；`calcBLogHDsIter`、`calcDsKeplerBinaryTree`（MUL_POT 与 min 两版）全部统一使用；
- 用户否决了 integrator 侧"重建后 ds 只减不增"的安全网（避免长期性能损失）。

**预期**：层级健康时行为不变（scale≈1、semi≈平均距离）；破碎/重构时 ds 收缩而非膨胀。

### BLogH ds 双曲修复（2026-08-15c）

**gdb 实证**（遭遇瞬间，unstable triple）：叶子 = 双曲飞行对 (0.9,1.0) `semi=-0.00722, ecc=1.164`，`ds_prod=0.189, P_eff_min=period_prod=2.80e-3`，根 `U_node=587`(a_root≈3.2e-4)。对照健康值 ds≈0.03，分解出两个缺陷：

1. **双曲叶丢失 8× 分辨率**：BLogH 原先椭圆/双曲共用 `coeff=2π`，末尾统一 /32 → 双曲只有 2π/32，而 LogH 约定双曲（近心穿越）需 2π/256。修正：双曲叶用 `2π/8`（组合后 = 2π/256）；
2. **P_eff_min 只扫叶子**：遭遇瞬间真实最快时标是**根层侵入轨道**（a_root≈3.2e-4 → P_root≈2.5e-5，比叶虚拟周期快 110 倍）却不在 min 里。修正：非叶节点也贡献其 `P·κ`（椭圆）或遭遇时标（双曲）——与叶子同一惯例，健康深层级下根 κ 大、P_eff 大，min 不变；破碎态 root κ=1 → P_eff=裸周期，正确驱动 ds 收缩。

组合效应（用 gdb 数值推算）：
| 组合 | ds 估计 |
|---|---|
| 现状 | ≈2.37（60×）|
| 修双曲 /8 | ≈0.30 |
| 修 P_min 锚点（×0.0091）| ≈0.021（≈健康值）|
| 两者 | ≈0.0027（遭遇期 10× 保守，合理）|

Tree 重建判据：`pert_out / pert_in > 1.0`（无量纲比值）。`pert_in = m1·m2 / apo³`，`pert_out` 为外部扰动。当任意内层 binary 的比值超过 1 时触发 `generateBinaryTree`。

已知局限：判据可能触发树结构不变的无效重建。在少数粒子系统中开销可忽略，后续可加拓扑变化检测优化。
