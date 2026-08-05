# Hierarchical BLogH 实现笔记

> **最后更新**: 2026-08-05
> - 重构：`hybrid_switch` → `g_func` / `g_func_user` / `g_func_switch` 三变量拆分
> - 统一：MUL_POT / MAX_POT / ADD_POT 三种方法共用 `--g-func` + `--g-func-switch` CLI
> - 重命名：`AR_HYBRID` → `AR_G_FUNC`，`AR_TIME_FUNCTION_*` → `AR_G_FUNC_*`

## Summary

1. **BTLogH (g_func=4)**: 树级层级 BLogH，`--g-func 4`，支持 B-B 四星及以上。✅ 已实现并验证。
2. **Auto-Switching**: `--g-func 4 --g-func-switch auto`，自动在 4↔0 间切换。✅ 已实现，待深度测试。

## Files Changed

- **`src/AR/symplectic_integrator.h`**:
  - `processOuterNode()` — fused traversal: U_node 乘入 gt_kick_inv + 梯度分发，一次递归完成
  - `addOuterGradientToMember()` — 叶子层梯度分发辅助函数
  - `calc_gt_cross` 门控扩展：`(g_func>0 && g_func<=2) || g_func==4`
  - `kickEtotAndGTDrift` 缩放扩展：`g_func==1 || g_func==3 || g_func==4`
  - `g_func` / `g_func_user` / `g_func_switch` — 三变量设计（`AR_G_FUNC`）
  - `switchGFuncAuto()` — 通用 auto-switch（`AR_G_FUNC`）
  - `checkGFuncCriterionIter()` — 扰动比判据（原 `checkHybridMethodCriterionIter`）
  - `GFUNC_CHECK_INTERVAL=100` — 检查频率

- **`src/AR/information.h`**:
  - `calcBLogHDsIter()` — BLogH 乘积 ds 公式
  - `multiplyDsByNodePotentials()` — hierarchical ds 缩放
  - 修复了所有 BLogH 模式 (1/2/3/4) 的 ds 计算

- **`sample/AR/ar.cxx`**: `--g-func 0-4` + `--g-func-switch fixed|auto` CLI（三种 POT 统一）
- **`src/Hermite/hermite_integrator.h`**: `g_func` 传递给 `calcDsAndStepOption`
- **`tools/ar.py`**: `g_func` → `g_func` 列名

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

## Auto-Switching — 2026-08-05 重构

### 概述

通过 `--g-func-switch auto` 启用自动切换。与 `--g-func` 指定的方法共同决定行为：

| `--g-func` | `--g-func-switch` | 行为 |
|------------|-------------------|------|
| 4 | fixed | 始终 BTLogH |
| 4 | auto | BTLogH ↔ LogH 自动切换 |

通用设计：`--g-func N --g-func-switch auto` 对任意 N=1-4 均自动在 N ↔ 0 间切换。

### 三变量设计

| 变量 | 含义 | 取值 |
|------|------|------|
| `g_func` | 当前生效的 g 函数 | 0-4（永远是合法分支选择器） |
| `g_func_user` | 用户通过 CLI 选择的方法 | 0-4 |
| `g_func_switch` | 自动切换模式 | `GFUNC_FIXED=0`, `GFUNC_AUTO=1` |

### 实现位置

- **`src/AR/symplectic_integrator.h`**:
  - `g_func` / `g_func_user` / `g_func_switch` — 三个独立变量（`AR_G_FUNC`）
  - `g_func_switch_pending_` / `g_func_switch_confirm_` — hysterisis 状态（`AR_G_FUNC_MUL_POT`）
  - `switchGFuncAuto()` — 通用切换函数，target = criterion ? `g_func_user` : 0
  - `checkGFuncCriterionIter()` — 扰动比判据（原 `checkHybridMethodCriterionIter`）
  - `initialIntegration()` — 三变量初始化
  - `integrateToTime()` — 每 100 步检查

- **`sample/AR/ar.cxx`**: `--g-func` + `--g-func-switch` CLI
- **`src/Hermite/hermite_integrator.h`**: `hybrid_switch` → `g_func`

### CLI 用法

```bash
--g-func 4 --g-func-switch fixed     # 始终 BTLogH
--g-func 4 --g-func-switch auto      # BTLogH ↔ LogH 自动切换
--g-func 1 --g-func-switch auto      # BLogH ↔ LogH 自动切换（通用）
```

### 判据

`checkGFuncCriterionIter()` 的扰动比判据：

$$\text{pert\_ratio} = \frac{M_{\text{out},1}M_{\text{out},2}}{M_{\text{in},1}M_{\text{in},2}} \cdot \left(\frac{a_{\text{in}}(1+e_{\text{in}})}{a_{\text{out}}(1-e_{\text{out}})}\right)^3$$

- 所有 inner binary `pert_ratio < 1` → `g_func = g_func_user`
- 任意 inner binary `pert_ratio >= 1` 或 hyperbolic → `g_func = 0`

### 切换时的关键操作

| 操作 | 说明 |
|------|------|
| `calcAccPotAndGTKickInv()` | 用新 `g_func` 重算力和 g |
| `gt_drift_inv_` 调整 | 变化 > 0.1% 时直接重置，否则 `+= diff` |
| `info.calcDsAndStepOption()` | 用新 `g_func` 重算 ds |
| `ds[0]`/`ds[1]`/`ds_init`/`ds_backup` 同步 | 取 `min(当前, 新 ds)` |

### Hysteresis

要求连续 3 次判据一致才切换。`g_func_switch_pending_` 记候选目标，`g_func_switch_confirm_` 记连续确认次数。

### 架构决策

- **内置于 `integrateToTime()`**（方案 A）：standalone AR、Hermite、PeTar 自动继承
- **不每步重建 binary tree**：判据依赖的 orbital 参数通过 `updateBinarySemiEccPeriodIter` 更新
- **检查频率**：每 100 步

### 编译状态

✅ 所有 AR binary variants 编译通过，零 warning。
✅ Hermite sample 编译通过。

### 测试状态

| 测试 | 状态 |
|------|:---:|
| 编译通过 | ✅ |
| 冒烟测试（固定 + auto 模式运行不崩溃） | ✅ |
| 层级三体 KL 循环（g_func 切换验证） | ❌ 待设计 |

### 设计文档

详见 [`hierarchical_blogh_plan.md` §2](./hierarchical_blogh_plan.md)。
