# Hierarchical BLogH 实现笔记

> **最后更新**: 2026-08-27
> - 2026-08-27：**g-func 宏体系重构**（详见下文"g-func 重构记录"节与"编译 Flag 速查"）——`AR_G_FUNC_MUL_POT` 拆分为每方法一宏（`AR_G_FUNC_BLOGH/NORM_BLOGH/MUL_ALL_POT/BTLOGH`），`ADD_POT` 改名 `ADD_INNER_POT`；`--g-func` 统一模板 0=LogH / 1=本方法 / 2=auto（BTLOGH 与 MUL_ALL_POT 无 2），`--g-func-switch` 删除；输出 `g_func` 列改打生效状态 0/1；宏体系集中定义于 `src/AR/g_func.h`。回归验收：全方法×{fixed,auto,gf0} 物理列逐位一致。
> - 2026-08-26b：**束缚段跳变归因与对照判决**（详见 `hierarchical_blogh_plan.md` §3.5）——分辨率减半（S512）、关闭 tree 重构（fixed）、切 LogH 三思路全部证伪：共振飞掠硬化事件（t=0.0368，内双星 a 9× 收缩）的 ~1e-7 跳变是事件型 floor，不随 ds 收敛、不因重构而变、LogH 反差 6 个量级；tree 重构+ds 重算实为保护机制（事件后窗口差 90×）。**使用警示**：standalone base 模式（无 `-m orbit`）不更新根数，交换后 stale semi>0 使 q-cap 永不触发 → dt 爆炸（fixed s512 终点 dE=218）；发生交换/逃逸的场景必须用 orbit/full 模式（生产路径 integrateToTime 不受影响）。
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

1. **BTLogH (新 `--g-func 1`，旧 g_func=4)**: 树级层级 BLogH，支持 B-B 四星及以上。✅ 已实现并验证。
2. **Auto-Switching**: 新 `--g-func 2`（旧 `--g-func 4 --g-func-switch auto`）。⚠️ 现状：BTLOGH 构建不提供选项 2（CLI 拒绝）——双曲层级由 q-cap 处理，无需 LogH 回退（logh 对照同段差 ~6 个量级，见数据档案）。

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

每二进制编译**且仅编译一个**方法宏（互斥，见 `src/AR/g_func.h`）：

| Flag | 含义 | `--g-func` 合法值 |
|------|------|:---:|
| `AR_G_FUNC_BLOGH` | 内层对势乘积（BLogH） | 0,1,2 |
| `AR_G_FUNC_NORM_BLOGH` | (内层乘积)^(1/N) | 0,1,2 |
| `AR_G_FUNC_MUL_ALL_POT` | 所有对势乘积（需 `--s`） | 0,1 |
| `AR_G_FUNC_BTLOGH` | 树级乘积（内层×外层节点） | 0,1 |
| `AR_G_FUNC_MAX_POT` | 最大内层 pair potential | 0,1,2 |
| `AR_G_FUNC_ADD_INNER_POT` | 内层 pair 之和（经 `calc_gt_cross` 间接实现） | 0,1,2 |

任一方法宏 → 自动定义 umbrella `AR_G_FUNC`；前四个方法（乘积家族）另自动定义
`AR_G_FUNC_MUL_POT_FAMILY`（共享 GtKickInv 乘积结构；MAX/ADD 结构不同，故家族宏
≠ 伞形宏，不可合并）。旧宏 `MUL_POT`/`ADD_POT` 已删除（保留 #error 提示）。

二进制目标（`sample/AR/Makefile`，命名 `ar.<method>[.ttl][.sd][.cm][.mpfrc]`）：
`ar.{blogh,normblogh,mulall,btlogh}.ttl.sd[.cm][.mpfrc]`、`ar.maxpot.ttl.sd.cm`、
`ar.addpot.ttl.sd.cm`。

## g-func 重构记录（2026-08-27）

运行期状态由两个成员承载：`int g_func`（CLI 意图 0/1/2，积分中不变）与
`bool g_func_on`（当前生效态，auto 每步可能翻转）。为何需要两个：单一 int 无法
同时携带"意图=auto"与"当前=0/1"两个正交信息，且生效态出现在 pair 热路径
（`calcAccPotAndGTKickInvTwo`）与切换驱动的能量修正/ds 重算中，必须缓存。
已删除：`g_func_user`、`g_func_switch`、`GFUNC_FIXED/GFUNC_AUTO` 枚举。

**`calc_gt_cross` 统一式**（`ADD_INNER_POT` 的实现核心——基础 sum 累积 + 跳过
cross 对 = 只对最内层对求和）：

```cpp
bool calc_gt_cross = true;              // LogH(0) 恒 true；MUL_ALL_POT 恒 true
#ifndef AR_G_FUNC_MUL_ALL_POT
if (g_func_on) calc_gt_cross = false;   // 其余五法生效时：仅 inner 对计入
#endif
```

旧条件 `(g_func>0&&g_func<=2)||g_func==4` 与此逐宏等价（MAX/ADD 的 g_func=1
命中 `<=2`）。

**分发点 → 宏映射**（维护参考）：

| 位置 | BLOGH | NORM_BLOGH | MUL_ALL_POT | BTLOGH | MAX_POT | ADD_INNER_POT |
|---|---|---|---|---|---|---|
| `GtKickInv` struct | family | family | family | family | 专属 | 基础 |
| pair 累积（`calcAccPotAndGTKickInvTwo`） | 乘积 | 乘积 | 乘积 | 乘积 | max+平滑 | sum |
| `pow(1/nbin)` 后处理 | — | ✅ | — | — | — | — |
| `processOuterNode` | — | — | — | ✅ | — | — |
| `dgt_drift_inv` 缩放 | `*=value` | `*=value/nbin` | `*=value` | `*=value` | —（专属 gtgrad 路径） | — |
| κ-capping（Step5，P_eff_min 计算） | — | — | — | ✅ | — | — |
| `checkGFuncCriterionIter` | ✅保守式 | ✅保守式 | ✅保守式 | 不编译 | ✅保守式 | ✅保守式 |
| auto-switch（Step6/initial） | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ |
| `calcDsAndStepOption` auto-ds | 乘积式 | 几何平均 | abort（需--s） | 乘积式+节点势 | min 式（忽略参数） | min 式（忽略参数） |

**新旧命令映射与验收**（基线为 git HEAD 重建的旧宏二进制；对比排除 wall-clock
的 profile 计时列 Total(s)/Int(s)，g_func 列按非零→1 映射）：

| 旧命令 | 新命令 | 验收 |
|---|---|---|
| mulpot `--g-func {1,2,4}` | blogh/normblogh/btlogh `--g-func 1` | 物理列逐位一致 |
| mulpot `--g-func k --g-func-switch auto` | 对应宏 `--g-func 2` | 逐位一致，翻转序列一致 |
| maxpot/addpot `--g-func 1`、auto | 同左 / `--g-func 2` | 逐位一致 |
| mulpot `--g-func 0` | 任一 g-func 宏 `--g-func 0` | 逐位一致（gf0 回归） |

实测 11 组全部逐位一致。gf3（MUL_ALL_POT）：旧代码 auto-ds 路径直接 abort
（`--s` 覆盖发生在其后，死路）；新代码 gf1 无 `--s` 时 CLI 报错，给了 `--s` 则
以占位 min-ds 初始化（仅用于设置 fix_step_option）后覆盖，gf2 直接拒绝；
`-m full` 自适应实测 dE≈2e-9。

顺带修复两个存量 bug：`-f`（快照输出）以 "r" 模式 fopen 新文件必然 abort
（改为 "w"）；`-l`（load 重启）跳过 `initialIntegration` 导致 `--g-func` 静默
失效（改为在 load 分支从 CLI 重建 `g_func/g_func_on`，auto 的历史翻转态不可
恢复，按 criterion 现场评估一次）。

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
--g-func 1                          # 始终本方法（按二进制宏）
--g-func 0                          # 标准 LogH

# 自动切换（blogh/normblogh/maxpot/addpot 提供；btlogh/mulall 拒绝）
--g-func 2                          # 本方法 ↔ LogH 自动切换
```

## g_func 模式速查（2026-08-27 统一模板）

| --g-func | 含义 |
|--------|----------|
| 0 | $\sum_{\rm all} U_{ij}$（标准 LogH） |
| 1 | 本二进制宏对应的方法（见编译 Flag 速查） |
| 2 | auto：1 ↔ 0 动态切换（btlogh/mulall 无） |

各方法 $g$ 定义：BLogH $\prod_{\rm inner} U_{ij}$；NORM $(\prod_{\rm inner} U_{ij})^{1/N}$；
MUL_ALL $\prod_{\rm all} U_{ij}$；BTLogH $\prod_{\rm n} U_{\rm n}$（树级，B-B 四星+）；
MAX $\max_{\rm inner} U_{ij}$；ADD_INNER $\sum_{\rm inner} U_{ij}$。

输出 `g_func` 列打印**当前生效状态**（0=LogH，1=本方法）；auto 模式下每步可翻转。
旧数据映射：1/2/3/4（含 -1/-2 auto 组合）→ 1。

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

> **N_particle / N_sd 必须与日志匹配**（列数敏感）：quad_sd2 用 `(4, N_sd=?)`、ustabtri 用 `(3, 2)`（见数据档案节）；N_sd 可从表头 SD 块数确认。

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

### ustabtri 数据档案（论文写作用，2026-08-26 定稿）

**初始条件**（输入 `ustabtri` 由 `keplertree ustabtri.orbit` 生成，3 体，G=1）：
- inner binary: m=(0.1, 0.9), a=1e-3, e=0.9
- outer: m=(1+1 等效), a=0.01, e=0.9，与 inner 反向（倾角 π）
- 轨道演化：t≈0.0022 交换 #1 → 0.0364 共振飞掠（内双星 a 从 5.2e-3 硬化至 5.7e-4，9×）→ t≈0.0598 交换 #2 + 双曲逃逸至 r≈1.38（v≈45）→ t_end=0.0889

**运行矩阵**（`ar.ttl.sd.t.mulpot.cm`，输出间隔 `-o 6.103515625e-05`，`-t 0.08885765876316733`，`--slowdown-ref 1e-20 --slowdown-timescale-max 1e10`）：

| 标签 | 命令差异 | 数据文件（`/home/lwang/localdata/SDAR_BLogH/`） |
|---|---|---|
| auto s256 | `--g-func 4 --g-func-switch auto -m orbit --ds-scale 0.125` | `fixg_test/ustabtri.btlogh_auto.qcap.log`（q-cap 后）；`ustabtri.btlogh_auto.s256.log` 为 8/26 用户重跑，位级一致 |
| fixed s256 | `--g-func 4 --ds-scale 0.125`（无 `-m orbit`：无 tree 重构/ds 重算） | `ustabtri.btlogh.s256.log` |
| logh s256 | `--g-func 0 --ds-scale 0.125` | `ustabtri.logh.s256.log` |
| s512 三模式 | 同上，`--ds-scale 0.0625` | `ustabtri.*.s512.log`（`ustabtri.sh` 当前版本即此批） |
| h4 参照 | `hermite --r-group 0.002 --r-neighbor-over-group 20 --dt-max-power 14 -e 1e-4 -o 14 -t ...` | `ustabtri.h4.log` |

**关键数值**（束缚段末 = t≈0.0597 交换 #2 前）：

| run | 束缚段末 \|dE\| | 终点 \|dE\| | Nstep | 备注 |
|---|---|---|---|---|
| auto s256 | 9.9e-7 | **9.9e-7** | 873238 | 推荐配置；逃逸尾零误差（终点=束缚段末） |
| auto s512 | 6.7e-6 | 6.75e-6 | 1788718 | 分辨率翻倍反差 7×：事件型 floor 证据 |
| fixed s256 | 4.7e-6 | 1.4e-4 | 77519 | 尾部亦有失配累积 |
| fixed s512 | 7.9e-7 | **218** | 119545 | stale semi>0 → q-cap 不触发 → dt 爆炸（§已知局限） |
| logh s256 | 2.4e-2 | 1.8e-2 | 248908 | |
| logh s512 | 8.4e-5 | 1.2e-4 | 568460 | 收敛快（~100×/减半）但同成本差 3–4 个量级 |
| h4 | ~0.33 全程 | 0.33 | — | 参照 |

**前 q-cap 基线**（8/15 原始 log 已被覆盖，数字以 plan §3.4 表为准）：auto s256 终点 3.8e-2、ds 冻结 0.0127135、Nstep 774265。

**归因结论**（详见 plan §3.5）：束缚段剩余误差为强相互作用**事件型 floor**（~1e-7–1e-6）——不随 ds 收敛（S512 判决）、不因 tree 重构而变（fixed 对照：事件瞬间增量同为 1e-7 级，但事件后窗口 auto 好 90×）、LogH 同段差 ~6 个量级。分析要点：跳变与共振飞掠硬化后的内双星极近心同时（t=0.03680，+2.2e-7）；通道排除——跳变处 ΔdE_SDC=0、生效 κ 不变、非重建行、q-cap 未触发。

**读回与对比注意事项**：`N_particle=3, N_sd=2`；计时列 18–20（Total/Int/Int_tsyn）每次运行不同，位级对比需排除；tree 拓扑变化用 SD 块的 (I1,I2) 列检测。

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
- **Auto-switching (4↔0)**: 判据已旁路（见 Summary 2），切换从不发生；mode 1/3 的 auto 路径未深度测试
- **standalone base 模式不适用于交换/逃逸场景**: 无 `-m orbit` 时不更新根数，交换后 stale semi>0 使 q-cap 永不触发 → dt 爆炸（fixed s512 终点 dE=218，见数据档案）。生产路径（Hermite/PeTar 经 `integrateToTime`）有根数更新，不受影响
- **强相互作用事件型误差 floor ~1e-7–1e-6**: 共振飞掠硬化事件的跳变不随 ds 收敛、不因重构而变（三组对照判决，见 plan §3.5）；突破需事件时刻相位精确穿越，未规划

---

## g_func 三变量设计

| 变量 | 含义 | 取值 |
|------|------|------|
| `g_func` | 当前生效的 g 函数 | 0-4 |
| `g_func_user` | 用户 CLI 选择的方法 | 0-4 |
| `g_func_switch` | 自动切换模式 | `GFUNC_FIXED=0`, `GFUNC_AUTO=1` |

`checkGFuncCriterionIter()` 遍历 tree 检查扰动比——所有内层 binary `pert_ratio < 1` 时可用 g_func，否则退为 0。提供 auto 的方法（blogh/normblogh/mulall/maxpot/addpot）用此保守公式；BTLogH 不编译此函数（CLI 拒绝选项 2，双曲外层由 `processOuterNode` 的 q-cap 处理）。

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
