# SDAR BTLogH ds 算法修复计划（v2）：双曲层级 ds 规范（gauge）统一

> 状态: 待执行（交给 PeTar Developer）；v2（2026-08-26）全文修订
> 迁移记录: 2026-08-26 由 `PeTar/plans/btlogh-ds-hyperbolic-escape-plan.md` 移入本目录（SDAR 开发计划统一存放处）
> 路径约定: 本计划位于 SDAR 仓库，文中文件路径均相对 multi-root workspace 根（`SDAR/...`、`PeTar/...` 前缀 = 各仓库根；分析数据为绝对路径），与本文档所在位置无关，可直接执行
> 分析数据: `/home/lwang/localdata/SDAR_BLogH/ustabtri.btlogh_auto.s256.log`（关键区段逐行复核）、`ustabtri.sh`
> 关联文档: `SDAR/docs/hierarchical_blogh_impl_notes.md`（2026-08-11/14/15a/b/c 条目）、`SDAR/docs/hierarchical_blogh_plan.md`
> v1 撤销项: ①"跳过双曲非叶节点 U_node"（特例补丁；dt=ds/g 少一个能量因子，逃逸时 dt∝1/U_out(t) 发散更严重）；②"破碎后切回 LogH"（g_func=4 不回退为刻意设计，LogH 并不比 BTLogH 精确）

## 1. 目标与动机

Unstable triple（inner: m=(0.1,0.9), a=1e-3, e=0.9; outer: m=(1,1), a=0.01, e=0.9）在 `--g-func 4 --g-func-switch auto -m orbit` 下发生两次双星交换。**v2 修正 v1 的误差归因**（复核 t=0.0596–0.0601 区段）：

- 交换 #2 的密近相遇本身被正确分辨：相遇期间 ds 收缩到 4.09e-5（Nstep 单输出区间 +2.5e4），dE 保持 ~1e-6 —— 2026-08-15c 的 P_eff_min 锚点有效，**不是**要修的对象；
- 误差在**其后的逃逸尾期**累积：交换后重建把 ds 设为 0.0127135，此后**逐位恒定**至终点；而运行时 g(t) 的外层因子 = G·M_L·M_R/r_esc(t) 按 1/r 单调衰减（r_esc: ~2e-3 → 1.31，无回归）→ dt=ds/g 无界增长 → 新双星分辨率从 ~40 步/轨道跌到 ~2.6 步/轨道（设计值 32），slowdown 因子 κ 同时被棘轮式抬升（5.5 → 53），dE 包络从 1e-6 单调爬到 3.8e-2。

**机制（一般表述）**：束缚层级的 U_node 有"轨道平均"（15b；有效性 = 一个外层周期，g 周期回归自校）；**双曲层级不存在轨道平均**——近心点口径 q=|a|(e−1) 是 Kepler 常数，ds 公式输出与时间无关，而 g 的对应因子是瞬时 1/r。ds 与 g 规范失配 → 相消破坏 → dt 无界。

**修复原则（一般化，非特例）**：ds 中每个 U 因子与运行时 g 的对应因子**同规范**。双曲非叶节点取 r_ref = max(q, r_inst)：接近段（r>q）U∝1/r → 接近遭遇时 ds 收缩（保守，兼修"大步长冲向 plunge"的反向隐患）；近心点退化为现行 15b/c 行为；离开段 ds 与 g 同步衰减 → dt=ds/g = ds_leaf/(32·U_in,inst)，**内双星分辨率严格恒定 32 步/轨道**（1/r 精确相消），本测试预估仅增加 ~6000 步（逃逸段现状 ~800 步，总量 +0.8%）。

修复目标：遭遇期行为不变（15c）；逃逸尾期分辨率恒定、终点 |dE| ≤ 1e-5（基线 3.8e-2）；健康层级（quad_sd2）位级不变。

## 2. 数据分析结论（基线记录，v2 复核修正）

### 2.1 btlogh_auto 运行的 ds / dE 时间线

| 时刻 t | ds | dE（局部） | 事件 |
|---|---|---|---|
| 0 | 0.004759 | 0 | 初始 |
| ≈0.002197 | 1.75e-5 | −2.9e-11 | 交换 #1：tree 重建，ds 收缩（P_eff_min 锚点生效） |
| ≈0.002258 | 3.42e-4 | −1.9e-10 | 重建后恢复 |
| ≈0.002319 | 6.75e-4 | +2.2e-9 | 稳定平台 |
| ≈0.0352–0.0401 | 7.39e-4 | 2.6e-7 | slowdown 上限调整，+9% |
| 0.059753 | **4.09e-5（瞬态）** | −9.9e-7 | **交换 #2 密近相遇中重建**（Nstep 单输出区间 +2.5e4）——相遇被正确分辨 |
| 0.059814 | **0.0127135** | −9.9e-7 | 交换后重建（新配对 (0.9,1.0) + 逃逸星 0.1）；**此后逐位恒定至终点** |
| 0.0600 | 0.0127135 | −9.9e-7 | 误差尚未增长（v1"误差在交换时刻产生"的归因有误） |
| 0.0627 / 0.065 / 0.0706 / 0.0889（终点） | 0.0127135 | 包络 1.5e-5 / 3e-4 / 1–4e-3 / 峰值 3.8e-2 | **逃逸尾期误差单调爬升**；r_esc≈0.14 / 0.26 / 0.55 / 1.31；κ≈5.5 / 10 / 20 / 53；步/轨道 ~40 → ~2.6 |

- `g_func` 全程 = 4（未切回 0，符合设计意图）。
- 终态：逃逸星 m=0.1 距离 1.31、v≈45；末端 Nstep 每输出区间仅 +1–4 → 单 ds 步 ≈ 6.1e-5 ≈ 0.4 P_in（新双星 a≈2.3e-3，P_in≈1.6e-4）。
- ds 恒定的直接原因：双曲根节点的 (a,e) 是 Kepler 常数 → U_node(q) 与伪周期均为常数 → `calcDsAndStepOption` 即使被每步调用也输出同一值——**公式级冻结，与调用频次无关**。
- 误差控制器（`time_error_max=2.5e-14`, `energy_error_relative_max=1e-10`, `ds_scale=0.125`）只门控 ds 增长，不阻碍也不修正物理性失配；实际步长另被输出间隔 6.1e-5 封顶。
- 独立 AR sample 无 group break（`r_break_crit` 仅打印；break 逻辑在 `H4::HermiteIntegrator::checkBreak`）→ 逃逸星被积分到 r=1.31（break 半径 0.001 的 650×）：standalone 测试即最坏情形压力测试；PeTar 生产中会在 r≈crit 处释放（见 Phase 3）。
- 对照：`h4.log` 同样交换+逃逸 → 交换是物理的；`logh` / 固定 `btlogh`（无 `-m orbit`）不重估 ds、轨迹不同，不构成对照组。

### 2.2 机制：ds 与 g 的规范失配（v2 修正版）

ds 公式（`calcDsAndStepOption`，g_func==4）：

```
ds = ds_prod · P_eff_min / period_prod · Π_o U_node · node_scale · (ds_scale/32)
```

运行时时间变换 g(t)=Π U_inst 由 `processOuterNode` 按**瞬时**成员位置求值（08-14 起即如此）。设计意图：U_node(ds) 与 U_node(g) 相消，使 dt=ds/g 给出恒定内层分辨率。三类层级：

| 层级 | U_node(ds) 现行口径 | 有效性 | 问题 |
|---|---|---|---|
| 束缚非叶 | Gm₁m₂/a（轨道平均，15b） | 一个外层周期（g 周期回归，近心自校） | 无（±(1±e) 振荡为已接受行为） |
| 双曲非叶·近心附近 | Gm₁m₂/q（15b 保守口径） | r≈q 时成立 | 无（15c 遭遇期已验证，dE~1e-6） |
| **双曲非叶·离开段** | Gm₁m₂/q = **常数** | **永不**（g 因子 1/r 单调衰减、无回归） | **ds 与 g 失配 → dt∝r_esc 无界** |

被否决的替代口径及失败原因（记录以免反复）：
- **跳过 U_node**（v1 Phase 1）：dt=ds/g 少一个能量因子 → dt∝1/U_out(t)，逃逸时发散更严重；对束缚群内部的双曲 flyby（chain 遭遇）直接错误。
- **仿椭圆用 |a| 口径**：仍是常数（同样冻结）；且对强双曲（e>2 时 q>|a|）在近心点**高估** U → 遭遇期 ds 过大，破坏 15c。
- **正确口径 r_ref = max(q, r_inst)**：与 g 自身、`calcPertRatio` 双曲分支（15c 的 `calcPertFromMR` 瞬时度量）同规范。代入后 dt=ds/g 的 1/r 因子精确相消 → 分辨率不变量 32 步/轨道；接近/近心/离开三阶段连续覆盖，无特例分支。

误差通道备注：dE 包络与 κ（5.5→53）同步增长，提示 slowdown 棘轮粒度随 dt_step 变粗是主要症状通道；即便 Phase 0 证明主导通道是 DKD 分裂误差而非 κ，两者都由"分辨率恢复 32/轨道"修复，Fix A 均为正确解。

### 2.3 历史一致性评审（修复必须满足的约束）

| 历史决策 | 位置/时间 | 对本修复的约束 |
|---|---|---|
| 节点势改轨道平均口径（semi），避免瞬时 r_sep 在**束缚**轨道近心点高估 U（1/(1-e)） | 15b | 束缚（semi>0）层级**保持** semi 口径，不得改回瞬时 |
| 双曲节点用近心点 q=|a|(e-1) 保守口径（声明"再退化为瞬时 r_sep"） | 15b | q 口径仅对遭遇期成立；尾期无轨道平均可依——瞬时口径正是 15b 留给 degenerate 情形的退化路径，与本修复方向一致 |
| P_eff_min 覆盖全部层级（非叶节点贡献 P·κ / 遭遇时标）；双曲叶 2π/8 | 15c | 遭遇期收缩机制**原样保留**；修复只改尾期行为 |
| 双曲层的 pert 比值改用活距离度量 `calcPertFromMR`（瞬时） | 15c `calcPertRatio` | 同一层级的 U gauge 用瞬时距离与之一致 |
| g 函数本身按瞬时位置求值（`processOuterNode`） | 08-14 | dt=ds/g 中的 g 是瞬时的 → ds 的 gauge 必须与 g 同口径才能让 dt 有界 |
| BTLogH 跳过 g_func 扰动比判据（"always use"，LogH 回退被刻意放弃） | 现行 `checkGFuncCriterionIter` | **不得**用 4→0 切换当修复手段；理论文档 `hierarchical_blogh_plan.md` §2（4↔0 自动切换设计）已过时，需在文档阶段标注 |

**综合结论**：双曲层级不存在"轨道平均"这一概念——它的 U 是瞬态的。因此 gauge 应取瞬时分离 r_inst（与 g 自身、`calcPertRatio` 双曲分支同口径），束缚层级维持半长轴口径。两处同属"口径统一"而非特例补丁。

### 2.4 r_inst 的定义（非用户指定，与 Fix C 完全正交）

**计算来源（零新增参数）**：

```
r_inst = | getMember(0)->pos − getMember(1)->pos |
```

- 在 `multiplyDsByNodePotentials` 双曲分支内计算：成员为叶子取粒子 `pos`，为子树取节点 CM `pos`。这与该函数 semi==0 退化回退分支、`calcPertRatio` 双曲分支（`calcPertFromMR`）、`processOuterNode` 的 U_inst 是**同一数据源、同一写法**——代码模式已存在，Fix A 只是复用。
- **无任何用户输入 / CLI 选项**。Fix C 的 `r_break_crit`（生产中由 PeTar 的 r_group 传入的群破判据，决定"何时不再积分逃逸星"）不进入 ds 公式，二者不共享任何量。

**为什么用瞬时位置而非由 (a,e) 解 Kepler 方程求 r**（对"应由轨道参数定"意见的回应）：
1. **相消要求同源**：dt=ds/g 的 1/r 相消只在 ds 的 U 因子与 g 自身的 U 因子是**同一个函数**时严格成立。`processOuterNode` 用瞬时成员位置求 U_inst（08-14 设计决策 #1/#2：`_bin.r` 过时、距离必须实时从位置算）——ds 若另走 Kepler 解算路径，两条实现间的任何偏差都会漏进 dt。
2. **位置本身即轨道状态**：Kepler 双曲轨道上 r 由状态唯一决定（r=|a|(e·cosh H−1)）；从位置取值与"根数+时间异常解算"数学等价，而后者需存 anomaly（BinaryTree 不存）+ 解双曲 Kepler 方程（成本 + 不一致风险），无额外收益。
3. **"轨道参数定标"已经体现**：r_ref = max(q, r_inst) 中 q=|a|(e−1) 是轨道根数给出的**保守下限**；真实 Kepler 轨道上 r≥q 恒成立（max 平凡退化为 r_inst），仅在扰动或根数陈旧的瞬态配置下 q 兜底。即"根数定标 + 状态取值"，衔接自动、无特例分支。
4. **普适性**：公式对任意层级、任意双曲节点（外逃逸 / 束缚群内 flyby / 链式遭遇）一致生效，无场景参数。

## 3. 相关文件与符号

| 文件 | 符号 | 角色 |
|---|---|---|
| `SDAR/src/AR/information.h` | `multiplyDsByNodePotentials`（~L185-230） | **Fix A 主修改点**：双曲分支改 r_ref=max(q, r_inst)（r_inst 复用 semi==0 回退分支的成员 CM 距离计算） |
| `SDAR/src/AR/information.h` | `calcPertRatio`（~L306） | 根节点 pert_out=0 → 返回 1.0，node_scale 失效 |
| `SDAR/src/AR/information.h` | `calcBLogHDsIter` / `calcEffectivePeriod` | P_eff_min 锚点（2026-08-15c 修复，已工作正常） |
| `SDAR/src/AR/symplectic_integrator.h` | `updateBinarySemiEccPeriodIter`（~L1481-1496，stab>1 门控） | 次要因素：双曲 (a,e) 本为 Kepler 常数，冻结根源在公式输入（§2.2）；Fix B 核验刷新链路时一并检查 |
| `SDAR/src/AR/symplectic_integrator.h` | `checkGFuncCriterionIter`（~L1741-1774） | `g_func_user==4` 跳过判据为**刻意设计**（见 §2.3），不作为修复对象 |
| `SDAR/src/AR/symplectic_integrator.h` | `syncTreeSlowDownAndDs`（~L1678-1720） | 触发 `calcDsAndStepOption` 的入口（orbit 模式每步调用） |
| `SDAR/sample/AR/ar.cxx` | `-m orbit` 循环（~L537-556） | 每步 `updateBinarySemiEccPeriodIter` → `syncTreeSlowDownAndDs`，步长 `integrateOneStep(info.ds, …)`，`dt=ds/g(t)`；Phase 3 增加可选 `--break-check`（生产样态：镜像 `H4::checkBreak` 双曲逃逸分支） |
| `SDAR/tools/ar.py` | `SDARData(g_func=True, N_particle=4, slowdown=True, N_sd=3, time_measure=True)` | 日志读回 |
| `SDAR/docs/hierarchical_blogh_impl_notes.md` | — | 需同步更新 |

依赖关系：PeTar 通过头文件复用 `information.h`/`symplectic_integrator.h`（`PeTar/src/.../ar_interaction.hpp`），修改会传播到 PeTar 构建。

## 4. 分阶段执行计划（v2）

### Phase 0 — 复现、仪器化与误差通道归因（无行为变更）

- **目标**：把 §2 的手工分析固化为可复现脚本；给 ds 公式加调试输出；确认误差主导通道（slowdown 棘轮粒度 vs 算子分裂）。
- **文件**：新建 `SDAR/sample/test/analysis_ustabtri_btlogh.py`（`sdar.SDARData(g_func=True, N_particle=4, slowdown=True, N_sd=3, time_measure=True)` 读回：ds 时间线表、dE 包络 vs r_esc、每输出区间 Nstep 增量 → 步/轨道估算、根轨道 (a,e,E_rel) 分解、κ 序列）；`SDAR/src/AR/information.h` 在 `AR_DEBUG_DUMP` 下打印 ds 分解（ds_prod、P_eff_min、period_prod、各节点 semi/ecc/U_node/node_scale/r_ref）。
- **检查**：脚本输出与 §2.1 表一致（7 个 ds 状态及时刻）；确认 t≥0.0598 后根节点 semi<0、U_node 逐位不变；步/轨道曲线从 ~40 跌至 ~2.6。
- **验收**：一键复现基线；通道归因结论写入脚本输出头部注释。
- **代理**：PeTar Researcher + PeTar Implementer。

### Phase 1 — Fix A：双曲非叶节点 U_node 规范统一（核心）

- **目标**：`multiplyDsByNodePotentials` 双曲分支改 r_ref = max(|a|(e−1), r_inst)，U_node = G·m1·m2/r_ref。
- **文件**：`SDAR/src/AR/information.h`（双曲分支；r_inst 计算复用 semi==0 回退分支；注释写明规范原则"ds 的每个 U 因子与运行时 g 同规范"；束缚分支、叶子系数、P_eff_min 锚点一律不动）。
- **检查**：重编 sample/AR（沿用 `ar.ttl.sd.t.mulpot.cm` 编译选项）；重跑 `ustabtri.sh`。
- **验收**：① 交换 #2 后 ds 随 r_esc 按 ~1/r 衰减（与 Gm₁m₂/r_inst 预测差 ×2 内）；② 逃逸全程内双星分辨率 ≥16 步/轨道（目标 32）；③ 终点 |dE| ≤ 1e-5（基线 3.8e-2）；④ 交换 #1 瞬态（1.75e-5）与相遇期 dE~1e-6 不回退；⑤ 总 Nstep ≤ 780k（基线 774k，预期 +0.8%）。
- **代理**：PeTar Implementer → PeTar Build and Test Maintainer（重编）→ PeTar Simulation Engineer（重跑）。

### Phase 2 — Fix B：ds 时变化的刷新链路核验 + 防御加固

**Fix B 指代**（验证型任务，非新功能）：Fix A 之前，双曲根节点的 (a,e) 是 Kepler 常数 → ds 是"公式级常数"，刷不刷新都一样；Fix A 之后 ds 经 r_inst 随时间变化（逃逸段 ~1/r 递减），其全部有效性取决于**每个积分周期把新鲜 r_inst 算进 ds 并传播到步长数组**。Fix B = 对这条链路的逐环验证与加固，四项：

- **B1 刷新节奏（核心）**：三条调用路径逐一确认"每周期重算"——① `ar.cxx -m orbit` 循环（~L537-556，每步 `syncTreeSlowDownAndDs(true,true)`）；② `integrateToTime` 更新块（symplectic_integrator.h ~L2930-2965，`syncTreeSlowDownAndDs` 在 `if (update_flag)` 之外、逐循环执行）；③ Hermite 直调点（hermite_integrator.h ~L2662）。关键论断：r_inst 取自成员位置（永远新鲜），**不受** `updateBinarySemiEccPeriodIter` 的 `stab>1` 门控影响——该门控只停 semi/ecc 刷新，而双曲 (a,e) 本是 Kepler 常数、陈旧无害。需实测证实（逃逸段打印 ds 序列，应见 ~1/r 递减而非恒定）。
- **B2 传播语义**：更新点 `ds[0]=min(ds[0],info.ds)` + `ds_backup.initial(info.ds)` + `ds_init=info.ds` 对**递减**方向即时生效；确认误差控制器（以 ds_init 为基的 `calcStepModifyFactorFromErrorRatio`）不会在两次刷新之间把 ds 增长越过新鲜物理值（刷新节奏=每周期时窗口≈0）。增大方向的一般语义见未决问题 #1。
- **B3 断言**：`r_ref>0`、`node_scale∈(0,1]`、`U_node` 有限正；双曲分支 `ecc>1`（保留现有）。
- **B4 编译矩阵**：`AR_G_FUNC_MUL_POT × AR_SLOWDOWN_TREE` 组合 + 非 SLOWDOWN_TREE 分支（~L2893）各编译运行一次。

- **检查**：ustabtri + quad_sd2 重跑；ds 序列在重建点无振荡/NaN。
- **验收**：B1 打印证实逃逸段 ds 按 ~1/r 递减；无新断言触发；编译变体间行为一致。
- **代理**：PeTar Implementer + PeTar Reviewer。

### Phase 3 — Fix C：standalone AR 的生产样态（可选 break 检查）

- **目标**：`sample/AR/ar.cxx` 增加可选 `--break-check`：根轨道 semi<0 且 r_inst>r_break_crit 时镜像 `H4::checkBreak` 的双曲逃逸分支，记录释放事件（默认仅记录不终止，"最坏情形压力测试"仍为默认模式）。注：`r_break_crit` 是生产群破判据（PeTar 中由 r_group 派生），只决定"何时停止积分逃逸星"，**不进入 ds 公式**——与 Fix A 的 r_inst 正交（见 §2.4）。
- **文件**：`SDAR/sample/AR/ar.cxx`（CLI + 事件打印；core 库不动）。
- **检查**：ustabtri 开/关 `--break-check` 各跑一遍。
- **验收**：开启时在 r≈(1–10)·r_break_crit 记录到释放事件；两模式 dE 曲线对比入档。
- **代理**：PeTar Implementer + PeTar Validation Analyst。

### Phase 4 — 验证：回归 + 对比矩阵

- **目标**：证明修复不破坏健康层级，并建立可复用基线。
- **场景**：
  1. ustabtri 四件套（logh / btlogh / btlogh_auto / h4）× {默认, `--break-check`}：记录 ds 时间线与终点 dE，对照 §2.1 基线；
  2. quad_sd2 B-B 四星（健康层级，`/home/lwang/localdata/SDAR_BLogH/quad_sd2*.log`）：ds 序列位级或 <1% 一致（外轨道 e=0.1 束缚，不触及双曲分支）；
  3. PeTar 侧：SDAR 头文件被 PeTar 复用，重编 PeTar + `petar.select` 后跑 `test/functional` 冒烟。
- **文件**：验证脚本/记录入 `SDAR/sample/test/`；基线数字写入 impl_notes（不依赖 `~/localdata` 临时文件）。
- **验收**：Phase 1 阈值全部满足；quad_sd2 不变；PeTar 冒烟通过。
- **代理**：PeTar Build and Test Maintainer + PeTar Simulation Engineer + PeTar Validation Analyst。

### Phase 5 — 文档同步

- **目标**：实现与文档一致。
- **文件**：`SDAR/docs/hierarchical_blogh_impl_notes.md` 新增 2026-08-26 条目（双曲 U_node 规范统一 r_ref=max(q,r_inst)、逃逸尾期 ds 跟随、ustabtri 基线表含 v1→v2 归因修正；条目中引用本计划文件 `SDAR/docs/hierarchical_blogh_ds_hyperbolic_gauge_plan.md`）；`hierarchical_blogh_plan.md` §2 的 4↔0 自动切换段落标注"历史设计，现行不使用 LogH 回退"；若 ar.cxx 增加 `--break-check`，同步 `SDAR/README.md`、sample 帮助文本与 `SDAR/.github/skills/sdar-fewbody-integration/SKILL.md`。
- **验收**：PeTar Reviewer 检查代码/文档/测试三者一致。
- **代理**：PeTar Documentation Maintainer + PeTar Reviewer。

## 5. 未决问题 / 决策点（v2）

1. **双向 ds 更新语义**：更新路径 `ds[0]=min(ds[0],info.ds)` 只对当前序列收缩（增大经 `ds_init` 后续生效）。逃逸段单调递减无碍；束缚群内 flyby"离开再接近"是否需要增大即时生效？注意与用户已否决的"重建后 ds 只减不增安全网"区分——那是误差控制器侧的单调钳制，此处是物理输入变化的正常重估，本计划不做任何钳制。是否补对称路径需用户拍板。
2. **孤立 group 根节点 node_scale**：pert_out=0 → `calcPertRatio` 返回 1.0。是否把子层 pert_out（逃逸星对内双星的潮汐）接入根节点 scale 作额外松弛（逃逸远处 ds 回升、步数再省）？默认**不做**（最小修改），Phase 4 若显示生产样态过保守再议。
3. **calcEffectivePeriod 双曲伪周期**（2π√(|a|³/μ) vs 真实近心穿越 q/v_peri，e≈1.7 时差 ~19×）：本数据遭遇期精度已足（dE~1e-6），**不动**，列为后续独立任务。
4. **`--break-check` 语义**：仅记录 vs 终止群积分；默认仅记录。
5. **性能预算**：生产中逃逸尾由 break 终止，理论无长尾；用 PeTar functional 冒烟确认。

## 6. 风险与缓解（v2）

- **混沌敏感性**：重跑轨迹漂移 → 验证用统计阈值（ds 跳变倍数、dE 上界、分辨率曲线、交换是否发生），Phase 0 脚本固化。
- **束缚群内双曲 flyby 行为变化**：接近段 ds 现按 1/r 收缩（v1 无此行为）。分辨率不变量论证表明每内层轨道步数恒定、成本平坦；以 quad_sd2 + ustabtri 交换 #1 前的接近段曲线实证。
- **重建点 ds 不连续**：现有 error-ratio 增长门 + min-clamp 已有处理；确认无 NaN/断言。
- **归因不确定性**：若 Phase 0 证明主导通道非 slowdown 棘轮，Fix A 仍正确（两通道均由分辨率恢复修复），但文档需如实记录最终归因。
- **PeTar 联动**：头文件变更需重编 PeTar；`petar.select` 选定 binary family 后跑 functional smoke。

## 7. 各阶段代理分工速查

- Phase 0: PeTar Researcher + PeTar Implementer
- Phase 1: PeTar Implementer → PeTar Build and Test Maintainer（重编）→ PeTar Simulation Engineer（重跑 ustabtri.sh）
- Phase 2: PeTar Implementer + PeTar Reviewer
- Phase 3: PeTar Implementer + PeTar Validation Analyst
- Phase 4: PeTar Build and Test Maintainer + PeTar Simulation Engineer + PeTar Validation Analyst
- Phase 5: PeTar Documentation Maintainer + PeTar Reviewer
- 全程协调: PeTar Developer
