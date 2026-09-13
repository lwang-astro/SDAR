# Lessons Learned — SDAR Development & Integration

This file records mistakes, gotchas, and pitfalls discovered during SDAR agent sessions
(integrator development, time synchronization, build/baseline methodology).
Each entry documents: what went wrong, why, and how to prevent recurrence.

**Lifecycle**: New entries are added by agents automatically after encountering issues.
Periodically reviewed → verified entries are elevated to `SKILL.md` as hard rules.

---

## Categories

- [Build & Baseline Methodology](#build--baseline-methodology)
- [AR Integrator & Time Synchronization](#ar-integrator--time-synchronization)

---

## Build & Baseline Methodology

### 2026-08-27: Derived preprocessor macros defined in an includer are invisible to included headers

**Mistake**: `AR_G_FUNC_MUL_POT_FAMILY`/`AR_G_FUNC` were defined in `symplectic_integrator.h`, but `information.h`
(included earlier by it) also gates code on them → all family branches silently compiled out, `ds` took the
wrong formula; the only symptom was last-digit float drift in the first output row.

**Root cause**: old macros came from the command line (`-D`), visible everywhere; refactored *derived* macros
inherit include order.

**Prevention rule**: shared derived macros live in their own header (`src/AR/g_func.h`) that every consumer includes
first; when verifying a pure refactor, insert one-shot stderr prints to confirm the expected branch is taken
before trusting bit-level diffs.

### 2026-08-27: Stale binaries in `build/` are not a valid regression baseline

**Mistake**: Phase-0 baselines were taken from existing `sample/AR/build/*` binaries; they differed from
HEAD-rebuilt binaries in 1000/1002 lines (hot fixes were never rebuilt), nearly misattributing a real bug
to binary staleness and hiding the include-order bug above.

**Prevention rule**: for refactors, generate baselines from binaries rebuilt via `git worktree add <tmp> HEAD` —
never trust pre-existing build outputs.

### 2026-08-27: Bit-level output diffs must exclude wall-clock columns

**Mistake**: `SDAR_TIME_MEASURE` prints `Total(s)/Int(s)` (per-step CPU time) which differ every run; they accounted
for 1000/1002 "differences" until physical columns were compared selectively.

**Prevention rule**: normalize/ignore timing/profile columns before declaring divergence (normalize/grep out
`Total(s)`/`Int(s)`-style columns; the SDAR g-func refactor acceptance log in
`SDAR/docs/hierarchical_blogh_impl_notes.md` used this approach).

### 2026-09-13: 自适应步长+混沌系统的步数对比跨构建旗标无效（-O0 vs -O2 差 3 倍）

**Mistake**: 对比 ustabtri `-m full` 步数时，把 debug 构建（-O0 -g，12:56 前的 `~/bin`）产生的日志（S32: 19407 步）与 release 构建（-O2）基线（60542 步）直接对比，得出"步数膨胀 3.8×"的量化结论；orbit 回归也因新旧构建不同被误判为"位级不一致"（实为混沌发散+输出行数不同）。另有两次把墙钟计时列（`Total(s)`/`Int(s)`）差异误当作物理差异。

**Root cause**: 步数受两个放大器控制：(1) 误差控制器阈值（1e-10）附近二值决策——微小的浮点舍入差（-O0/-O2 的 FMA 收缩等）翻转减步决策；(2) 系统混沌（unstable triple 的交换事件）放大一切初始差直到轨迹完全不同。两者叠加使步数对构建旗标、甚至二进制版本极端敏感；orbit 模式输出按步（非按区间）发生，尾部大步跨越多个区间时行数 < 区间数，行数差异又被误读。

**Prevention rule**:
1. 自适应步长+混沌系统的 A/B 对比必须**同一二进制、同一构建旗标**；改代码前先 `git stash` 构建基线二进制留存，改后重建对比（本次 ar.base 流程），不要复用历史日志当基线；
2. `cmp`/diff 判定位级回归时，先剔除计时列（标题行含 `Total(s)`/`Int(s)` 的列号），再比较物理列；
3. 报告此类系统的步数时注明二进制 md5 与构建旗标；能量误差（dE）在同一轨迹内可比，跨轨迹只能比量级；
4. 在混沌系统上验证机制假设（如"-o 影响步数"）必须加**非混沌对照**：ustabtri 上 regular 步随 cadence ±25% 无规律摆动，稳定三体对照显示平坦（±1.6%）——摆动是同步网格扰动引发的轨迹散布，不是机制；机制假设会被混沌散布伪造或掩盖；
5. 改动 ds 恢复/重置逻辑前，先确认它不是历史事故的防护再动手：`integrateToTime` 每次调用入口对工作 ds 的**无条件**向上重置，防的是混沌 dance 中误差长期高于阈值导致塌缩 ds 跨区间滞留（"ds 减小后回不去"）；committed 层的单向 valve 防的是 ustabquin S128 式停滞（24.4M 步 200× 过分辨）。条件重置/对称化这类"优化"会拆掉保险——本次收窄版 C 因此撤销。

---

## AR Integrator & Time Synchronization

### 2026-09-13: -m full 断言 `_ds>0` 根因是 live/backup H 评估器公式不一致，而非时间同步逻辑

**Mistake**: unstable triple 在 BTLogH+TTL 构建下 `-m orbit` 正常、切换 `-m full`（时间同步）后 `ASSERT(_ds>0)` 崩溃（`integrateOneStep` 收到 ds=0）。初步推断指向时间同步的 ds 缩放分支或负时间步分支（gt_drift_inv_ 变负）——gdb 条件断点 `dt<0` 证明负时间步分支**从未执行**，推断错误。真实机制：`integrateToTime` 的逐步误差检查 `|H - H_bk|` 中，live 端 `getHSlowDown()` 用 `calcH`（LogH 形式 log(pt)−log(−U)，误差除以 |U|~9e2），backup 端 `getHSlowDownFromBackup()` 在 `#ifdef AR_TTL` 分支遗留旧 TTL 形式 `(E−Eref)/gt_kick_inv`（误差除以 gt~2.4e5）。同一状态两公式相差 ~265 倍尺度，`|H−H_bk|` 含一个 **ds 减小永远无法消除的常数偏移**；当累积能量误差恰好使偏移压在 `energy_error_relative_max` 阈值之上时，能量误差分支每 2 次迭代无下限地减半 ds，~1074 次后 ds 下溢为精确 0 → 断言。`-m orbit` 不触发仅因该检查只存在于 `integrateToTime`。

**Root cause**: live 端 H 评估器曾整体切换为 `calcH`（TTL 形式全部被注释），但 backup 端两个孪生函数的 `#ifdef AR_TTL` 旧分支未同步删除——部分重构遗留。诊断被误导的原因：症状（仅 -m full 崩溃、崩溃前 ds 被 cost valve 抬升 11.5 倍、gt_drift 塌缩）都指向时间同步/valve，但那些只是把能量误差推到阈值边缘的背景条件，不是崩溃机制本身。

**Prevention rule**:
1. 改任何"成对评估器"（live / FromBackup、当前值 / 备份值）一侧的公式时，必须 grep 孪生实现并同步修改，二者用同一表达式；
2. 诊断"自适应步长死亡螺旋"类崩溃（ds 反复减半→下溢）时，先用 gdb 在修改点打印轨迹（每次 ds 修改输出 t/ds/H/H_bk/gt），确认误差信号是否随 ds 减小而下降——不下降即为常数偏移（公式不一致或状态函数混入非状态量），不要先怀疑步长调整逻辑本身；
3. gdb 条件断点是廉价证伪工具：先对初步假设设条件断点（如 `dt<0`），未命中即可排除，避免在错误分支上深挖。

### 2026-09-13: integrateToTime 着陆分支隐含"sorted cck ∈ (0,1]"假设；Yoshida-B(-8) 违反且会计口径随改动漂移

**Mistake**: 给 `integrateToTime` 加预测式截断时，用 `getSortCumSumCK(cd_pair_size-1)` 当作"ΣCK=1"使用——对正阶与 -6 恰好成立（偏和最大值=1），但 Yoshida 第二法（负系数，如 `-k -8`）的 sorted 表末项是最大**偏和** 1.406 而非总和，dt_pred 高估 40%。同轮 review 还发现：着陆分支 `ds[ds_switch] *= cck_prev` 在 -8 下直接产生负 ds（cck[0]=−0.406）——`-m full`+`-8` 在基线就崩（`ASSERT(_ds>0)`），属预存在缺陷。另两处自查教训：(1) ds 地板首版以 `ds_init`(=info.ds) 为基准，而估计器在密近相遇会给垃圾值（实测 2.6e11），健康步长被误报"下溢"——地板必须自洽（时间分辨率 `time_error*gt` 口径）；(2) CLI 校验首版用 `cck<=1.0`，把浮点累加出 `1+1ulp` 的正阶 +8 误拒——致命条件只是 cck≤0，边界比较必须给舍入留缝。

**Root cause**: 辛系数表有两种语义（时序偏和 vs 总和；比率乘子 vs 归一化权重），代码消费方各自隐含假设且无文档；`n_step_tsyn` 的语义还会随"完美着陆不进入 tsyn 计数分支"而改变（跨版本只能比总步数，不能比 tsyn 步数）。

**Prevention rule**:
1. 凡改 `integrateToTime` 的着陆/截断逻辑，必须用 `-k -8` 冒烟（暴露 cck 假设违规）+ `-k 8` 冒烟（暴露边界 ulp 问题）+ 两体输入（`integrateTwoOneStep` 路径无其他覆盖）；
2. ΣCK 按构造恒为 1（已核对生成器全部分支：正阶递归权重和为 1，-6/-8 硬编码表和为 1，偏差 ≤2e-15 纯舍入），代码直接用常数 1.0 即可；禁止把 sorted 表末项当总和（负系数阶末项 = 最大**偏和** ≠ 1，-8 实测 1.406）；
3. 自适应循环的防御地板锚定在**浮点推进下限**（单步 dt < ~4 ulp·max(|t_end|,1)，与 time_error clamp 第二项同口径），不锚定在可能为垃圾的估计器输出（ds_init/info.ds），也**不能**用 time_error 当地板（见下一条 ustabquin 误杀）；
4. 报告/论文引用 `n_step_tsyn` 时注明二进制版本——着陆算法改动会改变其会计口径。

### 2026-09-13: ds 地板用 time_error 口径误杀健康着陆的最后一步（ustabquin dance）

**Mistake**: ds 防御地板第二版用"单步 dt < time_error 即无法收敛"作 abort 条件。ustabquin `-m full` 在 dance 相（t=0.1006，gt_drift_inv~4e7）崩溃：gdb 轨迹显示 dt_end=2.5424e-14 仅比 time_error=2.5e-14 **大 1.7%**，landing 几何收敛（每步 ×0.098）从 ds=4.17e-6 直接跳到 4.07e-7、跳过目标 ds*≈1.06e-6——而那一步 dt=9.8e-15 < dt_end 属**欠冲**，time_ 推进后剩余 1.56e-14 < time_error，正好落入 finish window 完成着陆。地板把这条完全正常的收尾路径拦下 abort。同旗标 debug 基线对照确认基线在该点毫无困难。

**Root cause**: "dt < time_error ⇒ 窗口无法闭合"推理错误——finish window 判定是 time_ 落入 [end−err, end+err]，dt_end 略大于 time_error 时一步欠冲即可进窗；landing 收敛序列的几何步进（×0.1）远粗于 time_error 尺度，"跳过目标后欠冲进窗"是常态而非病态。真正的不可收敛只有浮点级：dt 小于 ulp(time) 时 time_ 不再推进。

**Prevention rule**:
1. 防御性 abort 的阈值必须取"物理上确实无法继续"的下限（浮点推进下限），不能取"看起来太慢"的工作阈值——后者会把合法收敛路径误判为卡死；
2. 新增 abort 守卫后，除常规回归外必须跑一个"极端 gt + 着陆"场景（ustabquin dance 即是最小复现：gt~4e7、dt_end/time_error≈1.02）；
3. 诊断此类误杀：gdb 在守卫处打印 dt_end 精确值与 time_error 比较——dt_end 仅比阈值大百分之几时，几乎必然是阈值设计错误而非系统卡死。

---
