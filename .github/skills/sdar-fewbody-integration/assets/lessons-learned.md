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
- [Documentation & Units](#documentation--units)
- [Python Tools](#python-tools)

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

### 2026-09-14: quad_sd2 slowdown B--B 的 Δa/a 漂移回归——LogH sum gauge 误乘 κ（85bb681）+ tsyn 截断放大（d86b1e2）；当日误诊两次

**Mistake**: 诊断链上连续犯两个错。(1) 第一轮：把 `Large_energy_error` 初始化爆发（全部消息 t<4.2e-4）+ 慢盘 I/O 节流误判为"运行卡死"，且用 **1e-3 量级的粗精度**比较各版本 a1/a0（打印 4 位有效数字），在 9 月各版本"彼此一致(≤3e-3)"后就宣布"物理健康、无回归"——**基线选错**：正确基线是参照数据的生产二进制（Aug 1-3, `hermite_old`/e55ff89），其 Δa/a 长期 ~4e-8，比 3e-3 小四个量级。(2) 第二轮全精度对比才确认真实回归：85bb681（08-29）给 LogH sum gauge 的椭圆有效周期**乘了 slowdown 因子 κ**（新 `calcLogHSumGaugeIter` 替换从不乘 κ 的 `calcDsKeplerBinaryTree`），ds 放大 κ 倍（quad_sd2 binary1 κ≈19、binary2 κ≈509）→ 每区间能量误差贴着 -e 检查阈值 → binary1 Δa/a 线性漂移（6.4e-6/t，t=4153 外推 2.7e-2）；d86b1e2 首步预测截断在此基础上再降 37% AR 步数 → 6 阶误差 ×16 → 漂移 9.5e-5/t（15×）。用户指出"老版 Δa/a<1e-7、新版>1e-3"后定位完成。

**Root cause**: (a) 85bb681 的验证集（ustabquin 等）全部无 slowdown（κ=1），P·κ 与 P 相同 → 回归只在 κ≫1 的 Hermite slowdown 组路径暴露；(b) 诊断时"同类版本互比"代替"参照基线比对"，且比较精度（4 位）掩盖了 1e-7 vs 1e-3 的四个量级差异；(c) 消息条数/推进速度是运行学症状，精度结论必须来自 .log 的轨道量全精度重算。

**Prevention rule**:
1. 精度回归判断必须满足两个条件：(i) 基线 = 参照数据的生产二进制（git worktree 重建或 `~/bin/hermite_old` 类存档），(ii) 比较用全浮点精度（%.2e）的轨道量（Δa/a、dE），禁止 4 位打印的"看起来一致"；
2. 改动 ds 估计器的任何公式，验证矩阵必须包含 **slowdown κ≫1 用例**（Hermite 组路径）而不只是无 slowdown 的 AR -m 路径；κ=1 时新旧公式逐位相同不代表 κ>1 时安全；
3. 2026-09-14 修复：`calcLogHSumGaugeIter` 椭圆层去掉 κ 乘子（恢复 08-29 前 LogH 语义；κ=1 逐位不变，ustabtri/quin/quad 论文数字不受影响；BTLogH 路径的 P·κ 是 Aug 论文数据既有基线，保持不动）。修复后 quad_sd2 h4：0 消息、Δa/a≤8.5e-8 至 t=600、比 AUG 省 26% AR 步（tsyn 截断收益保留）；κ=56000 的 sd1e2 配置 60 秒到 t=1024 正常；
4. "卡死"判定先分离 I/O（tmpfs 复跑）再看消息时间直方图（见同日 I/O 节流误判）——但 I/O 结论不能引申为精度结论，两者独立验证。

---

### 2026-09-14: BTLogH 乘积 g 在 2 粒子组（pair 路径）发散——验证矩阵系统性盲区

**Mistake**: BTLogH → Hermite/PeTar 迁移接线完成后做 gf1 冒烟，`--g-func 1` 在**所有含 2 粒子组的场景**立即发散：AR standalone 孤立双星 dE→1.4e15 且输出失控（26GB）；H4 孤立双星（κ=64000）、quad_sd2 双 2 粒子组（κ=21.8）、ustabtri 内双星组（**κ=1**）均 NaN abort。历史 BTLogH 验证矩阵（ustabtri/ustabquin/quad B--B--S、Scheme G 回归 R1-R7）**全部是 3+ 粒子单组**，走树路径（`integrateOneStep` 树版本）；2 粒子走独立的 pair 重载（`integrateOneStep` 2 体版本，KDK/non-KDK 两分支），其乘积形式分支从未被任何验证覆盖。"2 粒子时 BTLogH 的 g ≡ 单对势（=LogH）"的代数论断让人误以为无需专门验证——实际上 pair 路径的 drift 因子积分（`dgt_drift_inv *= gt_kick_inv_`，乘入**含 κ** 的 gt_kick_inv_=U/κ）与 gtgrad（乘积分支给**不含 κ** 的 ∇ln U）口径不一致，`gt_drift_inv_` 偏离 g_kick 量级并失控（dump 实测 2.5e-11/3.95e10 vs g_kick O(0.1)）；κ=1 亦失败，κ 大则恶化更快。

**Root cause**: (a) 组粒子数决定代码路径（2=pair 重载 / 3+=树路径），验证矩阵按"系统拓扑"（三体/四体/五体）设计而未按"代码路径"覆盖——全部多体测试恰好都绕开了 pair 重载的 g_func 分支；(b) 代数等价（g_product=U/κ=g_sum）只保证 g 值相同，不保证 kick/drift 分解与 gt 演变的数值实现一致；(c) H4/PeTar 中 2 粒子组恰是最常见组型（孤立双星），盲区落在生产最重的路径上。具体两缺陷：① `applyStableCheckAndSlowDown` factor=1 分支不同步 `sd_root.period`（2 粒子根=叶子永不进 `calcBinaryTreeSlowDown` 内循环）→ BLogH ds 估计器读到构造默认 DBL_MAX → ds=inf（LogH 口径在线算 P 从不读该成员故隐形）；② pair 路径 `dgt_drift_inv *= gt_kick_inv_` 平白多乘 U——pair 路径直接消费 interaction 层原始 ∇U gtgrad（÷U 的 ∇ln U 换算只在树路径包装层做），**圆轨道 v·∇U≡0 使该缺陷在两体圆轨道测试中完全不可见**。

**Prevention rule**:
1. g-func（或任何按宏分支的算法路径）的验证矩阵必须按**代码路径**而非系统拓扑设计：至少各包含一个 2 粒子（pair 重载）与 3+ 粒子（树路径）用例 × {κ=1, κ≫1} × {圆轨道, 偏心轨道}；"代数上等价"不能替代运行验证；
2. **数学上两模式恒等的路径不应有任何模式分支**：两体时乘积 g ≡ LogH 求和，`integrateTwoOneStep` 中的 `AR_G_FUNC` 分支本身即错误（最终修复=整段删除）。曾以 `*=κ⁻¹` 试图"对齐口径"——κ=1 时恰为无操作故通过两体验证，κ=21.8 仍失败；等价性应通过删除分支达成，而非补系数；
3. 圆轨道是危险的验证用例：v·∇U≡0 使 drift 梯度类缺陷完全静默；两体最小验证必须含偏心轨道（e=0.9 即可在数步内暴露）；
4. H4 场景的 gf1 冒烟最小集：孤立双星（`bin2.dat` 型输入）+ 一个含 2 粒子组的少体系统（quad_sd2 拆组型）；AR standalone 孤立双星是最小复现载体（`ar.btlogh.ttl.sd.cm --g-func 1 -t 1 <2粒子输入>`，几分钟内 dE 爆炸）；
5. 失控运行先设输出上限再诊断（`head -c` / `-n <nstep_max>` / timeout），本次 26GB 日志属可预防事故；
6. 诊断链：先查初始化行（t=0 的 ds/Gt_drift/SD 块）再查首步轨迹——本轮 ds=inf 在 t=0 行即肉眼可见；修复后验收标准=gf1≡gf0（AR 两体逐位；H4 两体末位差来自两条 ds 公式浮点路径，属预期）。

### 2026-09-14: `sample/input/fewbody_hermite.sh` 与当前 experiment 分支不兼容（存量）

**Mistake**: 按 sample 脚本注释（"Runtime: < 0.1 second"）直接跑 `hermite -t 1.0 -o 2 -G 1.0 triple.stable.lowm3`，plain 与安装版二进制均 NaN abort（`Assertion !std::isnan(integration_error_rel_abs)`），一度干扰新改动的冒烟判断。

**Root cause**: 脚本注释基于 2026-07-21 验证；experiment 分支后续演化（κ 修复、着陆截断等）改变了该输入在默认参数下的行为，脚本未同步。

**Prevention rule**: sample 脚本失效时先跑"已安装生产二进制 + 同参数"对照（本次 `~/bin/hermite` 同样 abort → 存量问题，与新改动无关），再寻找替代冒烟输入（ustabtri H4 命令见 plan §Phase 5）；skill/脚本文档在分支演化后需重新验证。

### 2026-09-15: ds 地板 `max(|t_end|,1)` 钳制在代码单位时间 <1 时退化成 time_error 口径，误杀 812-ulp 着陆步

**Mistake**: PeTar 生产运行（mcluster Kroupa IC，t=0.0039 Myr，二体组 id 34/35）触发 `Error! adaptive ds below time resolution in integrateToTime` 并 abort。dump 显示 H=-8.9e-16（扩展哈密顿量在舍入水平守恒）、gt_drift_inv=0.098（远非 pericenter 大 gt）——AR 积分本身完全健康，却因 ds 在 22 步内从 2.9e-5 塌缩到 6.9e-17 触发地板。直接原因：`dt_step_floor = 4·eps·max(|t_end|, 1.0)` 的 `max(,1)` 钳制——t_end=0.0039=2^-8 时地板虚高 256 倍，被杀的那步 dt=7.04e-16 = 812 个 ulp(time_)，完全可以推进浮点时间且大概率一步进窗收尾。更深一层：`max(,1)` 使地板在一切 |t_end|<1 的运行里恰好等于 time_error 的 clamp 项——2026-09-13 ustabquin 教训明令禁止的 "time_error 口径地板" 从后门回来了。塌缩驱动是着陆重标定 `i==0` 分支的绝对时间比 `_time_end/time_table[k]`（隐含 t_start=0；该组实际从 t≈0.00195 起积，比值错 ~2 倍），在着陆循环里逐次把 ds 压向地板。

**Root cause**: (a) 可表示性判据的尺度取错对象——`time_+dt` 是否推进取决于 ulp(|time_|)，与常数 1.0 无关；Myr 单位下星团早期演化整段都在 t<1 区间，正好落进钳制的放大区；(b) 着陆算术多处隐含 "组积分时间从 0 开始" 假设（`_time_end/time_table[k]`），而 Hermite 调用方传入的组起始时间一般非零；(c) 2026-09-13 的 Prevention rule 第 3 条把 "4 ulp·max(|t_end|,1)" 写成了规则本身，把当时发现的症状当成了正确口径固化下来。

**Prevention rule**:
1. 防御地板的正确口径是 `4·eps·max(|time_|, |t_end|)`（当前时刻与目标时刻的较大者，无 1.0 钳制）；修复后同 IC 完整跑完 100 Myr，能量误差 -7.5e-6；
2. 时间量之比必须以"步起点"为参考（`(t_end - t_step_start)/(table[k] - t_step_start)`），禁止绝对时间直除；循环顶部保存 `time_step_start` 供着陆分支使用；`time_table[k] - t_step_start <= 0`（首子步倒退）时保守减半而非套公式；
3. "X 口径不能用作判据"类教训在写成 Prevention rule 时，要检查新判据在全部参数区间（此处 |t|<1 与 |t|>1）是否真的避开了被禁口径——本条与 2026-09-13 第 3 条冲突即是固化症状的代价；
4. 负时间步缩步分支（pre-sync 与 time-sync 两处）已加 streak 上限（连续 2 次后停止缩减并显式报 gauge 问题）：ds 幅值改不了 gt 决定的符号，restore 后同状态重试必然同样结果——本 IC 实测 0 次触发，塌缩来自着陆比值而非负步长，但保护留存以防其它场景；
5. 诊断此族崩溃时先读 dump 里的 H 与 gt_drift_inv：H≈0 + gt 平缓 ⇒ 积分健康、控制器/着陆逻辑有病；gt 巨大 ⇒ 才是真正的 pericenter 正则化场景；
6. `-DAR_COLLECT_DS_MODIFY_INFO`（PeTar Makefile 中注释保留）一次重编即可区分 ds 修改来源（Large_energy_error / Negative_step / Negative_step_tsyn），先插桩再动控制器逻辑。

---

## Documentation & Units

### 2026-09-17: `G` 常量规则在 SKILL 内部自相矛盾——单位制静默出错风险

**Mistake**: `SKILL.md` 的 Non-Negotiable Rules 要求"使用 `G_MSUN_PC_MYR = 0.00449830997959438`"，而同文件的选项表与 Common Pitfalls 写"`-G` 默认 1.0（Henon），需确认"。两处并列时，agent 在 Msun/pc/Myr 场景可能直接套用 0.004498（正确），但在 Henon 场景也可能沿用它（错误）；更危险的是 `-G` 与 `-u` 不一致时 SDAR 不报错，只静默给出错误的能量与轨道。同一文件内的配置审计（2026-09-17）才发现。

**Root cause**: 规则按"写作时最方便的角度"分散落笔——Non-Negotiable 层从"该用什么值"表述（只列出一个值），选项表层从"默认值是什么"表述（列出另一个值），两个片段各自孤立看都正确，从未被并列对照过。

**Prevention rule**: 涉及单位制的常量必须写成"合法值集合 + 必须与运行一致"的单一声明，而不是若干各自正确的片段。现行口径：`-G` 只有两个合法值 —— `G_HENON = 1.0`（unscaled / Henon 单位，C++ 默认）与 `G_MSUN_PC_MYR = 0.00449830997959438`（`-u 4`，Msun/pc/Myr）；C++ 运行与 Python 分析必须一致，否则结果静默错误。配置文件审计时，应专门对照"同一常量在不同章节的取值"。

---

## Python Tools

### 2026-09-27: `readArray` 列数不匹配从 warning 改为抛错——列过剩静默丢列与列不足裸 IndexError 都是坑

**Mistake**: `DictNpArrayMix.readArray`（`tools/base.py`）在列数不匹配时只 `warnings.warn` 后继续执行：列不足时随后在 `_dat[:,icol]` 抛出无上下文的 `IndexError`（PeTar Pal5 会话读 2021 年 `data.status` 的直接崩溃点）；列过剩时静默丢弃多余列、不报任何错——错误的 reader kwargs（如 `external_mode` 选错导致中段插列）可产生静默错位数据，比崩溃更危险。

**Root cause**: 校验失败路径只照顾了"可能是有意的超宽读取"这一善意情形，未评估两种真实失败模式的下游后果；且 warning 文本不含类名与 reader 实际 kwargs，用户/agent 无法自助诊断（按文档穷举 kwargs 也试不到版本相关的 `spin_3d`）。

**Prevention rule**: 现行为：`ncol_check=True`（默认）且不匹配时抛 `ValueError`，消息含类名、文件列数、类列数+offset、reader initargs 与 legacy 提示（pre-2024-12 输出用 `spin_3d=False`）；有意读取列子集时先精确宽度切片或传 `ncol_check=False`（嵌套 readArray 本就传 False，不受影响）。通用规则：数据布局校验失败应默认 fail fast，"宽容继续"必须显式 opt-in。

### 2026-09-27: `fromfile` 字节错位默认抛错（实装 strict_mismatch）+ `loadtxt` 二进制检测——文档曾引用未实现的参数

**Mistake**: `fromfile` 对文件字节数与 dtype itemsize 不对齐只 `warnings.warn` 后继续：PeTar Pal5 管线漏 `-i bse` 时错位数据传导成 NaN → `np.histogram` "bins must increase monotonically"，崩溃点距根因三层；`loadtxt` 读二进制 `data.core` 抛裸 `UnicodeDecodeError`。且 PeTar patterns 文档早已写了 `strict_mismatch=False`（DSM interrupt 尾部填充场景），但该参数在代码中从未存在。

**Root cause**: 与 readArray 同类的"校验失败宽容继续"；文档先行描述了计划中的参数而未实现，形成 doc-code 脱节。

**Prevention rule**: `fromfile` 现默认 `strict_mismatch=True`：错位抛 `ValueError`（类名、字节数、itemsize、reader initargs、kwargs 提示）；已知填充/截断文件（DSM `data.interrupt`、petar.data 崩溃恢复 partial）显式传 `False` 读完整记录。`loadtxt` 捕获 `UnicodeDecodeError` 转为明确 `ValueError`（提示改用 fromfile）。新增文档参数必须同 change 实装，否则在文档显式标注"未实现"。
