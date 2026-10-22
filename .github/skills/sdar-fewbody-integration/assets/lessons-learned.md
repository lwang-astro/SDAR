# Lessons Learned — SDAR Development & Integration

This file records mistakes, gotchas, and pitfalls discovered during SDAR agent sessions
(integrator development, time synchronization, build/baseline methodology).
Each entry documents: what went wrong, why, and how to prevent recurrence.

**Lifecycle**: New entries are added by agents automatically after encountering issues.
Capture rules (grep-dedup before adding, extend in place, category + chronological placement,
Categories maintenance) are owned by `AGENTS.md` → "Lessons-Learned Capture".
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

### 2026-09-28: 字节级 diff 失败的第一分诊=同二进制复跑（重蹈 09-13 规则 2 + 新增分诊步）

**Mistake**: symplectic_integrator.h 清理后 diff 改动前后的 triple 冒烟日志报 DIFFERS，误以为控制流被改坏；实为计时列（本例第 17/18 列）墙上时间噪声——同二进制两次复跑也互相 DIFFERS，物理列逐位一致。09-13 规则 2 已写"先剔除计时列"，本次未先查教训文件直接 diff，重蹈覆辙。

**Root cause**: SDAR_TIME_MEASURE 默认开启且计时列列号随输出模式变化，肉眼 tail 看不出差异在哪列；教训存在但未进入 diff 前的例行动作。

**Prevention rule**: SDAR 日志 A/B 对比前的例行动作：(1) 同二进制复跑一次，区分"非确定性/计时噪声"与"真实行为改变"；(2) 按列剔除计时列后再 diff 物理列（09-13 规则 2）；(3) 判定"行为是否改变"只以物理列逐位一致为准。

### 2026-10-22: VERSION 的 `e` 后缀靠人工追加是陷阱——先归因脚本再读脚本源码（假 bug 记录）

**Mistake**: 提交前执行 `bash get_version.sh` 后又手工 `echo '469f' > VERSION` 覆盖（正确值 470e），随后把错误值归因为"脚本按十六进制自增了带 e 后缀的版本号"，并据此向用户断言脚本有 bug——脚本实际从不读 VERSION（只做提交计数十进制加法并写出裸数字），错误完全出自覆盖命令本身。诊断未读脚本源码，凭结果反推了不存在的机制。

**Root cause**: 版本号由"脚本产出裸数字 + 人工追加 `e` 实验后缀"两步拼成，人工步骤游离在脚本之外，拼写/进位错误直接写进文件；且当结果异常时，"脚本产出被人工覆盖"这一因果链肉眼不可见，极易归因到自动化环节。

**Prevention rule**:
1. 版本后缀已并入脚本（2026-10-22 起）：`get_version.sh` 写出 `<提交数+1>` 并自动保留既有 `e` 后缀，提交前只需运行脚本+stage，删除人工追加步骤（PeTar `tools/get_version.sh` 同步修复）；
2. 断言"工具 X 有 bug"之前必须先读 X 的源码并复现 X 的原始输出（本次脚本裸输出是正确的 `470`，读三行源码即可证伪）；
3. 检查"自动化输出是否被污染"时先查命令历史里的重定向/echo 覆盖，再怀疑工具本身。

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

### 2026-09-28: 负 dt 的 streak≥2 停止缩减/abort 属过度防御，已撤销（推翻 2026-09-15 规则 4）

**Mistake**: streak 上限（连续 2 次负 dt 后停止缩减并报 gauge 问题，d13b857 降级为 ASSERT）在 PeTar 实际模拟中多次触发 abort；生产构建忽略 ASSERT 后同一批运行全部完成、能量误差无异常——判据在杀健康运行。tsyn 分支的 streak<2 门还会在 streak≥2 后以同一 ds 确定性重放同一状态，直到 step_count_max abort，复刻同一误杀；HARD_DEBUG 构建里该 ASSERT 仍会在良性 streak≥3 上 dump+abort，调试构建比生产更脆。

**Root cause**: "restore 后同状态、g₀<0 ⇒ 缩减改不了符号"的极限论证只在步首 gauge g₀ = ekin − etot_ref < 0 时成立，而真实轨迹状态上 g₀ = U − ΔE > 0（无 bug 恒成立）；负 dt 是大 ds 把组合子步的 gauge 采样到偏离轨迹中间态的数学现象（hyperbolic pericenter + 高阶组合的大正负 CK 系数），缩 ds 使子步收敛回真实轨迹、dt 恢复为正，所需次数可以 >2。streak≥2 被误读为"卡死"，实为"还没缩够"。

**Prevention rule**: 负 dt（pre-sync 与 tsyn 两分支）一律无条件缩 ds；终止只交给既有兜底（ds 地板/下溢 abort、step_count_max abort——仅真 bug 类 g₀<0 触发，届时 dump 中 gt_drift_inv<0 可直接辨认）。从极限论证推出"不可能自愈"的结论，必须有一起实测卡死案例支撑才能固化为 abort 判据。

### 2026-09-28: correctSlowDownEnergy 快照取在调用方更新 epot_sd_ 之后——慢同步每事件漏记 Δepot，etot_sd_ref 漂移致 log() 越界 NaN（fewbody_hermite.sh abort 根因，已修复）

**Mistake**: `fewbody_hermite.sh`（`hermite -t 1.0 -o 2 -G 1.0 triple.stable.lowm3`）稳定复现 `ASSERT(!ISNAN(integration_error_rel_abs))` abort（t≈0.33；9 月 14 日安装版同炸）。-O0+gdb 诊断链：条件断点监视恒等式残差 R=(ekin_sd_+epot_sd_)−etot_sd_ref_（应≈0）→ 首次 <−1e-12 于 t=0.0617 → 残差出现时刻前停下 → `etot_sd_ref_` 硬件观察点+条件 → 单次写入泄漏 −9.8e-9，写入者为 correctSlowDownEnergy 全量路径 `etot_sd_ref_ += de_sd`（由步中 syncTreeSlowDownAndDs←integrateToTime:2961 触发）；数百事件后 R→−1.87e-4，g_sd=ekin_sd_−etot_sd_ref_ 过零，H=log(负数)=NaN。

**Root cause**: syncTreeSlowDownAndDs Step 7 先执行 calcAccPotAndGTKickInv()（更新 epot_sd_，含 κ_in 变化的嵌套标定跳变），之后才调 correctSlowDownEnergy——后者入口取 before 快照时 epot_sd_ 已是新值，Δepot 永远量不到，de_sd 只补了 Δekin（calcEKin 在函数内、快照之后）。触发频率放大器：hermite 组树根节点 stab 保留构造默认 DBL_MAX（该初始化路径未跑 stableCheckIter；物理真值≈0.015），updateBinarySemiEccPeriodIter 的 `stab>1.0` 门恒真 → 每个轨道周期执行一次步中 sync；ar 驱动 stab 正常（<1）无周期 sync，故修复前后 ar 输出逐位不变。诊断中被推翻的两个中间假设（避免重蹈）：(1) "ar 样例无 AR_COLLECT_DS_MODIFY_INFO 导致零事件是插桩假象"——补装插桩后确认 ar 确实无事件；(2) "缩放路径根-κ 与嵌套-κ 约定混合是根因"——条件断点证明该路径在整个失败运行中从未执行（`_force_full_recalc` 恒 true）。

**Prevention rule**:
1. 能量记账 de_sd 类增量的 before 快照必须在任何重算之前取——"调用方先改、被调方入口才量"的结构必然漏项；修复：syncTreeSlowDownAndDs 在 force sync 前快照 ekin_sd_/epot_sd_ 并作为必填参数传入 correctSlowDownEnergy；
2. 修复验收：hermite 冒烟跑通 t=1.0，Large_energy_error 从每块步爆发（err/max 1e6-1e7）降至 4 次边缘事件（1.7e-10~2e-9，一次减半即恢复）；ar.logh.ttl.sd.cm 与 ar.logh.sd.cm 修复前后逐位一致；
3. 慢速记账类 bug 的诊断模板：恒等式残差条件断点定位首发时刻 → 该时刻前停 → 成员观察点+条件继续 → bt 即写入者；比纯静态读码快一个量级（本次 6 次 gdb 会话定位）；
4. 遗留（未修）：hermite 组树创建路径 stab 未初始化（DBL_MAX）把稳定层级当永不稳定处理，每周期空跑 sync——性能损耗，且正是它把记账漏洞暴露成 NaN；后续应在组初始化处补 stableCheckIter。

### 2026-10-19: κ 时变系统的每步 H 误差是速率型(∝dt)——ds 控制器误判为截断误差塌 5 个量级并被 0.5×/1× 死区带钉死（PeTar "restart 磨洋工" 机制根因，gdb 实测）

**Mistake**: PeTar functional bse 的 10⁶ 步磨洋工（每 Hermite 调用 ~7e5–1e6 AR 子步,反复撞 `ar-max-nstep`）最初被归因为"restart 后 slowdown 冷态"。用 `petar.hard.debug` 位级重放 `data.dump_large_step_ar_n3_*`（一步 5.54e6 子步复现）+ gdb 在 `integrateToTime` 内 ignore-count 采样平衡态、两次采样差分能量分量后定案：slowdown 完全正常（κ_in≈1.3e6 顶在 timescale 帽上,SD_org~1e12）,塌的是 **ds**（自然值 8.3e-4→平衡值 1.3e-8）。

**Root cause**: 配置=超紧双星(P≈1.1e-8)+近抛物线慢速第三体(e≈1.001, r≈4e4·a_in 徘徊 ~1 时间单位)的三体组;第三体轨道相位驱动 ts_in 故 κ_in 以相对速率 ~50–100/单位时间持续漂移;sd 能量记账(ekin_sd_/epot_sd_ ∝1/κ, etot_sd_ref_ 追踪更新)每次更新残留 ~2% 增量 → 每步 |H−H_bk| = 漂移率×dt_step,**随 ds 线性缩放**(实测斜率恒定 ~0.6:err 3.3e-6@ds5.4e-6 → 6.5e-9@ds1.3e-8),不是截断型(ds^order)。控制器只会一招(缩 ds),缩到每步误差 <e_err_ar=1e-8 才停(此时 ds 已低自然值 5 个量级);而增长门(:3470)要求误差 <0.5×max、缩减门(:3332)在 1×max 触发——平衡态误差 6.5e-9 落在 (5e-9,1e-8) 死区带内永不回升(全部 19 次警告的每步误差都在带内)。物理能量交换比观测漂移小 ~10¹³ 倍(pert_out/pert_in≈6e-17),漂移纯属记账——与 2026-09-28 correctSlowDownEnergy 修复同族,κ 时变下仍有残差。

**Prevention rule**:
1. ds 控制器的误差信号必须区分速率型与截断型:对系统性 dH/dt 分量先行估计扣除,或按误差随 ds 的缩放指数在线判别;并给 ds 设自然步长下限(按 P_eff),速率型误差不允许把 ds 无限压低。
2. 增长门(0.5×max)与缩减门(1×max)不对齐构成死区带——任何误差≈rate×dt 的系统必被钉死在带内;至少对齐阈值或加逃逸计数。
3. 诊断模板(本次 6 次 gdb 会话提炼):`break symplectic_integrator.h:3231` + `ignore N` 跳到平衡态 → 打印 ds/err/κ/能量分量 → 两次采样(隔 ~200 步)差分定位携带漂移的分量;比读码推测快一个量级,且直接证伪"κ 冷态"假设(采样即见 κ 顶帽)。
4. 与 09-28 条联动:凡动 κ 的路径(syncTreeSlowDownAndDs/correctSlowDownEnergy/κ cap 变化)都要做"κ 时变注入多大 H 漂移"的账——本次表明残差足以主导误差信号。
5. **修复(2026-10-19 已实现,`integrateToTime` 速率型漂移分离)**:每步算完 H 后维护符号化漂移率 EMA `h_drift_rate_ema_`(前 10 个样本跑步均值、之后 α=0.1,每完成步更新),误差信号改为 `integration_error_eff_abs = max(0, |H−H_bk| − 1.5·|R̂|·dt)`(≥3 个样本后生效),缩减门与增长门都改用 eff(增长门内 ratio 用 eff 并加零除守卫);真尖峰(近心点)误差 ≫ rate·dt,行为不变;警告输出新增 `|Int_err_eff/E|` 与 `|H_drift_rate|` 字段。验证:确定性 dump 重放 5,540,018 → 16,393 子步(338×,0 警告,dE 同量级且略优);第二实现(34.87 双星)dump 整簇 183 子步;干净运行同 seed 重跑轨迹有微差(自适应路径分岔,预期)但守恒账目更干净(旧版一次 dE_SD/Etot_SD≈202 爆发,新版 ≤1.2e-4);PeTar validation T1–T4 全过。剩余开放:κ 时变记账残差本身(~2%/更新)未除,见 plan。**【2026-10-21 更新:本条的补偿架构已被根治方案取代,EMA 机制整体删除】**两个上游根因修复后,EMA 成员/更新块/eff 信号/打印列(`|Int_err_eff/E|`、`|H_drift_rate|`)全部从代码中移除:(a) SD 参考的扰动功记账改为 κ 加权(`kickEtotAndGTDriftTreeIter` 用 `vel_sd`,与 S 的梯度严格一致;逐步 A/B 实测速率信号中位数 1.06e-4→6.6e-8,~1600×);(b) 控制判据改用规范不变的每步能量误差 `control_error_abs = |ΔR|/max(|ref|,|epot|)`——κ 重标只重整 H 规范(实测 72% 步同号 0.97 的"速率"即此伪影),ΔR 按构造免疫。旧判据离开 EMA 会真实磨死(rp=0.35 重现 900s 超时);新判据无 EMA 下 4/4 完成、0 rebase。本条保留作为速率型诊断方法论与"补偿可被根治取代"的完整案例;后续演化见 docs/etot_sd_ref_corruption_fix_notes.md §Fix(授予层修订 + NaN 探测器删除)。
6. **κ 漂移机制与放大器边界(同日深挖,2026-10-20 完结)**:κ_in 顶在**漂移的相遇时标** ts_in(入侵体+场星相遇,簇内 0.0142<ts_max),dκ/κ/dt≈106/单位;E_sd 漂移 = 记账重标度恒等式 dE_sd/dt = −(dκ/κ/dt)·E_sd(数值精确:106×0.93≈98 vs 实测 94)。**最终根因(两环节闭合)**:(a) 帽身份——两侧上下文 inner stab=0、同为 t_min 帽、同一代码路径每步更新(硬件观察点回溯证实;standalone 瞬态期 κ 每步 +3.6% 漂移更快),无帽切换;(b) 放大器 = **perturber 表本身**:H 误差估计量只含内部能量、不含扰动势项,外部 perturber 与群的真实能量流 W 每步以 |ΔH|≈W·dt/|E| 出现——速率型误差,缩 ds 无解;无 perturber 时 W=0、误差为纯截断 ~2e-11;standalone 复现:同 IC 加大 `--r-neighbor-over-group`(10/30)填充邻居表后误差跳至 3.5e-8 且 Large_energy_error 出现——**此前 80 粒子仿真"干净"是因默认邻居半径下组邻居表为空(无效对照)**。量级全部自洽:簇内 W≈0.2-14/单位 ↔ 形成期 err 3.6e-5(dt 3.7e-4)/平衡期 6.5e-9(dt 5.6e-9)。诊断教训:(i) 对照实验必须核对**生效条件**(邻居表内容、调用窗口期——默认 dt-max 会跳过瞬态,早前 kdkpert/独立测试均因此无效);(ii) "平均做功太小"不构成排除——控制器对**每步增量**反应,不是均值;(iii) EMA 修复扣除的正是缺失的开系统能量流,属原理性修正;备选方案(守恒量中纳入 epert)会全局改变误差语义。

7. **自制 IC 验证闭环(2026-10-20,HardDump 克隆法)**:克隆 dump 二进制布局(PtclH4 步长 240B,mass@0/pos@8/vel@32/bin_state@56/r_in@224/r_out@232)替换粒子记录做"种子×激活子"分解:真三体单独=57 步干净;**真三体+16 颗自制合成场星 = 旧版 16,062,968 步(3 警告)vs 修复版 1,245 步(12,900×,0 警告)**——自制 IC 成功触发速率型事件且修复确证有效;合成三体(同标称几何,近抛物线)+真实场 = 764,891 步新旧**完全相同**(合法近心点成本,修复正确地不动它)。结论:种子=真实三体的刀尖态(入侵体勉强双曲、在 ~4e4·a_in 处以 v≈v_esc 徘徊——混沌俘获的产物,标称几何的合成三体落不进这个角,故 35 组合成配置零触发);激活子=任意足够的外场。**放大器=changeover 对扰动力的 r 依赖加权(中和实验定案)**:同种子同场、仅把成员 r_in/r_out 放大到 1.0/1.2 使全部场星处于 k=1 满权区(力更强、仅消除 r 调制)→ grind 消失(16,062,968→650 步,24,712×)——速率型误差地板由 `calcAcc0WTwo` 加权扰动力的逐步调制注入;standalone 样例交互层无 changeover,故**无法完整复现**(同种子+场只到 err≈2e-7 尖峰:AR_KDK_PERT 贡献 6×,TTL 在样例中无效),距 PeTar 级连续 3.6e-5 差 ~170×。诊断方法论:(i) HardDump 二进制可整段克隆重组(杂交实验比再生成 IC 更快定位载体);(ii) "合法昂贵"与"病理昂贵"必须用新控制器做判别(旧版成本高≠速率型;修复版不降的成本=合法);(iii) petar 完整求解器对紧凑合成系统的 r_out 自动推导会被紧双星塌缩(r_out~1e-6),CLI/par 覆盖又被启动重推导覆盖——PeTar 层的自制实验走 dump 重放而非求解器直跑;(iv) 中和实验要避免退化参数(r_in=r_out 会 abort;r_out>r_search 违反不变量),并确认"中和后力更强而病灶消失"以排除强度混淆。

9. **两类边界异常的定性(2026-10-20,4 体扫描;b 项 2026-10-20 晚研究修订)**:(a)"修复版更差"=混沌彩票非缺陷——修复改变早期 ds 决策→刀尖散射系统蝴蝶分岔→该实现内双星被瓦解(semi<0,e=2.15,树打乱,ds_init≈35 荒谬)→打乱的树击败**任何**控制器(rp=0.4 新旧逐位同 abort;rp=0.2 旧版也 grind);完成运行守恒同量级(3% vs 9%)。(b)"震荡型 ds 地板 abort"(rp=0.35-0.5、5.0,新旧逐位一致)=**预存控制器策略漏洞(全量写入归因实测定案,2026-10-20 深夜;推翻当晚早些时候的"de_sd 记账失准"假说)**:相遇把内双星抽到 e=0.9955(近心距 1.9e-9)→近心俯冲期间 κ(~5e4)剧变→**先说被推翻的**:条件观察点抓到的 `etot_sd_ref_ += de_sd`(`symplectic_integrator.h:1697`,−946→+988 过零写)曾被当作"腐化写入"——全量 701 次写入归因(bp+逐次打印)证明该写**良性**:所有写入皆 force-sync,`de_sd ≡ S_post−S_pre`(S=ekin_sd+epot_sd)精确到 1e-9 相对,每次写入**保持**残差 R=S−ref 不变;过零时 R **早已**是 −20407(kick 写入 Σ=−5.7e-4 亦排除)。**真实机制**:`etot_sd_ref≈S` 不是记账不变量而是积分误差规范(getDESlowDown/dE_cum);残差由**已提交的单步**注入——致命步 step_count=524(t=0.00128094584,dt=7.1e-12)一步 dR=−2.75e4(|S|的 0.4%),仅恢复 +4.9e3 后 R 冻结 −2.05e4。控制器没抓住的原因:`fix_step_option=later`(`information.h:516`,`bin_root.stab<1.0` 分支,d45b383 2021-05 引入)在组成型即设上(424 次重申)并**同时关闭两条误差缩减路径**(仅 `none` 缩减;`always` 按设计跳过;`later` 静默穿透)——近心期每步 |H−H_bk|≈1e-6..1e-5,超 ar-max-error(1e-8)百到千倍却全部提交;stab=1.35e-4 判"深稳"是 `stableCheckIter` 的 apo/peri 盲区(`peri_out=semi*(ecc+1)` 名为 peri 实算 apo,e_out 仅进 ^(1/6) 项,ecc=0.9995 的根近心 8.9e-6=22·a_in 完全不可见);ds 重估也救不了:静默门(×3 可复现性)在候选翻腾时全部 DEFER,成本阀需 g_func_on(此处 ar-g-func=0),安全阀因近心 gt 爆炸 n_live 巨大永不触发——ds 从成型到 abort 冻死在 0.0087518。**迟爆 NaN 链**:俯冲后 κ[0] 回帽(→49575),|epot_sd|=|E_in|/κ 随之缩小,t=0.0012895 时 20379 < |R|=20400 → num=ekin_sd−ref<0 → ln(negative)=NaN → NaN 分支减半 ds(误差已提交进历史,结构性)→地板 abort;log 形式还压缩信号(0.4% 注入只显 |dH|~1e-5)。新旧逐位一致另有解释:两版在 `later` 下都从不缩减。修复方向(plan 已更新):(i) `n>2` 分支给 `later` 加双曲/高 e 守卫(c3ce7b0 已给 2 体加过,n≥3 漏掉);(ii) 结构性 NaN 防御:有限分量下 num≤0 → 重置参考+记账+告警,不进减半循环;(iii) stab 的 apo/peri 缺陷单列决策(动它=全局改组管理语义);(iv) **与 2026-10-19 速率修复的互作用须在重现器上验证**——平滑俯冲的 dH 系统同号、dt 塌缩,可能被 EMA 当"速率"扣除(其"真尖峰≫rate·dt"前提对平滑欠解析失效)。判别手法沉淀:**新旧逐位一致的失败=预存问题(修复未激活);分岔后单侧失败=先查轨迹是否已物理分家(树结构/成员/守恒),再归因控制器**;NaN 定位三板斧:abort 帧成员打印(log 参数逐个验符号)→ 条件观察点(`watch var if var>0` 直达变号写入)→ 写入行回溯;**新增第四斧(方法修正)**:变号观察点抓到的是"过零穿越"不是"腐化写入"——残差先于变号时就已存在;归因必须**全量打印所有写入者**(含 R 快照)用相邻写入差分做预算分配(kick 份额=Δref−Σde_sd),或直接问"R 在哪个区间增长";另:`etot_sd_ref≈S` 违约≠记账 bug——它就是能量误差表,先查控制器门控(`later`)再查记账。

8. **最小系统定位(2026-10-20,单扰动星扫描)**:真三体+**单颗**静态扰动星(m=1,+x 方向,其余全部旁置),成本随 r_p 呈干净钟形:0.03/0.05(k=1 深区)=774/546;0.15/0.3/0.5/0.9(过渡带)=1.3e4/5.7e5/**2.07e6**(峰值在带中央)/4.5e5;1.5(k=0 外区)=53(=三体单独基线);m=4 与 m=1 成本几乎相同(与强度无关,纯调制);同位置中和 k(r)→43(基线)。修复版在带内 17-167× 治愈并保留合法成本梯度。触发条件最终版:**刀尖态三体种子 × ≥1 颗处于成员 changeover 过渡带(r_in<r_p<r_out)的扰动星**——k(r) 对扰动力的逐步调制通道。方法沉淀:最小系统(减到单变量)+ 分区扫描比多星场更能钉死条件;静态扰动星(窗口内自由落体可忽略)把调制源唯一化为组成员的运动。

### 2026-10-21: 防御性探测器的存在性用裕度扫描裁决——判据更换后被完全覆盖的旧防御即删;丢弃记账残差的修复必须物理量门控

1. **Mistake**: 控制判据更换(`control_error_abs`)后保留旧的 H-NaN 探测器,并在存在性争议中先构造注入实验证明"分支能触发"来辩护其保留——经过分类→门控→LogH-only→删除四轮评审才识别出冗余。注入实验只证明分支可工作,不证明事件自然可达,两者被混用为存在性证据。
2. **Root cause**: 探测器必要性随机制演化漂移。新控制器接管真实能量误差后,NaN 探测的全部用例被两条既有路径完全覆盖:controller(ΔR 判据,有限量,真实误差 transient 携带 |ΔR|=O(|epot_sd|) 必触发,restore 同时从完好 backup 治愈参考——注入实测 control_error_abs=1.68、ref 逐位恢复)与负 dt 行程(LogH 的 dt 符号=num 符号,域越界自带检测,同事后点同药方;committed 违规有界磨步 ~1e3 次到 ds 地板修复)。"能构造出来"≠"自然可达":裕度扫描(28,217 步逐步记录 num/|epot_sd|/R)实测健康运行离翻负 ≥4.6e6 倍,且 κ 封顶使分母塌缩有界。
3. **Prevention rules**: (a) 判据/机制更换后必须审计所有旧防御的覆盖关系,被完全覆盖即删,不以"保险"名义保留;(b) 存在性争议用裕度扫描量化自然可达距离,注入实验只用于验证分支行为,不得作为保留依据;(c) 丢弃记账残差类修复(rebase/重置参考)必须以物理能量误差门控(超阈 abort+dump,如 abortSdGaugeRepair),规范量(NaN/log 域)只做诊断打印;(d) 数值参考的再校准点应与 gt 重锚、时钟再锚定同族设计,由健康谓词触发而非定时。
4. 验证方法论沉淀:确定性 A/B 用日志逐行 diff 证明行为不变(比二进制 cmp 更抗构建元数据噪声);gdb 两级断点(冷锚点武装热断点+计数自退)在 -O3 -g 二进制上可做受控状态注入,每步打印的规范裕度扫描是"过度防御"争议的裁决手段。

### 2026-10-22: `Neighbor::r_min_mass` 构造与两处重置均缺初始化——NaN 堆垃圾使近邻选择永久中毒(GH issue #2,已修复)

1. **Mistake**: `H4::Neighbor` 的 `r_min_mass` 在构造函数初始化列表、`clearNoFreeMemNoResizeNeighborAdress()`、`resetNeighbor()` 三处全部缺失,而 `checkAndAddNeighborSingle/Group` 首次比较即读取(`_r2*r_min_mass < r_min_sq*mass`)。外部报告(GH issue #2)后用 placement-new 到受控垃圾内存的微测实证分档:垃圾为 0/有限值时无害(`r_min_sq` 每轮重置为 DBL_MAX 占绝对主导,首个正质量候选必被接受,argmin(r²/m) 正确;`neighbors.reserveMem` 只在初始化 `new[]` 一次,构造函数跑过但该成员留堆值);垃圾为 NaN 时比较永假 → `r_min_index` 永为 -1,且 `resetNeighbor()` 不清它使 NaN 跨轮存活——该粒子永远无法经近邻路径成团,静默、依赖堆布局(不可复现性来源)。下游 `checkNewGroup` 的 `j<0 continue` 保证不崩溃,纯正确性/可复现性风险。issue 第二项(`merge_mask` 按活跃组数定长)已在 ed6eb16 修复。
2. **Root cause**: 类成员部分初始化——比较用的质量字段加入时未同步全部构造/重置点;NaN 与一切比较为假的语义绕过了"首候选种子+DBL_MAX 哨兵"的自愈机制。稳态(stale 有限值)本无害的假象掩盖了首轮窗口。
3. **Prevention rules**: (a) 新增类成员时 grep 该类全部构造/clear/reset 函数同步初始化,比较表达式左端读到的成员必须有与输入无关的确定初值;(b) "未初始化读是否真问题"用分档微测定裁(0/有限/NaN 三档 × 修复前后),不以"通常堆是零"辩护;(c) 修复取 `r_min_mass(0.0)`(三处),0.0 恰是占主导的既有良性语义——零质量候选本就不会胜出(`_r2*0 < r_min_sq*0` 为假)且 `checkNewGroup` 另有 mass==0 跳过,故数值行为严格不变。

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

### 2026-09-28: NaN 能量误差断言把可自愈瞬态当致命错误——双曲俯冲 + 过期 ds 一步打出 15% 能量误差

**Mistake**: PeTar 开启 HARD_DEBUG 后 functional `std` 用例（双并合 IC，-u 1，-b 4）死于 `ASSERT(!ISNAN(integration_error_rel_abs))`（第 1 个区间内）。gdb 定位：组是 BSE 踢后的**双曲近距飞越**（etot>0，semi<0），`info.ds=0.0331` 是上个宽轨道时代的值、对该俯冲过大 ~40×：单个组合步穿越近心点，TTL gauge 反复变号（time_table 正负振荡 ±1.5e-8），单步能量误差 15% → `ekin_(1.98e7) < etot_ref_(2.14e7)` → LogH 型 `H = calcH(ekin-etot_ref, epot) = log(负) = NaN`。断言位于控制器回滚分支**之前**，且 NaN 参与的比较恒为 false——控制器本可 restore+减半自愈（负 dt 分支正是干这个的），却被断言抢先 abort。

**Root cause**: H 用 LogH 形式评估 gauge 域（要求 ekin−etot_ref>0），但单步大误差可把它推出定义域；"H 为 NaN"≠"状态损坏"——备份完好、误差是单步截断误差。另记两条 gdb/单位教训：(a) 函数**入口行**的条件断点在序言中评估 `this`，条件永假——本次四个 watchpoint 实验因此全部无效，必须在序言后（如 `dt_full` 行）设断点；(b) PeTar `-u 1` 是 galactic 单位 **G=4.498e-3** 而非 G=1——推 epot 数值/符号前先确认 G，本次因单位误判绕了远路（`epot=-1.6e6 @ r=2.8e-7` 在 G=4.5e-3 下完全正确）。

**Prevention rule**: NaN 能量误差改为恢复分支：restore 备份 + ds 减半重试（与负 dt 分支同构），终止仍交给 ds 地板/下溢与 step_count_max abort；真损坏（持续 NaN）会几何级减 ds 触发兜底。修复验证：该用例完整跑通 t=10、累计相对能量误差 4.4e-9、零断言；ar.logh.ttl.sd.cm 与 hermite 冒烟回归逐位不变。判据设计规则：**"检测量出定义域"类断言（NaN/log 域）必须先问状态是否可回滚**——备份完好时应走恢复路径而非 abort。

### 2026-09-28: SN 踢后双曲对的 ds 双重失配——中断路径丢失 ds 重估（重构回归）+ 能量标度 p_eff 低估飞越时标 2.4e5 倍

**Mistake**: 应用户要求补齐"BSE 事件后 ds 更新"。实施中发现两处独立缺陷叠加成步数爆炸（functional std 用例单区间 2.3M/11.3M 步，两个构建相同=既有问题）：(1) 881abf3 重构把中断路径的无条件 ds 重估（非 slowdown 分支保留的 `calcDsAndStepOption`）换成了 dg 门控——但 SN 速度踢不改变瞬时势能（dg≈0），门不触发，ds 留在事件前的宽轨道时代；(2) 双曲 p_eff 用能量标度 `2π√(|a|³/GM)`——高能掠射在 r≫|a| 处观测时低估穿越时标 ~(r/|a|)^1.5（实测 2.4e5 倍），osculating q 对近径向速度数值敏感（拟合 q=5e-12 vs 实际 r_min=2.8e-7）；再叠加 2 粒子组 `fix_step_option=later` 禁止 ds 恢复 → 近心点分辨率 ds（~1e-11）冻结整个区间，出腿以 dt~1e-12 爬行。

**Root cause**: ds 是"epoch 不变量"的设计对双曲一次性飞越不成立——近心点与出入腿的时标差可达 1e6；"两体一开始定步长然后固定"的语义只对椭圆双星正确。中断路径的 ds 重估在 slowdown 重构中回归丢失，说明该重构时未对照非 slowdown 分支的事件语义做等价审计。

**Prevention rule**: (1) 中断路径恢复 ds 重估：每次中断事件做廉价 O(N) 估算（刷新 ds_est_prev），提交只走单向 >2× 缩减门（`ds_candidate<0.5*ds_bk`），小幅变化/增长不动 committed ds——事件频繁也不抖动；(2) 双曲 p_eff 以当前成员间距的穿越时标 `2r/v` 为下限；(3) 双曲根的 2 粒子组 `fix_step_option=none`（椭圆双星保持 later）——允许近心点后恢复/增长 ds。验收：functional std 用例分钟级→1.27s、零步数爆炸警告、能量 4e-9；ar/hermite 椭圆冒烟逐位不变。性能问题先测再改：本例"每次 interrupt 检查"本身零可测开销（改动前后步数逐位相同），瓶颈全部在既有双曲失配。

### 2026-09-29: hermite 单机驱动 0 组初始化崩溃 + 输出时间分段重置——回文测试两坑

**Mistake**: 搭建组进出对称性回文测试时连环两个坑：(1) 输入组配置行写 `0`（无初始组）→ `COMM::List<int>::reserveMem(0)` 断言崩溃（`sample/Hermite/hermite.cxx` 对 `n_group_init==0` 无守卫）；(2) 前向日志末行 t=0.047 而积分到 T=0.29——一度误判运行提前中止；实际输出含**多段**，段边界处 `Time` 列重置（相邻行 dt 为负可作指纹），物理时间 = `Time + Time_offset`。

**Root cause**: 驱动器初始化假定至少一个组（sample 输入全带组）；输出分段由结构变化/重启逻辑触发，`Time` 是段内时间。另：`--print-width` 过小（20）时负数吃掉分隔符导致数值粘连不可解析（28 安全）。

**Prevention rule**: hermite 单机 IC 至少带一个组（`1 0 2 0 1` 格式：n_group、offsets、member 序号）；解析输出一律用 `Time+Time_offset`，并以负 dt 行识别段边界；快照解析用表头名→索引（`id<0` 滤除组 CM 块，`r_neighbor` 后列名 dt/SD 区分单粒子/成员）；`--print-width ≥ 28`。0 组守卫修复已在工作树（未提交）。

### 2026-09-29: 回文归因三明治——核心积分器时间对称、重启丢慢化状态、组事件是不对称之源

**Mistake**: 首轮回文测试（文件重启反转）测得 5–30% 回程差后一度归因于组判据不对称；分段重启对照显示**重启伪影本身（0.17–0.42）比回程差（0.05–0.10）大 2–6×**——文件重启只保留 (m,x,v)，慢化状态/步长表全部重导。

**Root cause**: 分层归因（新工具 `hermite --reverse-at`，输出边界内存内速度反转，零重启）：(a) **纯双星（零组事件）回程 5.7e-12 / 7.5e-8**——Hermite+AR+slowdown 核心近机器时间对称；(b) 任何含组进出的三体构型回程速度误差 30–50%（差 6–8 个量级），且事件序列非镜像（f/b=2/0、4/1、18/1）；(c) 误差幅度受混沌放大调制，与事件数无单调关系（rg=0.012 事件最少误差最大）。形成（距离判据）与释放（r_break+轨道根数）判据不配对是结构性的。

**Prevention rule**: 测时间对称性必须用内存内反转（--reverse-at），文件重启反转混入慢化状态丢失；对称化改进的靶点是组进出判据（单阈值状态函数 + 力平滑混合），核心积分器与重启机制是两个独立问题（后者值得单独修：checkpoint 补慢化状态）。

### 2026-09-30: 对称判据落地受阻于更深层——嵌套 slowdown 树自身不时间反演（回文分层测量的新发现）

**Mistake**: P2 对称判据（`--group-symmetric`：形成去 drdv<0 方向门、释放纯状态函数 r>R&&apo>R）实现后回文验收仍无改善（dpos 与 legacy 逐位相同），一度以为判据无效果；全组构型（rg=0.1，无单粒子无组事件）测试揭示真因。

**Root cause**: 分层回文测量：(a) 单层慢化双星（κ=32）回程 5.7e-12 ✓；(b) **嵌套慢化树三体（内对+外轨双层）全组无事件回程 mirror=9.6e-5、dvel_rt=1.3**——慢化树的内部状态（更新调度/子轨道相位）在速度反转后不对齐，位置近似回而内轨相位错（dvel 放大）；(c) 判据层在此之下被掩盖。另：驱动级 jerk（force_ 数组 acc1）在下次块更新前不被预测使用——翻转与否无差（假设被数据否定）。对称判据本身已生效（同轨迹事件 2→29，apo>R 对每近心一对进出）但当前过抖，需 β 平滑窗。

**Prevention rule**: 回文归因必须逐层：纯双星→全组多体→混合单粒子+组，每层一个反转语义问题。嵌套慢化反演语义（syncTreeSlowDownAndDs 调度的奇分量/慢相位重启）是判据对称化之前的先决修复。`reverseVelocities()`（vel+force.acc1+vcm_record+pred 副本）已加入 HermiteIntegrator 作为标准反转原语。

### 2026-09-30: 修正前条归因——组"事件注入"确证；判据对称(何时)≠过渡对称(如何)；回文注入地板 ~dt_max²

**Mistake**: 前条把全组三体回文缺口归因"嵌套 slowdown 树不可逆"——被数据否定：该构型与弱扰动 IC 的 κ 全程=1（`--slowdown-ref 1e-30` 与默认逐位相同，慢化从未激活）。真正因果链由同 IC 三配置对照闭合。

**Root cause**: (1) **回文注入地板**：反转点一次性局部误差（tsyn 着陆/子步非自伴），随 dt_max **平方**衰减（0.25→1e-2；0.0039→6e-6）——回文验收必须用细块；(2) **组事件注入**：同弱相遇动力学下，无事件=地板 5.6e-6，2 次事件=1.5e-3（270×）——每次组进出注入 O(1e-3) 不可逆误差，与判据模式无关（对称判据 7-11 次事件同幅度）；(3) 判据对称（事件**何时**触发已镜像）不充分——**过渡对称**（事件**如何**改写积分器状态）才是注入源：merge（addGroups：vcm_record/ds/epert 从新树重导）与 break（breakGroups）两条初始化路径不对称，AR 子步相位在树重构时重置。混沌放大只在地板×弱增长下起次要作用（无事件对照排除主因）。

**Prevention rule**: 回文验收配置=弱扰动 IC + κ=1 + dt-max-power 8（地板 6e-6）；判据层与过渡层分开验收：前者看事件时间镜像（对称判据已达成），后者需统一 merge/break 的状态重导路径（同一 state→state 函数，方向无关）——这是下一个实现目标。归因结论必须配 κ 与 SD_factor 列的实测（两次"slowdown 假设"都被 κ=1 否定）。

### 2026-10-01: 统一判据的 κ 分支在 break 侧是死代码——守卫条件与短路首分支矛盾；κ 估计约定两边不一致

**Mistake**: 统一判据重构后声称 "form=P / break=¬P 构造性互补"，但 break 侧保留了守卫 `(bin_root.r > r_break_crit) && n_member==2` 才算 κ，而 `groupedCriterion` 首分支 `d > r_crit → false` 直接短路——κ 只在 d>r_crit 时被计算，恰是它永不被咨询的区域；d≤r_crit 时 κ 取默认 1e-2（≥阈值自动通过）。break 实际退化为纯 `d > r_crit`，在宽 apo 区与 formation 不互补。同时 κ 估计约定两边不同：formation 喂力单位质量加权场（含内部互势项），break 喂 CM 比场；formation 的 κ_ref 带质量比缩放，break 用组内存储值（初始化时不带缩放，`symplectic_integrator.h:1874`）——`force²/pot` 不守恒于该缩放，1e-2 阈值两边含义不同（差 ~m_cm 倍）。

**Root cause**: 把 "κ 只在需要时才算" 的优化写成与谓词短路结构冲突的守卫；两处调用点各自内联 κ 配方导致约定漂移无编译期/测试期可见性（样本单位 G=1、m~1 下因子≈1，数值上不可见，PeTar 单位下 1-2 个量级）。用户质疑"判据真的对称吗"触发审计才发现。

**Prevention rule**: 统一谓词必须**吞掉**其全部输入的计算（估计器内聚进 `groupedCriterion`），调用点只传原始态量并统一到单一物理约定（这里是 CM 比场、external-only）；"互补性" 声明前必须逐分支核对可达性（守卫 × 短路路径），不能只看公式形态。验证：定向 IC（宽 apo e=0.9 双星 + m=20 掠过）触发二体 κ-break（d≪r_crit、κ≈4e-3）；回文不劣化（事件计数 3/1→4/3）。

### 2026-10-03: hermite 单机 checkpoint 用活跃组数直接索引组槽位——mask 复用后写错组并 abort

**Mistake**: `writeCheckpoint`（hermite.cxx）以 `getNGroup()`（活跃组数）为界直接索引 `groups[i]`；组槽位经 break→mask→复用后，前 n_g 个槽常是已 clear 的空槽：`members.reserveMem(0)` 断言 abort（Plummer N=20、rg=0.05，t≈246 首触发），即使不 abort 也写出错误的组配置（重启即静默错）。长事件 bench 被该 bug 阻断，一度疑为新插桩所致（plain 二进制同参复跑排除）。

**Root cause**: 活跃组身份在 `index_dt_sorted_group_`（排序索引表），`groups` 是含 mask 槽的存储数组——两个"第 i 个"含义不同；类无访问器，调用点自行猜测索引约定。

**Prevention rule**: `HermiteIntegrator::getGroupIndexSorted(i)` 已加；任何遍历活跃组的输出/后处理必须经它取槽位，禁止 `groups[i]`（i<getNGroup()）。诊断归属未知崩溃时先用未改动的二进制同参复跑（本次 40 s 排除插桩嫌疑）。

### 2026-10-03: 回文对照三坑——对比对象是初始态、slowdown-timescale-max 默认耦合 time-end、事件日志须用物理时间

**Mistake**: Phase-1 回文首轮全部 O(1) 无效：(1) pal(2T) 与 ref(T) 的 checkpoint 对比——正确对象是**初始态的镜像**；(2) ref(-t T) 与 pal(-t 2T) 未显式钉住 `--slowdown-timescale-max`（默认=time-end，两跑动力学不同）；(3) GEID 事件日志打 `time_`（内部时钟，shiftTimeOrigin 重锚后周期回绕）→"事件集中在初始相"误读 → 窗口选错；另 N=20 窗口 >数个时间单位即混沌饱和。

**Root cause**: 回文数学对象是 state(2T)≅镜像 state(0)，中间任何状态都不可作参照；time-end 身兼积分终点与 slowdown 上限默认两个语义；time_≠getTime()（物理=内部+offset）。

**Prevention rule**: 回文对照= pal(2T) checkpoint vs 质心平移后的 IC（`pal_vs_init.py` 模式，H4Particle 记录 184 B：Particle 112 + dt/time/acc0/acc1/pot 72）；成对跑必须显式同传 `--slowdown-timescale-max`；诊断日志一律打 `getTime()`；窗口选在事件簇后最近输出格点且总长 ≲2 时间单位。基准数：安静地板 5.1e-6（dt-max-power 8，与 2026-09-30 的 5.6e-6 一致）；+1 对 chatter 事件 5.3e-5（10×）；同相遇 Hermite-only（rg=1e-7）O(1)——混合方案比纯单粒子路径可逆性好 4 个量级。

### 2026-10-03: 切换位置主导缝合成本——早切换在孤立相遇好 400×、强扰动场仅 ~3×（回文 r-group 扫描）

**Mistake**: 归因链条上把"切换固有截断"当作不可调的常数残差收尾；用户质疑"是否切换时机不对"触发参数扫描，发现切换位置才是主导杠杆。

**Root cause**: 专用三体回文 IC（0.45+0.45 双曲相遇 peri=0.02 e=1.067 + 伴星）扫 `--r-group`：0.05（近心才切，55-60 AR 子步）dvel 3.0e-4；0.4（早切，460 子步，整个强场段含近心点由 AR 一体积分，两次切换都在 dr≈0.44 平滑场区，且反转点落在组存活期内不再触发新切换）dvel **7.5e-7**（400×）。伴星拉近加重扰动（pert 比进入生产量级）后：0.05→6.9e-4，0.4→2.4e-4，优势缩到 **3×**——AR 组内存活越长，冻结扰动误差累积越多，吃掉缝合节省。

**Prevention rule**: 组进出的误差评估必须报切换半径，不能只当固定残差；最优 r_group 是扰动场的函数（孤立系统宜大、密集场收益有限且伴随更多组/更长寿命/性能代价）。生产侧验证**不能**用回放：判据半径烧录在 dump 的逐粒子 changeover 态，data.par.hard 的 r-group 修改无效（实测 0.1–10× 因子下 r_crit 逐位相同，x30 触发 r_search 断言）——须用改参的目标生产运行段（或 dump 补丁），也别用孤立 IC 外推。另：r-group 经 r_neighbor→`calcAcc0OffsetSq` 间接影响块步控制器，"无组"对照在不同 r-group 下不是同一条轨迹。

### 2026-10-03: 【切换层注入归因三律】一次性记账看净值不看单项、审计表会过期、对照组必须同相遇——切换层终判"是解的一部分",成组门优化方向已死

**Mistake**: 切换层统一调查(Phase 1/2,2026-10-03 关闭)中途的三次归因偏差:① Phase-1 bench 报"每 break 3e-2 注入"并据此锁定 #4/#5 为待修项;② 审计表(§3)把 #4 kappa 采样调度、#7 重入初始化分裂列为待修不对称;③ pre-consolidation 时期的"无事件 5.6e-6 vs 2 事件 1.5e-3(270× 劣化)"测量被当作切换层有害的证据。

**Root cause**: ① FORM 时 AR 记 `de_sd=+E_bin(1-1/κ)`、BREAK 时 Hermite 记 `-E_bin(1-1/κ)`——**逐对相消**,单事件单项是定义性缩放项,每周期净注入 ~1e-3·E_bin 才是预算量;② 审计表描述的是 pre-consolidation 行为,#4 的"formation 锚定调度"(`time_update_`/`setUpdateTime`)在 HEAD 已是死代码,#7 重入路径已统一为 `initialIntegration` 同款;③ 270× 来自未配对对照——正确对象是**同一条相遇轨迹**的分组 vs 不分组(matched-control:分组 3.0e-4 vs 纯 Hermite 9.5e-3,**好 32×**)。

**Prevention rule**: ① 注入/误差预算先验证记账项是否成对相消,报净值;② 按审计表设计修复前,先在 HEAD 上确认审计项仍存在(grep 调用点/测死代码);③ 跨方案误差比较必须 matched-control(同 IC 同相遇,仅切换 r-group 门),非配对的历史测量不可复用;④ 附带裁定:fixed-ds(`-e 1e10`)与自适应位相同→自适应 ds 控制器(#6)对 chatter 类无贡献;"值不值得成组"门方向已死(单扰动源 68×..1.1e8×、12 伴星环 6.4×..7000× 分组全胜,无有害阈值);kappa gate 与 kappa_ref 约定强耦合(sample 默认 1e-6 在潮汐环境完全拒绝成组,生产用 1e-4);判据/切换工作的回归工具是 `sample/test/group_transition_bench.sh`(8 案例+`--scan` 杠杆方向门)。

### 2026-10-13: 【退化分支门控】crossing dt 界限只应在准则退化(零加速度/零 jerk)时介入——无条件 cap 干扰 Aarseth 准则,且门控版顺带修复 T1 最粗端点悬崖

**Mistake**: per-slot crossing 界限(459e-461e)在所有调用点无条件 `min()` 进 dt_limit,对一切有速度的粒子生效;用户指出这干扰了原算法——有真实力时 Aarseth 准则本身随加权力的增长自适应收缩,几何 cap 冗余且在弱力非零区产生不必要的过小步长。5765150 提交信息自证其职责本就是"force-free particle"场景(旧 HERMITE_MONOPOLE_DT_LIMIT 的替代)。

**Fix**: 把界限移入 calcBlockDt2nd/4th 内部:仅当 `s0==acc0_offset_sq || s1==0.0`(准则无力度尺度)时,dt_ref 用 `max(eta_4th·r_limit/|v|, dt_min)` 替代 NUMERIC_FLOAT_MAX;调用点恢复原算法签名(vel/r_limit 作模板参数传入,兼容 PS::F64vec)。力驱动的 sub-dt_min 需求仍按策略 abort。验证:note 复现器干净(−3e-12);T3 10/10 计数器逐位不变;T1 si00 与旧版逐位相同;**T1 crossfull_r015_si07(最粗 dt 端点,458e-461e 在 LWLaptop 全崩)两族均确定性通过**——机制:粗 dt 交界后的伪轨道在关键步出现精确零加权力,458e 无限制→dt_max 跳步落上碰撞航线,门控版以 eta·r_in/|v| 尺度解析接近段。

**Prevention rule**: ① 几何类 dt 界限只应填补充准则失效的退化分支,不得无条件叠加——先问"准则在该状态有没有力度尺度";② 安装二进制按 configure 族覆盖:一次 `make install` 只装当前族,跨族验证前必须逐族重装(本次 T1 non64b 首轮用了旧安装件假失败);③ 端点悬崖的 FP 敏感结论要保留"策略改动可移动悬崖位置"的开放性——语义修复需单变量对照(同机同编译器只改源码)确认。
