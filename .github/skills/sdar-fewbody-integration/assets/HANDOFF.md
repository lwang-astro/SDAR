# HANDOFF: BTLogH → Hermite/PeTar 集成（2026-09-14，SDAR experiment@444e）

> 交接对象：PeTar 端 agent（Developer / Implementer）。处理本任务前先读本文件与
> `lessons-learned.md`（2026-09-14 条目）、`data-readback-patterns.md`。

## 2026-09-14 执行记录（代码落地 + 关键阻塞发现）

**已完成（接线全部落地并验证）：**

1. SDAR 样例端：
   - `sample/Hermite/hermite.cxx`：`--g-func` CLI（默认 0=LogH；1=BTLogH；
     BTLogH 构建拒绝 2），`ar_manager.g_func` 传递。
   - `src/AR/symplectic_integrator.h`：`TimeTransformedSymplecticManager` 增加
     `g_func` 成员（`#ifdef AR_G_FUNC`，默认 0；注意改变 g-func 构建的
     manager 二进制 dump 布局，旧 dump 不兼容）。
   - `src/Hermite/hermite_integrator.h`：`addGroups` 中 `group_new.g_func =
     ar_manager->g_func`（`g_func_on` 由组 `initialIntegration` 解析）。
   - `sample/Hermite/Makefile`：`hermite.btlogh` target（默认 flags +
     `-DAR_G_FUNC_BTLOGH`）。
2. PeTar 端：
   - `src/hard.hpp`：IOParamsHard 新增 `--ar-g-func`（0/1，默认 0）；
     `HardManager::initial` 传递 `ar_manager.g_func`；孤立双星 sym_int 路径
     设 `g_func` 且 `calcDsAndStepOption` 传 `g_func_on`（`#ifdef AR_G_FUNC`）。
   - `configure.ac`/`Makefile.in`：`--with-sdar-g-func=logh|btlogh`
     （btlogh → `MT_FLAGS`/`HARD_MT_FLAGS` 加 `-DAR_G_FUNC_BTLOGH`，
     `PROG_NAME` 追加 `.btlogh` 后缀）；`configure` 已 autoreconf 重生成。
   - 文档：PeTar README（configure 小节）、`option-reference.md`（已用
     `generate_option_reference.py` 重生成，含 `--ar-g-func`）。
3. 验证：
   - SDAR Hermite 四目标编译 0 error。
   - **LogH 惰性逐位验证**：ustabtri H4（plain vs btlogh构建 --g-func 0），
     364 输出行 × 103 物理列（剔除 8 计时列）NaN-aware 全等。
   - PeTar `hard_test.cxx` 在 `-DAR_G_FUNC_BTLOGH` + 全硬积分器旗标下语法编译
     0 新增 error（唯一 error 为存量 hard_test.cxx 签名过期，基线同样存在）。
   - configure 管线验证（btlogh 配置旗标/后缀注入正确），用户原配置已恢复
     （`--with-interrupt=bseEmp --with-external=galpy`）。

**✅ 阻塞已解除（2026-09-14 晚，双修复落地并验证）：**

1. **Bug #1（ds=inf）**：`applyStableCheckAndSlowDown` 的 factor=1 分支不同步
   `sd_root.period`——2 粒子组（根=叶子）永不进 `calcBinaryTreeSlowDown`
   内循环，`SlowDown::period` 停留构造默认 DBL_MAX；BLogH ds 估计器经
   `getEffectivePeriod()=P·κ` 读到 DBL_MAX → `ds_prod·P_eff_min` 溢出 inf。
   LogH 口径在线算 P 从不读该成员（隐形），3+ 组叶子在内循环已设（历史
   验证全过）。修复：椭圆根无条件同步 period。
2. **Bug #2（gt_drift 失控）**：`integrateTwoOneStep`（pair 路径）内的
   `AR_G_FUNC` 分支**本不应存在**——两体时乘积 g ≡ 单对势 ≡ LogH 求和，
   drift 演化恒等。原 `dgt_drift_inv *= gt_kick_inv_` 平白多乘 U（pair 路径
   消费 interaction 层原始 ∇U gtgrad，÷U 换算只在树路径做）；圆轨道
   v·∇U≡0 完全掩盖、偏心轨道暴露（dE→1e-2..1e15）。修复：删除 pair 路径
   全部 g_func 分支（用户裁定；中间版 `*=κ⁻¹` 在 κ=1 恰为无操作故通过
   两体验证，κ=21.8 quad 仍失败——分支本身即错）。

**验收矩阵（全部通过）**：
| 场景 | 结果 |
|---|---|
| AR n=2 圆/e=0.9，κ=1 | gf1 ≡ gf0 逐位一致 |
| H4 n=2 双星 κ=64000（bin2） | 零错误，dE~1e-13（末位差=ds 公式浮点路径） |
| H4 quad_sd2（2×两体组 κ=21.8）t=20 | 零错误，终点 dE 与 gf0 四位一致（-0.001116） |
| H4 ustabtri gf1 全程 | 零错误，dE -0.0126（gf0 -0.018，同级） |
| 3+ 回归：ustabtri{auto,fixed,logh} s128 + ustabquin s128（59,806 行） | 对 HEAD 逐位一致（仅墙钟列差） |
| 惰性 plain vs gf0（ustabtri H4 364 行） | 0 差异 |

**Step 2/3（H4/PeTar A/B）已解除阻塞。** 检验系统重新设计仍需覆盖
2 粒子组 × 乘积 g × {κ=1, κ≫1}（本轮最脆弱面；两修复后已可用作回归项）。

证据目录：`/tmp/h4_gfunc_smoke/`、`/tmp/btlogh_n2_debug/`（易失，结论已录
本文件、impl notes 与 lessons）。

## 目标

SDAR Hermite (H4) 与 PeTar 硬积分器增加 **BTLogH 选项**（g-func tree-product，
`-DAR_G_FUNC_BTLOGH`），**默认保持 LogH**，BTLogH 显式开关开启。

## SDAR 端当前状态（全部实测，证据可复查）

1. **κ 修复已提交**（LogH sum gauge 椭圆周期去 κ 乘子，VERSION 444e）。
   quad_sd2 h4 全程验证：a1 终值 −1.15e-7，|dE|/|E| ≤ 2.24e-7，
   0 条 `Large_energy_error`，ar_step 8.35e8。
2. **接口完整**：`sample/Hermite` 以 `-DAR_G_FUNC_BTLOGH` 编译通过（0 error）；
   interaction 类**无需改动**（g 从 `info.binarytree` 求值，不要求新 API）。
3. **编译旗标对 LogH 惰性**（已证）：默认 `g_func=0` 时输出与无旗标构建
   逐位一致（120 行 × 126 列，剔除 `prof_*[s]` 墙钟列，NaN-aware 全等）。
   → 任何构建可以先加宏、零风险。
4. **H4 已有脚手架**：`src/Hermite/hermite_integrator.h:2664` 组循环内每步
   `calcDsAndStepOption(..., groups[k].g_func_on)`（`#ifdef AR_G_FUNC` 分支）
   → BTLogH-Adjust 的逐步 ds 自适应在 Hermite 架构下**天然成立**。
5. **κ 回归在 BTLogH 乘积形式结构性不存在**：
   `ds = Π(ds_i)·P_eff_min/Π(P_eff)`，κ 只出现在 P_eff——单层（双星组，
   Hermite 最常见）完全消去（ds = ds₁ = 旧 LogH 叶片规范）；多层时其他层
   κ 只进分母（更保守）。实测 Aug 基线 BTLogH+slowdown AR 长程
   （quad_sd2_1e6_tkz，t=4153）a₁ 采样 1e-8–1e-6。

## SDAR 样例端待做（✅ 2026-09-14 已全部完成，见顶部执行记录）

- ~~`sample/Hermite/hermite.cxx` 加 `--g-func` CLI~~（默认 0=LogH，1=BTLogH；
  auto(2) 在 `AR_G_FUNC_BTLOGH` 下拒绝）。
- ~~manager → 组的 `g_func`/`g_func_on` 传递~~（manager 新增 `g_func` 成员，
  `addGroups` 传递；`g_func_on` 由组 `initialIntegration` 解析）。
- ~~`sample/Hermite/Makefile` 加 `hermite.btlogh` target~~。

## 下一步（优先级序）

1. **修复 2 粒子组乘积 g 路径**（见顶部阻塞节；AR standalone 可独立复现，
   是最小调试载体——`ar.btlogh.ttl.sd.cm --g-func 1 -t 1 bin2 输入`）。
2. 修复后重跑 2 粒子 / 3 粒子 H4 冒烟，再按用户重新设计的检验系统执行
   A/B（原 quad_sd2 全程 A/B 作废）。
3. PeTar 侧构建一个 `.btlogh` 二进制族做小 cluster A/B（默认 LogH）。

## 风险与注意

1. **未验证组合**（中风险）：3+ 粒子组 × Hermite 块同步 + d86b1e2 着陆截断 +
   每块 ds 重估——ustabtri/quin 验证的是 AR `-m` 路径，此组合没跑过；
   且 4a0387b（09-13）刚修 "BTLogH -m full 模式 getH"（最接近 Hermite
   路径的代码年轻）。**切换后必须 A/B**。
2. **性能预期**（低风险）：Hermite/PeTar 的 AR 组多为孤立双星——2 粒子时
   BTLogH 的 g ≡ 单对势（=LogH），只剩 tree-g 开销、无收益；收益集中在
   交会期 3+/4 体组。论文的步数收益**不可直接外推**到 PeTar 生产。
   **默认 LogH、BTLogG 显式开关**。
3. **PeTar 侧接入点**：
   - SDAR 为同级目录引用 `-I../SDAR/src`（PeTar `Makefile:347`，**非
     submodule**）→ PeTar 重编即带上 κ 修复；LogH κ=1 逐位不变；
     slowdown 组行为变化 = 修复生效，**需重立能量误差基线**。
   - 需要：(a) 硬积分器编译单元加 `-DAR_G_FUNC_BTLOGH`（LogH 惰性已证）；
     (b) PeTar `src/ar_interaction.hpp` 在该宏下的编译测试（SDAR 样例通过，
     大概率 OK）；(c) PeTar 参数管线加 g-func 选项（`external_hard.hpp` 一带）；
     (d) `g_func_on` 为组内状态，线程安全模式与现有一致。

## 验证阶梯

1. **Step 1**：SDAR sample plumbing → `hermite.btlogh` 可运行（`--g-func 1`）。
2. **Step 2 A/B**（数据/命令见下）：quad_sd2 h4 BTLogH vs LogH 基线
   （Δa/a、ar_step、wall time、消息数）+ ustabtri H4（三体节点势能路径）。
3. **Step 3**：PeTar 小 cluster A/B（能量误差、wall time），默认 LogH。

## 测试数据与基线（路径均已核实存在）

- 输入：`/home/lwang/localdata/SDAR_BLogH/quad_sd2`（内容自 08-15 未变）
- h4 命令：`hermite -t 4153.283038753799 --r-group 0.1 --r-neighbor-over-group 20
  --dt-max-power 13 -o 3 -e 1e-4 <input>`（sd1e2 变体加 `--slowdown-ref 1e-2`）
- LogH 基线（修复版、全程）：`quad_sd2_h4.log` / `quad_sd2_h4-sd1e2.log`。
  注意 SD2 的 a₁ ±6e-4 游走 = 松 slowdown(ref 1e-2) 固有近似误差，非缺陷。
- 对照二进制：`~/bin/hermite`（修复版，生产）；`hermite.prefix-broken`
  （κ 回归演示）；`hermite.aug`（Aug-1 干净重建）；`~/bin/hermite_old`
  （07-31 原版）；`/tmp/hermite.btlogh`（btlogh 编译/LogH 运行，/tmp 易失，
  用下方命令重建）。
- A/B 沙盒：`/tmp/h4ab/`（含各版本 90s 测试的 out_*.log/.err）。
- 读回：`sdar.HermiteData(N_particle=4, N_sd=2, time_measure=True)` +
  数值行预过滤；a₁ 从 p0/p1 位置速度直接算（精度结论勿依赖消息条数）。
  详见 `data-readback-patterns.md`。

## 复现命令

```bash
# btlogh 编译（已验证 0 error；注意 Makefile 默认还有 -DSDAR_TIME_MEASURE）
cd SDAR/sample/Hermite && g++ -std=c++11 -I../../src -I./ -O2 -Wall \
  -DAR_TTL -DAR_SLOWDOWN_TREE -DAR_SLOWDOWN_TIMESCALE -DSDAR_TIME_MEASURE \
  -DAR_G_FUNC_BTLOGH hermite.cxx -o /tmp/hermite.btlogh

# 惰性验证法：同输入同 flags 跑两构建，比较 .log 剔除 prof_*[s] 列后逐位一致
```

## 论文约束（重要）

BTLogH 论文在审稿修改中（`/home/lwang/write/Astrophysics/BTLogH/paper.tex`）。
无 slowdown（κ=1）的所有论文数字与 κ 修复前**逐位一致**
（ustabtri/ustabquin/quad B--B--S 不受影响）；AR BTLogH+slowdown 的
Aug 基线数字维持不变（乘积路径未动）。集成工作**不得**改变这些数字。
