# LogH/TSI 步长效应(S 判据)与 ds 估计器双曲分支合规性

> **性质**:theory + implementation reference(Wang & Nitadori 2020, MNRAS 497,
> 4384 的 S 判据摘要;SDAR LogH 求和口径双曲分支对该判据的违反推导与 merger
> 案例实测闭合;修复设计依据与实施验证记录(§6);残留恢复注入极限环机制(§7)。
> 原文 PDF:PeTar/doc/Step_size_effect_in_LogH.pdf)。
> 最后更新:2026-10-23。

## 1. S 判据(Wang & Nitadori 2020 摘要)

对时间变换辛积分器(TSI/LogH,DKD 型)定义**无量纲步长**

$$S \equiv \frac{\delta s}{L},\qquad L = m_1 m_2 \sqrt{\frac{G|a|}{m_1+m_2}}$$

(δs 为一个 DKD 步的**半步**;L 为平近点角的共轭动量;注意
$L = $ `calcDsHyperbolic` 取 coff=1 时的 ds_i,即该函数的单位)。

- **椭圆**:$S=\tan(\delta E/2)$,任意实 S 都有解(δE<π)— 步长无上限,轨迹恒精确。
- **双曲**:$S=\tanh(\delta E/2)$,双曲正切值域 (−1,1) ⇒ **必须 S<1**。
  - S>1 时发生**相变**:r 跳到双曲轨迹的对称伪分支(排斥力解),
    守恒量从 $T+U=H$ 变为伪守恒 $T-U=H$,
    $\log(T_1-H_0)$ 中 $T_1-H_0<0$ ⇒ **H 为 NaN**,
    且 $dt/ds<0$ ⇒ **物理时间倒退**(δt<0)。
  - S=1 恰在拐点:T₁−H₀=0,第二 drift 除零,r→∞。
  - 该现象由超新星踢后双曲化的双星发现;BS 外推法中不收敛。

**这就是 merger 案例磨步的机制**:dt<0(96% 事件)、H=NaN、伪守恒 — 全部
是 S>1 相变的教科书症状。

## 2. LogH 求和口径双曲分支违反 S 判据的推导

`calcLogHSumGaugeIter`(g_func=0 路径,merger 案例所用)对双曲层取
u_min = 近心点规范 U(q)=Gm₁m₂/q,q=|a|(e−1);p_eff = 2π√(|a|³/GM)(后被
t_cross=2r/v 下限抬高)。单层贡献:

$$ds = \frac{1}{32}\cdot\frac{Gm_1m_2}{|a|(e{-}1)}\cdot 2\pi\sqrt{\frac{|a|^3}{GM}} = \frac{\pi}{16(e{-}1)}\,L$$

$$\Rightarrow\quad \boxed{S = \frac{ds}{2L} = \frac{\pi}{32\,(e{-}1)}}$$

**e→1(近抛物线)时 S 发散** — 违反 S<1。数值验证(merger 案例,首次拟合
|a|=1e-8, e=1.001):理论 S=98.2;直接算术 98.2(未加 floor)/ 1013(加
t_cross floor,×10.3 抬高 — floor 本意防出腿过分辨,此处反向放大违反);
实测 ds_init=0.288 与公式含 floor 输出 0.304 吻合;第二次拟合 S=246。

**对照**:BLogH 族 `calcDsHyperbolic`(coff=2π/8)单层 ds=(2π/8)L/32
⇒ S=(2π/8)/64≈0.012 — **构造上 S 安全**(1/256 相遇约定)。即:双曲安全
形式只存在于 BLogH 乘积路径,LogH 求和路径从未对齐。

## 3. E2(pre-etot 修复)幸存机制 = NaN 块即运行时 S 执行者

被 `6a8a357` 删除的 `ISNAN(H)` restore+halve 块:merger 案例中 H=NaN 正是
S>1 相变的 log 域症状(log(T₁−H₀<0))。实测该块触发 **11 次**减半:
0.288/2¹¹ = 1.4e-4 ⇒ S=0.469<1 — **恰好停在相变阈值之下**,随后收敛完成。
删除后(理由"H 仅为规范伪影、TTL 控制器接管")该防御消失,相变步被接受,
积分规范被腐蚀,进入不可恢复的负步磨步(864k 事件)。删除论证所用的裕度
扫描(bookfix/gauge 日志)未覆盖高速近正撞双曲类 — 该类中 NaN 是唯一前置
可恢复探测器。

## 4. 修复设计依据(plan item 6,2026-10-22)

1. **估计器(F1,主修)**:LogH 求和口径对每个双曲层施加 S 上限 — 该层对
   ds 的贡献不得超过 ds_Smax = 2L/c_max(c_max 为所用车贝雪夫/台形合成最
   大系数,留裕度取 S_max≈0.5)。实现即对该层 u_min 乘 min(1, ds_Smax/ds_layer)。
   t_cross floor 保留(出腿方向正确),S 上限在近抛物线方向接管。
2. **运行时防御(F2)**:恢复 NaN 瞬态闸门(restore+halve+retry,不带已删除
   的 rebase 部分),触发条件 ISNAN(integration_error_rel_abs) ||
   ISNAN(control_error_abs) — E2 已验证的收敛路径,兜底任何残余 S>1
   (如树打乱后 e≈1 拟合)。
3. 高阶合成(ar-sym-order −6 等)的负系数使子步半步为 c_k·ds — S 判据应
   对最大 |c_k| 取裕度;具体常数在实现时从 manager->step 系数表取。

## 5. 与 plan item 4 的关系(item 4 于 2026-10-23 关闭)

同属"双曲轨道 × ds 机制"家族、不同层:item 6 = 控制器/估计器 S 判据违反
(本文);item 4 = 估计器首估对垃圾双曲拟合无条件接受。LogH 侧由本修复闭合
(S-cap 使任何双曲拟合的首估 S-safe)。原计划的 g-func 侧补丁——给
`calcEffectivePeriod` 双曲分支加 t_cross 下限("镜像 LogH 侧约定")——经
设计审查**否决**(2026-10-23,用户相位无关性质疑触发):

1. **来源案例未隔离验证**:该 floor 出自 2026-09-28 的三缺陷捆绑修复
   (`881abf3` 重构丢失中断路径 ds 重估[主因] + p_eff 能量标度低估
   2.4e5× + 双粒子组 `later` 冻结),验收是捆绑验收;floor 的孤立贡献
   从未测量,且该案例的 osculating q 拟合本身是垃圾(q=5e-12 vs 实际
   r_min=2.8e--7)。承重修复是事件驱动重估+恢复许可,两者已在位。
2. **g-func 乘积结构免疫**:ds = Πds_i·P_eff_min/ΠP_i 中,孤立双曲对的
   周期在分子(P_eff_min)与分母(P_prod)中精确相消——ds = ds_i,纯
   轨道根数,S≈0.012 由 2π/8 相遇约定构造安全;长腿成本病理
   (ds ∝ p_eff 直接相乘)是求和口径特有。生产佐证:N=1000 BTLogH 与
   LogH 步数平价(g3 5076 vs g6 4459/树步)。
3. **相位无关性设计原则**(用户论点,采纳):元素法 ds 对椭圆层是
   原则性的(周期轨道,每轨道一 ds,轨道平均良定义);引入 r_now 使
   估计依赖事件发生时的相位——跨轨迹不可复现,且无已证的 g-func 失败
   案例为该代价背书。双曲层的"epoch 不变量"本无周期性支撑,但两害
   相权:g-func 现状(元素法+事件重估+构造性 S 安全)无 observed 病理,
   不为假想对称性引入相位依赖。

观察(不开条目):LogH 侧 2026-09-28 的 floor 同受相位依赖批评;其有害
方向已被本修复的 S-cap 钳住(merger 案例中 floor 单独把 S 从 98 抬到
1013,S-cap 现封顶 0.25);若要移除需在 std 用例重现器上做 floor on/off
受控 A/B。原"估计器需要统一的双曲合规层"合并结论作废——两族公式结构
不同,合规层不需要统一,各自已闭合。

## 6. 实施与验证(2026-10-23,F1 落地;F2 按决策不加)

**实现**(`src/AR/information.h` `calcLogHSumGaugeIter`,新增 `_ds_scale` 参数):

- 双曲层 u_min 改半长轴规范 U=Gm₁m₂/|a| —— 与 e 无关(§2 推导中 e−1 因子的
  唯一来源),与椭圆分支同口径;u·p_eff=2πL ⇒ 单层 S=π/32,构造性合规。
- p_eff 下限 t_cross=2r_now/v_now 保留(出腿方向正确,防飞越欠分辨)。
- S 上限实装 **S_cap=0.25**(设计稿 0.5,实装收紧,对高阶合成最大系数留
  裕度):`S_level=_ds_scale·p_eff·u_eff/(64·L_level)`,超限则
  `u_eff×(S_cap/S_level)`(除 64 = 半步口径)。
- 首估由此对垃圾双曲拟合构造性免疫 —— 任何 e→1 拟合的层贡献被 cap 在
  S 安全值(plan item 4 的 LogH 侧随之闭合)。
- F2(NaN 瞬态闸门)未实现 —— F1 后 merger 重现器 0 NaN/0 负步长,兜底
  暂无必要(用户决策)。

**验证**:

- 确定性:merger snapshot 重启重现器,F1 前逐位复现 NaN/负步长磨步,F1 后
  0 NaN、0 负 dt、中断正常触发。
- 生产受控 A/B(g5 vs g6;同 etot-fix 基座、同 LogH、唯一差异 F1):
  step 警告 168→6,墙钟 43.1→9.5 min,err/total 同量级(-4.4e-7 vs
  -5.0e-8)。单次全程对比受混沌发散影响,与确定性重现器互为印证。
- 规范对照:g3(BTLogH,item-5 修复)同 IC 全程 0 警告 / 10.0 min /
  1.15e-6 —— 磨步类为 LogH 路径现象;T8 r1/r2 "cascade" 类(同一颗
  105 Msun 致密双星)dump 重放 BTLogH 每事件便宜 38–934×。

## 7. 残留缺陷:控制器恢复注入极限环(LogH 路径,plan item 6 机制记录)

dump 重放(`g5.dump_large_step_h4_n35_g13_t1.90`,AR_COLLECT 构建)修正了早期
"ds 无关误差地板"的读法 —— 逐步误差**与 ds² 成正比**(级联本身健康),
并非地板:

- **循环机制**:备份恢复(`countAndRecover`/`Reuse_backup_ds`)把级联起点
  过大的 ds 放回 → 首步误差 1e2–1e7× 阈值 → 健康折半级联(err∝ds²)→
  渐近值钉在 1.00–1.02× 阈值 → 静默期 → 恢复再次注入。阈值扫描
  1e-6/1e-8/1e-10 → 710/127/118 事件(ds 5e-6/7.7e-8)—— 非 NaN、非负步长、
  非独立地板。
- **磨步组解剖**(g6):4 体层级四重系 —— 105 Msun 致密双星(a≈6.5e-7,
  κ≈4e5)+ 塌缩中轻对(semi 3.8e-4→3.7e-5,e 0.99→0.87)+ 宽根(e≈0.94)。
  slowdown 使时间尺度排序倒置:轻对塌缩(周期 ~1.6e-5)快于慢化内双星
  (有效 ~1.95e-3),估计器的"未慢化最小周期"结构无法表达(估计器 vs
  控制器 ~120× 分歧)。
- **修复方向(open,未批准)**:有界恢复(恢复 ds 上限 ×2 最后良好值)/
  生存窗 / 估计器按慢化有效周期重排序。生产影响已收窄:LogH+F1 残留
  6 警告,BTLogH 0。
