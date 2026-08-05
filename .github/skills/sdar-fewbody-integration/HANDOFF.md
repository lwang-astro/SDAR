# SDAR Skill 开发交接说明（2026-08-05 更新）

本文件用于在新电脑/新 VS Code 会话中快速恢复当前 SDAR skill 开发状态。

## 1) 当前完成状态

已完成：
- 建立 SDAR skill 主文件（`SKILL.md`），覆盖 AR、Hermite、Kepler 三大组件。
- 编写 sample run scripts（binary_logh.sh, triple_logh_sd.sh, fewbody_hermite.sh, triple_compare_methods.sh, build_kepler_tree.sh）。
- 创建 assets 文档：binary-scenario-map.md（方法选择指南）、data-readback-patterns.md（Python 数据读取范例）。

## 2) 核心文件

主规则文件：
- `.github/skills/sdar-fewbody-integration/SKILL.md`

能力资产文件：
- `.github/skills/sdar-fewbody-integration/assets/binary-scenario-map.md` — AR/Hermite 变体选择指南
- `.github/skills/sdar-fewbody-integration/assets/data-readback-patterns.md` — Python 数据读取验证范例

示例脚本：
- `sample/input/binary_logh.sh` — 双星 AR LogH 积分
- `sample/input/triple_logh_sd.sh` — 三星 AR LogH + slowdown
- `sample/input/fewbody_hermite.sh` — Hermite+AR 混合积分
- `sample/input/triple_compare_methods.sh` — 多方法对比
- `sample/input/build_kepler_tree.sh` — Kepler 树构建

## 3) 新电脑恢复步骤

```bash
cd /home/lwang/code/SDAR

# 编译全部示例代码
make -C sample/AR && make -C sample/AR install
make -C sample/Hermite && make -C sample/Hermite install
make -C sample/Kepler && make -C sample/Kepler install

# 安装 Python 工具
make -C tools

# 验证 Python 工具可用
python3 -c "import sys; sys.path.append('/home/lwang/include'); import sdar; print('SDAR tools OK')"
```

## 4) 在新会话里如何"接上上下文"

建议在 VS Code Chat 首条消息直接粘贴：

```text
请读取 SDAR/.github/skills/sdar-fewbody-integration/HANDOFF.md，
并基于 SKILL.md 继续维护 SDAR 相关工作。
先确认 sample/AR、sample/Hermite、sample/Kepler 的编译状态和 tools/ 的安装状态。
```

## 5) 当前约定（重要）

- SDAR 是 PeTar 的参考实现，PeTar 中关于 SDAR 的参数（如 `--r-group`、`--r-search-group`）映射到 SDAR 的组检测参数。
- AR 变体（logh/ttl/sd.t/kdk.pert）按场景选择，见 `assets/binary-scenario-map.md`。`.sd.a` 已废弃，不在当前 Makefile 中。
- Python 工具安装在 `/home/lwang/include/sdar/`（通过 `tools/Makefile`）。
- SDAR 不涉及 MPI、外势、恒星演化 —— 这些只在 PeTar 层面存在。
- SDAR 的 `-G` 默认值为 1.0（Henon 单位），物理单位下使用 0.00449830997959438。

## 5a) 已验证的 Python 读取约定（2026-07-21 验证修正）

以下约定经实际运行 5 个 sample 脚本 + Python 读回验证，SKILL 已同步：

1. **正确的 API 模式：先构造再 loadtxt。** `SDARData(N_particle=N, time_measure=True)` 然后 `data.loadtxt(file, skiprows=1)`。kwargs 不能传入 `loadtxt`（会报错）。
2. **Hermite 输出需预过滤。** 独立 `hermite` 的诊断消息混入数据流，过滤到数字起始行后可正确读取。
3. **列数随 AR 变体而异，N_sd 可匹配。** N=3 时：plain AR = 56 列，`.sd.t` = 75 列。`slowdown=True, N_sd=2` 可完美匹配 `.sd.t` 的 75 列。废弃的 `.sd.a` 不再使用。
4. **`findPair` 签名因 `use_kdtree` 而异。** `True`：`(kdt, singles, binary)`，`False`：`(singles, binary)`。
5. **HermiteData 访问器为 `data.energy_phy`，非 `data.energy`。**
6. **修改 Jupyter notebook 后必须验证 cell 完整性。** `edit_notebook_file` 的 `replace_string_in_file` 操作容易导致 cell 内容被错误截断（如丢失函数体后半段、缺失 for 循环头、留下孤立 try/except）。每次修改后必须执行被修改的 cell 确认无语法错误和运行时错误。此规则已写入 SKILL.md 的 Non-Negotiable Rules。

## 6) 下一步可继续做的增强

- 增加 `assets/minimal-question-sets.md` 用于最小提问清单（AR 场景、Hermite 场景、Kepler 场景）。
- 增加自动化编译状态检查脚本（检查所有 sample 子目录的 build 状态）。
- 增加 Python Jupyter Notebook 示例（参考 PeTar 的 `sample/data_analysis.ipynb`）。
- 增加 BlogH 混合方法的相关文档（`docs/hierarchical_blogh_impl_notes.md` 中已有实现笔记）。

## 7) PeTar ↔ SDAR 版本对应关系

PeTar 的 VERSION 文件格式为 `PeTar版本_SDAR版本`（如 `1745e_423e`）。
当调试 PeTar 中 SDAR 行为时，确认 PeTar 内嵌的 SDAR 版本与独立 SDAR 版本一致。
