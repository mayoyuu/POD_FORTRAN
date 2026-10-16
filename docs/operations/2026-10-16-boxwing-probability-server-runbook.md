# DRO-2 8D 概率传播：Git 迁移与无 agent 服务器操作手册

日期：2026-10-16。当前已经通过本机 50k/2k 两小时 smoke；新的 360 h 概率传播尚未启动。

本手册中的服务器命令需要由使用者执行。手册编写本身未代执行服务器操作；后续本地提交/运行状态以实际 Git 和 C2 输出为准，不代表已推送或部署服务器。

## 0. 先确认架构，不要把“能 Git clone”当作“能运行”

在服务器执行并保存输出：

```bash
uname -m
cat /etc/os-release
ldd --version | head -1
python3 --version
df -h "$HOME"
```

本机实际环境：Ubuntu 24.04 x86_64、gfortran 13.3、Python 3.12、fpm 0.12.0、NumPy 2.5.4、h5py 3.16.0。冻结的 P1 extractor 要求的符号包括 GLIBC_2.38、GLIBCXX_3.4.32、GFORTRAN_10。

选择执行环境：

| 服务器 | 采用路线 |
| --- | --- |
| Ubuntu 24.04 x86_64 | 按 §4A 原生环境继续，仍需实际运行审计/smoke |
| 其他 x86_64 Linux | 先确认 ABI；不兼容时用 §4B Ubuntu 24.04 amd64 容器 |
| ARM/aarch64 | 不能直接链接本机 `.a` 或执行本机 extractor；可按 §4B 在**管理员已配置 amd64 仿真**的 Docker 中尝试；必须先通过两小时 smoke |
| ARM 且无可用仿真/容器 | 本手册的已验证依赖复用路线不可用，停在此处；需要额外进行 ARM 原生移植验证，或换 x86_64 计算节点 |

ARM 原生不是“换个编译器”即可保证：DACE、项目 wrapper、Fortran SPICELIB 均需原生构建；既有 P1 的冻结 x86 extractor 也无法执行。不能把历史 metadata 的 executable hash 随意改成新 ARM 二进制来假装历史身份一致。需要单独定义和验证原生 P1/历元审计路线。本手册不宣称完成了该移植。

Docker 的 amd64 镜像在 ARM 上仍需要仿真；`--platform` 本身不会解决缺失的仿真。QEMU 对编译/数值计算可能明显更慢。不要在没有 smoke 结果时直接提交长算。[Docker 官方说明](https://docs.docker.com/build/building/multi-platform/)

## 1. 本机代码：必须从正确的 worktree 提交

本机主仓库和代码工作树不是同一个目录：

```text
主仓库/数据根：/home/songyu/fortran_test/POD_Fortran
代码工作树：  /home/songyu/fortran_test/POD_Fortran/.worktrees/boxwing-modes-p1
当前分支：    codex/boxwing-modes-p1
远程：        git@github.com:mayoyuu/POD_FORTRAN.git
```

手册最初编写时，代码尚未提交，HEAD 为 `8d999bfc317eee44e22f36e56e0a2e2e87f493e1`。若之后已执行本地提交，先用 `git log -1` / `git status` 确认，跳过重复 commit，按实际新 SHA 推送和部署。仅 `git push` 不会上传未提交的功能。下面在**本机 WSL Ubuntu/Bash** 执行，不是在 Windows PowerShell 直接执行：

```bash
cd /home/songyu/fortran_test/POD_Fortran/.worktrees/boxwing-modes-p1
set -eo pipefail
git branch --show-current
git status --short
git diff --check

# 包含 P1–P4 及概率功能的新源码。先人工检查，不用 git add .。
git add src app config tools test external fpm.toml setup_env.sh .gitignore
# docs/ 被忽略；只强制加入这份明确需要迁移的手册。
git add -f docs/operations/2026-10-16-boxwing-probability-server-runbook.md
git diff --cached --name-status
git diff --cached --check

# 确认暂存区没有 OPM/SRP/kernels/lib 的符号链接，没有历史数据、venv 或 build。
git commit -m "feat: add audited 8D modes probability output"
git push -u origin codex/boxwing-modes-p1
git rev-parse HEAD
```

记录最后的**新 commit SHA**。不要继续用旧的 `8d999...` 作为部署版本。不要为了省事把 `OPM`、`SRP`、`kernels`、`lib` 加入 Git：当前 worktree 中它们是指向本机位置的未跟踪链接；Git 即使记录链接也不会上传链接所指的文件。

至少确认下列文件已经被 Git 跟踪：

```bash
git ls-files --error-unmatch \
  app/run_srp_modes_ads.f90 \
  src/lib/uncertainty/io/pod_srp_modes_probability_io_module.f90 \
  tools/run_srp_modes_ads.py \
  tools/generate_srp_modes_probability_samples.py \
  tools/pack_srp_modes_probability_hdf5.py \
  tools/requirements-srp-modes-probability.txt
```

## 2. 本机导出独立运行资源包

### 2.1 为什么要单独导出

Git 不包含 OPM、SPICE 内核、重力场、模式 CSV、概率样本和静态库。P1 metadata 还有两个绝对路径：`opm_file`、`config_file`。这里在**新导出副本**中改为相对路径，并保留原 metadata 及其 SHA256；不修改本机已发布结果。

已有模式 CSV、OPM、P1 源配置和 extractor 的字节/哈希都保持不变。该重定位方式已在本机临时目录通过现有 `audit_table` 验证。迁移 metadata 自身的哈希会改变，服务器的新 C2 运行会如实记录新身份。

### 2.2 在本机 WSL 导出

```bash
set -eo pipefail
export POD_LOCAL_MAIN=/home/songyu/fortran_test/POD_Fortran
export POD_LOCAL_WORKTREE="$POD_LOCAL_MAIN/.worktrees/boxwing-modes-p1"
export POD_EXPORT_DIR=/home/songyu/fortran_test/transfer/boxwing_probability_20261016
export POD_ARCHIVE=/home/songyu/fortran_test/transfer/boxwing_probability_20261016.tar.gz

test ! -e "$POD_EXPORT_DIR"
test ! -e "$POD_ARCHIVE"
git -C "$POD_LOCAL_WORKTREE" diff --quiet
git -C "$POD_LOCAL_WORKTREE" diff --cached --quiet
```

下面的 Python 只创建新导出目录，复制实际文件而不是符号链接：

```bash
python3 - <<'PY'
import hashlib, json, os, shutil, subprocess
from pathlib import Path

main=Path(os.environ['POD_LOCAL_MAIN'])
work=Path(os.environ['POD_LOCAL_WORKTREE'])
out=Path(os.environ['POD_EXPORT_DIR'])
out.mkdir(parents=True,exist_ok=False)
sha=lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
def copy(src,rel):
    target=out/rel
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(src,target)  # follows links: export regular files
    assert sha(src)==sha(target)

for name in ('libdace_s.a','libdace_wrapper.a','libspicelib.a'):
    copy(main/'lib'/name,'lib/'+name)
copy(Path(shutil.which('fpm')),'.runtime/bin/fpm')
for rel in ('kernels/lsk/naif0012.tls',
            'kernels/pck/pck00010.tpc',
            'kernels/pck/earth_latest_high_prec.bpc',
            'kernels/spk/de440.bsp',
            'kernels/gravity_models/GGM05C.GEO',
            'kernels/gravity_models/gggrx_0660pm_sha.tab',
            'OPM/DRO-2/DRO-2_init.opm.json'):
    copy(main/rel,rel)

mode=main/'SRP/boxwing_modes_DRO-2_init'
rel=Path('SRP/boxwing_modes_DRO-2_init')
meta_path=mode/'boxwing_modes_meta.json'
meta=json.loads(meta_path.read_text())
for key in ('opm','config','executable'):
    source=Path(meta[key+'_file'])
    if not source.is_absolute(): source=mode/source
    assert sha(source)==meta[key+'_sha256'], key+' identity mismatch'
assert sha(mode/'boxwing_modes.csv')==meta['csv_sha256']
for name in ('boxwing_modes.csv','extractor'):
    copy(mode/name,rel/name)
copy(meta_path,rel/'boxwing_modes_meta.original.json')
copy(Path(meta['opm_file']),rel/'source_opm.opm.json')
copy(Path(meta['config_file']),rel/'source_config.txt')
portable=dict(meta)
portable.update(opm_file='source_opm.opm.json',config_file='source_config.txt')
portable['relocation']={
    'original_metadata_sha256':sha(meta_path),
    'original_opm_file':meta['opm_file'],
    'original_config_file':meta['config_file'],
    'changed_fields':['opm_file','config_file']}
(out/rel/'boxwing_modes_meta.json').write_text(json.dumps(portable,indent=2,allow_nan=False)+'\n')

samples=Path('SRP/boxwing_modes_DRO-2_init_probability_samples_20261016')
for name in ('probability_samples.csv','probability_samples.json'):
    copy(main/samples/name,samples/name)
commit=subprocess.check_output(['git','-C',str(work),'rev-parse','HEAD'],text=True).strip()
(out/'.runtime/CODE_COMMIT.txt').write_text(commit+'\n')
print('EXPORT=',out)
print('CODE_COMMIT=',commit)
print('Original data and metadata were not changed.')
PY

# 下载与本机验证一致的 Linux x86_64 / Python3.12 wheels；ARM 路线仅在 amd64 容器中使用。
"$POD_LOCAL_WORKTREE/.venv-probability/bin/python" -m pip download \
  --only-binary=:all: --dest "$POD_EXPORT_DIR/.runtime/wheels" \
  'numpy==2.5.4' 'h5py==3.16.0'

(
  cd "$POD_EXPORT_DIR"
  find . -type f ! -name SHA256SUMS.inputs -print0 \
    | sort -z | xargs -0 sha256sum > SHA256SUMS.inputs
)
tar -czf "$POD_ARCHIVE" -C "$POD_EXPORT_DIR" .
(
  cd /home/songyu/fortran_test/transfer
  sha256sum boxwing_probability_20261016.tar.gz > boxwing_probability_20261016.tar.gz.sha256
)
```

只需要这一组概率输入和 hourly P1 模式数据，不需要搬全部历史 SRP、C0/C1 结果、旧 build 或本机 venv。尤其不要漏两个 `gravity_models` 文件：仅有 de440 不足以运行当前高阶引力配置。

## 3. 上传资源、在服务器 clone 并固定 commit

在本机 WSL 中设置真实 SSH 地址（这是唯一无法替你填写的地址）：

```bash
export POD_SERVER_SSH=your_user@your_server
ssh "$POD_SERVER_SSH" 'mkdir -p "$HOME/pod_inbox"'
scp "$POD_ARCHIVE" "$POD_ARCHIVE.sha256" "$POD_SERVER_SSH:pod_inbox/"
```

服务器执行。使用新目录，不覆盖旧的工程或结果：

```bash
set -eo pipefail
export POD_HOST_ROOT="$HOME/pod_runs/POD_FORTRAN_20261016"
test ! -e "$POD_HOST_ROOT"
mkdir -p "$HOME/pod_runs"
git clone --single-branch --branch codex/boxwing-modes-p1 \
  git@github.com:mayoyuu/POD_FORTRAN.git "$POD_HOST_ROOT"

cd "$HOME/pod_inbox"
sha256sum -c boxwing_probability_20261016.tar.gz.sha256

# 本手册 §2 生成的可信自有压缩包，只解到新的 clone 中。
tar -xzf boxwing_probability_20261016.tar.gz -C "$POD_HOST_ROOT"
cd "$POD_HOST_ROOT"
sha256sum -c SHA256SUMS.inputs
export POD_RELEASE_COMMIT="$(tr -d '\r\n' < .runtime/CODE_COMMIT.txt)"
git cat-file -e "$POD_RELEASE_COMMIT^{commit}"
git switch --detach "$POD_RELEASE_COMMIT"
test "$(git rev-parse HEAD)" = "$POD_RELEASE_COMMIT"

# 只影响该 clone 的本地忽略规则，不修改提交。
printf '%s\n' '.runtime/' '.venv-probability/' 'SHA256SUMS.inputs' >> .git/info/exclude
```

若服务器没有该仓库的 SSH 读取权限，先由管理员/你配置密钥，或使用你有权限的 HTTPS clone URL。不要关闭主机密钥校验，不要把私钥放进 Git/输入包。

## 4. 建立运行环境（A/B 二选一）

### 4A. Ubuntu 24.04 x86_64 原生

由有权限的人安装系统依赖：

```bash
sudo apt-get update
sudo apt-get install -y git build-essential gfortran libopenblas-dev \
  libgfortran5 libstdc++6 python3.12 python3.12-venv tmux binutils
export POD_SERVER_ROOT="$POD_HOST_ROOT"
cd "$POD_SERVER_ROOT"
```

不要重新使用旧 README 的根目录 CMake 安装流程。本任务用 FPM；`external/dace_build/CMakeLists.txt` 还含本机 DACE 路径，本路线不运行它，而复用已验证 `.a`。

### 4B. ARM 或旧 Linux：Ubuntu 24.04 amd64 容器兼容路线

前提：管理员允许 Docker；ARM 主机已经配置可运行 amd64 的仿真。先测：

```bash
docker run --rm --platform linux/amd64 ubuntu:24.04 uname -m
```

必须成功得到 `x86_64`。若报 `exec format error`、权限不足或仿真缺失，**停下**，请管理员提供环境或换 x86_64 节点；不要直接运行后面的 `.a`。本手册不自动修改主机 binfmt 或启用特权容器。

服务器主机创建独立容器，主机数据通过 bind mount 保留：

```bash
export POD_CONTAINER=pod-boxwing-20261016
docker run -d --name "$POD_CONTAINER" --platform linux/amd64 \
  --mount "type=bind,source=$POD_HOST_ROOT,target=/work/POD_Fortran" \
  --workdir /work/POD_Fortran ubuntu:24.04 sleep infinity
docker exec --user root "$POD_CONTAINER" bash -lc \
  'apt-get update && apt-get install -y git build-essential gfortran libopenblas-dev libgfortran5 libstdc++6 python3.12 python3.12-venv tmux binutils'

# 在容器中选用与主机同 UID 的账号，防止产生不可访问的 root-owned 结果。
docker exec --user root "$POD_CONTAINER" bash -c '
  uid="$1"; gid="$2"
  if ! getent group "$gid" >/dev/null; then groupadd -g "$gid" podgroup; fi
  if ! getent passwd "$uid" >/dev/null; then
    useradd -m -u "$uid" -g "$gid" -s /bin/bash podrunner
  fi
' bash "$(id -u)" "$(id -g)"
export POD_CONTAINER_USER="$(docker exec "$POD_CONTAINER" getent passwd "$(id -u)" | cut -d: -f1)"
docker exec -it "$POD_CONTAINER" su - "$POD_CONTAINER_USER"
```

从这里起，在**容器内**执行后续命令：

```bash
export POD_SERVER_ROOT=/work/POD_Fortran
cd "$POD_SERVER_ROOT"
uname -m
```

容器内 `/work/POD_Fortran/...` 对应主机 `$POD_HOST_ROOT/...`。不要在主机 ARM Python 中创建/运行这份 amd64 venv。容器不能在任务结束前删除；stop/start 后 apt 环境保留，数据始终在主机目录。该路线的本项目实际结果仍需以服务器 smoke 为准，本机没有 ARM/QEMU 测试环境。

## 5. 运行前的目录与环境

此后所有命令都在**所选执行环境**内运行，且 cwd 始终为 `$POD_SERVER_ROOT`。

```text
$POD_SERVER_ROOT/
├── .git/                         固定提交的代码
├── app/ src/ tools/ config/ test/
├── .runtime/bin/fpm              已验证 Linux x86_64 静态 FPM
├── .runtime/wheels/              cp312 Linux x86_64 wheels
├── .runtime/CODE_COMMIT.txt       部署源码 commit
├── lib/
│   ├── libdace_s.a
│   ├── libdace_wrapper.a
│   └── libspicelib.a
├── kernels/
│   ├── lsk/naif0012.tls
│   ├── pck/pck00010.tpc
│   ├── pck/earth_latest_high_prec.bpc
│   ├── spk/de440.bsp
│   ├── gravity_models/GGM05C.GEO
│   └── gravity_models/gggrx_0660pm_sha.tab
├── OPM/DRO-2/DRO-2_init.opm.json
├── SRP/boxwing_modes_DRO-2_init/
│   ├── boxwing_modes.csv
│   ├── boxwing_modes_meta.json             迁移副本，两个输入路径为相对路径
│   ├── boxwing_modes_meta.original.json    未修改历史 metadata
│   ├── source_opm.opm.json                 P1 原始输入字节副本
│   ├── source_config.txt                  P1 原始配置字节副本
│   └── extractor                          原冻结 P1 程序，未改二进制/hash
├── SRP/boxwing_modes_DRO-2_init_probability_samples_20261016/
│   ├── probability_samples.csv            已验证 50000 个固定 8D 样本
│   └── probability_samples.json           先验/种子/映射/hash
├── .venv-probability/                      在当前执行环境新建
├── build/                                 当前执行环境重新编译
└── logs/                                  外部启动日志，不在新结果目录内
```

初始化。注意 `setup_env.sh` 引用一些可能尚未设置的变量，因此这里不启用 `set -u`：

```bash
cd "$POD_SERVER_ROOT"
set -eo pipefail
export PATH="$POD_SERVER_ROOT/.runtime/bin:$PATH"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
source setup_env.sh
mkdir -p logs

sha256sum -c SHA256SUMS.inputs
fpm --version
gfortran --version
python3.12 --version
ldd SRP/boxwing_modes_DRO-2_init/extractor
```

`ldd` 不得出现 `not found`；真正的 ABI 判据还包括下面的 P1 程序执行。不要因 fpm 帮助成功就认为 extractor 也兼容。

安装 Python 依赖到新 venv，不搬本机 venv：

```bash
test ! -e .venv-probability
python3.12 -m venv .venv-probability
.venv-probability/bin/python -m pip install --no-index \
  --find-links .runtime/wheels 'numpy==2.5.4' 'h5py==3.16.0'
.venv-probability/bin/python -c 'import numpy,h5py; print(numpy.__version__,h5py.__version__)'

# 验证原模式表、相对路径、冻结 extractor、UTC→TDB、CSV/物理门槛与输入 hash。
.venv-probability/bin/python tools/audit_boxwing_mode_table.py \
  --csv SRP/boxwing_modes_DRO-2_init/boxwing_modes.csv \
  --meta SRP/boxwing_modes_DRO-2_init/boxwing_modes_meta.json \
  > logs/p1_preflight.json

# 只编译/help，不开始长算；服务器上不要复用本机 build/runner。
fpm run --target run_srp_modes_ads -- --help
export POD_C2_EXE="$(find "$POD_SERVER_ROOT/build" -type f -path '*/app/run_srp_modes_ads' -print)"
test -x "$POD_C2_EXE"
printf 'C2 executable: %s\n' "$POD_C2_EXE"
df -h "$POD_SERVER_ROOT"
```

全 360 h、order4、max-depth8 的 wrapper 会保守计入最坏补丁数的每小时快照，当前估算预检约需 **47 GB 可用空间**，并在运行中维持 2 GB 余量。因此建议当前文件系统至少预留 **50 GB**；不能只按约 98 MB 的点状态主体估算磁盘。

建议为单个作业提供至少 8 GB 内存，运行中观察实际峰值；这是部署余量建议，不是已测得的完整新概率算例峰值。当前主要传播循环并未并行化，先按单作业/单核心运行，不要盲目同时开多个复制目录争用磁盘。

可选但建议的定向回归：

```bash
fpm test --target test_srp_modes_probability_io
.venv-probability/bin/python -m unittest discover -s test -p 'test_srp_modes_probability_*.py' -v
```

应有 15 项概率 Python 测试通过，原生 IO target 通过。这里不直接执行全仓 `fpm test`：历史有缺夹具/超时，而且部分旧测试固定引用原 worktree 的 `../../SRP` 布局，并不适合新的部署目录。不要为使测试表面通过而创建指向不存在本机位置的链接。

## 6. 先运行两小时 smoke，再决定是否进行 360 h

先进入持久终端：

```bash
tmux new -s pod-boxwing
```

如果新 shell 未继承变量，重新执行 §5 的环境导出（不要重复创建 venv）。在 tmux 内执行：

```bash
cd "$POD_SERVER_ROOT"
set -eo pipefail
export POD_SAMPLES_DIR="$POD_SERVER_ROOT/SRP/boxwing_modes_DRO-2_init_probability_samples_20261016"
export POD_RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
export POD_SMOKE_OUT="$POD_SERVER_ROOT/SRP/boxwing_modes_DRO-2_init_c2_probability_smoke_server_$POD_RUN_ID"
export POD_SMOKE_LOG="$POD_SERVER_ROOT/logs/probability_smoke_$POD_RUN_ID.wrapper.log"
test ! -e "$POD_SMOKE_OUT"

.venv-probability/bin/python -u tools/run_srp_modes_ads.py \
  --opm "$POD_SERVER_ROOT/OPM/DRO-2/DRO-2_init.opm.json" \
  --config "$POD_SERVER_ROOT/config/config_srp_modes_DRO-2.txt" \
  --executable "$POD_C2_EXE" \
  --output-dir "$POD_SMOKE_OUT" --hours 2 --save-hours 1 \
  --probability-samples "$POD_SAMPLES_DIR/probability_samples.csv" \
  --probability-manifest "$POD_SAMPLES_DIR/probability_samples.json" \
  --probability-count 50000 --probability-visualization-count 2000 \
  --probability-validation-count 64 \
  2>&1 | tee "$POD_SMOKE_LOG"
```

必须退出码为 0，且 `$POD_SMOKE_OUT/c2_report.json` 存在并满足：

```bash
.venv-probability/bin/python - <<'PY'
import json,os
from pathlib import Path
import h5py
p=Path(os.environ['POD_SMOKE_OUT'])
r=json.loads((p/'c2_report.json').read_text())
q=r['probability_outputs']
assert r['formal_result'] is True and r['duration_hours']==2
assert r['shape_count']==512 and r['validation_count']==337
assert q['statistics_sample_count']==50000 and q['visualization_sample_count']==2000
assert q['numerical_validation_passed'] and q['hourly_full_snapshot_count']==3
assert q['validation_count']==64 and q['snapshot_probe_checks']==192
assert q['sample_convergence_confirmed'] and q['probability_output_accepted']
with h5py.File(p/'c2_probability_points.h5') as h:
    assert h['visualization/state'].shape==(3,2000,6)
    assert h['keyframes/state'].shape==(2,50000,6)  # 无事件时
    assert h['moments/covariance'].shape==(3,3,6,6)
    assert h['keyframes/hour_index'][:].tolist()==[0,2]
print('PASS: server 50k/2k two-hour smoke')
print('Sample convergence:',q['sample_convergence_confirmed'])
PY
```

若服务器上实际出现事件，不能直接删除事件帧来凑该无事件断言；先查看 split/depth/probe 报告，查明与本机基线不同的原因，再决定长算。

## 7. 提交 360 h 完整概率传播

只在 §6 通过后，在同一 tmux shell/执行环境运行：

```bash
cd "$POD_SERVER_ROOT"
set -eo pipefail
export POD_RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
export POD_CASE_OUT="$POD_SERVER_ROOT/SRP/boxwing_modes_DRO-2_init_c2_probability_15d_server_$POD_RUN_ID"
export POD_JOB_LOG="$POD_SERVER_ROOT/logs/probability_15d_$POD_RUN_ID.wrapper.log"
test ! -e "$POD_CASE_OUT"
printf 'OUTPUT=%s\nCOMMIT=%s\n' "$POD_CASE_OUT" "$(git rev-parse HEAD)" | tee "$POD_JOB_LOG"

.venv-probability/bin/python -u tools/run_srp_modes_ads.py \
  --opm "$POD_SERVER_ROOT/OPM/DRO-2/DRO-2_init.opm.json" \
  --config "$POD_SERVER_ROOT/config/config_srp_modes_DRO-2.txt" \
  --executable "$POD_C2_EXE" \
  --output-dir "$POD_CASE_OUT" --hours 360 --save-hours 1 \
  --probability-samples "$POD_SAMPLES_DIR/probability_samples.csv" \
  --probability-manifest "$POD_SAMPLES_DIR/probability_samples.json" \
  --probability-count 50000 --probability-visualization-count 2000 \
  --probability-validation-count 64 \
  2>&1 | tee -a "$POD_JOB_LOG"
```

配置中的模式路径相对代码 ROOT，内核/重力场也依赖 ROOT 下的固定布局。不要挪走 ROOT 后单独拷一个 runner 到别的工作目录执行。不要绕过 Python wrapper 直接用原生 raw runner 作为正式结果。

`Ctrl+b` 然后 `d`：离开 tmux 而不终止任务。重新连接后 `tmux attach -t pod-boxwing`。容器路线先按 §4B 用 `docker exec ... su - ...` 进入原容器，再 attach。

若服务器是 Slurm/PBS 集群，应把相同环境/命令放进计算作业，不在登录节点长算；按管理员的队列、内存、walltime 和容器政策提交。本手册没有假设某个集群分区名或可用 Docker 权限。

## 8. 不依赖 agent 的监控与失败处理

原生传播写在隐藏 staging，全部审计成功才原子改名为 `$POD_CASE_OUT`。最终目录在前期不存在是正常现象。

```bash
# 启动器日志（tmux 内/另一个所选环境 shell）
tail -n 40 "$POD_JOB_LOG"

# 找当前 staging；用输出中的实际目录替换下一条路径。
find "$POD_SERVER_ROOT/SRP" -maxdepth 1 -type d \
  -name '.boxwing_modes_DRO-2_init_c2_probability_15d_server_*.staging-*' -print
tail -f /实际找到的staging目录/run.log
df -h "$POD_SERVER_ROOT"
```

注意：日志中的 `PASS: C2 raw run` 或 `hour=360` 只表示原生传播结束；随后还有 361 份快照、概率探针、HDF5 和身份审计。不能据此提前杀进程或开始改源码。

运行中禁止 `git pull/checkout`、改源码/输入/kernel、更新 Python 包、重编译/替换本次绑定 executable。输入或源码身份变化会使本次结果拒绝发布。

失败后保留 staging、run.log、wrapper.log 和 `.lock` 线索。不要盲目删除锁/半成品，不要手工把 raw report 改为 `formal_result=true`。当前没有经验证的中断续算功能；修好原因后用**新的 run ID / 输出目录**从头运行，旧 staging 保留诊断。

## 9. 完成判定：数值成功与概率充分性分别检查

将 `$POD_CASE_OUT` 设置为日志 `OUTPUT=` 记录的实际绝对目录，然后执行：

```bash
.venv-probability/bin/python - <<'PY'
import csv,hashlib,json,os
from pathlib import Path
import h5py
import numpy as np
p=Path(os.environ['POD_CASE_OUT'])
r=json.loads((p/'c2_report.json').read_text())
q=r['probability_outputs']
assert r['formal_result'] is True and r['duration_hours']==360 and r['save_hours']==1
assert r['shape_count']==512 and r['validation_count']==337
assert r['da_handles_baseline']==r['da_handles_final']==0
assert not r['validation_failed'] and not r['precision_unconfirmed']
assert q['statistics_sample_count']==50000 and q['visualization_sample_count']==2000
assert q['numerical_validation_passed'] and q['validation_count']==64
assert q['hour_count']==361 and q['statistics_row_count']==1083
assert q['hourly_full_snapshot_count']==361 and q['snapshot_probe_checks']==361*64
with (p/'c2_probability_hourly_stats.csv').open() as f:
    rows=list(csv.DictReader(f))
assert len(rows)==1083
with h5py.File(p/'c2_probability_points.h5') as h:
    assert np.array_equal(h['hour'][:],np.arange(361))
    assert h['visualization/state'].shape==(361,2000,6)
    assert h['moments/covariance'].shape==(361,3,6,6)
    assert np.array_equal(h['moments/sample_count'][:],[10000,20000,50000])
    keys=h['keyframes/hour_index'][:]
    assert {0,24,72,168,240,360}<=set(keys.tolist())
    assert h['keyframes/state'].shape==(len(keys),50000,6)
    assert np.all(h['validity/ads_accuracy_valid'][:])
    assert np.all(h['validity/precision_confirmed'][:])
for name,expected in r['artifacts_sha256'].items():
    with (p/name).open('rb') as f:
        assert hashlib.file_digest(f,'sha256').hexdigest()==expected,name
print('PASS: complete 360h numerical and artifact checks')
print('sample_convergence_confirmed=',q['sample_convergence_confirmed'])
print('probability_output_accepted=',q['probability_output_accepted'])
if not q['sample_convergence_confirmed'] or not q['probability_output_accepted']:
    raise SystemExit('NOT ACCEPTED as probability result: inspect convergence JSON; do not hide or auto-rerun')
print('PASS: probability output accepted for the specified prior and 360h arc')
PY
```

`formal_result=true` 不替代 `probability_output_accepted`。即使主统计收敛通过，P99.9 仍只有约 50 个尾部样本的量级，不是极罕见事件概率保证。

## 10. 输出清单、备份与回传

正式目录至少包含：

| 文件 | 内容 |
| --- | --- |
| `c2_report.json` | 正式状态、数值门槛、源码/输入/执行文件身份与所有最终文件哈希 |
| `c2_probability_points.h5` | 每小时 2k 状态/误差/patch_id、关键帧全量 50k、完整初始总体、每小时各前缀均值与完整协方差 |
| `c2_probability_hourly_stats.csv` | 1083 行，包含全 36 项协方差、偏度/峰度、分位数、有效性 |
| `c2_probability_probes.csv` | 每小时 64 个概率样本的 ADS/Real 验证 |
| `c2_probability_convergence.json` | 全弧 20k→50k 稳定性与充分性 |
| `probability_samples.csv/json` | 本次冻结的输入总体与先验身份 |
| `c2_history.csv`, `c2_shape.csv`, `c2_probes.csv`, `c2_split_events.csv` | 原有 512/337/零模式与事件输出 |
| `c2_patch_s*_meta/map/splits/coeff.csv` | 每小时完整 v2 8D ADS 快照，文件名由原程序定义 |
| `runner`, `run.log` | 冻结的执行文件、原生运行日志 |

点 staging CSV 默认在 HDF5 成功核验后删除，这是正常行为，不代表没有保存点云。没有每小时 50k 点大文件；每小时全量统计 N 始终为 50000。

建议备份**整个正式目录**和相应 wrapper 日志、`SHA256SUMS.inputs`、部署 commit，而非只取 HDF5。要保留未来离线重采样能力，不能删掉逐小时快照。

回传（本机 WSL，替换为日志记录的实际主机路径；容器 `/work/...` 要换回主机路径）：

```bash
export POD_RETURN_DIR=/home/songyu/fortran_test/server_results/boxwing_20261016
mkdir -p "$POD_RETURN_DIR"
rsync -a --info=progress2 \
  your_user@your_server:/home/你的服务器用户名/pod_runs/POD_FORTRAN_20261016/SRP/实际15d结果目录/ \
  "$POD_RETURN_DIR/"
```

不要把新服务器结果覆盖回本机旧的 `*_20261014`、`*_20261015` 或已发布 smoke。迁移不需要重跑 C0/C1，也不把轨道改成 DRO-3。

## 11. 遇错时按此表处理

| 现象 | 检查/处理 |
| --- | --- |
| `Exec format error` / `file in wrong format` | 架构不符；确认运行在 amd64 环境，不是直接在 ARM 上使用本机库 |
| `GLIBC_2.38` / `GLIBCXX_3.4.32 not found` | 旧系统 ABI；用兼容容器/节点，不要覆盖系统 glibc |
| `missing source opm/config` | 没有使用 §2 的重定位副本，或漏了两个 `source_*` 文件 |
| extractor/hash mismatch | 没带原冻结 extractor，或文件被替换；不能随意改 metadata hash |
| `cannot find -lspicelib/-ldace_wrapper/-ldace_s` | `lib` 缺三份静态库，或没在 ROOT source 环境 |
| `cannot find -lopenblas` | 所选运行环境缺开发包/链接库 |
| 内核/重力场警告 | 核对 §5 六个固定文件及 cwd；不要带警告启动长算 |
| `ModuleNotFoundError numpy/h5py` | 当前 shell 用错 Python；必须 `.venv-probability/bin/python` |
| `insufficient probability output disk space` | 最坏快照预检需约 47GB；选择 ≥50GB 可用的文件系统 |
| `output already exists` | 改新 run ID；不删、不覆盖已发布目录 |
| 最终目录还没有，但 run.log 到了 360 h | 审计/打包阶段可能还在进行；检查进程和 wrapper.log |
| `inputs/sources changed` | 运行期间发生改动；保留失败 staging，用冻结版本新目录重跑 |
| `probability_output_accepted=false` | 读 convergence JSON；数值完整不代表概率充分性通过 |

## 参考与验证边界

- [FPM 官方安装文档](https://fpm.fortran-lang.org/install/index.html)：二进制/源码安装方式。本路线复用本机 0.12.0 可执行文件，不要求猜最新版本。
- [NAIF Fortran Toolkit 平台列表](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)：当前列出的 Linux Fortran 分发为 PC 平台，不能据此宣称现有包已支持 ARM 原生。
- [Docker 多平台/仿真说明](https://docs.docker.com/build/building/multi-platform/)：amd64/arm64 区别、仿真前提与性能限制。

本机已验证：代码状态/忽略规则、输入清单、冻结二进制 ABI、相对路径重定位审计、50k/2k 两小时实际传播。本手册未在用户服务器执行；服务器的权限、架构、容器/仿真与完整长算需要按照上述门槛实测。
