#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
SHARC 轨迹 -> HDF5（含 expec.out 解析，列名+单位写入，ΔE=动能口径）
===============================================================================

【功能】
1) 遍历多个 ROOTS 下的状态/轨迹目录，解析每条轨迹的
   - geometry (Å)、velocity (a.u.)，来自 output.dat
   - expec.out 原文 + 表格，并把“表头+单位”写入 energy_table.attrs["columns"]
   - 轨迹属性（elements、atomic_masses、nstep、mtime、dir_path、xyz_path、state、
               multiplicity、碎片信息等）
   - delta_energy（单位 eV）：默认 ΔE = Ekin(final) − Ekin(initial)

2) 两种输出模式：
   - MODE=1（或输出名以“.1.hdf5”）：仅 geometry/velocity + 属性
   - MODE=2（或输出名以“.2.hdf5”）：再附 expec_out_raw 与 energy_table

3) 执行逻辑分离（Phase 1 & Phase 2）：
   - 在提取数据前，会先集中检查并生成缺失的 output_data 目录，避免日志混乱与写入冲突。

【用法】
- 直接运行：    python extrac_hdf5.py
- 常改项：见下方“配置区”——ROOTS、OUTPUT_PATH、MODE、N_WORKERS、DELTA_E_POLICY 等

===============================================================================
"""

# ============================== 配置区（只改这里） ==============================

# 1) 数据根目录（可多个）
ROOTS = [
    "/data/dliu/CO23+/07.newtest",
]

# 2) 输出 HDF5 文件
OUTPUT_PATH = "Data.2.hdf5"    # 建议：Data.1.hdf5 / Data.2.hdf5

# 3) 模式：1=轻量(geom+vel)，2=完整(含 expec)；None=按输出名后缀推断
MODE = None

# 4) 并行核数（0 或 1 表示单线程）
N_WORKERS = 32

# 5) ΔE 口径（与 demo.py 对齐）
#   "ekin_delta" : ΔE = Ekin(final) − Ekin(initial)  （默认，常用于近似 KER）
#   "ekin_final" : ΔE = Ekin(final)                  （若你把最终动能直接当 KER）
#   "auto"       : 先按 ekin_delta，不行再 ekin_final，再 epot，再 etot
DELTA_E_POLICY = "ekin_delta"

# 6) 若 MODE=2，需要 expec.out；如果没有，是否尝试调用外部程序生成？
GENERATE_EXPEC = True
DATA_EXTRACTOR = "$SHARC/data_extractor.x"
GEN_CMD = "{exe} output.dat"

# 7) 目录名匹配（不区分大小写）
STATE_DIR_REGEXES = [
    r"(?i)^state_\d+",
    r"(?i)^singlet_\d+",
    r"(?i)^triplet_\d+",
    r"(?i)^doublet_\d+",
    r"(?i)^duplet_\d+",
    r"(?i)^s_\d+",
    r"(?i)^t_\d+",
]
TRAJ_DIR_REGEX = r"(?i)^traj_\d+"

# 8) 其他常量
BOHR_TO_A = 0.529177210903
AU_TO_FS  = 0.02418884326505
FRAG_RCUT = 3.5   # 末帧聚类距离阈值（Å）

# ============================ 代码区（一般不动） ============================

import os, re, sys, time, json, subprocess
from pathlib import Path
import numpy as np
import h5py

def _compile_regexes(patterns):
    return [re.compile(p) for p in patterns]

RE_STATE = _compile_regexes(STATE_DIR_REGEXES)
RE_TRAJ  = re.compile(TRAJ_DIR_REGEX)

def f2(x: str) -> float:
    return float(x.replace("D","E").replace("d","E"))

def parse_output_dat(path: Path):
    """从 output.dat 解析 geometry/velocity、元素、质量、步长等。"""
    try:
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception as e:
        return {"ok": False, "msg": f"read-failed: {e}", "path": str(path)}

    dt_au, nsteps, natom = None, None, None
    for ln in lines:
        m = re.search(r"\bdtstep\s+([\-0-9\.EedD\+]+)", ln);   dt_au   = f2(m.group(1)) if m else dt_au
        m = re.search(r"\bnsteps\s+(\d+)", ln);                nsteps  = int(m.group(1)) if m else nsteps
        m = re.search(r"\bnatom\s+(\d+)", ln);                 natom   = int(m.group(1)) if m else natom
    dt_fs = (dt_au or 20.670686894780374) * AU_TO_FS

    elements, masses = [], []
    reading_el = reading_mass = False
    for ln in lines:
        s = ln.strip()
        if s.startswith("! Elements"):         reading_el, reading_mass = True, False;  continue
        if s.startswith("! Atomic masses"):    reading_el, reading_mass = False, True;  continue
        if s.startswith("********************************"): reading_el = reading_mass = False; continue
        if reading_el:
            toks = [t for t in re.split(r"[\s,]+", s) if t]
            elements += toks
        if reading_mass:
            toks = [t for t in re.split(r"[\s,]+", s) if t]
            for t in toks:
                try: masses.append(f2(t))
                except: pass

    R_list, V_list = [], []
    i, N = 0, len(lines)
    while i < N:
        s = lines[i].strip()
        if s.startswith("! 11 Geometry"):
            block = []
            for j in range(1, 1 + (natom or 3)):
                parts = [p for p in re.split(r"[\s,]+", lines[i+j].strip()) if p]
                block.append([f2(parts[0]), f2(parts[1]), f2(parts[2])])
            R_list.append(np.array(block, dtype=float))
            i += 1 + (natom or 3);  continue
        if s.startswith("! 12 Velocities"):
            block = []
            for j in range(1, 1 + (natom or 3)):
                parts = [p for p in re.split(r"[\s,]+", lines[i+j].strip()) if p]
                block.append([f2(parts[0]), f2(parts[1]), f2(parts[2])])
            V_list.append(np.array(block, dtype=float))
            i += 1 + (natom or 3);  continue
        i += 1

    if not R_list or not V_list:
        return {"ok": False, "msg": "no geometry/velocity blocks", "path": str(path)}

    R = np.array(R_list)            # bohr
    V = np.array(V_list)            # a.u.
    T = min(len(R), len(V))
    if natom is None: natom = R.shape[1]

    return {
        "ok": True,
        "path": str(path),
        "R_A": (R[:T] * BOHR_TO_A).astype(np.float64),
        "V_au": V[:T].astype(np.float64),
        "nstep": int(T),
        "natom": int(natom),
        "elements": np.array(elements[:natom], dtype=object) if elements else None,
        "atomic_masses": np.array(masses[:natom], dtype=np.float64) if masses else None,
        "dt_fs": float(dt_fs),
    }

# ---------- expec.out 解析（带表头+单位，写入 columns） ----------
def parse_expec_table(text: str):
    import re
    lines = [ln.rstrip() for ln in text.splitlines()]
    head = None
    for i in range(len(lines)-1):
        a = lines[i].lstrip("#! ").strip()
        b = lines[i+1].lstrip("#! ").strip()
        if "|" in a and "|" in b:
            head = i; break

    def split_bar_row(s):
        return [c.strip() for c in s.split("|") if c.strip()!=""]

    def join_name_unit(names, units):
        out=[]
        for n,u in zip(names, units):
            m = re.match(r"^(.*?)(?:\s+(\d+))?$", n.strip())
            base = m.group(1).strip()
            num  = m.group(2)
            if num: base = f"{base} {num}"
            u2 = u.strip()
            if u2.startswith("[") and u2.endswith("]"):
                out.append(f"{base} {u2}")
            else:
                out.append(base)
        return out

    if head is not None:
        names = split_bar_row(lines[head].lstrip("#! "))
        units = split_bar_row(lines[head+1].lstrip("#! "))
        K = min(len(names), len(units))
        columns = join_name_unit(names[:K], units[:K])
        data_start = head + 2
    else:
        columns, data_start = None, None
        for i, ln in enumerate(lines):
            if ln.startswith("#") or ln.startswith("!"):
                s = ln.lstrip("#! ").strip()
                if any(ch.isalpha() for ch in s):
                    columns = re.split(r"[,\s]+", s)
                    data_start = i+1
                    break
        if columns is None:
            columns=[]
        if data_start is None:
            data_start = 0

    rows=[]
    for ln in lines[data_start:]:
        if not ln or ln.lstrip().startswith("#") or ln.lstrip().startswith("!"):
            continue
        parts = [p for p in re.split(r"[,\s]+", ln.strip()) if p]
        try:
            rows.append([f2(p) for p in parts])
        except:
            pass
    if not rows:
        return None, None
    lens=[len(r) for r in rows]
    K=max(set(lens), key=lens.count)
    rows=[r for r in rows if len(r)==K]
    arr=np.array(rows, dtype=np.float64)

    if not columns or len(columns)!=K:
        columns=[f"col{i}" for i in range(K)]
    return arr, columns

def guess_expec_path(traj_dir: Path) -> Path:
    p1 = traj_dir / "output_data" / "expec.out"
    p2 = traj_dir / "expec.out"
    return p1 if p1.exists() else (p2 if p2.exists() else p1)

# ---------- 碎片聚类（末帧） ----------
def connected_fragments(last_R_A: np.ndarray, elems: np.ndarray, rcut=FRAG_RCUT):
    N = last_R_A.shape[0]
    adj = {i:set() for i in range(N)}
    for i in range(N):
        for j in range(i+1,N):
            if np.linalg.norm(last_R_A[i]-last_R_A[j]) <= rcut:
                adj[i].add(j); adj[j].add(i)
    seen=set(); frags=[]
    for i in range(N):
        if i in seen: continue
        stack=[i]; comp=[]
        while stack:
            u=stack.pop()
            if u in seen: continue
            seen.add(u); comp.append(u)
            for v in adj[u]:
                if v not in seen: stack.append(v)
        frags.append(sorted(comp))

    from collections import Counter
    def formula(idxs):
        syms=[str(elems[k]) for k in idxs]
        cnt=Counter(syms)
        out=[]
        if 'C' in cnt: out.append(('C', cnt.pop('C')))
        if 'H' in cnt: out.append(('H', cnt.pop('H')))
        for k in sorted(cnt.keys()): out.append((k, cnt[k]))
        def show(p): s,n=p; return s if n==1 else f"{s}{n}"
        return "".join(show(p) for p in out) if out else "X"

    frag_list=", ".join(formula(g) for g in frags)
    frag_json=json.dumps(
        [{"atoms":g, "elements":[str(elems[i]) for i in g]} for g in frags],
        ensure_ascii=False
    )
    return len(frags), frag_list, frag_json

# ---------- ΔE 计算帮助 ----------
def _find_col_idx(columns, regex):
    pat=re.compile(regex, re.IGNORECASE)
    for i,c in enumerate(columns):
        if pat.search(str(c)): return i
    return None

def compute_delta_e_from_table(arr, cols, policy="ekin_delta"):
    if arr is None or cols is None or arr.shape[0]<1: return np.nan
    i_ekin = _find_col_idx(cols, r"\bE\s*kin\b|\bKinetic\b")
    i_epot = _find_col_idx(cols, r"\bE\s*pot\b|\bPotential\b")
    i_etot = _find_col_idx(cols, r"\bE\s*tot\b|\bTotal\b")

    def safe_col(i): return arr[:, i] if i is not None else None

    ekin = safe_col(i_ekin)
    epot = safe_col(i_epot)
    etot = safe_col(i_etot)

    if policy == "ekin_final" and ekin is not None:
        return float(ekin[-1])
    if policy == "ekin_delta" and ekin is not None and ekin.size>=2:
        return float(ekin[-1] - ekin[0])

    if policy == "auto":
        if ekin is not None and ekin.size>=2: return float(ekin[-1]-ekin[0])
        if ekin is not None:                return float(ekin[-1])
        if epot is not None and epot.size>=2: return float(epot[0]-epot[-1])
        if etot is not None and etot.size>=2: return float(etot[-1]-etot[0])
        return np.nan

    if ekin is not None and ekin.size>=2: return float(ekin[-1]-ekin[0])
    if epot is not None and epot.size>=2: return float(epot[0]-epot[-1])
    if etot is not None and etot.size>=2: return float(etot[-1]-etot[0])
    return np.nan

# ---------- HDF5 写入 ----------
def ensure_group(h5: h5py.File, gpath: str) -> h5py.Group:
    cur=h5["/"]
    for p in [p for p in gpath.strip("/").split("/") if p]:
        cur=cur.require_group(p)
    return cur

def write_one(h5: h5py.File, state_name: str, traj_name: str, parsed: dict,
              traj_dir: Path, mode: int):
    g = ensure_group(h5, f"/{state_name}/{traj_name}")
    str_t = h5py.string_dtype("utf-8")

    for name,data in [("geometry", parsed["R_A"]), ("velocity", parsed["V_au"])]:
        if name in g: del g[name]
        g.create_dataset(name, data=data, compression="gzip", shuffle=True, fletcher32=True)

    g.attrs["nstep"] = int(parsed.get("nstep",0))
    try: g.attrs["mtime"] = int(os.path.getmtime(traj_dir/"output.dat"))
    except: g.attrs["mtime"] = int(time.time())
    g.attrs["dir_path"] = np.array(str(traj_dir), dtype=str_t)
    g.attrs["xyz_path"] = np.array(f"{state_name}/{traj_name}/output.xyz", dtype=str_t)

    els = parsed.get("elements")
    g.attrs["elements"] = np.asarray(els if els is not None else np.array([], dtype=str), dtype=str_t)
    g.attrs["atomic_masses"] = np.asarray(parsed.get("atomic_masses", []), dtype=np.float64)

    sname = state_name.lower()
    if   "singlet" in sname or re.match(r"^s[_\-]?\d+", sname): mult="Singlet"
    elif "triplet" in sname or re.match(r"^t[_\-]?\d+", sname): mult="Triplet"
    elif "doublet" in sname or "duplet" in sname:               mult="Doublet"
    else: mult="Unknown"
    g.attrs["state"] = np.array(state_name, dtype=str_t)
    g.attrs["multiplicity"] = np.array(mult, dtype=str_t)

    if parsed["R_A"].size and g.attrs["elements"].size:
        nfrag, frag_list, frag_json = connected_fragments(parsed["R_A"][-1], g.attrs["elements"])
    else:
        nfrag, frag_list, frag_json = 0, "", "[]"
    g.attrs["nfrag"]     = int(nfrag)
    g.attrs["frag_list"] = np.array(frag_list, dtype=str_t)
    g.attrs["frag"]      = np.array(frag_json, dtype=str_t)

    delta_e = np.nan
    if int(mode)==2:
        expec_path = guess_expec_path(traj_dir)
        if expec_path.exists():
            txt = expec_path.read_text(encoding="utf-8", errors="ignore")
            arr, cols = parse_expec_table(txt)
            if "expec_out_raw" in g: del g["expec_out_raw"]
            g.create_dataset("expec_out_raw", data=np.array(txt, dtype=str_t))
            if arr is not None:
                if "energy_table" in g: del g["energy_table"]
                dset = g.create_dataset("energy_table", data=arr, compression="gzip", shuffle=True, fletcher32=True)
                dset.attrs["columns"] = np.array(cols, dtype=str_t)
                policy = DELTA_E_POLICY if DELTA_E_POLICY in ("ekin_delta","ekin_final","auto") else "ekin_delta"
                delta_e = compute_delta_e_from_table(arr, cols, policy=policy)

    g.attrs["delta_energy"] = np.float64(delta_e)

# ---------- 枚举所有“状态/轨迹”并去重 ----------
def iter_all_trajs(roots):
    trajs=[]
    for root in roots:
        root=Path(root)
        if not root.exists(): continue
        state_dirs=set()
        for p in root.rglob("*"):
            if p.is_dir() and any(rx.match(p.name) for rx in RE_STATE):
                state_dirs.add(p)
        for sdir in sorted(state_dirs):
            for tdir in sorted(p for p in sdir.iterdir() if p.is_dir() and RE_TRAJ.match(p.name)):
                odat = tdir/"output.dat"
                if odat.exists():
                    trajs.append( (sdir.name, tdir.name, tdir.resolve()) )
    uniq={}
    for s,t,d in trajs:
        key=str(d)
        if key not in uniq: uniq[key]=(s,t,d)
    return list(uniq.values())


def update_progress(phase_name, current_dir, current, total):
    """在终端原地刷新双行进度条，不滚动屏幕"""
    bar_len = 40
    frac = current / total if total > 0 else 1
    filled = int(bar_len * frac)
    bar = '=' * filled + '-' * (bar_len - filled)
    
    # \033[K 用于清空当前行，防止旧的较长字符串残留
    sys.stdout.write(f"\033[K[{phase_name}] 当前处理: {current_dir}\n")
    sys.stdout.write(f"\033[K总进度: [{bar}] {current}/{total} ({frac*100:.1f}%)\n")
    # \033[2A 将光标上移两行，为下一次覆盖做准备
    sys.stdout.write("\033[2A")
    sys.stdout.flush()
    
    
# ---------- 主程序 ----------
def main():
    outp = Path(OUTPUT_PATH)
    mode = MODE if MODE is not None else (2 if outp.name.endswith(".2.hdf5") else 1)

    trajs = iter_all_trajs(ROOTS)
    total_trajs = len(trajs)
    print(f"[INFO] unique trajectories found: {total_trajs}")

    # ================= Phase 1: 预先执行外部程序生成数据 =================
    if mode == 2 and GENERATE_EXPEC:
        print(f"[INFO] Phase 1: Checking and generating output_data using {DATA_EXTRACTOR} ...\n")
        
        def pre_process(item):
            s, t, d = item
            expec_path = guess_expec_path(d)
            if not expec_path.exists():
                cmd = GEN_CMD.format(exe=DATA_EXTRACTOR)
                subprocess.run(cmd, cwd=str(d), shell=True, 
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return f"{s}/{t}" # 返回目录名用于进度条显示

        if N_WORKERS and N_WORKERS > 1:
            from concurrent.futures import ThreadPoolExecutor, as_completed
            with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
                futures = {ex.submit(pre_process, item): item for item in trajs}
                completed = 0
                for fu in as_completed(futures):
                    dir_name = fu.result()
                    completed += 1
                    update_progress("Phase 1", dir_name, completed, total_trajs)
        else:
            completed = 0
            for item in trajs:
                s, t, d = item
                pre_process(item)
                completed += 1
                update_progress("Phase 1", f"{s}/{t}", completed, total_trajs)
                
        sys.stdout.write("\n\n") # 循环结束后下移两行，防止结尾语覆盖进度条
        print("[INFO] Phase 1 Complete: All data directories checked/generated.")
    
    print("\n[INFO] Phase 2: Starting data extraction and writing to HDF5 ...\n")
    # =====================================================================

    outp.parent.mkdir(parents=True, exist_ok=True)
    str_t = h5py.string_dtype("utf-8")

    if N_WORKERS and N_WORKERS>1:
        from concurrent.futures import ThreadPoolExecutor, as_completed
        with h5py.File(outp, "w") as h5:
            h5.attrs["created_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
            h5.attrs["source_roots"] = np.array([str(r) for r in ROOTS], dtype=str_t)
            h5.attrs["units_geometry"] = "Angstrom"
            h5.attrs["units_velocity"] = "a.u."
            h5.attrs["mode"] = np.array(f"{mode}", dtype=str_t)

            n_ok = n_fail = completed = 0
            def work(item):
                s,t,d = item
                parsed = parse_output_dat(d/"output.dat")
                return (s,t,d,parsed)

            with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
                futures = [ex.submit(work, it) for it in trajs]
                for fu in as_completed(futures):
                    s,t,d,parsed = fu.result()
                    completed += 1
                    update_progress("Phase 2", f"{s}/{t}", completed, total_trajs)
                    
                    if not parsed.get("ok", False): 
                        n_fail+=1
                        continue
                    write_one(h5, s, t, parsed, d, mode)
                    n_ok+=1
                    
            sys.stdout.write("\n\n")
            print(f"[DONE] mode={mode} wrote {n_ok} trajectories to {outp} (failed: {n_fail})")
    else:
        with h5py.File(outp, "w") as h5:
            # 写入单线程属性同上... (为节省篇幅略去冗余 attrs 写入，直接补齐进度条逻辑)
            h5.attrs["created_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
            h5.attrs["source_roots"] = np.array([str(r) for r in ROOTS], dtype=str_t)
            h5.attrs["units_geometry"] = "Angstrom"
            h5.attrs["units_velocity"] = "a.u."
            h5.attrs["mode"] = np.array(f"{mode}", dtype=str_t)

            n_ok = n_fail = completed = 0
            for s,t,d in trajs:
                parsed = parse_output_dat(d/"output.dat")
                completed += 1
                update_progress("Phase 2", f"{s}/{t}", completed, total_trajs)
                
                if not parsed.get("ok", False): 
                    n_fail+=1
                    continue
                write_one(h5, s, t, parsed, d, mode)
                n_ok+=1
                
            sys.stdout.write("\n\n")
            print(f"[DONE] mode={mode} wrote {n_ok} trajectories to {outp} (failed: {n_fail})")

if __name__ == "__main__":
    main()