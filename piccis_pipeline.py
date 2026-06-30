#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#stdlib
import argparse
import csv
import itertools
import multiprocessing
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
#bioinformatics
from Bio import Phylo, SeqIO
#data / viz
import matplotlib
matplotlib.use("Agg")   # backend no interactivo guarda SVG sin abrir ventanas
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Polygon, Rectangle
from scipy.stats import mannwhitneyu
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
# geopandas / geopy / distinctipy
try:
    import geopandas as gpd
    from geopy.exc import GeocoderTimedOut
    from geopy.geocoders import Nominatim
    from distinctipy import distinctipy
    _GEO_OK = True
except ImportError:
    _GEO_OK = False
#################################################
#banner
##############################################
def banner(color: bool = None) -> str:
    if color is None:
        color = sys.stdout.isatty()

    letters = {
        'P': ["██████╗ ", "██╔══██╗", "██████╔╝", "██╔═══╝ ", "██║     ", "╚═╝     "],
        'I': ["██╗", "██║", "██║", "██║", "██║", "╚═╝"],
        'C': [" ██████╗", "██╔════╝", "██║     ", "██║     ", "╚██████╗", " ╚═════╝"],
        'S': ["███████╗", "██╔════╝", "███████╗", "╚════██║", "███████║", "╚══════╝"],
    }
    helix_l = ["o   o", " \\ / ", "  X  ", " / \\ ", "o   o", " \\ / "]
    helix_r = ["o   o", " \\ / ", "  X  ", " / \\ ", "o   o", " \\ / "]
    grad = [45, 39, 63, 99, 135, 171]   # cian → azul → morado → rosa

    def col(code, s):
        return f"\033[38;5;{code}m{s}\033[0m" if color else s

    title = ["" for _ in range(6)]
    for ch in "PICCIS":
        for i in range(6):
            title[i] += letters[ch][i] + " "

    lines = []
    for i in range(6):
        lines.append("  " + col(grad[i], helix_l[i]) + "   "
                     + col(grad[i], title[i]) + "  "
                     + col(grad[i], helix_r[i]))

    bar = col(135, "═" * 74)
    out = ["", bar] + lines + [
        "",
        col(141, "  Plasmid Identification, Clustering and Comparative Integrated Score"),
        bar, "",
    ]
    return "\n".join(out)


# ═══════════════════════════════════════════════════════════════════
#  COMANDOS DEL SISTEMA
# ═══════════════════════════════════════════════════════════════════

def run_cmd(cmd: list, cwd: Path = None, check: bool = True,
            fatal: bool = True) -> int:
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"\n[{ts}] ▶ {' '.join(str(c) for c in cmd)}")
    r = subprocess.run([str(c) for c in cmd], cwd=cwd)
    if check and r.returncode != 0:
        msg = f"✖  Error en: {' '.join(str(c) for c in cmd)}"
        if fatal:
            sys.exit(msg)
        else:
            raise RuntimeError(msg)
    return r.returncode

#ejecutables en path
def check_exe(*exes: str):
    missing = [e for e in exes if not shutil.which(e)]
    if missing:
        sys.exit(
            "✖  Executables not found: " + ", ".join(missing) +
            "\n   Activate environment: conda activate plasmidos_env"
        )

#verificar los entornos (con wich por que algunos no tienen version
def check_env(env_name: str, exe: str) -> bool:
    r = subprocess.run(
        ["conda", "run", "-n", env_name, "which", exe],
        capture_output=True
    )
    if r.returncode != 0:
        print(f"⚠  Environment '{env_name}' not found or '{exe}' not available.")
        print(f"   Instal with: mamba create -n {env_name} "
              f"-c bioconda -c conda-forge {env_name.replace('_env','')} -y")
        return False
    path = r.stdout.decode().strip()
    print(f"   ✔  {exe} → {path} ({env_name})")
    return True

#ejecutar dentro del entorno especifico
def conda_run(env_name: str, cmd: list, **kwargs):
    full_cmd = ["conda", "run", "--no-capture-output", "-n", env_name] + \
               [str(c) for c in cmd]
    return run_cmd(full_cmd, **kwargs)

#input interactivo
def ask(prompt: str, valid: list = None) -> str:
    """Interactive input with validation."""
    while True:
        val = input(prompt).strip()
        if valid is None or val in valid:
            return val
        print(f"Invalid option. Choose between: {valid}")

#piccis.config que tiene los directorios de BD
def leer_conf() -> dict:
    conf_path = Path(__file__).parent / "piccis.conf"
    conf: dict[str, str] = {}
    if not conf_path.exists():
        return conf
    with open(conf_path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                key, _, val = line.partition("=")
                val = val.strip()
                if val:
                    conf[key.strip()] = val
    return conf

#ruta de las BD
def pedir_db(conf: dict, key: str, prompt: str) -> Path:
    if key in conf and Path(conf[key]).exists():
        print(f"   [conf] {key} = {conf[key]}")
        return Path(conf[key]).resolve()
    ruta = input(prompt).strip()
    return Path(ruta).resolve()


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path

#longitud de fastas
def fasta_len(path: Path) -> int:
    return sum(len(l.strip()) for l in open(path) if not l.startswith(">"))

#figura como SVG
def guardar_fig(fig, nombre: str, out_dir: Path):
    ruta = out_dir / f"{nombre}.svg"
    fig.savefig(ruta, format="svg", bbox_inches="tight")
    print(f"   💾  {ruta.name}")


#Detección de cores disponibles sea por la entrada de --cores N en linea de comando, PICCIS_WORKERS=N  (variable de entorno) o cpu_count() - 1  (automático, mínimo 1)
def detectar_workers(cli_cores: int = None) -> int:
    total = multiprocessing.cpu_count()

    if cli_cores is not None:
        n = max(1, cli_cores)
        print(f"   [parallel] --cores {n} (CLI argument)")
        return n

    env_val = os.environ.get("PICCIS_WORKERS")
    if env_val and env_val.isdigit():
        n = max(1, int(env_val))
        print(f"   [parallel] PICCIS_WORKERS={n} (environment variable)")
        return n

    n = max(1, total - 1)
    print(f"   [parallel] {total} cores detected → using {n} workers (automatic)")
    return n


# Ejecutor paralelo

def ejecutar_en_paralelo(fn, trabajos: list, n_workers: int,
                          desc: str = "") -> list:

    resultados = []

    if n_workers == 1 or len(trabajos) <= 1:
        for t in trabajos:
            resultados.append(fn(t))
        return resultados

    label = f" [{desc}]" if desc else ""
    total = len(trabajos)
    print(f"\n   [parallel]{label} {total} tasks · {n_workers} workers")

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futuros = {pool.submit(fn, t): t for t in trabajos}
        completados = 0
        for fut in as_completed(futuros):
            completados += 1
            try:
                res = fut.result()
                resultados.append(res)
            except Exception as exc:
                trabajo = futuros[fut]
                print(f"   ✖  Error in {trabajo}: {exc}")
            ts = datetime.now().strftime("%H:%M:%S")
            print(f"   [{ts}] {label} {completados}/{total} completed")

    return resultados


# colores global (Niche)
NICHE_COLORS = {
    "Environment":  "darkmagenta",
    "Clinic":       "yellow",
    "Undetermined": "linen",
}


def _fmt_x(x, _):
    """Formats large numbers as K / M for axis ticks."""
    if x >= 1e6:  return f"{x*1e-6:.1f}M"
    if x >= 1e3:  return f"{x*1e-3:.0f}K"
    return f"{x:.0f}"


# ═══════════════════════════════════════════════════════════════════
#  MÓDULO 1 – ENTRADA FASTQ  →  SPAdes + MOB-recon
# ═══════════════════════════════════════════════════════════════════
#spades entorno propio por que colisiona con platon, eggnog-mapper y blast=2.12 en plasmidos_env
def _spades_single(args: tuple):
    library_type, fasta_or_pair, out_dir, threads_per_job, meta = args
    mkdir(out_dir)

    mode_flag = "--metaplasmid" if meta else "--plasmid"
    base = ["spades.py", mode_flag, "-t", str(threads_per_job), "-o", out_dir]

    if library_type == "iontorrent":
        cmd = ["spades.py", "--iontorrent", "-s", fasta_or_pair,
               "--plasmid", "-t", str(threads_per_job), "-o", out_dir]
    elif library_type == "single":
        cmd = base + ["-s", fasta_or_pair]
    elif library_type == "paired":
        r1, r2 = fasta_or_pair
        cmd = base + ["-1", r1, "-2", r2]
    elif library_type == "paired+pacbio":
        r1, r2, pb = fasta_or_pair
        cmd = base + ["-1", r1, "-2", r2, "--pacbio", pb]
    elif library_type == "paired+nanopore":
        r1, r2, np_ = fasta_or_pair
        cmd = base + ["-1", r1, "-2", r2, "--nanopore", np_]
    else:
        return

    conda_run("spades_env", cmd)

#funcion q ejecuta unycicler _indica que es interna
def _unicycler_single(args: tuple):

    library_type, fasta_or_pair, out_dir, threads_per_job = args
    mkdir(out_dir)

    base = ["unicycler", "--mode", "conservative",
            "-t", str(threads_per_job),
            "--spades_options", "-m 6",
            "-o", out_dir]

    if library_type in ("single", "iontorrent"):
        cmd = base + ["-s", fasta_or_pair]
    elif library_type == "paired":
        r1, r2 = fasta_or_pair
        cmd = base + ["-1", r1, "-2", r2]
    elif library_type == "paired+pacbio":
        r1, r2, pb = fasta_or_pair
        cmd = base + ["-1", r1, "-2", r2, "--existing_long_read_assembly", pb]
    elif library_type == "paired+nanopore":
        r1, r2, np_ = fasta_or_pair
        cmd = base + ["-1", r1, "-2", r2, "-l", np_]
    else:
        return

    conda_run("unicycler_env", cmd)

#recibe fastq calcula cuantos threads calcular por muestra
def run_unicycler(fastq_dir: Path, out_dir: Path, library_type: str,
                  n_workers: int = 1, total_threads: int = None):
    mkdir(out_dir)
    total_threads   = total_threads or multiprocessing.cpu_count()
    threads_per_job = max(1, total_threads // max(1, n_workers))

    exts   = ("*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz")
    all_fq = [f for ext in exts for f in fastq_dir.glob(ext)]
    trabajos = []

    if library_type == "iontorrent":
        for fq in all_fq:
            sample = fq.name.split(".fastq")[0].split(".fq")[0]
            if ya_corrido(out_dir, sample, "assembly.fasta"):
                print(f"   [Unicycler] {sample} already assembled, skipping.")
                continue
            trabajos.append(("iontorrent", fq,
                              mkdir(out_dir / sample), threads_per_job))

    elif library_type == "single":
        for fq in all_fq:
            sample = fq.name.split(".fastq")[0].split(".fq")[0]
            if ya_corrido(out_dir, sample, "assembly.fasta"):
                print(f"   [Unicycler] {sample} already assembled, skipping.")
                continue
            trabajos.append(("single", fq,
                              mkdir(out_dir / sample), threads_per_job))

    elif library_type in ("paired", "paired+pacbio", "paired+nanopore"):
        pairs = {}
        long_tag = {"paired+pacbio": "_pb", "paired+nanopore": "_np"}.get(library_type)
        for fq in all_fq:
            name = fq.name
            for tag in ("_R1", "_1"):
                if tag in name:
                    pairs.setdefault(name.split(tag)[0], {})["R1"] = fq
            for tag in ("_R2", "_2"):
                if tag in name:
                    pairs.setdefault(name.split(tag)[0], {})["R2"] = fq
            if long_tag and long_tag in name:
                pairs.setdefault(name.split(long_tag)[0], {})["LONG"] = fq

        for sample, reads in pairs.items():
            if "R1" not in reads or "R2" not in reads:
                print(f"⚠  Unicycler: incomplete pair {sample}, skipping.")
                continue
            if ya_corrido(out_dir, sample, "assembly.fasta"):
                print(f"   [Unicycler] {sample} already assembled, skipping.")
                continue
            if library_type == "paired":
                payload = (reads["R1"], reads["R2"])
            else:
                payload = (reads["R1"], reads["R2"],
                           reads.get("LONG", reads["R1"]))
            trabajos.append((library_type, payload,
                              mkdir(out_dir / sample), threads_per_job))
    else:
        print(f"⚠  Unicycler: not supported ({library_type}), skipping.")
        return

    if not trabajos:
        print("   ⚠  Unicycler: no new FASTQs to process.")
        return

    print(f"\n   [Unicycler] {len(trabajos)} samples · "
          f"{threads_per_job} threads/sample · conservative")
    ejecutar_en_paralelo(_unicycler_single, trabajos, n_workers, "Unicycler")
    print(f"✔  Unicycler → {out_dir}")

#correr mobrecon sobre fastas de ensamblados
def run_mob_recon_from_dir(assembly_dir: Path, out_dir: Path,
                            fasta_name: str = "contigs.fasta",
                            n_workers: int = 1):
    """
    Runs mob_recon on all assembled FASTAs in assembly_dir.
    Skips samples where MOB-recon already ran (chromosome.fasta exists).
    """
    mkdir(out_dir)
    trabajos = []
    for sample_dir in sorted(assembly_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        fasta = sample_dir / fasta_name
        if not fasta.exists():
            print(f"   ⚠  mob_recon: no {fasta_name} in "
                  f"{sample_dir.name}, skipping")
            continue
        # checkpoint: mob_recon genera chromosome.fasta al terminar
        mob_done = out_dir / sample_dir.name / "chromosome.fasta"
        if mob_done.exists() and mob_done.stat().st_size > 0:
            print(f"   [mob_recon] {sample_dir.name} already processed, skipping.")
            continue
        trabajos.append((fasta, mkdir(out_dir / sample_dir.name)))
    if not trabajos:
        print(f"   [mob_recon] All samples already processed.")
        return
    ejecutar_en_paralelo(_mob_recon_single, trabajos, n_workers, "MOB-recon")
    print(f"✔  MOB-recon → {out_dir}")

#abricate
def extraer_contigs_con_replicon(fuentes: list, home: Path,
                                  n_workers: int = 1) -> Path:
    out_dir  = mkdir(home / "unicycler-replicons")
    ab_dir   = mkdir(home / "unicycler-abricate")
    from Bio import SeqIO as _SeqIO
    import subprocess as _sp

    def _procesar_fasta(fasta: Path, sample: str):
        """Runs Abricate on a FASTA and extracts the replicon regions."""
        ab_sample = mkdir(ab_dir / sample)
        ab_tab    = ab_sample / "plasmidfinder.tab"

        # Checkpoint
        if ab_tab.exists() and ab_tab.stat().st_size > 100:
            pass  # ya corrió
        else:
            cmd = ["conda", "run", "--no-capture-output", "-n", "abricate_env",
                   "abricate", str(fasta),
                   "--db", "plasmidfinder", "--minid", "80", "--mincov", "80"]
            with open(ab_tab, "w") as fh:
                _sp.run(cmd, stdout=fh, check=False)

        try:
            df = pd.read_csv(ab_tab, sep="\t")
            if df.empty or "SEQUENCE" not in df.columns:
                return
        except Exception:
            return

        if "START" not in df.columns or "END" not in df.columns:
            return

        n_hits = len(df)
        if n_hits == 0:
            return

        print(f"   [Abricate-replicon] {sample}: {n_hits} hit(s) with replicon")
        contig_index = {rec.id: rec for rec in _SeqIO.parse(fasta, "fasta")}

        for _, hit in df.iterrows():
            seq_id = str(hit.get("SEQUENCE", "")).strip()
            start  = int(hit.get("START", 1)) - 1
            end    = int(hit.get("END", 0))
            gene   = str(hit.get("GENE", "replicon")).strip().replace("/", "_")

            if seq_id not in contig_index:
                continue

            rec    = contig_index[seq_id]
            region = rec[start:end]
            region.id          = f"{sample}__{seq_id}__{gene}"
            region.description = f"replicon={gene} coords={start+1}-{end}"

            out_fasta = out_dir / f"{region.id}.fasta"
            if not out_fasta.exists():
                with open(out_fasta, "w") as fh:
                    _SeqIO.write(region, fh, "fasta")
                print(f"   → extracted: {out_fasta.name} ({len(region)} bp)")

    # Procesar cada fuente
    for src in fuentes:
        if not src.exists():
            continue

        # Caso 1: directorio con subdirectorios de muestra (Unicycler)
        subdirs = [d for d in src.iterdir() if d.is_dir()]
        if subdirs:
            for sample_dir in sorted(subdirs):
                assembly = sample_dir / "assembly.fasta"
                if assembly.exists():
                    _procesar_fasta(assembly, sample_dir.name)
        else:
            # Caso 2: directorio plano con FASTAs individuales (GBK-FASTA-plasmid)
            for fasta in sorted(src.glob("*.fasta")):
                _procesar_fasta(fasta, fasta.stem)

    n = len(list(out_dir.glob("*.fasta")))
    if n > 0:
        print(f"✔  Abricate-replicon: {n} regions extracted → {out_dir}")
    else:
        print(f"   [Abricate-replicon] No new replicons detected")

    return out_dir

#descomprimir fastq a partir de .tar.gz / .tar.bz2 / .tgz / .zip
def descomprimir_fastqs(fastq_dir: Path, home: Path) -> Path:
    import tarfile, zipfile

    patrones = (list(fastq_dir.glob("*.tar.gz")) +
                list(fastq_dir.glob("*.tar.bz2")) +
                list(fastq_dir.glob("*.tgz")) +
                list(fastq_dir.glob("*.zip")))

    if not patrones:
        return fastq_dir   # nada que descomprimir, usar carpeta original

    dest_dir = mkdir(home / "FASTQ-descomprimidos")
    print(f"\n   [decompress] {len(patrones)} compressed file(s) → {dest_dir}")

    for archivo in patrones:
        try:
            if archivo.suffix == ".zip":
                with zipfile.ZipFile(archivo) as z:
                    for m in z.namelist():
                        if m.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
                            dest = dest_dir / Path(m).name
                            if not dest.exists():
                                print(f"   → {archivo.name}  →  {dest.name}")
                                z.extract(m, dest_dir)
                            else:
                                print(f"   → {dest.name} already exists, skipping")
            else:
                # Intentar como tar primero
                try:
                    with tarfile.open(archivo) as t:
                        for m in t.getnames():
                            if m.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
                                dest = dest_dir / Path(m).name
                                if not dest.exists():
                                    print(f"   → {archivo.name}  →  {dest.name}")
                                    obj = t.getmember(m)
                                    obj.name = Path(m).name
                                    t.extract(obj, dest_dir)
                                else:
                                    print(f"   → {dest.name} already exists, skipping")
                except tarfile.TarError:
                    # No es un tar real puede ser un .fastq.gz con extensión .tar.gz
                    import gzip
                    try:
                        with gzip.open(archivo, 'rb') as g:
                            g.read(10)  # test
                        # si gzip válido copiarlo como .fastq.gz
                        nuevo_nombre = archivo.name.replace(".tar.gz", ".fastq.gz").replace(".tgz", ".fastq.gz")
                        dest = dest_dir / nuevo_nombre
                        if not dest.exists():
                            print(f"   → {archivo.name} is not tar, copying as {nuevo_nombre}")
                            shutil.copy(archivo, dest)
                    except Exception:
                        print(f"   ⚠  {archivo.name} is not a valid tar or gzip — skipping.")
        except Exception as e:
            print(f"   ⚠  Could not decompress {archivo.name}: {e}")

    # Copiar también los .fastq.gz que ya estaban sueltos en fastq_dir
    for f in fastq_dir.glob("*.fastq.gz"):
        dest = dest_dir / f.name
        if not dest.exists():
            shutil.copy(f, dest)
    for f in fastq_dir.glob("*.fq.gz"):
        dest = dest_dir / f.name
        if not dest.exists():
            shutil.copy(f, dest)

    print(f"   [decompress] done → using {dest_dir}\n")
    return dest_dir

#
def ya_corrido(out_dir: Path, sample: str, archivo_esperado: str) -> bool:
    """
    Checks whether a tool already ran for a specific sample.
    Looks for out_dir/<sample>/<expected_file>.
    If it exists and is not empty, returns True (skip).

    Usage:
        if ya_corrido(spades_raw, "HST74", "contigs.fasta"):
            print("HST74 already has SPAdes contigs, skipping")
    """
    target = out_dir / sample / archivo_esperado
    if target.exists() and target.stat().st_size > 0:
        return True
    return False

#devolver la lista de nombres de las muestras
def samples_desde_fastqs(fastq_dir: Path) -> list:

    exts = ("*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz")
    muestras = set()
    for ext in exts:
        for fq in fastq_dir.glob(ext):
            name = fq.stem.replace(".fastq", "").replace(".fq", "")
            for sufijo in ("_R1", "_R2", "_1", "_2", "_pb", "_np"):
                if name.endswith(sufijo):
                    name = name[: -len(sufijo)]
                    break
            muestras.add(name)
    return sorted(muestras)

#correr spaled sobre fastq
def run_spades(fastq_dir: Path, out_dir: Path, library_type: str,
               n_workers: int = 1, total_threads: int = None,
               meta: bool = False):
    mkdir(out_dir)
    if total_threads is None:
        total_threads = multiprocessing.cpu_count()
    threads_per_job = max(1, total_threads // max(1, n_workers))

    exts   = ("*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz")
    all_fq = [f for ext in exts for f in fastq_dir.glob(ext)]
    trabajos = []

    if library_type == "iontorrent":
        for fq in all_fq:
            sample = fq.name.split(".fastq")[0].split(".fq")[0]
            if ya_corrido(out_dir, sample, "contigs.fasta"):
                print(f"   [SPAdes] {sample} already assembled, skipping.")
                continue
            trabajos.append(("iontorrent", fq,
                              mkdir(out_dir / sample), threads_per_job, meta))

    elif library_type == "single":
        for fq in all_fq:
            sample = fq.name.split(".fastq")[0].split(".fq")[0]
            if ya_corrido(out_dir, sample, "contigs.fasta"):
                print(f"   [SPAdes] {sample} already assembled, skipping.")
                continue
            trabajos.append(("single", fq,
                              mkdir(out_dir / sample), threads_per_job, meta))

    elif library_type in ("paired", "paired+pacbio", "paired+nanopore"):
        pairs = {}
        long_tag = {"paired+pacbio": "_pb", "paired+nanopore": "_np"}.get(library_type)
        for fq in all_fq:
            name = fq.name
            for tag in ("_R1", "_1"):
                if tag in name:
                    pairs.setdefault(name.split(tag)[0], {})["R1"] = fq
            for tag in ("_R2", "_2"):
                if tag in name:
                    pairs.setdefault(name.split(tag)[0], {})["R2"] = fq
            if long_tag and long_tag in name:
                pairs.setdefault(name.split(long_tag)[0], {})["LONG"] = fq

        for sample, reads in pairs.items():
            if "R1" not in reads or "R2" not in reads:
                print(f"⚠  Incomplete pair: {sample}, skipping.")
                continue
            if ya_corrido(out_dir, sample, "contigs.fasta"):
                print(f"   [SPAdes] {sample} already assembled, skipping.")
                continue
            if library_type == "paired":
                payload = (reads["R1"], reads["R2"])
            else:
                payload = (reads["R1"], reads["R2"],
                           reads.get("LONG", reads["R1"]))
            trabajos.append((library_type, payload,
                              mkdir(out_dir / sample), threads_per_job, meta))
    else:
        sys.exit(f"✖  Unknown library type: {library_type}")

    print(f"   [SPAdes] {len(trabajos)} samples · "
          f"{threads_per_job} threads/sample · "
          f"{'--metaplasmid' if meta else '--plasmid'}")
    ejecutar_en_paralelo(_spades_single, trabajos, n_workers, "SPAdes")

#agrupa contigs de spades por component_x y concatena con NN
def process_spades_output(spades_dir: Path, processed_dir: Path):

    mkdir(processed_dir)
    for contigs in spades_dir.glob("*/contigs.fasta"):
        sample = contigs.parent.name
        comps: dict[str, list[str]] = {}
        with open(contigs) as fh:
            comp_id, seq = None, []
            for line in fh:
                if line.startswith(">"):
                    if comp_id is not None:
                        comps.setdefault(comp_id, []).append("".join(seq))
                    hdr = line.strip()
                    comp_id = (hdr.split("component_")[-1].split()[0]
                               if "component_" in hdr else hdr[1:].split()[0])
                    seq = []
                else:
                    seq.append(line.strip())
            if comp_id is not None:
                comps.setdefault(comp_id, []).append("".join(seq))
        for cid, frags in comps.items():
            out = processed_dir / f"{sample}-plasmid_component_{cid}.fasta"
            with open(out, "w") as w:
                w.write(f">{sample}_component_{cid}\n")
                w.write("NNNN".join(frags) + "\n")
    print(f"✔  SPAdes processed → {processed_dir}")

#corre mob_recon sobre fastq
def _mob_recon_single(args: tuple):
    """Runs mob_recon on a FASTQ using mob_env."""
    fq, out_dir = args
    conda_run("mob_env",
              ["mob_recon", "--infile", fq, "--outdir", out_dir, "--force"])


#correr platon
def _platon_single(args: tuple):
    fasta, platon_db, sample_out, mode, threads_per_job = args
    # Verificar longitud mínima (Platon necesita ≥1000 bp para predecir ORFs)
    try:
        total_len = sum(len(r.seq) for r in SeqIO.parse(fasta, "fasta"))
        if total_len < 1000:
            print(f"   [Platon] {fasta.name} too short ({total_len} bp) — skipping.")
            return
    except Exception:
        pass
    try:
        conda_run("platon_env", [
            "platon", "--db", platon_db,
            "--output", sample_out,
            "--mode", mode,
            "--threads", str(threads_per_job),
            fasta,
        ], fatal=False)
    except RuntimeError as e:
        print(f"   ⚠  Platon failed for {fasta.name} — skipping. ({e})")

#correr platon en paralelo si ya existe saltea el proceso
def run_platon(fasta_dir: Path, home: Path, platon_db: Path,
               mode: str = "accuracy", threads: int = 4,
               n_workers: int = 1) -> Path:

    platon_out      = mkdir(home / "platon_out")
    threads_por_job = max(1, threads // max(1, n_workers))
    trabajos = []

    # Buscar FASTAs: tanto en el raíz (GBK input) como en subdirectorios (Unicycler)
    fastas = list(fasta_dir.glob("*.fasta"))
    for subdir in fasta_dir.iterdir():
        if subdir.is_dir():
            asm = subdir / "assembly.fasta"
            if asm.exists():
                fastas.append(asm)

    for fasta in sorted(set(fastas)):
        # Si es un assembly.fasta de Unicycler, usar el nombre del directorio padre
        if fasta.name == "assembly.fasta":
            sample_name = fasta.parent.name
        else:
            sample_name = fasta.stem
        sample_out = mkdir(platon_out / sample_name)
        # Checkpoint: platon genera <sample>.plasmid.fasta al terminar
        done = list(sample_out.glob("*.plasmid.fasta"))
        if done:
            print(f"   [Platon] {sample_name} already processed, skipping.")
            continue
        trabajos.append((fasta, platon_db, sample_out, mode, threads_por_job))
    if not trabajos:
        print(f"   [Platon] All samples already processed.")
    else:
        ejecutar_en_paralelo(_platon_single, trabajos, n_workers, "Platon")
    print(f"✔  Platon → {platon_out}")
    return platon_out


# ═══════════════════════════════════════════════════════════════════
#  MÓDULO 2 – ENTRADA GBK, multiFASTA  →  FASTAs individuales
# ═══════════════════════════════════════════════════════════════════

def _gbk_id(record) -> str:
    locus      = record.name.replace(" ", "_")
    defin      = record.description.replace(" ", "_")
    acc        = record.annotations.get("accessions", ["accession_unknown"])[0]
    ver        = record.annotations.get("sequence_version", "1")
    length_bp  = f"{len(record.seq)}_bp"
    bioproject = biosample = assembly = "unknown"
    for dbx in record.dbxrefs:
        if "BioProject" in dbx:
            bioproject = dbx.split(":")[1].strip().replace(" ", "_")
        elif "BioSample" in dbx:
            biosample  = dbx.split(":")[1].strip().replace(" ", "_")
        elif "Assembly" in dbx:
            assembly   = dbx.split(":")[1].strip().replace(" ", "_")
    return f"{locus}|{length_bp}|{defin}|{acc}|{ver}|{bioproject}|{biosample}|{assembly}"


def gbk_to_fastas(gbk_path: Path, out_dir: Path):
    mkdir(out_dir)
    prefix = _safe_name(gbk_path.stem)
    for i, rec in enumerate(SeqIO.parse(gbk_path, "gb"), 1):
        rec_id = _gbk_id(rec)
        rec.id = rec_id
        rec.description = ""
        fname = f"{prefix}_contig_{i}.fasta"
        with open(out_dir / fname, "w") as fout:
            SeqIO.write(rec, fout, "fasta")
    print(f"✔  Individual FASTAs → {out_dir}")


def split_multifasta(fasta_path: Path, out_dir: Path):
    mkdir(out_dir)
    prefix = _safe_name(fasta_path.stem)
    for i, rec in enumerate(SeqIO.parse(fasta_path, "fasta"), 1):
        fname = f"{prefix}_contig_{i}.fasta"
        with open(out_dir / fname, "w") as fout:
            SeqIO.write(rec, fout, "fasta")
    print(f"✔  Individual FASTAs → {out_dir}")


# ═══════════════════════════════════════════════════════════════════
#  MÓDULO 3 – BLAST por muestra,fastas-unicos/
# ═══════════════════════════════════════════════════════════════════


def _extraer_sample(fasta: Path) -> str:
    
    import re
    stem = fasta.stem

    # Esquema con doble guion bajo → tomar todo antes del primer __
    if "__" in stem:
        sample = stem.split("__")[0]
        # Eliminar sufijo _contig_N si existe (viene de gbk_to_fastas)
        sample = re.sub(r'_contig_\d+$', '', sample)
        return sample

    # Esquemas legacy con guion simple
    for sep in ("_spades_", "_mob_recon_", "_platon_", "_genomad_",
                "_abricate_", "_plasflow_"):
        if sep in stem:
            sample = stem.split(sep)[0]
            sample = re.sub(r'_contig_\d+$', '', sample)
            return sample

    for tag in ("-plasmid_component_", "-platon-plasmid", "-plasmid"):
        if tag in stem:
            return stem.split(tag)[0]

    return stem


def _safe_name(name: str) -> str:
   
    import re
    # Reemplazar puntos por guiones (excepto la extensión .fasta)
    name = re.sub(r'\.(?!fasta$|fa$|fna$|gbk$|gff3?$)', '-', name)
    # Reemplazar espacios
    name = name.replace(" ", "_")
    return name


def _blast_par(q: Path, s: Path, blast_dir: Path,
               pident_thr: float, cov_thr: float):

    bout = blast_dir / f"{q.stem}_vs_{s.stem}.tsv"
    run_cmd(
        ["blastn", "-query", q, "-subject", s,
         "-outfmt", "6 qseqid sseqid pident length qlen slen",
         "-out", bout],
        check=False,
    )
    if not (bout.exists() and bout.stat().st_size > 0):
        return False, 0.0, 0.0

    best_pid = best_cov = 0.0
    with open(bout) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                pid    = float(parts[2])
                length = float(parts[3])
                qlen   = float(parts[4])
                slen   = float(parts[5])
            except ValueError:
                continue
            cov_q = length / qlen * 100 if qlen else 0
            cov_s = length / slen * 100 if slen else 0
            cov   = max(cov_q, cov_s)
            if pid >= pident_thr and cov >= cov_thr:
                return True, pid, cov
    return False, best_pid, best_cov


def blast_por_muestra(source_dirs: list, home: Path,
                      pident_thr: float = 95.0,
                      cov_thr:    float = 95.0) -> Path:

    unique_dir   = mkdir(home / "fastas-unicos")
    blast_dir    = mkdir(unique_dir / "blast")
    resumen_path = home / "resultados-blast-resumido.tsv"

    print(f"\n   [BLAST] Thresholds: identity ≥ {pident_thr}%  "
          f"| coverage ≥ {cov_thr}%")
    print( "   [BLAST] Comparison ONLY within each biological sample")

    por_sample: dict[str, list[Path]] = {}
    origen:     dict[Path, str]       = {}

    for d in source_dirs:
        d = Path(d)
        if not d.exists():
            continue
        for f in d.glob("*.fasta"):
            sid = _extraer_sample(f)
            por_sample.setdefault(sid, []).append(f)
            origen[f] = d.name

    if not por_sample:
        print("⚠  No FASTA files to compare.")
        return unique_dir

    total_fastas  = sum(len(v) for v in por_sample.values())
    total_pares   = sum(len(v)*(len(v)-1)//2 for v in por_sample.values())
    print(f"   [BLAST] {len(por_sample)} samples · "
          f"{total_fastas} FASTAs · {total_pares} intra-sample pairs")

    rows: list = []
    n_total_descartados = 0

    for sid, fastas in sorted(por_sample.items()):
        if len(fastas) == 1:
            # Sample con una sola herramienta: pasa directamente
            f = fastas[0]
            dest = unique_dir / f.name
            if not dest.exists():
                shutil.copy(f, dest)
            continue

        # Ordenar por longitud desc para que q >= s siempre
        fastas.sort(key=fasta_len, reverse=True)
        descartadas: set[Path] = set()

        for i, q in enumerate(fastas):
            for s in fastas[i + 1:]:
                if s in descartadas:
                    continue

                match, best_pid, best_cov = _blast_par(
                    q, s, blast_dir, pident_thr, cov_thr
                )

                if match:
                    # Conservar q (más largo), descartar s
                    descartadas.add(s)
                    rows.append([
                        sid,
                        q.name, s.name,
                        f"{fasta_len(q)}bp", f"{fasta_len(s)}bp",
                        f"{best_pid:.1f}%", f"{best_cov:.1f}%",
                        origen.get(q, "?"), origen.get(s, "?"),
                        "keep-longest", "discard",
                    ])
                else:
                    rows.append([
                        sid,
                        q.name, s.name,
                        f"{fasta_len(q)}bp", f"{fasta_len(s)}bp",
                        f"{best_pid:.1f}%", f"{best_cov:.1f}%",
                        origen.get(q, "?"), origen.get(s, "?"),
                        "keep", "keep",
                    ])

        # Copiar los supervivientes a fastas-unicos/
        for f in fastas:
            if f not in descartadas:
                dest = unique_dir / f.name
                if not dest.exists():
                    shutil.copy(f, dest)

        n_total_descartados += len(descartadas)

    #resumen
    with open(resumen_path, "w", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow([
            "sample",
            "query", "subject",
            "len_query", "len_subject",
            "pident_best", "cov_best",
            "origen_query", "origen_subject",
            "accion_query", "accion_subject",
        ])
        w.writerows(rows)

    n_unicos = len(list(unique_dir.glob("*.fasta")))
    print(f"\n✔  {n_unicos} FASTAs in fastas-unicos/")
    print(f"   {n_total_descartados} discarded by intra-sample deduplication")
    print(f"✔  Summary → {resumen_path}")
    return unique_dir


# Alias para compatibilidad con llamadas existentes en el pipeline
blast_all_vs_all = blast_por_muestra


# ═══════════════════════════════════════════════════════════════════
#  MÓDULO 4 – ANÁLISIS SOBRE fastas-unicos/
# ═══════════════════════════════════════════════════════════════════


def run_genomad(fasta: Path, out_dir: Path, genomad_db: Path,
                sample_name: str = None) -> Path:
    name       = sample_name or fasta.stem
    sample_out = mkdir(out_dir / name)
    # Checkpoint
    if list(sample_out.glob("*_plasmid_summary.tsv")):
        print(f"   [geNomad] {name} already processed, skipping.")
        return sample_out
    try:
        conda_run("genomad_env",
                  ["genomad", "end-to-end",
                   "--conservative",
                   "--cleanup",
                   "--threads", "1",
                   "--splits", "8",
                   fasta, sample_out, genomad_db],
                  fatal=False)
        print(f"✔  geNomad → {sample_out}")
    except (RuntimeError, SystemExit) as e:
        print(f"   ⚠  geNomad failed on {name}: {e}")
    return sample_out


def run_genomad_dir(fasta_dir: Path, home: Path, genomad_db: Path,
                    n_workers: int = 1) -> Path:
    genomad_out = mkdir(home / "genomad_out")
    tareas = []

    subdirs = [d for d in fasta_dir.iterdir() if d.is_dir()]
    if subdirs:
        for subdir in sorted(subdirs):
            for fn in ("assembly.fasta", "assembly.fa"):
                assembly = subdir / fn
                if assembly.exists():
                    tareas.append((assembly, subdir.name))
                    break
    else:
        for fasta in sorted(fasta_dir.glob("*.fasta")):
            tareas.append((fasta, fasta.stem))

    if not tareas:
        print(f"   ⚠  geNomad: no FASTAs found in {fasta_dir}")
        return genomad_out

    print(f"\n   [geNomad] {len(tareas)} samples...")
    for fasta, name in tareas:
        run_genomad(fasta, genomad_out, genomad_db, sample_name=name)

    return genomad_out


def collect_genomad_plasmids(genomad_out: Path, todos_plasmidos: Path,
                              sample_prefix: str = "") -> int:

    from Bio import SeqIO as _SeqIO
    n_extraidos = 0

    for sample_dir in sorted(genomad_out.iterdir()):
        if not sample_dir.is_dir():
            continue
        # geNomad genera los summaries en <sample>/<prefix>_summary/
        summaries = list(sample_dir.rglob("*_plasmid_summary.tsv"))
        if not summaries:
            continue

        import re as _re
        prefix_raw = sample_prefix or sample_dir.name
        # Quitar _contig_N si existe (viene de gbk_to_fastas)
        prefix = _re.sub(r'_contig_\d+$', '', prefix_raw)

        for summary in summaries:
            try:
                # geNomad genera el FASTA de plásmidos junto al summary
                plasmid_fna = summary.parent / summary.name.replace(
                    "_plasmid_summary.tsv", "_plasmid.fna")
                if not plasmid_fna.exists():
                    candidates = list(summary.parent.glob("*_plasmid.fna"))
                    plasmid_fna = candidates[0] if candidates else None
                if plasmid_fna is None or not plasmid_fna.exists():
                    continue

                for rec in _SeqIO.parse(plasmid_fna, "fasta"):
                    # Filtrar secuencias muy grandes (>500 kb) — probablemente
                    # son contigs cromosómicos que geNomad clasificó como plásmidos
                    if len(rec.seq) > 500000:
                        print(f"   [geNomad] {rec.id} ({len(rec.seq):,} bp) exceeds 500 kb → discarded")
                        continue
                    safe_id = _safe_name(rec.id)
                    dest = todos_plasmidos / f"{prefix}__genomad__{safe_id}.fasta"
                    if not dest.exists():
                        with open(dest, "w") as fh:
                            _SeqIO.write(rec, fh, "fasta")
                        n_extraidos += 1
            except Exception as e:
                print(f"   ⚠  geNomad collect error: {e}")

    if n_extraidos:
        print(f"✔  geNomad: {n_extraidos} plasmids → todos_plasmidos/")
    return n_extraidos


def run_mob_typer(fu: Path, home: Path):

    out_dir = mkdir(home / "mobtyper_out")
    fastas  = list(fu.glob("*.fasta")) + list(fu.glob("*.fa")) + \
              list(fu.glob("*.fna"))

    if not fastas:
        print(f"   ⚠  MOB-typer: no FASTAs found in {fu}")
        return

    # Deduplicar por nombre normalizado (guion → guion bajo)
    seen_norm = {}
    for fasta in sorted(fastas):
        norm = fasta.stem.replace("-", "_")
        if norm not in seen_norm:
            seen_norm[norm] = fasta
        else:
            # Preferir el que tiene guion bajo (esquema nuevo __)
            if "__" in fasta.stem:
                seen_norm[norm] = fasta

    for norm, fasta in sorted(seen_norm.items()):
        out_file = out_dir / f"{norm}_mobtyper.txt"
        if out_file.exists() and out_file.stat().st_size > 0:
            print(f"   [mob_typer] {norm} already processed, skipping.")
            continue
        try:
            conda_run("mob_env",
                      ["mob_typer", "-i", fasta, "-o", out_file],
                      fatal=False)
        except RuntimeError as e:
            print(f"   ⚠  mob_typer failed on {fasta.stem}: {e}")

    print(f"✔  MOB-typer → {out_dir}")


def _abricate_single(args: tuple):

    fu, db, outfile = args
    fastas = sorted([f for f in fu.iterdir()
                     if f.is_file() and f.suffix in (".fasta", ".fa", ".fna")])
    if not fastas:
        return
    import subprocess
    cmd = ["conda", "run", "--no-capture-output", "-n", "abricate_env",
           "abricate", "--db", db] + [str(f) for f in fastas]
    with open(outfile, "w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)


def run_abricate(fu: Path, home: Path, n_workers: int = 1):

    abdir    = mkdir(home / "Abricate")
    trabajos = [
        (fu, db, abdir / f"{db}.tab")
        for db in ["resfinder", "card", "vfdb", "plasmidfinder"]
    ]
    ejecutar_en_paralelo(_abricate_single, trabajos,
                          min(n_workers, 4), "Abricate")
    print(f"✔  Abricate → {abdir}")


def _downstream_independiente(args: tuple):

    herramienta, fu, home = args
    if herramienta == "mob_typer":
        run_mob_typer(fu, home)


def run_downstream_fase1_paralelo(fu: Path, home: Path,
                                   n_workers: int = 1):

    trabajos = [("mob_typer", fu, home)]
    print(f"\n   [Downstream phase 1] MOB-typer ({n_workers} workers)")
    ejecutar_en_paralelo(_downstream_independiente, trabajos,
                          1, "Downstream-F1")


def _bakta_single(args: tuple):

    fasta, bakta_db, sample_out = args
    # Checkpoint: si ya tiene .gff3 no vuelve a correr
    if list(sample_out.glob("*.gff3")) or list(sample_out.glob("*.gff")):
        print(f"   [Bakta] {fasta.stem} already annotated, skipping.")
        return
    # Si el directorio existe pero está vacío/incompleto, borrarlo
    if sample_out.exists():
        shutil.rmtree(sample_out)
    sample_out.mkdir(parents=True, exist_ok=True)
    conda_run("bakta_env",
              ["bakta", "--db", bakta_db,
               "--output", sample_out,
               "--prefix", fasta.stem,
               "--force",
               fasta])


def run_bakta(fu: Path, home: Path, bakta_db: Path,
              n_workers: int = 1) -> Path:

    bakta_out = mkdir(home / "bakta_out")
    fastas    = sorted(fu.glob("*.fasta"))

    # Hash de los IDs actuales de fastas-unicos/
    import hashlib
    ids_actual = "|".join(f.stem for f in fastas)
    hash_actual = hashlib.md5(ids_actual.encode()).hexdigest()[:8]
    hash_file   = bakta_out / ".input_hash"

    # Si el hash cambió, borrar bakta_out y empezar de nuevo
    if hash_file.exists():
        hash_previo = hash_file.read_text().strip()
        if hash_previo != hash_actual:
            print(f"   [Bakta] input FASTAs changed → re-annotating everything")
            shutil.rmtree(bakta_out)
            bakta_out.mkdir(parents=True, exist_ok=True)

    # Guardar hash actual
    hash_file.write_text(hash_actual)

    trabajos = [
        (fasta, bakta_db, bakta_out / fasta.stem)
        for fasta in fastas
    ]
    ejecutar_en_paralelo(_bakta_single, trabajos, n_workers, "Bakta")
    print(f"✔  Bakta → {bakta_out}")
    return bakta_out


def _sanitizar_gff3(gff_path: Path, fasta_path: Path = None) -> Path:

    import tempfile

    seq_index = {}
    if fasta_path and fasta_path.exists():
        from Bio import SeqIO as _SeqIO
        seq_index = {r.id: r.seq for r in _SeqIO.parse(fasta_path, "fasta")}

    lines_out = []
    genes_removidos = 0

    with open(gff_path) as fh:
        lines = fh.readlines()

    fasta_section = False
    for line in lines:
        if line.startswith("##FASTA"):
            fasta_section = True
        if fasta_section:
            lines_out.append(line)
            continue
        if line.startswith("#") or not line.strip():
            lines_out.append(line)
            continue

        parts = line.split("\t")
        if len(parts) < 9:
            lines_out.append(line)
            continue

        feature_type = parts[2]
        if feature_type not in ("CDS", "gene"):
            lines_out.append(line)
            continue

        if seq_index and feature_type == "CDS":
            seqid  = parts[0]
            start  = int(parts[3]) - 1
            end    = int(parts[4])
            strand = parts[6]
            if seqid in seq_index:
                seq = seq_index[seqid][start:end]
                if strand == "-":
                    seq = seq.reverse_complement()
                # Traducir y buscar stop interno
                aa = seq.translate()
                aa_str = str(aa)
                # Si tiene stop antes del final → gen inválido
                if "*" in aa_str[:-1]:
                    genes_removidos += 1
                    continue

        lines_out.append(line)

    if genes_removidos == 0:
        return gff_path  # sin cambios, usar original

    tmp = Path(tempfile.mktemp(suffix=".gff3", prefix="piccis_"))
    tmp.write_text("".join(lines_out))
    print(f"   [Panaroo] {gff_path.parent.name}: {genes_removidos} gene(s) with internal stop removed")
    return tmp


def run_panaroo(bakta_out: Path, home: Path, n_workers: int = 1) -> Path:

    pan_out = mkdir(home / "panaroo_out")

    todos_gffs = list(bakta_out.rglob("*.gff3")) + list(bakta_out.rglob("*.gff"))
    if not todos_gffs:
        print("⚠  No Bakta GFF; skipping Panaroo.")
        return pan_out

    gffs_validos = []
    for gff in todos_gffs:
        try:
            content = gff.read_text(errors="ignore")
            if "\tCDS\t" in content or "\tgene\t" in content:
                gffs_validos.append(gff)
            else:
                print(f"   [Panaroo] {gff.parent.name}: no CDS → excluded "
                      f"(secuencia demasiado corta)")
        except Exception:
            pass

    if len(gffs_validos) < 2:
        print(f"⚠  Panaroo requires at least 2 GFF3 with genes "
              f"({len(gffs_validos)} disponibles). Se omite.")
        return pan_out

    # Checkpoint
    import hashlib
    ids_actual  = "|".join(sorted(g.stem for g in gffs_validos))
    hash_actual = hashlib.md5(ids_actual.encode()).hexdigest()[:8]
    hash_file   = pan_out / ".input_hash"

    csv_out = pan_out / "gene_presence_absence_roary.csv"
    if csv_out.exists() and hash_file.exists():
        if hash_file.read_text().strip() == hash_actual:
            print(f"   [Panaroo] already run (no changes), skipping.")
            return pan_out
        print(f"   [Panaroo] input GFF3 changed → recomputing pangenome")
        shutil.rmtree(pan_out)
        pan_out.mkdir(parents=True, exist_ok=True)
    elif csv_out.exists() and not hash_file.exists():
        print(f"   [Panaroo] output without hash record → recomputing")
        shutil.rmtree(pan_out)
        pan_out.mkdir(parents=True, exist_ok=True)

    gffs_sanitizados = []
    for gff in gffs_validos:

        fasta = gff.parent / gff.name.replace(".gff3", ".fna").replace(".gff", ".fna")
        if not fasta.exists():
            fasta = next(gff.parent.glob("*.fna"), None)
        gffs_sanitizados.append(_sanitizar_gff3(gff, fasta))

    print(f"   [Panaroo] {len(gffs_sanitizados)} valid GFF3 out of {len(todos_gffs)} total")
    conda_run("panaroo_env",
              ["panaroo", "-i", *gffs_sanitizados, "-o", pan_out,
               "--clean-mode", "sensitive",
               "-t", str(max(1, n_workers))])
    hash_file.write_text(hash_actual)
    print(f"✔  Panaroo → {pan_out}")
    return pan_out

###eggnog mapper que corre remoto si no esta la base de datos descargada
def _eggnog_single(args: tuple):
    faa, outp, eggnog_db = args
    cmd = ["emapper.py", "-i", faa, "--output", outp, "--cpu", "1",
           "--override", "--sensmode", "fast"]
    if eggnog_db:
        cmd += ["--data_dir", eggnog_db]
    else:
        cmd += ["--server-mode"]
    conda_run("eggnog_env", cmd)


def run_eggnog(bakta_out: Path, home: Path, n_workers: int = 1,
               eggnog_db: Path = None):

    if eggnog_db is None:
        print("   [EggNOG] No DB configured — EggNOG skipped.")
        print("   To enable it, run install_databases.sh and choose option 1 or 2.")
        return

    egg_dir  = mkdir(home / "EggNOGmapper")

    stems_validos = {
        faa.stem for faa in bakta_out.rglob("*.faa")
        if "hypotheticals" not in faa.name
    }

    
    n_obsoletos = 0
    for f in list(egg_dir.glob("*.emapper.*")):
        stem_real = f.name.split(".emapper.")[0]
        if stem_real not in stems_validos:
            f.unlink()
            n_obsoletos += 1
    if n_obsoletos:
        print(f"   [sync] {n_obsoletos} stale EggNOG file(s) "
              f"removed (no longer present in bakta_out/)")

    trabajos = []
    for faa in bakta_out.rglob("*.faa"):
        if "hypotheticals" in faa.name:
            continue
        # Saltear FASTAs vacíos (secuencias muy cortas sin ORFs anotados)
        if faa.stat().st_size < 50:
            print(f"   [EggNOG] {faa.stem}: empty .faa, skipping.")
            continue
        outp = egg_dir / faa.stem
        if (egg_dir / f"{faa.stem}.emapper.annotations").exists():
            print(f"   [EggNOG] {faa.stem} already annotated, skipping.")
            continue
        trabajos.append((faa, outp, eggnog_db))
    if not trabajos:
        print(f"   [EggNOG] All samples already annotated.")
        return
    ejecutar_en_paralelo(_eggnog_single, trabajos, n_workers, "EggNOG")
    print(f"✔  EggNOG-mapper → {egg_dir}")


def run_tani(fu: Path, home: Path, n_threads: int = None,
             identity: float = 0.7, coverage: float = 0.7,
             bootstraps: int = 100):

    tani_repo = None
    for cand in [
        Path(__file__).parent / "tANI_tool",
        Path.home() / "Documentos" / "PICCIS" / "tANI_tool",
        Path.home() / "tANI_tool",
        home / "tANI_tool",
        fu.parent / "tANI_tool",
    ]:
        if cand.exists() and (cand / "tANI_tool.pl").exists():
            tani_repo = cand
            break

    if tani_repo is None:
        tani_repo = Path(__file__).parent / "tANI_tool"
        print(f"   [tANI] Cloning repository...")
        run_cmd(["git", "clone",
                 "https://github.com/sophiagosselin/tANI_tool.git",
                 tani_repo])

    tani_pl = tani_repo / "tANI_tool.pl"

    if not tani_pl.exists():
        run_cmd(["git", "clone",
                 "https://github.com/sophiagosselin/tANI_tool.git",
                 tani_repo])

    tani_local = fu / "tANI_tool.pl"
    if not tani_local.exists():
        shutil.copy(tani_pl, tani_local)

    # Copiar buildtree_w_support.R a fastas-unicos/ para que el script R
    # corra en el mismo directorio que la matriz tANI_original.matrix
    r_script_src = tani_repo / "buildtree_w_support.R"
    r_script_local = fu / "buildtree_w_support.R"
    if r_script_src.exists() and not r_script_local.exists():
        shutil.copy(r_script_src, r_script_local)

    if n_threads is None:
        n_threads = max(1, multiprocessing.cpu_count() - 1)

    # Checkpoint: comparar la cantidad EXACTA de líneas esperadas
    # (antes había un margen de tolerancia de 1.5x que no detectaba
    # diferencias chicas, ej. 21 vs 15 muestras)
    tani_matrix = fu / "outputs" / "tANI" / "tANI_original.matrix"

    def _limpiar_checkpoints_tani_pl():

        for nombre in ["setup.log", "original_calculations.log",
                       "blast_database.log"]:
            f = fu / nombre
            if f.exists():
                f.unlink()
        shutil.rmtree(fu / "intermediates", ignore_errors=True)
        shutil.rmtree(fu / "blast", ignore_errors=True)

    if tani_matrix.exists():
        n_lines   = sum(1 for _ in open(tani_matrix))
        n_samples = len(list(fu.glob("*.fasta")))
        expected  = n_samples + 1  # header + filas
        if n_lines != expected:
            print(f"   [tANI] Matrix does not match current fastas-unicos/ "
                  f"({n_lines} lines, exactly {expected} expected) — regenerating.")
            shutil.rmtree(fu / "outputs" / "tANI", ignore_errors=True)
            _limpiar_checkpoints_tani_pl()
        else:
            print(f"   [tANI] already run ({tani_matrix}), skipping.")
    elif (fu / "setup.log").exists():
        print(f"   [tANI] stale setup.log found without matrix — cleaning up.")
        _limpiar_checkpoints_tani_pl()

    if not tani_matrix.exists():
        run_cmd([
            "perl", tani_local,
            "-id", str(identity),
            "-cv", str(coverage),
            "-bt", str(bootstraps),
            "-t",  str(n_threads),
            "-v",  "1",
        ], cwd=fu)
    else:
        pass  # usar matriz existente

    # Copiar buildtree_w_support.R al directorio donde está la matriz
    # tANI genera los resultados en outputs/tANI/
    tani_output_dir = fu / "outputs" / "tANI"
    r_script_src = tani_repo / "buildtree_w_support.R"
    if not r_script_src.exists():
        r_script_src = fu / "buildtree_w_support.R"

    if r_script_src.exists() and tani_output_dir.exists():
        r_script = tani_output_dir / "buildtree_w_support.R"
        if not r_script.exists():
            shutil.copy(r_script_src, r_script)

        # Checkpoint: si el árbol ya existe, no regenerar
        best_tree = tani_output_dir / "BestTree_wSupport.tre"
        if best_tree.exists():
            print(f"   [tANI] tree already generated → {best_tree}")
        else:
            check = subprocess.run(
                ["conda", "run", "-n", "tani_env", "which", "Rscript"],
                capture_output=True
            )
            if check.returncode == 0:
                conda_run("tani_env", ["Rscript", r_script], cwd=tani_output_dir)
            elif shutil.which("Rscript"):
                run_cmd(["Rscript", r_script], cwd=tani_output_dir)
            else:
                print("⚠  Rscript not found.")
    else:
        print("⚠  buildtree_w_support.R or outputs/tANI/ not found.")


# ═══════════════════════════════════════════════════════════════════
#  MÓDULO 5 – TABLA UNIFICADA
# ═══════════════════════════════════════════════════════════════════

def _load_mob_typer(home: Path) -> pd.DataFrame:

    mob_dir = home / "mobtyper_out"
    if not mob_dir.exists():
        return pd.DataFrame()

    frames = []
    seen_ids = set()

    for f in sorted(mob_dir.glob("*.txt")):
        try:
            df = pd.read_csv(f, sep="\t")
            if df.empty:
                continue
            if "sample_id" in df.columns:
                df = df.rename(columns={"sample_id": "ID"})
            if "ID" not in df.columns:
                continue
            # Normalizar IDs: guion → guion bajo para deduplicar
            df["ID"] = df["ID"].astype(str).str.strip()
            df["_id_norm"] = df["ID"].str.replace("-", "_")
            for _, row in df.iterrows():
                norm = row["_id_norm"]
                if norm not in seen_ids:
                    seen_ids.add(norm)
                    frames.append(df[df["_id_norm"] == norm].drop(columns=["_id_norm"]))
        except Exception:
            pass

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def _load_plasmidfinder(home: Path) -> pd.DataFrame:
    tsv = home / "plasmidfinder_out" / "results_tab.tsv"
    if not tsv.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(tsv, sep="\t")
        if df.empty:
            return pd.DataFrame()
        if "Query ID" in df.columns:
            df = df.rename(columns={"Query ID": "ID"})
        return df
    except (pd.errors.EmptyDataError, Exception):
        return pd.DataFrame()


def _load_abricate(home: Path) -> pd.DataFrame:

    frames = []
    for db in ["resfinder", "card", "vfdb", "plasmidfinder"]:
        f = home / "Abricate" / f"{db}.tab"
        if f.exists():
            try:
                tmp = pd.read_csv(f, sep="\t")
                if not tmp.empty:
                    tmp["DB"] = db
                    frames.append(tmp)
            except (pd.errors.EmptyDataError, Exception):
                pass
    if not frames:
        return pd.DataFrame()

    df = pd.concat(frames, ignore_index=True)
    if "#FILE" not in df.columns:
        return pd.DataFrame()

    df["ID"] = df["#FILE"].apply(lambda x: Path(x).stem)

    # Agregar por ID: un resumen por plásmido, no una fila por hit
    rows = []
    for _id, grp in df.groupby("ID"):
        genes_resf = sorted(grp.loc[grp["DB"]=="resfinder", "GENE"].dropna().unique())
        genes_card = sorted(grp.loc[grp["DB"]=="card",      "GENE"].dropna().unique())
        genes_vfdb = sorted(grp.loc[grp["DB"]=="vfdb",       "GENE"].dropna().unique())
        rows.append({
            "ID":                     _id,
            "n_resistance_genes":     len(genes_resf) + len(genes_card),
            "resistance_genes":       ";".join(genes_resf + genes_card) or "-",
            "n_virulence_genes":      len(genes_vfdb),
            "virulence_genes":        ";".join(genes_vfdb) or "-",
        })
    return pd.DataFrame(rows)


def _load_genomad(home: Path) -> pd.DataFrame:
    """
    Reads geNomad results (*_plasmid_summary.tsv).
    Returns a DataFrame with ID and genomad_score.
    """
    genomad_out = home / "genomad_out"
    if not genomad_out.exists():
        return pd.DataFrame()
    rows = []
    for sample_dir in sorted(genomad_out.iterdir()):
        if not sample_dir.is_dir():
            continue
        for summary in sample_dir.glob("*_plasmid_summary.tsv"):
            try:
                df = pd.read_csv(summary, sep="\t")
                if "seq_name" not in df.columns:
                    continue
                for _, row in df.iterrows():
                    rows.append({
                        "ID":            str(row["seq_name"]).strip(),
                        "genomad_score": row.get("virus_score", row.get("plasmid_score", None)),
                        "genomad_topology": row.get("topology", "-"),
                    })
            except Exception:
                pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def _load_platon(home: Path) -> pd.DataFrame:
    """
    Reads all <sample>.tsv produced by Platon and concatenates them.
    Key columns: contig, length, score, rds, circularity, plasmid_hits,
    chromosome_hits, replication, mobilization, conjugation, orf_types.
    Adds column 'ID' = stem of the sample directory (= FASTA name).
    """
    platon_out = home / "platon_out"
    frames = []
    for tsv in platon_out.rglob("*.tsv"):
        try:
            tmp = pd.read_csv(tsv, sep="\t")
            tmp["ID"] = tsv.parent.name          # nombre de la muestra
            frames.append(tmp)
        except Exception:
            pass
    if frames:
        return pd.concat(frames, ignore_index=True)
    return pd.DataFrame()


PANAROO_METADATA_COLS = {
    "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
    "Avg sequences per isolate", "Genome Fragment", "Order within Fragment",
    "Accessory Fragment", "Accessory Order with Fragment", "QC",
    "Min group size nuc", "Max group size nuc", "Avg group size nuc",
}


def _load_panaroo(home: Path) -> pd.DataFrame:

    pan_dir   = home / "panaroo_out"
    rtab_path = pan_dir / "gene_presence_absence.Rtab"
    try:
        if rtab_path.exists():
            df_rtab = pd.read_csv(rtab_path, sep="\t", index_col=0)
            mat_t   = df_rtab.astype(int).T  # filas = muestras, columnas = genes
            resumen = pd.DataFrame({
                "ID":              mat_t.index,
                "n_genes_panaroo": (mat_t != 0).sum(axis=1).values,
            })
            return resumen

        pan_csv = pan_dir / "gene_presence_absence_roary.csv"
        if not pan_csv.exists():
            pan_csv = next(pan_dir.glob("gene_presence_absence*.csv"), None)
        if not (pan_csv and pan_csv.exists()):
            return pd.DataFrame()

        df = pd.read_csv(pan_csv, index_col=0, low_memory=False)
        sample_cols = [c for c in df.columns if c not in PANAROO_METADATA_COLS]
        mat   = df[sample_cols]
        mat_t = mat.T.notna() & (mat.T != "")
        resumen = pd.DataFrame({
            "ID":              mat_t.index,
            "n_genes_panaroo": mat_t.sum(axis=1).values,
        })
        return resumen
    except Exception:
        return pd.DataFrame()


def _normalizar_id(raw_id: str, ids_tabla: list) -> str:

    import re
    raw = str(raw_id).strip()
    candidatos = [raw]
    candidatos.append(Path(raw).stem)
    candidatos.append(raw.replace(".", "_"))
    dot = re.sub(r'_(\d+)$', r'.\1', raw)
    if dot != raw:
        candidatos.append(dot)
    candidatos.append(raw.replace("-", "_"))
    candidatos.append(raw.replace("_", "-"))

    for cand in candidatos:
        for tid in ids_tabla:
            if cand == tid:
                return tid
    for cand in candidatos:
        for tid in ids_tabla:
            if cand and cand in tid:
                return tid
            if tid and tid in cand:
                return tid
    return raw


def build_unified_table(home: Path, fu: Path, add_meta: bool = False,
                        metadata_file: Path = None) -> Path:

    table_path = home / "tabla_unificada.tsv"
    ids = [f.stem for f in sorted(fu.glob("*.fasta"))]
    df  = pd.DataFrame({"ID": ids})

    # Tamaño y GC desde los FASTA
    sizes, gcs = [], []
    for fid in ids:
        fasta = fu / f"{fid}.fasta"
        seqs  = list(SeqIO.parse(fasta, "fasta")) if fasta.exists() else []
        if seqs:
            seq_str = "".join(str(r.seq) for r in seqs)
            sizes.append(len(seq_str))
            gcs.append(round((seq_str.count("G") + seq_str.count("C"))
                             / len(seq_str) * 100, 2))
        else:
            sizes.append(None)
            gcs.append(None)
    df["Tamaño"] = sizes
    df["gc"]     = gcs

    # Metadatos opcionales ingresados por el usuario
    if metadata_file and Path(metadata_file).exists():
        sep_meta = "," if str(metadata_file).endswith(".csv") else "\t"
        df_meta = pd.read_csv(metadata_file, sep=sep_meta)
        # Normalizar nombre de columna de sample
        for col in ["sample", "Sample", "SAMPLE", "cepa", "strain"]:
            if col in df_meta.columns:
                df_meta = df_meta.rename(columns={col: "sample"})
                break
        if "sample" in df_meta.columns:
            df_meta["sample_norm_dash"] = df_meta["sample"].str.replace(".", "-", regex=False)
            df_meta["sample_norm_us"]   = df_meta["sample"].str.replace(".", "_", regex=False)

            def get_meta(plasmid_id, col):
                import re as _re3
                if "__" in plasmid_id:
                    sample = plasmid_id.split("__")[0]
                    sample = _re3.sub(r'_contig_\d+$', '', sample)
                else:
                    sample = plasmid_id
                row = df_meta[df_meta["sample"] == sample]
                if row.empty:
                    row = df_meta[df_meta["sample_norm_dash"] == sample]
                if row.empty:
                    row = df_meta[df_meta["sample_norm_us"] == sample]
                if row.empty:
                    candidatos_validos = [c for c in df_meta["sample"]
                                          if sample.startswith(str(c))]
                    if candidatos_validos:
                        mejor = max(candidatos_validos, key=len)
                        row = df_meta[df_meta["sample"] == mejor]
                if not row.empty and col in row.columns:
                    return row[col].iloc[0]
                return ""
            cols_meta_agregadas = []
            for col in ["Niche", "Especie", "Pais"]:
                if col in df_meta.columns:
                    df[col] = df["ID"].apply(lambda x: get_meta(x, col))
                    cols_meta_agregadas.append(col)
            n_match = (df[cols_meta_agregadas[0]] != "").sum() if cols_meta_agregadas else 0
            print(f"   [metadata] {metadata_file.name}: columns {cols_meta_agregadas} "
                  f"added — {n_match}/{len(df)} plasmids matched.")
            if n_match == 0:
                print(f"   ⚠  NO plasmid matched the metadata. Check that the "
                      f"'sample' column of the CSV matches the input file "
                      f"name (without extension).")
        else:
            print(f"   ⚠  Metadata: file {metadata_file} has no "
                  f"'sample'/'Sample'/'cepa'/'strain' column — skipped.")
    elif metadata_file:
        print(f"   ⚠  Metadata: file not found → {metadata_file} — skipped.")
    elif add_meta:
        print("\n   Metadata per plasmid (Enter to leave blank):")
        niches, envs, paises, especies = [], [], [], []
        for fid in ids:
            print(f"\n   ── {fid}")
            niches.append( input("     Niche [Environment/Clinic/Undetermined]: ").strip())
            envs.append(   input("     Environment [clinic/environment/undetermined]: ").strip())
            paises.append( input("     Country/region of origin: ").strip())
            especies.append(input("     Species (e.g.: Klebsiella pneumoniae): ").strip())
        df["Niche"]       = niches
        df["Environment"] = envs
        df["Pais"]        = paises
        df["Especie"]     = especies

    # Fusionar resultados de cada herramienta normalizando IDs
    for loader, tag in [
        (_load_mob_typer,     "mob"),
        (_load_plasmidfinder, "pf"),
        (_load_genomad,       "genomad"),
        (_load_platon,        "platon"),
        (_load_abricate,      "abr"),
        (_load_panaroo,       "panaroo"),
    ]:
        sub = loader(home)
        if sub.empty or "ID" not in sub.columns:
            continue
        sub = sub.copy()
        sub["ID"] = sub["ID"].astype(str).str.strip().apply(
            lambda x: _normalizar_id(x, ids)
        )
        if sub["ID"].duplicated().any():
            n_dup = sub["ID"].duplicated().sum()
            print(f"   ⚠  [{tag}] {n_dup} duplicate IDs — keeping only the first.")
            sub = sub.drop_duplicates(subset="ID", keep="first")
        # columnas nuevas
        cols_nuevas = ["ID"] + [c for c in sub.columns
                                 if c != "ID" and c not in df.columns]
        sub = sub[cols_nuevas]
        # Agregar por ID normalizado
        df = df.merge(sub, on="ID", how="left")

    # Salvaguarda final: la tabla unificada nunca debe tener más filas
    # que IDs únicos (uno por plásmido)
    if df["ID"].duplicated().any():
        df = df.drop_duplicates(subset="ID", keep="first")

    df.to_csv(table_path, sep="\t", index=False)
    print(f"✔  Unified table → {table_path}")
    return table_path


# ═══════════════════════════════════════════════════════════════════
#  MÓDULO 5 – PICCIS SCORE
# ═══════════════════════════════════════════════════════════════════

def calcular_piccis_score(home: Path, fu: Path) -> pd.DataFrame:

    ids_unicos_raw = [
        f.stem for f in sorted(fu.glob("*.fasta"))
        if f.is_file() and "__" in f.stem
    ]

    # Leer BLAST resumido para saber qué genomad matcheó con otro programa
    blast_tsv = home / "resultados-blast-resumido.tsv"
    genomad_con_match: set[str] = set()
    if blast_tsv.exists():
        try:
            with open(blast_tsv) as fh:
                for line in fh:
                    if line.startswith("sample") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 3:
                        continue
                    q = Path(parts[1]).stem  # quitar .fasta
                    s = Path(parts[2]).stem
                    q_genomad = "__genomad__" in q
                    s_genomad = "__genomad__" in s
                    # Si matcheó con algo que no es genomad → confirmado
                    if q_genomad and not s_genomad:
                        genomad_con_match.add(q)
                    if s_genomad and not q_genomad:
                        genomad_con_match.add(s)
                    # Si dos genomad matchean entre sí no cuenta como confirmación
        except Exception as e:
            print(f"   ⚠  Error reading summarized BLAST: {e}")

    # Filtro si genomad no matcheó con ningún otro programa
    ids_unicos = []
    for fid in ids_unicos_raw:
        if "__genomad__" in fid and fid not in genomad_con_match:
            continue  # genomad sin confirmación → excluir
        ids_unicos.append(fid)

    # Cepas con plásmidos — usando el nombre del archivo fuente (sin _contig_N)
    import re as _re2
    cepas_con_plasmido: set[str] = set()
    for fid in ids_unicos:
        if "__" in fid:
            sample = fid.split("__")[0]
            # Eliminar _contig_N para obtener el nombre del archivo fuente
            sample = _re2.sub(r'_contig_\d+$', '', sample)
            cepas_con_plasmido.add(sample)

    # Todas las cepas de entrada (FASTQ + GBK/FASTA)
    cepas_input_conocidas: set[str] = set()

    # Cepas FASTQ: de FASTQ-unicycler/
    unicycler_dir = home / "FASTQ-unicycler"
    if unicycler_dir.exists():
        for d in unicycler_dir.iterdir():
            if d.is_dir():
                cepas_input_conocidas.add(d.name)

    # Cepas GBK/FASTA: de GBK-FASTA-plasmid/
    gbk_dir = home / "GBK-FASTA-plasmid"
    if gbk_dir.exists():
        for f in gbk_dir.glob("*.fasta"):
            fuente = _re2.sub(r'_contig_\d+$', '', f.stem)
            cepas_input_conocidas.add(fuente)


    def _variantes(nombre: str) -> set:
        return {nombre, nombre.replace("-", "_"), nombre.replace("_", "-")}

    cepas_con_plasmido_norm: set[str] = set()
    for c in cepas_con_plasmido:
        cepas_con_plasmido_norm |= _variantes(c)

    cepas_sin_plasmidos = sorted(
        c for c in cepas_input_conocidas
        if not (_variantes(c) & cepas_con_plasmido_norm)
    )

    # IDs finales: plásmidos + cepas sin plásmidos
    ids = ids_unicos + cepas_sin_plasmidos

    # Leer BLAST resumido para transferir detectores 
    blast_aliases: dict[str, set] = {fid: {fid} for fid in ids}
    blast_tsv = home / "resultados-blast-resumido.tsv"
    if blast_tsv.exists():
        try:
            df_blast = pd.read_csv(blast_tsv, sep="\t")

            def _registrar_alias(winner: str, descartado: str):
                for w in {winner, winner.replace("-", "_"), winner.replace("_", "-")}:
                    if w not in blast_aliases:
                        blast_aliases[w] = {w}
                    blast_aliases[w].add(descartado)

            for _, row in df_blast.iterrows():
                q = Path(str(row.get("query", ""))).stem
                s = Path(str(row.get("subject", ""))).stem
                aq = str(row.get("accion_query", "")).strip().lower()
                as_ = str(row.get("accion_subject", "")).strip().lower()
                if "discard" in as_:
                    # q es el winner, s es el descartado
                    _registrar_alias(q, s)
                elif "discard" in aq:
                    # s es el winner, q es el descartado
                    _registrar_alias(s, q)
        except Exception:
            pass

    def _match_any(candidatos: list, mapa: dict) -> int:
        for c in candidatos:
            if mapa.get(c, 0):
                return 1
        return 0

    def _candidatos(fid: str) -> list:

        import re
        base = [fid]
        if "__" in fid:
            sample = fid.split("__")[0]
            base.append(sample)
            parts = fid.split("__", 2)
            if len(parts) == 3:
                nombre = parts[2]
                base.append(nombre)
                base.append(Path(nombre).stem)
                # Variante con guion en lugar de guion bajo
                nombre_h = nombre.replace("_plasmid_", "-plasmid_")
                base.append(nombre_h)
                base.append(Path(nombre_h).stem)

        expandidos = []
        for c in base:
            expandidos.append(c)
            # Punto a guion bajo 
            dot = re.sub(r'_(\d+)$', r'.\1', c)
            if dot != c:
                expandidos.append(dot)
            # Guion a guion bajo
            expandidos.append(c.replace("-", "_"))
            expandidos.append(c.replace("_", "-"))

        # Agregar alias a la secuencia
        for alias in blast_aliases.get(fid, set()):
            if alias != fid:
                expandidos += _candidatos_simple(alias)

        return list(dict.fromkeys(expandidos))

    def _candidatos_simple(fid: str) -> list:
        import re
        base = [fid]
        if "__" in fid:
            base.append(fid.split("__")[0])
            parts = fid.split("__", 2)
            if len(parts) == 3:
                nombre = parts[2]
                base.append(nombre)
                base.append(Path(nombre).stem)
                base.append(nombre.replace("_plasmid_", "-plasmid_"))
                base.append(Path(nombre.replace("_plasmid_", "-plasmid_")).stem)
        expandidos = []
        for c in base:
            expandidos.append(c)
            dot = re.sub(r'_(\d+)$', r'.\1', c)
            if dot != c:
                expandidos.append(dot)
            expandidos.append(c.replace("-", "_"))
            expandidos.append(c.replace("_", "-"))
        return list(dict.fromkeys(expandidos))
###spades
    spades_map: dict[str, int] = {}
    todos_plas = home / "todos_plasmidos"
    if todos_plas.exists():
        for f in todos_plas.glob("*_spades_*.fasta"):
            spades_map[f.stem] = 1
            spades_map[f.stem.replace("_spades_", "-spades_")] = 1
    spades_proc = home / "FASTQ-plasmid-procesados"
    if spades_proc.exists():
        for f in spades_proc.glob("*.fasta"):
            if "component_" in f.stem or "plasmid" in f.stem:
                spades_map[f.stem] = 1
                spades_map[f.stem.replace("-", "_")] = 1
                spades_map[f.stem.replace("_", "-")] = 1
#MOB-recon
    mob_map: dict[str, int] = {}
    # Leer desde todos_plasmidos/ (nombre: <cepa>_mob_recon_plasmid_*)
    if todos_plas.exists():
        for f in todos_plas.glob("*_mob_recon_*.fasta"):
            mob_map[f.stem] = 1
            mob_map[f.stem.replace("_mob_recon_", "-mob_recon_")] = 1
    for mob_dir in [home / "MOB-recon-plasmid", home / "MOB-recon-gbk"]:
        if not mob_dir.exists():
            continue
        for sample_dir in mob_dir.iterdir():
            if not sample_dir.is_dir():
                continue
            for f in sample_dir.glob("plasmid_*.fasta"):
                mob_map[f.stem] = 1
                key = f"{sample_dir.name}_mob_recon_{f.stem}"
                mob_map[key] = 1
                mob_map[key.replace("-", "_")] = 1
#geNomad
    genomad_map: dict[str, int] = {}
    genomad_out = home / "genomad_out"
    if genomad_out.exists():
        for sample_dir in genomad_out.iterdir():
            if not sample_dir.is_dir():
                continue
            for summary in sample_dir.glob("*_plasmid_summary.tsv"):
                try:
                    df_g = pd.read_csv(summary, sep="\t")
                    if "seq_name" in df_g.columns:
                        for seq_id in df_g["seq_name"].tolist():
                            genomad_map[str(seq_id).strip()] = 1
                            # También con prefijo de cepa
                            genomad_map[f"{sample_dir.name}_genomad_{seq_id}"] = 1
                except Exception:
                    pass

    if todos_plas.exists():
        for f in todos_plas.glob("*_genomad_*.fasta"):
            genomad_map[f.stem] = 1

#Platon
    platon_map: dict[str, int] = {}
    for platon_pl in [home / "platon_plasmids", home / "platon_plasmids_gbk",
                      home / "platon_out"]:
        if not platon_pl.exists():
            continue
        # Archivos planos (platon_plasmids/)
        for f in platon_pl.glob("*.fasta"):
            safe = _safe_name(f.stem)
            platon_map[safe] = 1
            # Con y sin _platon al final
            platon_map[safe.rstrip("_platon")] = 1
            if not safe.endswith("_platon"):
                platon_map[safe + "_platon"] = 1
        # Subdirectorios (platon_out/)
        for sample_dir in platon_pl.iterdir():
            if not sample_dir.is_dir():
                continue
            for pf in sample_dir.glob("*.plasmid.fasta"):
                from Bio import SeqIO as _SeqIO2
                try:
                    for rec in _SeqIO2.parse(pf, "fasta"):
                        safe = _safe_name(rec.id)
                        platon_map[safe] = 1
                        platon_map[sample_dir.name] = 1
                        platon_map[_safe_name(sample_dir.name)] = 1
                except Exception:
                    pass
    # También desde todos_plasmidos/
    if todos_plas.exists():
        for f in todos_plas.glob("*_platon*.fasta"):
            platon_map[f.stem] = 1
            platon_map[_safe_name(f.stem)] = 1

#PlasmidFinder
    pf_finder_map: dict[str, int] = {}

    # Abricate/plasmidfinder sobre fastas-unicos/
    ab_pf = home / "Abricate" / "plasmidfinder.tab"
    if ab_pf.exists():
        try:
            df_abpf = pd.read_csv(ab_pf, sep="\t")
            if not df_abpf.empty and "#FILE" in df_abpf.columns:
                for _id in df_abpf["#FILE"].apply(lambda x: Path(x).stem).unique():
                    pf_finder_map[str(_id).strip()] = 1
        except Exception:
            pass

    # Abricate/plasmidfinder sobre ensamblados completos (unicycler-abricate/)
    ab_repl_dir = home / "unicycler-abricate"
    if ab_repl_dir.exists():
        for sample_dir in ab_repl_dir.iterdir():
            if not sample_dir.is_dir():
                continue
            tab = sample_dir / "plasmidfinder.tab"
            if not tab.exists():
                continue
            try:
                df_ar = pd.read_csv(tab, sep="\t")
                if not df_ar.empty and "SEQUENCE" in df_ar.columns:
                    for seq_id in df_ar["SEQUENCE"].unique():
                        key = f"{sample_dir.name}__abricate__{_safe_name(seq_id)}"
                        pf_finder_map[key] = 1
                        pf_finder_map[seq_id] = 1
            except Exception:
                pass

    # Fragmentos abricate en todos_plasmidos/
    if todos_plas.exists():
        for f in todos_plas.glob("*__abricate__*.fasta"):
            pf_finder_map[f.stem] = 1

##############################################
    #Construir tabla
#############################################

    DETECTORES_FASTQ = {
        "det_spades":        spades_map,
        "det_mob_recon":     mob_map,
        "det_genomad":       genomad_map,
        "det_platon":        platon_map,
        "det_plasmidfinder": pf_finder_map,  # Abricate/plasmidfinder
    }

    DETECTORES_GBK = {
        "det_platon":        platon_map,
        "det_mob_recon":     mob_map,
        "det_genomad":       genomad_map,
        "det_plasmidfinder": pf_finder_map,  # Abricate/plasmidfinder
    }

    def _es_gbk(fid: str) -> bool:

        if "__" in fid:
            sample = fid.split("__")[0]
        else:
            sample = fid
        # Los samples FASTQ vienen de FASTQ-unicycler/
        fastq_samples = set()
        unicycler_dir = home / "FASTQ-unicycler"
        if unicycler_dir.exists():
            fastq_samples = {d.name for d in unicycler_dir.iterdir()
                             if d.is_dir()}
        return sample not in fastq_samples

    rows = []
    for fid in ids:
        row: dict = {"ID": fid}
        candidatos = _candidatos(fid)
        es_sin_plasmido = fid in cepas_sin_plasmidos
        es_gbk = _es_gbk(fid)

        # Si es una cepa sin plásmidos → todos NaN excepto score=0
        if es_sin_plasmido:
            for col_name in DETECTORES_FASTQ:
                row[col_name] = float("nan")
            row["n_herramientas"] = 0
            row["n_detectan"]     = 0
            row["piccis_score"]   = 0.0
            row["confiabilidad"]  = "No plasmids"
            rows.append(row)
            continue

        DETECTORES = DETECTORES_GBK if es_gbk else DETECTORES_FASTQ
        todos_cols  = list(DETECTORES_FASTQ.keys())
        valores_disponibles = []

        for col_name in todos_cols:
            mapa = DETECTORES_FASTQ[col_name]
            if col_name not in DETECTORES:
                # No aplica para este origen → NaN
                row[col_name] = float("nan")
            elif mapa:
                val = _match_any(candidatos, mapa)
                row[col_name] = val
                valores_disponibles.append(val)
            else:
                # Herramienta corrió pero mapa vacío → 0
                row[col_name] = 0
                valores_disponibles.append(0)

        n_herramientas = len(valores_disponibles)
        n_detectan     = sum(v for v in valores_disponibles if v == 1)
        score          = round(n_detectan / n_herramientas, 4) if n_herramientas else 0.0
        confiabilidad  = ("High" if score >= 0.60 else
                          "Medium" if score >= 0.40 else "Low")

        row["n_herramientas"] = n_herramientas
        row["n_detectan"]     = n_detectan
        row["piccis_score"]   = score
        row["confiabilidad"]  = confiabilidad
        rows.append(row)

    df_score = pd.DataFrame(rows)
    out = home / "piccis_score.tsv"
    df_score.to_csv(out, sep="\t", index=False)
    print(f"✔  PICCIS Score → {out}  ({len(df_score)} plasmids)")
    return df_score


def grafico_piccis_score(df_score: pd.DataFrame, gdir: Path):

    if df_score.empty:
        print("⚠  PICCIS Score plot skipped: no data.")
        return

    det_cols = [c for c in df_score.columns if c.startswith("det_")]
    labels   = [c.replace("det_", "").replace("_", "\n") for c in det_cols]

    CONF_COLORS = {
        "High":           "#2ecc71",
        "Medium":         "#f39c12",
        "Low":            "#e74c3c",
        "No plasmids":    "#bdc3c7",
    }

    df_sorted = df_score.sort_values("piccis_score", ascending=True).reset_index(drop=True)
    n = len(df_sorted)

    fig, axes = plt.subplots(
        1, 2,
        figsize=(max(14, len(det_cols) * 2 + 4), max(6, n * 0.45 + 2)),
        gridspec_kw={"width_ratios": [len(det_cols), 5]},
    )

    # ── Panel A: heatmap de detección ───────────────────────────────
    ax0 = axes[0]
    mat = df_sorted[det_cols].values.astype(float)

    # colormap personalizado: NaN = gris, 0 = rojo, 1 = verde
    cmap_det = plt.cm.RdYlGn
    ax0.imshow(mat, aspect="auto", cmap=cmap_det,
                           vmin=0, vmax=1,
                           interpolation="nearest")

    # marcar NaN en gris
    nan_mask = np.isnan(mat)
    if nan_mask.any():
        ax0.imshow(
            np.where(nan_mask, 0.5, np.nan),
            aspect="auto", cmap=plt.cm.Greys,
            vmin=0, vmax=1, interpolation="nearest", alpha=0.6,
        )

    ax0.set_xticks(range(len(det_cols)))
    ax0.set_xticklabels(labels, fontsize=8, ha="center")
    ax0.set_yticks(range(n))
    ax0.set_yticklabels(df_sorted["ID"], fontsize=7)
    ax0.set_title("Detection per tool\n(green=yes, red=no, gray=N/A)",
                  fontsize=9)

    # anotaciones numéricas en cada celda
    for i in range(n):
        for j in range(len(det_cols)):
            val = mat[i, j]
            if not np.isnan(val):
                ax0.text(j, i, str(int(val)), ha="center", va="center",
                          fontsize=7, color="white" if val == 0 else "black",
                          fontweight="bold")
            else:
                ax0.text(j, i, "N/D", ha="center", va="center",
                          fontsize=6, color="grey")

    # ── Panel B: barras horizontales por score ───────────────────────
    ax1 = axes[1]
    colors_bar = [CONF_COLORS[c] for c in df_sorted["confiabilidad"]]
    scores     = df_sorted["piccis_score"].fillna(0).values

    bars = ax1.barh(range(n), scores, color=colors_bar, edgecolor="white", height=0.7)
    ax1.set_xlim(0, 1.15)
    ax1.set_yticks(range(n))
    ax1.set_yticklabels(df_sorted["ID"], fontsize=7)
    ax1.set_xlabel("PICCIS Score  (0 = none detect · 1 = all detect)")
    ax1.set_title("PICCIS plasmid reliability Score", fontsize=9)
    ax1.axvline(x=0.6, color="green",  linestyle="--", linewidth=0.8, alpha=0.6)
    ax1.axvline(x=0.4, color="orange", linestyle="--", linewidth=0.8, alpha=0.6)
    ax1.text(0.81, n - 0.5, "High",   fontsize=7, color="green",  alpha=0.8)
    ax1.text(0.51, n - 0.5, "Medium", fontsize=7, color="orange", alpha=0.8)

    # etiquetas de score en cada barra
    for i, (bar, sc, nd, nh) in enumerate(
        zip(bars, scores,
            df_sorted["n_detectan"], df_sorted["n_herramientas"])
    ):
        ax1.text(sc + 0.02, i,
                  f"{sc:.2f}  ({nd}/{nh})",
                  va="center", fontsize=7)

    # leyenda de confiabilidad
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, color=v, label=k)
        for k, v in CONF_COLORS.items()
    ]
    ax1.legend(handles=legend_handles, title="Reliability",
               loc="lower right", fontsize=7, title_fontsize=8)

    plt.suptitle("PICCIS – Plasmid detection reliability",
                 fontsize=11, fontweight="bold", y=1.01)
    plt.tight_layout()
    guardar_fig(fig, "PICCIS_score", gdir)
    


# ═══════════════════════════════════════════════════════════════════
#  MÓDULO 6 – GRÁFICOS
# ═══════════════════════════════════════════════════════════════════

def _encontrar_codo(valores) -> int:
    arr    = np.array(valores, dtype=float)
    n      = len(arr)
    puntos = np.vstack((range(n), arr)).T
    linea  = puntos[-1] - puntos[0]
    norma  = np.linalg.norm(linea)
    if norma == 0:
        return 0
    linea  /= norma
    vecs    = puntos - puntos[0]
    proyec  = np.outer(np.dot(vecs, linea), linea)
    dist    = np.linalg.norm(vecs - proyec, axis=1)
    return int(np.argmax(dist))


def _color_subtree(clade, color):
    clade.color = color
    for sub in clade.clades:
        _color_subtree(sub, color)

#g1 torta por nicho
def grafico_ambiente(df: pd.DataFrame, gdir: Path):
    col  = "Niche" if "Niche" in df.columns else df.columns[1]
    dist = df[col].fillna("Undetermined").value_counts()
    colors  = [NICHE_COLORS.get(k, "grey") for k in dist.index]
    explode = [0.05] * len(dist)

    fig, ax = plt.subplots(figsize=(8, 6))
    wedges, _, _ = ax.pie(dist.values, labels=None, autopct="%1.1f%%",
                          colors=colors, startangle=90,
                          wedgeprops={"edgecolor": "black"},
                          pctdistance=0.75, explode=explode)
    for i, (lbl, sz) in enumerate(zip(dist.index, dist.values)):
        ax.text(1.5, 0.5 - i*0.2, f"{lbl}: {sz}", fontsize=11, fontweight="bold")
    ax.legend(wedges, dist.index, title="Categories",
              loc="upper right", bbox_to_anchor=(1.35, 1))
    ax.set_title("Distribution by Environment")
    guardar_fig(fig, "G1_ambiente", gdir)


#G2 Barras apiladas plásmidos por especie y Niche. 
def grafico_descripcion(df: pd.DataFrame, gdir: Path):
    if "Especie" not in df.columns or "Niche" not in df.columns:
        print("⚠  G2 skipped: missing 'Especie' or 'Niche'.");  return
    data = df.copy()
    data["Especie2"] = data["Especie"].str.split().str[1].fillna(data["Especie"])
    filtrada = data[data["Niche"].isin(NICHE_COLORS)]
    tabla    = filtrada.groupby(["Especie2","Niche"]).size().unstack(fill_value=0)
    colors   = [NICHE_COLORS[c] for c in tabla.columns if c in NICHE_COLORS]

    fig, ax = plt.subplots(figsize=(12, 6))
    tabla.plot(kind="bar", stacked=True, ax=ax, color=colors[:len(tabla.columns)])
    ax.set_title("Plasmids by Species and Niche")
    ax.set_xlabel("Species");  ax.set_ylabel("Count")
    plt.xticks(rotation=45, ha="right")
    ax.legend(title="Niche")
    texto = ""
    for esp, row in tabla.iterrows():
        partes = "  ".join(f"{n}={row.get(n,0)}" for n in NICHE_COLORS)
        texto += f"{esp}: {partes}  Total={row.sum()}\n"
    plt.figtext(0.92, 0.5, texto, fontsize=7, va="center")
    plt.subplots_adjust(right=0.82)
    guardar_fig(fig, "G2_descripcion", gdir)


#G3 histograma de tamaños
def grafico_tamano(df: pd.DataFrame, gdir: Path):
    if "Tamaño" not in df.columns:
        print("⚠  G3 skipped: missing 'Tamaño'.");  return
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.histplot(df["Tamaño"].dropna(), bins=30, kde=True, color="purple", ax=ax)
    ax.set_title("Size distribution")
    ax.set_xlabel("Size (bp)");  ax.set_ylabel("Plasmid count")
    for p in ax.patches:
        if p.get_height() > 0:
            ax.annotate(f"{int(p.get_height())}",
                        (p.get_x()+p.get_width()/2, p.get_height()),
                        ha="center", va="bottom", fontsize=7)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(_fmt_x))
    plt.tight_layout()
    guardar_fig(fig, "G3_tamano", gdir)


#G4 Dispersión Tamaño vs %GC
def grafico_tamano_gc(df: pd.DataFrame, gdir: Path):
    if "Tamaño" not in df.columns or "gc" not in df.columns:
        print("⚠  G4 skipped.");  return
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(x="Tamaño", y="gc", data=df,
                    alpha=0.7, color="mediumpurple", ax=ax)
    ax.set_title("Size vs %GC")
    ax.set_xlabel("Size (bp)");  ax.set_ylabel("%GC")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(_fmt_x))
    plt.xticks(rotation=45);  plt.tight_layout()
    guardar_fig(fig, "G4_tamano_gc", gdir)


#G5 Boxplot %GC por Niche + test Mann-Whitney
def grafico_gc_ambiente(df: pd.DataFrame, gdir: Path):
    if "gc" not in df.columns or "Niche" not in df.columns:
        print("⚠  G5 skipped.");  return
    order   = [n for n in NICHE_COLORS if n in df["Niche"].unique()]
    palette = {n: NICHE_COLORS[n] for n in order}
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(x="Niche", y="gc", data=df,
                order=order, palette=palette, ax=ax)
    ax.set_title("%GC distribution by Environment")
    ax.set_xlabel("Environment");  ax.set_ylabel("%GC")
    plt.xticks(rotation=45)
    print("\n   Mann-Whitney U test (%GC):")
    for a, b in itertools.combinations(order, 2):
        g1 = df.loc[df["Niche"]==a, "gc"].dropna()
        g2 = df.loc[df["Niche"]==b, "gc"].dropna()
        if len(g1) and len(g2):
            stat, p = mannwhitneyu(g1, g2, alternative="two-sided")
            print(f"     {a} vs {b}: U={stat:.2f}, p={p:.4f}")
    plt.tight_layout()
    guardar_fig(fig, "G5_gc_ambiente", gdir)


#G6: Torta de predicted_mobility.
def grafico_mobility_torta(df: pd.DataFrame, gdir: Path):
    """G6 – """
    col = "predicted_mobility"
    if col not in df.columns:
        print("⚠  G6 skipped.");  return
    data    = df[col].replace(["-",""], "Undetermined").fillna("Undetermined")
    dist    = data.value_counts()
    colors  = ["mediumblue","thistle","lavender","lightgray"]
    explode = [0.05] * len(dist)
    fig, ax = plt.subplots(figsize=(8, 8))
    porciones, _, _ = ax.pie(dist.values, colors=colors[:len(dist)],
                              autopct="%1.1f%%", startangle=140,
                              explode=explode[:len(dist)],
                              wedgeprops={"edgecolor":"black"})
    ax.set_title("predicted_mobility distribution")
    ax.axis("equal")
    ax.legend(porciones, dist.index, title="Categories",
              loc="upper left", bbox_to_anchor=(1, 1))
    for i, (cat, cant) in enumerate(dist.items()):
        ax.text(1.2, 0.5-i*0.1, f"{cat}: {cant}",
                transform=ax.transAxes, va="top")
    guardar_fig(fig, "G6_mobility_torta", gdir)


# G7 / G8 (helper compartido) 
def _barras_relaxasa(df: pd.DataFrame, mobility: str,
                     titulo: str, nombre: str, gdir: Path):
    needed = {"predicted_mobility", "relaxase_type(s)"}
    if not needed.issubset(df.columns):
        print(f"⚠  {nombre} skipped: missing columns.");  return
    sub = df[df["predicted_mobility"] == mobility].copy()
    if sub.empty:
        print(f"⚠  No records with mobility='{mobility}'.");  return
    sub["genero"] = (sub["Especie"].str.split().str[1].fillna("Desconocido")
                     if "Especie" in sub.columns else "Desconocido")
    expanded = (sub.assign(
        relaxase_type=sub["relaxase_type(s)"].fillna("Undetermined").str.split(","))
        .explode("relaxase_type"))
    expanded["relaxase_type"] = expanded["relaxase_type"].str.strip()
    tabla = pd.crosstab(expanded["genero"], expanded["relaxase_type"])
    fig, ax = plt.subplots(figsize=(12, 8))
    tabla.plot(kind="bar", stacked=True, colormap="viridis", ax=ax)
    for p in ax.patches:
        h = p.get_height()
        if h > 0:
            ax.annotate(f"{int(h)}",
                        (p.get_x()+p.get_width()/2, p.get_y()+h/2),
                        ha="center", va="center", fontsize=8)
    ax.set_title(titulo);  ax.set_xlabel("Genus");  ax.set_ylabel("Frequency")
    plt.xticks(rotation=45, ha="right");  plt.tight_layout()
    guardar_fig(fig, nombre, gdir)


def grafico_tipos_movilizables(df: pd.DataFrame, gdir: Path):
    _barras_relaxasa(df, "mobilizable",
                     "Mobilizable plasmids – Genus × Relaxase",
                     "G7_movilizables", gdir)


def grafico_tipos_conjugativos(df: pd.DataFrame, gdir: Path):
    _barras_relaxasa(df, "conjugative",
                     "Conjugative plasmids – Genus × Relaxase",
                     "G8_conjugativos", gdir)


# G9 archivos para genrar planisferio y el planisferio con especie y cantidades

def obtener_shapefile(home: Path) -> Path:
    
    ne_dir = home / "ne_countries"
    shp    = ne_dir / "ne_110m_admin_0_countries.shp"

    if shp.exists():
        print(f"   [G9] Shapefile already available → {shp}")
        return shp

    mkdir(ne_dir)
    zip_path = ne_dir / "ne_110m_admin_0_countries.zip"
    url = ("https://naturalearth.s3.amazonaws.com/110m_cultural/"
           "ne_110m_admin_0_countries.zip")

    print(f"   [G9] Downloading Natural Earth shapefile (~60 KB)...")
    try:
        import urllib.request
        urllib.request.urlretrieve(url, zip_path)
    except Exception as e:
        print(f"   ⚠  Could not download the shapefile: {e}")
        print(f"       Download it manually from:")
        print(f"       {url}")
        print(f"       and unzip into: {ne_dir}/")
        return None

    import zipfile
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(ne_dir)
    zip_path.unlink()   # borrar el zip una vez descomprimido

    print(f"   [G9] Shapefile ready → {shp}")
    return shp


def grafico_planisferio(df: pd.DataFrame, shp_path: Path, gdir: Path):
    if not _GEO_OK:
        print("⚠  G9 skipped: install geopandas, geopy and distinctipy.");  return
    if not shp_path.exists():
        print(f"⚠  G9 skipped: shapefile not found.");  return
    col_pais = next((c for c in ["Pais","pais","Country","country","Columna12"]
                     if c in df.columns), None)
    if col_pais is None:
        print("⚠  G9 skipped: no country column.");  return
    df = df.copy()
    df["species"] = (df["Especie"].fillna("Unknown")
                     if "Especie" in df.columns else "Unknown")

    geolocator = Nominatim(user_agent="piccis_geo")
    coords_pais = {}
    for country in df[col_pais].dropna().unique():
        try:
            loc = geolocator.geocode(str(country), timeout=15)
            if loc:
                coords_pais[country] = (loc.latitude, loc.longitude)
            else:
                print(f"   ⚠  Country not found: {country}")
        except GeocoderTimedOut:
            print(f"   ⚠  Timeout: {country}")
    if not coords_pais:
        print("⚠  G9 skipped: no geolocated countries.");  return

    # Agrupar por (país, especie) — un punto por combinación
    grp = (df[df[col_pais].isin(coords_pais)]
           .groupby([col_pais, "species"]).size()
           .reset_index(name="count"))

    especies = sorted(grp["species"].unique())
    colores  = distinctipy.get_colors(len(especies))
    color_map = dict(zip(especies, colores))

    world = gpd.read_file(shp_path)
    fig, ax = plt.subplots(figsize=(20, 12))
    world.boundary.plot(ax=ax, edgecolor="black", linewidth=0.8)
    world.plot(color="lightyellow", ax=ax, edgecolor="black")
    ax.set_facecolor("lightblue")

    # Jitter: si un país tiene varias especies, separar los puntos
    jitter_step = 1.5
    for country, sub in grp.groupby(col_pais):
        lat, lon = coords_pais[country]
        n_sp = len(sub)
        for i, (_, row) in enumerate(sub.iterrows()):
            offset = (i - (n_sp - 1) / 2) * jitter_step
            ax.scatter(lon + offset, lat,
                       s=100 + row["count"] * 60,
                       color=color_map[row["species"]],
                       alpha=0.8, zorder=5,
                       edgecolors="black", linewidths=0.5)
        # Una sola anotación de país por grupo (encima de todos los puntos)
        ax.annotate(f"{country}", (lon, lat),
                    textcoords="offset points", xytext=(0, 12),
                    fontsize=8, ha="center")

    handles = [plt.Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=color_map[sp], markersize=10,
                          label=sp)
               for sp in especies]
    ax.legend(handles=handles, title="Species",
              bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    ax.set_title("Geographic distribution of plasmids by species", fontsize=16)
    ax.set_xlabel("Longitud");  ax.set_ylabel("Latitud")
    plt.tight_layout()
    guardar_fig(fig, "G9_planisferio", gdir)


# COG (EggNOG-mapper + Bakta) 
def grafico_cog(home: Path, df_meta: pd.DataFrame, gdir: Path):
    egg_dir = home / "EggNOGmapper"
    frames  = []
    for ann in egg_dir.rglob("*.annotations"):
        try:
            lines = ann.read_text(errors="ignore").splitlines()
            header_idx = next(
                (i for i, l in enumerate(lines) if l.startswith("#query")), None
            )
            if header_idx is None:
                continue
            header = lines[header_idx].lstrip("#").split("\t")
            data_lines = [l for l in lines[header_idx + 1:]
                          if l.strip() and not l.startswith("#")]
            if not data_lines:
                continue
            from io import StringIO
            tmp = pd.read_csv(StringIO("\n".join(data_lines)),
                              sep="\t", names=header)
            tmp["ID"] = ann.stem.split(".emapper")[0]
            frames.append(tmp)
        except Exception:
            pass
    if not frames:
        print("⚠  COG skipped: no EggNOG files.");  return

    df_egg = pd.concat(frames, ignore_index=True)
    # Buscar columna COG de forma case-insensitive
    cog_col = next((c for c in df_egg.columns if c.lower() == "cog_category"), None)
    if cog_col is None:
        print("⚠  COG skipped: no 'COG_category' column.");  return
    df_egg["COG_cat_simple"] = df_egg[cog_col].astype(str).str[0].replace("n", "S").fillna("S")

    if not df_meta.empty and "ID" in df_meta.columns:
        cols = ["ID"] + [c for c in ["Especie","Niche"] if c in df_meta.columns]
        df_egg = df_egg.merge(df_meta[cols], on="ID", how="left")

    cats  = sorted(df_egg["COG_cat_simple"].unique())
    cmap  = plt.get_cmap("tab20")
    cdict = {c: cmap(i % 20) for i, c in enumerate(cats)}

    #Barras generales de distribución COG 
    counts_total = df_egg["COG_cat_simple"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.bar(counts_total.index,
                  counts_total.values,
                  color=[cdict[c] for c in counts_total.index])
    ax.set_xlabel("COG category")
    ax.set_ylabel("Number of genes")
    ax.set_title("COG functional categories — all plasmids")
    for bar, val in zip(bars, counts_total.values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                str(val), ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    guardar_fig(fig, "COG_categories", gdir)

    # Heatmap COG por plásmido 
    pivot = df_egg.groupby(["ID","COG_cat_simple"]).size().unstack(fill_value=0)
    if not pivot.empty:
        fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns)*0.6),
                                        max(4, len(pivot)*0.4)))
        sns.heatmap(pivot, cmap="YlOrRd", linewidths=0.3,
                    linecolor="lightgray", ax=ax,
                    cbar_kws={"label":"Gene count"})
        ax.set_title("COG categories per plasmid")
        ax.set_xlabel("COG category")
        ax.set_ylabel("Plasmid")
        plt.tight_layout()
        guardar_fig(fig, "COG_heatmap", gdir)

    #Tortas por Especie y por Niche (si hay metadatos)
    def _pie(ax, sub, title):
        counts = sub["COG_cat_simple"].value_counts(normalize=True)*100
        colors = [cdict[l] for l in counts.index]
        wedges, _ = ax.pie(counts, colors=colors, labels=None,
                            startangle=90, wedgeprops={"width":0.5})
        ax.legend(wedges,
                  [f"{l} – {p:.1f}%" for l,p in zip(counts.index, counts)],
                  title="COG", loc="center left",
                  bbox_to_anchor=(1,0,0.5,1), fontsize=7)
        ax.set_title(title, fontsize=10)

    for col_g, prefix in [("Especie","COG_especie"), ("Niche","COG_niche")]:
        if col_g not in df_egg.columns:
            continue
        for val in df_egg[col_g].dropna().unique():
            fig, ax = plt.subplots(figsize=(7,7))
            _pie(ax, df_egg[df_egg[col_g]==val], f"COG – {val}")
            guardar_fig(fig, f"{prefix}_{val.replace(' ','_')}", gdir)
            


#SVD + KMeans 
def grafico_svd_kmeans(home: Path, gdir: Path):

    pan_dir   = home / "panaroo_out"
    rtab_path = pan_dir / "gene_presence_absence.Rtab"

    if rtab_path.exists():
        # Matriz binaria desde el .Rtab 
        df_pan = pd.read_csv(rtab_path, sep="\t", index_col=0)
        df_bin = df_pan.astype(int).T  # filas = plásmidos, columnas = genes
        df_bin.index.name = "ID"
        print(f"   [SVD] Using {rtab_path.name} ({df_bin.shape[1]} genes)")
    else:
        csv_path = pan_dir / "gene_presence_absence_roary.csv"
        if not csv_path.exists():
            csv_path = next(pan_dir.glob("gene_presence_absence*.csv"), None)
        if csv_path is None or not csv_path.exists():
            print("⚠  SVD/KMeans skipped: no gene_presence_absence from Panaroo.")
            return None

        df_pan = pd.read_csv(csv_path, index_col=0, low_memory=False)
        sample_cols = [c for c in df_pan.columns if c not in PANAROO_METADATA_COLS]
        n_excluidas = len(df_pan.columns) - len(sample_cols)
        if n_excluidas:
            print(f"   [SVD] {n_excluidas} Roary metadata columns "
                  f"excluded from {csv_path.name}.")
        df_bin = df_pan[sample_cols].notna().astype(int).T  # filas = plásmidos
        df_bin.index.name = "ID"

    # Variables de tabla_unificada
    tabla_path = home / "tabla_unificada.tsv"
    df_extra = pd.DataFrame(index=df_bin.index)

    if tabla_path.exists():
        df_meta = pd.read_csv(tabla_path, sep="\t")
        # Eliminar duplicados de ID antes de set_index
        df_meta = df_meta.drop_duplicates(subset="ID")
        df_meta = df_meta.set_index("ID")

        # Tamaño (continua, escalada)
        if "Tamaño" in df_meta.columns:
            df_extra["Tamaño"] = df_meta["Tamaño"].reindex(df_bin.index).fillna(0)

        cat_cols = [c for c in
                    ["rep_type(s)", "relaxase_type(s)", "predicted_mobility"]
                    if c in df_meta.columns]

        if cat_cols:
            df_cat = df_meta[cat_cols].reindex(df_bin.index).fillna("Unknown")
            df_cat = df_cat.replace(["-", "", "nan"], "Unknown")
            df_ohe = pd.get_dummies(df_cat, prefix=cat_cols)
            df_extra = pd.concat([df_extra, df_ohe], axis=1)

        # Resistencias y virulencia agregadas de Abricate (numéricas)
        for col in ["n_resistance_genes", "n_virulence_genes"]:
            if col in df_meta.columns:
                df_extra[col] = pd.to_numeric(
                    df_meta[col].reindex(df_bin.index), errors="coerce"
                ).fillna(0)

    # Combinar Panaroo + variables extra
    if not df_extra.empty:
        df_extra = df_extra.fillna(0)
        # Escalar variables extra para que no dominen sobre Panaroo
        scaler_extra = StandardScaler()
        extra_sc = pd.DataFrame(
            scaler_extra.fit_transform(df_extra),
            index=df_extra.index,
            columns=df_extra.columns
        )
        df_combined = pd.concat([df_bin, extra_sc], axis=1)
        print(f"   [SVD] Panaroo: {df_bin.shape[1]} genes + "
              f"tabla_unificada: {df_extra.shape[1]} variables")
    else:
        df_combined = df_bin
        print(f"   [SVD] Panaroo only: {df_bin.shape[1]} genes")

    datos_esc = StandardScaler().fit_transform(df_combined.fillna(0))
    n_comp    = min(min(datos_esc.shape)-1, 100)
    svd       = TruncatedSVD(n_components=n_comp, random_state=42)
    reducida  = svd.fit_transform(datos_esc)

    var_acum = np.cumsum(svd.explained_variance_ratio_)
    fig, ax  = plt.subplots(figsize=(10, 5))
    ax.plot(range(1, len(var_acum)+1), var_acum, marker="o")
    ax.set_title("Varianza explicada acumulada – SVD (Panaroo + tabla_unificada)")
    ax.set_xlabel("Componentes");  ax.set_ylabel("Varianza acumulada")
    ax.grid(True);  plt.tight_layout()
    guardar_fig(fig, "SVD_varianza", gdir)

    codo_var = max(_encontrar_codo(var_acum), 2)
    reducida  = reducida[:, :codo_var]
    print(f"   Components selected by elbow: {codo_var}")

    k_range = range(2, min(11, len(df_combined)))
    inertias, silhouettes = [], []
    for k in k_range:
        km = KMeans(n_clusters=k, random_state=42, n_init="auto")
        cl = km.fit_predict(reducida)
        inertias.append(km.inertia_)
        silhouettes.append(silhouette_score(reducida, cl))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].plot(k_range, inertias, marker="o")
    axes[0].set_title("Codo");  axes[0].set_xlabel("k");  axes[0].grid(True)
    axes[1].plot(k_range, silhouettes, marker="o", color="green")
    axes[1].set_title("Silhouette");  axes[1].set_xlabel("k");  axes[1].grid(True)
    plt.tight_layout()
    guardar_fig(fig, "SVD_codo_silhouette", gdir)

    k_sil   = list(k_range)[int(np.argmax(silhouettes))]
    k_elbow = _encontrar_codo(inertias) + 2
    print(f"\n   1 – Silhouette (k={k_sil})   2 – Elbow (k={k_elbow})   3 – Manual")
    op      = ask("   Option [1/2/3]: ", ["1","2","3"])
    mejor_k = {"1":k_sil, "2":k_elbow,
                "3":int(input("   k value: ").strip())}[op]

    km_fin  = KMeans(n_clusters=mejor_k, random_state=42, n_init="auto")
    cl_fin  = km_fin.fit_predict(reducida)
    cnames  = [f"Comp_{i+1}" for i in range(reducida.shape[1])]
    df_comp = pd.DataFrame(reducida, index=df_combined.index, columns=cnames)
    df_comp["Cluster"] = cl_fin + 1

    if tabla_path.exists():
        todos_ids = pd.read_csv(tabla_path, sep="\t")["ID"].dropna().unique()
        faltantes = [i for i in todos_ids if i not in df_comp.index]
        filas_heredadas = []
        for fid in faltantes:
            cand = fid.replace("-", "_")
            if cand not in df_comp.index:
                cand = fid.replace("_", "-")
            if cand in df_comp.index:
                fila = df_comp.loc[[cand]].copy()
                fila.index = [fid]
                filas_heredadas.append(fila)
        if filas_heredadas:
            df_comp = pd.concat([df_comp] + filas_heredadas)
            print(f"   [SVD] {len(filas_heredadas)} duplicate ID(s) from "
                  f"tabla_unificada inherited the cluster from their Panaroo variant.")

    df_comp.to_csv(home / "componentes_y_clusters.csv")
    pd.DataFrame(svd.components_[:codo_var], columns=df_combined.columns,
                 index=cnames).to_csv(home / "cargas_componentes.csv")
    print(f"✔  componentes_y_clusters.csv and cargas_componentes.csv → {home}")

    fig, ax = plt.subplots(figsize=(10, 8))
    sc = ax.scatter(df_comp["Comp_1"], df_comp["Comp_2"],
                    c=df_comp["Cluster"], cmap="tab10", s=60)
    for idx in df_comp.index:
        ax.annotate(idx.split("__")[0],
                    (df_comp.loc[idx,"Comp_1"], df_comp.loc[idx,"Comp_2"]),
                    fontsize=7, alpha=0.7)
    plt.colorbar(sc, ax=ax, label="Cluster")
    ax.set_title(f"KMeans k={mejor_k} – Panaroo + tabla_unificada (SVD)")
    ax.set_xlabel("Componente 1");  ax.set_ylabel("Componente 2")
    ax.grid(True);  plt.tight_layout()
    guardar_fig(fig, "SVD_clusters", gdir)

    return df_comp


#Árbol tANI + Heatmap Panaroo 
def grafico_arbol_heatmap(home: Path, df_meta: pd.DataFrame,
                           df_clusters, gdir: Path):

    fu             = home / "fastas-unicos"
    tree_candidates = (list(fu.glob("BestTree*.tre")) + 
                       list(fu.glob("*.tre")) +
                       list((fu / "outputs" / "tANI").glob("*.tre")) +
                       list(fu.rglob("BestTree*.tre")))
    if not tree_candidates:
        print("⚠  Tree+Heatmap skipped: no .tre in fastas-unicos/");  return

    tree     = Phylo.read(str(tree_candidates[0]), "newick")
    tree_ids = [t.name for t in tree.get_terminals()]

    pan_dir   = home / "panaroo_out"
    rtab_path = pan_dir / "gene_presence_absence.Rtab"

    if rtab_path.exists():
        df_pan = pd.read_csv(rtab_path, sep="\t", index_col=0)
        df_pa  = df_pan.astype(int).T  # filas = plásmidos, columnas = genes
    else:
        csv_path = pan_dir / "gene_presence_absence_roary.csv"
        if not csv_path.exists():
            csv_path = next(pan_dir.glob("gene_presence_absence*.csv"), None)
        if csv_path is None:
            print("⚠  Tree+Heatmap skipped: no Panaroo table.");  return

        df_pan = pd.read_csv(csv_path, index_col=0, low_memory=False)
        sample_cols = [c for c in df_pan.columns if c not in PANAROO_METADATA_COLS]
        df_pa = df_pan[sample_cols].notna().astype(int).T

    df_pa.index = df_pa.index.str.replace(r"\.", "_", regex=True)
    filtered = [s for s in tree_ids if s in df_pa.index]
    df_pa    = df_pa.loc[filtered]

    # Mapas de Niche y Cluster
    niche_map, cluster_map = {}, {}
    if not df_meta.empty and "ID" in df_meta.columns:
        df_m = df_meta.copy()
        df_m["clean"] = df_m["ID"].str.replace(r"\.", "_", regex=True)
        if "Niche" in df_m.columns:
            niche_map = dict(zip(df_m["clean"], df_m["Niche"]))
    if df_clusters is not None and not df_clusters.empty:
        idx = df_clusters.index.astype(str).str.replace(r"\.", "_", regex=True)
        cluster_map = dict(zip(idx, df_clusters["Cluster"].astype(str)))

    # Colorear primeras bifurcaciones del árbol
    BRANCH_COLORS = ["#006400","#228B22","#32CD32","#3CB371","#6B8E23",
                     "#800080","#FF69B4","#FFD700","#FF6347","#4169E1"]
    internals = [c for c in tree.find_clades()
                 if not c.is_terminal() and len(c.clades) > 1]
    for i, clade in enumerate(internals[:len(BRANCH_COLORS)]):
        _color_subtree(clade, BRANCH_COLORS[i])

    unique_cl = sorted(set(cluster_map.values()) - {""})
    cl_pal    = dict(zip(unique_cl, sns.color_palette("hsv", len(unique_cl))))

    n_rows = len(df_pa)
    fig    = plt.figure(figsize=(30, n_rows*0.3 + 3))
    gs     = gridspec.GridSpec(2, 2,
                               height_ratios=[n_rows*0.3, 1.5],
                               width_ratios=[6, 4],
                               hspace=0.2, wspace=0.05)

    # Panel árbol
    ax0 = fig.add_subplot(gs[0, 0])
    Phylo.draw(tree, axes=ax0, branch_labels=False, do_show=False)
    ax0.set_title("Phylogenetic tree (tANI)", fontsize=12)
    for sp in ax0.spines.values():
        sp.set_visible(False)

    xlim  = ax0.get_xlim()
    leaf_y = {t.name: y
               for y, t in enumerate(reversed(tree.get_terminals()))
               if t.name in niche_map or t.name in cluster_map}
    x_niche   = xlim[1] + 1
    x_cluster = x_niche + 1.2

    for strain, y in leaf_y.items():
        niche   = niche_map.get(strain, "Undetermined")
        cluster = str(cluster_map.get(strain, ""))
        color_n = NICHE_COLORS.get(niche, "grey")
        if niche == "Clinic":
            ax0.add_patch(Polygon([[x_niche,y],[x_niche+0.5,y+0.25],
                                    [x_niche+0.5,y-0.25]],
                                   closed=True, color=color_n))
        elif niche == "Environment":
            ax0.add_patch(Rectangle((x_niche, y-0.25), 0.5, 0.5, color=color_n))
        else:
            ax0.add_patch(Circle((x_niche+0.25, y), 0.25, color=color_n))
        if cluster and cluster in cl_pal:
            ax0.add_patch(Circle((x_cluster+0.25, y), 0.25, color=cl_pal[cluster]))
    ax0.set_xlim(xlim[0], x_cluster + 1.5)

    gene_freq = df_pa.sum(axis=0)  # cuántas muestras tienen cada gen
    gene_order = gene_freq.sort_values(ascending=False).index
    df_pa_sorted = df_pa[gene_order]

    ax1 = fig.add_subplot(gs[0, 1])
    sns.heatmap(df_pa_sorted, annot=False,
                cmap=ListedColormap(["seashell","indigo"]),
                cbar=False, yticklabels=False, xticklabels=False,
                ax=ax1, linecolor="lightgray")
    ax1.set_xlabel("Genes (Panaroo) — sorted by frequency ↓")
    ax1.set_ylabel("Plasmids")
    ax1.set_title("Gene presence / absence", fontsize=12)

    prev_freq = None
    for i, gene in enumerate(gene_order):
        freq = gene_freq[gene]
        if prev_freq is not None and freq != prev_freq:
            ax1.axvline(x=i, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)
        prev_freq = freq

    # Leyenda
    ax_leg = fig.add_subplot(gs[1, :])
    ax_leg.axis("off")
    niche_handles = [
        Line2D([0],[0], marker="s", color="w", label="Environment",
               markerfacecolor=NICHE_COLORS["Environment"], markersize=9),
        Line2D([0],[0], marker="^", color="w", label="Clinic",
               markerfacecolor=NICHE_COLORS["Clinic"], markersize=9),
        Line2D([0],[0], marker="o", color="w", label="Undetermined",
               markerfacecolor=NICHE_COLORS["Undetermined"], markersize=9),
    ]
    cl_handles = [
        Line2D([0],[0], marker="o", color="w", label=f"Cluster {cl}",
               markerfacecolor=cl_pal[cl], markersize=8,
               markeredgecolor="black", lw=0.2)
        for cl in unique_cl
    ]
    ax_leg.legend(handles=niche_handles + cl_handles,
                  title="Niche & Cluster", loc="center",
                  ncol=6, fontsize=8, title_fontsize=9, frameon=False)
    plt.tight_layout()
    guardar_fig(fig, "Arbol_Heatmap_Panaroo", gdir)
    


# Tablas de resumen 
COLS_MULTI = {
    "relaxase_type(s)":              "resumen_relaxase_type",
    "mpf_type":                      "resumen_mpf_type",
    "predicted_mobility":            "resumen_predicted_mobility",
    "rep_type(s)":                   "resumen_rep_type",
    "primary_cluster_id":            "resumen_primary_cluster_id",
    "observed_host_range_ncbi_rank": "resumen_observed_host_range",
}

def generar_tablas_resumen(df: pd.DataFrame, out_dir: Path):
    """Generates statistical summary CSVs from the unified table."""
    df = df.copy()
    if "Especie" in df.columns:
        df["Especie2"] = df["Especie"].str.split().str[1].fillna(df["Especie"])

    if "Especie2" in df.columns:
        t = df["Especie2"].value_counts().reset_index()
        t.columns = ["Especie","Cantidad"]
        t.to_csv(out_dir / "resumen_especies.csv", index=False)

    if "Niche" in df.columns:
        t = df["Niche"].value_counts().reset_index()
        t.columns = ["Niche","Cantidad"]
        t.to_csv(out_dir / "resumen_niche.csv", index=False)

    if "Especie2" in df.columns and "Niche" in df.columns:
        t = df.groupby(["Especie2","Niche"]).size().reset_index(name="Cantidad")
        t.to_csv(out_dir / "resumen_especie_niche.csv", index=False)

    grp_base = [c for c in ["Especie2","Niche"] if c in df.columns]
    for col, nombre in COLS_MULTI.items():
        if col not in df.columns:
            continue
        tmp = df.dropna(subset=[col]).copy()
        tmp[col] = tmp[col].astype(str).str.replace(" ","").str.split(",")
        exp  = tmp.explode(col)
        gcols = [c for c in grp_base if c in exp.columns] + [col]
        t    = exp.groupby(gcols).size().reset_index(name="Cantidad")
        t.to_csv(out_dir / f"{nombre}.csv", index=False)

    print(f"✔  Summary tables → {out_dir}")


# ═══════════════════════════════════════════════════════════════════
#  PROGRAMA PRINCIPAL
# ═══════════════════════════════════════════════════════════════════

def main():
    #Argumentos de línea de comandos
    parser = argparse.ArgumentParser(
        prog="piccis_pipeline.py",
        description="PICCIS v2.0 – Plasmid Identification, Clustering "
                    "and Comparative Integrated Score",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python piccis_pipeline.py
  python piccis_pipeline.py --cores 8
  python piccis_pipeline.py --cores 1

Core selection priority:
  1. --cores N          (this flag)
  2. PICCIS_WORKERS=N   (environment variable)
  3. cpu_count() - 1    (automatic default)
        """,
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=None,
        metavar="N",
        help="Number of parallel workers to use. "
             "Default: cpu_count() - 1. "
             "Use 1 to disable parallelism (easier to debug).",
    )
    parser.add_argument("--workdir", type=str, default=None, metavar="DIR",
        help="Working directory.")
    parser.add_argument("--tipo", type=str, default=None, choices=["1","2","3"],
        help="Data type: 1=FASTQ, 2=GBK/FASTA, 3=both.")
    parser.add_argument("--fastq-dir", type=str, default=None, metavar="DIR",
        help="Folder with the FASTQ files.")
    parser.add_argument("--gbk-dir", type=str, default=None, metavar="DIR",
        help="Folder with GBK/FASTA files.")
    parser.add_argument("--library", type=str, default=None,
        choices=["1","2","3","4","5"],
        help="SPAdes library type: 1=PE, 2=SE, 3=IonTorrent, 4=PE+PacBio, 5=PE+Nanopore.")
    parser.add_argument("--metagenomic", action="store_true", default=False,
        help="Metagenomic data (uses --metaplasmid in SPAdes).")
    parser.add_argument("--metadata", type=str, default=None, metavar="TSV",
        help="TSV file with metadata: columns sample, Niche, Especie, Pais.")
    args = parser.parse_args()

    print(banner())

    #Directorio de trabajo
    if args.workdir:
        home = Path(args.workdir).resolve()
    else:
        home = Path(
            input("📁 Working directory (e.g.: home-plasmid): ").strip()
            or "home-plasmid"
        ).resolve()
    mkdir(home)
    gdir = mkdir(home / "graficos")
    print(f"   Directory: {home}")
    print(f"   Plots     : {gdir}\n")

    #Cores
    n_workers = detectar_workers(cli_cores=args.cores)
    print(f"   Parallel workers: {n_workers}  "
          f"(change with: --cores N  or  export PICCIS_WORKERS=N)\n")

    # Leer configuración de bases de datos
    conf = leer_conf()
    if conf:
        print(f"   [conf] piccis.conf found — DB paths loaded automatically")
    else:
        print(f"   [conf] piccis.conf not found — paths will be requested manually")
        print(f"          (run install_databases.sh to configure automatically)\n")

    # Recolectar todos los inputs de una vez 
    print("\n" + "─"*60)
    print("  INPUT CONFIGURATION")
    print("─"*60)

    # Tipo de entrada
    if args.tipo:
        tipo = args.tipo
    else:
        print("\nWhat input data do you have?")
        print("  1 – FASTQ only")
        print("  2 – Assembled GBK / FASTA only")
        print("  3 – FASTQ + GBK / FASTA")
        tipo = ask("Option [1/2/3]: ", ["1","2","3"])

    # Carpeta FASTQ
    fastq_dir    = None
    library_type = None
    meta_mode    = False
    if tipo in ("1","3"):
        if args.fastq_dir:
            fastq_dir = Path(args.fastq_dir).resolve()
        else:
            fastq_dir = Path(input("\n📁 FASTQ folder: ").strip()).resolve()

        lib_map = {
            "1": "paired", "2": "single", "3": "iontorrent",
            "4": "paired+pacbio", "5": "paired+nanopore",
        }
        if args.library:
            library_type = lib_map[args.library]
        else:
            print("\nSPAdes library type:")
            print("  1 – Paired-end Illumina")
            print("  2 – Single-end / SE")
            print("  3 – IonTorrent")
            print("  4 – Paired-end + PacBio CLR")
            print("  5 – Paired-end + Nanopore")
            library_type = lib_map[ask("Option [1-5]: ", list(lib_map.keys()))]

        if getattr(args, 'metagenomic', False):
            meta_mode = True
        elif getattr(args, 'no_meta', False):
            meta_mode = False
        else:
            print("\nIs the data metagenomic?")
            print("  y → --metaplasmid   n → --plasmid")
            meta_mode = ask("Option [y/n]: ", ["y","n"]) == "y"

    # Carpeta GBK/FASTA
    gbk_input = None
    if tipo in ("2","3"):
        if args.gbk_dir:
            gbk_input = Path(args.gbk_dir).resolve()
        else:
            gbk_input = Path(input(
                "\n📁 GBK/FASTA folder (or single file): "
            ).strip()).resolve()

    # Archivo de metadatos
    meta_file = Path(args.metadata).resolve() if getattr(args, 'metadata', None) else None
    if meta_file and not meta_file.exists():
        print(f"\n✖  ERROR: the metadata file passed with --metadata does not exist:")
        print(f"   {meta_file}")
        print(f"   Check the path and run the command again.")
        sys.exit(1)
    if not meta_file:
        meta_input = input(
            "\n📋 Metadata CSV/TSV file (Enter to skip): "
        ).strip()
        if meta_input:
            meta_file = Path(meta_input).expanduser().resolve()
            if not meta_file.exists():
                print(f"   ⚠  Not found: {meta_file} — metadata skipped.")
                meta_file = None

    tipo_labels = {"1":"FASTQ","2":"GBK/FASTA","3":"FASTQ + GBK/FASTA"}
    print("\n" + "─"*60)
    print("  CONFIGURATION SUMMARY")
    print("─"*60)
    print(f"  Input type      : {tipo_labels[tipo]}")
    if fastq_dir:
        print(f"  FASTQ           : {fastq_dir}")
        print(f"  Library         : {library_type}")
        print(f"  Metagenomic     : {'yes' if meta_mode else 'no'}")
    if gbk_input:
        print(f"  GBK/FASTA       : {gbk_input}")
    if meta_file:
        print(f"  Metadata        : {meta_file}")
    else:
        print(f"  Metadata        : (none)")
    print("─"*60 + "\n")

    # Carpeta central donde se acumulan todos los plásmidos reconstruidos
    todos_plasmidos = mkdir(home / "todos_plasmidos")

    # Registrar todas las cepas de entrada para el PICCIS Score
    cepas_input: list[str] = []

    
    if tipo in ("1","3"):
        check_env("spades_env",     "spades.py")
        check_env("unicycler_env",  "unicycler")
        check_env("mob_env",        "mob_recon")
        check_env("platon_env",     "platon")
        platon_db_fastq = pedir_db(conf, "PLATON_DB",
                                    "📁 Path to the Platon database: ")

        spades_raw    = mkdir(home / "FASTQ-plasmid")
        unicycler_raw = mkdir(home / "FASTQ-unicycler")
        mob_out_fq    = mkdir(home / "MOB-recon-plasmid")

        fastq_dir = descomprimir_fastqs(fastq_dir, home)

        # Registrar cepas desde los FASTQs
        cepas_input += samples_desde_fastqs(fastq_dir)

        
        print("\n   [assembly] Running SPAdes (plasmid mode)...")
        run_spades(fastq_dir, spades_raw, library_type,
                   n_workers, multiprocessing.cpu_count(), meta_mode)
        
        spades_proc = mkdir(home / "FASTQ-plasmid-procesados")
        process_spades_output(spades_raw, spades_proc)
        for f in spades_proc.glob("*.fasta"):
            
            stem = f.stem
            for tag in ("-plasmid_component_", "_plasmid_component_"):
                if tag in stem:
                    parts = stem.split(tag)
                    stem = f"{_safe_name(parts[0])}__spades__component_{parts[1]}"
                    break
            else:
                stem = _safe_name(stem)
            dest = todos_plasmidos / f"{stem}.fasta"
            if not dest.exists():
                shutil.copy(f, dest)

        print("\n   [assembly] Running Unicycler (conservative)...")
        run_unicycler(fastq_dir, unicycler_raw, library_type,
                      n_workers, multiprocessing.cpu_count())

        print("\n   [mob_recon] Running on Unicycler assembly...")
        run_mob_recon_from_dir(unicycler_raw, mob_out_fq, "assembly.fasta", n_workers)
        for sample_dir in mob_out_fq.iterdir():
            if not sample_dir.is_dir():
                continue
            for f in sample_dir.glob("plasmid_*.fasta"):
                stem = f"{_safe_name(sample_dir.name)}__mob_recon__{f.stem}"
                dest = todos_plasmidos / f"{stem}.fasta"
                if not dest.exists():
                    shutil.copy(f, dest)

        replicon_fq = extraer_contigs_con_replicon([unicycler_raw], home, n_workers)
        for f in replicon_fq.glob("*.fasta"):
            dest = todos_plasmidos / _safe_name(f.name)
            if not dest.exists():
                shutil.copy(f, dest)

        print("\n   [Platon] Running on Unicycler assemblies...")
        platon_fq_out = run_platon(unicycler_raw, home, platon_db_fastq,
                                    threads=multiprocessing.cpu_count(),
                                    n_workers=n_workers)
        for sample_dir in platon_fq_out.iterdir():
            if not sample_dir.is_dir():
                continue
            for pf in sample_dir.glob("*.plasmid.fasta"):
                stem = f"{_safe_name(sample_dir.name)}__platon__{_safe_name(sample_dir.name)}-platon-plasmid"
                dest = todos_plasmidos / f"{stem}.fasta"
                if not dest.exists():
                    shutil.copy(pf, dest)


        genomad_db_path = pedir_db(conf, "GENOMAD_DB",
                                    "📁 Path to the geNomad database: ")
        genomad_fq_out = run_genomad_dir(unicycler_raw, home, genomad_db_path, n_workers)
        collect_genomad_plasmids(genomad_fq_out, todos_plasmidos)

#RAMA GBK/FASTA
    if tipo in ("2","3"):
        check_env("platon_env", "platon")
        gbk_out = mkdir(home / "GBK-FASTA-plasmid")

        if gbk_input.is_dir():
            gbk_exts   = (".gb", ".gbk", ".genbank")
            fasta_exts = (".fasta", ".fa", ".fna", ".fas")
            archivos   = [f for f in gbk_input.iterdir()
                          if f.suffix.lower() in gbk_exts + fasta_exts]
            if not archivos:
                print(f"   ⚠  No GBK/FASTA found in {gbk_input}")
            for arch in sorted(archivos):
                print(f"   → {arch.name}")
                cepas_input.append(arch.stem)
                if arch.suffix.lower() in gbk_exts:
                    gbk_to_fastas(arch, gbk_out)
                else:
                    split_multifasta(arch, gbk_out)
        elif gbk_input.is_file():
            cepas_input.append(gbk_input.stem)
            if gbk_input.suffix.lower() in (".gb", ".gbk", ".genbank"):
                gbk_to_fastas(gbk_input, gbk_out)
            else:
                split_multifasta(gbk_input, gbk_out)

        platon_db_gbk = pedir_db(conf, "PLATON_DB",
                                   "📁 Path to the Platon database: ")

        mob_out_gbk = mkdir(home / "MOB-recon-gbk")
        trabajos_mob_gbk = []
        for fasta in sorted(gbk_out.glob("*.fasta")):
            mob_done = mob_out_gbk / fasta.stem / "chromosome.fasta"
            if mob_done.exists() and mob_done.stat().st_size > 0:
                print(f"   [mob_recon] {fasta.stem} already processed, skipping.")
            else:
                trabajos_mob_gbk.append((fasta, mkdir(mob_out_gbk / fasta.stem)))
        if trabajos_mob_gbk:
            ejecutar_en_paralelo(_mob_recon_single, trabajos_mob_gbk,
                                 n_workers, "MOB-recon GBK")
        for sample_dir in mob_out_gbk.iterdir():
            if not sample_dir.is_dir():
                continue
            import re as _re

            fuente = _re.sub(r'_contig_\d+$', '', sample_dir.name)
            for f in sample_dir.glob("plasmid_*.fasta"):
                stem = f"{_safe_name(fuente)}__mob_recon__{f.stem}"
                dest = todos_plasmidos / f"{stem}.fasta"
                if not dest.exists():
                    shutil.copy(f, dest)


        replicon_gbk = extraer_contigs_con_replicon([gbk_out], home, n_workers)
        for f in replicon_gbk.glob("*.fasta"):
            dest = todos_plasmidos / _safe_name(f.name)
            if not dest.exists():
                shutil.copy(f, dest)

        platon_gbk_out = run_platon(gbk_out, home, platon_db_gbk,
                                     threads=multiprocessing.cpu_count(),
                                     n_workers=n_workers)
        import re as _re
        for sample_dir in platon_gbk_out.iterdir():
            if not sample_dir.is_dir():
                continue
            for pf in sample_dir.glob("*.plasmid.fasta"):
                # sample_dir.name = LMG23361_contig_6
                # Extraer fuente: quitar _contig_N
                fuente = _re.sub(r'_contig_\d+$', '', sample_dir.name)
                contig_part = sample_dir.name
                stem = f"{_safe_name(fuente)}__platon__{_safe_name(contig_part)}-platon-plasmid"
                dest = todos_plasmidos / f"{stem}.fasta"
                if not dest.exists():
                    shutil.copy(pf, dest)

        genomad_db_path = conf.get("GENOMAD_DB") or pedir_db(
            conf, "GENOMAD_DB", "📁 Path to the geNomad database: ")
        genomad_gbk_out = run_genomad_dir(gbk_out, home, Path(genomad_db_path), n_workers)
        collect_genomad_plasmids(genomad_gbk_out, todos_plasmidos)

#BLAST ALL VS ALL
    unique_dir = home / "fastas-unicos"
    if unique_dir.exists() and any(f for f in unique_dir.glob("*.fasta")
                                   if "__" in f.stem):
        n = len([f for f in unique_dir.glob("*.fasta") if "__" in f.stem])
        print(f"\n   [checkpoint] fastas-unicos/ already exists ({n} FASTAs), skipping BLAST.")
    else:
        n_todos = len(list(todos_plasmidos.glob("*.fasta")))
        if n_todos == 0:
            print("\n   ⚠  No reconstructed plasmids found in todos_plasmidos/")
            unique_dir.mkdir(parents=True, exist_ok=True)
        else:
            check_exe("blastn")
            print(f"\n── BLAST intra-sample ({n_todos} FASTAs in todos_plasmidos/) ──")
            unique_dir = blast_all_vs_all([todos_plasmidos], home)

    blast_tsv_path = home / "resultados-blast-resumido.tsv"
    genomad_confirmado: set[str] = set()
    if blast_tsv_path.exists():
        try:
            with open(blast_tsv_path) as fh:
                for line in fh:
                    if line.startswith("sample") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 3:
                        continue
                    q = Path(parts[1]).stem
                    s = Path(parts[2]).stem
                    if "__genomad__" in q and "__genomad__" not in s:
                        genomad_confirmado.add(q)
                    if "__genomad__" in s and "__genomad__" not in q:
                        genomad_confirmado.add(s)
        except Exception as e:
            print(f"   ⚠  Error reading summarized BLAST: {e}")

    plasmidos_analizados = mkdir(home / "plasmidos_analizados")
    n_backup = 0
    n_genomad_filtrado = 0
    nombres_actuales = set()
    for f in sorted(unique_dir.glob("*.fasta")):
        if "__" not in f.stem:
            continue
        # Filtrar genomad-only (no confirmado por otro programa)
        if "__genomad__" in f.stem and f.stem not in genomad_confirmado:
            n_genomad_filtrado += 1
            # eliminar de fastas-unicos/
            f.unlink(missing_ok=True)
            continue
        nombres_actuales.add(f.name)
        dest_bk = plasmidos_analizados / f.name
        if not dest_bk.exists():
            shutil.copy(f, dest_bk)
            n_backup += 1

    # Sincronizar: eliminar de plasmidos_analizados/ los archivos que ya
    # no existen en fastas-unicos/ 
    n_obsoletos = 0
    for f in list(plasmidos_analizados.glob("*.fasta")):
        if f.name not in nombres_actuales:
            f.unlink()
            n_obsoletos += 1
    if n_obsoletos:
        print(f"   [sync] {n_obsoletos} stale backup(s) removed from "
              f"plasmidos_analizados/ (no longer present in fastas-unicos/)")

    if n_genomad_filtrado:
        print(f"   [filter] {n_genomad_filtrado} genomad-only plasmids removed")
    if n_backup:
        print(f"✔  Backup → plasmidos_analizados/ ({n_backup} FASTAs)")

    # ─── 8. Downstream pesado ─────────────────────────────────────────
    print("\n── Downstream analysis ───────────────────────────────────")
    print("   PlasmidFinder · PlasFlow · MOB-typer · Abricate")
    print("   Bakta · Panaroo · EggNOG-mapper · tANI")
    check_exe("plasmidfinder.py", "perl", "git")
    check_env("eggnog_env",   "emapper.py")
    check_env("mob_env",      "mob_recon")
    check_env("panaroo_env",  "panaroo")
    check_env("bakta_env",    "bakta")
    check_env("abricate_env", "abricate")
    check_env("genomad_env",  "genomad")
    check_env("tani_env",     "Rscript")
    bakta_db = pedir_db(conf, "BAKTA_DB",
                         "📁 Path to the Bakta database: ")


    run_downstream_fase1_paralelo(unique_dir, home, n_workers=n_workers)
    run_abricate(unique_dir, home, n_workers=n_workers)

    # PICCIS Score 
    print("\n── PICCIS plasmid reliability Score ───────────────────────")
    score_path = home / "piccis_score.tsv"
    if score_path.exists() and score_path.stat().st_size > 0:
        print(f"   [checkpoint] piccis_score.tsv already exists, skipping.")
        df_score = pd.read_csv(score_path, sep="\t")
    else:
        df_score = calcular_piccis_score(home, plasmidos_analizados)


    det_cols = [c for c in df_score.columns
                if c.startswith("det_") and c != "det_genomad"]
    if det_cols:
        mask_genomad     = df_score.get("det_genomad", pd.Series(0)) == 1
        mask_others_zero = df_score[det_cols].fillna(0).eq(0).all(axis=1)
        genomad_only_ids = set(df_score.loc[mask_genomad & mask_others_zero, "ID"])
        if genomad_only_ids:
            for fid in genomad_only_ids:
                for carpeta in [unique_dir, plasmidos_analizados]:
                    f = carpeta / f"{fid}.fasta"
                    if f.exists():
                        f.unlink()
            print(f"   [genomad-only filter] {len(genomad_only_ids)} sequences removed")

            score_path.unlink(missing_ok=True)
            df_score = calcular_piccis_score(home, plasmidos_analizados)


    print("\n── Unified table ──────────────────────────────────────────")

    table_path = home / "tabla_unificada.tsv"
    if meta_file and table_path.exists():
        table_path.unlink()

    if table_path.exists() and table_path.stat().st_size > 0:
        print(f"   [checkpoint] tabla_unificada.tsv already exists, skipping.")
        df_meta = pd.read_csv(table_path, sep="\t")
    else:
        table_path = build_unified_table(home, plasmidos_analizados,
                                         add_meta=False,
                                         metadata_file=meta_file)
        df_meta = pd.read_csv(table_path, sep="\t")
        df_meta = df_meta.merge(
            df_score[["ID", "piccis_score", "confiabilidad",
                      "n_detectan", "n_herramientas"]],
            on="ID", how="left"
        )
        df_meta.to_csv(table_path, sep="\t", index=False)
        print(f"✔  Unified table → {table_path}")

    bakta_out = run_bakta(unique_dir, home, bakta_db, n_workers=n_workers)
    run_panaroo(bakta_out, home, n_workers=n_workers)
    eggnog_db = conf.get("EGGNOG_DB") or None
    if eggnog_db:
        eggnog_db = Path(eggnog_db)
        print(f"   [EggNOG] using local DB: {eggnog_db}")
    else:
        print(f"   [EggNOG] local DB not configured → remote mode (requires internet)")

    from concurrent.futures import ProcessPoolExecutor as _PPE
    with _PPE(max_workers=2) as pool:
        fut_egg   = pool.submit(run_eggnog, bakta_out, home, n_workers, eggnog_db)
        fut_tani  = pool.submit(run_tani, unique_dir, home, n_workers)
        fut_egg.result()
        fut_tani.result()

    # Gráficos 
    print("\n── Generating plots ───────────────────────────────────────")

    grafico_piccis_score(df_score, gdir)
    grafico_ambiente(df_meta, gdir)
    grafico_descripcion(df_meta, gdir)
    grafico_tamano(df_meta, gdir)
    grafico_tamano_gc(df_meta, gdir)
    grafico_gc_ambiente(df_meta, gdir)
    grafico_mobility_torta(df_meta, gdir)
    grafico_tipos_movilizables(df_meta, gdir)
    grafico_tipos_conjugativos(df_meta, gdir)

    shp = obtener_shapefile(home)
    if shp is not None:
        grafico_planisferio(df_meta, shp, gdir)

    grafico_cog(home, df_meta, gdir)

    df_clusters = None
    if (home / "panaroo_out").exists():
        df_clusters = grafico_svd_kmeans(home, gdir)
    elif (home / "componentes_y_clusters.csv").exists():
        df_clusters = pd.read_csv(home / "componentes_y_clusters.csv", index_col=0)

    grafico_arbol_heatmap(home, df_meta, df_clusters, gdir)

    # Tablas de resumen 
    generar_tablas_resumen(df_meta, mkdir(home / "tablas_resumen"))

    print(f"""
╔══════════════════════════════════════════════════════════════════════╗
║  ✔  PICCIS v2.0 completed                                            ║
║     Results    : {str(home):<52}║
║     Plots      : {str(gdir):<52}║
╚══════════════════════════════════════════════════════════════════════╝
""")


if __name__ == "__main__":
    main()
