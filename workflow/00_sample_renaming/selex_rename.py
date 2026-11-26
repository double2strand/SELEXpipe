#!/usr/bin/env python3
"""
selex_rename.py

One-pass SELEX renamer: Sanger sample-ID folders/FASTQs -> final cycle layout.

Combines the two earlier scripts:
  1) map renamer    : <base-dir>/<SangerID>  -> sample name from the map
  2) cycle renamer  : C1HTS200SOX2R1         -> SOX2R1/SOX2R1c1
so it goes straight from the map IDs to the final layout:
  <base-dir>/8030STDY17011914/8030STDY17011914_L1_R1.fastq.gz
      -> <base-dir>/HTS/SOX2R1/SOX2R1c1/SOX2R1c1_L1_R1.fastq.gz       (default)
      -> <base-dir>/SOX2R1/SOX2R1c1/SOX2R1c1_L1_R1.fastq.gz           (--no-assay-subdir)

Map file (--map): two whitespace-separated columns, TAG and NAME, one per line.
Blank lines and lines starting with '#' are ignored; unparseable lines and
duplicate tags are warned about and skipped (first occurrence kept).

Name rules (as in the map renamer):
  * a trailing lowercase 'hts' is stripped   (MAFG_TDY1c4hts -> MAFG_TDY1c4)
  * characters outside [A-Za-z0-9._-] -> '_', runs of '_' collapsed, '_' stripped
Then (as in the cycle renamer) the name is parsed as
  ^(c\\d+)_?<ASSAY>_?(.+)$   (case-insensitive; ASSAY from --assay-pattern)
and the final folder name is <root><cycle> with the cycle lower-cased.
Names that do not parse (e.g. MAFG_TDY1c4) are used as-is, in place.
Duplicate names in the map get unique final names in map order (_2, _3, ...).

Placement (cycle renamer's rule, plus the assay level):
  * parent already named <root>  -> parent/<root><cycle>
  * otherwise                    -> parent/<root>/<root><cycle>
  * --assay-subdir (default on): parent/<ASSAY>/<root>/<root><cycle>, where ASSAY
    is the assay token without trailing digits (HTS200 -> HTS). An existing HTS/
    or NCAP/ folder (any case) is reused, and nothing is added when the parent
    is already HTS/ or NCAP/. This matches the HTS/SOX2R1/SOX2R1c1 layout and
    keeps HTS and NCAP folders of the same sample/cycle from colliding.

Files renamed inside each folder (recursive by default):
  * FASTQs named <TAG>_... ending .fastq.gz / .fq.gz -> <final>_...  (tag matched
    case-insensitively; other files are left alone)
  * any file whose name contains the old map/folder name (step-1 output or a
    cycle-renamer input folder) -> that part replaced with the final name

Also still supports trees that are half done: folders already named by the map
or any folder matching the cycle pattern anywhere under --base-dir, and leftover
tag-named FASTQs inside an already-final folder.

Safety: everything is planned and checked before anything is touched (targets
on disk, two items -> same target, folders blocked by files, file-name clashes).
On any clash nothing changes and the script exits with code 2, unless
--skip-existing is given (clashing items are skipped, the first in map/sorted
order wins). Never overwrites; re-runs skip what is already done.
Requires Python 3.9+.
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

DEFAULT_ASSAY_PATTERNS = ("HTS200", "NCAP200")
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz")


# ----------------------------------------------------------------- map / names
def parse_mapping_file(mapping_path: Path) -> list[tuple[str, str]]:
    """Return ordered (tag, raw_name) pairs; warn and skip bad/duplicate lines."""
    items: list[tuple[str, str]] = []
    seen: set[str] = set()
    with open(mapping_path, encoding="utf-8-sig") as fh:
        for i, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line, maxsplit=1)
            if len(parts) != 2:
                print(f"[WARN] Skipping line {i} (cannot parse into 2 columns): {line}",
                      file=sys.stderr)
                continue
            tag, raw_name = parts[0], parts[1].strip()
            if tag in seen:
                print(f"[WARN] Duplicate tag on line {i}: {tag} (keeping first occurrence only)",
                      file=sys.stderr)
                continue
            seen.add(tag)
            items.append((tag, raw_name))
    return items


def core_name(raw_name: str) -> str:
    """Strip a trailing lowercase 'hts' (MAFG_TDY1c4hts -> MAFG_TDY1c4)."""
    return re.sub(r"hts$", "", raw_name)


def sanitize_name(name: str) -> str:
    name = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    name = re.sub(r"_+", "_", name)
    return name.strip("_")


def parse_folder_name(name: str, assay_patterns) -> tuple[str, str, str] | None:
    """C1HTS200SOX2R1 -> ('SOX2R1c1', 'SOX2R1', 'HTS200'), or None."""
    for p in assay_patterns:
        m = re.match(rf"^(c\d+)_?{re.escape(p)}_?(.+)$", name, re.IGNORECASE)
        if m:
            cycle, root_name = m.group(1).lower(), m.group(2)
            return f"{root_name}{cycle}", root_name, p
    return None


def assay_label(pattern: str) -> str:
    return (re.sub(r"\d+$", "", pattern) or pattern).upper()


class Entry:
    """One map line resolved to its final name."""

    def __init__(self, tag, raw, name, new, root, label):
        self.tag, self.raw, self.name = tag, raw, name   # name = sanitised map name
        self.new, self.root, self.label = new, root, label


def resolve_entries(mapping_items, assay_patterns) -> list[Entry]:
    """Final names in map order; duplicate map names get _2, _3, ... (as before)."""
    entries, counts, used = [], defaultdict(int), set()
    for tag, raw in mapping_items:
        name = sanitize_name(core_name(raw))
        parsed = parse_folder_name(name, assay_patterns)
        new, root, label = (parsed[0], parsed[1], assay_label(parsed[2])) if parsed \
            else (name, None, "")
        counts[name] += 1
        if counts[name] > 1:
            k = counts[name]
            while (root, label, f"{new}_{k}") in used:
                k += 1
            new = f"{new}_{k}"
        used.add((root, label, new))
        entries.append(Entry(tag, raw, name, new, root, label))
    return entries


# ----------------------------------------------------------------- planning
class Op:
    """One planned folder move plus the file renames inside it."""

    def __init__(self, kind, src, rel, new, label, tag, olds):
        self.kind, self.src, self.rel, self.new = kind, src, rel, new
        self.label, self.tag, self.olds = label, tag, olds
        self.dest: Path | None = None
        self.files: list[tuple[Path, str]] = []   # (path relative to src, new name)
        self.skip = False


def word_rx(olds):
    """Case-insensitive regex for any of olds, not glued to letters/digits."""
    olds = sorted({o for o in olds if o}, key=len, reverse=True)
    if not olds:
        return None
    alt = "|".join(re.escape(o) for o in olds)
    return re.compile(rf"(?<![A-Za-z0-9])(?:{alt})(?![A-Za-z0-9])", re.IGNORECASE)


def new_file_name(name: str, tag: str | None, rx, new: str) -> str:
    """Tag-prefixed FASTQs (map renamer rule), else old-stem replacement."""
    if tag and name.upper().startswith(tag.upper() + "_") and name.lower().endswith(FASTQ_SUFFIXES):
        return new + name[len(tag):]
    if rx is not None:
        return rx.sub(new, name)
    return name


_assay_dir_cache: dict[tuple[Path, str], str] = {}


def assay_dir_name(parent: Path, label: str) -> str:
    """Reuse an existing HTS/ or NCAP/ folder (any case) under parent."""
    key = (parent, label)
    if key not in _assay_dir_cache:
        name = label
        if parent.is_dir():
            for c in sorted(parent.iterdir()):
                if c.is_dir() and c.name.upper() == label:
                    name = c.name
                    break
        _assay_dir_cache[key] = name
    return _assay_dir_cache[key]


def destination_for(parent: Path, root_name, new, label, assay_subdir) -> tuple[str, ...]:
    """Destination relative to parent (cycle renamer's rule + assay level)."""
    if root_name is None:                 # name without cycle/assay: rename in place
        return (new,)
    if parent.name == root_name:
        return (new,)
    if assay_subdir and label:
        if parent.name.upper() == label:
            return (root_name, new)
        return (assay_dir_name(parent, label), root_name, new)
    return (root_name, new)


def build_plan(base, entries, assay_patterns, assay_subdir, recursive, ids_only, dry_run):
    ops: list[Op] = []
    stats = defaultdict(int)
    fixups: list[tuple[Path, Path]] = []
    children = {c.name.upper(): c for c in base.iterdir()} if base.is_dir() else {}
    tag_srcs = set()

    # 1) Sanger-ID folders, exactly as the map renamer looked for them
    for e in entries:
        src = base / e.tag
        if not src.is_dir() and e.tag.upper() in children and children[e.tag.upper()].is_dir():
            src = children[e.tag.upper()]          # tolerate case differences
        rel = destination_for(base, e.root, e.new, e.label, assay_subdir)
        dst = base.joinpath(*rel)
        if not src.is_dir():
            if dst.is_dir():
                if dry_run:
                    print(f"[DRY] SKIP already-renamed folder: {dst.name}")
                stats["already"] += 1
                # finish leftovers: tag-named FASTQs still inside the final folder
                for f in sorted(dst.iterdir()):
                    if f.is_file():
                        nn = new_file_name(f.name, e.tag, None, e.new)
                        if nn != f.name:
                            fixups.append((f, f.with_name(nn)))
            else:
                print(f"[WARN] Folder not found, skipping: {src}")
                stats["not_found"] += 1
            continue
        tag_srcs.add(src)
        ops.append(Op("tag", src, rel, e.new, e.label, e.tag, {e.name}))

    # 2) folders already named by the map / matching the cycle pattern (cycle renamer)
    by_name = {e.name.upper(): e for e in entries}
    for d in sorted(p for p in base.rglob("*") if p.is_dir() and not p.is_symlink()):
        if d in tag_srcs or any(a in tag_srcs for a in d.parents):
            continue
        e = by_name.get(d.name.upper())
        if e is not None and e.root is not None:
            new, root_name, label, tag = e.new, e.root, e.label, e.tag
        elif ids_only:
            continue
        else:
            parsed = parse_folder_name(d.name, assay_patterns)
            if parsed is None:
                continue
            new, root_name, label, tag = parsed[0], parsed[1], assay_label(parsed[2]), None
        rel = destination_for(d.parent, root_name, new, label, assay_subdir)
        if d.parent.joinpath(*rel) == d:
            if dry_run:
                print(f"[DRY] SKIP already correct: {d}")
            stats["already"] += 1
            continue
        ops.append(Op("cycle", d, rel, new, label, tag, {d.name}))

    # final paths (children of moved folders follow their parent)
    by_src = {op.src: op for op in ops}
    memo: dict[Path, Path] = {}

    def final(p: Path) -> Path:
        if p not in memo:
            if p == base or base not in p.parents:
                memo[p] = p
            elif p in by_src:
                memo[p] = final(p.parent).joinpath(*by_src[p].rel)
            else:
                memo[p] = final(p.parent) / p.name
        return memo[p]

    for op in ops:
        op.dest = final(op.src)
        rx = word_rx(op.olds)
        it = op.src.rglob("*") if recursive else op.src.iterdir()
        for f in sorted(it):
            if not f.is_file():
                continue
            if any(a in by_src for a in f.parents if a != op.src and op.src in a.parents):
                continue                       # belongs to a nested planned folder
            op.files.append((f.relative_to(op.src), new_file_name(f.name, op.tag, rx, op.new)))

    mkdirs: dict[Path, Op] = {}
    for op in ops:
        pf = final(op.src.parent)
        for i in range(1, len(op.rel)):
            mkdirs.setdefault(pf.joinpath(*op.rel[:i]), op)
    return ops, fixups, mkdirs, final, stats


# ----------------------------------------------------------------- preflight
def preflight(base, ops, fixups, mkdirs, skip_existing):
    """Return list of clash messages; mark skipped ops and files."""
    msgs: list[str] = []
    skipped_files: set[Path] = set()
    tail = " (skipping)" if skip_existing else ""

    def warn_exists(op, dst):
        if op.kind == "tag":
            msgs.append(f"[WARN] Target folder already exists, cannot rename {op.src} -> {dst}{tail}")
        else:
            msgs.append(f"[WARN] Destination exists, cannot move {op.src} -> {dst}{tail}")

    groups: dict[Path, list[Op]] = defaultdict(list)
    for op in ops:                         # map order first, then sorted folders
        groups[op.dest].append(op)
    for dst, group in groups.items():
        if dst.exists() or dst.is_symlink():
            for op in group:
                warn_exists(op, dst)
                op.skip = True
            continue
        if dst in mkdirs:
            for op in group:
                msgs.append(f"[WARN] Planned clash: {dst} is needed as a folder for "
                            f"{mkdirs[dst].src}, cannot move {op.src} there{tail}")
                op.skip = True
            continue
        if len(group) > 1:
            labels = sorted({op.label for op in group if op.label})
            note = ""
            if len(labels) > 1:
                note = (f"\n       ({' and '.join(labels)} of the same sample and cycle map to "
                        f"the same folder; drop --no-assay-subdir to keep assays apart)")
            for op in group[1:]:
                msgs.append(f"[WARN] Planned clash: {op.src} and {group[0].src} would both "
                            f"move to {dst}{tail}{note}")
                op.skip = True
    for d, op in mkdirs.items():
        if (d.exists() or d.is_symlink()) and not d.is_dir():
            for o in ops:
                if not o.skip and d in o.dest.parents:
                    msgs.append(f"[WARN] Cannot create folder {d} (a file of that name exists), "
                                f"cannot move {o.src} -> {o.dest}{tail}")
                    o.skip = True

    for op in ops:
        if op.skip:
            continue
        seen: dict[Path, list[Path]] = defaultdict(list)
        for rel, new in op.files:
            seen[rel.parent / new].append(rel)
        for tgt, rels in seen.items():
            if len(rels) > 1:
                keep = tgt if tgt in rels else rels[0]     # an existing name always wins
                for r in rels:
                    if r != keep:
                        msgs.append(f"[WARN] File target already exists, cannot rename "
                                    f"{op.src / r} -> {op.dest / tgt}{tail}")
                        skipped_files.add(op.src / r)
    tgts: dict[Path, list[Path]] = defaultdict(list)
    for f, tgt in fixups:
        tgts[tgt].append(f)
    for tgt, fs in tgts.items():
        if tgt.exists() or tgt.is_symlink():
            bad = fs
        else:
            bad = fs[1:]
        for f in bad:
            msgs.append(f"[WARN] File target already exists, cannot rename {f} -> {tgt}{tail}")
            skipped_files.add(f)
    return msgs, skipped_files


# ----------------------------------------------------------------- actions
def rel_to(p: Path, base: Path) -> str:
    try:
        return str(p.relative_to(base))
    except ValueError:
        return str(p)


def safe_rename(src: Path, dst: Path, dry_run: bool, base: Path, final_dst: Path | None = None) -> None:
    shown = rel_to(final_dst or dst, base)
    print(f"{'[DRY] ' if dry_run else '      '}rename {rel_to(src, base)} -> {shown}")
    if dry_run:
        return
    if dst.exists() or dst.is_symlink():
        print(f"[ERROR] Target appeared during run, stopping: {dst}", file=sys.stderr)
        sys.exit(2)
    src.rename(dst)


def safe_mkdir(path: Path, dry_run: bool, base: Path) -> None:
    print(f"{'[DRY] ' if dry_run else '      '}mkdir {rel_to(path, base)}")
    if not dry_run:
        path.mkdir(parents=True, exist_ok=True)


def safe_move(src: Path, dst: Path, dry_run: bool, base: Path) -> None:
    print(f"{'[DRY] ' if dry_run else '      '}move {rel_to(src, base)} -> {rel_to(dst, base)}")
    if dry_run:
        return
    if dst.exists() or dst.is_symlink():
        print(f"[ERROR] Target appeared during run, stopping: {dst}", file=sys.stderr)
        sys.exit(2)
    shutil.move(str(src), str(dst))


# ----------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser(
        description="Rename Sanger-ID SELEX folders and FASTQs straight to the final "
                    "<sample>/<sample><cycle> layout using a TAG -> NAME map file.")
    ap.add_argument("--base-dir", required=True, type=Path,
                    help="Folder containing the Sanger-ID folders (searched recursively "
                         "for cycle folders).")
    ap.add_argument("--map", dest="map_path", required=True, type=Path,
                    help="Map file with two columns: TAG NAME "
                         "(e.g. 8030STDY17011914 C1HTS200SOX2R1).")
    ap.add_argument("--assay-pattern", action="append", dest="assay_patterns", default=None,
                    help="Middle pattern between cycle and sample name. May be given multiple "
                         f"times. Default: {', '.join(DEFAULT_ASSAY_PATTERNS)}.")
    ap.add_argument("--assay-subdir", action=argparse.BooleanOptionalAction, default=True,
                    help="Put outputs under HTS/ and NCAP/ (reusing existing ones) so the two "
                         "assays never collide. Default: on. Use --no-assay-subdir for "
                         "<sample>/<sample><cycle> directly under the parent.")
    ap.add_argument("--no-recursive-file-rename", action="store_true",
                    help="Only rename files directly inside each folder, not in sub-folders.")
    ap.add_argument("--ids-only", action="store_true",
                    help="Only act on Sanger-ID folders and folders named by the map, not on "
                         "other folders that merely match the cycle pattern.")
    ap.add_argument("--skip-existing", action="store_true",
                    help="Skip items whose target already exists or clashes instead of aborting.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Show what would be done without changing anything.")
    args = ap.parse_args()
    sys.stdout.reconfigure(line_buffering=True)     # keep [WARN]/[INFO] in order

    base_dir: Path = args.base_dir.resolve()
    map_path: Path = args.map_path
    if not base_dir.is_dir():
        print(f"[ERROR] Base directory not found: {base_dir}", file=sys.stderr)
        sys.exit(1)
    if not map_path.is_file():
        print(f"[ERROR] Map file not found: {map_path}", file=sys.stderr)
        sys.exit(1)
    assay_patterns = tuple(args.assay_patterns) if args.assay_patterns else DEFAULT_ASSAY_PATTERNS
    recursive = not args.no_recursive_file_rename
    dry = args.dry_run

    print(f"[INFO] Base directory: {base_dir}")
    print(f"[INFO] Map file: {map_path}")
    mapping_items = parse_mapping_file(map_path)
    print(f"[INFO] Mappings loaded: {len(mapping_items)}")
    print(f"[INFO] Assay sub-folders: {'on' if args.assay_subdir else 'off'}"
          f"{'  (DRY RUN)' if dry else ''}")
    print()

    entries = resolve_entries(mapping_items, assay_patterns)
    ops, fixups, mkdirs, final, stats = build_plan(
        base_dir, entries, assay_patterns, args.assay_subdir, recursive, args.ids_only, dry)
    msgs, skipped_files = preflight(base_dir, ops, fixups, mkdirs, args.skip_existing)
    if msgs:
        out = sys.stdout if args.skip_existing else sys.stderr
        for m in msgs:
            print(m, file=out)
        if not args.skip_existing:
            print(f"\n[ERROR] {len(msgs)} clash(es) found before starting; nothing was changed. "
                  f"Fix them or re-run with --skip-existing.", file=sys.stderr)
            sys.exit(2)
        print()

    folders_done = files_done = created = skipped_existing = files_skipped = 0
    made: set[Path] = set()
    order = sorted(ops, key=lambda o: -len(o.src.parts))      # deepest first (stable)
    for op in order:
        if op.skip:
            skipped_existing += 1
            continue
        dest_now = op.src.parent.joinpath(*op.rel)               # ancestors not moved yet
        old = op.tag if op.kind == "tag" else op.src.name
        if op.kind == "tag":
            print(f"[INFO] Processing {old}  ->  {rel_to(op.dest, base_dir)}")
        else:
            print(f"[INFO] {old}  ->  {rel_to(op.dest, base_dir)}")
        # FASTQs/files first, then the folder (map renamer's order)
        for rel, new in op.files:
            if new == rel.name:
                continue
            if op.src / rel in skipped_files:
                files_skipped += 1
                continue
            cur = op.src / rel
            safe_rename(cur, cur.with_name(new), dry, base_dir, op.dest / rel.parent / new)
            files_done += 1
        for i in range(1, len(op.rel)):
            d = op.src.parent.joinpath(*op.rel[:i])
            if d not in made and not d.exists():
                safe_mkdir(d, dry, base_dir)
                made.add(d)
                created += 1
        safe_move(op.src, dest_now, dry, base_dir)
        folders_done += 1
        print()
    for f, tgt in fixups:
        if f in skipped_files:
            files_skipped += 1
            continue
        print(f"[INFO] Finishing already-renamed folder {rel_to(f.parent, base_dir)}")
        safe_rename(f, tgt, dry, base_dir)
        files_done += 1

    print()
    print("========== SUMMARY ==========")
    print(f"Folders moved/renamed: {folders_done}")
    print(f"Files renamed:         {files_done}")
    print(f"Folders created:       {created}")
    if stats["already"]:
        print(f"Skipped already-renamed (idempotent): {stats['already']}")
    if skipped_existing:
        print(f"Skipped (existing targets via --skip-existing): {skipped_existing}")
    if files_skipped:
        print(f"Files skipped (existing targets via --skip-existing): {files_skipped}")
    if stats["not_found"]:
        print(f"Folders not found:     {stats['not_found']}")
    print("=============================")
    if dry:
        print("[DRY] Dry run only - nothing was changed.")


if __name__ == "__main__":
    main()
