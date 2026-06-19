#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re, gzip, random, copy, hashlib
from collections import OrderedDict, defaultdict
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple, Sequence, Optional
import pandas as pd

# ---------- FASTA I/O ----------
def read_fasta(path: str) -> OrderedDict:
    op = gzip.open if path.endswith(".gz") else open
    fa = OrderedDict()
    name, buff = None, []
    with op(path, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith(">"):
                if name:
                    fa[name] = "".join(buff).upper()
                name = line[1:].strip().split()[0]
                buff = []
            else:
                buff.append(line.strip())
        if name:
            fa[name] = "".join(buff).upper()
    return fa

def write_fasta(path: str, fa: Dict[str,str], width: int = 60):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "wt") as out:
        for name, seq in fa.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i+width] + "\n")

def subset_fasta(fa: Dict[str,str], keep_names: Sequence[str]) -> OrderedDict:
    return OrderedDict((k, v) for k, v in fa.items() if k in keep_names)

# ---------- RC ----------
RCMAP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s: str) -> str:
    return s.translate(RCMAP)[::-1]

# ---------- HPV circular helpers ----------
def normalize(n: int, L: int) -> int:
    """Normalize any integer to 1..L for circular coordinates."""
    return ((n - 1) % L) + 1

def extract_circular(seq: str, a: int, b: int) -> str:
    """Extract [a,b] 1-based inclusive on a circular sequence (wrap if a>b)."""
    L = len(seq)
    a = normalize(a, L)
    b = normalize(b, L)
    a0, b0 = a-1, b-1
    if a <= b:
        return seq[a0:b0+1]
    else:
        return seq[a0:] + seq[:b0+1]

def build_hpv_unit(hpv_seq: str, paths: List[Tuple[int,int]], hpv_strand: str, copies: int = 1) -> str:
    segs = [extract_circular(hpv_seq, a, b) for (a,b) in paths]
    unit = "".join(segs)
    unit = rc(unit) if hpv_strand == "-" else unit
    return unit * copies

# ---------- composite insert (hpv / host / literal) ----------
def build_composite_insert(hpv_seq: str, blocks: List[Dict], host_ref: Dict[str,str]) -> Tuple[str, str]:
    """
    Build insert sequence by concatenating blocks.
    Returns: (sequence, description_string)
    """
    pieces, descs = [], []
    for blk in blocks:
        t = blk.get("type")
        if t == "hpv":
            paths  = blk["paths"]
            strand = blk.get("strand", "+")
            copies = int(blk.get("copies", 1))
            unit   = build_hpv_unit(hpv_seq, paths, strand, copies=copies)
            pieces.append(unit)
            pstr = ",".join([f"{a}-{b}" for a,b in paths])
            rep  = f"x{copies}" if copies != 1 else ""
            descs.append(f"HPV:{pstr}({strand}){rep}")

        elif t == "host":
            chrom  = blk["chrom"]
            a = int(blk["start"])
            b = int(blk["end"])
            strand = blk.get("strand", "+")
            if chrom not in host_ref:
                raise ValueError(f"Host block chrom not loaded: {chrom}")
            ref = host_ref[chrom]
            if not (1 <= a <= len(ref) and 1 <= b <= len(ref)):
                raise ValueError(f"Host block coords out of range: {chrom}:{a}-{b} (len={len(ref)})")
            if a > b:
                raise ValueError(f"Host block expects linear coords a<=b, got {a}-{b}")
            sub = ref[a-1:b]
            sub = rc(sub) if strand == "-" else sub
            pieces.append(sub)
            descs.append(f"{chrom}:{a}-{b}({strand})")

        elif t == "literal":
            seq = blk["seq"].upper()
            pieces.append(seq)
            descs.append(f"LIT:{len(seq)}")

        else:
            raise ValueError(f"Unknown block type: {t}")

    return "".join(pieces), "+".join(descs)

# ---------- intervals & placement ----------
def overlaps(intervals: List[Tuple[int,int]], new_iv: Tuple[int,int]) -> bool:
    a,b = new_iv
    for x,y in intervals:
        if not (b < x or a > y):  # overlap
            return True
    return False

def read_avoid_bed(path: Optional[str]) -> Dict[str, List[Tuple[int,int]]]:
    """
    Read BED file of regions to avoid.

    BED is 0-based half-open: chrom start end.
    This function converts it to 1-based inclusive intervals
    to match the coordinate logic used in this script.
    """
    avoid = defaultdict(list)
    if not path:
        return avoid

    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split()
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start0 = int(fields[1])
            end0 = int(fields[2])
            if end0 <= start0:
                continue
            # BED [start0, end0) -> 1-based inclusive [start0+1, end0]
            avoid[chrom].append((start0 + 1, end0))

    return dict(avoid)

def site_overlaps_avoid(chrom: str, interval: Tuple[int,int], avoid_regions: Optional[Dict[str, List[Tuple[int,int]]]]) -> bool:
    if not avoid_regions:
        return False
    return overlaps(avoid_regions.get(chrom, []), interval)

def choose_random_site(chrom_sizes: Dict[str,int],
                       allow_chrs: List[str],
                       del_len: int,
                       occupied: Dict[str,List[Tuple[int,int]]],
                       rng: random.Random,
                       min_gap: int = 50000,
                       max_retry: int = 10000,
                       avoid_regions: Optional[Dict[str, List[Tuple[int,int]]]] = None) -> Tuple[str,int]:
    """
    Choose random insertion anchor `pos` in allowed chromosomes.
    We place a protected window [pos+1 - min_gap, pos+del_len + min_gap] to avoid nearby events.
    """
    weights = [(c, chrom_sizes[c]) for c in allow_chrs]
    tot = sum(w for _,w in weights)
    if tot <= 0:
        raise ValueError("Total chromosome length is zero in allow_chrs.")

    for _ in range(max_retry):
        r = rng.randint(1, tot)
        acc = 0; chrom = None
        for c,w in weights:
            acc += w
            if r <= acc:
                chrom = c
                break

        L = chrom_sizes[chrom]
        max_pos = max(1, L - del_len)
        pos = rng.randint(1, max_pos)

        core_start = pos + 1
        core_end   = pos + del_len if del_len > 0 else pos

        start = max(1, core_start - min_gap)
        end   = min(L, core_end + min_gap)
        iv = (start, end)

        if overlaps(occupied[chrom], iv):
            continue

        if site_overlaps_avoid(chrom, iv, avoid_regions):
            continue

        occupied[chrom].append(iv)
        return chrom, pos

    raise RuntimeError("Failed to place event without overlap")

# ---------- event schema ----------
@dataclass
class Event:
    event_id: str
    chrom: str
    host_pos: int
    del_len: int
    ins_len: int
    hpv_ref: str
    hpv_path: str
    hpv_strand: str
    copies: int
    source: str

# ---------- site parsing ----------
def parse_site(s: str) -> Tuple[str, int]:
    s = s.strip()
    m = re.match(r'^(chr[^\s:]+)\s*[:\s]\s*(\d+)$', s)
    if not m:
        raise ValueError(f"Bad site format: {s}")
    return m.group(1), int(m.group(2))

def normalize_pos_within_chr(pos: int, L: int) -> int:
    return ((pos - 1) % L) + 1

def pick_sites_from_list(cands: List[str],
                         allow_chrs: List[str],
                         chrom_sizes: Dict[str,int],
                         occupied: Dict[str,List[Tuple[int,int]]],
                         k: int,
                         del_len: int,
                         rng: random.Random,
                         min_gap: int = 0,
                         avoid_regions: Optional[Dict[str, List[Tuple[int,int]]]] = None) -> List[Tuple[str,int]]:
    """
    Pick sites from provided list (chr:pos), skipping overlaps.
    If min_gap>0, will also protect +/- min_gap around the deletion window to avoid clustering.
    """
    out: List[Tuple[str,int]] = []
    tried = set()
    order = list(cands)
    rng.shuffle(order)

    for s in order:
        if len(out) >= k:
            break
        if s in tried:
            continue
        tried.add(s)
        try:
            chrom, pos = parse_site(s)
        except Exception:
            continue
        if chrom not in allow_chrs:
            continue
        L = chrom_sizes.get(chrom)
        if not L:
            continue
        pos = normalize_pos_within_chr(pos, L)
        if pos + del_len > L:
            continue

        core_start = pos + 1
        core_end   = pos + del_len if del_len > 0 else pos
        start = max(1, core_start - min_gap)
        end   = min(L, core_end + min_gap)
        iv = (start, end)

        if overlaps(occupied[chrom], iv):
            continue

        if site_overlaps_avoid(chrom, iv, avoid_regions):
            continue

        occupied[chrom].append(iv)
        out.append((chrom, pos))

    return out

# ---------- per-event perturbation that preserves length & preserves repeats inside one event ----------
def _hpv_key(blk: Dict) -> Tuple:
    paths = tuple((int(a), int(b)) for a,b in blk["paths"])
    strand = blk.get("strand", "+")
    copies = int(blk.get("copies", 1))
    return ("hpv", strand, copies, paths)

def _host_key(blk: Dict) -> Tuple:
    chrom = blk["chrom"]
    a = int(blk["start"]); b = int(blk["end"])
    strand = blk.get("strand", "+")
    return ("host", chrom, a, b, strand)

HPV16_BIASED_REGIONS = [
    (1, 858, 6.0),        # E6/E7
    (866, 3853, 1.0),     # E1/E2
    (3854, 7906, 3.0)
]

def circular_len(a, b, hpv_len):
    a = normalize(a, hpv_len)
    b = normalize(b, hpv_len)

    if a <= b:
        return b - a + 1

    return hpv_len - a + 1 + b


def pick_hpv_segment_biased(hpv_len, seg_len, rng):
    starts = []
    weights = []

    for start, end, weight in HPV16_BIASED_REGIONS:
        for s in range(start, end + 1):
            starts.append(s)
            weights.append(weight)

    a = rng.choices(starts, weights=weights, k=1)[0]
    b = normalize(a + seg_len - 1, hpv_len)

    return a, b
def relocate_host_blocks_near_integration(blocks, event_chrom, event_pos, chrom_sizes,
                                          upstream_offset=1000):
    """
    For host blocks:
      - keep original block length
      - replace chrom with event_chrom
      - set start near upstream of integration site
      - identical host templates reuse same relocated interval within one event
    """
    blks = copy.deepcopy(blocks)
    host_map = {}

    for blk in blks:
        if blk.get("type") != "host":
            continue

        key = _host_key(blk)

        old_start = int(blk["start"])
        old_end = int(blk["end"])
        seg_len = old_end - old_start + 1

        if key not in host_map:
            L = chrom_sizes[event_chrom]

            new_start = event_pos - upstream_offset - seg_len + 1
            new_end = event_pos - upstream_offset

            if new_start < 1:
                new_start = event_pos + upstream_offset
                new_end = new_start + seg_len - 1

            if new_end > L:
                raise RuntimeError(
                    f"Cannot relocate host block near {event_chrom}:{event_pos}, "
                    f"seg_len={seg_len}, chr_len={L}"
                )

            host_map[key] = (event_chrom, new_start, new_end)

        new_chrom, new_start, new_end = host_map[key]
        blk["chrom"] = new_chrom
        blk["start"] = new_start
        blk["end"] = new_end

    return blks



def perturb_blocks_hpv_and_host_shared(
    blocks: List[Dict],
    hpv_len: int,
    chrom_sizes: Dict[str,int],
    rng: random.Random,
    hpv_shift_max: int = 2000,
    host_shift_max: int = 50000,
    min_hpv_shift: int = 200,
    min_host_shift: int = 2000,
    max_tries: int = 100
) -> List[Dict]:
    """
    Perturb blocks *per event*:
      - HPV blocks: circular shift all paths in that block by Delta (length preserved)
      - HOST blocks: linear shift interval by delta in same chrom (length preserved)
    IMPORTANT:
      - If the SAME block template repeats inside a single event (like UD2),
        we reuse the SAME Delta/delta so the repeated content stays identical within that event.
    """
    blks = copy.deepcopy(blocks)
    hpv_delta_map: Dict[Tuple, int] = {}
    host_delta_map: Dict[Tuple, int] = {}

    for blk in blks:
        t = blk.get("type")

        if t == "hpv":
            key = _hpv_key(blk)
            if key not in hpv_delta_map:
                if hpv_shift_max <= 0:
                    Delta = 0
                else:
                    Delta = rng.randint(-hpv_shift_max, hpv_shift_max)
                    if abs(Delta) < min_hpv_shift:
                        Delta = min_hpv_shift if Delta >= 0 else -min_hpv_shift
                hpv_delta_map[key] = Delta
            Delta = hpv_delta_map[key]

            new_paths = []

            for (a, b) in blk["paths"]:

                seg_len = circular_len(int(a), int(b), hpv_len)

                new_a, new_b = pick_hpv_segment_biased(
                    hpv_len,
                    seg_len,
                    rng
                )

                new_paths.append((new_a, new_b))

            blk["paths"] = new_paths

        elif t == "host":
            key = _host_key(blk)
            chrom = blk["chrom"]
            if chrom not in chrom_sizes:
                continue
            L = chrom_sizes[chrom]
            a = int(blk["start"]); b = int(blk["end"])
            if a > b:
                raise ValueError(f"Host block expects start<=end, got {a}-{b}")
            seg_len = b - a + 1

            if key not in host_delta_map:
                if host_shift_max <= 0:
                    host_delta_map[key] = 0
                else:
                    ok = False
                    for _ in range(max_tries):
                        delta = rng.randint(-host_shift_max, host_shift_max)
                        if abs(delta) < min_host_shift:
                            delta = min_host_shift if delta >= 0 else -min_host_shift
                        a2 = a + delta
                        b2 = a2 + seg_len - 1
                        if 1 <= a2 and b2 <= L:
                            host_delta_map[key] = delta
                            ok = True
                            break
                    if not ok:
                        host_delta_map[key] = 0

            delta = host_delta_map[key]
            a2 = a + delta
            b2 = a2 + seg_len - 1
            if 1 <= a2 and b2 <= L:
                blk["start"], blk["end"] = a2, b2
            # else keep original

        # literal: keep as-is

    return blks

def md5_seq(s: str) -> str:
    return hashlib.md5(s.encode("utf-8")).hexdigest()

# ---------- core ----------
def insert_events_from_site_lists(host_fa_path: str,
                                  hpv_fa_path: str,
                                  out_fa_path: str,
                                  site_lists: Dict[str, List[str]],
                                  allow_chrs: Tuple[str,...],
                                  recipes: Optional[List[Dict]] = None,
                                  out_prefix: str = "run_sitelist",
                                  seed: int = 20251004,
                                  extra_sites: Optional[List[str]] = None,
                                  output_only_modified: bool = True,
                                  target_total: int = 100,
                                  # NEW knobs:
                                  min_gap: int = 50000,
                                  ensure_unique_insert: bool = True,
                                  unique_retry: int = 20,
                                  avoid_bed: Optional[str] = None):
    """
    Build mutated host reference with multiple HPV integration events.
    Main differences vs your original:
      - per-event perturbation (HPV shift + host shift), length-preserving
      - repeated blocks inside an event stay identical (shared shift)
      - optional md5 de-dup to avoid identical inserts across events
    """
    rng = random.Random(seed)

    # load & subset host / hpv
    host = read_fasta(host_fa_path)
    host = subset_fasta(host, [c for c in allow_chrs if c in host])
    if not host:
        raise ValueError("No allowed chromosomes found in host FASTA.")

    hpv = read_fasta(hpv_fa_path)
    if len(hpv) != 1:
        raise ValueError("HPV FASTA should contain exactly one sequence.")
    hpv_name, hpv_seq = list(hpv.items())[0]

    chrom_sizes = {c: len(s) for c,s in host.items()}
    mutable = OrderedDict((c, list(seq)) for c,seq in host.items())
    occupied = {c: [] for c in host.keys()}
    avoid_regions = read_avoid_bed(avoid_bed) if avoid_bed else None

    if avoid_regions:
        n_avoid = sum(len(v) for v in avoid_regions.values())
        print(f"Avoid BED loaded: {avoid_bed} ({n_avoid} intervals)")

    # default recipes if none
    if recipes is None:
        recipes = [
            {"name":"C", "paths":[(1471,7906),(1,2727)], "strand": "+", "del_len": 33063, "copies": 1},
            {"name":"D", "paths":[(5983,7906),(1,7906),(1,7906),(1,3406)], "strand": "+", "del_len": 342, "copies": 1},
        ]

    def recipe_max_del(recs: List[Dict]) -> int:
        return max(int(r.get("del_len", 0)) for r in recs) if recs else 0

    max_dlen = recipe_max_del(recipes)

    chosen_sites: List[Tuple[str,int,str]] = []

    # 1) fixed lists
    for tag, lst in site_lists.items():
        if not lst:
            continue
        picks = pick_sites_from_list(
            lst,
            list(allow_chrs),
            chrom_sizes,
            occupied,
            k=len(lst),
            del_len=max_dlen,
            rng=rng,
            min_gap=min_gap,
            avoid_regions=avoid_regions
        )
        for chrom, pos in picks:
            chosen_sites.append((chrom, pos, tag))

    # 2) extra_sites
    if extra_sites:
        extra_picks = pick_sites_from_list(
            extra_sites,
            list(allow_chrs),
            chrom_sizes,
            occupied,
            k=len(extra_sites),
            del_len=max_dlen,
            rng=rng,
            min_gap=min_gap,
            avoid_regions=avoid_regions
        )
        for chrom, pos in extra_picks:
            chosen_sites.append((chrom, pos, "EXTRA"))

    # 3) random sites until target_total
    while len(chosen_sites) < target_total:
        chrom, pos = choose_random_site(
            chrom_sizes=chrom_sizes,
            allow_chrs=list(allow_chrs),
            del_len=max_dlen,
            occupied=occupied,
            rng=rng,
            min_gap=min_gap,
            avoid_regions=avoid_regions
        )
        chosen_sites.append((chrom, pos, f"RAND{len(chosen_sites)+1}"))

    # Apply from RIGHT to LEFT per chromosome to keep `pos` referring to original coordinates
    events: List[Event] = []
    bedpe_rows = []
    recipe_idx = 0
    modified_chroms = set()

    chr_order = {c:i for i, c in enumerate(allow_chrs)}
    chosen_sites = sorted(chosen_sites, key=lambda x: (chr_order.get(x[0], 10**9), -x[1]))

    seen_insert_md5 = set()

    for idx, (chrom, pos, src) in enumerate(chosen_sites, start=1):
        r = recipes[recipe_idx % len(recipes)]
        recipe_idx += 1

        # build insert with per-event perturbation if blocks exist
        if "blocks" in r:
            # allow per-recipe overrides
            hpv_shift_max   = int(r.get("hpv_shift_max", 2000))
            host_shift_max  = int(r.get("host_shift_max", 50000))
            min_hpv_shift   = int(r.get("min_hpv_shift", 200))
            min_host_shift  = int(r.get("min_host_shift", 2000))

            # optional: retry until insert becomes unique
            for attempt in range(unique_retry if ensure_unique_insert else 1):
                # pert_blocks = perturb_blocks_hpv_and_host_shared(
                #     blocks=r["blocks"],
                #     hpv_len=len(hpv_seq),
                #     chrom_sizes=chrom_sizes,
                #     rng=rng,
                #     hpv_shift_max=hpv_shift_max,
                #     host_shift_max=host_shift_max,
                #     min_hpv_shift=min_hpv_shift,
                #     min_host_shift=min_host_shift
                # )
                pert_blocks = perturb_blocks_hpv_and_host_shared(
                    blocks=r["blocks"],
                    hpv_len=len(hpv_seq),
                    chrom_sizes=chrom_sizes,
                    rng=rng,
                    hpv_shift_max=hpv_shift_max,
                    host_shift_max=0,  
                    min_hpv_shift=min_hpv_shift,
                    min_host_shift=0
                )
                pert_blocks = relocate_host_blocks_near_integration(
                    blocks=pert_blocks,
                    event_chrom=chrom,
                    event_pos=pos,
                    chrom_sizes=chrom_sizes,
                    upstream_offset=1000
                )
                insert_seq, hpv_path_desc = build_composite_insert(hpv_seq, pert_blocks, host_ref=host)
                h = md5_seq(insert_seq)
                if (not ensure_unique_insert) or (h not in seen_insert_md5):
                    seen_insert_md5.add(h)
                    break
            else:
                # after retries still duplicated; accept but warn via hash table size
                seen_insert_md5.add(h)

            dlen   = int(r.get("del_len", 0))
            copies = 1
            strand = "."
        else:
            # legacy recipe using HPV paths only (optionally you could also perturb these similarly)
            paths   = r["paths"]
            strand  = r.get("strand", "+")
            dlen    = int(r.get("del_len", 0))
            copies  = int(r.get("copies", 1))
            insert_seq = build_hpv_unit(hpv_seq, paths, strand, copies=copies)
            hpv_path_desc = ";".join([f"{a}-{b}" for a,b in paths]) + f"({strand})"

        ins_len = len(insert_seq)
        L = chrom_sizes[chrom]
        if pos + dlen > L:
            raise RuntimeError(f"Chosen site {chrom}:{pos} with del_len={dlen} exceeds chromosome length {L}.")

        # apply edit: delete [pos, pos+dlen) and insert insert_seq
        seq = mutable[chrom]
        s0 = pos
        e0 = pos + dlen
        new_seq = seq[:s0] + list(insert_seq) + seq[e0:]
        mutable[chrom] = new_seq
        modified_chroms.add(chrom)

        ev_id = f"E{idx}_{src}_{r.get('name', f'R{recipe_idx}')}"
        events.append(Event(
            event_id=ev_id,
            chrom=chrom,
            host_pos=pos,
            del_len=dlen,
            ins_len=ins_len,
            hpv_ref=hpv_name,
            hpv_path=hpv_path_desc,
            hpv_strand=strand,
            copies=copies,
            source=src
        ))

        bedpe_rows.append({
            "chrom1": chrom, "start1": pos, "end1": pos+1,
            "chrom2": hpv_name, "start2": 1, "end2": 2,
            "name": ev_id, "score": ".",
            "strand1": "+", "strand2": strand
        })

    # write mutated FASTA
    if output_only_modified:
        to_write = [(c, "".join(mutable[c])) for c in mutable.keys() if c in modified_chroms]
        if not to_write:
            raise RuntimeError("No chromosomes had insertions; nothing to write.")
    else:
        to_write = [(c, "".join(seq)) for c,seq in mutable.items()]

    mutated = OrderedDict((f"{c}|HPVsiteLists", seq) for c, seq in to_write)
    write_fasta(out_fa_path, mutated)

    ev_df = pd.DataFrame([asdict(e) for e in events])
    ev_path = out_prefix + ".events.tsv"
    ev_df.to_csv(ev_path, sep="\t", index=False)

    bedpe_df = pd.DataFrame(bedpe_rows)
    bedpe_path = out_prefix + ".bedpe.tsv"
    bedpe_df.to_csv(bedpe_path, sep="\t", index=False)

    print(f"FASTA written: {out_fa_path}")
    print(f"Events table : {ev_path}")
    print(f"BEDPE table  : {bedpe_path}")
    print(f"Unique inserts (md5): {len(seen_insert_md5)} / {len(events)}")
    return ev_path, bedpe_path, out_fa_path


# ===================== Example usage =====================
if __name__ == "__main__":
    LINE = [
        "chr11:131789342","chr1:16783209","chr18:22565190","chr2:151136063",
        "chr3:81992132","chr3:82226566","chr7:83120494","chr9:26645008",
        "chr9:97813335","chr9:97814224","chr9:97902954","chr9:97908821",
        "chr9:97914603","chr9:97944420","chr9:97945714","chr9:97952091"
    ]
    SINE = [
        "chr12:66057682","chr3:189879690","chr6:36910796","chr6:37167442",
        "chr6:37176440","chr6:79484344","chr9:116921587"
    ]
    Dup = ["chr1:16783209", "chr1:234783386"]
    lowgc = ["chr4:103351588","chr9:99956760"]

    extra_sites = ["chr1:11651759", "chr1:16783209", "chr1:18033714", "chr1:22412818",
    "chr1:125976080", "chr1:13226313", "chr1:13370676", "chr1:13372446",
    "chr1:13373838", "chr1:13382314", "chr1:13382507", "chr1:13390569",
    "chr2:32916236", "chr2:32916253", "chr2:100978863", "chr2:151136063",
    "chr2:239803038", "chr3:51502997", "chr3:51503437", "chr3:51512714",
    "chr3:51527285", "chr3:51533589", "chr3:51543376", "chr3:81785783",
    "chr3:81992132", "chr3:82226566", "chr3:93470597", "chr3:93470625",
    "chr3:93470710", "chr3:93470734", "chr3:141401612", "chr3:146257083",
    "chr3:163137594", "chr3:189879690", "chr3:189889702", "chr3:189895061",
    "chr4:103351588", "chr4:147855664", "chr5:124471917", "chr6:36910796",
    "chr6:37161250", "chr6:37167442", "chr6:37175775", "chr6:37176440",
    "chr6:79484344", "chr7:83120494", "chr8:46619611", "chr8:33292629",
    "chr8:52372259", "chr8:122107037", "chr9:26645008", "chr9:97813335",
    "chr9:97814224", "chr9:97833766", "chr9:97857784", "chr9:97860564",
    "chr9:97869989", "chr9:97891262", "chr9:97899374", "chr9:97901339",
    "chr9:97902954", "chr9:97908821", "chr9:97914603", "chr9:97940224",
    "chr9:97942608", "chr9:97942781", "chr9:97945714", "chr9:97952091",
    "chr9:99956760", "chr9:116921587", "chr9:128338917", "chr10:5092886",
    "chr10:43440220", "chr11:107488160", "chr11:19907776", "chr11:38609871",
    "chr11:70269302", "chr11:70283619", "chr11:131789342", "chr12:66057678",
    "chr12:66057682", "chr12:98557895", "chr12:124889086", "chr13:73257354",
    "chr13:85257545", "chr15:31335022", "chr15:31335186", "chr17:38267232",
    "chr17:38269138", "chr17:38323632", "chr18:22565190", "chr18:30317896",
    "chr19:4782478", "chr19:56233676", "chr20:21968549", "chr20:23622855",
    "chr21:15273221", "chr21:15315693", "chrX:27024471", "chrX:33461718",
    "chrX:96960123", "chrX:96983771", "chrX:97114880", "chrX:97118472",
    "chrX:97120637", "chrX:97125896"]

    site_lists = {"LINE": LINE, "SINE": SINE, "Dup": Dup, "lowgc": lowgc}
    allow = (
        "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
        "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
        "chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"
    )

    recipes = [

    # =========================
    # Group AB (HPV-only)
    # =========================

    # # A: HPV:(1706-6292)(+), del_len=0
    # {
    #     "name": "A",
    #     "blocks": [
    #         {"type": "hpv", "paths": [(1706, 6292)], "strand": "+", "copies": 1}
    #     ],
    #     "del_len": 0,

    #     # 可选：让 A 的 HPV shift 更强一点
    #     "hpv_shift_max": 2500,
    #     "min_hpv_shift": 200,
    #     "host_shift_max": 0,
    # },

    # # B: HPV:(5305-7906, 1-3372)(-), del_len=134237
    # {
    #     "name": "B",
    #     "blocks": [
    #         {"type": "hpv", "paths": [(5305, 7906), (1, 3372)], "strand": "-", "copies": 1}
    #     ],
    #     "del_len": 134_237,

    #     "hpv_shift_max": 2500,
    #     "min_hpv_shift": 200,
    #     "host_shift_max": 0,
    # },


    # =========================
    # Group CD (HPV-only, longer / multi-copy)
    # =========================

    # C: HPV:(1471-7906, 1-2727)(+), del_len=33063
    {
        "name": "C",
        "blocks": [
            {"type": "hpv", "paths": [(1471, 7906), (1, 2727)], "strand": "+", "copies": 1}
        ],
        "del_len": 33_063,

        "hpv_shift_max": 3000,
        "min_hpv_shift": 200,
        "host_shift_max": 0,
    },

    # D: HPV:(5983-7906, 1-7906, 1-7906, 1-3406)(+), del_len=342
    {
        "name": "D",
        "blocks": [
            {"type": "hpv", "paths": [(5983, 7906), (1, 7906), (1, 7906), (1, 3406)], "strand": "+", "copies": 1}
        ],
        "del_len": 342,

        # D，shift
        "hpv_shift_max": 6000,
        "min_hpv_shift": 500,
        "host_shift_max": 0,
    },


    # =========================
    # Group EFG (HPV+host composite)
    # =========================

    # {"name":"E",
    #  "blocks":[
    #     {"type":"hpv","paths":[(2599,4572)],"strand":"+","copies":1},
    #     {"type":"host","chrom":"chr17","start":36094082,"end":36099602,"strand":"+"},
    #     {"type":"hpv","paths":[(5979,7906),(1,4572)],"strand":"+","copies":1},
    #     {"type":"host","chrom":"chr17","start":36094082,"end":36099602,"strand":"+"},
    #     {"type":"hpv","paths":[(5979,7906),(1,4572)],"strand":"+","copies":1},
    #  ],
    #  "del_len":0,
    #  "hpv_shift_max": 2000,
    #  "host_shift_max": 60000,
    #  "min_hpv_shift": 200,
    #  "min_host_shift": 2000
    # },

    # {"name":"F",
    #  "blocks":[
    #     {"type":"hpv","paths":[(4100,7906),(1,3572)],"strand":"-","copies":1},
    #     {"type":"literal","seq":"TATTA"},
    #     {"type":"hpv","paths":[(6000,7906),(1,2114)],"strand":"+","copies":1},
    #  ],
    #  "del_len":651,
    #  "hpv_shift_max": 2000,
    #  "min_hpv_shift": 200,
    #  "host_shift_max": 0
    # },

    # {"name":"UD2",
    #  "blocks":[
    #     {"type":"hpv","paths":[(1,3106)],"strand":"+","copies":1},
    #     {"type":"host","chrom":"chrX","start":97114879,"end":97125896,"strand":"+"},
    #     {"type":"hpv","paths":[(7710,7892)],"strand":"+","copies":1},
    #     {"type":"hpv","paths":[(1,3106)],"strand":"+","copies":1},
    #     {"type":"host","chrom":"chrX","start":97114879,"end":97125896,"strand":"+"},
    #     {"type":"hpv","paths":[(7710,7892)],"strand":"+","copies":1},
    #     {"type":"hpv","paths":[(1,3106)],"strand":"+","copies":1},
    #     {"type":"host","chrom":"chrX","start":97114879,"end":97125896,"strand":"+"},
    #     {"type":"hpv","paths":[(7710,7892)],"strand":"+","copies":1},
    #     {"type":"hpv","paths":[(1,3106)],"strand":"+","copies":1},
    #     {"type":"host","chrom":"chrX","start":97114879,"end":97125896,"strand":"+"},
    #     {"type":"hpv","paths":[(7710,7892)],"strand":"+","copies":1},
    #     {"type":"hpv","paths":[(1,3106)],"strand":"+","copies":1},
    #     {"type":"host","chrom":"chrX","start":97114879,"end":97125896,"strand":"+"},
    #     {"type":"hpv","paths":[(7710,7892)],"strand":"+","copies":1},
    #     {"type":"hpv","paths":[(1,3106)],"strand":"+","copies":1},
    #     {"type":"host","chrom":"chrX","start":97114879,"end":97125896,"strand":"+"},
    #     {"type":"hpv","paths":[(7710,7892)],"strand":"+","copies":1},
    #  ],
    #  "del_len":0,
    #  "hpv_shift_max": 3000,
    #  "host_shift_max": 80000,
    #  "min_hpv_shift": 300,
    #  "min_host_shift": 5000
    # },
]

    host_fa_path = "hg38.p14.fa"
    hpv_fa_path  = "HPV16pave.fa"
    out_fa_path  = "LevelCD.withHPV.fa.gz"

    _ = insert_events_from_site_lists(
        host_fa_path=host_fa_path,
        hpv_fa_path=hpv_fa_path,
        out_fa_path=out_fa_path,
        site_lists = {"EXTRA": extra_sites},
        extra_sites = None,
        allow_chrs=allow,
        recipes=recipes,
        out_prefix="100_0521planCD_shifted",
        seed=20251110,
        output_only_modified=True,
        target_total=100,
        min_gap=33000,
        ensure_unique_insert=True,
        unique_retry=30,
        avoid_bed="hg38.avoid.bed"
    )
