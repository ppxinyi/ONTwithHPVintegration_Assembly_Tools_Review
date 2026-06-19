#!/usr/bin/env python3
import argparse
import csv
import re
import sys


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def revcomp(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def fetch_linear(seqs, ref_name, start, end, strand="+"):
    if ref_name not in seqs:
        raise KeyError(f"Reference '{ref_name}' not found in FASTA")
    seq = seqs[ref_name]
    start = max(1, start)
    end = min(len(seq), end)
    if end < start:
        out = ""
    else:
        out = seq[start - 1:end]
    return revcomp(out) if strand == "-" else out


def fetch_circular(seq, start, end):
    """Fetch 1-based inclusive interval from a circular sequence in + orientation."""
    n = len(seq)
    if n == 0:
        return ""
    start = ((start - 1) % n) + 1
    end = ((end - 1) % n) + 1
    if start <= end:
        return seq[start - 1:end]
    return seq[start - 1:] + seq[:end]


def split_top_level_path(path):
    """Split hpv_path on + signs that join path tokens, not strand '(+)' marks."""
    parts = []
    buf = []
    depth = 0
    for ch in path:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "+" and depth == 0:
            parts.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    if buf:
        parts.append("".join(buf))
    return [p.strip() for p in parts if p.strip()]


def get_hpv_seq(hpv_seqs, hpv_ref):
    if hpv_ref in hpv_seqs:
        return hpv_seqs[hpv_ref]
    if len(hpv_seqs) == 1:
        return next(iter(hpv_seqs.values()))
    candidates = [name for name in hpv_seqs if hpv_ref in name or name in hpv_ref]
    if len(candidates) == 1:
        return hpv_seqs[candidates[0]]
    raise KeyError(
        f"Could not choose HPV sequence for hpv_ref='{hpv_ref}'. "
        f"Available FASTA names: {', '.join(hpv_seqs)[:500]}"
    )


def build_insert_sequence(path, hpv_seq, host_seqs, literal_base="N"):
    out = []
    for token in split_top_level_path(path):
        if token.startswith("LIT:"):
            length = int(token.split(":", 1)[1])
            out.append(literal_base.upper() * length)
            continue

        m = re.fullmatch(r"HPV:(\d+)-(\d+)(?:,(\d+)-(\d+))*\(([+-])\)", token)
        if m:
            strand = m.group(5)
            pieces = []
            for s, e in re.findall(r"(\d+)-(\d+)", token.split("(", 1)[0]):
                pieces.append(fetch_circular(hpv_seq, int(s), int(e)))
            seq = "".join(pieces)
            out.append(revcomp(seq) if strand == "-" else seq)
            continue

        m = re.fullmatch(r"([^:]+):(\d+)-(\d+)\(([+-])\)", token)
        if m:
            chrom, start, end, strand = m.groups()
            out.append(fetch_linear(host_seqs, chrom, int(start), int(end), strand))
            continue

        raise ValueError(f"Could not parse hpv_path token: {token}")
    return "".join(out).upper()


def wrap(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    parser = argparse.ArgumentParser(
        description="Create event-specific truth FASTA from HPV integration event TSV."
    )
    parser.add_argument("--events", required=True, help="Input event TSV with hpv_path column")
    parser.add_argument("--hg38", required=True, help="Host genome FASTA")
    parser.add_argument("--hpv", required=True, help="HPV FASTA")
    parser.add_argument("--out", required=True, help="Output event truth FASTA")
    parser.add_argument("--flank", type=int, default=5000, help="Host flank length on each side")
    parser.add_argument(
        "--literal-base",
        default="N",
        help="Base to use for LIT:length path tokens. Default: N",
    )
    args = parser.parse_args()

    host_seqs = read_fasta(args.hg38)
    hpv_seqs = read_fasta(args.hpv)

    n = 0
    with open(args.events, newline="") as handle, open(args.out, "w") as out:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"event_id", "chrom", "host_pos", "del_len", "ins_len", "hpv_ref", "hpv_path"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Missing required columns: {', '.join(sorted(missing))}")

        for row in reader:
            event_id = row["event_id"]
            chrom = row["chrom"]
            host_pos = int(row["host_pos"])
            del_len = int(row["del_len"])
            expected_ins_len = int(row["ins_len"])
            hpv_seq = get_hpv_seq(hpv_seqs, row["hpv_ref"])

            left = fetch_linear(host_seqs, chrom, host_pos - args.flank, host_pos - 1)
            insert = build_insert_sequence(
                row["hpv_path"], hpv_seq, host_seqs, literal_base=args.literal_base
            )
            right_start = host_pos + del_len
            right = fetch_linear(host_seqs, chrom, right_start, right_start + args.flank - 1)
            truth = left + insert + right

            if expected_ins_len and len(insert) != expected_ins_len:
                print(
                    f"WARNING {event_id}: built insert length {len(insert)} "
                    f"!= ins_len {expected_ins_len}",
                    file=sys.stderr,
                )

            header = (
                f">{event_id}|{chrom}:{host_pos}|del={del_len}|ins={len(insert)}|"
                f"flank={args.flank}|left_junction={len(left)}|"
                f"right_junction={len(left) + len(insert)}"
            )
            out.write(header + "\n")
            out.write(wrap(truth) + "\n")
            n += 1

    print(f"Wrote {n} event truth sequences to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
