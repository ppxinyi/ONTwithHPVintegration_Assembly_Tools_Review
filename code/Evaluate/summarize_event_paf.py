#!/usr/bin/env python3
import argparse
import re


def parse_truth_fasta(path):
    truth = {}
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            name = header.split()[0]
            event_id = name.split("|", 1)[0]
            fields = dict(re.findall(r"([A-Za-z_]+)=([^|]+)", header))
            truth[event_id] = {
                "target": name,
                "left_junction": int(fields.get("left_junction", 0)),
                "right_junction": int(fields.get("right_junction", 0)),
            }
    return truth


def paf_tags(fields):
    tags = {}
    for item in fields[12:]:
        parts = item.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def aln_identity(fields):
    tags = paf_tags(fields)
    matches = int(fields[9])
    block = int(fields[10])
    return matches / block if block else 0.0


def covers(pos_start, pos_end, junction, pad):
    return pos_start <= junction + pad and pos_end >= junction - pad


def classify(truth_cov, left, right, identity, min_identity, min_full_cov):
    if left and right and identity >= min_identity and truth_cov >= min_full_cov:
        return "PASS"
    if left or right or truth_cov > 0:
        return "PARTIAL"
    return "MISS"


def main():
    parser = argparse.ArgumentParser(
        description="Summarize assembly-vs-event-truth PAF into event-level calls."
    )
    parser.add_argument("--truth-fa", required=True, help="truth_events.fa")
    parser.add_argument("--paf", required=True, help="minimap2 PAF")
    parser.add_argument("--out", required=True, help="Output TSV")
    parser.add_argument(
        "--junction-pad",
        type=int,
        default=50,
        help="Allow this many bp around junction when checking coverage. Default: 50",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.95,
        help="Minimum identity for PASS. Default: 0.95",
    )
    parser.add_argument(
        "--min-full-cov",
        type=float,
        default=0.80,
        help="Minimum truth coverage for PASS. Default: 0.80",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum MAPQ to consider. Default: 0",
    )
    args = parser.parse_args()

    truth = parse_truth_fasta(args.truth_fa)
    best = {}

    with open(args.paf) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query = fields[0]
            target = fields[5]
            event_id = target.split("|", 1)[0]
            if event_id not in truth:
                continue

            mapq = int(fields[11])
            if mapq < args.min_mapq:
                continue

            tlen = int(fields[6])
            tstart = int(fields[7])
            tend = int(fields[8])
            identity = aln_identity(fields)
            truth_cov = (tend - tstart) / tlen if tlen else 0.0
            left = covers(tstart, tend, truth[event_id]["left_junction"], args.junction_pad)
            right = covers(tstart, tend, truth[event_id]["right_junction"], args.junction_pad)

            # Prefer alignments that cover junctions, then coverage, then identity, then MAPQ.
            score = (
                int(left) + int(right),
                truth_cov,
                identity,
                mapq,
                int(fields[10]),
            )
            if event_id not in best or score > best[event_id]["score"]:
                best[event_id] = {
                    "score": score,
                    "best_contig": query,
                    "identity": identity,
                    "truth_cov": truth_cov,
                    "left_junc": "yes" if left else "no",
                    "right_junc": "yes" if right else "no",
                    "mapq": mapq,
                    "target_start": tstart,
                    "target_end": tend,
                }

    with open(args.out, "w") as out:
        out.write(
            "event_id\tbest_contig\tidentity\ttruth_cov\tleft_junc\t"
            "right_junc\tstatus\tmapq\ttarget_start\ttarget_end\n"
        )
        for event_id in truth:
            row = best.get(event_id)
            if row is None:
                out.write(f"{event_id}\t.\t.\t0\tno\tno\tMISS\t.\t.\t.\n")
                continue
            status = classify(
                row["truth_cov"],
                row["left_junc"] == "yes",
                row["right_junc"] == "yes",
                row["identity"],
                args.min_identity,
                args.min_full_cov,
            )
            out.write(
                f"{event_id}\t{row['best_contig']}\t{row['identity']:.4f}\t"
                f"{row['truth_cov']:.4f}\t{row['left_junc']}\t"
                f"{row['right_junc']}\t{status}\t{row['mapq']}\t"
                f"{row['target_start']}\t{row['target_end']}\n"
            )


if __name__ == "__main__":
    main()
