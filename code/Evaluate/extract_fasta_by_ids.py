#!/usr/bin/env python3
import argparse


def read_ids(path, column, skip_header):
    ids = set()
    with open(path) as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if skip_header and line_no == 1:
                continue
            if len(fields) < column:
                continue
            ids.add(fields[column - 1])
    return ids


def fasta_records(path):
    name = None
    header = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(chunks)
                header = line[1:]
                name = header.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, header, "".join(chunks)


def wrap(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    parser = argparse.ArgumentParser(
        description="Extract FASTA records whose IDs appear in a text file."
    )
    parser.add_argument("--ids", required=True, help="Text file containing contig IDs")
    parser.add_argument("--fasta", required=True, help="Input FASTA")
    parser.add_argument("--out", required=True, help="Output FASTA")
    parser.add_argument(
        "--column",
        type=int,
        default=1,
        help="1-based column in --ids containing contig IDs. Default: 1",
    )
    parser.add_argument(
        "--skip-header",
        action="store_true",
        help="Skip the first non-comment line in --ids",
    )
    args = parser.parse_args()

    keep = read_ids(args.ids, args.column, args.skip_header)
    found = 0
    with open(args.out, "w") as out:
        for name, header, seq in fasta_records(args.fasta):
            if name in keep:
                out.write(f">{header}\n{wrap(seq)}\n")
                found += 1

    missing = len(keep) - found
    print(f"Requested IDs: {len(keep)}")
    print(f"Extracted FASTA records: {found}")
    print(f"Missing IDs: {missing}")


if __name__ == "__main__":
    main()
