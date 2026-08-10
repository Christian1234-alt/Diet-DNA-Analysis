#!/usr/bin/env python3

import argparse
import csv


# IUPAC nucleotide definitions
IUPAC = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"}
}


COMPLEMENT = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
    "N": "N"
}


def reverse_complement(seq):
    return "".join(COMPLEMENT[b] for b in reversed(seq.upper()))


def read_fasta(filename):

    header = None
    sequence = []

    with open(filename) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header:
                    yield header, "".join(sequence).upper()

                header = line[1:]
                sequence = []

            else:
                sequence.append(line)

        if header:
            yield header, "".join(sequence).upper()



def bases_match(primer_base, sequence_base):

    return not IUPAC[primer_base].isdisjoint(
        IUPAC.get(sequence_base, set())
    )



def count_mismatches(primer, sequence):

    mismatches = 0

    for p, s in zip(primer, sequence):

        if not bases_match(p, s):
            mismatches += 1

    return mismatches



def find_primer_hits(sequence, primer, max_mismatches):

    hits = []

    length = len(primer)

    for i in range(len(sequence) - length + 1):

        window = sequence[i:i+length]

        mismatches = count_mismatches(primer, window)

        if mismatches <= max_mismatches:
            hits.append(
                {
                    "position": i,
                    "mismatches": mismatches
                }
            )

    return hits



def main():

    parser = argparse.ArgumentParser(
        description="In silico PCR with IUPAC ambiguity and mismatch support"
    )

    parser.add_argument("--fasta", required=True)
    parser.add_argument("--forward", required=True)
    parser.add_argument("--reverse", required=True)
    parser.add_argument("--mismatches", type=int, default=0)
    parser.add_argument("--output", default="primer_results.csv")

    args = parser.parse_args()


    forward_primer = args.forward.upper()
    reverse_primer = reverse_complement(
        args.reverse.upper()
    )


    total = 0
    amplified = 0


    with open(args.output, "w", newline="") as outfile:

        writer = csv.writer(outfile)

        writer.writerow([
            "Accession",
            "Description",
            "Forward_position",
            "Forward_mismatches",
            "Reverse_position",
            "Reverse_mismatches",
            "Amplicon_length",
            "Amplifies"
        ])


        for header, sequence in read_fasta(args.fasta):

            total += 1

            accession = header.split()[0]


            forward_hits = find_primer_hits(
                sequence,
                forward_primer,
                args.mismatches
            )


            reverse_hits = find_primer_hits(
                sequence,
                reverse_primer,
                args.mismatches
            )


            valid_amplicons = []


            for f in forward_hits:

                for r in reverse_hits:

                    if r["position"] > f["position"]:

                        length = (
                            r["position"]
                            + len(reverse_primer)
                            - f["position"]
                        )

                        valid_amplicons.append(
                            {
                                "forward": f,
                                "reverse": r,
                                "length": length
                            }
                        )


            if valid_amplicons:

                # Choose shortest PCR product
                best = min(
                    valid_amplicons,
                    key=lambda x: x["length"]
                )


                writer.writerow([
                    accession,
                    header,
                    best["forward"]["position"] + 1,
                    best["forward"]["mismatches"],
                    best["reverse"]["position"] + 1,
                    best["reverse"]["mismatches"],
                    best["length"],
                    "Yes"
                ])

                amplified += 1


            else:

                writer.writerow([
                    accession,
                    header,
                    "",
                    "",
                    "",
                    "",
                    "",
                    "No"
                ])


    print()
    print("Primer analysis complete")
    print("------------------------")
    print(f"Sequences analyzed : {total}")
    print(f"Amplified          : {amplified}")

    if total > 0:
        print(f"Coverage           : {100*amplified/total:.2f}%")

    print(f"Results file       : {args.output}")



if __name__ == "__main__":
    main()
