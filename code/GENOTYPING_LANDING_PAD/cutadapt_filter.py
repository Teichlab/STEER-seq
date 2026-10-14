# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 13:43:39 2026

@author: ig8
"""

import argparse

from Bio import SeqIO
from cutadapt.adapters import FrontAdapter, BackAdapter


def revcomp(seq):
    return seq.translate(
        str.maketrans("ACGTacgt", "TGCAtgca")
    )[::-1]


def load_primers(filename):
    records = {
        record.id.lower(): str(record.seq).upper()
        for record in SeqIO.parse(filename, "fasta")
    }

    if "fw" not in records:
        raise ValueError("Primer FASTA must contain >fw")

    if "rev" not in records:
        raise ValueError("Primer FASTA must contain >rev")

    return records["fw"], records["rev"]


def make_adaptors(primer_fasta):
    fw_primer, rev_primer = load_primers(primer_fasta)

    # Forward cassette:
    # A1 ... A2
    A1 = fw_primer
    A2 = revcomp(rev_primer)

    # Reverse cassette:
    # A2_RC ... A1_RC
    A1_RC = revcomp(A1)
    A2_RC = rev_primer

    f1 = FrontAdapter(
        A1,
        max_errors=0.15,
        indels=True,
    )

    f2 = BackAdapter(
        A2,
        max_errors=0.15,
        indels=True,
    )

    r1 = FrontAdapter(
        A2_RC,
        max_errors=0.15,
        indels=True,
    )

    r2 = BackAdapter(
        A1_RC,
        max_errors=0.15,
        indels=True,
    )

    return [f1, f2, r1, r2]


def find_all_matches(seq, adapter):
    """
    Yield successive adapter matches by repeatedly searching
    the remaining suffix of the sequence.
    """
    offset = 0

    while offset < len(seq):
        m = adapter.match_to(seq[offset:])

        if m is None:
            break

        # Convert coordinates back to the original sequence
        m.rstart += offset
        m.rstop += offset

        yield m

        # Continue searching after this match
        offset = m.rstop


def find_cassettes(seq, adaptors):
    f1, f2, r1, r2 = adaptors

    candidates = []

    # Forward cassettes: A1 ... A2
    for m_f1 in find_all_matches(seq, f1):
        tail = seq[m_f1.rstart:]

        m_f2 = f2.match_to(tail)

        if m_f2 is None:
            continue

        candidates.append(
            (
                "F",
                m_f1.rstart,
                m_f1.rstart + m_f2.rstop,
            )
        )

    # Reverse cassettes: A2RC ... A1RC
    for m_r1 in find_all_matches(seq, r1):
        tail = seq[m_r1.rstart:]

        m_r2 = r2.match_to(tail)

        if m_r2 is None:
            continue

        candidates.append(
            (
                "R",
                m_r1.rstart,
                m_r1.rstart + m_r2.rstop,
            )
        )

    # Return cassettes in sequence order
    candidates.sort(key=lambda x: x[1])

    yield from candidates


VALID_PAIRS = {
    "A1_A2",
    "A2RC_A1RC",
}


INVALID_PAIRS = {
    "A1_A1",
    "A1_A1RC",
    "A1_A2RC",

    "A2_A1",
    "A2_A2",
    "A2_A1RC",
    "A2_A2RC",

    "A1RC_A1",
    "A1RC_A2",
    "A1RC_A1RC",
    "A1RC_A2RC",

    "A2RC_A1",
    "A2RC_A2",
    "A2RC_A2RC",
}


def classify_subreads_extra(seq, adaptors):
    """
    Classify adjacent adapter pairs as VALID, BOUNDARY,
    MALFORMED, or UNKNOWN.

    This is diagnostic only. It does not determine which
    subreads are written to the output.
    """

    valid = list(find_cassettes(seq, adaptors))
    valid_intervals = [(s, e) for _, s, e in valid]

    def containing_interval(pos):
        for i, (s, e) in enumerate(valid_intervals):
            if s <= pos <= e:
                return i

        return None

    f1, f2, r1, r2 = adaptors

    events = []

    for m in find_all_matches(seq, f1):
        events.append(("A1", m))

    for m in find_all_matches(seq, f2):
        events.append(("A2", m))

    for m in find_all_matches(seq, r1):
        events.append(("A2RC", m))

    for m in find_all_matches(seq, r2):
        events.append(("A1RC", m))

    events.sort(key=lambda x: x[1].rstart)

    for (name1, m1), (name2, m2) in zip(events, events[1:]):

        pair = f"{name1}_{name2}"

        if pair in VALID_PAIRS:

            status = "VALID"
            start = m1.rstart
            stop = m2.rstop

        elif pair in INVALID_PAIRS:

            left = containing_interval(m1.rstart)
            right = containing_interval(m2.rstop)

            if (
                left is not None
                and right is not None
                and left != right
            ):
                status = "BOUNDARY"
                start = m1.rstart
                stop = m2.rstop

            else:
                status = "MALFORMED"
                start = m1.rstart
                stop = m2.rstop

                # Clip malformed region to the start
                # of the next valid cassette.
                for _, s, e in valid:
                    if s > m1.rstart:
                        stop = s
                        break

        else:

            status = "UNKNOWN"
            start = m1.rstart
            stop = m2.rstop

        yield pair, status, start, stop


def write_valid_subreads(input_fastq, output_fastq, adaptors):
    """
    Write every complete valid cassette as a separate FASTQ record.

    Output IDs are:
        original_id_1
        original_id_2
        original_id_3
        ...
    """

    with open(output_fastq, "w") as out:

        for record in SeqIO.parse(input_fastq, "fastq"):

            seq = str(record.seq)

            subread_num = 0

            for orient, start, stop in find_cassettes(
                seq,
                adaptors,
            ):

                subread_num += 1

                subread = record[start:stop]

                subread.id = f"{record.id}_{subread_num}"
                subread.name = subread.id
                subread.description = ""

                SeqIO.write(
                    subread,
                    out,
                    "fastq",
                )


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Extract complete forward and reverse adapter "
            "cassettes from FASTQ reads."
        )
    )

    parser.add_argument(
        "input_fastq",
        help="Input FASTQ file",
    )

    parser.add_argument(
        "output_fastq",
        help="Output FASTQ file containing valid subreads",
    )

    parser.add_argument(
        "--primers",
        required=True,
        help=(
            "FASTA containing >fw and >rev primer sequences"
        ),
    )

    args = parser.parse_args()

    adaptors = make_adaptors(args.primers)

    write_valid_subreads(
        args.input_fastq,
        args.output_fastq,
        adaptors,
    )


if __name__ == "__main__":
    main()
    
