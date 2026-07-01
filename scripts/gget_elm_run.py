#!/usr/bin/env python
"""Run gget elm (ortholog + regex) on each region sequence.

Reads data/elm_gget/region.fasta, runs gget.elm() per sequence in ortholog
mode, and writes combined results:
  data/elm_gget/elm_ortho.tsv  - ortholog-transferred motif hits
  data/elm_gget/elm_regex.tsv  - direct regex motif matches
Both carry a `query` column = the FASTA sequence name.
"""
import os
import sys
import warnings

import pandas as pd
import gget

warnings.filterwarnings("ignore")

DATA = "/mnt/c/Users/chris/R_projects/project-hail-mary/data/elm_gget"
FASTA = os.path.join(DATA, "region.fasta")


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def main():
    records = list(read_fasta(FASTA))
    print(f"Loaded {len(records)} sequences", flush=True)

    ortho_all, regex_all = [], []
    for i, (name, seq) in enumerate(records, 1):
        try:
            ortho_df, regex_df = gget.elm(seq, uniprot=False, verbose=False)
            if ortho_df is not None and len(ortho_df):
                ortho_df = ortho_df.copy()
                ortho_df.insert(0, "query", name)
                ortho_all.append(ortho_df)
            if regex_df is not None and len(regex_df):
                regex_df = regex_df.copy()
                regex_df.insert(0, "query", name)
                regex_all.append(regex_df)
            n_o = 0 if ortho_df is None else len(ortho_df)
            n_r = 0 if regex_df is None else len(regex_df)
            print(f"[{i}/{len(records)}] {name}: {n_o} ortho, {n_r} regex", flush=True)
        except Exception as e:  # keep going on per-seq failures
            print(f"[{i}/{len(records)}] {name}: ERROR {e}", flush=True)

    if ortho_all:
        pd.concat(ortho_all, ignore_index=True).to_csv(
            os.path.join(DATA, "elm_ortho.tsv"), sep="\t", index=False)
    if regex_all:
        pd.concat(regex_all, ignore_index=True).to_csv(
            os.path.join(DATA, "elm_regex.tsv"), sep="\t", index=False)

    print(f"DONE ortho_rows={sum(len(d) for d in ortho_all)} "
          f"regex_rows={sum(len(d) for d in regex_all)}", flush=True)


if __name__ == "__main__":
    main()
