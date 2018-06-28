#!/usr/bin/env python
# coding=utf-8
"""
save only fungi from euks reference before joining with prokaryotes
"""
import click
import sys
from tqdm import tqdm


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.argument("misc_taxonomy")
@click.argument("family_euks_faa")
def main(misc_taxonomy, family_euks_faa):
    g1 = ""
    g2 = ""
    g3 = ""
    # ag : addendum categories
    exclude = ["ag"]
    with open(misc_taxonomy) as fh:
        for line in fh:
            line = line.strip("\r\n")
            if line.startswith("# "):
                g1 = line.partition(" ")[-1]
                continue
            elif line.startswith("## "):
                g2 = line.partition(" ")[-1]
                continue
            elif line.startswith("### "):
                g3 = line.partition(" ")[-1]
                continue
            elif line.startswith("#### "):
                continue
            toks = line.split("\t")
            if g1 == "Eukaryotes":
                if g2 == "Fungi":
                    continue
                else:
                    exclude.append(toks[1])
    with open(family_euks_faa) as fh, tqdm(desc="Records processed") as pbar:
        seen_and_excluded = 0
        exclude_sequence = False
        for line in fh:
            line = line.strip("\r\n")
            if line.startswith(">"):
                pbar.update()
                exclude_sequence = False
                for code in exclude:
                    if line.startswith(">%s" % code):
                        exclude_sequence = True
                        break
                if not exclude_sequence:
                    print(line)
                else:
                    seen_and_excluded += 1
            else:
                if not exclude_sequence:
                    print(line)
        print(
            "Excluded",
            seen_and_excluded,
            "genes from",
            len(exclude),
            "species",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
