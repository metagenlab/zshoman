import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='Merge COG hit files')
    parser.add_argument('cdd_to_cog', help='CDD to COG mapping')
    parser.add_argument('cog_hit_files', nargs='*', help='hit files')
    args = parser.parse_args()

    hsh_cdd_to_cog = {}
    with open(args.cdd_to_cog, "r") as cdd_to_cog_file:
        for line in cdd_to_cog_file:
            cog, cdd = line.split()
            hsh_cdd_to_cog[int(cdd)] = int(cog)

    data = []
    for chunk in args.cog_hit_files:
        cogs_hits = pd.read_csv(chunk, sep="\t", header=None,
                                names=["seq_hsh", "cdd", "pident", "length", "mismatch", "gapopen", "qstart",
                                       "qend", "sstart", "send", "evalue", "bitscore"])

        # Select only the best hits: using pandas clearly is an overkill here
        cogs_hits = cogs_hits[["seq_hsh", "cdd", "evalue", "pident"]]
        min_hits = cogs_hits.sort_values(
            ["evalue", "pident"],
            ascending=[True, False]).drop_duplicates("seq_hsh")
        for index, row in min_hits.iterrows():
            cog = hsh_cdd_to_cog[int(row["cdd"].split(":")[1])]
            evalue = float(row["evalue"])
            data.append([row["seq_hsh"], cog, evalue])
    df = pd.DataFrame(data, columns=["gene", "cog", "evalue"])
    df.to_csv("cogs.csv", index=False)


if __name__ == '__main__':
    main()
