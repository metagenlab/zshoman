import argparse

import pandas as pd


def get_ko_id(string):
    return int(string[1:])


def prepare_ko_table(ko_files):
    data = []
    for ko_file in ko_files:
        curr_gene = None
        with open(ko_file, "r") as ko_fh:
            for ko_line in ko_fh:
                tokens = ko_line.split()
                # ignore all but the best hits
                if tokens[0] != "*":
                    continue
                gene, ko_str, thrs_str, score_str, evalue_str, * \
                    descr = tokens[1:]
                if gene == curr_gene:
                    # skip the entries that were classified as significant, but
                    # with a higher e-value
                    continue
                else:
                    curr_gene = gene
                ko = get_ko_id(ko_str)
                thrs = float(thrs_str)
                score = float(score_str)
                evalue = float(evalue_str)
                entry = [gene, ko, thrs, score, evalue]
                data.append(entry)
    df = pd.DataFrame(data, columns=["gene", "ko", "threshold", "score", "evalue"])
    return df


def main():
    parser = argparse.ArgumentParser(description='Merge COG hit files')
    parser.add_argument('ko_hit_files', nargs='*', help='hit files')
    args = parser.parse_args()

    df = prepare_ko_table(args.ko_hit_files)
    df.to_csv("kos.csv", index=False)


if __name__ == '__main__':
    main()
