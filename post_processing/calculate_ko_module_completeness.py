"""
This script will calculate the Ko module completeness and add it to the corresponding table.
It needs access to the KEGG modules database which can be downloaded with the
download_kegg_db.py script. It also needs the corrected table of module abundances
obtained with the correct_module_abundances.py script.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).parent.parent))

import subprocess
import tempfile

from utils.utils import parse_arguments

if __name__ == "__main__":
    args = parse_arguments(
        samples_file=False,
        postprocessed_dir=True,
    )

    kos = pd.read_csv(Path(args.postprocessed_dir, "KEGG_ko.csv"), index_col=0)
    kos.index = kos.index.str.replace("ko:", "")

    # Create a table with samples followed by tab-delimited list of KOs present
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir=Path(tmpdir)
        input_file = tmpdir / "kos.csv"
        with input_file.open("w") as fout:
            for sample, row in kos.T.iterrows():
                fout.write(f"{sample}\t{'\t'.join(row[row != 0].index)}\n")

        subprocess.call(["give_completeness", "-i", input_file, "-o", tmpdir, "--add-per-contig"])
        completeness = pd.read_csv(tmpdir / "summary.kegg_contigs.tsv", sep="\t")

        # The output table has one line for each module for each sample. We will reshape this
        # to just keep one line per module with samples as columns
        new_columns = ["contig"]
        new_indices = ["module_accession", "pathway_name", "pathway_class"]
        completeness = completeness[new_columns + new_indices + ["completeness"]].pivot(
            columns=new_columns, index=new_indices).fillna(0)
        # The result here has a multiindex columns ("compelteness", "sample"), so we change this
        # to keep only the sample names
        completeness.columns = [el[1] for el in completeness.columns]
        completeness.to_csv(args.postprocessed_dir / "KEGG_modules_completeness.csv", index=True)