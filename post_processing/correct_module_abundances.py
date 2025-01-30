"""
The current eggnog database uses an older version of the KEGG database, which includes
modules that do not exist anymore. We therefore rewrite the modules table, using the
modules associated to KOs to write the table
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
from Bio.KEGG import REST

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Post-processing")


def load_ko_to_module(db_dir):
    return pd.read_csv(Path(db_dir, "ko_to_modules.csv"), index_col=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("correct_module_abundances.py")
    parser.add_argument("--db_dir", default="db", help="path to the database directory")

    parser.add_argument(
        "-o",
        "--output_dir",
        default="output/post_processed",
        help="path to the post_processed annotation tables",
    )

    args = parser.parse_args()
    db_dir = Path(args.db_dir)
    output_dir = Path(args.output_dir)

    ko_to_modules = load_ko_to_module(db_dir)

    kos = pd.read_csv(Path(output_dir, "KEGG_ko.csv"), index_col=0)
    kos.index = kos.index.str.replace("ko:", "")

    ko_to_modules.index.rename("KEGG_ko", inplace=True)
    # Lists are loaded as strings
    ko_to_modules["module"] = ko_to_modules["module"].apply(eval)

    modules = kos.join(ko_to_modules).explode("module").groupby("module").sum()
    modules.to_csv(Path(output_dir, "KEGG_modules_corrected.csv"))
