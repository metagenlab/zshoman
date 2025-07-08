"""
This script will calculate the Ko module completeness and add it to the corresponding table.
It needs access to the KEGG modules database which can be downloaded with the
download_kegg_db.py script. It also needs the corrected table of module abundances
obtained with the correct_module_abundances.py script.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
from Bio.KEGG import REST

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Post-processing")


def load_module_data(db_dir):
    return pd.read_csv(Path(db_dir, "kegg_modules.csv"), index_col=0)


def parse_definitions(module_data):
    """
    We prepare the definitions to be evaluable as boolean expressions.
    The definitions include optional parts indicated by the - sign. + and
    space represent AND, while coma represent OR.
    Some definitions have -- which is undefined, so we simply skip them.
    """
    module_data["definition"] = (
        module_data["definition"]
        .str.replace(r"--", "", regex=True)
        .str.replace(r"-K\d{5}", "", regex=True)
        .str.replace(r"-\(.*?\)", "", regex=True)
        .str.strip()
        .str.replace(" +", " and ", regex=True)
        .str.replace("+", " and ")
        .str.replace(",", " or ")
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser("annotations.py")
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

    module_data = load_module_data(db_dir)
    parse_definitions(module_data)

    modules = pd.read_csv(Path(output_dir, "KEGG_modules_corrected.csv"), index_col=0)
    modules = modules.join(module_data["definition"])
    kos = pd.read_csv(Path(output_dir, "KEGG_ko.csv"), index_col=0)
    kos.index = kos.index.str.replace("ko:", "")
    for i, sample in enumerate(kos, 1):
        if i % 100 == 0:
            logger.info(f"Done {i}/{len(kos.columns)}")
        kos_present = "|".join(kos.index[kos[sample] > 0])
        if not kos_present:
            # not a single KO present, which means all modules are already
            # set to 0 as well
            modules[sample] = 0
            continue
        complete_modules = (
            modules["definition"]
            .str.replace(kos_present, "True", regex=True)
            .replace(r"K\d{5}", "False", regex=True)
            .replace(r"M\d{5}", "False", regex=True)
            .map(eval, na_action="ignore")
        )
        modules[sample] = modules[sample].where(complete_modules, 0)
    modules.drop(columns="definition", inplace=True)
    modules.to_csv(Path(output_dir, "KEGG_complete_modules.csv"))
