"""
This script will download the KEGG module database and save it as a csv file.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
from Bio.KEGG import REST

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Post-processing")


# from REST documentation, can get a max of 10 queries
# in kegg_get
MAX_N_QUERIES = 10


def parse_modules(handle):
    resp = handle.read()
    modules = []
    for block in resp.split("///"):
        if not block.strip():
            continue
        entry = name = definition = module_type = module_category = (
            module_subcategory
        ) = None
        for line in block.split("\n"):
            token = line[:12].strip()
            value = line[12:]
            if not token:
                continue
            if token == "ENTRY":
                entry = value.split(None, 1)[0].strip()
            elif token == "NAME":
                name = value.strip()
            elif token == "DEFINITION":
                definition = value.strip()
            elif token == "CLASS":
                module_type, module_category, module_subcategory = (
                    el.strip() for el in value.split(";")
                )
            if all(
                [
                    entry,
                    name,
                    definition,
                    module_type,
                    module_category,
                    module_subcategory,
                ]
            ):
                break
        modules.append(
            [entry, name, module_type, module_category, module_subcategory, definition]
        )
    return modules


if __name__ == "__main__":
    args = argparse.ArgumentParser("annotations.py")
    args.add_argument("--db_dir", default="db", help="path to the database directory")

    args = args.parse_args()
    db_dir = Path(args.db_dir)
    if not db_dir.exists():
        db_dir.mkdir()

    res = REST.kegg_list("module")
    modules = pd.read_csv(res, sep="\t", names=["module", "description"])

    modules_data = []
    logger.info("Downloading module data")
    for i in range(0, len(modules), MAX_N_QUERIES):
        resp = REST.kegg_get(list(modules["module"].iloc[i : i + MAX_N_QUERIES]))
        modules_data.extend(parse_modules(resp))

    df = pd.DataFrame(
        modules_data,
        columns=[
            "entry",
            "name",
            "module_type",
            "module_category",
            "module_subcategory",
            "definition",
        ],
    )
    df.to_csv(Path(db_dir, "kegg_modules.csv"), index=False)
