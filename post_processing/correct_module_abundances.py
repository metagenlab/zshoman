"""
The current eggnog database uses an older version of the KEGG database, which includes
modules that do not exist anymore. We therefore rewrite the modules table, using the
modules associated to KOs to write the table
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).parent.parent))

from utils.utils import parse_arguments


def load_ko_to_module(db_dir):
    return pd.read_csv(Path(db_dir, "ko_to_modules.csv"), index_col=0)


if __name__ == "__main__":
    args = parse_arguments(
        samples_file=False,
        postprocessed_dir=True,
        db_dir=True,
    )

    ko_to_modules = load_ko_to_module(args.db_dir)

    kos = pd.read_csv(Path(args.postprocessed_dir, "KEGG_ko.csv"), index_col=0)
    kos.index = kos.index.str.replace("ko:", "")

    ko_to_modules.index.rename("KEGG_ko", inplace=True)
    # Lists are loaded as strings
    ko_to_modules["module"] = ko_to_modules["module"].apply(eval)

    modules = kos.join(ko_to_modules).explode("module").groupby("module").sum()
    modules.to_csv(Path(args.postprocessed_dir, "KEGG_modules_corrected.csv"))
