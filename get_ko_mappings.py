import argparse
from collections import defaultdict

import pandas as pd
from Bio.KEGG import REST

# from REST documentation, can get a max of 10 queries
# in kegg_get
MAX_N_QUERIES = 10


class Pathway(object):
    def __init__(self, entry, descr):
        self.entry = entry
        self.descr = descr

    def simplified_entry(self):
        return int(self.entry[len("map"):])

    def __hash__(self):
        return hash(self.entry)

    def __eq__(self, other):
        return self.entry == other.entry

    def __str__(self):
        return self.entry + ": " + self.descr


class Record(object):
    def __init__(self, entry, name, modules, pathways):
        self.entry = entry
        self.name = name
        self.modules = modules
        self.pathways = pathways

    def simplified_entry(self):
        return int(self.entry[len("K"):])

    def __str__(self):
        acc = self.entry + " " + self.definition
        return acc


class Module(object):
    def __init__(self):
        self.entry = None
        self.descr = None
        self.classes = []
        self.definition = []

    def simplified_entry(self):
        return int(self.entry[len("M"):])

    def get_definition(self):
        if len(self.definition) == 1:
            return self.definition[0]
        # Module definitions may come on multiple lines, join them
        # with an AND
        return " ".join(f"{part}" for part in self.definition)

    def sub_category(self):
        return self.classes[-1]

    def category(self):
        return self.classes[-2]

    def is_signature(self):
        return self.classes[0] == "Signature modules"

    def __str__(self):
        return self.entry + ": " + self.descr

    def __hash__(self):
        return hash(self.entry)

    def __eq__(self, other):
        return self.entry == other.entry


def get_ko_modules_mapping():
    data = REST.kegg_link("module", "ko")
    mapping = defaultdict(list)
    for line in data:
        entry, module = line.split("\t")
        entry = entry.split("ko:")[1].strip()
        module = module.split("md:")[1].strip()
        # it seems there are a few dupplicate
        if module not in mapping[entry]:
            mapping[entry].append(module)
    return mapping


def get_ko_pathways_mapping():
    pathway_mapping = {}
    for line in REST.kegg_list("pathway"):
        entry, descr = line.split("\t")
        entry = entry.strip()
        descr = descr.strip()
        pathway_mapping[entry] = Pathway(entry, descr)

    data = REST.kegg_link("pathway", "ko")
    mapping = defaultdict(list)
    for line in data:
        entry, pathway_id = line.split("\t")
        entry = entry.split("ko:")[1].strip()
        pathway_id = pathway_id.split("path:")[1].strip()
        if pathway_id.startswith("map"):
            mapping[entry].append(pathway_mapping[pathway_id])
    return mapping


def chunk_modules(modules):
    curr = []
    for i, module in enumerate(modules):
        curr.append(module)
        if i % MAX_N_QUERIES == 0:
            yield curr
            curr = []
    if len(curr) > 0:
        yield curr


def parse_module(handle):
    module = Module()
    for line in handle:
        if line[:3] == "///":
            yield module
            module = Module()
            continue
        if line[:12] != "            ":
            keyword = line[:12].strip()
        data = line[12:].strip()
        if keyword == "ENTRY":
            module.entry = data.split()[0]
        if keyword == "NAME":
            module.descr = data
        if keyword == "CLASS":
            tokens = data.split(";")
            module.classes = [token.strip() for token in tokens]
        if keyword == "DEFINITION":
            module.definition.append(data)


def load_KO_references():
    ko_modules = get_ko_modules_mapping()
    ko_pathways = get_ko_pathways_mapping()

    genes = []
    for line in REST.kegg_list("KO"):
        entry, rest = line.split("\t", 1)
        if ";" in rest:
            name = rest.split(";", 1)[1]  # Record for K17548 has two semicola.
        else:
            # record K23479 does not have a semicolon
            name = rest.split("/")[1]
        entry = entry.strip()
        name = name.strip()
        genes.append(
            Record(entry, name, ko_modules[entry], ko_pathways[entry])
            )

    # get only the module present in the genes
    modules_set = {module for gene in genes for module in gene.modules}
    hsh_modules = {}
    print("Downloading module data")
    for i, module_list in enumerate(chunk_modules(modules_set)):
        mod_data = REST.kegg_get(module_list)
        for module in parse_module(mod_data):
            hsh_modules[module.entry] = module
    print("Done")

    pathway_set = set()
    category_set = set()
    for entry, module in hsh_modules.items():
        if len(module.classes) > 0:
            category_set.add(module.sub_category())
            category_set.add(module.category())
    for gene in genes:
        for pathway in gene.pathways:
            pathway_set.add(pathway)

    # really ugly, but spares some web requests
    module_classes = []
    hsh_class_to_id = {}
    for cat_id, cat in enumerate(category_set):
        module_classes.append((cat_id, cat))
        hsh_class_to_id[cat] = cat_id
    for entry, module in hsh_modules.items():
        cat_id = hsh_class_to_id[module.category()]
        subcat_id = hsh_class_to_id[module.sub_category()]
        module.cat_id = cat_id
        module.subcat_id = subcat_id

    module_classes_df = pd.DataFrame(module_classes, columns=["class_id", "description"])

    ko_modules = [
        (m.simplified_entry(),
         m.descr, m.get_definition(),
         m.is_signature(),
         m.cat_id, m.subcat_id)
        for m in hsh_modules.values()]
    modules_df = pd.DataFrame(
        ko_modules,
        columns=["module_id", "description", "definition",
                 "is_signature_module", "class_id", "subclass_id"])

    pathways_df = pd.DataFrame(
        [(p.simplified_entry(), p.descr) for p in pathway_set],
        columns=["pathway_id", "description"])

    kos_df = pd.DataFrame(
        [(gene.simplified_entry(), gene.name) for gene in genes],
        columns=["ko_id", "description"])

    ko_to_path = []
    ko_to_module = []
    for gene in genes:
        simp = gene.simplified_entry()
        ko_to_path.extend((simp, path.simplified_entry())
                          for path in gene.pathways)
        ko_to_module.extend((simp, hsh_modules[mod].simplified_entry())
                            for mod in gene.modules)

    ko_to_path_df = pd.DataFrame(ko_to_path, columns=["ko_id", "pathway_id"])
    ko_to_module_df = pd.DataFrame(ko_to_module, columns=["ko_id", "module_id"])
    return module_classes_df, modules_df, pathways_df, kos_df, ko_to_path_df, ko_to_module_df


def main():
    parser = argparse.ArgumentParser(description='Prepare KO reference tables')
    parser.parse_args()

    module_classes_df, modules_df, pathways_df, kos_df, ko_to_path_df, ko_to_module_df = load_KO_references()
    module_classes_df.to_csv("ko_module_classes.csv", index=False)
    modules_df.to_csv("ko_modules.csv", index=False)
    pathways_df.to_csv("ko_pathways.csv", index=False)
    kos_df.to_csv("ko_definitions.csv", index=False)
    ko_to_path_df.to_csv("ko_to_pathway.csv", index=False)
    ko_to_module_df.to_csv("ko_to_module.csv", index=False)


if __name__ == '__main__':
    main()
