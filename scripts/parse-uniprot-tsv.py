import re
import csv
import argparse
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


HDR_ENTRY = "Entry"
HDR_NAME = "Entry Name"
HDR_PROTEIN_NAME = "Protein names"
HDR_ORGANISM = "Organism"
HDR_EC_NUMBER = "EC number"
HDR_TAX_LINEAGE = "Taxonomic lineage"
HDR_TAXID_LINEAGE = "Taxonomic lineage (Ids)"

HEADERS = [
  HDR_ENTRY,
  HDR_NAME,
  HDR_PROTEIN_NAME,
  HDR_ORGANISM,
  HDR_EC_NUMBER,
  HDR_TAX_LINEAGE,
  HDR_TAXID_LINEAGE
]


@dataclass
class Node(object):
    taxon_id: int
    parent_taxon_id: int
    name: str

    _parent: Optional[object] = None
    _children: Optional[List[object]] = None


class Tree(object):

    def __init__(self):
        self.nodes = {}
        self.root = []

    def _load_uniprot_rows(self, path):
        with Path(path).open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for hdr in HEADERS:
                assert hdr in reader.fieldnames
            return [row for row in reader]

    def _add_lineage(self, tax_lineage, tax_id_lineage):

        tax_names = [n.strip() for n in tax_lineage.split("),") if n.strip()]
        tax_ids = [int(re.match(r'(\d+)', n.strip())[1]) for n in tax_id_lineage.split("),") if n.strip()]
        assert len(tax_names) == len(tax_ids)
     
        last_id = None 
        for name, id in zip(tax_names, tax_ids):
            if id not in self.nodes:
                node = Node(taxon_id=id, parent_taxon_id=last_id, name=name)
                self.nodes[id] = node
                if last_id:
                    self.nodes[id]._parent = self.nodes[node.parent_taxon_id]
                    if self.nodes[node.parent_taxon_id]._children is None:
                        self.nodes[node.parent_taxon_id]._children = []
                    self.nodes[node.parent_taxon_id]._children.append(node)
                else:
                    self.root.append(node)
            last_id = id

    def load_uniprot_tsv(self, path):
        rows = self._load_uniprot_rows(path)
        for i,row in enumerate(rows):
            self._add_lineage(row[HDR_TAX_LINEAGE], row[HDR_TAXID_LINEAGE])
            if i%10000 == 0:
                print(i)


if __name__ == "__main__":

    ap = argparse.ArgumentParser()
    ap.add_argument("uniprot_tsv")
    args = ap.parse_args()

    tree = Tree()
    tree.load_uniprot_tsv(args.uniprot_tsv)

    print("root", len(tree.root))
    print("size", len(tree.nodes.keys()))
