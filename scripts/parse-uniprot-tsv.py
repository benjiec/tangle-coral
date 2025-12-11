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
        self.counts = {}

    def _load_uniprot_rows(self, path):
        with Path(path).open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for hdr in HEADERS:
                assert hdr in reader.fieldnames
            return [row for row in reader]

    def _add_lineage(self, tax_lineage, tax_id_lineage):

        tax_names = [n.strip() for n in tax_lineage.split("),") if n.strip()]
        tax_names = [n if n.endswith(')') else n+')' for n in tax_names]
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

        if last_id:
            if last_id not in self.counts:
                self.counts[last_id] = 0
            self.counts[last_id] += 1

    def load_uniprot_tsv(self, path):
        rows = self._load_uniprot_rows(path)
        for i,row in enumerate(rows):
            self._add_lineage(row[HDR_TAX_LINEAGE], row[HDR_TAXID_LINEAGE])
            if i%10000 == 0:
                print(i)

    def tally(self):

        def tally_at_node(node):
            n = self.counts[node.taxon_id] if node.taxon_id in self.counts else 0
            if node._children is None:
                return n
            n = n + sum([tally_at_node(c) for c in node._children])
            self.counts[node.taxon_id] = n
            return n

        for n in self.root:
            tally_at_node(n)


if __name__ == "__main__":

    ap = argparse.ArgumentParser()
    ap.add_argument("uniprot_tsv")
    ap.add_argument("output_tsv")
    args = ap.parse_args()

    tree = Tree()
    tree.load_uniprot_tsv(args.uniprot_tsv)

    print("root", len(tree.root))
    print("size", len(tree.nodes.keys()))

    if len(tree.root):
        assert 0 not in tree.nodes
        real_root = Node(taxon_id=0, parent_taxon_id=None, name="root")
        for n in tree.root:
            n.parent_taxon_id = 0
            n._parent = real_root
        real_root._children = tree.root
        tree.nodes[0] = real_root
        tree.root = [real_root]

    tree.tally()

    with Path(args.output_tsv).open("w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=["node", "parent", "nproteins"])
        writer.writeheader()
        for taxon_id, count in tree.counts.items():
            node = tree.nodes[taxon_id]
            if node.parent_taxon_id is None:
                parent_name = None
            else:
                parent_name = tree.nodes[node.parent_taxon_id].name+f" ({node.parent_taxon_id})"
            d = dict(node=node.name+f" ({taxon_id})",
                     parent=parent_name,
                     nproteins=count)
            writer.writerow(d)
