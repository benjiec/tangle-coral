import requests
from dataclasses import dataclass
from typing import List, Optional, Tuple
from parsimonious.grammar import Grammar
from parsimonious.nodes import NodeVisitor


grammar = Grammar(
    r"""
    pathway    = step (ws step)*
    step       = options
    options    = option (comma option)*
    option     = plusminus? component (plusminus component)*
    component  = grouped / atom
    grouped    = "(" (pathway / options) ")"

    atom       = ~"[KRM][0-9]{5}"
    comma      = ","
    plusminus  = plus / minus
    plus       = "+"
    minus      = "-"
    ws         = ~r"\s+"
    """
)


@dataclass
class Pathway:
    steps: List

@dataclass
class Step:
   options: List

@dataclass
class Option:
   components: List


class DefVisitor(NodeVisitor):
    def visit_pathway(self, node, visited_children):
        steps = []
        step, more_steps = visited_children
        steps.append(step)
        for _, step in more_steps:
            steps.append(step)
        if len(steps) > 1:
            return Pathway(steps=steps)
        else:
            return steps[0]

    def visit_step(self, node, visited_children):
        return visited_children[0]

    def visit_grouped(self, node, visited_children):
        _1, options_or_pathway, _2 = visited_children
        return options_or_pathway[0]

    def visit_options(self, node, visited_children):
        options = []
        option, more_options = visited_children
        options.append(option)
        for _, option in more_options:
            options.append(option)
        if len(options) == 1 and len(options[0].components) == 1 and type(options[0].components[0][1]) is Step:
            return options[0].components[0][1]
        return Step(options=options)

    def visit_option(self, node, visited_children):
        components = []
        plus_minus, component, more_components = visited_children
        if len(plus_minus) == 0:
            components.append(("+", component))
        else:
            components.append((plus_minus[0], component))
        for modifier, component in more_components:
            components.append((modifier[0], component))
        return Option(components=components)

    def visit_component(self, node, visited_children):
        return visited_children[0]

    def generic_visit(self, node, visited_children):
        return visited_children or node.text


def fetch_module_data(module_id):
    url = f"https://rest.kegg.jp/get/{module_id}"
    module_data = requests.get(url).text
    return module_data


class DefFormatter(object):

    def __init__(self):
        self.row_module = []
        self.row_coords = []
        self.row_values = []

    def gather_option(self, module_id, step_num, option_num, par_coords, tree, spacing):
        for i, component in enumerate(tree.components):
            coords = par_coords+[[step_num, option_num, i+1]]

            if type(component[1]) == type(""):
                # print(" "*spacing, "component", component[0], component[1])
                self.row_module.append(module_id)
                self.row_coords.append(coords)
                self.row_values.append([1 if component[0] == "+" else 0, component[1]])

            else:
                # print(" "*spacing, "component", component[0], type(component[1]))
                if type(component[1]) is Pathway:
                    self.gather_pathway(module_id, coords, component[1], spacing+4)
                elif type(component[1]) is Step:
                    self.gather_step(module_id, 1, coords, component[1], spacing+4)

    def gather_step(self, module_id, step_num, par_coords, tree, spacing):
        for i, option in enumerate(tree.options):
            # print(" "*spacing, "option", i+1)
            self.gather_option(module_id, step_num, i+1, par_coords, option, spacing+4)

    def gather_pathway(self, module_id, par_coords, tree, spacing=None):
        spacing = 0 if spacing is None else spacing
        for i, step in enumerate(tree.steps):
            # print(" "*spacing, "step", i+1)
            self.gather_step(module_id, i+1, par_coords, step, spacing+4)

    def to_csv(self, fn):
        nlevels = max([len(x) for x in self.row_coords])
       
        headers = []
        headers.append("module_id")
        for i in range(nlevels):
            headers.append("sub_"*i+"step")
            headers.append("sub_"*i+"step_option")
            headers.append("sub_"*i+"step_option_component")
        headers.append("essential")
        headers.append("identifier")

        with open(fn, "w") as f:
            f.write(",".join(headers)+"\n")
            for module_id,coords,values in zip(self.row_module, self.row_coords, self.row_values):
                row = [module_id]
                for i in range(nlevels):
                    if i < len(coords):
                        row.extend(coords[i])
                    else:
                        row.extend(["", "", ""])
                row.extend(values)
                f.write(",".join([str(x) for x in row])+"\n")


def parse_module_definition(formatter, module_id, def_line):
    print(def_line)
    def_line = def_line.replace("--", "").strip()
    tree = grammar.parse(def_line)
    visitor = DefVisitor()
    output = visitor.visit(tree)
    if type(output) is Step:
        if len(output.options) == 1 and len(output.options[0].components) == 1 and type(output.options[0].components[0][1]) is Pathway:
            output = output.options[0].components[0][1]
        else:
            output = Pathway(steps=[output])
    formatter.gather_pathway(module_id, [], output)


def parse_module_data(formatter, module_id, module_data):
    keyword = "DEFINITION"
    lines = module_data.split("\n")
    def_lines = [line for line in lines if line.startswith(keyword)]
    if len(def_lines) != 1:
        return None
    def_line = def_lines[0][len(keyword):]
    def_line = def_line.strip()
    return parse_module_definition(formatter, module_id, def_line)


formatter = DefFormatter()

with open('data/modules.tsv') as f:
    for line in f.readlines():
        if line.startswith("Module ID"):
            continue
        module_id, module_name = line.split("\t")
        if module_id[0] != "M":
            continue
        print(module_id)
        module_data = fetch_module_data(module_id)
        parse_module_data(formatter, module_id, module_data)

"""
parse_module_definition(formatter, "x1", "((K00134,K00150) K00927,K11389)")
parse_module_definition(formatter, "x2", "(K00134,K00150) K00927,K11389")
parse_module_definition(formatter, "x3", "(K00134,K00150,K00927,K11389)")
parse_module_definition(formatter, "x4", "(K01647,K05942,K01659) (K01681,K27802,K01682) (K00031,K00030) ((K00164+K00658,K01616)+K00382,K00174+K00175-K00177-K00176) (K01902+K01903,K01899+K01900,K18118) (K00234+K00235+K00236+(K00237,K25801),K00239+K00240+K00241-(K00242,K18859,K18860),K00244+K00245+K00246-K00247) (K01676,K01679,K01677+K01678) (K00026,K00025,K00024,K00116)")
parse_module_definition(formatter, "x5", "(K00844,K12407,K00845,K25026,K00886,K08074,K00918) (K01810,K06859,K13810,K15916) (K00850,K16370,K21071,K24182,K00918) (K01623,K01624,K11645,K16305,K16306) K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) (K01689,K27394) (K00873,K12406)")
parse_module_definition(formatter, "x", "(K01596,K01610) (K01689,K27394) (K01834,K15633,K15634,K15635) K00927 (K00134,K00150) K01803 ((K01623,K01624,K11645) (K03841,K02446,K11532,K01086,K04041),K01622)")
"""

formatter.to_csv('data/module_defs.tsv')
