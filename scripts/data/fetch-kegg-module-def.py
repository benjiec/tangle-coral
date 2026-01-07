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
    option     = component (plusminus component)*
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
        component, more_components = visited_children
        components.append(("+", component))
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


def print_option(tree, spacing):
    for component in tree.components:
        if type(component[1]) == type(""):
            print(" "*spacing, "component", component[0], component[1])
        else:
            print(" "*spacing, "component", component[0])
            if type(component[1]) is Pathway:
                print_pathway(component[1], spacing+4)
            elif type(component[1]) is Step:
                print_step(component[1], spacing+4)


def print_step(tree, spacing):
    for option in tree.options:
        print(" "*spacing, "option")
        print_option(option, spacing+4)


def print_pathway(tree, spacing=None):
    spacing = 0 if spacing is None else spacing
    for step in tree.steps:
        print(" "*spacing, "step")
        print_step(step, spacing+4)


def parse_module_definition(def_line):
    print(def_line)
    tree = grammar.parse(def_line)
    visitor = DefVisitor()
    output = visitor.visit(tree)
    if type(output) is Step:
        if len(output.options) == 1 and len(output.options[0].components) == 1 and type(output.options[0].components[0][1]) is Pathway:
            output = output.options[0].components[0][1]
        else:
            output = Pathway(steps=[output])
    print_pathway(output)


def parse_module_data(module_data):
    keyword = "DEFINITION"
    lines = module_data.split("\n")
    def_lines = [line for line in lines if line.startswith(keyword)]
    if len(def_lines) != 1:
        return None
    def_line = def_lines[0][len(keyword):]
    def_line = def_line.strip()
    return parse_module_definition(def_line)


"""
with open('data/modules.tsv') as f:
    for line in f.readlines():
        if line.startswith("Module ID"):
            continue
        module_id, module_name = line.split("\t")
        if module_id not in ("M00009"):
            continue
        print(module_id)
        module_data = fetch_module_data(module_id)
        parse_module_data(module_data)
"""

parse_module_definition("((K00134,K00150) K00927,K11389)")

parse_module_definition("(K00134,K00150) K00927,K11389")

parse_module_definition("(K00134,K00150,K00927,K11389)")

parse_module_definition("(K01647,K05942,K01659) (K01681,K27802,K01682) (K00031,K00030) ((K00164+K00658,K01616)+K00382,K00174+K00175-K00177-K00176) (K01902+K01903,K01899+K01900,K18118) (K00234+K00235+K00236+(K00237,K25801),K00239+K00240+K00241-(K00242,K18859,K18860),K00244+K00245+K00246-K00247) (K01676,K01679,K01677+K01678) (K00026,K00025,K00024,K00116)")

parse_module_definition("(K00844,K12407,K00845,K25026,K00886,K08074,K00918) (K01810,K06859,K13810,K15916) (K00850,K16370,K21071,K24182,K00918) (K01623,K01624,K11645,K16305,K16306) K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) (K01689,K27394) (K00873,K12406)")
