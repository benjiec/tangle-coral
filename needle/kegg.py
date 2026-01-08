from dataclasses import dataclass
from typing import List
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


class DefinitionVisitor(NodeVisitor):
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


class ModuleDefinitionFormatter(object):

    def __init__(self):
        self.row_module = []
        self.row_coords = []
        self.row_values = []

    def _gather_option(self, module_id, step_num, option_num, par_coords, tree, spacing):
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
                    self._gather_step(module_id, 1, coords, component[1], spacing+4)

    def _gather_step(self, module_id, step_num, par_coords, tree, spacing):
        for i, option in enumerate(tree.options):
            # print(" "*spacing, "option", i+1)
            self._gather_option(module_id, step_num, i+1, par_coords, option, spacing+4)

    def gather_pathway(self, module_id, par_coords, tree, spacing=None):
        spacing = 0 if spacing is None else spacing
        for i, step in enumerate(tree.steps):
            # print(" "*spacing, "step", i+1)
            self._gather_step(module_id, i+1, par_coords, step, spacing+4)

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
    # print(def_line)
    def_line = def_line.replace("--", "").strip()
    tree = grammar.parse(def_line)
    visitor = DefinitionVisitor()
    output = visitor.visit(tree)
    if type(output) is Step:
        if len(output.options) == 1 and len(output.options[0].components) == 1 and type(output.options[0].components[0][1]) is Pathway:
            output = output.options[0].components[0][1]
        else:
            output = Pathway(steps=[output])
    formatter.gather_pathway(module_id, [], output)
