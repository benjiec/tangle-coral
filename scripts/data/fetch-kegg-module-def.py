import requests
from needle.kegg import parse_module_definition, ModuleDefinitionFormatter


def fetch_module_data(module_id):
    url = f"https://rest.kegg.jp/get/{module_id}"
    module_data = requests.get(url).text
    return module_data

def parse_module_data(formatter, module_id, module_data):
    keyword = "DEFINITION"
    lines = module_data.split("\n")
    def_lines = [line for line in lines if line.startswith(keyword)]
    if len(def_lines) != 1:
        return None
    def_line = def_lines[0][len(keyword):]
    def_line = def_line.strip()
    return parse_module_definition(formatter, module_id, def_line)

formatter = ModuleDefinitionFormatter()

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

formatter.to_csv('data/module_defs.csv')
