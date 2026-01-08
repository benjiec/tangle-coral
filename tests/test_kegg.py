import os
import tempfile
import unittest

from needle.kegg import parse_module_definition, ModuleDefinitionFormatter

class TestKEGGModuleDefinition(unittest.TestCase):

    def setUp(self):
        tmpf = tempfile.NamedTemporaryFile(delete=False, suffix=".csv", mode="w")
        tmpf.close()
        self.tmpfn = tmpf.name

    def tmpData(self):
        with open(self.tmpfn, "r") as f:
            data = f.read().split("\n")
            lines = [l.strip().split(",") for l in data if l.strip()]
            return lines

    def tearDown(self):
        os.remove(self.tmpfn)

    def test_single_level_steps(self):
        fmt = ModuleDefinitionFormatter()
        parse_module_definition(fmt, "x1", "K00134 K00150")
        fmt.to_csv(self.tmpfn)

        lines = self.tmpData()
        self.assertEqual(len(lines), 3)
        # just single level, each line has module ID, step ID, option ID, component ID, then essential and identifier
        self.assertEqual(len(lines[0]), 6)
        self.assertEqual(lines[1], ["x1", "1", "1", "1", "1", "K00134"])
        self.assertEqual(lines[2], ["x1", "2", "1", "1", "1", "K00150"])

    def test_single_level_options(self):
        fmt = ModuleDefinitionFormatter()
        parse_module_definition(fmt, "x1", "K00134,K00150")
        fmt.to_csv(self.tmpfn)

        lines = self.tmpData()
        self.assertEqual(len(lines), 3)
        # just single level, each line has module ID, step ID, option ID, component ID, then essential and identifier
        self.assertEqual(len(lines[0]), 6)
        self.assertEqual(lines[1], ["x1", "1", "1", "1", "1", "K00134"])
        self.assertEqual(lines[2], ["x1", "1", "2", "1", "1", "K00150"])

    def test_single_level_components(self):
        fmt = ModuleDefinitionFormatter()
        parse_module_definition(fmt, "x1", "K00134+K00150")
        fmt.to_csv(self.tmpfn)
        lines = self.tmpData()
        self.assertEqual(len(lines), 3)
        # just single level, each line has module ID, step ID, option ID, component ID, then essential and identifier
        self.assertEqual(len(lines[0]), 6)
        self.assertEqual(lines[1], ["x1", "1", "1", "1", "1", "K00134"])
        self.assertEqual(lines[2], ["x1", "1", "1", "2", "1", "K00150"])

        fmt = ModuleDefinitionFormatter()
        parse_module_definition(fmt, "x1", "K00134-K00150")
        fmt.to_csv(self.tmpfn)
        lines = self.tmpData()
        self.assertEqual(len(lines), 3)
        # just single level, each line has module ID, step ID, option ID, component ID, then essential and identifier
        self.assertEqual(len(lines[0]), 6)
        self.assertEqual(lines[1], ["x1", "1", "1", "1", "1", "K00134"])
        self.assertEqual(lines[2], ["x1", "1", "1", "2", "0", "K00150"])

        fmt = ModuleDefinitionFormatter()
        parse_module_definition(fmt, "x1", "-K00134-K00150")
        fmt.to_csv(self.tmpfn)
        lines = self.tmpData()
        self.assertEqual(len(lines), 3)
        # just single level, each line has module ID, step ID, option ID, component ID, then essential and identifier
        self.assertEqual(len(lines[0]), 6)
        self.assertEqual(lines[1], ["x1", "1", "1", "1", "0", "K00134"])
        self.assertEqual(lines[2], ["x1", "1", "1", "2", "0", "K00150"])

    def test_multiple_steps_and_options(self):
        fmt = ModuleDefinitionFormatter()
        parse_module_definition(fmt, "x2", "(K00134,K00150) K00927,K11389-K11400")
        fmt.to_csv(self.tmpfn)

        lines = self.tmpData()
        self.assertEqual(len(lines), 6)
        # just single level, each line has module ID, step ID, option ID, component ID, then essential and identifier
        self.assertEqual(len(lines[0]), 6)
        self.assertEqual(lines[1], ["x2", "1", "1", "1", "1", "K00134"])
        self.assertEqual(lines[2], ["x2", "1", "2", "1", "1", "K00150"])
        self.assertEqual(lines[3], ["x2", "2", "1", "1", "1", "K00927"])
        self.assertEqual(lines[4], ["x2", "2", "2", "1", "1", "K11389"])
        self.assertEqual(lines[5], ["x2", "2", "2", "2", "0", "K11400"])

    def test_extra_grouping_ignored(self):
        fmt = ModuleDefinitionFormatter()
        parse_module_definition(fmt, "x2", "K00134,K00150 K00927,K11389")
        fmt.to_csv(self.tmpfn)
        lines_simple = self.tmpData()

        fmt = ModuleDefinitionFormatter()
        parse_module_definition(fmt, "x2", "(K00134,K00150) K00927,K11389")
        fmt.to_csv(self.tmpfn)
        self.assertEqual(self.tmpData(), lines_simple)

        fmt = ModuleDefinitionFormatter()
        parse_module_definition(fmt, "x2", "(K00134,K00150 K00927,K11389)")
        fmt.to_csv(self.tmpfn)
        self.assertEqual(self.tmpData(), lines_simple)

    def test_two_levels(self):
        fmt = ModuleDefinitionFormatter()
        parse_module_definition(fmt, "x3", "(K00134+K00150),K00161 K00927,(K11389-K11400,K11401 K11432)")
        fmt.to_csv(self.tmpfn)

        lines = self.tmpData()
        self.assertEqual(len(lines), 9)
        # two levels, each line has module ID, two of "step ID, option ID, component ID", then essential and identifier
        self.assertEqual(len(lines[0]), 9)
        # first step, option 1 is split into two components
        self.assertEqual(lines[1], ["x3", "1", "1", "1", "1", "1", "1", "1", "K00134"])
        self.assertEqual(lines[2], ["x3", "1", "1", "1", "1", "1", "2", "1", "K00150"])
        self.assertEqual(lines[3], ["x3", "1", "2", "1", "",  "",  "",  "1", "K00161"])
        self.assertEqual(lines[4], ["x3", "2", "1", "1", "",  "",  "",  "1", "K00927"])
        # second step, option 2 is itself a two step process
        self.assertEqual(lines[5], ["x3", "2", "2", "1", "1", "1", "1", "1", "K11389"])
        self.assertEqual(lines[6], ["x3", "2", "2", "1", "1", "1", "2", "0", "K11400"])
        #  .. and this two step process has two options in the first step
        self.assertEqual(lines[7], ["x3", "2", "2", "1", "1", "2", "1", "1", "K11401"])
        self.assertEqual(lines[8], ["x3", "2", "2", "1", "2", "1", "1", "1", "K11432"])

    def test_three_levels(self):
        fmt = ModuleDefinitionFormatter()
        parse_module_definition(fmt, "x3", "(K00134+K00150),K00161 K00927,((K11389 -K11400),K11401 K11432)")
        fmt.to_csv(self.tmpfn)

        lines = self.tmpData()
        self.assertEqual(len(lines), 9)
        # two levels, each line has module ID, three of "step ID, option ID, component ID", then essential and identifier
        self.assertEqual(len(lines[0]), 12)
        # first step, option 1 is split into two components
        self.assertEqual(lines[1], ["x3", "1", "1", "1", "1", "1", "1", "", "", "", "1", "K00134"])
        self.assertEqual(lines[2], ["x3", "1", "1", "1", "1", "1", "2", "", "", "", "1", "K00150"])
        self.assertEqual(lines[3], ["x3", "1", "2", "1", "",  "",  "",  "", "", "", "1", "K00161"])
        self.assertEqual(lines[4], ["x3", "2", "1", "1", "",  "",  "",  "", "", "", "1", "K00927"])
        # second step, option 2 is itself a two step process
        #  .. and first sub step, first option, is itself a two step process
        self.assertEqual(lines[5], ["x3", "2", "2", "1", "1", "1", "1", "1", "1", "1", "1", "K11389"])
        self.assertEqual(lines[6], ["x3", "2", "2", "1", "1", "1", "1", "2", "1", "1", "0", "K11400"])
        self.assertEqual(lines[7], ["x3", "2", "2", "1", "1", "2", "1", "",  "",  "",  "1", "K11401"])
        self.assertEqual(lines[8], ["x3", "2", "2", "1", "2", "1", "1", "",  "",  "",  "1", "K11432"])
