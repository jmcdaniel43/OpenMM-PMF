import re

class DiffusionSystem:
    def __init__(self, systemString):
        output_name_regex = re.compile("\.?/?([0-9])sheet/([a-z0-9]+)/([0-9]+)pore/output_bf4_(.*)_([A-Z0-9]+)diff(_.*)?")
        sim_path_regex = re.compile("bf4_(.*)_(dce|acn|h2o)_(1|2)_(7|10|14)_([A-Z0-9]+)diff(_.*)?")

        match_out = output_name_regex.match(systemString)
        match_path = sim_path_regex.match(systemString)

        if match_out is not None:
            self.sheetCount = match_out.group(1)     
            self.solvent = match_out.group(2)        
            self.poreSize = match_out.group(3)       
            self.ion_pair = match_out.group(4)
            self.diffusingIon = match_out.group(5)   
            self.tag = match_out.group(6)   

        elif match_path is not None:
            self.ion_pair = match_path.group(1)
            self.solvent = match_path.group(2)
            self.sheetCount = match_path.group(3)
            self.poreSize = match_path.group(4)
            self.diffusingIon = match_path.group(5)
            self.tag = match_path.group(6)

        else:
            raise ValueError("invalid system string", systemString)

        if self.tag is None:
            self.tag = ""

    def __hash__(self):
        return hash(self.output_name())

    def __eq__(self, other):
        if isinstance(other, DiffusionSystem):
            return False

        return NotImplemented
    
    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __lt__(self, other):
        return self.output_name() < other.output_name()

    def __gt__(self, other):
        return self.output_name() > other.output_name()

    def __le__(self, other):
        return self.output_name() <= other.output_name()

    def __ge__(self, other):
        return self.output_name() >= other.output_name()

    def __repr__(self):
        return self.output_name()

    def __format__(self, spec):
        return self.__repr__().__format__(spec)

    def sim_path(self):
        return "{0:s}sheet/{1:s}/{2:s}pore/output_bf4_{3:s}_{4:s}diff{5:s}".format(self.sheetCount, self.solvent, self.poreSize, self.ion_pair, self.diffusingIon, self.tag)

    def output_name(self):
        return "bf4_{0:s}_{1:s}_{2:s}_{3:s}_{4:s}diff{5:s}".format(self.ion_pair, self.solvent, self.sheetCount, self.poreSize, self.diffusingIon, self.tag)

    def pandas_tuple(self):
        return (int(self.sheetCount), int(self.poreSize), self.solvent, self.ion_pair, self.diffusingIon)

    """
    Checks for a single difference between systems. For comparing systems based on one variable.

    Returns True if a single difference exists; False otherwise.
    """
    def is_analogue(self, other): # what a gross method impl
        differences_count = 0

        if self.sheetCount != other.sheetCount:
            differences_count += 1
        if self.ion_pair != other.ion_pair:
            differences_count += 1
        if self.solvent != other.solvent:
            differences_count += 1
        if self.poreSize != other.poreSize:
            differences_count += 1
        if self.diffusingIon != other.diffusingIon:
            differences_count += 1

        return differences_count == 1
