import re

class DiffusionSystem:
    def __init__(self, systemString):
        output_name_regex = re.compile("\.?/?([0-9])sheet/([a-z]+)/([0-9]+)pore/output_bf4(_tma)?_tmea_([A-Z0-9]+)diff(_.*)?") 
        sim_path_regex = re.compile("bf4(_tma)?_tmea_(dce|acn)_(1|2)_(7|10|14)_([A-Z0-9]+)diff(_.*)?") 

        match_out = output_name_regex.match(systemString)
        match_path = sim_path_regex.match(systemString)

        if match_out is not None:
            self.sheetCount = match_out.group(1)     
            self.solvent = match_out.group(2)        
            self.poreSize = match_out.group(3)       
            self.tmaPresent = match_out.group(4)
            self.diffusingIon = match_out.group(5)   
            self.tag = match_out.group(6)   

        elif match_path is not None:
            self.tmaPresent = match_path.group(1)
            self.solvent = match_path.group(2)
            self.sheetCount = match_path.group(3)
            self.poreSize = match_path.group(4)
            self.diffusingIon = match_path.group(5)
            self.tag = match_path.group(6)

        else:
            raise ValueError("invalid system string", systemString)

        if self.tmaPresent is None:
            self.tmaPresent = ""
        if self.tag is None:
            self.tag = ""

    def __hash__(self):
        return hash(self.output_name())

    def __eq__(self, other):
        print(dir(self))
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

    def sim_path(self):
        return "{0:s}sheet/{1:s}/{2:s}pore/output_bf4{3:s}_tmea_{4:s}diff{5:s}".format(self.sheetCount, self.solvent, self.poreSize, self.tmaPresent, self.diffusingIon, self.tag)

    def output_name(self):
        return "bf4{0:s}_tmea_{1:s}_{2:s}_{3:s}_{4:s}diff{5:s}".format(self.tmaPresent, self.solvent, self.sheetCount, self.poreSize, self.diffusingIon, self.tag)
