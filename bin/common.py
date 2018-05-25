import re

class DiffusionSystem:
    def __init__(self, systemString):
        output_name_regex = re.compile("\.?/?([0-9])sheet/([a-z]+)/([0-9]+)pore/output_bf4(_tma)?_tmea_([A-Z0-9]+)diff.*") 
        sim_path_regex = re.compile("bf4(_tma)?_tmea_(dce|acn)_(1|2)_(7|10|14)_([A-Z0-9]+)diff.*") 

        match_out = output_name_regex.match(systemString)
        match_path = sim_path_regex.match(systemString)

        if match_out is not None:
            self.sheetCount = match_out.group(1)     
            self.solvent = match_out.group(2)        
            self.poreSize = match_out.group(3)       
            self.tmaPresent = match_out.group(4)
            self.diffusingIon = match_out.group(5)   

        elif match_path is not None:
            self.tmaPresent = match_path.group(1)
            self.solvent = match_path.group(2)
            self.sheetCount = match_path.group(3)
            self.poreSize = match_path.group(4)
            self.diffusingIon = match_path.group(5)

        else:
            raise ValueError("invalid system string", systemString)

        if self.tmaPresent is None:
            self.tmaPresent = ""


    def sim_path(self):
       return "{0:s}sheet/{1:s}/{2:s}pore/output_bf4{3:s}_tmea_{4:s}diff".format(self.sheetCount, self.solvent, self.poreSize, self.tmaPresent, self.diffusingIon)

    def output_name(self):
        return "bf4{0:s}_tmea_{1:s}_{2:s}_{3:s}_{4:s}diff".format(self.tmaPresent, self.solvent, self.sheetCount, self.poreSize, self.diffusingIon)
