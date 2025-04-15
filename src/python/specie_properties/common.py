import cantera as ct
from math import isclose

def formulaToCompositionDict(formula):
    if formula == 'AR':
        return dict(Ar=1)
    import re
    l = re.findall('[A-Z][a-z]?|\d+|.', formula)
    d = dict()
    i=0
    while i < len(l):
        if (i+1<len(l) and l[i+1].isdigit()):
            d[l[i]] = int(l[i+1])
            i+=1
        else:
            d[l[i]] = 1
        i+=1
    return d

assert formulaToCompositionDict('FeO') == dict(Fe=1,O=1)
assert formulaToCompositionDict('Fe2O3') == dict(Fe=2,O=3)
assert formulaToCompositionDict('CH4') == dict(C=1,H=4)
assert formulaToCompositionDict('O2') == dict(O=2)
assert formulaToCompositionDict('AR') == dict(Ar=1)
assert formulaToCompositionDict('Ar') == dict(Ar=1)

def AW(element):
    """Atomic weight in kg/mol"""
    return ct.Element(element).weight/1000
assert isclose(AW('Ar'), 39.95e-3)

def MW(formula):
    """Molar weight of a specie in kg/mol"""
    composition = formulaToCompositionDict(formula)
    s = ct.Solution()
    test = ct.Species.from_dict(dict(name=formula, composition=composition))
    s.add_species(test)
    return s.molecular_weights[0]/1000

assert isclose(MW('Fe2O3'), 159.687e-3)
assert isclose(MW('FeO'), 71.844e-3)
assert isclose(MW('N2'), 28.014e-3)
assert isclose(MW('O2'), 31.998e-3)
assert isclose(MW('CO'), 28.01e-3)
assert isclose(MW('CO2'), 44.009e-3)
assert isclose(MW('H2O'), 18.015e-3)
assert isclose(MW('H2'), 2.016e-3)
assert isclose(MW('CH4'), 16.043e-3)
assert isclose(MW('AR'), 39.95e-3)
