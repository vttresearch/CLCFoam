import os
from foam.common import read_fo_df

def parseParamString(params_string: str) -> dict:
    return {k:v for k,v in (kv.split('=') for kv in params_string.split(','))}

def parseFunctionObjectString(fo_string):
    if '(' in fo_string:
        splitted = fo_string.split('(',1)
        name = splitted[0]
        params_string = splitted[1][:-1]
        return name, parseParamString(params_string)
    else:
        return fo_string, dict()

def find_function_object(ofcase, foname, foparams):
    found = list()
    none_keys = []
    needed_params = foparams.copy()
    for k,v in foparams.items():
        if v is None:
            none_keys.append(k)
            del needed_params[k]
    for fo in ofcase.function_objects:
        # if foname == fo.name:
        #     print(fo.params, needed_params, none_keys)
        #     # print(fo)
        if foname == fo.name \
            and all(item in fo.params.items() for item in needed_params.items()):
            none_pass = True
            for none_key in none_keys:
                if none_key in fo.params.keys():
                    none_pass = False
            if none_pass:
                found.append(fo)
    return found

class FoamCase:
    """OpenFOAM case wrapper, mainly meant for post-processing.
    """
    def __init__(self, path):
        self.path = path

        self._read_mesh()

        self._read_function_objects()
        self._read_faceZones()
        self._read_cellSets()
        self._read_probes()

        self.domain = Region('domain', dict(), self)
    
    def _read_mesh(self):
        self._read_patches()
    
    def _read_patches(self):
        path_mesh = os.path.join(self.path, 'constant', 'polyMesh')
        patches = dict()
        with open(os.path.join(path_mesh, 'boundary'), 'r') as f:
            # TODO: general openfoam file parser 
            # or make use of foamDictionary either by direct calls
            # or using cython compile only the part that is needed
            inside_list = False
            inside_dict = False
            currentEntry = None
            for line in f.readlines():
                if line.startswith(')'):
                    inside_list = False
                if line.strip().startswith('}'):
                    inside_dict = False
                if inside_list and not inside_dict and '}' not in line and '{' not in line:
                    patches[line.strip()] = dict()
                    currentEntry = line.strip()
                if inside_list and inside_dict:
                    splitted = line.strip()[:-1].split(' ', 1)
                    patches[currentEntry][splitted[0]] = splitted[1].strip()
                if line.startswith('('):
                    inside_list = True
                if line.strip().startswith('{'):
                    inside_dict = True
        self.patches = {k: Patch(k, v, self) for k,v in patches.items()}
    
    def _read_faceZones(self):
        path_mesh = os.path.join(self.path, 'constant', 'polyMesh')
        path_faceZones = os.path.join(path_mesh, 'faceZones')
        faceZones = dict()
        if not os.path.isfile(path_faceZones):
            self.faceZones = dict()
            return
        with open(path_faceZones, 'rb') as f:
            # TODO: general openfoam file parser 
            # or make use of foamDictionary either by direct calls
            # or using cython compile only the part that is needed
            inside_main_list = False
            inside_dict = False
            currentEntry = None
            last_byte = None
            while(last_byte != b''):
                line = f.readline().decode('ascii').strip()
                if line.startswith('faceLabels'):
                    arrlen = int(f.readline().decode('ascii').strip())
                    assert f.read(1).decode('ascii') == '('
                    data = f.read(4*arrlen)
                    assert f.read(1).decode('ascii') == ')'
                if line.startswith('flipMap'):
                    arrlen = int(f.readline().decode('ascii').strip())
                    assert f.read(1).decode('ascii') == '('
                    data = f.read(1*arrlen)
                    assert f.read(1).decode('ascii') == ')'

                if line.startswith(')'):
                    if not inside_dict:
                        inside_main_list = False
                        # break
                if line.strip().startswith('}'):
                    inside_dict = False
                if inside_main_list and not inside_dict and '}' not in line and '{' not in line and len(line.strip())>0:
                    faceZones[line.strip()] = FaceZone(line.strip(), dict(), self)
                    currentEntry = line.strip()
                if line.startswith('('):
                    if not inside_main_list:
                        inside_main_list = True
                if line.strip().startswith('{'):
                    inside_dict = True
                if f.tell() == os.fstat(f.fileno()).st_size:
                    break
        self.faceZones = faceZones
    
    def _read_cellSets(self):
        path = os.path.join(self.path,'constant', 'polyMesh', 'sets')
        cellSets = dict()
        for fname in os.listdir(path):
            fpath = os.path.join(path, fname)
            if os.path.isfile(fpath):
                with open(fpath, 'r') as f:
                    for l in f:
                        if 'class' in l:
                            stype = l.split()[1][:-1]
                            if stype == 'cellSet':
                                cellSets[fname] = CellSet(fname, dict(), self)
        self.cellSets = cellSets

        # class       cellSet;
    def _read_probes(self):
        self.probes = dict()
        for fo in self.function_objects:
            if fo.op_type == 'probe':
                self.probes[fo.name] = fo

    def _read_function_objects(self):
        self.function_objects = list()
        for d in os.scandir(os.path.join(self.path,'postProcessing')):
            f = FunctionObject(d.name, self)
            self.function_objects.append(f)
            # print(f.params)
            # if 'patch' in f.params.keys() and f.params['patch'] == 'inlet_FR':
            #     print(d.name)

def get_function_object(ofcase:FoamCase, foname:str, foparamsd:dict):
    l = find_function_object(ofcase, foname, foparamsd)
    if len(l) > 1:
        print(f'Warning: {foname} function object with params {foparamsd} has multiple matches')
        for fo in l:
            print('Search: ' + foname + ' with ' + str(foparamsd))
            print(' '+fo.dirname)
        raise NotImplementedError
    elif len(l) < 1:
        print(f"Warning: {foname} function object with params {foparamsd} not found")
        return None
    else:
        return l[0]

class FunctionObject:
    def __init__(self, dirname:str, ofcase: FoamCase):
        self.dirname = dirname
        self.name, self.params = parseFunctionObjectString(dirname)
        self.ofcase = ofcase
        self.op_type = None
        
        self._determine_data_type()

    def _determine_data_type(self): # TODO: should be childrem, generated by factory
        # foam_path = '/u1/common/OpenFOAM/OpenFOAM-TASIFNG-20240119/'
        # foam_etc_pp_path = foam_path + 'etc/caseDicts/postProcessing/'
        # surfaceFO_path = foam_etc_pp_path + 'surfaceFieldValue/'
        # volFO_path     = foam_etc_pp_path + 'volFieldValue/'

        # surfaceFOs = [f for f in os.listdir(surfaceFO_path) if os.path.isfile(os.path.join(surfaceFO_path, f)) and '.cfg' not in f]
        # volFOs     = [f for f in os.listdir(volFO_path) if os.path.isfile(os.path.join(volFO_path, f)) and '.cfg' not in f]

        surfaceFOs = ['faceZoneAverage', 'faceZoneFlowRate', 'patchAverage', \
            'patchDifference', 'patchFlowRate', 'patchIntegrate', \
            'triSurfaceAverage', 'triSurfaceDifference', \
            'triSurfaceVolumetricFlowRate']
        volFOs = ['cellMax', 'cellMaxMag', 'cellMin', 'cellMinMag', 'volAverage', 'volIntegrate']

        if self.name in surfaceFOs:
            self.op_type = 'surfaceFieldValue'
            self.filename = 'surfaceFieldValue.dat'
        elif self.name in volFOs:
            self.op_type = 'volFieldValue'
            self.filename = 'volFieldValue.dat'
        elif 'probe' in self.name.lower():
            self.op_type = 'probe'
        #     # TODO: determine filename

        self.op_name = None
        path = os.path.join(self.ofcase.path, 'postProcessing', self.dirname, self.times()[0])
        files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
        if len(files) == 1:
            self.filename = files[0]
        else:
            self.filename = None
    #     if 'probe' in self.name.lower():
    #         self._configure_probe()
        
    # def _configure_probe(self):
    #     path = os.path.join(self.ofcase.path, 'postProcessing', self.dirname, self.times()[0])
    #     files = [f for f in os.listdir(path) if os.path.isfile(join(path, f))]
    #         if len(files) == 1:
    #             self.filename = files[0]
            

    def times(self):
        path_pp_fo = os.listdir(os.path.join(self.ofcase.path, 'postProcessing', self.dirname))
        times = sorted(path_pp_fo, key=float)
        return times

    def get_data(self):
        return read_fo_df(
            os.path.join(self.ofcase.path),
            self.dirname,
            self.filename,
            self.times()
            )
    
    def __repr__(self):
        return f'FunctionObject({self.name}, {self.params})'

class MeshObject:
    def __init__(self, name, params, ofcase):
        self.name = name
        self.params = params
        self.ofcase = ofcase
        self.hefield = 'h'
        self.heinterpfield = '{0}.{1}' # first param is hefield, second phase
        # # TODO: generalize for other solvers (without phases)
        # # self.typename - TO BE SET BY CHILD CLASSES
    
    def register_custom_function_object(self, method_name, fo_name, params={}):
        def method(self):
            d = dict() if self.typename is None or self.typename in params.keys() else {self.typename: self.name}
            d = d | params
            if get_function_object(self.ofcase, fo_name, d) is None:
                return None
            else:
                return get_function_object(self.ofcase, fo_name, d).get_data()
        setattr(self, method_name, method.__get__(self, self.__class__))

class Surface(MeshObject):
    """Abstract class for all surfaces"""
    def __init__(self, name, params, ofcase):
        super().__init__(name, params, ofcase)

    def mdot(self, phase=None):
        d = {self.typename: self.name, 'weightFields': None}
        if phase:
            d['fields'] = f'(alphaRhoPhi.{phase})'
        return get_function_object(self.ofcase, f'{self.typename}FlowRate', d).get_data()
    
    def hdot(self, phase=None):
        d = {self.typename: self.name}
        if phase:
            d['fields'] = f'(alphaRhoPhi.{phase})'
            d['weightFields'] = '('+self.heinterpfield.format(self.hefield, phase)+')'
        return get_function_object(self.ofcase, f'{self.typename}FlowRate', d).get_data() 

class Patch(Surface):
    def __init__(self, name, params, ofcase):
        super().__init__(name, params, ofcase)
        self.typename = 'patch'

class FaceZone(Surface):
    def __init__(self, name, params, ofcase):
        super().__init__(name, params, ofcase)
        self.typename = 'faceZone'
        self.heinterpfield = 'surfaceInterpolate({0}.{1})'
        # self.he = (surfaceInterpolate(h.gas))

class Volume(MeshObject):
    def __init__(self, name, params, ofcase):
        super().__init__(name, params, ofcase)
    
    def m(self, phase):
        d = {'select': self.typename, 'fields': f'(rho.{phase})', 'weightField': f'alpha.{phase}'}
        if not self.typename is None:
            d[self.typename] = self.name
        return get_function_object(self.ofcase, 'volIntegrate', d).get_data()
    
    def Tavgest(self, phase):
        d = {'select': self.typename, 'fields': f'(T.{phase})', 'weightFields': f'(rho.{phase}alpha.{phase})'}
        if not self.typename is None:
            d[self.typename] = self.name
        return get_function_object(self.ofcase, 'volAverage', d).get_data()

class Region(Volume):
    def __init__(self, name, params, ofcase):
        super().__init__(name, params, ofcase)
        self.typename = None
    
    def register_custom_function_object(self, method_name, fo_name, params={}):
        super().register_custom_function_object(method_name, fo_name, params | {'select':None})
    
class CellSet(Volume):
    def __init__(self, name, params, ofcase):
        super().__init__(name, params, ofcase)
        self.typename = 'cellSet'
