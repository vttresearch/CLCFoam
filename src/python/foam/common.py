from os.path import join
from math import isclose
import pandas as pd
import cantera as ct

def is_type(s: str, typ):
    try:
        typ(s)
        return True
    except ValueError:
        return False

def is_int(s: str):
    return is_type(s, int)

def is_float(s: str):
    return is_type(s, float)

def integrate_time_series(ts):
    """same as `integrate.cumtrapz(ts, ts.index, initial=0)`, but preserves index"""
    return ((ts.shift(1) + ts) / 2 * ts.index.to_series().diff() ).cumsum().fillna(0.)

def read_case_params(case_path, params):
    # TODO: create or take somewhere OpenFOAM file reader, using some automated
    # file parsing library preferrably. Or alternatilvely use system call to 
    # foamDictionary -expand
    # The difficulty is that it will be system-dependent, because path to
    # OpenFOAM needs to be provided
    d = dict()
    for param_name in params:
        param = getCaseParam(param_name, case_path)
        if is_int(param):
            d[param_name] = int(param)
        elif is_float(param):
            d[param_name] = float(param)
        else:
            d[param_name] = param
    return d

def get_first_line_and_names(file_path):
    line_number = 0
    with open(file_path, 'r') as f:
        for line_number, line in enumerate(f):
            if not line.strip().startswith('#'):
                break
            last_line = line
    # Last commented line contains them
    # wallHeatFlux specifies also dimensions as "min [W/m^2]"
    colnames = [s.strip() for s in last_line[1:].replace(' [', '[').split()] 
    return line_number, colnames

def read_fo_df(case_path, fo_name, file_name, time:str|list[str], ignore_index=False):
    if type(time) is str:
        time = [time]
    file_path = join(case_path, 'postProcessing', fo_name, time[0], file_name)
    line_number, colnames = get_first_line_and_names(file_path)
    kwargs = dict(skiprows=line_number, sep='\s+', index_col=0, names=colnames)
    df = pd.read_csv(file_path, **kwargs)
    for i in range(1,len(time)):
        file_path = join(case_path, 'postProcessing', fo_name, time[i], file_name)
        next_df = pd.read_csv(file_path, **kwargs)
        if not df.index[-1] <= next_df.index[0]:
            print(f"Warning: trimming {file_name} at time {time[i]}")
            next_df = next_df[next_df.index > df.index[-1]]
        df = pd.concat([df, next_df], ignore_index=ignore_index)
    return df

# TODO: make time arg variable:
# none - attempt to read all times
# str  - read one time
# list[str] - read multiple and join

# TODO: make consistent with __read
def read_multiple_df(pp_path, fo_name, op_name, times_list, read_csv_kwargs={}, ignore_index=False):
    def myread(time):
        path = join(pp_path,fo_name,str(time),op_name)
        with open(path) as f:
            for i, l in enumerate(f):
                if not l.startswith('#'):
                    break
        read_csv_kwargs['header'] = i - 1
        return pd.read_csv(path, **read_csv_kwargs)
    df = myread(times_list[0])
    for i in range(1,len(times_list)):
        # df = df.append(myread(times_list[i]))
        df = pd.concat([df, pd.DataFrame(myread(times_list[i]))], ignore_index=ignore_index)
    return df

def get_wpd_columns(csv_path):
    """Read webplotdigitizer column names"""
    with open(csv_path, 'r') as f:
        s1 = f.readline().strip().split(',')
        last = ''
        for i, el in enumerate(s1):
            if el == '':
                s1[i] = last
            last = el
        s2 = f.readline().strip().split(',')
    return ['_'.join(str(i) for i in s) for s in zip(s1, s2)]

def getCaseParam(param, casePath, caseParamsDict='caseParamsDict'):
    # getCaseParam('T', 'labReactor_act_H2_3g_1050C_dt0.000025_kN3')
    with open(join(casePath, caseParamsDict), 'r') as f:
        while line := f.readline():
            if line.strip().startswith(param + ' '):
                return line.strip().split()[1][:-1]
