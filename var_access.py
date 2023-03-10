
from_access = {}

def var_from_access(access):
    access = access.strip() # This really shouldn't be necessary
    v = from_access.get(access, None)
    if v is not None:
        return v
    g = re.match(r'^(\w+)(\[\w+\])*$', access)
    if g:
        return g.group(1)
    g = re.match(r'in_gfs\w*\[IDX4S\((\w+),i0,i1,i2\)\]', access)
    if g:
        return g.group(1)
    g = re.match(r'^const\s+(\w+)\s+(\w+)', access)
    if g:
        return g.group(2)
    raise Exception(f("Could not identify a variable name from the access string '{access}'"))

def gfaccess(gfarrayname = "", varname = "", ijklstring = "", context = "DECL"):
    ret = _gfaccess(gfarrayname, varname, ijklstring, context)
    from_access[ret] = varname
    return ret
