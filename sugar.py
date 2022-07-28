import inspect, re
import sympy as sp
import safewrite
import indexedexp as ixp
from here import here
from outputC import lhrh
from colored import colored
from here import here, herecc

properties = {}
variants = {}
definitions = {}

verbose = False

sim_params = {}

def gfparams(**kw):
    global sim_params
    sim_params = {}
    for k in kw:
        if k == "symmetries":
            sim_params["symmetry_option"] = kw[k]
        else:
            sim_params[k] = kw[k]
    if "DIM" not in sim_params:
        sim_params["DIM"] = 3

def gfdecl(*args):
    g = inspect.stack()[1].frame.f_globals
    namelist = []
    for arg in args:
        if type(arg) == str:
            namelist += [arg]
            #sim_params["gf_basename"] = name
            #g[name] = ixp.register_gridfunctions_for_single_rankN(**sim_params)
        elif type(arg) == list:
            #gfparams(rank=len(arg))
            rank = len(arg)
            suffix = ""
            for k in arg:
                assert type(k) == sp.tensor.indexed.Idx
                if str(k)[0] == "u":
                    suffix += "U"
                elif str(k)[0] == "l":
                    suffix += "D"
                else:
                    assert False,f"{k} {type(k)}"
            for basename in namelist:
                assert not re.match(r'^.*[DU]$',basename), f"Bad declaration for '{basename}'. Basenames should not end in D or U."
                fullname = basename + suffix
                if basename not in g:
                    if rank>0:
                        g[basename] = sp.IndexedBase(basename,shape=tuple([sim_params["DIM"]]*len((arg))))
                name = basename + suffix
                copy = {}
                for k in sim_params:
                    copy[k] = sim_params[k]
                copy["gf_basename"] = name
                copy["rank"] = rank
                if rank < 2:
                    copy["symmetry_option"] = None
                if copy["gf_type"] != "EXTERNAL":
                    copy["external_module"] = None

                assert fullname not in properties, f"Redefinition of {fullname}"

                base_variants = variants.get(basename,set())
                for k in base_variants:
                    variant_copy = properties.get(k)
                    sym1 = copy.get("symmetry_option", "")
                    sym2 = variant_copy.get("symmetry_option", "")
                    assert sym1 == sym2, \
                        f"Inconsistent declaration of {basename}. Variant {k} has {sym2}, and {fullname} has {sym1}"

                if verbose:
                    print(colored("Adding Definition for:","cyan"),basename)
                    for k in copy:
                        print("  ",colored(k+":","yellow"),copy[k])
                    if len(base_variants) > 0:
                        print("  ",colored("Previous Definitions:","yellow"),base_variants)
                    print()

                base_variants.add(fullname)
                variants[basename] = base_variants

                properties[fullname] = copy

                gf = ixp.register_gridfunctions_for_single_rankN(**copy)
                g[name] = gf
                definitions[name] = gf
            namelist = []
    assert len(namelist)==0, "Missing final index args"

def declIndexes():
    g = inspect.stack()[1].frame.f_globals
    for c in range(ord('a'),ord('z')+1):
        letter = chr(c)
        dn, up = sp.symbols(f"l{letter} u{letter}", cls=sp.Idx)
        g[f"l{letter}"] = dn
        g[f"u{letter}"] = up

def ixnam(i):
    return ["x", "y", "z"][i]

def namefun(symbol, index, shape, prefix):
    symbol = prefix
    result = [sp.Symbol(symbol + ''.join(ixnam(n) for n in index + [i]))
            if symbol else sp.sympify(0) for i in range(shape[0])]
    return result

def name_xyz(sym,ind,shape):
    symbase = re.sub("[UD]+$","",sym)
    return namefun(sym,ind,shape,symbase)

def matchindex(s):
    return re.match(r'^([ul])([a-z])$', str(s))

def getindexes(expr):
    indexes = {}
    for sym in expr.free_symbols:
        ssym = str(sym)
        g = matchindex(ssym)
        if g:
            updn = g.group(1)
            let = g.group(2)
            if updn == "u":
                mask = 1
            else:
                mask = 2
            indexes[let] = indexes.get(let,0) | mask
    return indexes

def incrindexes(indexes_input,dim,symmetries):
    indexes = [x for x in indexes_input]
    yield indexes
    while True:
        for i in range(len(indexes)):
            indexes[i] += 1
            max_val = dim
            for symmetry in symmetries:
                if symmetry[0] == i:
                    j = symmetry[1]
                    max_val = min(max_val,indexes[j]+1)
            if indexes[i] >= max_val:
                if i+1 == len(indexes):
                    return
                indexes[i] = 0
            else:
                result = [x for x in indexes[::-1]]
                yield result
                break

def lookup(array, indexes, i=0):
    if i >= len(indexes):
        return array
    else:
        return lookup(array[indexes[i]], indexes, i+1)
    
def getsyms(syms):
    li = []
    if syms == "":
        return li
    for sym in syms.split('_'):
        g = re.match(r'^(a?sym)(\d)(\d)', sym)
        assert g, f"Bad symmetry: '{sym}'"
        if g.group(1) == "sym":
            sign = 1
        else:
            sign = -1
        li += [(int(g.group(2)),int(g.group(3)),sign)]
    return li

def getsuffix(expr):
    suffix = ""
    for sym in expr.free_symbols:
        g = matchindex(sym) #re.match(r'^([ul])([a-z])$', str(sym))
        if g:
            if g.group(1) == "u":
                suffix += "U"
            else:
                suffix += "D"
    return suffix

def getname(expr):
    assert type(sym) == sp.tensor.indexed.Indexed
    nm = str(sym.base)
    indexes = []
    for k in sym.args[1:]:
        ks = str(k)
        g = re.match(r'([ul]).$', ks)
        if g.group(1) == "u":
            nm += "U"
        else:
            nm += "D"
    return nm

def makesum(expr, dim=3):
    #globs = inspect.stack()[2].frame.f_globals

    if type(expr) == sp.core.add.Add:
        sume = sp.sympify(0)
        for a in expr.args:
            sume += makesum(expr,dim)
        return sume
    elif expr.is_Function:
        assert len(expr.args)==1
        expr = expr.func(makesum(expr.args[0]))
        return expr

    indexes = {}
    for f in expr.free_symbols:
        fs = str(f)
        g = matchindex(fs)
        if g:
            # This is an index
            updown = g.group(1)
            letter = g.group(2)
            if letter not in indexes:
                indexes[letter] = 0
            if updown == "u":
                indexes[letter] |= 1
            else:
                indexes[letter] |= 2
    count = 0
    for index in indexes:
        if indexes[index] == 3:
            new_expr = sp.sympify(0)
            for d in range(dim):
                un, dn = sp.symbols(f"u{index} l{index}", cls=sp.Idx)
                u1, d1 = sp.symbols(f"u{d} l{d}", cls=sp.Idx)
                new_expr += expr.subs(un,u1).subs(dn,d1)
            expr = new_expr
    subs = {}
    for sym in expr.free_symbols:
        if type(sym) == sp.tensor.indexed.Indexed:
            nm = str(sym.base)
            indexes = []
            for k in sym.args[1:]:
                ks = str(k)
                g = re.match(r'([ul])(\d)$', ks)
                if g.group(1) == "u":
                    nm += "U"
                else:
                    nm += "D"
                indexes += [int(g.group(2))]
            #here("sym:",sym,"->",lookup(globs[nm],indexes))
            subs[sym] = lookup(definitions[nm],indexes)
    return expr.subs(subs)

def geneqns(lhs,rhs=None,values=None,DIM=3):
    globs = inspect.stack()[1].frame.f_globals
    if values is None:
        assert rhs is not None, "Must supply either values or rhs to geneqns"
        lhs_indexes = getindexes(lhs)
        rhs_indexes = getindexes(rhs)
        for k in lhs_indexes:
            assert lhs_indexes[k] != 3, f"Contracted indexes are not allowed on the left hand side: '{lhs}'"
            assert lhs_indexes.get(k,-1) == rhs_indexes.get(k,-1), f"Free index '{k}' does not match on the lhs and rhs."
        for k in rhs_indexes:
            if rhs_indexes[k] == 3:
                continue
            assert lhs_indexes.get(k,-1) == rhs_indexes.get(k,-1), f"Free index '{k}' does not match on the rhs and lhs."
        if len(lhs_indexes)==0:
            return [lhrh(lhs=lhs, rhs=makesum(rhs))]
        else:
            assert False
    elif rhs is None:
        result = []
        assert values is not None, "Must supply either values or rhs to geneqns"
        indexes = getindexes(lhs)
        for index in indexes:
            assert indexes[index] != 3, f"Error, contracted index in lhs: '{index}'"
        nm = str(lhs.base)+getsuffix(lhs)
        props = properties[nm]
        symmetries = props.get("symmetry_option","")
        for index in incrindexes([0] * len(indexes),3, getsyms(symmetries)):
            result += [lhrh(lhs=lookup(globs[nm],index),rhs=lookup(values, index))]
        return result #generatevalues(lhs,rhs,[0]*len(indexes),symmetries)
    else:
        assert False, "Must supply either values or rhs to geneqns"
