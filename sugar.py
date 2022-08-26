import grid
import sys
import inspect, re
import sympy as sp
import safewrite
import indexedexp as ixp
from here import here
from outputC import lhrh
from colored import colored
from here import here, herecc
from traceback import print_exc

import nrpylatex
from nrpylatex import parse_latex as parse_latex_

import warnings
warnings.filterwarnings('ignore',r'some variable\(s\) in the namespace were overridden')

coords = None

class excepthook:
    def __init__(self):
        self.excepthook = sys.excepthook
    def __enter__(self):
        pass
    def __exit__(self, ty, val, tb):
        sys.excepthook = self.excepthook

def parse_latex(*args,**kwargs):
    with excepthook():
        return parse_latex_(*args,**kwargs)

def set_coords(*c):
    global coords
    latex = "% coord ["+", ".join(c)+"]"
    parse_latex(latex)
    coords = c

###
class Seq:
    """
    Convert a sequence of D's and U's signifying "down" and "up" indexes
    into a latex equivalent which, when combined with specific index values
    becomes a format sequence for latex. Thus, "DD" becomes "_{%s %s}",
    "U" becomes "^{%s}" and "UUD" becomes "^{%s %s}_{%s}".
    """
    def __init__(self,seq):
        self.fmt = ""
        self.prev = ""
        for letter in seq:
            self.add(letter)
        if self.fmt != "":
            self.fmt += "}"
    def add(self, letter):
        if letter in ["l","D"]:
            if self.prev == letter:
                self.fmt += " %s"
            elif self.fmt == "":
                self.fmt += "_{%s"
            else:
                self.fmt += "}_{%s"
        elif letter in ["u","U"]:
            if self.prev == letter:
                self.fmt += " %s"
            elif self.fmt == "":
                self.fmt += "^{%s"
            else:
                self.fmt += "}^{%s"
        else:
            assert False
        self.prev = letter
            

def _latex_def(basename, fmt, gf, args):
    assert coords is not None, "Please call set_coords()"
    if type(gf) == list:
        for i in range(len(gf)):
            _latex_def(basename, fmt, gf[i], args + [coords[i]])
    else:
        indexes = fmt % tuple(args)
        latex = f"% \\text{{{basename}}}{indexes} = \\text{{{gf}}}"
        parse_latex(latex)


def latex_def(basename, sig, gf):
    s = Seq(sig)
    _latex_def(basename, s.fmt, gf, [])
###
def matchexpr(expr):
    result = []
    while True:
        g = re.match(r'\s+|%.*|\\text{([a-z]+)}|([a-z]+)|([\^_]){ *((?:[a-z] )*[a-z]) *}|([\^_])([a-z])|=',expr)
        if g:
            s = g.group(0)
            if g.group(1) is not None:
                result += [("symbol",g.group(1))]
            elif g.group(2) is not None:
                result += [("symbol",g.group(2))]
            elif g.group(3) is not None:
                if g.group(3) == "^":
                    typestr = "upindex"
                else:
                    typestr = "downindex"
                answer = re.split(r'\s+',g.group(4))
                result += [(typestr,answer)] 
            elif g.group(5) is not None:
                if g.group(5) == "^":
                    typestr = "upindex"
                else:
                    typestr = "downindex"
                answer = [g.group(6)]
                result += [(typestr,answer)] 
            elif g.group(0) == "=":
                result += [("equals","=")]
            expr = expr[g.end():]
        else:
            result += [("end",expr)]
            break
    return result
###

def getTypes():
    a,b = sp.symbols("a b")
    return type(a+b)

# Need to get type type of two symbols added together
# comparing with sp.core.add.Add seems to sometimes
# not work.
Add = getTypes()

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
    if len(args)>0 and type(args[-1]) == dict:
        g = args[-1]
        args = args[:-1]
    else:
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
                sym1 = copy.get("symmetry_option", "")
                for k in base_variants:
                    variant_copy = properties.get(k)
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
                latex = f"% define {fullname} --dim {properties.get('DIM',3)}"
                if sym1 not in ["", None]:
                    latex += f" --sym {sym1}"
                parse_latex(latex)

                base_variants.add(fullname)
                variants[basename] = base_variants

                properties[fullname] = copy

                gf = ixp.register_gridfunctions_for_single_rankN(**copy)

                namefun = copy.get("namefun",None)
                if namefun is not None:
                    latex_def(basename, suffix, gf)

                g[name] = gf
                definitions[name] = gf
            namelist = []
    assert len(namelist)==0, "Missing final index args"

indexdefs = {}

def gflatex(inp):
    globs = inspect.stack()[1].frame.f_globals
    assert "," not in inp, f"Commas are not valid in input: '{inp}'"
    args = matchexpr(inp)
    assert len(args) > 0, f"Failed to parse {inp}"
    assert args[-1][0] == "end", args[-1]
    end_len = len(args[-1][1])
    assert end_len == 0, f"Parse failure in input: '{inp[:end_len]}{colored('|','red')}{inp[end_len:]}"
    symbol = None
    indexes = []
    for a in args:
        if a[0] == "symbol":
            if symbol is not None:
                gfdecl(symbol, indexes, globs)
                symbol = None
                indexes = []
            symbol = a[1]
        elif a[0] == "upindex":
            for letter in a[1]:
                pair = indexdefs[letter]
                up = pair[1]
                indexes += [ up ]
        elif a[0] == "downindex":
            for letter in a[1]:
                pair = indexdefs[letter]
                down = pair[0]
                indexes += [ down ]
    if symbol is not None:
        gfdecl(symbol, indexes, globs)

def declIndexes():
    g = inspect.stack()[1].frame.f_globals
    for c in range(ord('a'),ord('z')+1):
        letter = chr(c)
        dn, up = sp.symbols(f"l{letter} u{letter}", cls=sp.Idx)
        g[f"l{letter}"] = dn
        g[f"u{letter}"] = up
        indexdefs[letter] = (dn, up)

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

UP_INDEX         = 1
DOWN_INDEX       = 2
CONTRACTED_INDEX = UP_INDEX | DOWN_INDEX

def getindexes(expr):
    """
    getindexes(expr) finds all the symbols in expr that
    represent up or down indexes.
    """
    indexes = {}
    for sym in expr.free_symbols:
        ssym = str(sym)
        g = matchindex(ssym)
        if g:
            updn = g.group(1)
            let = g.group(2)
            if updn == "u":
                mask = UP_INDEX
            else:
                mask = DOWN_INDEX
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
    if syms in [None,""]:
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

    if type(expr) == Add:
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
            subs[sym] = lookup(definitions[nm],indexes)
    return expr.subs(subs)

def geneqns2(lhs,rhs):
    lhs_str = r"\text{result}"
    last = ""
    suffix = ""
    # We expect a tensor expression of the form Foo[la,lb,ua,ub]
    assert type(lhs) == sp.Indexed, f"Type of lhs was '{type(lhs)}', not Indexed"
    for a in lhs.args[1:]:
        sa = str(a)
        # Ensure that we have an index, ua, ub,... or la, lb, ...
        assert len(sa) == 2
        if sa[0] == 'u':
            if last != "u":
                if last != "":
                    lhs_str += "}"
                lhs_str += "^{"
            suffix += "U"
        else:
            assert sa[0] == 'l'
            if last != "l":
                if last != "":
                    lhs_str += "}"
                lhs_str += "_{"
            suffix += "D"
        lhs_str += sa[1] +  " "
        last = sa[0]
    latex = lhs_str + "}=" + rhs
    r=parse_latex(latex,verbose=True)

    return geneqns(lhs=lhs, values=globals()["result"+suffix])

def geneqns(lhs,rhs=None,values=None,DIM=3):
    globs = inspect.stack()[1].frame.f_globals
    if values is None:
        assert rhs is not None, "Must supply either values or rhs to geneqns"
        lhs_indexes = getindexes(lhs)
        rhs_indexes = getindexes(rhs)
        for k in lhs_indexes:
            assert lhs_indexes[k] != CONTRACTED_INDEX, f"Contracted indexes are not allowed on the left hand side: '{lhs}'"
            assert lhs_indexes.get(k,-1) == rhs_indexes.get(k,-1), f"Free index '{k}' does not match on the lhs and rhs."
        for k in rhs_indexes:
            if rhs_indexes[k] == CONTRACTED_INDEX:
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
            assert indexes[index] != CONTRACTED_INDEX, f"Error, contracted index in lhs: '{index}'"
        nm = str(lhs.base)+getsuffix(lhs)
        props = properties[nm]
        symmetries = props.get("symmetry_option","")
        for index in incrindexes([0] * len(indexes), DIM, getsyms(symmetries)):
            result += [lhrh(lhs=lookup(definitions[nm],index),rhs=lookup(values, index))]
        for r in result:
            pass #here(r)
        return result #generatevalues(lhs,rhs,[0]*len(indexes),symmetries)
    else:
        assert False, "Must supply either values or rhs to geneqns"
