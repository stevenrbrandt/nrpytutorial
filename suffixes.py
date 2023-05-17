import sympy as sp

subtable = {}

def setsuffix(f,t):
    assert type(f) == str or type(f) == unicode
    assert type(t) == str
    subtable[f] = t

def getsuffix(f):
    assert type(f) == str or type(f) == unicode
    return subtable.get(f,"")

def dosubs(expr):
    for sym in [s for s in expr.free_symbols]:
        ss = str(sym)
        if ss in subtable:
            expr = expr.subs(sym, sp.symbols(ss+subtable[ss]))
    return expr
