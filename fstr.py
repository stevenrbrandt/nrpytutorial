import re
import sys

def f(s):
    """
    Mimic the functionality of formatted strings in Python3. Convert curly brackets in s
    to expressions.
    >>> f('3+2={3+2}.')
    '3+2=5.'
    >>> f('3+2={3+2}')
    '3+2=5'
    >>> f('{"="*3} test {"="*3}')
    '=== test ==='
    >>> f('{{hello}}')
    '{hello}'
    """
    count = 0
    ns = ''
    w = ''
    for g in re.finditer(r'\\(.)|.', s):
        if g.group(1) is not None:
            if count > 0:
                w += g.group(0)
            else:
                ns += g.group(0)
        elif g.group(0) == '{':
            if count > 0:
                w += '{' 
            count += 1
        elif g.group(0) == '}':
            assert count > 0
            if count > 1:
                w += "}"
            count -= 1
            if w == '':
                pass
            elif w.startswith("{"):
                ns += w
            else:
                ns += str(eval(w))
            w = ''
        elif count > 0:
            w += g.group(0)
        else:
            ns += g.group(0)
    assert w == '', f("Unclosed curly bracket in expression: '{s}'")
    return ns

if __name__ == "__main__":
    import doctest
    doctest.testmod()
