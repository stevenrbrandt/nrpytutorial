def here(*args):
    """
    Here is a debugging tool. It works like print(), except it prefixes any
    output with "HERE file:line ".
    """
    import inspect
    stack = inspect.stack()
    frame = stack[1]
    print("HERE:","%s:%d" % (frame.filename, frame.lineno), *args, flush=True)
    frame = None
    stack = None
