import sys

def openOrFail(name, mode, msg = 'cannot open "%s".\n'):
    try:
        return open(name, mode)
    except:
        try:
            sys.stderr.write(msg % name)
        except:
            sys.stderr.write(msg)
        sys.exit(1)
        
