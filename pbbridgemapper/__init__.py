VERSION = (0, 4)

_changelist = "$Change: 138434 $"

def get_version():
    """Return the version as a string. "O.7"

    This uses a major.minor
    """
    return ".".join([str(i) for i in VERSION])

def get_full_version():
    """Return major.minor.changelist."""

    version = get_version()
    cl = str(_get_changelist(_changelist))
    return '.'.join([version, cl])

def _get_changelist(perforce_str):
    import re
    rx = re.compile(r'Change: (\d+)')
    match = rx.search(perforce_str)
    if match is None:
        v = 'UnknownChangelist'
    else:
        try:
            v = int(match.group(1))
        except (TypeError, IndexError):
            v = "UnknownChangelist"
    return v
