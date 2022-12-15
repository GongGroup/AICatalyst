import itertools
import re
import tokenize
from collections.abc import Iterable
from io import BytesIO


def sort_defaultdict(ddict):
    return {key: value for key, value in sorted(ddict.items(), key=lambda x: len(x[1]), reverse=True)}


def get_combinations(sequence1, sequence2):
    return set(tuple([tuple(sorted(item)) for item in itertools.product(sequence1, sequence2) if item[0] != item[1]]))


def get_tokens(lines: list):
    tokens = []
    for line in lines:
        _tokens = tokenize.tokenize(BytesIO(line.encode('utf-8')).readline)
        tokens.append([token for token in _tokens])
    return tokens


def flatten(items, ignore_types=(str, bytes)):
    for x in items:
        if isinstance(x, Iterable) and not isinstance(x, ignore_types):
            yield from flatten(x)
        else:
            yield x


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False


def is_ratio(s):
    if re.search(r"[0-9]:[0-9]", s) is not None:
        return True
    else:
        return False
