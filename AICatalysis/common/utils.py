import itertools


def sort_defaultdict(ddict):
    return {key: value for key, value in sorted(ddict.items(), key=lambda x: len(x[1]), reverse=True)}


def get_combinations(sequence1, sequence2):
    return set(tuple([tuple(sorted(item)) for item in itertools.product(sequence1, sequence2) if item[0] != item[1]]))
