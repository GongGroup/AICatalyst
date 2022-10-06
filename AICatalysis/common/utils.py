def sort_defaultdict(ddict):
    return {key: value for key, value in sorted(ddict.items(), key=lambda x: len(x[1]), reverse=True)}
