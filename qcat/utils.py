try:
    from string import maketrans
except ImportError:
    maketrans = str.maketrans

LOOKUP = []

for q in range(100):
    LOOKUP.append(pow(10,-.1*q))


def qstring_to_phred(quality):
    """ Compute standard phred scores from a quality string. """
    qscores = []
    if quality is not None:
        qscores = [ord(q) - 33 for q in quality]
    return qscores


def revcomp(seq):
    return seq.translate(maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]


def mean_error_prob(scores):
    """ Returns the phred score corresponding to the mean of the probabilities
    associated with the phred scores provided.

    :param scores: Iterable of phred scores.

    :returns: Phred score corresponding to the average error rate, as
        estimated from the input phred scores.
    """
    if not scores:
        return -1.0
    sum_prob = 0.0
    for val in scores:
        sum_prob += LOOKUP[val]
    mean_prob = sum_prob / len(scores)
    # return -10.0 * log10(mean_prob)
    return mean_prob


