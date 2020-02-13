import numpy as np
from scipy.stats import ortho_group


def get_four_norm_maximization_iterates(size):
    """
    Runs the matching, stretching and projection (MSP) algorithm to iteratively
    maximize the (entry-wise) four norm of a SIZE by SIZE matrix, starting at a
    random orthogonal matrix, drawn from the Haar distribution. Returns an
    infinite generator which yields the current matrix at each iteration,
    beginning with random starting point.

    :param size: (int) The size of the matrix.
    :return: (generator) An infinite generator which yields all intermediate
        iterates of the algorithm.
    """
    current = get_random_orthogonal(size)
    yield current
    while True:
        delta_current = current * current * current
        left, _, right = np.linalg.svd(delta_current, compute_uv=True)
        projection = left @ right
        current = projection
        yield current


def get_dictionary_learning_iterates(observations):
    """
    Runs the matching, stretching and projection (MSP) algorithm to iteratively
    learn an orthogonal dictionary given a matrix of observations. Similarly
    starts at a random orthogonal matrix and iteratively improves the
    dictionary. Returns an infinite generator which yields the current
    dictionary at each iteration. Note that for practical purposes, this should
    be truncated.

    :param observations: (numpy.ndarray) A numpy array, where each column
        represents a new observation.
    :return: (generator) An infinite generator which yields all intermediate
        iterates of the algorithm.
    """
    current = get_random_orthogonal(len(observations))
    yield current
    while True:
        matched = current @ observations
        delta_current = (matched * matched * matched) @ observations.T
        left, _, right = np.linalg.svd(delta_current, compute_uv=True)
        projection = left @ right
        current = projection
        yield current


def random_dictionary_learning_instance(features, samples, theta):
    """
    :param features: (int) The number of features for each sample, equivalently
        the length of the signal.
    :param samples: (int) The number of samples or signals.
    :param theta: (float) The probability a particular entry in the decoded
        samples is non-zero. The smaller THETA is, the sparser the signals
        are in the optimal, intended, basis.
    :return: (tuple)
    """
    dictionary = get_random_orthogonal(features)
    samples = get_bernoulli_gaussian(theta, (features, samples))
    observations = dictionary @ samples
    return observations, dictionary, samples


def sum_of_fourth_powers(matrix):
    """
    :param matrix: (numpy.ndarray) A numpy array.
    :return: The fourth power of the four-norm of the matrix. In other words,
        the sum of the fourth power of all of its entries.
    """
    squared_entries = matrix * matrix
    return np.sum(squared_entries * squared_entries)


def get_random_orthogonal(size):
    """
    :param size: (int) The dimension of the matrix.
    :return: (numpy.ndarray) Returns a random orthogonal matrix from O(SIZE),
        the orthogonal group of dimension SIZE, drawn from the Haar
        distribution. The matrix has size (SIZE, SIZE).
    """
    return ortho_group.rvs(size)


def get_bernoulli_gaussian(theta, size):
    """
    :param theta: (float) The probability a particular entry is non-zero: must
        be between 0 and 1 inclusive. The smaller THETA is, the more sparse the
        output will be in expectation.
    :param size: (int or tuple) The shape of the output.
    :return: (numpy.ndarray) A random numpy array where each entry is from
        independently and identically distributed according  to a bernoulli-
        gaussian
    """
    bernoulli = np.random.binomial(1, theta, size)
    gaussian = np.random.standard_normal(size)
    result = bernoulli * gaussian
    return result
