from scipy.stats.distributions import binom
import sys

def pbinom(successes, fail, prob):
    total = successes + fail
    return binom.cdf(successes, total, prob)


if __name__ == "__main__":
    ref_read, alt_read = int(sys.argv[1]),  int(sys.argv[2])
    print pbinom(ref_read, alt_read, 0.5)




