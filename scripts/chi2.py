"""
Chi squared test for 2 categories with Yates correction as desctibed here:
http://archive.bio.ed.ac.uk/jdeacon/statistics/tress9.html#Chi-squared test
"""

from scipy.stats import chisqprob
from math import fabs

A = 1427
B = 1031

expected = (A + B)/2.0  # the null hypothesis is that the numbers are the same

A_stat = ((fabs(A-expected)-0.5)**2)/expected
B_stat = ((fabs(B-expected)-0.5)**2)/expected

chi = A_stat + B_stat

probability = chisqprob(chi, 1)

print A_stat, B_stat, probability
