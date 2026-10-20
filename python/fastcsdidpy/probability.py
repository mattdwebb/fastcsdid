"""MIT 2026 Matt Webb; vendored from didsuitepy's independent math-only tails."""
import math


def normal_two_sided(statistic):
    return math.erfc(abs(statistic)/math.sqrt(2))


def _fraction(a, b, x):
    tiny, tolerance = 1e-300, 3e-14
    c = 1.0
    d = 1.0 - (a + b) * x / (a + 1.0)
    d = 1.0 / (d if abs(d) > tiny else tiny)
    result = d
    for m in range(1, 10001):
        for numerator in (m * (b - m) * x / ((a + 2*m - 1) * (a + 2*m)),
                          -(a + m) * (a + b + m) * x / ((a + 2*m) * (a + 2*m + 1))):
            d = 1 + numerator * d
            c = 1 + numerator / c
            if abs(d) < tiny:
                d = tiny
            if abs(c) < tiny:
                c = tiny
            d = 1 / d
            delta = c * d
            result *= delta
        if abs(delta - 1) < tolerance:
            return result
    raise ArithmeticError("incomplete beta continued fraction did not converge")


def student_two_sided(statistic, df):
    if df <= 0 or not math.isfinite(df) or math.isnan(statistic):
        raise ValueError("finite positive df and non-NaN statistic required")
    if math.isinf(statistic):
        return 0.0
    if statistic == 0:
        return 1.0
    x = df / (df + statistic * statistic)
    if x == 0:
        return 0.0
    a, b = df / 2, 0.5
    front = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                     + a * math.log(x) + b * math.log1p(-x)) if x < 1 else 0.0
    if x == 1:
        return 1.0
    if x < (a + 1) / (a + b + 2):
        p = front * _fraction(a, b, x) / a
    else:
        p = 1 - front * _fraction(b, a, 1-x) / b
    return min(1.0, max(0.0, p))
