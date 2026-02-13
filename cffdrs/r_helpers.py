import math


def safe_div(numerator, denominator):
    """
    Divide two numbers in a way that mimics R's numeric behavior.

    Returns Inf or NaN instead of raising exceptions for division by zero:
      - If denominator is 0 and numerator is 0, returns NaN.
      - If denominator is 0 and numerator is non-zero, returns Inf with the correct sign.
      - Otherwise returns numerator / denominator normally.

    Useful for porting R code where divisions by zero occur but NaN/Inf are valid results.

    :param numerator: The numerator of the division.
    :param denominator: The denominator of the division.
    :return: The result of the division, or Inf/NaN if the denominator is zero.
    """
    if denominator == 0:
        if numerator == 0:
            return math.nan
        return math.copysign(math.inf, numerator)
    return numerator / denominator
