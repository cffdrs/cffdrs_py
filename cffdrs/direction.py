def direction(bearingT1T2, bearingT1T3, ThetaAdeg):
    """
    Direction definer

    :param bearingT1T2: Bearing between T1 and T2
    :param bearingT1T3: Bearing between T1 and T3
    :param ThetaAdeg: Direction

    :returns: dir - a direction in degrees

    Note: This implementation normalizes the result using 360-degree adjustments,
    fixing what I believe is a bug in the original R version which used 10-degree adjustments.
    The result is that every direction returned will be between +180° and -180°, which wasn't always the
    case with the R implementation.
    """
    dir = None
    if bearingT1T2 > 0 and bearingT1T3 > 0 and bearingT1T2 > bearingT1T3:
        dir = bearingT1T2 - ThetaAdeg
    elif bearingT1T2 > 0 and bearingT1T3 > 0 and bearingT1T2 < bearingT1T3:
        dir = bearingT1T2 + ThetaAdeg
    elif bearingT1T2 < 0 and bearingT1T3 < 0 and bearingT1T2 > bearingT1T3:
        dir = bearingT1T2 - ThetaAdeg
    elif bearingT1T2 < 0 and bearingT1T3 < 0 and bearingT1T2 < bearingT1T3:
        dir = bearingT1T2 + ThetaAdeg
    elif bearingT1T2 > 0 and bearingT1T2 < 90 and bearingT1T3 < 0 and bearingT1T3 > -90:
        dir = bearingT1T2 - ThetaAdeg
    elif bearingT1T2 < 0 and bearingT1T2 > -90 and bearingT1T3 > 0 and bearingT1T3 < 90:
        dir = bearingT1T2 + ThetaAdeg
    elif bearingT1T2 > 90 and bearingT1T3 < -90 and bearingT1T2 + ThetaAdeg > 180:
        dir = bearingT1T2 + ThetaAdeg - 360
    elif bearingT1T2 > 90 and bearingT1T3 < -90 and bearingT1T2 + ThetaAdeg < 180:
        dir = bearingT1T2 + ThetaAdeg
    elif bearingT1T2 < -90 and bearingT1T3 > 90 and bearingT1T2 - ThetaAdeg < -180:
        dir = bearingT1T2 - ThetaAdeg + 360
    elif bearingT1T2 < -90 and bearingT1T3 > 90 and bearingT1T2 - ThetaAdeg > -180:
        dir = bearingT1T2 - ThetaAdeg
    if dir is not None:
        if dir < -180:
            dir += 360
        if dir > 180:
            dir -= 360
    return dir
