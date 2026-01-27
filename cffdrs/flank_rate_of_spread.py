def flank_rate_of_spread(ros, bros, lb):
    """
    Flank Fire Rate of Spread Calculator

    Calculate the Flank Fire Spread Rate.

    All variables names are laid out in the same manner as Forestry Canada
    Fire Danger Group (FCFDG) (1992). Development and Structure of the
    Canadian Forest Fire Behavior Prediction System." Technical Report
    ST-X-3, Forestry Canada, Ottawa, Ontario.

    :param ros: Fire Rate of Spread (m/min)
    :param bros: Back Fire Rate of Spread (m/min)
    :param lb: Length to breadth ratio

    :returns: FROS Flank Fire Spread Rate (m/min) value
    """
    # Eq. 89 (FCFDG 1992)
    FROS = (ros + bros) / lb / 2
    return FROS
