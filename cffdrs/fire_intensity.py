def fire_intensity(fc, ros):
    """
    Fire Intensity Calculator

    Calculate the Predicted Fire Intensity

    All variables names are laid out in the same manner as Forestry Canada Fire
    Danger Group (FCFDG) (1992). Development and Structure of the Canadian Forest
    Fire Behavior Prediction System." Technical Report ST-X-3, Forestry Canada,
    Ottawa, Ontario.

    :param fc: Fuel Consumption (kg/m^2)
    :param ros: Rate of Spread (m/min)

    :returns: FI Fire Intensity (kW/m)
    """
    # Eq. 69 (FCFDG 1992) Fire Intensity (kW/m)
    FI = 300 * fc * ros
    return FI
