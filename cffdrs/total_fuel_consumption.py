from typing import Literal
from cffdrs.constants import FuelType


def crown_fuel_consumption(fuel_type: FuelType, cfl, cfb, pc, pdf):
    # Eq. 66a (Wotton 2009) - Crown Fuel Consumption (CFC)
    CFC = cfl * cfb
    if fuel_type in ["M1", "M2"]:
        # Eq. 66b (Wotton 2009) - CFC for M1/M2 types
        CFC = pc / 100 * CFC
    elif fuel_type in ["M3", "M4"]:
        # Eq. 66c (Wotton 2009) - CFC for M3/M4 types
        CFC = pdf / 100 * CFC
    return CFC


def total_fuel_consumption(
    fuel_type: FuelType, cfl, cfb, sfc, pc, pdf, option: Literal["TFC", "CFC"] = "TFC"
):
    """
    Total Fuel Consumption calculation

    Computes the Total (Surface + Crown) Fuel Consumption by Fuel
    Type.
    All variables names are laid out in the same manner as FCFDG (1992) or
    Wotton et. al (2009)

    Forestry Canada Fire Danger Group (FCFDG) (1992). "Development and
    Structure of the Canadian Forest Fire Behavior Prediction System."
    Technical Report ST-X-3, Forestry Canada, Ottawa, Ontario.

    Wotton, B.M., Alexander, M.E., Taylor, S.W. 2009. Updates and revisions to
    the 1992 Canadian forest fire behavior prediction system. Nat. Resour.
    Can., Can. For. Serv., Great Lakes For. Cent., Sault Ste. Marie, Ontario,
    Canada. Information Report GLC-X-10, 45p.

    :param fuel_type: The Fire Behaviour Prediction FuelType
    :param cfl: Crown Fuel Load (kg/m^2)
    :param cfb: Crown Fraction Burned (0-1)
    :param sfc: Surface Fuel Consumption (kg/m^2)
    :param pc: Percent Conifer (%)
    :param pdf: Percent Dead Balsam Fir (%)
    :param option: Type of output (TFC, CFC, default="TFC")

    :returns: TFC Total (Surface + Crown) Fuel Consumption (kg/m^2) OR CFC Crown Fuel Consumption (kg/m^2)
    """
    cfc = crown_fuel_consumption(fuel_type, cfl, cfb, pc, pdf)
    # Return CFC if requested
    if option == "CFC":
        return cfc
    # Eq. 67 (FCFDG 1992) - Total Fuel Consumption
    tfc = sfc + cfc
    return tfc
