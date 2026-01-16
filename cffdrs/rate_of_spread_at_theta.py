import math

from cffdrs.r_helpers import safe_div


def rate_of_spread_at_theta(ros, fros, bros, theta):
    """
    Rate of spread at a point along the perimeter calculator

    Computes the Rate of Spread at any point along the perimeter of
    an elliptically shaped fire. Equations are from Wotton et. al. (2009).

    Wotton, B.M., Alexander, M.E., Taylor, S.W. 2009. Updates and revisions to
    the 1992 Canadian forest fire behavior prediction system. Nat. Resour.
    Can., Can. For. Serv., Great Lakes For. Cent., Sault Ste. Marie, Ontario,
    Canada. Information Report GLC-X-10, 45p.

    :param ros: Rate of Spread (m/min)
    :param fros: Flank Fire Rate of Spread (m/min)
    :param bros: Back Fire Rate of Spread (m/min)
    :param theta: THETA

    :returns: ROSTHETA - Rate of spread at point theta(m/min)
    """
    c1 = math.cos(theta)
    s1 = math.sin(theta)
    c1 = math.cos(theta + 0.001) if c1 == 0 else c1
    # Eq. 94 - Calculate the Rate of Spread at point THETA
    # large equation, view the paper to see a better representation
    ROStheta = (
        ((ros - bros) / (2 * c1) + (ros + bros) / (2 * c1))
        * (
            safe_div((fros * c1 * math.sqrt(fros * fros * c1 * c1 + (ros * bros) * s1 * s1)
              - ((ros * ros - bros * bros) / 4) * s1 * s1),
            (fros * fros * c1 * c1
                + ((ros + bros) / 2) * ((ros + bros) / 2) * s1 * s1))
          )
    )
    return ROStheta