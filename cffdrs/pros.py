from dataclasses import dataclass
import math
from cffdrs.direction import direction


@dataclass
class ProsOutput:
    ros: float
    direction: float


def haversine_distance(lon1, lat1, lon2, lat2):
    "R uses the geosphere library, this is a manual implementation of haversine distance"
    R = 6371000  # radius of Earth in meters
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (
        math.sin(dlat / 2) ** 2
        + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2) ** 2
    )
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c


def bearing(lon1, lat1, lon2, lat2):
    "R uses the geosphere library, this is a manual implementation of bearing"
    dlon = math.radians(lon2 - lon1)
    lat1_rad = math.radians(lat1)
    lat2_rad = math.radians(lat2)
    x = math.sin(dlon) * math.cos(lat2_rad)
    y = math.cos(lat1_rad) * math.sin(lat2_rad) - math.sin(lat1_rad) * math.cos(
        lat2_rad
    ) * math.cos(dlon)
    bearing_rad = math.atan2(x, y)
    return math.degrees(bearing_rad)


def pros(T1, Long1, Lat1, T2, Long2, Lat2, T3, Long3, Lat3):
    """
    Point-based input for Simard Rate of Spread and Direction

    pros is used to calculate the rate of spread and direction given one set of
    three point-based observations of fire arrival time. The function requires
    that the user specify the time that the fire crossed each point, along with
    the latitude and longitude of each observational point. The function allows
    quick input for a single triangle.

    pros allows users to calculate the rate of spread and direction of a fire
    across a triangle, given three time measurements and details about the
    orientation and distance between observational points. The algorithm is
    based on the description from Simard et al. (1984).

    Rate of spread and direction of spread are primary variables of interest
    when observing wildfire growth over time. Observations might be recorded
    during normal fire management operations (e.g., by a Fire Behaviour
    Analyst), during prescribed fire treatments, and during experimental
    research burns. Rate of spread is especially important for estimating
    Byram's fireline intensity, fireline intensity = heat constant of fuel ×
    weight of fuel consumed × forward rate of spread (Byram 1959).

    Rate of spread is difficult to measure and highly variable in the field.
    Many techniques were proposed over the years, but most were based on
    observations collected from a pre-placed reference grid and stopwatch
    (Curry and Fons 1938; Simard et al. 1982). Early approaches required that
    observers be in visual contact with the reference grid, but later,
    thermocouples and dataloggers were employed to measure the onset of the
    heat pulse at each point.

    Simard et al. (1982) proposed calculations for spread based on an
    equilateral triangle layout. Simard et al. (1984) proposed calculations for
    spread based on any type of triangle. Both articles also discussed field
    sampling design and layout, with special attention to the size of the
    triangles (large enough that the fire traverses the triangle in one to two
    minutes) and even using triangles of varying size within one field plot
    (but no triangle larger than one fourth of the site's total area).

    The underlying algorithms use trigonometry to solve for rate of spread and
    direction of spread. One important assumption is that the spread rate and
    direction is uniform across one triangular plot, and that the fire front
    is spreading as a straight line; Simard et al. (1982, 1984) acknowledge
    that these assumption are likely broken to some degree during fire spread
    events.

    :author: Tom Schiks, Xianli Wang, Alan Cantin

    :references:
        1. Simard, A.J., Eenigenburg, J.E., Adams, K.B., Nissen, R.L., Deacon,
           and Deacon, A.G. 1984. A general procedure for sampling and
           analyzing wildland fire spread.
        2. Byram, G.M. 1959. Combustion of forest fuels. In: Davis, K.P. Forest
           Fire Control and Use. McGraw-Hill, New York.
        3. Curry, J.R., and Fons, W.L. 1938. Rate of spread of surface fires in
           the Ponderosa Pine Type of California. Journal of Agricultural
           Research 57(4): 239-267.
        4. Simard, A.J., Deacon, A.G., and Adams, K.B. 1982. Nondirectional
           sampling wildland fire spread. Fire Technology: 221-228.

    :param T1: (required) Time that the fire front crossed point 1. Time
        entered in fractional format. Output ROS will depend on the level of
        precision entered (minute, second, decisecond)
    :param Long1: (required) Longitude for datalogger 1. (decimal degrees)
    :param Lat1: (required) Latitude for datalogger 1. (decimal degrees)
    :param T2: (required) Time that the fire front crossed point 2. Time
        entered in fractional format. Output ROS will depend on the level of
        precision entered (minute, second, decisecond)
    :param Long2: (required) Longitude for datalogger 2. (decimal degrees)
    :param Lat2: (required) Latitude for datalogger 2. (decimal degrees)
    :param T3: (required) Time that the fire front crossed point 3. Time
        entered in fractional format. Output ROS will depend on the level of
        precision entered (minute, second, decisecond)
    :param Long3: (required) Longitude for datalogger 3. (decimal degrees)
    :param Lat3: (required) Latitude for datalogger 3. (decimal degrees)

    :returns: Dict with Ros and Direction. Output units depend on the user's
        inputs for distance (typically meters) and time (seconds or minutes).

    """
    # Compute lengths
    length_t1t2 = haversine_distance(Long1, Lat1, Long2, Lat2)
    length_t1t3 = haversine_distance(Long1, Lat1, Long3, Lat3)
    length_t2t3 = haversine_distance(Long2, Lat2, Long3, Lat3)

    # Compute bearings
    bearing_t1t2 = bearing(Long1, Lat1, Long2, Lat2)
    bearing_t1t3 = bearing(Long1, Lat1, Long3, Lat3)

    angle_arad = math.acos(
        (length_t1t3**2 + length_t1t2**2 - length_t2t3**2) / (2 * length_t1t3 * length_t1t2)
    )
    theta_arad = math.atan(
        (T3 - T1) / (T2 - T1) * (length_t1t2 / (length_t1t3 * math.sin(angle_arad)))
        - (1 / math.tan(angle_arad))
    )
    theta_adeg = (theta_arad * 180) / math.pi
    dir = direction(bearing_t1t2, bearing_t1t3, theta_adeg)
    ros = (length_t1t2 * math.cos(theta_arad)) / (T2 - T1)
    return ProsOutput(ros=ros, direction=dir)
