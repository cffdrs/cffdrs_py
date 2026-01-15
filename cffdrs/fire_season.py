import datetime


def fire_season(input, fs_start=12, fs_end=5, method="WF93", consistent_snow=False, multi_year=False):
    """
    Fire Season Start and End

    Calculates the start and end fire season dates for given weather data.
    The current method used is based on three consecutive daily maximum temperature thresholds.

    :param input: List of dicts with keys yr, mon, day, tmax, snow_depth (optional)
    :param fs_start: Temperature threshold to start the fire season (default=12)
    :param fs_end: Temperature threshold to end the fire season (default=5)
    :param method: Method of fire season calculation ("WF93" or "LA08", default="WF93")
    :param consistent_snow: Is consistent snow data in the input? (default=False)
    :param multi_year: Should the fire season span multiple years? (default=False)

    :returns: List of dicts with yr, mon, day, fsdatetype, date
    """
    # Normalize input keys to lowercase
    normalized_input = []
    for row in input:
        normalized_input.append({k.lower(): v for k, v in row.items()})

    # Extract variables
    yr = [row.get('yr') for row in normalized_input]
    mon = [row.get('mon') for row in normalized_input]
    day = [row.get('day') for row in normalized_input]
    tmax = [row.get('tmax') for row in normalized_input]
    snow_depth = [row.get('snow_depth', 0) for row in normalized_input]

    # Validate required fields
    if not all(yr):
        raise ValueError("Year (yr) is required for this function.")
    if not all(mon):
        raise ValueError("Month (mon) is required for this function.")
    if not all(x in range(1, 13) for x in mon if x is not None):
        raise ValueError("mon value is out of bounds (1-12).")
    if not all(day):
        raise ValueError("Day (day) is required for this function.")
    if not all(x in range(1, 32) for x in day if x is not None):
        raise ValueError("day value is out of bounds (1-31).")
    if not all(tmax):
        raise ValueError("Maximum Daily Temperature (tmax) is required for this function.")

    method = method.upper()
    if method not in ["WF93", "LA08"]:
        raise ValueError(f"Selected method '{method}' is unavailable.")

    n0 = len(tmax)
    season_active = False
    season_start_end = []

    if method == "WF93":
        for k in range(n0):
            if k > 2:
                # Check for start
                if not season_active and all(t > fs_start for t in tmax[k-3:k]):
                    season_active = True
                    the_day = day[k]
                    if not multi_year and mon[k] == 1 and day[k] == 4:
                        the_day = day[k-3]
                    season_start_end.append({
                        'yr': yr[k],
                        'mon': mon[k],
                        'day': the_day,
                        'fsdatetype': 'start'
                    })
                # Check for end
                if season_active and all(t < fs_end for t in tmax[k-3:k]):
                    season_active = False
                    season_start_end.append({
                        'yr': yr[k],
                        'mon': mon[k],
                        'day': day[k],
                        'fsdatetype': 'end'
                    })
    elif method == "LA08":
        if consistent_snow:
            if not any('snow_depth' in row for row in normalized_input):
                raise ValueError("Snow depth is required for the selected method 'LA08'.")
            for k in range(n0):
                if k > 2:
                    # Check for start
                    if not season_active and all(s <= 0 for s in snow_depth[k-2:k+1]):
                        season_active = True
                        the_day = day[k]
                        if not multi_year and mon[k] == 1 and day[k] == 4:
                            the_day = day[k-3]
                        season_start_end.append({
                            'yr': yr[k],
                            'mon': mon[k],
                            'day': the_day,
                            'fsdatetype': 'start'
                        })
                    # Check for end
                    if season_active and (snow_depth[k] > 0 or (mon[k] == 12 and all(t < fs_end for t in tmax[k-2:k+1]))):
                        season_active = False
                        season_start_end.append({
                            'yr': yr[k],
                            'mon': mon[k],
                            'day': day[k],
                            'fsdatetype': 'end'
                        })
        else:
            # Fall back to WF93
            return fire_season(input, fs_start, fs_end, "WF93", consistent_snow, multi_year)

    # Add date field and remove duplicates
    for entry in season_start_end:
        entry['date'] = datetime.date(entry['yr'], entry['mon'], entry['day'])

    # Remove duplicates (same date)
    seen_dates = set()
    unique_entries = []
    for entry in season_start_end:
        date_str = entry['date']
        if date_str not in seen_dates:
            seen_dates.add(date_str)
            unique_entries.append(entry)

    return unique_entries