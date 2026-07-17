import math
from typing import Literal, NamedTuple
from cffdrs.constants import D1, S1, S2, S3, O1A, O1B, C6, NF, WA
from cffdrs.fwi import initial_spread_index, initial_spread_index_core
from cffdrs.rate_of_spread import rate_of_spread_extended, rate_of_spread_extended_core
from cffdrs.slope_calc import slope_adjustment, slope_adjustment_core
from cffdrs.surface_fuel_consumption import surface_fuel_consumption, surface_fuel_consumption_core
from cffdrs.models import FBPInput, FBPPrimaryOutput, FBPSecondaryOutput, FBPAllOutput
from cffdrs.back_rate_of_spread import back_rate_of_spread, back_rate_of_spread_core
from cffdrs.buildup_effect import buildup_effect, buildup_effect_core
from cffdrs.cfb_calc import crown_fraction_burned
from cffdrs.crown_base_height import crown_base_height, crown_base_height_core
from cffdrs.crown_fuel_load import crown_fuel_load, crown_fuel_load_core
from cffdrs.distance_at_time import distance_at_time, distance_at_time_core
from cffdrs.fire_intensity import fire_intensity
from cffdrs.flank_rate_of_spread import flank_rate_of_spread
from cffdrs.foliar_moisture_content import foliar_moisture_content
from cffdrs.length_to_breadth import length_to_breadth, length_to_breadth_core
from cffdrs.length_to_breadth_at_time import (
    length_to_breadth_at_time,
    length_to_breadth_at_time_core,
)
from cffdrs.rate_of_spread_at_time import rate_of_spread_at_time, rate_of_spread_at_time_core
from cffdrs.total_fuel_consumption import (
    total_fuel_consumption,
    total_fuel_consumption_core,
    crown_fuel_consumption_core,
)


def fire_behaviour_prediction(
    input: FBPInput, output: Literal["Primary", "Secondary", "All"] = "Primary"
):
    """
    Fire Behaviour Prediction System Calculation

    Fire Behavior Prediction System calculations. This is the
    primary function for calculating FBP for a single timestep. Not all
    equations are calculated within this function, but have been broken down
    further.

    :param input: FBPInput dataclass with input parameters
    :param output: What fbp outputs to return. Options are "Primary", "Secondary" and "All". Default: "Primary"

    :returns: FBPPrimaryOutput, FBPSecondaryOutput, or FBPAllOutput depending on output parameter
    """
    if not isinstance(input, FBPInput):
        input = FBPInput()
    # Unpack input (already validated and converted in FBPInput.__post_init__)
    fuel_type = input.fuel_type
    ffmc = input.ffmc
    bui = input.bui
    ws = input.ws
    wd_rad = input.wd  # Already in radians
    gs = input.gs
    aspect_rad = input.aspect  # Already in radians
    pc = input.pc
    pdf = input.pdf
    cc = input.cc
    gfl = input.gfl
    cbh = input.cbh
    cfl = input.cfl
    fmc = input.fmc
    isi = input.isi
    lat = input.lat
    long = input.lon  # Already adjusted for negative longitudes
    elv = input.elv
    dj = input.dj
    d0 = input.d0
    sd = input.sd
    sh = input.sh
    hr = input.hr
    theta_rad = input.theta  # Already in radians
    accel = input.accel
    buieff = input.bui_eff
    id = input.id or "1"
    output = output.upper()

    # Convert hours to minutes
    hr_min = hr * 60
    # Corrections to reorient Wind Azimuth(WAZ) and Uphill slope azimuth(SAZ)
    waz = wd_rad + math.pi
    waz = waz - 2 * math.pi if waz > 2 * math.pi else waz
    saz = aspect_rad + math.pi
    saz = saz - 2 * math.pi if saz > 2 * math.pi else saz

    # Initializing variables
    sfc = tfc = hfi = cfb = ros = 0
    raz = -999
    validOutTypes = ["SECONDARY", "ALL", "S", "A", "RAZ0", "WSV0"]
    if output in validOutTypes:
        fros = bros = tros = hrost = frost = brost = trost = fcfb = bcfb = tcfb = ffi = bfi = (
            tfi
        ) = ftfc = btfc = ttfc = 0
        ti = fti = bti = tti = lb = wsv = -999

    cbh = crown_base_height(fuel_type, cbh, sd, sh)
    cfl = crown_fuel_load(fuel_type, cfl)
    fmc = (
        foliar_moisture_content(lat, long, elv, dj, d0)
        if (fmc <= 0 or fmc > 120 or math.isnan(fmc))
        else fmc
    )
    fmc = 0 if fuel_type in ["D1", "S1", "S2", "S3", "O1A", "O1B"] else fmc

    # Calculate Surface fuel consumption (SFC)
    sfc = surface_fuel_consumption(fuel_type, ffmc, bui, pc, gfl)
    # Disable BUI Effect if necessary
    bui_eff = 0 if buieff != 1 else bui
    slope_values = slope_adjustment(
        fuel_type, ffmc, bui_eff, ws, waz, gs, saz, fmc, sfc, pc, pdf, cc, cbh, isi
    )
    # Calculate the net effective windspeed (WSV)
    wsv0 = slope_values.wsv
    if output == "WSV0":
        return wsv0
    wsv = wsv0 if gs > 0 and ffmc > 0 else ws
    # Calculate the net effective wind direction (RAZ)
    raz0 = slope_values.raz
    if output == "RAZ0":
        return raz0
    raz = raz0 if gs > 0 and ffmc > 0 else waz
    # Calculate or keep Initial Spread Index (ISI)
    isi = isi if isi > 0 else initial_spread_index(ffmc, wsv, True)
    # HACK: C6 ROS depends on CFB so do this to not repeat calculations
    ros_vars = rate_of_spread_extended(fuel_type, isi, bui_eff, fmc, sfc, pc, pdf, cc, cbh)
    ros = ros_vars.ros
    cfb = ros_vars.cfb if cfl > 0 else 0
    csi = ros_vars.csi
    rso = ros_vars.rso
    # Calculate Total Fuel Consumption (TFC)
    tfc = total_fuel_consumption(fuel_type, cfl, cfb, sfc, pc, pdf)
    # Calculate Head Fire Intensity(HFI)
    hfi = fire_intensity(tfc, ros)
    # Adjust Crown Fraction Burned
    cfb = -cfb if hr < 0 else cfb
    # Adjust RAZ
    raz_deg = math.degrees(raz)
    raz = 0 if raz_deg == 360 else raz_deg
    # Calculate Fire Type (S = Surface, C = Crowning, I = Intermittent Crowning)
    fd = "I"
    if cfb < 0.1:
        fd = "S"
    elif cfb >= 0.9:
        fd = "C"
    # Calculate Crown Fuel Consumption(CFC)
    cfc = total_fuel_consumption(fuel_type, cfl, cfb, sfc, pc, pdf, option="CFC")
    # Calculate the Secondary Outputs
    if output in ["SECONDARY", "ALL", "S", "A"]:
        # Eq. 39 (FCFDG 1992) Calculate Spread Factor (GS is group slope)
        sf = 10 if gs >= 70 else math.exp(3.533 * (gs / 100) ** 1.2)
        # Calculate The Buildup Effect
        be = buildup_effect(fuel_type, bui_eff)
        # Calculate length to breadth ratio
        lb = length_to_breadth(fuel_type, wsv)
        lbt = lb if accel == 0 else length_to_breadth_at_time(fuel_type, lb, hr_min, cfb)
        # Calculate Back fire rate of spread (BROS)
        bros = back_rate_of_spread(fuel_type, ffmc, bui_eff, wsv, fmc, sfc, pc, pdf, cc, cbh)
        # Calculate Flank fire rate of spread (FROS)
        fros = flank_rate_of_spread(ros, bros, lb)
        # Calculate the eccentricity
        e = math.sqrt(1 - 1 / lb / lb)
        # Calculate the rate of spread towards angle theta (TROS)
        tros = ros * (1 - e) / (1 - e * math.cos(theta_rad - raz))
        # Calculate rate of spread at time t for Flank, Back of fire and at angle
        # theta.
        rost = ros if accel == 0 else rate_of_spread_at_time(fuel_type, ros, hr_min, cfb)
        brost = bros if accel == 0 else rate_of_spread_at_time(fuel_type, bros, hr_min, cfb)
        frost = fros if accel == 0 else flank_rate_of_spread(rost, brost, lbt)
        # Calculate rate of spread towards angle theta at time t (TROSt)
        if accel == 0:
            trost = tros
        else:
            trost = (
                rost
                * (1 - math.sqrt(1 - 1 / lbt / lbt))
                / (1 - math.sqrt(1 - 1 / lbt / lbt) * math.cos(theta_rad - raz))
            )
        # Calculate Crown Fraction Burned for Flank, Back of fire and angle theta.
        fcfb = 0 if cfl == 0 else (0 if fuel_type == "C6" else crown_fraction_burned(fros, rso))
        bcfb = 0 if cfl == 0 else (0 if fuel_type == "C6" else crown_fraction_burned(bros, rso))
        tcfb = 0 if cfl == 0 else (0 if fuel_type == "C6" else crown_fraction_burned(tros, rso))
        # Calculate Total fuel consumption for the Flank fire, Back fire and at
        #  angle theta
        ftfc = total_fuel_consumption(fuel_type, cfl, fcfb, sfc, pc, pdf)
        btfc = total_fuel_consumption(fuel_type, cfl, bcfb, sfc, pc, pdf)
        ttfc = total_fuel_consumption(fuel_type, cfl, tcfb, sfc, pc, pdf)
        # Calculate the Fire Intensity at the Flank, Back and at angle theta fire
        ffi = fire_intensity(ftfc, fros)
        bfi = fire_intensity(btfc, bros)
        tfi = fire_intensity(ttfc, tros)
        # Calculate Rate of spread at time t for the Head, Flank, Back of fire and
        #  at angle theta.
        hrost = -rost if hr < 0 else rost
        frost = -frost if hr < 0 else frost
        brost = -brost if hr < 0 else brost
        trost = -trost if hr < 0 else trost

        # Calculate the elapsed time to crown fire initiation for Head, Flank, Back
        # fire and at angle theta. The (a# variable is a constant for Head, Flank,
        # Back and at angle theta used in the *TI equations)
        a1 = 0.115 - (18.8 * cfb**2.5 * math.exp(-8 * cfb))
        ti = math.log(1 if 1 - rso / ros <= 0 else 1 - rso / ros) / (-a1)
        a2 = 0.115 - (18.8 * fcfb**2.5 * math.exp(-8 * fcfb))
        fti = math.log(1 if 1 - rso / fros <= 0 else 1 - rso / fros) / (-a2)
        a3 = 0.115 - (18.8 * bcfb**2.5 * math.exp(-8 * bcfb))
        bti = math.log(1 if 1 - rso / bros <= 0 else 1 - rso / bros) / (-a3)
        a4 = 0.115 - (18.8 * tcfb**2.5 * math.exp(-8 * tcfb))
        tti = math.log(1 if 1 - rso / tros <= 0 else 1 - rso / tros) / (-a4)

        # Fire spread distance for Head, Back, and Flank of fire
        dh = distance_at_time(fuel_type, ros, hr_min, cfb) if accel == 1 else ros * hr_min
        db = distance_at_time(fuel_type, bros, hr_min, cfb) if accel == 1 else bros * hr_min
        df = (dh + db) / (lbt * 2) if accel == 1 else (dh + db) / (lb * 2)

    # Prepare output
    if output in ["PRIMARY", "P"]:
        fbp = FBPPrimaryOutput(
            id=id, cfb=cfb, cfc=cfc, fd=fd, hfi=hfi, raz=raz, ros=ros, sfc=sfc, tfc=tfc
        )
    elif output in ["SECONDARY", "S"]:
        fbp = FBPSecondaryOutput(
            id=id,
            be=be,
            sf=sf,
            isi=isi,
            ffmc=ffmc,
            fmc=fmc,
            d0=d0,
            rso=rso,
            csi=csi,
            fros=fros,
            bros=bros,
            hrost=hrost,
            frost=frost,
            brost=brost,
            fcfb=fcfb,
            bcfb=bcfb,
            ffi=ffi,
            bfi=bfi,
            ftfc=ftfc,
            btfc=btfc,
            ti=ti,
            fti=fti,
            bti=bti,
            lb=lb,
            lbt=lbt,
            wsv=wsv,
            dh=dh,
            db=db,
            df=df,
            tros=tros,
            trost=trost,
            tcfb=tcfb,
            tfi=tfi,
            ttfc=ttfc,
            tti=tti,
        )
    elif output in ["ALL", "A"]:
        fbp = FBPAllOutput(
            id=id,
            cfb=cfb,
            cfc=cfc,
            fd=fd,
            hfi=hfi,
            raz=raz,
            ros=ros,
            sfc=sfc,
            tfc=tfc,
            be=be,
            sf=sf,
            isi=isi,
            ffmc=ffmc,
            fmc=fmc,
            d0=d0,
            rso=rso,
            csi=csi,
            fros=fros,
            bros=bros,
            hrost=hrost,
            frost=frost,
            brost=brost,
            fcfb=fcfb,
            bcfb=bcfb,
            ffi=ffi,
            bfi=bfi,
            ftfc=ftfc,
            btfc=btfc,
            ti=ti,
            fti=fti,
            bti=bti,
            lb=lb,
            lbt=lbt,
            wsv=wsv,
            dh=dh,
            db=db,
            df=df,
            tros=tros,
            trost=trost,
            tcfb=tcfb,
            tfi=tfi,
            ttfc=ttfc,
            tti=tti,
        )
    if fuel_type in ["WA", "NF"]:
        for field in fbp.__dataclass_fields__:
            if field == "id":
                continue
            elif field == "fd":
                fbp.fd = None
            else:
                setattr(fbp, field, 0)
    return fbp


# Fire Type (fd) codes for FBPCoreOutput.fd_code - a plain int in place of the
# "S"/"I"/"C"/None string fbp() returns, since a string field has no place in
# a numeric, vectorization-ready output.
FD_SURFACE = 0
FD_INTERMITTENT = 1
FD_CROWN = 2
FD_NONE = -1  # non-fuel or unknown fuel types (NF, WA) - mirrors fd=None today


class FBPCoreOutput(NamedTuple):
    cfb: float
    cfc: float
    fd_code: int
    hfi: float
    raz: float
    ros: float
    sfc: float
    tfc: float
    be: float
    sf: float
    isi: float
    ffmc: float
    fmc: float
    d0: float
    rso: float
    csi: float
    fros: float
    bros: float
    hrost: float
    frost: float
    brost: float
    fcfb: float
    bcfb: float
    ffi: float
    bfi: float
    ftfc: float
    btfc: float
    ti: float
    fti: float
    bti: float
    lb: float
    lbt: float
    wsv: float
    dh: float
    db: float
    df: float
    tros: float
    trost: float
    tcfb: float
    tfi: float
    ttfc: float
    tti: float


def fire_behaviour_prediction_core(
    fuel_type_code: int,
    ffmc,
    bui,
    ws,
    wd_rad,
    gs,
    aspect_rad,
    pc,
    pdf,
    cc,
    gfl,
    cbh,
    cfl,
    fmc,
    isi,
    lat,
    lon,
    elv,
    dj,
    d0,
    sd,
    sh,
    hr,
    theta_rad,
    accel,
    buieff,
):
    """
    Vectorization-ready Fire Behaviour Prediction System Calculation.

    Same computation as fire_behaviour_prediction(input, "All"), but:
      - takes an int fuel_type_code (see cffdrs.constants.FUEL_TYPE_CODES)
        instead of a fuel type string, with every fuel-type-branching callee
        swapped for its "_core" counterpart
      - takes already-validated/range-clamped scalar inputs directly, rather
        than an FBPInput dataclass - FBPInput.__post_init__ does that
        clamping (and unit conversion to radians for wd/aspect/theta) via
        plain Python (isinstance checks, warnings.warn, string parsing) that
        can't be traced/compiled by array-vectorization tools; a caller
        vectorizing over arrays of inputs is expected to do the equivalent
        clamping once, up front, over the whole array
      - has no `id` field (not numeric)
      - always computes the full secondary-output set - a batched kernel
        can't cheaply take a different code path per element the way the
        output="Primary" early-exit does today, so this always returns
        everything "All" would; a caller just ignores fields it doesn't need
      - represents Fire Type as an int fd_code (FD_SURFACE/FD_INTERMITTENT/
        FD_CROWN/FD_NONE) instead of a "S"/"I"/"C"/None string
      - does not reproduce the output="WSV0"/"RAZ0" shortcuts (call
        slope_adjustment_core directly for those two raw values)

    :returns: FBPCoreOutput with all Primary + Secondary fields (fd as fd_code)
    """
    # Convert hours to minutes
    hr_min = hr * 60
    # Corrections to reorient Wind Azimuth(WAZ) and Uphill slope azimuth(SAZ)
    waz = wd_rad + math.pi
    waz = waz - 2 * math.pi if waz > 2 * math.pi else waz
    saz = aspect_rad + math.pi
    saz = saz - 2 * math.pi if saz > 2 * math.pi else saz

    cbh = crown_base_height_core(fuel_type_code, cbh, sd, sh)
    cfl = crown_fuel_load_core(fuel_type_code, cfl)
    fmc = (
        foliar_moisture_content(lat, lon, elv, dj, d0)
        if (fmc <= 0 or fmc > 120 or math.isnan(fmc))
        else fmc
    )
    fmc = 0 if fuel_type_code in (D1, S1, S2, S3, O1A, O1B) else fmc

    # Calculate Surface fuel consumption (SFC)
    sfc = surface_fuel_consumption_core(fuel_type_code, ffmc, bui, pc, gfl)
    # Disable BUI Effect if necessary
    bui_eff = 0 if buieff != 1 else bui
    slope_values = slope_adjustment_core(
        fuel_type_code, ffmc, bui_eff, ws, waz, gs, saz, fmc, sfc, pc, pdf, cc, cbh, isi
    )
    # Calculate the net effective windspeed (WSV)
    wsv0 = slope_values.wsv
    wsv = wsv0 if gs > 0 and ffmc > 0 else ws
    # Calculate the net effective wind direction (RAZ)
    raz0 = slope_values.raz
    raz = raz0 if gs > 0 and ffmc > 0 else waz
    # Calculate or keep Initial Spread Index (ISI)
    isi = isi if isi > 0 else initial_spread_index_core(ffmc, wsv, True)
    # HACK: C6 ROS depends on CFB so do this to not repeat calculations
    ros_vars = rate_of_spread_extended_core(
        fuel_type_code, isi, bui_eff, fmc, sfc, pc, pdf, cc, cbh
    )
    ros = ros_vars.ros
    cfb = ros_vars.cfb if cfl > 0 else 0
    csi = ros_vars.csi
    rso = ros_vars.rso
    # Calculate Total Fuel Consumption (TFC)
    tfc = total_fuel_consumption_core(fuel_type_code, cfl, cfb, sfc, pc, pdf)
    # Calculate Head Fire Intensity(HFI)
    hfi = fire_intensity(tfc, ros)
    # Adjust Crown Fraction Burned
    cfb = -cfb if hr < 0 else cfb
    # Adjust RAZ
    raz_deg = math.degrees(raz)
    raz = 0 if raz_deg == 360 else raz_deg
    # Calculate Fire Type (S = Surface, C = Crowning, I = Intermittent Crowning)
    fd_code = FD_INTERMITTENT
    if cfb < 0.1:
        fd_code = FD_SURFACE
    elif cfb >= 0.9:
        fd_code = FD_CROWN
    # Calculate Crown Fuel Consumption(CFC)
    cfc = crown_fuel_consumption_core(fuel_type_code, cfl, cfb, pc, pdf)

    # Calculate the Secondary Outputs (always computed - see docstring)
    # Eq. 39 (FCFDG 1992) Calculate Spread Factor (GS is group slope)
    sf = 10 if gs >= 70 else math.exp(3.533 * (gs / 100) ** 1.2)
    # Calculate The Buildup Effect
    be = buildup_effect_core(fuel_type_code, bui_eff)
    # Calculate length to breadth ratio
    lb = length_to_breadth_core(fuel_type_code, wsv)
    lbt = lb if accel == 0 else length_to_breadth_at_time_core(fuel_type_code, lb, hr_min, cfb)
    # Calculate Back fire rate of spread (BROS)
    bros = back_rate_of_spread_core(fuel_type_code, ffmc, bui_eff, wsv, fmc, sfc, pc, pdf, cc, cbh)
    # Calculate Flank fire rate of spread (FROS)
    fros = flank_rate_of_spread(ros, bros, lb)
    # Calculate the eccentricity
    e = math.sqrt(1 - 1 / lb / lb)
    # Calculate the rate of spread towards angle theta (TROS)
    tros = ros * (1 - e) / (1 - e * math.cos(theta_rad - raz))
    # Calculate rate of spread at time t for Flank, Back of fire and at angle
    # theta.
    rost = ros if accel == 0 else rate_of_spread_at_time_core(fuel_type_code, ros, hr_min, cfb)
    brost = bros if accel == 0 else rate_of_spread_at_time_core(fuel_type_code, bros, hr_min, cfb)
    frost = fros if accel == 0 else flank_rate_of_spread(rost, brost, lbt)
    # Calculate rate of spread towards angle theta at time t (TROSt)
    if accel == 0:
        trost = tros
    else:
        trost = (
            rost
            * (1 - math.sqrt(1 - 1 / lbt / lbt))
            / (1 - math.sqrt(1 - 1 / lbt / lbt) * math.cos(theta_rad - raz))
        )
    # Calculate Crown Fraction Burned for Flank, Back of fire and angle theta.
    fcfb = 0 if cfl == 0 else (0 if fuel_type_code == C6 else crown_fraction_burned(fros, rso))
    bcfb = 0 if cfl == 0 else (0 if fuel_type_code == C6 else crown_fraction_burned(bros, rso))
    tcfb = 0 if cfl == 0 else (0 if fuel_type_code == C6 else crown_fraction_burned(tros, rso))
    # Calculate Total fuel consumption for the Flank fire, Back fire and at
    #  angle theta
    ftfc = total_fuel_consumption_core(fuel_type_code, cfl, fcfb, sfc, pc, pdf)
    btfc = total_fuel_consumption_core(fuel_type_code, cfl, bcfb, sfc, pc, pdf)
    ttfc = total_fuel_consumption_core(fuel_type_code, cfl, tcfb, sfc, pc, pdf)
    # Calculate the Fire Intensity at the Flank, Back and at angle theta fire
    ffi = fire_intensity(ftfc, fros)
    bfi = fire_intensity(btfc, bros)
    tfi = fire_intensity(ttfc, tros)
    # Calculate Rate of spread at time t for the Head, Flank, Back of fire and
    #  at angle theta.
    hrost = -rost if hr < 0 else rost
    frost = -frost if hr < 0 else frost
    brost = -brost if hr < 0 else brost
    trost = -trost if hr < 0 else trost

    # Calculate the elapsed time to crown fire initiation for Head, Flank, Back
    # fire and at angle theta. The (a# variable is a constant for Head, Flank,
    # Back and at angle theta used in the *TI equations)
    a1 = 0.115 - (18.8 * cfb**2.5 * math.exp(-8 * cfb))
    ti = math.log(1 if 1 - rso / ros <= 0 else 1 - rso / ros) / (-a1)
    a2 = 0.115 - (18.8 * fcfb**2.5 * math.exp(-8 * fcfb))
    fti = math.log(1 if 1 - rso / fros <= 0 else 1 - rso / fros) / (-a2)
    a3 = 0.115 - (18.8 * bcfb**2.5 * math.exp(-8 * bcfb))
    bti = math.log(1 if 1 - rso / bros <= 0 else 1 - rso / bros) / (-a3)
    a4 = 0.115 - (18.8 * tcfb**2.5 * math.exp(-8 * tcfb))
    tti = math.log(1 if 1 - rso / tros <= 0 else 1 - rso / tros) / (-a4)

    # Fire spread distance for Head, Back, and Flank of fire
    dh = distance_at_time_core(fuel_type_code, ros, hr_min, cfb) if accel == 1 else ros * hr_min
    db = distance_at_time_core(fuel_type_code, bros, hr_min, cfb) if accel == 1 else bros * hr_min
    df = (dh + db) / (lbt * 2) if accel == 1 else (dh + db) / (lb * 2)

    if fuel_type_code in (NF, WA):
        return FBPCoreOutput(
            cfb=0,
            cfc=0,
            fd_code=FD_NONE,
            hfi=0,
            raz=0,
            ros=0,
            sfc=0,
            tfc=0,
            be=0,
            sf=0,
            isi=0,
            ffmc=0,
            fmc=0,
            d0=0,
            rso=0,
            csi=0,
            fros=0,
            bros=0,
            hrost=0,
            frost=0,
            brost=0,
            fcfb=0,
            bcfb=0,
            ffi=0,
            bfi=0,
            ftfc=0,
            btfc=0,
            ti=0,
            fti=0,
            bti=0,
            lb=0,
            lbt=0,
            wsv=0,
            dh=0,
            db=0,
            df=0,
            tros=0,
            trost=0,
            tcfb=0,
            tfi=0,
            ttfc=0,
            tti=0,
        )

    return FBPCoreOutput(
        cfb=cfb,
        cfc=cfc,
        fd_code=fd_code,
        hfi=hfi,
        raz=raz,
        ros=ros,
        sfc=sfc,
        tfc=tfc,
        be=be,
        sf=sf,
        isi=isi,
        ffmc=ffmc,
        fmc=fmc,
        d0=d0,
        rso=rso,
        csi=csi,
        fros=fros,
        bros=bros,
        hrost=hrost,
        frost=frost,
        brost=brost,
        fcfb=fcfb,
        bcfb=bcfb,
        ffi=ffi,
        bfi=bfi,
        ftfc=ftfc,
        btfc=btfc,
        ti=ti,
        fti=fti,
        bti=bti,
        lb=lb,
        lbt=lbt,
        wsv=wsv,
        dh=dh,
        db=db,
        df=df,
        tros=tros,
        trost=trost,
        tcfb=tcfb,
        tfi=tfi,
        ttfc=ttfc,
        tti=tti,
    )
