import math
from typing import Literal, NamedTuple
from cffdrs.constants import D1, S1, S2, S3, O1A, O1B, C6, NF, WA, FUEL_TYPE_CODES
from cffdrs.fwi import _initial_spread_index
from cffdrs.rate_of_spread import _rate_of_spread_extended
from cffdrs.slope_calc import _slope_adjustment
from cffdrs.surface_fuel_consumption import _surface_fuel_consumption
from cffdrs.models import FBPInput, FBPPrimaryOutput, FBPSecondaryOutput, FBPAllOutput
from cffdrs.back_rate_of_spread import _back_rate_of_spread
from cffdrs.buildup_effect import _buildup_effect
from cffdrs.cfb_calc import crown_fraction_burned
from cffdrs.crown_base_height import _crown_base_height
from cffdrs.crown_fuel_load import _crown_fuel_load
from cffdrs.distance_at_time import _distance_at_time
from cffdrs.fire_intensity import fire_intensity
from cffdrs.flank_rate_of_spread import flank_rate_of_spread
from cffdrs.foliar_moisture_content import foliar_moisture_content
from cffdrs.length_to_breadth import _length_to_breadth
from cffdrs.length_to_breadth_at_time import _length_to_breadth_at_time
from cffdrs.rate_of_spread_at_time import _rate_of_spread_at_time
from cffdrs.total_fuel_consumption import _total_fuel_consumption, _crown_fuel_consumption


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
    id = input.id or "1"
    output = output.upper()

    result = _fire_behaviour_prediction(
        FUEL_TYPE_CODES[input.fuel_type],
        input.ffmc,
        input.bui,
        input.ws,
        input.wd,  # Already in radians
        input.gs,
        input.aspect,  # Already in radians
        input.pc,
        input.pdf,
        input.cc,
        input.gfl,
        input.cbh,
        input.cfl,
        input.fmc,
        input.isi,
        input.lat,
        input.lon,  # Already adjusted for negative longitudes
        input.elv,
        input.dj,
        input.d0,
        input.sd,
        input.sh,
        input.hr,
        input.theta,  # Already in radians
        input.accel,
        input.bui_eff,
    )

    if output == "WSV0":
        return result.wsv0
    if output == "RAZ0":
        return result.raz0

    fd = {FD_SURFACE: "S", FD_INTERMITTENT: "I", FD_CROWN: "C", FD_NONE: None}[result.fd_code]

    if output in ["PRIMARY", "P"]:
        return FBPPrimaryOutput(
            id=id,
            cfb=result.cfb,
            cfc=result.cfc,
            fd=fd,
            hfi=result.hfi,
            raz=result.raz,
            ros=result.ros,
            sfc=result.sfc,
            tfc=result.tfc,
        )
    elif output in ["SECONDARY", "S"]:
        return FBPSecondaryOutput(
            id=id,
            be=result.be,
            sf=result.sf,
            isi=result.isi,
            ffmc=result.ffmc,
            fmc=result.fmc,
            d0=result.d0,
            rso=result.rso,
            csi=result.csi,
            fros=result.fros,
            bros=result.bros,
            hrost=result.hrost,
            frost=result.frost,
            brost=result.brost,
            fcfb=result.fcfb,
            bcfb=result.bcfb,
            ffi=result.ffi,
            bfi=result.bfi,
            ftfc=result.ftfc,
            btfc=result.btfc,
            ti=result.ti,
            fti=result.fti,
            bti=result.bti,
            lb=result.lb,
            lbt=result.lbt,
            wsv=result.wsv,
            dh=result.dh,
            db=result.db,
            df=result.df,
            tros=result.tros,
            trost=result.trost,
            tcfb=result.tcfb,
            tfi=result.tfi,
            ttfc=result.ttfc,
            tti=result.tti,
        )
    elif output in ["ALL", "A"]:
        return FBPAllOutput(
            id=id,
            cfb=result.cfb,
            cfc=result.cfc,
            fd=fd,
            hfi=result.hfi,
            raz=result.raz,
            ros=result.ros,
            sfc=result.sfc,
            tfc=result.tfc,
            be=result.be,
            sf=result.sf,
            isi=result.isi,
            ffmc=result.ffmc,
            fmc=result.fmc,
            d0=result.d0,
            rso=result.rso,
            csi=result.csi,
            fros=result.fros,
            bros=result.bros,
            hrost=result.hrost,
            frost=result.frost,
            brost=result.brost,
            fcfb=result.fcfb,
            bcfb=result.bcfb,
            ffi=result.ffi,
            bfi=result.bfi,
            ftfc=result.ftfc,
            btfc=result.btfc,
            ti=result.ti,
            fti=result.fti,
            bti=result.bti,
            lb=result.lb,
            lbt=result.lbt,
            wsv=result.wsv,
            dh=result.dh,
            db=result.db,
            df=result.df,
            tros=result.tros,
            trost=result.trost,
            tcfb=result.tcfb,
            tfi=result.tfi,
            ttfc=result.ttfc,
            tti=result.tti,
        )
    # Matches the original implementation's behavior for an unrecognized
    # output value (previously an accidental UnboundLocalError).
    raise ValueError(f"Invalid output type: {output}")


# Fire Type (fd) codes for _FBPOutput.fd_code - a plain int in place of the
# "S"/"I"/"C"/None string fbp() returns, since a string field has no place in
# a numeric, vectorization-ready output.
FD_SURFACE = 0
FD_INTERMITTENT = 1
FD_CROWN = 2
FD_NONE = -1  # non-fuel or unknown fuel types (NF, WA) - mirrors fd=None today


class _FBPOutput(NamedTuple):
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
    # Raw slope_adjustment() outputs, ahead of the wsv/raz fallback-to-ws/waz
    # logic and the raz degrees/mod-360 conversion below. Not part of
    # Primary/Secondary/All - exposed so fire_behaviour_prediction()'s
    # output="WSV0"/"RAZ0" shortcuts can read them straight off this result
    # instead of recomputing.
    wsv0: float
    raz0: float


def _fire_behaviour_prediction(
    fuel_type_code: int,
    ffmc: float,
    bui: float,
    ws: float,
    wd_rad: float,
    gs: float,
    aspect_rad: float,
    pc: float,
    pdf: float,
    cc: float,
    gfl: float,
    cbh: float,
    cfl: float,
    fmc: float,
    isi: float,
    lat: float,
    lon: float,
    elv: float,
    dj: float,
    d0: float,
    sd: float,
    sh: float,
    hr: float,
    theta_rad: float,
    accel: int,
    buieff: int,
) -> _FBPOutput:
    """
    Vectorization-ready Fire Behaviour Prediction System Calculation.

    Same computation as fire_behaviour_prediction(input, "All"), but:
      - takes an int fuel_type_code (see cffdrs.constants.FUEL_TYPE_CODES)
        instead of a fuel type string, with every fuel-type-branching callee
        swapped for its leading-underscore counterpart
      - takes already-validated/range-clamped scalar inputs directly, rather
        than an FBPInput dataclass - FBPInput.__post_init__ does that
        clamping (and unit conversion to radians for wd/aspect/theta) via
        plain Python (isinstance checks, warnings.warn, string parsing) that
        can't be traced/compiled by array-vectorization tools; a caller
        vectorizing over arrays of inputs is expected to do the equivalent
        clamping once, up front, over the whole array
      - has no `id` field (not numeric)
      - always computes and returns every field, since a batched kernel can't
        cheaply take a different code path per element - fire_behaviour_prediction()
        below delegates here unconditionally too, then just selects which
        fields go into the Primary/Secondary/All dataclass it returns
      - represents Fire Type as an int fd_code (FD_SURFACE/FD_INTERMITTENT/
        FD_CROWN/FD_NONE) instead of a "S"/"I"/"C"/None string
      - includes wsv0/raz0 (the raw slope_adjustment() outputs) so that
        fire_behaviour_prediction()'s output="WSV0"/"RAZ0" shortcuts can read
        them straight off this result instead of recomputing

    :returns: _FBPOutput with all Primary + Secondary fields (fd as fd_code)
    """
    # Convert hours to minutes
    hr_min = hr * 60
    # Corrections to reorient Wind Azimuth(WAZ) and Uphill slope azimuth(SAZ)
    waz = wd_rad + math.pi
    waz = waz - 2 * math.pi if waz > 2 * math.pi else waz
    saz = aspect_rad + math.pi
    saz = saz - 2 * math.pi if saz > 2 * math.pi else saz

    cbh = _crown_base_height(fuel_type_code, cbh, sd, sh)
    cfl = _crown_fuel_load(fuel_type_code, cfl)
    fmc = (
        foliar_moisture_content(lat, lon, elv, dj, d0)
        if (fmc <= 0 or fmc > 120 or math.isnan(fmc))
        else fmc
    )
    fmc = 0 if fuel_type_code in (D1, S1, S2, S3, O1A, O1B) else fmc

    # Calculate Surface fuel consumption (SFC)
    sfc = _surface_fuel_consumption(fuel_type_code, ffmc, bui, pc, gfl)
    # Disable BUI Effect if necessary
    bui_eff = 0 if buieff != 1 else bui
    slope_values = _slope_adjustment(
        fuel_type_code, ffmc, bui_eff, ws, waz, gs, saz, fmc, sfc, pc, pdf, cc, cbh, isi
    )
    # Calculate the net effective windspeed (WSV)
    wsv0 = slope_values.wsv
    wsv = wsv0 if gs > 0 and ffmc > 0 else ws
    # Calculate the net effective wind direction (RAZ)
    raz0 = slope_values.raz
    raz = raz0 if gs > 0 and ffmc > 0 else waz
    # Calculate or keep Initial Spread Index (ISI)
    isi = isi if isi > 0 else _initial_spread_index(ffmc, wsv, True)
    # HACK: C6 ROS depends on CFB so do this to not repeat calculations
    ros_vars = _rate_of_spread_extended(fuel_type_code, isi, bui_eff, fmc, sfc, pc, pdf, cc, cbh)
    ros = ros_vars.ros
    cfb = ros_vars.cfb if cfl > 0 else 0
    csi = ros_vars.csi
    rso = ros_vars.rso
    # Calculate Total Fuel Consumption (TFC)
    tfc = _total_fuel_consumption(fuel_type_code, cfl, cfb, sfc, pc, pdf)
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
    cfc = _crown_fuel_consumption(fuel_type_code, cfl, cfb, pc, pdf)

    # Calculate the Secondary Outputs (always computed - see docstring)
    # Eq. 39 (FCFDG 1992) Calculate Spread Factor (GS is group slope)
    sf = 10 if gs >= 70 else math.exp(3.533 * (gs / 100) ** 1.2)
    # Calculate The Buildup Effect
    be = _buildup_effect(fuel_type_code, bui_eff)
    # Calculate length to breadth ratio
    lb = _length_to_breadth(fuel_type_code, wsv)
    lbt = lb if accel == 0 else _length_to_breadth_at_time(fuel_type_code, lb, hr_min, cfb)
    # Calculate Back fire rate of spread (BROS)
    bros = _back_rate_of_spread(fuel_type_code, ffmc, bui_eff, wsv, fmc, sfc, pc, pdf, cc, cbh)
    # Calculate Flank fire rate of spread (FROS)
    fros = flank_rate_of_spread(ros, bros, lb)
    # Calculate the eccentricity
    e = math.sqrt(1 - 1 / lb / lb)
    # Calculate the rate of spread towards angle theta (TROS)
    tros = ros * (1 - e) / (1 - e * math.cos(theta_rad - raz))
    # Calculate rate of spread at time t for Flank, Back of fire and at angle
    # theta.
    rost = ros if accel == 0 else _rate_of_spread_at_time(fuel_type_code, ros, hr_min, cfb)
    brost = bros if accel == 0 else _rate_of_spread_at_time(fuel_type_code, bros, hr_min, cfb)
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
    ftfc = _total_fuel_consumption(fuel_type_code, cfl, fcfb, sfc, pc, pdf)
    btfc = _total_fuel_consumption(fuel_type_code, cfl, bcfb, sfc, pc, pdf)
    ttfc = _total_fuel_consumption(fuel_type_code, cfl, tcfb, sfc, pc, pdf)
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
    dh = _distance_at_time(fuel_type_code, ros, hr_min, cfb) if accel == 1 else ros * hr_min
    db = _distance_at_time(fuel_type_code, bros, hr_min, cfb) if accel == 1 else bros * hr_min
    df = (dh + db) / (lbt * 2) if accel == 1 else (dh + db) / (lb * 2)

    if fuel_type_code in (NF, WA):
        return _FBPOutput(
            cfb=0.0,
            cfc=0.0,
            fd_code=FD_NONE,
            hfi=0.0,
            raz=0.0,
            ros=0.0,
            sfc=0.0,
            tfc=0.0,
            be=0.0,
            sf=0.0,
            isi=0.0,
            ffmc=0.0,
            fmc=0.0,
            d0=0.0,
            rso=0.0,
            csi=0.0,
            fros=0.0,
            bros=0.0,
            hrost=0.0,
            frost=0.0,
            brost=0.0,
            fcfb=0.0,
            bcfb=0.0,
            ffi=0.0,
            bfi=0.0,
            ftfc=0.0,
            btfc=0.0,
            ti=0.0,
            fti=0.0,
            bti=0.0,
            lb=0.0,
            lbt=0.0,
            wsv=0.0,
            dh=0.0,
            db=0.0,
            df=0.0,
            tros=0.0,
            trost=0.0,
            tcfb=0.0,
            tfi=0.0,
            ttfc=0.0,
            tti=0.0,
            # wsv0/raz0 bypass the NF/WA zero-out, same as the WSV0/RAZ0
            # shortcuts in fire_behaviour_prediction() returning early before
            # that zero-out logic ever runs.
            wsv0=wsv0,
            raz0=raz0,
        )

    return _FBPOutput(
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
        wsv0=wsv0,
        raz0=raz0,
    )
