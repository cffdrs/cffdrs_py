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

import math
from typing import Literal
import warnings
from cffdrs.fwi import initial_spread_index
from cffdrs.rate_of_spread import rate_of_spread_extended
from cffdrs.slope_calc import slope_adjustment
from cffdrs.surface_fuel_consumption import surface_fuel_consumption
from cffdrs.models import FBPInput, FBPPrimaryOutput, FBPSecondaryOutput, FBPAllOutput
from cffdrs.back_rate_of_spread import back_rate_of_spread
from cffdrs.buildup_effect import buildup_effect
from cffdrs.cfb_calc import crown_fraction_burned 
from cffdrs.crown_base_height import crown_base_height
from cffdrs.crown_fuel_load import crown_fuel_load
from cffdrs.distance_at_time import distance_at_time
from cffdrs.fire_intensity import fire_intensity
from cffdrs.flank_rate_of_spread import flank_rate_of_spread
from cffdrs.foliar_moisture_content import foliar_moisture_content
from cffdrs.length_to_breadth import length_to_breadth
from cffdrs.length_to_breadth_at_time import length_to_breadth_at_time
from cffdrs.rate_of_spread_at_time import rate_of_spread_at_time
from cffdrs.total_fuel_consumption import total_fuel_consumption


def fire_behaviour_prediction(input: FBPInput, output: Literal["Primary", "Secondary", "All"] = "Primary"):
    if not isinstance(input, FBPInput):
        input = FBPInput()
    # Unpack input
    fuel_type = input.fuel_type
    ffmc = input.ffmc
    bui = input.bui
    ws = input.ws
    wd = input.wd
    gs = input.gs
    aspect = input.aspect
    pc = input.pc
    pdf = input.pdf
    cc = input.cc
    gfl = input.gfl
    cbh = input.cbh
    cfl = input.cfl
    fmc = input.fmc
    isi = input.isi
    lat = input.lat
    long = input.lon
    elv = input.elv
    dj = input.dj
    d0 = input.d0
    sd = input.sd
    sh = input.sh
    hr = input.hr
    theta = input.theta
    accel = input.accel
    buieff = input.bui_eff
    id = input.id or 1
    output = output.upper()
    # Convert Wind Direction from degrees to radians
    wd_rad = math.radians(wd)
    # Convert Theta from degrees to radians
    theta_rad = math.radians(theta)
    aspect = 0 if math.isnan(aspect) else aspect
    aspect = aspect + 360 if aspect < 0 else aspect
    # Convert Aspect from degrees to radians
    aspect_rad = math.radians(aspect)
    accel = 0 if math.isnan(accel) or accel < 0 else accel
    if accel not in [0, 1]:
        warnings.warn("Input variable Accel is out of range, will be assigned to 1")
        accel = 1
    dj = 0 if dj < 0 or dj > 366 else dj
    dj = 180 if math.isnan(dj) else dj
    d0 = 0 if math.isnan(d0) or d0 < 0 or d0 > 366 else d0
    elv = 0 if elv < 0 or elv > 10000 else elv
    elv = 0 if math.isnan(elv) else elv
    buieff = 0 if buieff <= 0 else 1
    buieff = 1 if math.isnan(buieff) else buieff
    hr = -hr if hr < 0 else hr
    hr = 24 if hr > 366 * 24 else hr
    hr = 0 if math.isnan(hr) else hr
    ffmc = 0 if ffmc < 0 or ffmc > 101 else ffmc
    ffmc = 90 if math.isnan(ffmc) else ffmc
    isi = 0 if math.isnan(isi) or isi < 0 or isi > 300 else isi
    bui = 0 if bui < 0 or bui > 1000 else bui
    bui = 60 if math.isnan(bui) else bui
    ws = 0 if ws < 0 or ws > 300 else ws
    ws = 10 if math.isnan(ws) else ws
    wd_rad = 0 if math.isnan(wd_rad) or wd_rad < -2 * math.pi or wd_rad > 2 * math.pi else wd_rad
    gs = 0 if math.isnan(gs) or gs < 0 or gs > 200 else gs
    gs = 0 if aspect_rad < -2 * math.pi or aspect_rad > 2 * math.pi else gs
    pc = 50 if math.isnan(pc) or pc < 0 or pc > 100 else pc
    pdf = 35 if math.isnan(pdf) or pdf < 0 or pdf > 100 else pdf
    cc = 95 if cc <= 0 or cc > 100 else cc
    cc = 80 if math.isnan(cc) else cc
    gfl = 0.35 if math.isnan(gfl) or gfl <= 0 or gfl > 100 else gfl
    lat = 0 if lat < -90 or lat > 90 else lat
    lat = 55 if math.isnan(lat) else lat
    long = 0 if long < -180 or long > 360 else long
    long = -120 if math.isnan(long) else long
    theta_rad = 0 if math.isnan(theta_rad) or theta_rad < -2 * math.pi or theta_rad > 2 * math.pi else theta_rad
    sd = -999 if sd < 0 or sd > 1e5 else sd
    sd = 0 if math.isnan(sd) else sd
    sh = -999 if sh < 0 or sh > 100 else sh
    sh = 0 if math.isnan(sh) else sh

    # Convert hours to minutes
    hr_min = hr * 60
    # Corrections to reorient Wind Azimuth(WAZ) and Uphill slope azimuth(SAZ)
    waz = wd_rad + math.pi
    waz = waz - 2 * math.pi if waz > 2 * math.pi else waz
    saz = aspect_rad + math.pi
    saz = saz - 2 * math.pi if saz > 2 * math.pi else saz
    # Any negative longitudes (western hemisphere) are translated to positive
    #  longitudes
    long = -long if long < 0 else long

    # Initializing variables
    sfc = tfc = hfi = cfb = ros = 0
    raz = -999
    validOutTypes = ["SECONDARY", "ALL", "S", "A", "RAZ0", "WSV0"]
    if output in validOutTypes:
        fros = bros = tros = hrost = frost = brost = trost = fcfb = bcfb = tcfb = ffi = bfi = tfi = ftfc = btfc = ttfc = 0
        ti = fti = bti = tti = lb = wsv = -999

    cbh = crown_base_height(fuel_type, cbh, sd, sh)
    cfl = crown_fuel_load(fuel_type, cfl)
    fmc = foliar_moisture_content(lat, long, elv, dj, d0) if (fmc <= 0 or fmc > 120 or math.isnan(fmc)) else fmc
    fmc = 0 if fuel_type in ["D1", "S1", "S2", "S3", "O1A", "O1B"] else fmc

    # Calculate Surface fuel consumption (SFC)
    sfc = surface_fuel_consumption(fuel_type, ffmc, bui, pc, gfl)
    # Disable BUI Effect if necessary
    bui_eff = 0 if buieff != 1 else bui
    slope_values = slope_adjustment(fuel_type, ffmc, bui_eff, ws, waz, gs, saz, fmc, sfc, pc, pdf, cc, cbh, isi)
    # Calculate the net effective windspeed (WSV)
    wsv0 = slope_values["WSV"]
    if output == "WSV0":
        return wsv0
    wsv = wsv0 if gs > 0 and ffmc > 0 else ws
    # Calculate the net effective wind direction (RAZ)
    raz0 = slope_values["RAZ"]
    if output == "RAZ0":
        return raz0
    raz = raz0 if gs > 0 and ffmc > 0 else waz
    # Calculate or keep Initial Spread Index (ISI)
    isi = isi if isi > 0 else initial_spread_index(ffmc, wsv, True)
    # HACK: C6 ROS depends on CFB so do this to not repeat calculations
    ros_vars = rate_of_spread_extended(fuel_type, isi, bui_eff, fmc, sfc, pc, pdf, cc, cbh)
    ros = ros_vars["ROS"]
    cfb = ros_vars["CFB"] if cfl > 0 else 0
    csi = ros_vars["CSI"]
    rso = ros_vars["RSO"]
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
        sf = 10 if gs >= 70 else math.exp(3.533 * (gs / 100)**1.2)
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
            trost = (rost * (1 - math.sqrt(1 - 1 / lbt / lbt)) / (1 - math.sqrt(1 - 1 / lbt / lbt) * math.cos(theta_rad - raz)))
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
        fbp = FBPPrimaryOutput(id=id, cfb=cfb, cfc=cfc, fd=fd, hfi=hfi, raz=raz, ros=ros, sfc=sfc, tfc=tfc)
        if fuel_type in ["WA", "NF"]:
            fbp.cfb = 0
            fbp.cfc = 0
            fbp.hfi = 0
            fbp.raz = 0
            fbp.ros = 0
            fbp.sfc = 0
            fbp.tfc = 0
            fbp.fd = "NA"
    elif output in ["SECONDARY", "S"]:
        fbp = FBPSecondaryOutput(id=id, be=be, sf=sf, isi=isi, ffmc=ffmc, fmc=fmc, d0=d0, rso=rso,
            csi=csi, fros=fros, bros=bros, hrost=hrost, frost=frost, brost=brost, fcfb=fcfb, bcfb=bcfb,
            ffi=ffi, bfi=bfi, ftfc=ftfc, btfc=btfc, ti=ti, fti=fti, bti=bti, lb=lb, lbt=lbt, wsv=wsv,
            dh=dh, db=db, df=df, tros=tros, trost=trost, tcfb=tcfb, tfi=tfi, ttfc=ttfc, tti=tti)
        if fuel_type in ["WA", "NF"]:
            for field in fbp.__dataclass_fields__:
                if field != 'id':
                    setattr(fbp, field, 0)
    elif output in ["ALL", "A"]:
        fbp = FBPAllOutput(id=id, cfb=cfb, cfc=cfc, fd=fd, hfi=hfi, raz=raz, ros=ros, sfc=sfc,
            tfc=tfc, be=be, sf=sf, isi=isi, ffmc=ffmc, fmc=fmc, d0=d0, rso=rso, csi=csi, fros=fros,
            bros=bros, hrost=hrost, frost=frost, brost=brost, fcfb=fcfb, bcfb=bcfb, ffi=ffi, bfi=bfi,
            ftfc=ftfc, btfc=btfc, ti=ti, fti=fti, bti=bti, lb=lb, lbt=lbt, wsv=wsv, dh=dh, db=db, df=df,
            tros=tros, trost=trost, tcfb=tcfb, tfi=tfi, ttfc=ttfc, tti=tti)
        if fuel_type in ["WA", "NF"]:
            for field in ["cfb", "cfc", "hfi", "raz", "ros", "sfc", "tfc", "be", "sf", "isi", "ffmc", "fmc", "d0", "rso", "csi", "fros",
                        "bros", "hrost", "frost", "brost", "fcfb", "bcfb", "ffi", "bfi", "ftfc", "btfc", "ti", "fti", "bti", "lb", "lbt", "wsv", "dh", "db", "df",
                        "tros", "trost", "tcfb", "tfi", "ttfc", "tti"]:
                setattr(fbp, field, 0)
            fbp.fd = "NA"
    return fbp
