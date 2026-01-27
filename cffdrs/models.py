from dataclasses import dataclass
from typing import Optional, get_args
import warnings
import math

from cffdrs.constants import FuelType


@dataclass
class FBPInput:
    fuel_type: FuelType = "C2"
    ffmc: float = 90.0
    bui: float = 60.0
    ws: float = 10.0
    wd: float = 0.0
    gs: float = 0.0
    lat: float = 55.0
    lon: float = -120.0
    elv: float = 0.0
    dj: float = 180.0
    d0: float = 0.0
    hr: float = 1.0
    pc: float = 50.0
    pdf: float = 35.0
    cc: float = 80.0
    gfl: float = 0.35
    cbh: float = 0.0
    cfl: float = 0.0
    isi: float = 0.0
    fmc: float = 0.0
    theta: float = 0.0
    accel: int = 0
    aspect: float = 0.0
    bui_eff: int = 1
    id: Optional[str] = None
    sd: float = 0.0
    sh: float = 0.0

    def __post_init__(self):
        # Fuel type normalization and validation
        original = self.fuel_type
        self.fuel_type = (
            self.fuel_type.upper().replace("-", "").replace(" ", "") if self.fuel_type else "C2"
        )
        if not self.fuel_type:
            warnings.warn("FuelType not provided, using C2 (default) in the calculation")
            self.fuel_type = "C2"

        # Validate against allowed fuel types
        allowed_fuel_types = get_args(FuelType)
        if self.fuel_type not in allowed_fuel_types:
            raise ValueError(
                f"Invalid fuel type '{original}' (normalized to '{self.fuel_type}'). Must be one of: {', '.join(allowed_fuel_types)}"
            )

        # Apply validations and defaults with warnings for changes
        original_accel = self.accel
        self.accel = 0 if math.isnan(self.accel) or self.accel < 0 else self.accel
        if self.accel not in [0, 1]:
            warnings.warn(
                f"Input variable Accel ({original_accel}) is out of range, will be assigned to 1"
            )
            self.accel = 1

        original_dj = self.dj
        self.dj = 0 if self.dj < 0 or self.dj > 366 else self.dj
        if self.dj != original_dj:
            warnings.warn(f"DJ ({original_dj}) out of range, clamped to {self.dj}")
        self.dj = 180 if math.isnan(self.dj) else self.dj

        original_d0 = self.d0
        self.d0 = 0 if math.isnan(self.d0) or self.d0 < 0 or self.d0 > 366 else self.d0
        if self.d0 != original_d0:
            warnings.warn(f"D0 ({original_d0}) out of range, clamped to {self.d0}")

        original_elv = self.elv
        self.elv = 0 if self.elv < 0 or self.elv > 10000 else self.elv
        if self.elv != original_elv:
            warnings.warn(f"ELV ({original_elv}) out of range, clamped to {self.elv}")
        self.elv = 0 if math.isnan(self.elv) else self.elv

        original_hr = self.hr
        self.hr = -self.hr if self.hr < 0 else self.hr
        if self.hr != original_hr:
            warnings.warn(f"HR ({original_hr}) negative, negated to {self.hr}")
        original_hr = self.hr
        self.hr = 24 if self.hr > 366 * 24 else self.hr
        if self.hr != original_hr:
            warnings.warn(f"HR ({original_hr}) too large, clamped to {self.hr}")
        self.hr = 0 if math.isnan(self.hr) else self.hr

        original_ffmc = self.ffmc
        self.ffmc = 0 if self.ffmc < 0 or self.ffmc > 101 else self.ffmc
        if self.ffmc != original_ffmc:
            warnings.warn(f"FFMC ({original_ffmc}) out of range, clamped to {self.ffmc}")
        self.ffmc = 90 if math.isnan(self.ffmc) else self.ffmc

        original_isi = self.isi
        self.isi = 0 if math.isnan(self.isi) or self.isi < 0 or self.isi > 300 else self.isi
        if self.isi != original_isi:
            warnings.warn(f"ISI ({original_isi}) out of range, clamped to {self.isi}")

        original_bui = self.bui
        self.bui = 0 if self.bui < 0 or self.bui > 1000 else self.bui
        if self.bui != original_bui:
            warnings.warn(f"BUI ({original_bui}) out of range, clamped to {self.bui}")
        self.bui = 60 if math.isnan(self.bui) else self.bui

        original_ws = self.ws
        self.ws = 0 if self.ws < 0 or self.ws > 300 else self.ws
        if self.ws != original_ws:
            warnings.warn(f"WS ({original_ws}) out of range, clamped to {self.ws}")
        self.ws = 10 if math.isnan(self.ws) else self.ws

        # Convert WD to radians and clamp
        self.wd = math.radians(self.wd)
        original_wd = self.wd
        self.wd = (
            0 if math.isnan(self.wd) or self.wd < -2 * math.pi or self.wd > 2 * math.pi else self.wd
        )
        if self.wd != original_wd:
            warnings.warn(f"WD (in radians) ({original_wd}) out of range, clamped to {self.wd}")

        original_gs = self.gs
        self.gs = 0 if math.isnan(self.gs) or self.gs < 0 or self.gs > 200 else self.gs
        if self.gs != original_gs:
            warnings.warn(f"GS ({original_gs}) out of range, clamped to {self.gs}")

        # Handle aspect
        self.aspect = 0 if math.isnan(self.aspect) else self.aspect
        original_aspect = self.aspect
        self.aspect = self.aspect + 360 if self.aspect < 0 else self.aspect
        if self.aspect != original_aspect:
            warnings.warn(f"ASPECT ({original_aspect}) negative, adjusted to {self.aspect}")
        self.aspect = math.radians(self.aspect)
        original_aspect_rad = self.aspect
        self.gs = 0 if self.aspect < -2 * math.pi or self.aspect > 2 * math.pi else self.gs
        if self.aspect != original_aspect_rad:
            warnings.warn(f"ASPECT (in radians) ({original_aspect_rad}) out of range, GS set to 0")

        original_pc = self.pc
        self.pc = 50 if math.isnan(self.pc) or self.pc < 0 or self.pc > 100 else self.pc
        if self.pc != original_pc:
            warnings.warn(f"PC ({original_pc}) out of range, clamped to {self.pc}")

        original_pdf = self.pdf
        self.pdf = 35 if math.isnan(self.pdf) or self.pdf < 0 or self.pdf > 100 else self.pdf
        if self.pdf != original_pdf:
            warnings.warn(f"PDF ({original_pdf}) out of range, clamped to {self.pdf}")

        original_cc = self.cc
        self.cc = 95 if self.cc <= 0 or self.cc > 100 else self.cc
        if self.cc != original_cc:
            warnings.warn(f"CC ({original_cc}) out of range, clamped to {self.cc}")
        self.cc = 80 if math.isnan(self.cc) else self.cc

        original_gfl = self.gfl
        self.gfl = 0.35 if math.isnan(self.gfl) or self.gfl <= 0 or self.gfl > 100 else self.gfl
        if self.gfl != original_gfl:
            warnings.warn(f"GFL ({original_gfl}) out of range, clamped to {self.gfl}")

        original_lat = self.lat
        self.lat = 0 if self.lat < -90 or self.lat > 90 else self.lat
        if self.lat != original_lat:
            warnings.warn(f"LAT ({original_lat}) out of range, clamped to {self.lat}")
        self.lat = 55 if math.isnan(self.lat) else self.lat

        original_lon = self.lon
        self.lon = 0 if self.lon < -180 or self.lon > 360 else self.lon
        if self.lon != original_lon:
            warnings.warn(f"LON ({original_lon}) out of range, clamped to {self.lon}")
        self.lon = -120 if math.isnan(self.lon) else self.lon

        # Convert theta to radians and clamp
        self.theta = math.radians(self.theta)
        original_theta = self.theta
        self.theta = (
            0
            if math.isnan(self.theta) or self.theta < -2 * math.pi or self.theta > 2 * math.pi
            else self.theta
        )
        if self.theta != original_theta:
            warnings.warn(
                f"THETA (in radians) ({original_theta}) out of range, clamped to {self.theta}"
            )

        original_sd = self.sd
        self.sd = -999 if self.sd < 0 or self.sd > 1e5 else self.sd
        if self.sd != original_sd:
            warnings.warn(f"SD ({original_sd}) out of range, clamped to {self.sd}")
        self.sd = 0 if math.isnan(self.sd) else self.sd

        original_sh = self.sh
        self.sh = -999 if self.sh < 0 or self.sh > 100 else self.sh
        if self.sh != original_sh:
            warnings.warn(f"SH ({original_sh}) out of range, clamped to {self.sh}")
        self.sh = 0 if math.isnan(self.sh) else self.sh

        # Corrections for longitudes
        self.lon = -self.lon if self.lon < 0 else self.lon


@dataclass
class FBPPrimaryOutput:
    id: str
    cfb: float
    cfc: float
    fd: str
    hfi: float
    raz: float
    ros: float
    sfc: float
    tfc: float


@dataclass
class FBPSecondaryOutput:
    id: str
    be: float
    sf: float
    isi: float
    ffmc: float
    fmc: float
    d0: int
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


@dataclass
class FBPAllOutput:
    id: str
    cfb: float
    cfc: float
    fd: str
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
    d0: int
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
