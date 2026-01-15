from dataclasses import dataclass
from typing import Optional

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
    dj: int = 180
    d0: int = 0
    hr: float = 1.0
    pc: float = 50.0
    pdf: float = 35.0
    cc: float = 80.0
    gfl: float = 0.35
    cbh: float = 0.0
    cfl: float = 0.0
    isi: float = 0.0
    fmc: Optional[float] = None
    theta: float = 0.0
    accel: int = 0
    aspect: float = 0.0
    bui_eff: int = 1
    id: Optional[str] = None
    sd: float = 0.0
    sh: float = 0.0

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