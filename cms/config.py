"""
Configuration parameters for the Cloud Modeling System (CMS).
Physical constants and grid parameters as defined in IMPLEMENTATION_GUIDE.md.
"""

from dataclasses import dataclass

# Physical Constants (Section 13)
VON_KARMAN = 0.4         # k [-]
RHO_A_STP = 1.225        # rho_0 [kg/m^3]
RHO_W = 1000.0           # rho_w [kg/m^3]
RHO_I = 917.0            # rho_i [kg/m^3]
L_V = 2.5e6              # Latent heat vaporization [J/kg]
L_F = 3.34e5             # Latent heat fusion [J/kg]
L_S = 2.834e6            # Latent heat sublimation [J/kg]
K_A = 2.4e-2             # Thermal conductivity air [W/m/K]
G = 9.81                 # Gravity [m/s^2]
CP = 1004.0              # Specific heat capacity [J/kg/K]
RD = 287.05              # Gas constant for dry air [J/kg/K]

# Grid Parameters (Section 1.2)
DEFAULT_NX = 101
DEFAULT_NY = 101
DEFAULT_NZ = 101
DEFAULT_DX = 1000.0      # [m]
DEFAULT_DY = 1000.0      # [m]
DEFAULT_DZ = 100.0       # [m]

@dataclass
class PhysicsConfig:
    von_karman: float = VON_KARMAN
    rho_a_stp: float = RHO_A_STP
    rho_w: float = RHO_W
    rho_i: float = RHO_I
    l_v: float = L_V
    l_f: float = L_F
    l_s: float = L_S
    k_a: float = K_A
    g: float = G
    cp: float = CP
    rd: float = RD

@dataclass
class GridConfig:
    nx: int = DEFAULT_NX
    ny: int = DEFAULT_NY
    nz: int = DEFAULT_NZ
    dx: float = DEFAULT_DX
    dy: float = DEFAULT_DY
    dz: float = DEFAULT_DZ
