import numpy as np
from typing import Optional
from cms.config import PhysicsConfig

def compute_radar_reflectivity(
    rho: np.ndarray,
    qr: np.ndarray, nr: np.ndarray,
    qi: Optional[np.ndarray] = None, ni: Optional[np.ndarray] = None,
    qs: Optional[np.ndarray] = None, ns: Optional[np.ndarray] = None,
    qg: Optional[np.ndarray] = None, ng: Optional[np.ndarray] = None,
    config: Optional[PhysicsConfig] = None
) -> np.ndarray:
    """
    Computes the equivalent radar reflectivity factor (Ze) in dBZ.
    Ref: IMPLEMENTATION_GUIDE.md Eq 8.2
    
    The calculation assumes a Gamma distribution with alpha=1 (exponential)
    for all hydrometeor species and uses the Rayleigh scattering approximation.
    
    Args:
        rho: Air density (kg/m^3)
        qr, nr: Rain mixing ratio (kg/kg) and number concentration (m^-3)
        qi, ni: Cloud ice mixing ratio and number concentration
        qs, ns: Snow mixing ratio and number concentration
        qg, ng: Graupel mixing ratio and number concentration
        config: Physics configuration containing densities
        
    Returns:
        dbz: Reflectivity in dBZ (3D array)
    """
    if config is None:
        config = PhysicsConfig()
        
    # Dielectric factors |K|^2
    kw2 = 0.93  # Water
    ki2 = 0.176 # Ice
    
    z_total = np.zeros_like(rho)
    
    def calculate_moment_6(q, n, rho_h, rho_a):
        """Calculates the 6th moment of exponential distribution."""
        mask = (q > 1e-12) & (n > 1e-3)
        m6 = np.zeros_like(q)
        # For exponential: M6 = 720 / lambda^6
        # where lambda = (pi * rho_h * n / (rho_a * q))^(1/3)
        # M6 = 720 * (rho_a * q / (pi * rho_h))**2 / n
        m6[mask] = (720.0 * (rho_a[mask] * q[mask])**2) / (
            (np.pi**2) * (rho_h**2) * n[mask]
        )
        return m6

    # 1. Rain
    z_total += calculate_moment_6(qr, nr, config.rho_w, rho)
    
    # 2. Ice
    if qi is not None and ni is not None:
        z_total += (ki2 / kw2) * calculate_moment_6(qi, ni, config.rho_i, rho)
        
    # 3. Snow (Assuming rho_s = 100 kg/m^3 if not in config)
    if qs is not None and ns is not None:
        rho_s = getattr(config, 'rho_s', 100.0)
        z_total += (ki2 / kw2) * calculate_moment_6(qs, ns, rho_s, rho)
        
    # 4. Graupel (Assuming rho_g = 400 kg/m^3 if not in config)
    if qg is not None and ng is not None:
        rho_g = getattr(config, 'rho_g', 400.0)
        z_total += (ki2 / kw2) * calculate_moment_6(qg, ng, rho_g, rho)
        
    # Convert to dBZ: 10 * log10( Z_mm6_m3 )
    # Z here is in m^6/m^3, so multiply by 10^18 to get mm^6/m^3
    z_total_mm6 = z_total * 1e18
    # Clipping to -30 dBZ floor to avoid log(0)
    z_total_mm6 = np.maximum(z_total_mm6, 1e-3) 
    
    return 10.0 * np.log10(z_total_mm6)

def compute_seeding_metrics(
    precip_seeded: np.ndarray,
    precip_control: np.ndarray,
    ni_seeded: np.ndarray,
    ni_control: np.ndarray,
    reagent_mass: float
) -> dict:
    """
    Computes seeding efficiency metrics.
    Ref: IMPLEMENTATION_GUIDE.md Section 8.3
    """
    metrics = {}
    
    # 1. Enhanced precipitation (%)
    total_p_seeded = np.sum(precip_seeded)
    total_p_control = np.sum(precip_control)
    if total_p_control > 0:
        metrics['enhanced_precipitation_pct'] = ((total_p_seeded - total_p_control) / total_p_control) * 100.0
    else:
        metrics['enhanced_precipitation_pct'] = 0.0
        
    # 2. Ice Enhancement Ratio (IER)
    avg_ni_seeded = np.mean(ni_seeded)
    avg_ni_control = np.mean(ni_control)
    if avg_ni_control > 0:
        metrics['ier'] = avg_ni_seeded / avg_ni_control
    else:
        metrics['ier'] = 1.0
        
    # 3. Reagent Utilization (eta)
    if reagent_mass > 0:
        metrics['utilization_eta'] = (total_p_seeded - total_p_control) / reagent_mass
    else:
        metrics['utilization_eta'] = 0.0
        
    return metrics
