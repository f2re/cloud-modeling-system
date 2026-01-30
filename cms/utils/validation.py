import numpy as np
from typing import Dict, Optional
from cms.utils.diagnostics import compute_radar_reflectivity
from cms.config import PhysicsConfig

def compute_statistics(obs: np.ndarray, model: np.ndarray) -> Dict[str, float]:
    """
    Computes RMSE, Bias, and Correlation between observations and model.
    Ref: IMPLEMENTATION_GUIDE.md Section 9.2
    """
    # Flatten and remove NaNs
    obs_f = obs.flatten()
    mod_f = model.flatten()
    mask = ~np.isnan(obs_f) & ~np.isnan(mod_f)
    
    o = obs_f[mask]
    m = mod_f[mask]
    
    if len(o) == 0:
        return {"rmse": np.nan, "bias": np.nan, "corr": np.nan}
        
    rmse = np.sqrt(np.mean((o - m)**2))
    bias = np.mean(m - o)
    
    std_o = np.std(o)
    std_m = np.std(m)
    if std_o > 1e-10 and std_m > 1e-10:
        corr = np.mean((o - np.mean(o)) * (m - np.mean(m))) / (std_o * std_m)
    else:
        corr = 0.0
        
    return {
        "rmse": rmse,
        "bias": bias,
        "correlation": corr
    }

def compute_lwp(qc: np.ndarray, qr: np.ndarray, rho: np.ndarray, dz: float) -> np.ndarray:
    """
    Computes Liquid Water Path (kg/m^2).
    Assumes vertical axis is the last dimension (2).
    """
    return np.sum((qc + qr) * rho, axis=2) * dz

def compute_iwp(qi: np.ndarray, qs: np.ndarray, qg: np.ndarray, rho: np.ndarray, dz: float) -> np.ndarray:
    """
    Computes Ice Water Path (kg/m^2).
    Assumes vertical axis is the last dimension (2).
    """
    return np.sum((qi + qs + qg) * rho, axis=2) * dz

def run_validation_suite(model_output: Dict[str, np.ndarray], 
                         observations: Dict[str, np.ndarray],
                         dz: float,
                         config: Optional[PhysicsConfig] = None) -> Dict[str, Dict[str, float]]:
    """
    Runs a full comparison between model and observations for multiple sensors.
    
    Expected model_output keys: 'qc', 'qr', 'qi', 'qs', 'qg', 'nc', 'nr', 'ni', 'ns', 'ng', 'rho'
    Expected observations keys: 'Z' (radar), 'LWP' (radiometer), 'IWP' (radiometer/sat)
    """
    results = {}
    
    # 1. Radar Reflectivity
    if 'Z' in observations:
        z_model = compute_radar_reflectivity(
            model_output['rho'],
            model_output['qr'], model_output['nr'],
            model_output.get('qi'), model_output.get('ni'),
            model_output.get('qs'), model_output.get('ns'),
            model_output.get('qg'), model_output.get('ng'),
            config=config
        )
        results['radar_Z'] = compute_statistics(observations['Z'], z_model)
        
    # 2. Liquid Water Path
    if 'LWP' in observations:
        lwp_model = compute_lwp(model_output['qc'], model_output['qr'], model_output['rho'], dz)
        results['radiometer_LWP'] = compute_statistics(observations['LWP'], lwp_model)
        
    # 3. Ice Water Path
    if 'IWP' in observations:
        iwp_model = compute_iwp(model_output['qi'], model_output['qs'], model_output['qg'], model_output['rho'], dz)
        results['radiometer_IWP'] = compute_statistics(observations['IWP'], iwp_model)
        
    return results
