import numpy as np
from typing import Tuple
from cms.config import PhysicsConfig

class IceMicrophysics:
    """
    Double-moment ice microphysics including AgI interaction.
    Ref: IMPLEMENTATION_GUIDE.md Section 3.5.2 - 3.5.5
    """
    def __init__(self, config: PhysicsConfig):
        self.config = config

    def compute_nucleation_agi(self, T: np.ndarray, N_agi: np.ndarray, N_inp_nat: np.ndarray) -> np.ndarray:
        """
        Primary heterogeneous nucleation by AgI with competition.
        Ref: Eq 3.5.2 & 3.6
        """
        # INF(T)
        T_c = T - 273.15
        inf = np.zeros_like(T)
        mask_nuc = T_c < -5.0
        
        # 0.0007 * exp(0.28 * (T + 15))
        inf[mask_nuc] = 0.0007 * np.exp(0.28 * (T_c[mask_nuc] + 15.0))
        
        # Saturate INF at 1.0 (100%)
        inf = np.minimum(inf, 1.0)

        # Assuming optimal size range (50-200nm) -> f_size = 1.0 for simplified bulk AgI
        f_size = 1.0
        
        N_inp_agi = N_agi * inf * f_size
        
        # Competition Effect (Eq 3.6)
        # beta = N_nat / (N_nat + 10^4)
        beta = N_inp_nat / (N_inp_nat + 1.0e4)
        
        # Total INP available to form ice
        # N_total = N_nat + N_agi * (1 - beta)
        # But this function returns the Target N_i from AgI specifically? 
        # Or total target?
        # The equation in guide says N_i,total.
        # Let's return the Total Target Ni.
        
        N_target = N_inp_nat + N_inp_agi * (1.0 - beta)
        
        return N_target

    def compute_secondary_production(self, 
                                   T: np.ndarray, 
                                   qc: np.ndarray, qs: np.ndarray, qi: np.ndarray,
                                   ni: np.ndarray, ns: np.ndarray) -> np.ndarray:
        """
        Secondary Ice Production (SIP).
        Ref: Eq 3.5.3
        """
        T_c = T - 273.15
        dni_dt = np.zeros_like(T)
        
        # 1. Hallett-Mossop (-3 to -8 C)
        # dNi/dt_HM = C_HM * R_rime * f(T)
        # Simplified R_rime ~ qc * qs assumption for prototype
        c_hm = 3.5e8
        mask_hm = (T_c > -8.0) & (T_c < -3.0)
        
        if np.any(mask_hm):
            f_t = (T_c[mask_hm] + 8.0) / 5.0
            # Rough approx of riming rate for SIP
            r_rime = 1.0 * qc[mask_hm] * qs[mask_hm] 
            dni_dt[mask_hm] += c_hm * r_rime * f_t

        # 2. Collisional Breakup (Ice-Ice)
        # dNi/dt_coll = C_coll * Ni * Ns * |dv| * phi(T)
        c_coll = 9e-5
        mask_coll = (T_c > -27.0) & (T_c < -3.0)
        
        if np.any(mask_coll):
            phi = 0.5
            dv = 1.0 # Placeholder for relative velocity diff
            dni_dt[mask_coll] += c_coll * ni[mask_coll] * ns[mask_coll] * dv * phi
            
        return dni_dt

    def _compute_mean_diameter(self, qx: np.ndarray, nx: np.ndarray, rho_x: float) -> np.ndarray:
        """Computes mean diameter from mass and number concentration."""
        # Avoid division by zero
        nx_safe = np.maximum(nx, 1e-12)
        qx_safe = np.maximum(qx, 0.0)
        
        # Mean mass per particle
        m_x = qx_safe / nx_safe
        
        # Diameter from mass, assuming spherical particles
        D_x = (6.0 * m_x / (np.pi * rho_x))**(1.0/3.0)
        return D_x

    def compute_rates(self,
                      T: np.ndarray, rho: np.ndarray,
                      qc: np.ndarray, qr: np.ndarray, nc: np.ndarray,
                      qi: np.ndarray, qs: np.ndarray, qg: np.ndarray,
                      ni: np.ndarray, ns: np.ndarray, ng: np.ndarray,
                      N_agi: np.ndarray, N_inp_nat: np.ndarray) -> Tuple[np.ndarray, ...]:
        """
        Основной интерфейс для расчета скоростей изменения в ледяной микрофизике.
        Возвращает тенденции для всех гидрометеоров.
        """
        T_c = T - 273.15
        # Ограничиваем температуру для стабильности расчетов.
        T_c_clipped = np.clip(T_c, -80.0, 20.0)

        dqc = np.zeros_like(qc)
        dqr = np.zeros_like(qr)
        dqi = np.zeros_like(qi)
        dqs = np.zeros_like(qs)
        dqg = np.zeros_like(qg)
        dni = np.zeros_like(ni)
        dns = np.zeros_like(ns)
        dng = np.zeros_like(ng)

        # 1. Первичная нуклеация (AgI + Естественные INP) -> Ледяные кристаллы
        n_inp_target = self.compute_nucleation_agi(T, N_agi, N_inp_nat)
        nucleation_rate = np.maximum(n_inp_target - ni, 0) / 10.0  # Релаксация к целевому значению за 10с
        dni += nucleation_rate
        dqi += nucleation_rate * 1e-12  # Добавляем массу для новых кристаллов

        # 2. Вторичное образование льда
        sip_rate = self.compute_secondary_production(T, qc, qs, qi, ni, ns)
        dni += sip_rate
        dqi += sip_rate * 1e-12

        # 3. Намерзание (сбор облачной воды льдом/градом)
        # Ref: Eq 3.5.4
        E_ci = np.exp(-0.09 * (T_c_clipped + 10.0)) # Экспонента с затуханием при высоких T
        mask_rime = (T < 273.15) & (qc > 1e-9) & (qi > 1e-9)
        if np.any(mask_rime):
            D_i = self._compute_mean_diameter(qi[mask_rime], ni[mask_rime], self.config.rho_i)
            # Упрощенная разница скоростей (кристаллы < капли)
            v_diff = 1.0  # м/с
            rate_rime = (np.pi/4.0) * E_ci[mask_rime] * ni[mask_rime] * qc[mask_rime] * D_i**2 * v_diff
            
            # Конвертируем облачную воду в град
            transfer_mass = np.minimum(rate_rime, qc[mask_rime] / 1.0) # Ограничиваем скоростью исчерпания qc
            dqg[mask_rime] += transfer_mass
            dqc[mask_rime] -= transfer_mass
            #TODO: учесть изменение ni -> ng

        # 4. Агрегация (сбор ледяных кристаллов снегом)
        # Ref: Eq 3.5.4
        E_is = 0.1 * np.exp(0.025 * T_c_clipped)  # Эффективность зависит от температуры
        mask_agg = (T < 273.15) & (qi > 1e-9) & (qs > 1e-9)
        if np.any(mask_agg):
            D_i = self._compute_mean_diameter(qi[mask_agg], ni[mask_agg], self.config.rho_i)
            D_s = self._compute_mean_diameter(qs[mask_agg], ns[mask_agg], self.config.rho_i)
            v_diff_is = 0.5  # Упрощенная разница скоростей
            rate_agg = (np.pi/4.0) * E_is[mask_agg] * (ni[mask_agg]*D_i**2 + ns[mask_agg]*D_s**2) * qi[mask_agg] * v_diff_is
            
            # Конвертируем лед в снег
            transfer_mass_agg = np.minimum(rate_agg, qi[mask_agg] / 1.0)
            dqs[mask_agg] += transfer_mass_agg
            dqi[mask_agg] -= transfer_mass_agg
            #TODO: учесть изменение ni, ns

        # 5. Таяние (T > 0C)
        mask_melt = (T > 273.15)
        if np.any(mask_melt):
            # Простая релаксация к таянию
            melt_rate_i = qi[mask_melt] / 20.0 # Время таяния 20с
            melt_rate_s = qs[mask_melt] / 20.0
            melt_rate_g = qg[mask_melt] / 20.0
            
            total_melt_to_rain = melt_rate_i + melt_rate_s + melt_rate_g
            dqi[mask_melt] -= melt_rate_i
            dqs[mask_melt] -= melt_rate_s
            dqg[mask_melt] -= melt_rate_g
            dqr[mask_melt] += total_melt_to_rain
            #TODO: учесть изменение ni, ns, ng -> nr

        # 6. Гомогенное замерзание (T < -40C / 233.15K)
        mask_fz = (T < 233.15)
        if np.any(mask_fz):
            # Мгновенная конвертация всей облачной воды в лед
            fz_rate = qc[mask_fz] / 60.0  # Релаксация за 60с для стабильности
            dqc[mask_fz] -= fz_rate
            dqi[mask_fz] += fz_rate
            dni[mask_fz] += nc[mask_fz] # Каждая капля становится кристаллом
            # dnc должен быть обнулен в warm.py или здесь

        return dqc, dqr, dqi, dqs, dqg, dni, dns, dng
