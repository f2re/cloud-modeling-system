import numpy as np
from typing import Tuple
from cms.config import PhysicsConfig

class IceMicrophysics:
    """
    Модуль, описывающий микрофизические процессы с участием ледяной фазы.
    Включает гетерогенную нуклеацию (в том числе под действием AgI),
    вторичную генерацию льда, намерзание, агрегацию, таяние и гомогенное замерзание.
    
    Ref: IMPLEMENTATION_GUIDE.md, Section 3.5.2 - 3.5.5
    """
    def __init__(self, config: PhysicsConfig):
        """Инициализирует модуль с физическими константами."""
        self.config = config

    def compute_nucleation_agi(self, T: np.ndarray, N_agi: np.ndarray, N_inp_nat: np.ndarray) -> np.ndarray:
        """
        Расчет первичной гетерогенной нуклеации на частицах AgI с учетом
        конкуренции с естественными ядрами замерзания (INP).
        
        Args:
            T: Температура (K).
            N_agi: Концентрация частиц AgI (1/м^3).
            N_inp_nat: Концентрация естественных ядер замерзания (1/м^3).

        Returns:
            Целевая концентрация ледяных кристаллов (1/м^3), которая должна образоваться.
            
        Ref: IMPLEMENTATION_GUIDE.md, Eq 3.5.2 & 3.6
        """
        T_c = T - 273.15
        
        # Расчет доли активных ядер AgI в зависимости от температуры (Ice Nucleation Fraction).
        inf = np.zeros_like(T)
        mask_nuc = T_c < -5.0 # Нуклеация на AgI активна при T < -5 C
        inf[mask_nuc] = 0.0007 * np.exp(0.28 * (T_c[mask_nuc] + 15.0))
        inf = np.minimum(inf, 1.0) # Доля не может быть больше 1

        # Упрощенно считаем, что все частицы AgI имеют оптимальный размер (f_size = 1.0)
        f_size = 1.0
        N_inp_agi = N_agi * inf * f_size
        
        # Эффект конкуренции: наличие естественных ядер подавляет эффективность AgI
        # Ref: IMPLEMENTATION_GUIDE.md, Eq 3.6
        beta = N_inp_nat / (N_inp_nat + 1.0e4)
        N_target = N_inp_nat + N_inp_agi * (1.0 - beta)
        
        return N_target

    def compute_secondary_production(self, 
                                   T: np.ndarray, 
                                   qc: np.ndarray, qs: np.ndarray, qi: np.ndarray,
                                   ni: np.ndarray, ns: np.ndarray) -> np.ndarray:
        """
        Расчет вторичной генерации льда (Secondary Ice Production, SIP).
        Включает процессы Hallett-Mossop и дробление кристаллов при столкновениях.
        
        Ref: IMPLEMENTATION_GUIDE.md, Eq 3.5.3
        """
        T_c = T - 273.15
        dni_dt = np.zeros_like(T)
        
        # 1. Процесс Hallett-Mossop (активен в диапазоне -3..-8 C)
        c_hm = 3.5e8
        mask_hm = (T_c > -8.0) & (T_c < -3.0)
        if np.any(mask_hm):
            f_t = (T_c[mask_hm] + 8.0) / 5.0
            # Упрощенная аппроксимация скорости намерзания (riming rate)
            r_rime = 1.0 * qc[mask_hm] * qs[mask_hm] 
            dni_dt[mask_hm] += c_hm * r_rime * f_t

        # 2. Дробление при столкновениях ледяных частиц
        c_coll = 9e-5
        mask_coll = (T_c > -27.0) & (T_c < -3.0)
        if np.any(mask_coll):
            phi = 0.5
            dv = 1.0 # Заглушка для разности скоростей
            dni_dt[mask_coll] += c_coll * ni[mask_coll] * ns[mask_coll] * dv * phi
            
        return dni_dt

    def _compute_mean_diameter(self, qx: np.ndarray, nx: np.ndarray, rho_x: float) -> np.ndarray:
        """Вспомогательная функция для расчета среднего диаметра частиц."""
        nx_safe = np.maximum(nx, 1e-12)
        qx_safe = np.maximum(qx, 0.0)
        
        m_x = qx_safe / nx_safe # Средняя масса одной частицы
        D_x = (6.0 * m_x / (np.pi * rho_x))**(1.0/3.0) # Диаметр из массы для сферы
        return D_x

    def compute_rates(self,
                      T: np.ndarray, rho: np.ndarray,
                      qc: np.ndarray, qr: np.ndarray, nc: np.ndarray,
                      qi: np.ndarray, qs: np.ndarray, qg: np.ndarray,
                      ni: np.ndarray, ns: np.ndarray, ng: np.ndarray,
                      N_agi: np.ndarray, N_inp_nat: np.ndarray) -> Tuple[np.ndarray, ...]:
        """
        Основная функция, вычисляющая скорости всех процессов ледяной микрофизики.
        """
        T_c = T - 273.15
        T_c_clipped = np.clip(T_c, -80.0, 20.0) # Отсечка для стабильности

        # Инициализация массивов тенденций
        dqc, dqr, dqi, dqs, dqg = np.zeros_like(qc), np.zeros_like(qr), np.zeros_like(qi), np.zeros_like(qs), np.zeros_like(qg)
        dni, dns, dng = np.zeros_like(ni), np.zeros_like(ns), np.zeros_like(ng)

        # --- 1. Первичная нуклеация (AgI + Естественные INP) ---
        n_inp_target = self.compute_nucleation_agi(T, N_agi, N_inp_nat)
        nucleation_rate = np.maximum(n_inp_target - ni, 0) / 10.0  # Релаксация к цели за 10с
        dni += nucleation_rate
        dqi += nucleation_rate * 1e-12 # Добавляем начальную массу кристаллам

        # --- 2. Вторичная генерация льда ---
        sip_rate = self.compute_secondary_production(T, qc, qs, qi, ni, ns)
        dni += sip_rate
        dqi += sip_rate * 1e-12

        # --- 3. Намерзание (riming) ---
        # Сбор переохлажденных облачных капель ледяными кристаллами/градом.
        # Ref: IMPLEMENTATION_GUIDE.md, Eq 3.5.4
        E_ci = np.exp(-0.09 * (T_c_clipped + 10.0))
        mask_rime = (T < 273.15) & (qc > 1e-9) & (qi > 1e-9)
        if np.any(mask_rime):
            D_i = self._compute_mean_diameter(qi[mask_rime], ni[mask_rime], self.config.rho_i)
            v_diff = 1.0  # Упрощенная разница скоростей (заглушка)
            rate_rime = (np.pi/4.0) * E_ci[mask_rime] * ni[mask_rime] * qc[mask_rime] * D_i**2 * v_diff
            
            # Ограничиваем скорость, чтобы не "съесть" больше воды, чем есть
            transfer_mass = np.minimum(rate_rime, qc[mask_rime] / 1.0)
            dqg[mask_rime] += transfer_mass # Масса переходит в град/крупу
            dqc[mask_rime] -= transfer_mass # Масса уходит из облачных капель

        # --- 4. Агрегация ---
        # Слипание ледяных кристаллов друг с другом с образованием снега.
        # Ref: IMPLEMENTATION_GUIDE.md, Eq 3.5.4
        E_is = 0.1 * np.exp(0.025 * T_c_clipped)
        mask_agg = (T < 273.15) & (qi > 1e-9) & (qs > 1e-9)
        if np.any(mask_agg):
            D_i = self._compute_mean_diameter(qi[mask_agg], ni[mask_agg], self.config.rho_i)
            D_s = self._compute_mean_diameter(qs[mask_agg], ns[mask_agg], self.config.rho_i)
            v_diff_is = 0.5  # Упрощенная разница скоростей (заглушка)
            rate_agg = (np.pi/4.0) * E_is[mask_agg] * (ni[mask_agg]*D_i**2 + ns[mask_agg]*D_s**2) * qi[mask_agg] * v_diff_is
            
            # Ограничиваем скорость, чтобы не "съесть" больше льда, чем есть
            transfer_mass_agg = np.minimum(rate_agg, qi[mask_agg] / 1.0)
            dqs[mask_agg] += transfer_mass_agg # Масса переходит в снег
            dqi[mask_agg] -= transfer_mass_agg # Масса уходит из кристаллов

        # --- 5. Таяние (при T > 0 C) ---
        mask_melt = (T > 273.15)
        if np.any(mask_melt):
            # Простая схема релаксации с характерным временем таяния 20с
            melt_rate_i = qi[mask_melt] / 20.0
            melt_rate_s = qs[mask_melt] / 20.0
            melt_rate_g = qg[mask_melt] / 20.0
            
            total_melt_to_rain = melt_rate_i + melt_rate_s + melt_rate_g
            dqi[mask_melt] -= melt_rate_i
            dqs[mask_melt] -= melt_rate_s
            dqg[mask_melt] -= melt_rate_g
            dqr[mask_melt] += total_melt_to_rain

        # --- 6. Гомогенное замерзание (при T < -40 C) ---
        mask_fz = (T < 233.15)
        if np.any(mask_fz):
            # Быстрая, но стабильная конвертация всей облачной воды в лед
            fz_rate = qc[mask_fz] / 60.0  # Релаксация за 60с для стабильности
            dqc[mask_fz] -= fz_rate
            dqi[mask_fz] += fz_rate
            dni[mask_fz] += nc[mask_fz] # Каждая капля становится кристаллом

        return dqc, dqr, dqi, dqs, dqg, dni, dns, dng