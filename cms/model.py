import numpy as np
from cms.core.grid import Grid
from cms.core.time_integration import IMEXIntegrator
from cms.dynamics.navier_stokes import NavierStokesSolver
from cms.dynamics.boundary import BoundaryConditions
from cms.microphysics.warm import WarmMicrophysics
from cms.microphysics.ice import IceMicrophysics
from cms.microphysics.activation import Activation
from cms.microphysics.condensation import CondensationEvaporation
from cms.diffusion.turbulence import Turbulence
from cms.config import PhysicsConfig, GridConfig, ComputeConfig

class CMSModel:
    """
    Основной класс модели, который связывает воедино все компоненты:
    динамику, микрофизику, турбулентность и аэрозоли.
    
    Использует полунеявную схему интегрирования по времени (IMEX) для
    обеспечения численной стабильности при работе с жесткими (stiff)
    уравнениями микрофизики.
    """
    def __init__(self, grid_config: GridConfig, physics_config: PhysicsConfig, compute_config: ComputeConfig = ComputeConfig()):
        """
        Инициализирует модель, сетку, физические модули и состояние модели.
        
        Args:
            grid_config: Конфигурация расчетной сетки.
            physics_config: Конфигурация физических констант.
            compute_config: Конфигурация вычислительных параметров (GPU/Numba).
        """
        self.grid = Grid(grid_config)
        self.physics = physics_config
        self.compute_config = compute_config
        self.integrator = IMEXIntegrator()  # Используем стабильный IMEX интегратор
        self.time = 0.0  # Время моделирования
        
        # --- Инициализация физических модулей ---
        self.dynamics = NavierStokesSolver(self.grid, self.physics, use_gpu=self.compute_config.use_gpu)
        self.warm_micro = WarmMicrophysics(self.physics)
        self.ice_micro = IceMicrophysics(self.physics)
        self.activation = Activation(self.physics)
        self.condensation = CondensationEvaporation(self.physics)
        self.turbulence = Turbulence(self.grid, self.physics)
        self.bc = BoundaryConditions(self.grid)

        # --- Инициализация полей состояния модели (прогностические переменные) ---
        # Динамические поля
        self.u = self.grid.create_field()  # Скорость ветра по X, м/с
        self.v = self.grid.create_field()  # Скорость ветра по Y, м/с
        self.w = self.grid.create_field()  # Скорость ветра по Z, м/с
        self.rho = self.grid.create_field() + self.physics.rho_a_stp  # Плотность воздуха, кг/м^3
        self.theta = self.grid.create_field() + 300.0  # Потенциальная температура, K
        self.p = self.grid.create_field() + 101325.0 # Давление, Па
        
        # Поля гидрометеоров: массовые доли (q) и числовые концентрации (N)
        self.qv = self.grid.create_field() + 0.01  # Водяной пар, кг/кг
        self.qc = self.grid.create_field()         # Облачные капли (масса), кг/кг
        self.qr = self.grid.create_field()         # Дождевые капли (масса), кг/кг
        self.nc = self.grid.create_field()         # Облачные капли (концентрация), 1/м^3
        self.nr = self.grid.create_field()         # Дождевые капли (концентрация), 1/м^3
        
        self.qi = self.grid.create_field()         # Кристаллы льда (масса), кг/кг
        self.qs = self.grid.create_field()         # Снег (масса), кг/кг
        self.qg = self.grid.create_field()         # Град/крупа (масса), кг/кг
        self.ni = self.grid.create_field()         # Кристаллы льда (концентрация), 1/м^3
        self.ns = self.grid.create_field()         # Снег (концентрация), 1/м^3
        self.ng = self.grid.create_field()         # Град/крупа (концентрация), 1/м^3
        
        # Поля аэрозолей и реагентов
        self.c_agi = self.grid.create_field()      # Концентрация реагента AgI, кг/кг
        self.n_ccn = self.grid.create_field() + 1e8 # Ядра конденсации, 1/м^3
        self.n_inp_nat = self.grid.create_field() + 1e3 # Естественные ядра замерзания, 1/м^3
        
        # Предварительное выделение памяти для массивов тенденций
        self.num_vars = 19
        self.tendencies = [self.grid.create_field() for _ in range(self.num_vars)]

    def _unpack_state(self, state_vector: list) -> tuple:
        """Вспомогательная функция для распаковки вектора состояния в именованные переменные."""
        return tuple(state_vector)

    def explicit_rhs(self, state_vector: list) -> list:
        """
        Вычисляет тенденции для явной части системы (быстрые процессы).
        Сюда входят адвекция и динамика (градиент давления, плавучесть).
        """
        # Гарантия положительной определенности на входе от интегратора
        for field in state_vector[3:]:  # rho и все скаляры
            np.maximum(field, 0.0, out=field)
            
        u, v, w, rho, theta, qv, qc, qr, nc, nr, qi, qs, qg, ni, ns, ng, c_agi, n_ccn, n_inp_nat = self._unpack_state(state_vector)
        
        # Сброс массива тенденций
        for t in self.tendencies: t.fill(0)
        du, dv, dw, drho, dtheta, dqv, dqc, dqr, dnc, dnr, dqi, dqs, dqg, dni, dns, dng, dc_agi, dn_ccn, dn_inp_nat = self._unpack_state(self.tendencies)

        # --- Динамика (адвекция + градиент давления) ---
        du_dyn, dv_dyn, dw_dyn, drho_dyn, dtheta_dyn = self.dynamics.compute_tendencies(u, v, w, rho, theta, self.p)
        du += du_dyn; dv += dv_dyn; dw += dw_dyn; drho += drho_dyn; dtheta += dtheta_dyn

        # --- Адвекция всех скалярных полей ---
        weno = self.dynamics.weno
        scalars = [qv, qc, qr, nc, nr, qi, qs, qg, ni, ns, ng, c_agi, n_ccn, n_inp_nat]
        d_scalars = [weno.advect(s, u, v, w) for s in scalars]
        
        (dqv_adv, dqc_adv, dqr_adv, dnc_adv, dnr_adv, dqi_adv, dqs_adv, dqg_adv, 
         dni_adv, dns_adv, dng_adv, dc_agi_adv, dn_ccn_adv, dn_inp_nat_adv) = d_scalars
        
        dqv += dqv_adv; dqc += dqc_adv; dqr += dqr_adv; dnc += dnc_adv; dnr += dnr_adv
        dqi += dqi_adv; dqs += dqs_adv; dqg += dqg_adv; dni += dni_adv; dns += dns_adv; dng += dng_adv
        dc_agi += dc_agi_adv; dn_ccn += dn_ccn_adv; dn_inp_nat += dn_inp_nat_adv
        
        return self.tendencies

    def implicit_rhs(self, state_vector: list) -> list:
        """
        Вычисляет тенденции для неявной части системы (жесткие процессы).
        Сюда входят микрофизика, турбулентность и термодинамические обратные связи.
        """
        # Гарантия положительной определенности на входе
        for field in state_vector[3:]: # rho и все скаляры
            np.maximum(field, 0.0, out=field)

        u, v, w, rho, theta, qv, qc, qr, nc, nr, qi, qs, qg, ni, ns, ng, c_agi, n_ccn, n_inp_nat = self._unpack_state(state_vector)

        for t in self.tendencies: t.fill(0)
        du, dv, dw, drho, dtheta, dqv, dqc, dqr, dnc, dnr, dqi, dqs, dqg, dni, dns, dng, dc_agi, dn_ccn, dn_inp_nat = self._unpack_state(self.tendencies)
        
        # --- Термодинамика: вычисление абсолютной температуры ---
        T = theta * np.power(self.p / self.physics.p0, self.physics.rd / self.physics.cp)

        # --- Турбулентность (обрабатывается неявно) ---
        nu_t = self.turbulence.compute_eddy_viscosity(u, v, w)
        du += self.turbulence.compute_diffusion(u, nu_t)
        dv += self.turbulence.compute_diffusion(v, nu_t)
        dw += self.turbulence.compute_diffusion(w, nu_t)
        # Диффузия всех скаляров
        all_scalars = state_vector[5:]
        all_d_scalars = self.tendencies[5:]
        for i in range(len(all_scalars)):
            all_d_scalars[i] += self.turbulence.compute_diffusion(all_scalars[i], nu_t)

        # --- Микрофизика (обрабатывается неявно) ---
        # 1. Конденсация/Испарение
        dqc_cond, dqv_cond = self.condensation.compute_condensation_rate(qv, qc, T, self.p, rho)
        dqc += dqc_cond; dqv += dqv_cond
        
        # 2. Активация ядер конденсации
        dnc_act, dn_ccn_act = self.activation.compute_activation(w, n_ccn, nc)
        dnc += dnc_act; dn_ccn += dn_ccn_act
        # Увеличение массы облачных капель при активации
        dqc += dnc_act * (4/3 * np.pi * (1e-6)**3 * self.physics.rho_w)

        # 3. Физика теплого облака
        w_dqc, w_dqr, w_dnc, w_dnr = self.warm_micro.compute_rates(qc, qr, nc, nr, rho)
        dqc += w_dqc; dqr += w_dqr; dnc += w_dnc; dnr += w_dnr
        
        # 4. Физика льда
        i_dqc, i_dqr, i_dqi, i_dqs, i_dqg, i_dni, i_dns, i_dng = self.ice_micro.compute_rates(
            T, rho, qc, qr, nc, qi, qs, qg, ni, ns, ng, c_agi, n_inp_nat
        )
        dqc += i_dqc; dqr += i_dqr; dqi += i_dqi; dqs += i_dqs; dqg += i_dqg
        dni += i_dni; dns += i_dns; dng += i_dng
        
        # --- Термодинамическая обратная связь от неявных процессов ---
        dtheta -= (self.physics.l_v / self.physics.cp) * dqv_cond # от конденсации/испарения
        dtheta += (self.physics.l_f / self.physics.cp) * (i_dqi)   # от замерзания/таяния (упрощенно)
        
        return self.tendencies

    def step(self, dt: float):
        """Выполняет один шаг моделирования по времени с использованием схемы IMEX."""
        self._step_internal(dt)
        self.time += dt
        
    def _step_internal(self, dt: float):
        """Внутренняя функция-обертка для шага по времени."""
        # Упаковка всех полей модели в единый вектор состояния
        state = [
            self.u, self.v, self.w, self.rho, self.theta, self.qv, self.qc, self.qr, 
            self.nc, self.nr, self.qi, self.qs, self.qg, self.ni, self.ns, self.ng, 
            self.c_agi, self.n_ccn, self.n_inp_nat
        ]
                 
        # Вызов интегратора с разделенными явной и неявной частями
        new_state = self.integrator.step(state, dt, self.explicit_rhs, self.implicit_rhs)
        
        # Распаковка и обновление полей модели
        (self.u, self.v, self.w, self.rho, self.theta, self.qv, self.qc, self.qr, 
         self.nc, self.nr, self.qi, self.qs, self.qg, self.ni, self.ns, self.ng, 
         self.c_agi, self.n_ccn, self.n_inp_nat) = new_state
         
        # Применение гарантии положительной определенности после полного шага
        np.maximum(self.rho, 1e-9, out=self.rho)
        for field in state[5:]: # Только скалярные величины
            np.maximum(field, 0.0, out=field)
         
        # Применение граничных условий
        for field in new_state:
            self.bc.apply_lateral_periodic(field)
        self.bc.apply_bottom_noslip(self.u, self.v, self.w)