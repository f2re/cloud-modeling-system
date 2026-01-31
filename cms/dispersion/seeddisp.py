import numpy as np
from cms.core.grid import Grid
from cms.core.time_integration import RK5Integrator
from cms.dynamics.advection import WENO5
from cms.diffusion.turbulence import Turbulence
from cms.dynamics.boundary import BoundaryConditions
from cms.dynamics.sedimentation import compute_terminal_velocity
from cms.microphysics.seeding_source import compute_reagent_emission
from cms.config import GridConfig, PhysicsConfig, ComputeConfig

class SeedDispModel:
    """
    Модель переноса и рассеяния реагентов (Глава 5, монография).
    Эта модель является автономной и управляет полем концентрации реагента.
    """
    def __init__(self, grid_config: GridConfig, physics_config: PhysicsConfig, compute_config: ComputeConfig = ComputeConfig()):
        self.grid = Grid(grid_config)
        self.physics = physics_config
        self.compute_config = compute_config
        self.integrator = RK5Integrator()
        
        # Physics and Dynamics Modules
        self.advection = WENO5(self.grid, use_gpu=self.compute_config.use_gpu)
        self.diffusion = Turbulence(self.grid, self.physics)
        self.bc = BoundaryConditions(self.grid)
        
        # State variables
        self.c_reagent = self.grid.create_field()
        self.time = 0.0

    def _vertical_sedimentation(self, field: np.ndarray, v_fall: float) -> np.ndarray:
        """
        Вычисляет тенденцию от гравитационного оседания с постоянной скоростью.
        Используется простая схема против потока первого порядка.
        """
        # v_fall положительна для движения вниз (в сторону уменьшения индекса k)
        # Поток F = v_fall * field
        # d(field)/dt = -dF/dz
        flux_k_p1 = v_fall * np.roll(field, -1, axis=2) # Поток на верхней границе ячейки
        flux_k = v_fall * field                          # Поток на нижней границе ячейки
        
        # Граничное условие на верхней границе (k=nz-1): нулевой поток
        flux_k_p1[:, :, -1] = 0
        
        tendency = -(flux_k_p1 - flux_k) / self.grid.dz
        return tendency

    def _rhs(self, c_reagent: np.ndarray, u: np.ndarray, v: np.ndarray, w: np.ndarray, source_params: dict = None) -> np.ndarray:
        """
        Вычисляет полную правую часть (тенденцию) для уравнения переноса концентрации.
        """
        # 1. Адвекция фоновым ветром
        dc_advection = self.advection.advect(c_reagent, u, v, w)
        
        # 2. Турбулентная диффузия
        # nu_t должен быть вычислен извне или здесь
        nu_t = self.diffusion.compute_eddy_viscosity(u, v, w)
        dc_diffusion = self.diffusion.compute_diffusion(c_reagent, nu_t)
        
        # 3. Гравитационное оседание
        v_fall = compute_terminal_velocity(
            self.physics.r_agi, 
            self.physics.rho_agi, 
            self.physics.rho_a_stp, # Используем плотность воздуха СТП для простоты
            self.physics.eta_air
        )
        dc_sedimentation = self._vertical_sedimentation(c_reagent, v_fall)

        # 4. Источники
        dc_source = self.grid.create_field()
        if source_params:
            dc_source = compute_reagent_emission(self.grid, **source_params)

        # Суммируем все тенденции
        return dc_advection + dc_diffusion + dc_sedimentation + dc_source

    def step(self, dt: float, u: np.ndarray, v: np.ndarray, w: np.ndarray, source_params: dict = None):
        """
        Выполняет один шаг моделирования по времени.
        
        Args:
            dt (float): Шаг по времени.
            u, v, w (np.ndarray): 3D поля компонент скорости ветра.
            source_params (dict, optional): Параметры для источника реагента.
        """
        # Создаем обертку для передачи в интегратор
        def rhs_wrapper(state_vector):
            # source_params должны быть переданы в _rhs
            return [self._rhs(state_vector[0], u, v, w, source_params)]

        # Интегрируем
        new_state = self.integrator.step([self.c_reagent], dt, rhs_wrapper)
        self.c_reagent = new_state[0]
        
        # Применяем граничные условия и физические ограничения
        np.maximum(self.c_reagent, 0.0, out=self.c_reagent)
        self.bc.apply_lateral_periodic(self.c_reagent)
        
        self.time += dt
