import unittest
import numpy as np
from cms.dispersion.seeddisp import SeedDispModel
from cms.config import GridConfig, PhysicsConfig

class TestSeedDisp(unittest.TestCase):
    def test_seeddisp_concentration(self):
        """
        Тест на основе экспериментальных данных 2015-2017 гг.
        (Глава 6.1, монография).
        
        Проверяет порядок величины максимальной концентрации после
        моделирования выброса.
        """
        # 1. Параметры из монографии и настройки симуляции
        dose_g = 220  # граммов AgI
        height_m = 500  # метров
        wind_speed_ms = 5  # м/с
        
        # Конвертируем дозу в кг
        dose_kg = dose_g / 1000.0
        
        # Предполагаем, что выброс происходит в течение 10 секунд
        injection_duration_s = 10.0
        Q_kgs = dose_kg / injection_duration_s
        
        # Ожидаемые значения концентрации
        expected_conc_min_m3 = 1e4  # частиц/м^3
        expected_conc_max_m3 = 1e5  # частиц/м^3

        # 2. Настройка модели
        # Используем небольшую сетку для скорости теста
        g_config = GridConfig(nx=50, ny=50, nz=50, dx=100, dy=100, dz=50)
        p_config = PhysicsConfig()
        model = SeedDispModel(g_config, p_config)

        # 3. Создание начальных условий (поля ветра)
        u = np.full(model.grid.shape, wind_speed_ms)
        v = np.zeros(model.grid.shape)
        w = np.zeros(model.grid.shape)
        
        # 4. Параметры источника
        # Выброс в центре горизонтальной области
        source_params = {
            "x_source": g_config.dx * g_config.nx / 2,
            "y_source": g_config.dy * g_config.ny / 2,
            "h_source": height_m,
            "Q": Q_kgs,
            "sigma_x": 300, # Настраиваемые параметры для соответствия
            "sigma_y": 300,
            "sigma_z": 75,
        }

        # 5. Запуск симуляции
        dt_s = 1.0
        sim_duration_s = 60 # Симулируем 1 минуту
        num_steps = int(sim_duration_s / dt_s)

        for i in range(num_steps):
            # Включаем источник только в начале
            current_source_params = source_params if i * dt_s < injection_duration_s else None
            model.step(dt_s, u, v, w, source_params=current_source_params)

        # 6. Проверка результата
        # Концентрация в модели - это массовая концентрация (кг/м^3)
        max_mass_conc_kg_m3 = np.max(model.c_reagent)
        
        # Конвертируем в концентрацию частиц (частиц/м^3)
        mass_per_particle_kg = p_config.rho_agi * (4/3) * np.pi * p_config.r_agi**3
        max_particle_conc_m3 = max_mass_conc_kg_m3 / mass_per_particle_kg
        
        print(f"\n[TestSeedDisp] Max particle concentration reached: {max_particle_conc_m3:.2e} particles/m^3")
        
        # Проверяем, что результат находится в ожидаемом диапазоне
        self.assertTrue(
            expected_conc_min_m3 <= max_particle_conc_m3 <= expected_conc_max_m3,
            f"Max concentration {max_particle_conc_m3:.2e} is outside the expected range [{expected_conc_min_m3:.2e}, {expected_conc_max_m3:.2e}]"
        )

if __name__ == '__main__':
    unittest.main()
