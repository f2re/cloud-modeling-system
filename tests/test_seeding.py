import unittest
import numpy as np
from cms.model import CMSModel
from cms.config import GridConfig, PhysicsConfig

class TestSeedingModel(unittest.TestCase):
    def test_seeding_ice_enhancement(self):
        """
        Тест на основе данных главы 6 монографии.
        Проверка увеличения концентрации льда после засева.
        """
        # 1. Настройка сетки и физики в соответствии с монографией
        g_config = GridConfig(nx=20, ny=20, nz=40, # Уменьшенная сетка для скорости
                             dx=500, dy=500, dz=50)
        p_config = PhysicsConfig()
        model = CMSModel(g_config, p_config)
        
        # 2. Установка начальных условий
        # Создаем переохлажденное облако
        z_cloud_base = 10 # индекс высоты 500м
        z_cloud_top = 30  # индекс высоты 1500м
        
        # Температура от 0C на 0м до -10C на 2000м
        # T(z) = T0 - gamma * z
        T0 = 273.15
        gamma = 0.005 # 5 K/km
        for k in range(model.grid.nz):
            model.theta[:, :, k] = T0 - gamma * (k * model.grid.dz)
        
        # Насыщаем слой, где будет облако
        model.qv[:, :, z_cloud_base:z_cloud_top] = 7e-3 # ~7 g/kg

        # Прогоняем несколько шагов для образования облака
        for _ in range(200): # 20s total
            model.step(dt=0.1)
        
        # Убедимся, что облако сформировалось
        self.assertGreater(np.max(model.qc), 1e-5, "Облако не сформировалось перед засевом")

        # 3. Фиксируем начальную концентрацию льда
        ni_initial = np.max(model.ni)
        
        # 4. Добавляем реагент AgI в центр облака
        mid_x, mid_y = g_config.nx // 2, g_config.ny // 2
        z_seed_level = (z_cloud_base + z_cloud_top) // 2
        
        # Концентрация AgI в кг/кг
        # 1e11 частиц/м^3 * (масса одной частицы) / (плотность воздуха)
        mass_per_particle = p_config.rho_agi * (4/3) * np.pi * p_config.r_agi**3
        c_agi_conc = 1e11 * mass_per_particle / model.rho[mid_x, mid_y, z_seed_level]
        model.c_agi[mid_x, mid_y, z_seed_level] = c_agi_conc
        
        # 5. Прогоняем модель после засева
        sim_duration_minutes = 1 # Уменьшено для скорости теста
        dt_s = 0.5 # Возвращено к исходному значению
        num_steps = int(sim_duration_minutes * 60 / dt_s)
        
        for step in range(num_steps):
            model.step(dt=dt_s)
        
        # 6. Проверяем усиление ледяной фазы
        ni_final = np.max(model.ni)
        
        # Избегаем деления на ноль, если льда не было совсем
        ice_enhancement_ratio = ni_final / max(ni_initial, 1.0) # Используем 1/м3 как минимум
        
        print(f"\n[TestSeeding] Initial Ni: {ni_initial:.2e}, Final Ni: {ni_final:.2e}, Enhancement: {ice_enhancement_ratio:.2f}x")
        
        # Из главы 6: IER (Ice Enhancement Ratio) должен быть > 10
        self.assertGreater(ice_enhancement_ratio, 10, "Засев не привел к значительному увеличению концентрации льда")
        self.assertLess(ice_enhancement_ratio, 1000, "Усиление слишком велико, возможна численная неустойчивость")

if __name__ == '__main__':
    unittest.main()
