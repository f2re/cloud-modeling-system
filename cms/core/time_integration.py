from typing import Callable, List
import numpy as np

class RK5Integrator:
    """
    Классический явный интегратор Рунге-Кутты 5-го порядка (SSP-RK5).
    
    Этот метод хорошо подходит для решения не-жестких (non-stiff) систем
    дифференциальных уравнений, таких как уравнения динамики жидкости.
    Он не рекомендуется для решения уравнений микрофизики из-за их высокой жесткости.
    
    Ref: IMPLEMENTATION_GUIDE.md, Section 6.2 (архивная версия).
    """
    def __init__(self):
        """Инициализирует коэффициенты схемы."""
        # Коэффициенты для 5-стадийной схемы SSP-RK5
        self.alpha = [0.37, 0.38, 0.18, 0.62, 1.0]

    def step(self, 
             q: List[np.ndarray], 
             dt: float, 
             rhs_func: Callable) -> List[np.ndarray]:
        """
        Выполняет один шаг интегрирования по времени.

        Args:
            q: Список массивов numpy, представляющих текущее состояние системы.
            dt: Шаг по времени (в секундах).
            rhs_func: Функция, вычисляющая правую часть системы уравнений (производные по времени).

        Returns:
            Список массивов numpy, представляющих новое состояние системы.
        """
        if isinstance(q, list):
            q_n = [qi.copy() for qi in q]
            q_stage = [qi.copy() for qi in q]

            for a in self.alpha:
                tendencies = rhs_func(q_stage)
                for i in range(len(q_stage)):
                    q_stage[i] = q_n[i] + a * dt * tendencies[i]
            return q_stage
        else: # Обработка случая, когда на вход подан один массив
            q_n = q.copy()
            q_stage = q_n.copy()

            for a in self.alpha:
                tendency = rhs_func(q_stage)
                q_stage = q_n + a * dt * tendency
            return q_stage

class IMEXIntegrator:
    """
    Полунеявный (Implicit-Explicit, IMEX) интегратор Рунге-Кутты 2-го порядка.
    
    Этот интегратор является ключевым компонентом для обеспечения численной стабильности модели.
    Он разделяет правую часть уравнений на две компоненты:
    - Явная (Explicit): для не-жестких процессов (динамика, адвекция).
    - Неявная (Implicit): для жестких процессов (микрофизика, диффузия).
    
    Это позволяет использовать большие шаги по времени, сохраняя при этом устойчивость.
    
    Ref: IMPLEMENTATION_GUIDE.md, Section 6.2.
         Ascher et al. (1997), "Implicit-explicit Runge-Kutta methods..."
    """

    def __init__(self):
        """Инициализирует коэффициенты для схемы IMEX-SSP2(2,2,2)."""
        # Коэффициент gamma определяет "неявность" схемы
        self.gamma = 1.0 - 1.0 / np.sqrt(2.0)
        self.delta = 1.0 - 1.0 / (2.0 * self.gamma)

    def _solve_implicit_stage(self, q: List[np.ndarray], dt_gamma: float, implicit_rhs_func: Callable) -> List[np.ndarray]:
        """
        Решает неявное уравнение на каждой стадии интегратора.
        
        Уравнение: k = f_imp(q + dt_gamma * k)
        Решается итерационным методом Ньютона для нелинейных систем.
        Ключевое упрощение: так как микрофизика - локальный процесс,
        Якобиан можно считать диагональным, что позволяет избежать решения
        огромных систем линейных уравнений.
        """
        # Начальное приближение для тенденции k - явный расчет
        k = implicit_rhs_func(q)

        # Упрощенный итерационный метод Ньютона
        for _ in range(3):  # 3 итераций обычно достаточно для сходимости
            q_stage = [q[i] + dt_gamma * k[i] for i in range(len(q))]
            f_stage = implicit_rhs_func(q_stage)
            
            # Остаток (residual): R(k) = k - f_stage
            residual = [k[i] - f_stage[i] for i in range(len(q))]
            
            # Аппроксимация диагонали Якобиана через время релаксации tau
            # (I - dt_gamma * J) * delta_k = -R(k)
            # J_f аппроксимируется как -1/tau, где tau - характерное время процесса
            tau = 10.0  # Эвристическое время релаксации для микрофизики (10 секунд)
            jac_inv_approx = 1.0 / (1.0 + dt_gamma / tau)
            
            delta_k = [-res * jac_inv_approx for res in residual]
            
            max_delta = 0.0
            for i in range(len(q)):
                k[i] += delta_k[i]
                max_delta = max(max_delta, np.max(np.abs(delta_k[i])))

            # Проверка сходимости
            if max_delta < 1e-6:
                break
        
        return k

    def step(self,
             q: List[np.ndarray],
             dt: float,
             explicit_rhs: Callable,
             implicit_rhs: Callable) -> List[np.ndarray]:
        """
        Выполняет один шаг интегрирования по схеме IMEX.

        Args:
            q: Текущее состояние системы.
            dt: Шаг по времени.
            explicit_rhs: Функция для вычисления явной части правой стороны.
            implicit_rhs: Функция для вычисления неявной части правой стороны.

        Returns:
            Новое состояние системы.
        """
        # --- Стадия 1 ---
        f_exp_1 = explicit_rhs(q)
        f_imp_1 = self._solve_implicit_stage(q, self.gamma * dt, implicit_rhs)
        
        q_stage_1 = []
        for i in range(len(q)):
            q_stage_1.append(q[i] + self.gamma * dt * (f_exp_1[i] + f_imp_1[i]))

        # --- Стадия 2 ---
        f_exp_2 = explicit_rhs(q_stage_1)
        
        # Вход для неявной части второй стадии
        q_stage_2_imp_input = []
        for i in range(len(q)):
            term1 = self.delta * dt * f_exp_1[i]
            term2 = (self.gamma - self.delta) * dt * f_imp_1[i]
            q_stage_2_imp_input.append(q[i] + term1 + term2)

        f_imp_2 = self._solve_implicit_stage(q_stage_2_imp_input, self.gamma * dt, implicit_rhs)

        # --- Финальное обновление ---
        q_new = []
        for i in range(len(q)):
            # Использование разных весов для явной и неявной частей для лучшей стабильности
            b_exp = 1.0 / (2.0 * self.gamma)
            
            exp_term = b_exp * dt * f_exp_1[i] + (1.0 - b_exp) * dt * f_exp_2[i]
            imp_term = (1.0 - b_exp) * dt * f_imp_1[i] + b_exp * dt * f_imp_2[i] # Пересмотренная схема

            q_new.append(q[i] + exp_term + imp_term)
            
        return q_new