import numpy as np
from typing import Tuple
from cms.config import GridConfig

class Grid:
    """
    Класс, описывающий 3D эйлерову сетку для модели CMS.
    
    Отвечает за хранение параметров сетки (размеры, шаг), создание полей
    и предоставление удобных интерфейсов для работы с основной расчетной
    областью и "призрачными" ячейками (ghost cells), которые используются
    для корректной постановки граничных условий.
    
    Ref: IMPLEMENTATION_GUIDE.md, Section 1.2 & 7.2
    """
    def __init__(self, config: GridConfig, ghost_cells: int = 3):
        """
        Инициализирует сетку.
        
        Args:
            config: Конфигурация сетки с размерами и шагами.
            ghost_cells: Количество призрачных ячеек с каждой стороны оси.
                         Необходимо для численных схем высоких порядков (например, WENO5).
        """
        # Размеры внутренней (физической) области
        self.nx = config.nx
        self.ny = config.ny
        self.nz = config.nz
        
        # Шаги сетки по осям
        self.dx = config.dx
        self.dy = config.dy
        self.dz = config.dz
        
        # Количество призрачных ячеек
        self.ng = ghost_cells

        # Полные размеры сетки, включая призрачные ячейки
        self.gnx = self.nx + 2 * self.ng
        self.gny = self.ny + 2 * self.ng
        self.gnz = self.nz + 2 * self.ng

        # Координатные массивы для внутренней области
        self.x = np.linspace(0, (self.nx - 1) * self.dx, self.nx)
        self.y = np.linspace(0, (self.ny - 1) * self.dy, self.ny)
        self.z = np.linspace(0, (self.nz - 1) * self.dz, self.nz)

        # Срезы (slices) для удобного доступа к внутренней области поля
        self.inner = (
            slice(self.ng, -self.ng),
            slice(self.ng, -self.ng),
            slice(self.ng, -self.ng)
        )

        # Кэш для 3D координатных сеток (meshgrid)
        self._X = None
        self._Y = None
        self._Z = None

    def create_field(self) -> np.ndarray:
        """
        Создает новый 3D массив (поле) с размерами полной сетки (включая призрачные ячейки),
        заполненный нулями и с нужным типом данных.
        """
        return np.zeros((self.gnx, self.gny, self.gnz), dtype=np.float64)

    def get_inner(self, field: np.ndarray) -> np.ndarray:
        """Возвращает срез поля, соответствующий внутренней расчетной области (без призрачных ячеек)."""
        return field[self.inner]

    @property
    def shape(self) -> Tuple[int, int, int]:
        """Возвращает полные размеры массивов, включая призрачные ячейки."""
        return (self.gnx, self.gny, self.gnz)

    @property
    def gx(self) -> np.ndarray:
        """Возвращает 1D массив координат по оси X для полной сетки."""
        return np.linspace(
            -(self.ng * self.dx),
            (self.nx - 1 + self.ng) * self.dx,
            self.gnx
        )

    @property
    def gy(self) -> np.ndarray:
        """Возвращает 1D массив координат по оси Y для полной сетки."""
        return np.linspace(
            -(self.ng * self.dy),
            (self.ny - 1 + self.ng) * self.dy,
            self.gny
        )

    @property
    def gz(self) -> np.ndarray:
        """Возвращает 1D массив координат по оси Z для полной сетки."""
        return np.linspace(
            -(self.ng * self.dz),
            (self.nz - 1 + self.ng) * self.dz,
            self.gnz
        )

    def _create_meshgrid(self):
        """Создает и кэширует 3D координатные сетки."""
        if self._X is None:
            self._X, self._Y, self._Z = np.meshgrid(self.gx, self.gy, self.gz, indexing='ij')

    @property
    def X(self) -> np.ndarray:
        """Возвращает 3D координатную сетку для X."""
        self._create_meshgrid()
        return self._X

    @property
    def Y(self) -> np.ndarray:
        """Возвращает 3D координатную сетку для Y."""
        self._create_meshgrid()
        return self._Y

    @property
    def Z(self) -> np.ndarray:
        """Возвращает 3D координатную сетку для Z."""
        self._create_meshgrid()
        return self._Z