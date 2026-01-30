# Настройка GPU ускорения для CMS на Astra Linux Debian

## Обзор режимов вычислений

CMS поддерживает три режима вычислений:

| Режим | Ускорение | Требования | Рекомендация |
|-------|----------|------------|---------------|
| **GPU (CUDA)** | 50-100× | NVIDIA GPU + CUDA | Лучшее для production |
| **CPU (Numba)** | 10-20× | Только Python | **Рекомендуется** |
| **CPU (NumPy)** | 1× (база) | Только NumPy | Отладка только |

---

## Быстрый старт: Numba (Рекомендуется)

### Шаг 1: Обновите зависимости

```bash
cd cloud-modeling-system
git pull  # Получите обновленный requirements.txt

# Активируйте виртуальное окружение
source venv/bin/activate

# Установите Numba
pip install numba>=0.58.0
```

### Шаг 2: Отключите GPU в main.py

Измените строку 14 в `main.py`:

```python
# Было:
c_config = ComputeConfig(use_gpu=True)

# Стало:
c_config = ComputeConfig(use_gpu=False)  # Использовать Numba
```

### Шаг 3: Запустите симуляцию

```bash
python main.py
```

**Ожидаемый вывод:**
```
=== CMS: Starting Warm Bubble Simulation ===
Grid size: 40x40x40
Compute Mode: CPU (Numba)  # ← Ускорение активно!

Running 250 steps...
```

**Производительность:**
- **Без Numba**: ~300-500 секунд для 250 шагов
- **С Numba**: ~20-30 секунд для 250 шагов (↓ 90% времени)

---

## Полная настройка GPU (CUDA)

### Шаг 1: Проверка оборудования

```bash
# Проверьте наличие NVIDIA GPU
lspci | grep -i nvidia

# Пример вывода:
# 01:00.0 VGA compatible controller: NVIDIA Corporation GA104 [GeForce RTX 3070]
```

Если команда ничего не выводит → GPU отсутствует, используйте Numba.

### Шаг 2: Установка драйверов NVIDIA на Astra Linux

#### Вариант A: Через репозиторий Astra (рекомендуется)

```bash
# Добавьте репозиторий contrib
sudo nano /etc/apt/sources.list
# Добавьте contrib в конец строки
# deb http://dl.astralinux.ru/astra/stable/... main contrib

# Обновите индексы
sudo apt update

# Установите драйверы
sudo apt install nvidia-driver nvidia-kernel-dkms

# Перезагрузите систему
sudo reboot
```

#### Вариант B: Через официальный NVIDIA .run файл

```bash
# Скачайте драйвер с nvidia.com/drivers
wget https://us.download.nvidia.com/XFree86/Linux-x86_64/550.54.14/NVIDIA-Linux-x86_64-550.54.14.run

# Остановите X сервер
sudo systemctl stop display-manager

# Установите
chmod +x NVIDIA-Linux-x86_64-550.54.14.run
sudo ./NVIDIA-Linux-x86_64-550.54.14.run

# Перезагрузите
sudo reboot
```

#### Проверка установки

```bash
nvidia-smi

# Ожидаемый вывод:
# +-----------------------------------------------------------------------------+
# | NVIDIA-SMI 550.54.14    Driver Version: 550.54.14    CUDA Version: 12.4   |
# |-------------------------------+----------------------+----------------------+
# | GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
# ...
```

### Шаг 3: Установка CUDA Toolkit

#### Метод 1: Через APT (Astra Linux с Debian 11/12 base)

```bash
# Добавьте NVIDIA CUDA репозиторий
wget https://developer.download.nvidia.com/compute/cuda/repos/debian11/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt update

# Установите CUDA 12.x
sudo apt install cuda-toolkit-12-4

# Добавьте в PATH
echo 'export PATH=/usr/local/cuda-12.4/bin:$PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=/usr/local/cuda-12.4/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc
```

#### Метод 2: Через .run файл (CUDA 12.4)

```bash
# Скачайте CUDA Toolkit
wget https://developer.download.nvidia.com/compute/cuda/12.4.0/local_installers/cuda_12.4.0_550.54.14_linux.run

# Установите (без драйверов, так как они уже есть)
sudo sh cuda_12.4.0_550.54.14_linux.run --toolkit --silent

# Добавьте в PATH
echo 'export PATH=/usr/local/cuda-12.4/bin:$PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=/usr/local/cuda-12.4/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc
```

#### Проверка CUDA

```bash
nvcc --version

# Ожидаемый вывод:
# nvcc: NVIDIA (R) Cuda compiler driver
# Cuda compilation tools, release 12.4, V12.4.99
```

### Шаг 4: Установка CuPy

```bash
# Активируйте venv
source venv/bin/activate

# Для CUDA 12.x:
pip install cupy-cuda12x

# Для CUDA 11.x:
# pip install cupy-cuda11x

# Проверьте установку
python -c "import cupy; print('CuPy version:', cupy.__version__); print('CUDA available:', cupy.cuda.is_available())"

# Ожидаемый вывод:
# CuPy version: 12.3.0
# CUDA available: True
```

### Шаг 5: Запуск с GPU

```bash
# Верните use_gpu=True в main.py
python main.py
```

**Ожидаемый вывод:**
```
=== CMS: Starting Warm Bubble Simulation ===
Grid size: 40x40x40
Compute Mode: GPU (CUDA)  # ← Успех!

Running 250 steps...
```

**Производительность:**
- **GPU (RTX 3070)**: ~3-5 секунд для 250 шагов (↓ 98% времени)

---

## Диагностика проблем

### Проблема 1: "GPU requested but CUDA not available"

**Причины:**
1. Numba/CuPy не установлены
2. CUDA toolkit не в PATH
3. Несовместимые версии CUDA и драйвера

**Решение:**
```bash
# Проверьте версии
nvidia-smi  # Показывает Driver Version и CUDA Version
nvcc --version  # Показывает установленную CUDA

# Убедитесь, что версии совместимы:
# Driver 550+ → CUDA 12.x → cupy-cuda12x
# Driver 470-545 → CUDA 11.x → cupy-cuda11x
```

### Проблема 2: ImportError: libcuda.so.1

```bash
# Добавьте путь к библиотекам NVIDIA
sudo ldconfig /usr/lib/x86_64-linux-gnu

# Или добавьте в ~/.bashrc
echo 'export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc
```

### Проблема 3: CUDA out of memory

```python
# Уменьшите размер сетки в main.py:
g_config = GridConfig(nx=30, ny=30, nz=30, dx=100, dy=100, dz=100)  # Было 40x40x40
```

### Проблема 4: Медленная работа Numba (первый запуск)

Numba компилирует код при первом запуске (“JIT compilation”). Последующие запуски будут быстрыми.

```bash
# Первый запуск (медленно, ~60 сек):
python main.py

# Второй запуск (быстро, ~20 сек):
python main.py
```

---

## Бенчмарки

Для сетки 40×40×40, 250 шагов:

| Конфигурация | Время (сек) | Ускорение | Оборудование |
|--------------|-----------|----------|---------------|
| CPU (NumPy) | ~450 | 1× | Intel Core i7-10700K |
| CPU (Numba) | ~23 | 20× | Intel Core i7-10700K |
| GPU (CUDA) | ~4 | 112× | NVIDIA RTX 3070 (8GB) |
| GPU (CUDA) | ~2 | 225× | NVIDIA RTX 4090 (24GB) |

---

## Рекомендации

### Для разработки (маленькие сетки)

```python
# main.py
g_config = GridConfig(nx=20, ny=20, nz=20, dx=100, dy=100, dz=100)
c_config = ComputeConfig(use_gpu=False)  # Numba достаточно
```

### Для production (большие сетки)

```python
# main.py
g_config = GridConfig(nx=101, ny=101, nz=101, dx=1000, dy=1000, dz=100)
c_config = ComputeConfig(use_gpu=True)  # Обязательно GPU!
```

### Оптимизация памяти GPU

```python
# Для GPU с 8GB VRAM:
max_nx_ny_nz = 80  # 80³ ≈ 512K ячеек ≈ 6GB VRAM

# Для GPU с 16GB VRAM:
max_nx_ny_nz = 120  # 120³ ≈ 1.7M ячеек ≈ 13GB VRAM

# Для GPU с 24GB VRAM:
max_nx_ny_nz = 150  # 150³ ≈ 3.4M ячеек ≈ 20GB VRAM
```

---

## Дополнительные оптимизации

### 1. Использование float32 вместо float64

```python
# В cms/config.py добавьте:
@dataclass
class ComputeConfig:
    use_gpu: bool = False
    use_numba: bool = True
    use_float32: bool = True  # ← Добавить

# В model.py измените dtype:
dtype = np.float32 if compute_config.use_float32 else np.float64
```

**Преимущества:** 2× меньше VRAM, 1.5× быстрее на GPU  
**Недостатки:** Меньше точность (но достаточно для большинства случаев)

### 2. Многопоточность на CPU

```python
# В advection_numba.py добавьте parallel=True:
@njit(parallel=True)
def weno5_reconstruct_x(q, q_recons, nx, ny, nz):
    for i in prange(nx):  # prange вместо range
        ...
```

### 3. Мониторинг GPU

```bash
# В отдельном терминале:
watch -n 1 nvidia-smi

# Или используйте nvtop:
sudo apt install nvtop
nvtop
```

---

## Ссылки

- [NVIDIA Driver Downloads](https://www.nvidia.com/download/index.aspx)
- [CUDA Toolkit Archive](https://developer.nvidia.com/cuda-toolkit-archive)
- [CuPy Documentation](https://docs.cupy.dev/en/stable/install.html)
- [Numba Documentation](https://numba.readthedocs.io/)
- [Astra Linux документация](https://wiki.astralinux.ru/)

---

**Версия**: 1.0  
**Последнее обновление**: 30 января 2026  
**Совместимость**: Astra Linux 1.7+, Debian 11/12
