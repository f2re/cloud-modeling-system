# Cloud Seeding Numerical Model: Implementation Guide

## 1. Model Architecture Overview

### 1.1 Core Components
Implement four integrated modules:
- **SeedDisp**: Reagent dispersion in atmosphere (3D Eulerian)
- **Seeding**: Cloud microphysics with seeding agents (warm+ice physics)
- **FogSeeding**: Fog dissipation specialized module
- **Cloud**: Deep convective cloud dynamics (compressible)

### 1.2 Computational Domain
- Grid: 101×101×101 nodes (adjustable)
- Horizontal resolution: 500-2000 m
- Vertical resolution: 50-100 m
- Time step: 5-10 s (adaptive CFL)
- Integration scheme: 5-stage Runge-Kutta or 3-stage for efficiency

---

## 2. SeedDisp Module: Reagent Transport

### 2.1 Governing Equations

**Advection-diffusion for reagent concentration C (m⁻³):**

\[
\frac{\partial C}{\partial t} = -\nabla \cdot (C\vec{u}) + \nabla \cdot (K_h \nabla_h C) + \frac{\partial}{\partial z}\left(K_z \frac{\partial C}{\partial z}\right) + S_C
\]

where:
- \(\vec{u} = (u, v, w)\): wind velocity (m/s)
- \(K_h\): horizontal diffusivity (m²/s)
- \(K_z\): vertical diffusivity (m²/s)
- \(S_C\): source/sink terms (nucleation, deposition)

**Vertical velocity with buoyancy:**

\[
w = w_{\text{mean}} + w_{\text{turb}} + \frac{g}{\theta_0}(\theta - \theta_0)z
\]

### 2.2 Turbulent Diffusivity

**Stable conditions (Ri > 0):**

\[
K_z = \frac{k u_* z}{(1 + 5Ri)^2}
\]

**Unstable conditions (Ri < 0):**

\[
K_z = k u_* z \left(1 - 16\frac{z}{L}\right)^{1/2}
\]

where:
- \(k = 0.4\): von Kármán constant
- \(u_*\): friction velocity
- \(L\): Monin-Obukhov length
- \(Ri\): Richardson number

**Compute L iteratively:**

\[
L = -\frac{\theta_v u_*^3}{k g (w'\theta_v')_s}
\]

### 2.3 Surface Deposition

\[
v_d = \frac{1}{r_a + r_s}
\]

where:
- \(r_a = \frac{1}{k u_*} \ln\frac{z}{z_0}\): aerodynamic resistance
- \(r_s\): surface resistance (depends on particle size)

### 2.4 Particle Settling

**For AgI particles (50-200 nm optimal size):**

\[
v_s = \frac{2gr^2(\rho_p - \rho_a)}{9\mu} \cdot \text{Cc}(r)
\]

Cunningham slip correction:

\[
\text{Cc}(r) = 1 + \frac{\lambda}{r}\left(1.257 + 0.4e^{-1.1r/\lambda}\right)
\]

where \(\lambda = 65\) nm (mean free path at STP).

---

## 3. Seeding Module: Double-Moment Microphysics

### 3.1 Core Philosophy
**Use double-moment scheme**: predict both mass \(q_x\) (kg/kg) and number concentration \(N_x\) (m⁻³) for all hydrometeor categories.

### 3.2 Hydrometeor Categories
1. Cloud droplets (c)
2. Rain drops (r)
3. Cloud ice (i)
4. Snow (s)
5. Graupel/hail (g)
6. AgI aerosol (a)

### 3.3 Prognostic Equations

**Mass mixing ratio:**

\[
\frac{\partial q_x}{\partial t} = -\nabla \cdot (q_x \vec{u}) + \frac{\partial}{\partial z}(v_{tx} q_x) + S_{q_x}
\]

**Number concentration (CRITICAL - adds Morrison scheme):**

\[
\frac{\partial N_x}{\partial t} = -\nabla \cdot (N_x \vec{u}) + \frac{\partial}{\partial z}(v_{tn} N_x) + S_{N_x}
\]

where:
- \(v_{tx}\): mass-weighted terminal velocity
- \(v_{tn}\): number-weighted terminal velocity
- \(S_{q_x}, S_{N_x}\): microphysical source/sink terms

### 3.4 Terminal Velocity Parameterization

**General power law:**

\[
v_t(D) = a_v D^{b_v} \left(\frac{\rho_0}{\rho}\right)^{0.5}
\]

**Mass-weighted velocity:**

\[
v_{tx} = a_v \left(\frac{\rho q_x}{\pi \rho_x N_x}\right)^{b_v/3} \Gamma\left(\frac{b_v + 4}{\alpha}\right) / \Gamma\left(\frac{4}{\alpha}\right)
\]

**Number-weighted velocity:**

\[
v_{tn} = a_v \left(\frac{\rho q_x}{\pi \rho_x N_x}\right)^{b_v/3} \Gamma\left(\frac{b_v + 1}{\alpha}\right) / \Gamma\left(\frac{4}{\alpha}\right)
\]

where \(\alpha\): shape parameter of gamma distribution.

### 3.5 Microphysical Process Rates

#### 3.5.1 Warm-Rain Processes

**Autoconversion (cloud → rain):**

\[
\left(\frac{dq_r}{dt}\right)_{\text{auto}} = \frac{1350 q_c^{2.47} N_c^{-1.79}}{\rho^{1.47}}
\]

**Self-collection (rain):**

\[
\left(\frac{dN_r}{dt}\right)_{\text{self}} = -5.78 N_r^2 \bar{D}_r^3
\]

**Accretion (cloud by rain):**

\[
\left(\frac{dq_r}{dt}\right)_{\text{accr}} = \frac{\pi}{4} E_{cr} N_r q_c \bar{D}_r^2 v_r
\]

where \(E_{cr} = 1.0\) for rain-cloud collisions.

#### 3.5.2 Ice Nucleation (UPDATED 2024)

**Primary heterogeneous nucleation by AgI:**

\[
N_{\text{INP}}^{\text{AgI}} = N_{\text{AgI}} \cdot \text{INF}(T) \cdot f_{\text{size}}(d_{\text{AgI}})
\]

**Ice-nucleated fraction (Jiang et al. 2025):**

\[
\text{INF}(T) = \begin{cases}
0.0007 \cdot \exp(0.28 \cdot (T + 15)) & T < -5°C \\
0 & T \geq -5°C
\end{cases}
\]

**Size efficiency function:**

\[
f_{\text{size}}(d) = \begin{cases}
0.1 & d < 50\text{ nm} \\
1.0 & 50 \leq d \leq 200\text{ nm} \\
0.5 & d > 200\text{ nm}
\end{cases}
\]

**Depositional nucleation rate:**

\[
\frac{dN_i}{dt}\bigg|_{\text{dep}} = N_{\text{INP}}^{\text{AgI}} \cdot \max(S_i - 1, 0)
\]

where \(S_i = e/e_{si}\): ice saturation ratio.

#### 3.5.3 Secondary Ice Production (NEW 2024)

**Hallett-Mossop (-3 to -8°C):**

\[
\left(\frac{dN_i}{dt}\right)_{\text{HM}} = C_{\text{HM}} \cdot R_{\text{rime}} \cdot f(T)
\]

\[
f(T) = \begin{cases}
\frac{T + 8}{5} & -8 < T < -3 \\
0 & \text{otherwise}
\end{cases}
\]

where \(C_{\text{HM}} = 3.5 \times 10^8\) kg⁻¹, \(R_{\text{rime}}\): riming rate.

**Ice-ice collisional breakup (NEW):**

\[
\left(\frac{dN_i}{dt}\right)_{\text{coll}} = C_{\text{coll}} \cdot N_i N_s \cdot |\Delta v| \cdot \phi(T)
\]

\[
\phi(T) = \begin{cases}
0.5 & -27 < T < -3°C \\
0 & \text{otherwise}
\end{cases}
\]

where \(C_{\text{coll}} = 9 \times 10^{-5}\) m³ kg⁻¹.

**Sublimational breakup (NEW):**

\[
\left(\frac{dN_i}{dt}\right)_{\text{subl}} = C_{\text{subl}} \cdot \left|\frac{dq_i}{dt}\bigg|_{\text{subl}}\right| \cdot (1 - S_i)
\]

where \(C_{\text{subl}} = 5 \times 10^6\) kg⁻¹.

#### 3.5.4 Riming and Aggregation

**Riming (supercooled water + ice):**

\[
\left(\frac{dq_g}{dt}\right)_{\text{rime}} = \frac{\pi}{4} E_{ci} N_i q_c \bar{D}_i^2 |v_i - v_c|
\]

**Collection efficiency:**

\[
E_{ci} = \exp\left[-0.09(T + 273.15)\right]
\]

**Aggregation (ice + snow):**

\[
\left(\frac{dq_s}{dt}\right)_{\text{agg}} = \frac{\pi}{4} E_{is} (N_i \bar{D}_i^2 + N_s \bar{D}_s^2) q_i |v_i - v_s|
\]

\[
E_{is} = 0.1 \cdot \exp(0.025 T)
\]

#### 3.5.5 Melting and Freezing

**Melting rate (T > 0°C):**

\[
\left(\frac{dq_x}{dt}\right)_{\text{melt}} = \frac{2\pi N_x \bar{D}_x K_a (T - T_0)}{L_f}
\]

where:
- \(K_a = 2.4 \times 10^{-2}\) W m⁻¹ K⁻¹: thermal conductivity of air
- \(L_f = 3.34 \times 10^5\) J kg⁻¹: latent heat of fusion

**Homogeneous freezing (T < -40°C):**

\[
\left(\frac{dq_i}{dt}\right)_{\text{frz}} = \frac{q_c}{\Delta t} \quad \text{(instant)}
\]

**Heterogeneous freezing (immersion):**

\[
\left(\frac{dN_i}{dt}\right)_{\text{imm}} = N_c \cdot J_{\text{het}}(T) \cdot V_{\text{drop}}
\]

\[
J_{\text{het}}(T) = 10^{-5} \exp\left[-0.6(T + 12.96)\right] \text{ m}^{-3}\text{ s}^{-1}
\]

### 3.6 Aerosol-Cloud Interactions (NEW 2024)

**CCN activation:**

\[
S_{\text{max}} = C_1 w^{3/4} N_{\text{CCN}}^{-1/2}
\]

**Activated droplet number:**

\[
N_{c,\text{act}} = C_2 N_{\text{CCN}}^{k_1} S_{\text{max}}^{k_2}
\]

where \(k_1 = 0.7\), \(k_2 = 1.5\) for continental aerosol.

**Competition between natural INP and AgI:**

\[
N_{i,\text{total}} = N_{\text{INP}}^{\text{natural}} + N_{\text{INP}}^{\text{AgI}} \cdot (1 - \beta)
\]

\[
\beta = \frac{N_{\text{INP}}^{\text{natural}}}{N_{\text{INP}}^{\text{natural}} + 10^4} \quad \text{(saturation effect)}
\]

---

## 4. Cloud Module: Dynamics

### 4.1 Compressible Navier-Stokes

**Momentum:**

\[
\frac{\partial u_i}{\partial t} = -u_j \frac{\partial u_i}{\partial x_j} - \frac{1}{\rho}\frac{\partial p}{\partial x_i} + g\delta_{i3}\frac{\theta_v - \theta_{v0}}{\theta_{v0}} + F_i
\]

where \(F_i\): subgrid turbulent forcing (LES/DNS).

**Continuity:**

\[
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) = 0
\]

**Thermodynamic energy:**

\[
\frac{\partial \theta}{\partial t} = -\vec{u} \cdot \nabla \theta + \frac{1}{\rho c_p}\left(L_v \frac{dq_v}{dt} + L_s \frac{dq_i}{dt}\right) + Q_{\text{rad}}
\]

where:
- \(L_v = 2.5 \times 10^6\) J kg⁻¹: latent heat of vaporization
- \(L_s = 2.834 \times 10^6\) J kg⁻¹: latent heat of sublimation

### 4.2 Subgrid Turbulence (LES)

**Smagorinsky model:**

\[
\tau_{ij} = -2\rho \nu_t S_{ij}
\]

\[
\nu_t = (C_s \Delta)^2 \sqrt{2S_{ij}S_{ij}}
\]

where \(C_s = 0.18\), \(\Delta = (dx \cdot dy \cdot dz)^{1/3}\).

**Turbulent enhancement of coagulation (NEW):**

\[
\left(\frac{dN}{dt}\right)_{\text{turb}} = -K_{\text{turb}} \epsilon^{1/2} N^2 r^3
\]

where:
- \(K_{\text{turb}} = 1.3\): turbulent coagulation constant
- \(\epsilon\): turbulent dissipation rate (m² s⁻³)

---

## 5. FogSeeding Module: Specialized

### 5.1 Fog Microphysics

**Droplet distribution (modified gamma):**

\[
n(r) = N_0 r^\alpha \exp(-\lambda r)
\]

**Liquid water content:**

\[
L = \frac{4\pi \rho_w}{3} \int_0^\infty r^3 n(r) dr = \frac{4\pi \rho_w N_0}{3\lambda^{3+\alpha}} \Gamma(4 + \alpha)
\]

**Visibility:**

\[
V = \frac{3.912}{Q_{\text{ext}}}
\]

\[
Q_{\text{ext}} = \pi \int_0^\infty r^2 Q(r, \lambda_{\text{light}}) n(r) dr
\]

where \(Q(r, \lambda)\): Mie efficiency factor.

### 5.2 Hygroscopic Seeding

**Droplet growth on hygroscopic nuclei (CaCl₂, NaCl):**

\[
\frac{dr}{dt} = \frac{S - S_{\text{eq}}(r, m_s)}{F_k + F_d}
\]

**Köhler equation:**

\[
S_{\text{eq}} = \frac{a}{r} - \frac{b}{r^3}
\]

\[
a = \frac{2\sigma M_w}{\rho_w R T}, \quad b = \frac{3i \nu m_s M_w}{4\pi \rho_w M_s}
\]

where:
- \(m_s\): mass of solute (kg)
- \(\nu\): van't Hoff factor (2.8 for CaCl₂)
- \(i\): dissociation factor

**Growth rate includes diffusion and kinetics:**

\[
F_k = \frac{L_v^2 M_w}{R T^2 K_a}, \quad F_d = \frac{R T}{\rho_w D_v e_{sw} M_w}
\]

---

## 6. Numerical Methods

### 6.1 Spatial Discretization

**Advection (5th-order WENO):**

\[
\frac{\partial q}{\partial x}\bigg|_i \approx \frac{1}{\Delta x}\sum_{k=-2}^{2} w_k q_{i+k}
\]

**Benefits**: minimizes numerical diffusion, captures steep gradients.

### 6.2 Time Integration

**5-stage Runge-Kutta (SSP-RK5):**

\[
\begin{aligned}
q^{(1)} &= q^n + \alpha_1 \Delta t \mathcal{L}(q^n) \\
q^{(2)} &= q^n + \alpha_2 \Delta t \mathcal{L}(q^{(1)}) \\
&\vdots \\
q^{n+1} &= q^{(5)}
\end{aligned}
\]

where \(\alpha = [0.37, 0.38, 0.18, 0.62, 1.0]\).

### 6.3 CFL Condition

\[
\Delta t \leq \text{CFL} \cdot \min\left(\frac{\Delta x}{|u|}, \frac{\Delta y}{|v|}, \frac{\Delta z}{|w|}\right)
\]

Use \(\text{CFL} = 0.5\) for stability.

---

## 7. Initialization and Boundary Conditions

### 7.1 Atmospheric Sounding

**Input required:**
- Temperature profile \(T(z)\)
- Dew point \(T_d(z)\)
- Wind \(u(z), v(z)\)
- Pressure \(p(z)\)

**Compute stability indices:**

\[
\text{CAPE} = \int_{LFC}^{EL} g\frac{T_p - T_e}{T_e} dz
\]

\[
\text{CIN} = \int_{z_0}^{LFC} g\frac{T_p - T_e}{T_e} dz
\]

### 7.2 Boundary Conditions

**Top (z = z_top): absorbing (sponge layer)**

\[
\frac{\partial \phi}{\partial t}\bigg|_{\text{sponge}} = -\nu_{\text{sponge}}(z) (\phi - \phi_{\text{ref}})
\]

\[
\nu_{\text{sponge}}(z) = \nu_0 \sin^2\left(\frac{\pi(z - z_s)}{2(z_{\text{top}} - z_s)}\right)
\]

**Bottom (z = 0): no-slip or slip**

\[
u = v = w = 0 \quad \text{(no-slip)}
\]

**Lateral: periodic or open (radiative)**

\[
\frac{\partial \phi}{\partial t} + c_{\text{phase}} \frac{\partial \phi}{\partial x} = 0
\]

### 7.3 Seeding Source Initialization

**Aircraft line source:**

\[
S_C(x,y,z,t) = \frac{Q_{\text{rate}}}{\sqrt{2\pi\sigma_y\sigma_z}} \exp\left[-\frac{(y-y_0)^2}{2\sigma_y^2} - \frac{(z-z_0)^2}{2\sigma_z^2}\right] \delta(x - x_0 - v_a t)
\]

where:
- \(Q_{\text{rate}}\): emission rate (particles/s)
- \(\sigma_y, \sigma_z\): plume dispersion (50-100 m initially)
- \(v_a\): aircraft speed

**Ground generator (vertical plume):**

\[
S_C(x,y,z,t) = \frac{Q_{\text{rate}}}{2\pi\sigma_x\sigma_y u} \exp\left[-\frac{(x-x_0)^2}{2\sigma_x^2} - \frac{(y-y_0)^2}{2\sigma_y^2}\right] f_z(z)
\]

\[
f_z(z) = \begin{cases}
\frac{z^2}{\sigma_z^3}\exp(-z/\sigma_z) & z > 0 \\
0 & z \leq 0
\end{cases}
\]

---

## 8. Model Output and Diagnostics

### 8.1 Primary Variables (every time step)

- 3D fields: \(u, v, w, T, p, \rho\)
- Hydrometeors: \(q_c, q_r, q_i, q_s, q_g, N_c, N_r, N_i, N_s, N_g\)
- Reagents: \(C_{\text{AgI}}\)
- Aerosol: \(N_{\text{CCN}}, N_{\text{INP}}\)

### 8.2 Derived Diagnostics

**Radar reflectivity (dBZ):**

\[
Z = 10\log_{10}\left(\sum_x \frac{N_x |K_w|^2}{\rho \lambda^4} \int D^6 n_x(D) dD\right)
\]

where \(|K_w|^2 = 0.93\) for water.

**Precipitation rate (mm/h):**

\[
R = 3600 \times \rho_w^{-1} \int_{D_{\min}}^{D_{\max}} \frac{\pi D^3}{6} v_t(D) n(D) dD
\]

**Ice water path (kg/m²):**

\[
\text{IWP} = \int_{z_{\text{cloud base}}}^{z_{\text{cloud top}}} \rho (q_i + q_s + q_g) dz
\]

**Cloud optical depth:**

\[
\tau = \int Q_{\text{ext}}(z) dz
\]

### 8.3 Seeding Efficiency Metrics

**Enhanced precipitation (target - control):**

\[
\Delta P = \frac{P_{\text{seeded}} - P_{\text{control}}}{P_{\text{control}}} \times 100\%
\]

**Ice enhancement ratio:**

\[
\text{IER} = \frac{N_{i,\text{seeded}}}{N_{i,\text{unseeded}}}
\]

**Optimal IER: 10-100** (too high = overseeding effect).

**Reagent utilization:**

\[
\eta = \frac{\text{mass of precipitation increase}}{\text{mass of reagent used}}
\]

---

## 9. Validation Strategy

### 9.1 Ensemble Approach

**Generate N=30 members varying:**
- Initial conditions (±5% perturbation)
- CCN concentration (50-2000 cm⁻³)
- AgI emission rate (±20%)
- Microphysics parameters (within uncertainty)

**Compute ensemble mean and spread:**

\[
\bar{P} = \frac{1}{N}\sum_{i=1}^N P_i, \quad \sigma_P = \sqrt{\frac{1}{N-1}\sum_{i=1}^N (P_i - \bar{P})^2}
\]

### 9.2 Multi-Sensor Validation

**Required observations:**
1. Weather radar: Z, ZDR, KDP, ρhv (polarimetric)
2. Microwave radiometer: liquid water path, temperature profile
3. In-situ aircraft: N_c, N_i, LWC, IWC, PSD
4. Surface: precipitation gauges, disdrometer
5. Satellite: cloud top temperature, optical depth

**Statistical metrics:**

\[
\text{RMSE} = \sqrt{\frac{1}{N}\sum_{i=1}^N (P_{\text{obs},i} - P_{\text{model},i})^2}
\]

\[
\text{Bias} = \frac{1}{N}\sum_{i=1}^N (P_{\text{model},i} - P_{\text{obs},i})
\]

\[
\text{Correlation} = \frac{\sum(P_{\text{obs}} - \bar{P}_{\text{obs}})(P_{\text{model}} - \bar{P}_{\text{model}})}{\sigma_{\text{obs}} \sigma_{\text{model}}}
\]

**Target: RMSE < 20%, |Bias| < 10%, r > 0.75**

---

## 10. Machine Learning Enhancement (Optional 2024+)

### 10.1 Physics-Informed Neural Network Emulator

**Train NN to emulate microphysics:**

\[
\frac{d\vec{q}}{dt} = f_{\text{NN}}(\vec{q}, T, p, w; \vec{\theta}) + \text{Physics Constraint}
\]

**Loss function:**

\[
\mathcal{L} = \underbrace{\|\vec{q}_{\text{pred}} - \vec{q}_{\text{true}}\|^2}_{\text{Data loss}} + \lambda_1 \underbrace{\left\|\frac{\partial \vec{q}}{\partial t} - \mathcal{L}_{\text{phys}}(\vec{q})\right\|^2}_{\text{Physics loss}} + \lambda_2 \|\vec{\theta}\|^2
\]

**Benefits**: 50-80% speedup, retain physical consistency.

### 10.2 Conditional GAN for Subgrid Variability

**Generator:**

\[
C_{\text{fine}} = G(C_{\text{coarse}}, z | \theta_G)
\]

where \(z \sim \mathcal{N}(0, I)\): latent noise.

**Discriminator loss:**

\[
\mathcal{L}_D = -\mathbb{E}[\log D(C_{\text{real}})] - \mathbb{E}[\log(1 - D(G(C_{\text{coarse}}, z)))]
\]

---

## 11. Computational Considerations

### 11.1 Parallelization

**Domain decomposition (MPI):**
- Split domain into subdomains (e.g., 8×8×4 for 256 cores)
- Halo exchange: 3 ghost cells per boundary
- Communication pattern: non-blocking (MPI_Isend/Irecv)

**Shared memory (OpenMP):**
- Parallelize loops over grid points
- Use `#pragma omp parallel for collapse(3)` for 3D loops

### 11.2 GPU Acceleration (CUDA/HIP)

**Offload kernels:**
1. Advection operator: 90% of compute time
2. Microphysics rates: column-parallel
3. Diffusion: tridiagonal solvers

**Expected speedup: 10-50× vs. CPU**

### 11.3 Adaptive Mesh Refinement (AMR)

**Refine where:**

\[
\epsilon_{\text{local}} = \max\left(\left|\nabla w\right|, \left|\nabla q_c\right|\right) > \epsilon_{\text{threshold}}
\]

**Benefits**: 3-5× reduction in grid points for similar accuracy.

---

## 12. Implementation Checklist

### Phase 1: Core Infrastructure
- [ ] Grid setup and domain decomposition
- [ ] Time integration framework (RK5)
- [ ] Advection scheme (WENO-5)
- [ ] Diffusion solver (implicit or explicit)
- [ ] I/O system (NetCDF format)

### Phase 2: Basic Microphysics
- [ ] Single-moment warm rain (Kessler)
- [ ] Terminal velocity calculations
- [ ] Saturation adjustment
- [ ] Condensation/evaporation

### Phase 3: Double-Moment Upgrade
- [ ] Add N_x prognostic equations
- [ ] Morrison scheme for ice
- [ ] Self-collection and breakup
- [ ] Gamma distribution moments

### Phase 4: Seeding Physics
- [ ] AgI dispersion module
- [ ] Ice nucleation parameterization (2024 update)
- [ ] Secondary ice production (3 mechanisms)
- [ ] Competition effects

### Phase 5: Advanced Features
- [ ] Aerosol-cloud interactions
- [ ] LES turbulence model
- [ ] Radiation coupling
- [ ] Ensemble framework

### Phase 6: Validation
- [ ] Radar forward operator
- [ ] Multi-sensor comparison
- [ ] Statistical evaluation (30+ cases)
- [ ] Sensitivity analysis

### Phase 7: Optimization (Optional)
- [ ] MPI parallelization
- [ ] GPU porting (CUDA/HIP)
- [ ] ML emulator training
- [ ] AMR implementation

---

## 13. Key Parameters Summary

| Parameter | Symbol | Value | Unit | Source |
|-----------|--------|-------|------|--------|
| von Kármán constant | k | 0.4 | - | Standard |
| Air density (STP) | ρ₀ | 1.225 | kg/m³ | Standard |
| Water density | ρ_w | 1000 | kg/m³ | Standard |
| Ice density | ρ_i | 917 | kg/m³ | Standard |
| Latent heat vaporization | L_v | 2.5×10⁶ | J/kg | Standard |
| Latent heat fusion | L_f | 3.34×10⁵ | J/kg | Standard |
| Latent heat sublimation | L_s | 2.834×10⁶ | J/kg | Standard |
| Thermal conductivity air | K_a | 2.4×10⁻² | W/m/K | Standard |
| AgI optimal size | d_AgI | 50-200 | nm | Jiang 2025 |
| AgI INF at -8°C | INF | 0.07-1.6 | % | Jiang 2025 |
| Hallett-Mossop coeff | C_HM | 3.5×10⁸ | kg⁻¹ | Korolev 2022 |
| Collisional breakup | C_coll | 9×10⁻⁵ | m³/kg | Phillips 2024 |
| Sublimational breakup | C_subl | 5×10⁶ | kg⁻¹ | Morrison 2024 |
| Smagorinsky constant | C_s | 0.18 | - | Standard |
| CFL number | CFL | 0.5 | - | Stability |

---

## 14. Critical Implementation Notes

### ⚠️ Common Pitfalls
1. **Underestimating secondary ice**: Can increase N_i by 10-100×
2. **Ignoring AgI size distribution**: Only 50-200 nm particles effective
3. **Overseeding**: IER > 1000 reduces precipitation (competition)
4. **Mass-number inconsistency**: Always enforce \(q_x = \rho \bar{m}_x N_x\)
5. **Negative values**: Apply flux limiters to prevent \(q_x < 0\) or \(N_x < 0\)

### ✅ Best Practices
1. **Conservation checks**: Verify total water, energy, momentum at each step
2. **Diagnostic output**: Save 3D fields every 5-10 min for analysis
3. **Restart capability**: Write checkpoint every 30 min of simulation
4. **Unit testing**: Test each microphysical process in isolation
5. **Code documentation**: Document all constants and parameter choices

---

## 15. References (Key 2023-2025)

- Morrison et al. (2024): Updated double-moment scheme with secondary ice
- Jiang et al. (2025): AgI ice nucleation efficiency (ACP)
- French et al. (2023): SNOWIE field campaign validation (JAMC)
- Korolev et al. (2024): Secondary ice production mechanisms
- Xue et al. (2024): Machine learning microphysics emulators (EDS)

---

## End of Implementation Guide

**Version**: 1.0  
**Last Updated**: 2026-01-30  
**Status**: Ready for implementation

**Next Steps**:
1. Set up development environment (Python/Fortran + MPI)
2. Implement Phase 1 (core infrastructure)
3. Validate with idealized test cases (rising bubble, squall line)
4. Progress to real case studies

---
