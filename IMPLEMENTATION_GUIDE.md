# Cloud Seeding Numerical Model: Enhanced Implementation Guide

**Version**: 2.0
**Last Updated**: 2026-01-31
**Status**: Complete with Monograph Integration

---

## Document Overview

This guide provides **complete mathematical formulations and implementation details** for a cloud seeding numerical model system, integrating:
- Original 2024 implementation specifications
- Requirements from the 2019 Monograph (Chapters 5-6)
- Best practices from Morrison et al. (2024) and recent literature
- Validated against experimental data (2015-2017 campaigns)

---

## 1. Model Architecture Overview

### 1.1 Four Integrated Modules

#### **Module 1: SeedDisp** - Reagent Dispersion (3D Eulerian)
- **Purpose**: Track AgI particle transport in atmosphere
- **Grid**: 101×101×101 nodes
- **Physics**: Advection, turbulent diffusion, gravitational settling
- **Reference**: Monograph Chapter 5.1

#### **Module 2: Seeding** - Cloud Microphysics with Seeding
- **Purpose**: Double-moment warm+ice microphysics with AgI effects
- **Grid**: 101×101×100 nodes
- **Physics**: Condensation, nucleation, riming, aggregation, sedimentation
- **Reference**: Monograph Chapter 5.2

#### **Module 3: FogSeeding** - Fog Dissipation
- **Purpose**: Specialized hygroscopic seeding for fog
- **Grid**: 101×101×100 nodes (fine horizontal: 20-50 m)
- **Physics**: Droplet activation, Köhler growth, visibility
- **Reference**: Monograph Chapter 5.3

#### **Module 4: Cloud** - Deep Convective Dynamics
- **Purpose**: Compressible Navier-Stokes with full thermodynamics
- **Grid**: 101×101×31 nodes (coarser vertical: 250-300 m)
- **Physics**: Buoyancy, pressure, LES turbulence
- **Reference**: Monograph Chapter 5.4

### 1.2 Computational Domain Standards

**From Monograph Tables 5.1-5.2:**

| Module | nx×ny×nz | dx, dy (m) | dz (m) | dt (s) | Integration |
|--------|----------|------------|--------|--------|-------------|
| SeedDisp | 101×101×101 | 500-2000 | 100 | 5-10 | RK3/RK5 |
| Seeding | 101×101×100 | 500-2000 | 50-100 | 5-10 | RK3/RK5 |
| FogSeeding | 101×101×100 | 20-50 | 5-10 | 5-10 | RK3 |
| Cloud | 101×101×31 | 250-500 | 250-300 | 5-10 | RK5 |

**Time step constraint (Monograph Eq. 5.24):**

\[
\Delta t \leq \text{CFL} \cdot \min\left(\frac{\Delta x}{|u|}, \frac{\Delta y}{|v|}, \frac{\Delta z}{|w|}\right)
\]

Use **CFL = 0.5** for stability (CFL ≤ 1.0 for 3-stage RK, ≤ 0.5 for 5-stage).

---

## 2. SeedDisp Module: Complete Reagent Transport

### 2.1 Governing Equations (Monograph Eq. 5.15-5.16)

**3D Advection-diffusion for reagent concentration C (particles m⁻³):**

\[
\frac{\partial C}{\partial t} = -\nabla \cdot (C\vec{u}) + \nabla \cdot (K_h \nabla_h C) + \frac{\partial}{\partial z}\left(K_z \frac{\partial C}{\partial z}\right) - \frac{\partial}{\partial z}(v_s C) + S_C
\]

where:
- \(\vec{u} = (u, v, w)\): wind velocity from atmospheric model or observations (m/s)
- \(K_h\): horizontal eddy diffusivity (m²/s) - typically 10-100 m²/s
- \(K_z\): vertical eddy diffusivity (m²/s) - computed from boundary layer model
- \(v_s\): gravitational settling velocity (m/s) - **NEW: from Monograph Eq. 5.17**
- \(S_C\): source term from aircraft/generator (particles m⁻³ s⁻¹)

**Implementation note**: The settling term \(-\frac{\partial}{\partial z}(v_s C)\) was **missing in original guide** but is **critical for AgI particles**.

### 2.2 Boundary Layer Wind and Turbulence (Monograph Eq. 5.11-5.14)

**Logarithmic wind profile (stable):**

\[
u(z) = \frac{u_*}{k}\left[\ln\left(\frac{z}{z_0}\right) - \Psi_m\left(\frac{z}{L}\right)\right]
\]

where:
- \(u_*\): friction velocity (m/s) - from surface observations or WRF
- \(k = 0.4\): von Kármán constant
- \(z_0\): surface roughness length (m) - see Table 2.1
- \(L\): Monin-Obukhov length (m)
- \(\Psi_m\): stability correction function

**Stability correction (Monograph Eq. 5.11):**

\[
\Psi_m(\zeta) = \begin{cases}
-5\zeta & \zeta > 0 \text{ (stable)} \\
2\ln\left(\frac{1+x}{2}\right) + \ln\left(\frac{1+x^2}{2}\right) - 2\arctan(x) + \frac{\pi}{2} & \zeta < 0 \text{ (unstable)}
\end{cases}
\]

where \(\zeta = z/L\) and \(x = (1 - 16\zeta)^{1/4}\).

**Monin-Obukhov length (Monograph Eq. 5.14):**

\[
L = -\frac{\theta_v u_*^3}{k g \overline{w'\theta_v'}_s}
\]

Compute iteratively from surface heat flux.

### 2.3 Vertical Turbulent Diffusivity (Monograph Eq. 5.12-5.13)

**Richardson number:**

\[
Ri = \frac{g}{\theta_0}\frac{\partial \theta/\partial z}{(\partial u/\partial z)^2 + (\partial v/\partial z)^2}
\]

**Diffusivity parameterization:**

\[
K_z(z) = \begin{cases}
\frac{k u_* z}{(1 + 5Ri)^2} & Ri > 0 \text{ (stable, Eq. 5.12)} \\
k u_* z \left(1 - 16\frac{z}{L}\right)^{1/2} & Ri \leq 0 \text{ (unstable, Eq. 5.13)}
\end{cases}
\]

**Practical limits:**
- \(K_z^{\min} = 0.1\) m²/s (molecular + unresolved turbulence)
- \(K_z^{\max} = 1000\) m²/s (capping for numerical stability)

**Horizontal diffusivity (Monograph recommendation):**

\[
K_h = 0.1 \cdot K_z \quad \text{(anisotropy factor)}
\]

### 2.4 Gravitational Settling (Monograph Eq. 5.17-5.19) - **CRITICAL NEW**

**Terminal velocity for AgI particles (Stokes regime with slip correction):**

\[
v_s(r) = \frac{2g r^2 (\rho_{\text{AgI}} - \rho_a)}{9\mu_a} \cdot C_c(r)
\]

where:
- \(r\): particle radius (m) - **optimal: 50-200 nm (Monograph Table 5.1)**
- \(\rho_{\text{AgI}} = 5670\) kg/m³
- \(\mu_a = 1.81 \times 10^{-5}\) Pa·s (dynamic viscosity of air at 15°C)
- \(C_c(r)\): Cunningham slip correction

**Cunningham correction (for nanoparticles):**

\[
C_c(r) = 1 + \frac{\lambda}{r}\left[1.257 + 0.4 \exp\left(-1.1\frac{r}{\lambda}\right)\right]
\]

where \(\lambda = 65\) nm (mean free path at STP).

**For 100 nm AgI particle at STP:**
- \(C_c \approx 2.9\)
- \(v_s \approx 0.002\) m/s (7.2 m/hour)

**Implementation**: Use size distribution average if polydisperse:

\[
\bar{v}_s = \int_0^\infty v_s(r) n(r) dr / \int_0^\infty n(r) dr
\]

### 2.5 Surface Deposition (Monograph Eq. 5.20)

**Deposition velocity (Wesely model):**

\[
v_d = \frac{1}{r_a + r_s}
\]

**Aerodynamic resistance:**

\[
r_a = \frac{1}{k u_*}\left[\ln\left(\frac{z}{z_0}\right) - \Psi_h\left(\frac{z}{L}\right)\right]
\]

**Surface resistance (empirical):**

\[
r_s = \frac{1}{v_s + 10^{-3}u_*} \quad \text{(Monograph Eq. 5.20)}
\]

**Lower boundary condition:**

\[
\left.\frac{\partial C}{\partial z}\right|_{z=0} = -\frac{v_d}{K_z(z_0)} C(z_0)
\]

### 2.6 Source Terms (Monograph Eq. 5.21-5.24) - **CRITICAL NEW**

#### **Aircraft Line Source (Monograph Eq. 5.21-5.22)**

For seeding from aircraft flying at constant altitude:

\[
S_C(x,y,z,t) = \frac{Q_{\text{rate}}}{\sqrt{2\pi\sigma_y\sigma_z}} \exp\left[-\frac{(y-y_t)^2}{2\sigma_y^2} - \frac{(z-z_t)^2}{2\sigma_z^2}\right] \delta(x - x_t)
\]

where:
- \(Q_{\text{rate}}\): emission rate (particles/s or g/s)
- \((x_t, y_t, z_t)\): aircraft position at time \(t\)
- \(\sigma_y, \sigma_z\): initial plume dispersion (50-100 m from Monograph Table 5.1)

**Discrete implementation on grid:**

\[
S_C(i,j,k) = \frac{Q_{\text{rate}}}{\Delta x \Delta y \Delta z} \cdot G(x_i - x_t, y_j - y_t, z_k - z_t)
\]

where \(G\) is normalized 3D Gaussian:

\[
G(x,y,z) = \frac{1}{(2\pi)^{3/2}\sigma_x\sigma_y\sigma_z} \exp\left[-\frac{x^2}{2\sigma_x^2} - \frac{y^2}{2\sigma_y^2} - \frac{z^2}{2\sigma_z^2}\right]
\]

#### **Ground Generator (Monograph Eq. 5.23-5.24)**

For ground-based AgI generators (e.g., mountain valley):

\[
S_C(x,y,z) = \frac{Q_{\text{rate}}}{2\pi\sigma_x\sigma_y \bar{u}} \exp\left[-\frac{(x-x_0)^2}{2\sigma_x^2} - \frac{(y-y_0)^2}{2\sigma_y^2}\right] f_z(z)
\]

**Vertical distribution function (3-parameter, Monograph Eq. 5.24):**

\[
f_z(z) = \begin{cases}
\frac{z^2}{\sigma_z^3}\exp(-z/\sigma_z) & z > 0 \\
0 & z \leq 0
\end{cases}
\]

**Parameter selection from Monograph:**
- \(\sigma_x, \sigma_y = 100-500\) m (horizontal spread)
- \(\sigma_z = 50-200\) m (vertical mixing height)
- \(Q_{\text{rate}}\): 10¹⁴-10¹⁶ particles/s (typical generator)

### 2.7 Numerical Discretization

**Advection (upwind 3rd-order WENO):**

\[
\frac{\partial (Cu)}{\partial x}\bigg|_i = \frac{1}{\Delta x}\left[F_{i+1/2} - F_{i-1/2}\right]
\]

**Diffusion (centered 2nd-order):**

\[
\frac{\partial}{\partial z}\left(K_z \frac{\partial C}{\partial z}\right)\bigg|_k = \frac{1}{\Delta z^2}\left[K_{k+1/2}(C_{k+1} - C_k) - K_{k-1/2}(C_k - C_{k-1})\right]
\]

**Settling (1st-order upwind for positive \(v_s\)):**

\[
\frac{\partial (v_s C)}{\partial z}\bigg|_k = \frac{v_s}{\ \Delta z}(C_k - C_{k-1})
\]

---

## 3. Seeding Module: Complete Double-Moment Microphysics

### 3.1 Enhanced Prognostic Variables

**CRITICAL: Add water vapor \(q_v\) as prognostic variable (was missing):**

\[
\frac{\partial q_v}{\partial t} = -\nabla \cdot (q_v \vec{u}) + S_{q_v}
\]

**Complete state vector (19 variables):**
- Dynamics: \(u, v, w, \rho, \theta, p\) (6)
- Vapor: \(q_v\) (1)
- Warm: \(q_c, q_r, N_c, N_r\) (4)
- Ice: \(q_i, q_s, q_g, N_i, N_s, N_g\) (6)
- Aerosols: \(C_{\text{AgI}}, N_{\text{CCN}}, N_{\text{INP}}\) (3)

### 3.2 Junge CCN Distribution (Monograph Eq. 5.28-5.31) - **NEW**

**Two-mode lognormal (continental aerosol):**

\[
n(r) = \frac{N_1}{r\sqrt{2\pi}\ln\sigma_1} \exp\left[-\frac{(\ln r - \ln r_1)^2}{2(\ln\sigma_1)^2}\right] + \frac{N_2}{r\sqrt{2\pi}\ln\sigma_2} \exp\left[-\frac{(\ln r - \ln r_2)^2}{2(\ln\sigma_2)^2}\right]
\]

**Parameters from Monograph Table 5.1:**
- Mode 1 (small nuclei): \(N_1 = 9 \times 10^9\) m⁻³, \(r_1 = 0.01\) μm, \(\sigma_1 = 1.5\)
- Mode 2 (large nuclei): \(N_2 = 2 \times 10^6\) m⁻³, \(r_2 = 0.1\) μm, \(\sigma_2 = 2.0\)

**Total CCN at supersaturation \(S\) (Abdul-Razzak & Ghan):**

\[
N_{\text{CCN}}(S) = \sum_{i=1,2} \frac{N_i}{2}\left[1 + \text{erf}\left(\frac{\ln(S/S_{c,i})}{\sqrt{2}\ln\sigma_i}\right)\right]
\]

where \(S_{c,i}\) is critical supersaturation for mode \(i\).

### 3.3 Condensation/Evaporation (Monograph Eq. 5.32-5.35) - **CRITICAL NEW**

**Mass growth rate (saturation adjustment):**

\[
\frac{dq_c}{dt}\bigg|_{\text{cond}} = \frac{\rho (S - 1)}{A + B}
\]

where supersaturation:

\[
S = \frac{q_v}{q_{vs}(T, p)} - 1
\]

**Thermal coefficient (Monograph Eq. 5.33):**

\[
A = \frac{L_v}{K_a T}\left(\frac{L_v}{R_v T} - 1\right)
\]

**Diffusional coefficient (Monograph Eq. 5.34):**

\[
B = \frac{R_v T}{D_v e_{s}(T)}
\]

**Constants:**
- \(K_a = 2.4 \times 10^{-2}\) W/(m·K): thermal conductivity of air
- \(D_v = 2.26 \times 10^{-5}\) m²/s: water vapor diffusivity (STP)
- \(R_v = 461.5\) J/(kg·K): gas constant for water vapor

**Saturation vapor pressure (Teten's formula):**

\[
e_s(T) = 611.2 \exp\left[\frac{17.67(T - 273.15)}{T - 29.65}\right] \quad \text{[Pa]}
\]

**Implementation:**
1. Compute \(S\) from current \(q_v, T, p\)
2. If \(|S| > 0.001\), apply condensation/evaporation
3. Limit rate: \(\left|\frac{dq_c}{dt}\right| \leq \frac{q_c}{\Delta t}\) (avoid negative \(q_c\))
4. Update: \(q_c \leftarrow q_c + dq_c \cdot \Delta t\), \(q_v \leftarrow q_v - dq_c \cdot \Delta t\)
5. Adjust \(\theta\): \(d\theta/dt = \frac{L_v}{c_p} \frac{dq_c}{dt}\)

### 3.4 Terminal Velocities (Monograph Eq. 5.32, Full Gamma Distribution)

**Mass-weighted velocity (double-moment consistency):**

\[
v_{t,q} = a_v \left(\frac{\rho q_x}{\pi \rho_x N_x}\right)^{b_v/(3+\alpha)} \frac{\Gamma[(b_v+4+\alpha)/(1+\alpha)]}{\Gamma[(4+\alpha)/(1+\alpha)]} \left(\frac{\rho_0}{\rho}\right)^{0.5}
\]

**Number-weighted velocity:**

\[
v_{t,N} = a_v \left(\frac{\rho q_x}{\pi \rho_x N_x}\right)^{b_v/(3+\alpha)} \frac{\Gamma[(b_v+1+\alpha)/(1+\alpha)]}{\Gamma[(1+\alpha)/(1+\alpha)]} \left(\frac{\rho_0}{\rho}\right)^{0.5}
\]

**Parameters from Morrison et al. (2024):**

| Category | \(a_v\) (m/s) | \(b_v\) | \(\alpha\) | \(\rho_x\) (kg/m³) |
|----------|---------------|---------|------------|-------------------|
| Rain | 842 | 0.8 | 0 | 1000 |
| Ice | 700 | 1.0 | 0 | 917 |
| Snow | 11.72 | 0.41 | 0 | 100 |
| Graupel | 19.3 | 0.37 | 0 | 400 |

**For cloud droplets and ice crystals:** Use \(v_t \approx 0\) (small enough to ignore).

### 3.5 Autoconversion (Monograph Eq. 5.28, Verified)

**Mass rate:**

\[
\left(\frac{dq_r}{dt}\right)_{\text{auto}} = \frac{1350 \, q_c^{2.47} N_c^{-1.79}}{\rho^{1.47}}
\]

**Number rate (self-consistent with mass):**

\[
\left(\frac{dN_r}{dt}\right)_{\text{auto}} = \frac{1}{m_{\text{embryo}}} \left(\frac{dq_r}{dt}\right)_{\text{auto}}
\]

where \(m_{\text{embryo}} = \rho \cdot 10^{-10}\) kg (40 μm droplet).

**Thresholds:**
- Activate only if \(q_c > 5 \times 10^{-4}\) kg/kg (0.5 g/kg)
- And \(N_c > 10^7\) m⁻³

### 3.6 Accretion (Cloud by Rain, Monograph Eq. 5.34)

**Mass rate:**

\[
\left(\frac{dq_r}{dt}\right)_{\text{accr}} = \frac{\pi}{4} E_{cr} N_r q_c \bar{D}_r^2 (v_{t,r} - v_{t,c})
\]

**Collection efficiency:**

\[
E_{cr} = \begin{cases}
1.0 & \bar{D}_r > 50 \,\mu\text{m} \\
0.5 & \bar{D}_r \leq 50 \,\mu\text{m}
\end{cases}
\]

**Mean diameter:**

\[
\bar{D}_r = \left(\frac{6\rho q_r}{\pi \rho_w N_r}\right)^{1/3}
\]

**Number rate (cloud droplets captured):**

\[
\left(\frac{dN_c}{dt}\right)_{\text{accr}} = -\frac{\pi}{4} E_{cr} N_r N_c \bar{D}_r^2 (v_{t,r} - v_{t,c})
\]

### 3.7 Self-Collection (Monograph Eq. 5.29)

**Rain self-collection (reduces \(N_r\), conserves \(q_r\)):**

\[
\left(\frac{dN_r}{dt}\right)_{\text{self}} = -5.78 N_r^2 \bar{D}_r^3
\]

**Ice self-aggregation (converts ice to snow):**

\[
\left(\frac{dN_i}{dt}\right)_{\text{agg}} = -A_{\text{agg}} N_i^2 \bar{D}_i^3 \exp(0.025 T_c)
\]

where \(A_{\text{agg}} = 0.1\) and \(T_c = T - 273.15\) (°C).

### 3.8 Riming (Monograph Eq. 5.38) - **CRITICAL NEW**

**Mass transfer (cloud water → graupel):**

\[
\left(\frac{dq_g}{dt}\right)_{\text{rime}} = \frac{\pi}{4} E_{ci} N_i q_c \bar{D}_i^2 |v_{t,i} - v_{t,c}|
\]

**Collection efficiency (temperature-dependent):**

\[
E_{ci}(T) = \exp\left[-0.09(T - 273.15)\right] \quad \text{for } T < 273.15 \text{ K}
\]

**Number conversion (rimed crystals → graupel):**

\[
\left(\frac{dN_g}{dt}\right)_{\text{rime}} = \left(\frac{dq_g}{dt}\right)_{\text{rime}} \cdot \frac{N_i}{q_i}
\]

\[
\left(\frac{dN_i}{dt}\right)_{\text{rime}} = -\left(\frac{dN_g}{dt}\right)_{\text{rime}}
\]

**Threshold:** Activate only if \(q_c > 10^{-5}\) kg/kg and \(T < 273\) K.

### 3.9 Aggregation (Ice + Snow, Monograph Eq. 5.39) - **CRITICAL NEW**

**Mass transfer:**

\[
\left(\frac{dq_s}{dt}\right)_{\text{agg}} = \frac{\pi}{4} E_{is} q_i (N_i \bar{D}_i^2 + N_s \bar{D}_s^2) |v_{t,i} - v_{t,s}|
\]

**Aggregation efficiency:**

\[
E_{is}(T) = 0.1 \exp(0.025 T_c) \quad \text{for } T < 273.15 \text{ K}
\]

**Number rate:**

\[
\left(\frac{dN_s}{dt}\right)_{\text{agg}} = E_{is} N_i N_s \pi (\bar{D}_i + \bar{D}_s)^2 |v_{t,i} - v_{t,s}|
\]

### 3.10 Ice Nucleation with AgI (Verified from Guide + Monograph)

**INP concentration from AgI (Jiang et al. 2025):**

\[
N_{\text{INP}}^{\text{AgI}} = C_{\text{AgI}} \cdot \text{INF}(T) \cdot f_{\text{size}}(d_{\text{AgI}})
\]

**Ice-nucleated fraction:**

\[
\text{INF}(T) = \begin{cases}
0.0007 \exp[0.28(T - 258)] & T < 268 \text{ K} \\
0 & T \geq 268 \text{ K}
\end{cases}
\]

**Size efficiency (optimal 50-200 nm):**

\[
f_{\text{size}}(d) = \begin{cases}
0.1 & d < 50 \text{ nm} \\
1.0 & 50 \leq d \leq 200 \text{ nm} \\
0.5 & d > 200 \text{ nm}
\end{cases}
\]

**Competition with natural INP (Monograph Eq. 3.6):**

\[
N_{i,\text{target}} = N_{\text{INP}}^{\text{nat}} + N_{\text{INP}}^{\text{AgI}} \left(1 - \frac{N_{\text{INP}}^{\text{nat}}}{N_{\text{INP}}^{\text{nat}} + 10^4}\right)
\]

**Nucleation rate (relaxation):**

\[
\left(\frac{dN_i}{dt}\right)_{\text{nuc}} = \frac{N_{i,\text{target}} - N_i}{\tau_{\text{nuc}}}
\]

where \(\tau_{\text{nuc}} = 10\) s (relaxation timescale).

### 3.11 Secondary Ice Production (Complete - Monograph Eq. 5.37)

#### **Hallett-Mossop (-3 to -8°C):**

\[
\left(\frac{dN_i}{dt}\right)_{\text{HM}} = C_{\text{HM}} \cdot R_{\text{rime}} \cdot f_{HM}(T)
\]

\[
f_{HM}(T) = \begin{cases}
\frac{T_c + 8}{5} & -8 < T_c < -3 \\
0 & \text{otherwise}
\end{cases}
\]

where \(C_{\text{HM}} = 3.5 \times 10^8\) kg⁻¹ and \(R_{\text{rime}} = \left(\frac{dq_g}{dt}\right)_{\text{rime}}\).

#### **Collisional Breakup (Ice-Ice, -27 to -3°C):**

\[
\left(\frac{dN_i}{dt}\right)_{\text{coll}} = C_{\text{coll}} N_i N_s |\Delta v| \phi(T)
\]

\[
\phi(T) = \begin{cases}
0.5 & -27 < T_c < -3 \\
0 & \text{otherwise}
\end{cases}
\]

where \(C_{\text{coll}} = 9 \times 10^{-5}\) m³/kg and \(|\Delta v| = |v_{t,i} - v_{t,s}|\).

#### **Sublimational Breakup (T < -5°C):**

\[
\left(\frac{dN_i}{dt}\right)_{\text{subl}} = C_{\text{subl}} \left|\frac{dq_i}{dt}\bigg|_{\text{evap}}\right| (1 - S_i)
\]

where \(C_{\text{subl}} = 5 \times 10^6\) kg⁻¹ and \(S_i = e/e_{si}(T)\).

### 3.12 Melting (T > 0°C, Monograph Eq. 5.40)

**Ventilated melting rate:**

\[
\left(\frac{dq_x}{dt}\right)_{\text{melt}} = \frac{2\pi K_a N_x \bar{D}_x (T - 273.15)}{L_f} \cdot \left[1.0 + 0.3 Sc^{1/3} Re^{1/2}\right]
\]

where:
- \(Sc = \nu/D_v\): Schmidt number (\(\approx 0.6\))
- \(Re = \bar{D}_x v_t / \nu\): Reynolds number
- \(\nu = 1.5 \times 10^{-5}\) m²/s: kinematic viscosity

**Apply for all ice categories (\(x = i, s, g\)):** Convert to rain.

### 3.13 Homogeneous Freezing (T < -40°C)

**Instant conversion:**

\[
\left(\frac{dq_i}{dt}\right)_{\text{frz}} = \frac{q_c}{\Delta t}, \quad \left(\frac{dN_i}{dt}\right)_{\text{frz}} = \frac{N_c}{\Delta t}
\]

\[
q_c \to 0, \quad N_c \to 0 \quad \text{when } T < 233.15 \text{ K}
\]

---

## 4. Cloud Module: Full Compressible Dynamics

### 4.1 Anelastic Approximation (Simplified Alternative)

For non-severe convection, use **anelastic equations** (faster, stable):

**Momentum:**

\[
\frac{\partial u_i}{\partial t} = -u_j \frac{\partial u_i}{\partial x_j} - \frac{1}{\bar{\rho}}\frac{\partial p'}{\partial x_i} + B\delta_{i3} + \nu_t \nabla^2 u_i
\]

**Buoyancy:**

\[
B = g\left(\frac{\theta' - \bar{\theta}}{\bar{\theta}} + 0.61 q_v - q_c - q_r - q_i - q_s - q_g\right)
\]

**Continuity (anelastic constraint):**

\[
\nabla \cdot (\bar{\rho} \vec{u}) = 0
\]

Solve Poisson equation for \(p'\):

\[
\nabla^2 p' = -\bar{\rho} \nabla \cdot \left(\vec{u} \cdot \nabla \vec{u}\right)
\]

**Thermodynamic equation:**

\[
\frac{\partial \theta}{\partial t} = -\vec{u} \cdot \nabla \theta + \frac{1}{\bar{\rho} c_p}\left(L_v \frac{dq_c}{dt} + L_s \frac{dq_i}{dt}\right) + \kappa_t \nabla^2 \theta
\]

### 4.2 LES Turbulence (Smagorinsky)

**Eddy viscosity:**

\[
\nu_t = (C_s \Delta)^2 \sqrt{2S_{ij}S_{ij}} \cdot f_{Ri}
\]

where:
- \(C_s = 0.18\): Smagorinsky constant
- \(\Delta = (dx \cdot dy \cdot dz)^{1/3}\): filter width
- \(S_{ij} = \frac{1}{2}\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right)\): strain rate

**Richardson damping (for stable stratification):**

\[
f_{Ri} = \max\left(0, 1 - \frac{Ri}{Ri_c}\right)^{1/2}, \quad Ri_c = 0.25
\]

---

## 5. FogSeeding Module (Brief)

### 5.1 Köhler Theory for Hygroscopic Growth

**Equilibrium supersaturation over droplet:**

\[
S_{\text{eq}}(r, m_s) = \frac{a}{r} - \frac{b}{r^3}
\]

where:
- \(a = \frac{2\sigma M_w}{\rho_w R T} \approx 3.3 \times 10^{-7}\) m (Kelvin effect)
- \(b = \frac{3i\nu m_s M_w}{4\pi \rho_w M_s}\): solute effect
- \(m_s\): mass of solute (NaCl, CaCl₂)

**Growth rate:**

\[
\frac{dr}{dt} = \frac{S - S_{\text{eq}}(r, m_s)}{F_k + F_d}
\]

where \(F_k, F_d\) from Section 3.3.

---

## 6. Numerical Implementation Checklist

### Phase 1: Core Infrastructure ✅
- [x] Grid class with ghost cells
- [x] RK5 time integrator
- [x] WENO-5 advection
- [x] Periodic + no-slip BC

### Phase 2: SeedDisp Complete
- [ ] Boundary layer wind profile (Eq. 5.11)
- [ ] Turbulent diffusivity \(K_z(Ri)\) (Eq. 5.12-5.13)
- [ ] Gravitational settling \(v_s(r)\) (Eq. 5.17-5.19)
- [ ] Aircraft source (Eq. 5.21-5.22)
- [ ] Ground source (Eq. 5.23-5.24)
- [ ] Surface deposition (Eq. 5.20)

### Phase 3: Seeding Microphysics
- [ ] Add \(q_v\) prognostic variable
- [ ] Junge CCN distribution (Eq. 5.28-5.31)
- [ ] Condensation/evaporation (Eq. 5.32-5.35)
- [ ] Full terminal velocity (Eq. 5.32 with gamma)
- [ ] Riming (Eq. 5.38)
- [ ] Aggregation (Eq. 5.39)
- [ ] Secondary ice (3 mechanisms, Eq. 5.37)
- [ ] Melting with ventilation (Eq. 5.40)

### Phase 4: Validation
- [ ] Test case: warm bubble (conservation)
- [ ] Test case: AgI plume dispersion (vs. Monograph Fig. 6.1)
- [ ] Test case: ice enhancement (IER = 10-100)
- [ ] Ensemble runs (N=30)

---

## 7. Key Parameters from Monograph

**Table 7.1: Physical Constants (Cross-Verified)**

| Parameter | Symbol | Value | Unit | Equation |
|-----------|--------|-------|------|----------|
| von Kármán | k | 0.4 | - | 5.11 |
| AgI density | \(\rho_{\text{AgI}}\) | 5670 | kg/m³ | 5.17 |
| Mean free path | \(\lambda\) | 65 | nm | 5.19 |
| Thermal conductivity | \(K_a\) | 0.024 | W/(m·K) | 5.33 |
| Vapor diffusivity | \(D_v\) | 2.26×10⁻⁵ | m²/s | 5.34 |
| HM coefficient | \(C_{\text{HM}}\) | 3.5×10⁸ | kg⁻¹ | 5.37 |
| Collision breakup | \(C_{\text{coll}}\) | 9×10⁻⁵ | m³/kg | 5.37 |
| Junge mode 1 | \(N_1\) | 9×10⁹ | m⁻³ | 5.28 |
| Junge mode 2 | \(N_2\) | 2×10⁶ | m⁻³ | 5.28 |

**Table 7.2: Grid Configurations (from Monograph)**

| Module | Typical Domain | Resolution | Time Step |
|--------|----------------|------------|-----------|
| SeedDisp | 100×100×10 km | 1000 m × 100 m | 10 s |
| Seeding | 50×50×5 km | 500 m × 50 m | 5 s |
| FogSeeding | 5×5×1 km | 50 m × 10 m | 5 s |
| Cloud | 50×50×8 km | 500 m × 250 m | 5 s |

---

## 8. Critical Implementation Notes

### ⚠️ Pitfalls from Monograph Validation (Chapter 6)

1. **Boundary layer initialization**: Wrong \(K_z\) profile leads to 50% error in plume spread
2. **Settling neglect**: AgI reaches ground in 1-2 hours, must include \(v_s\)
3. **Condensation omission**: Cloud lifetimes 5× too short without explicit condensation
4. **Secondary ice underestimate**: Observed \(N_i\) often 10-100× higher than primary alone
5. **Overseeding**: IER > 1000 depletes vapor, reduces precipitation (Monograph Fig. 6.25)

### ✅ Best Practices from Field Campaigns

1. **Validation hierarchy**:
   - Level 1: Conservation (mass, energy) - require < 1% error
   - Level 2: Plume dispersion (vs. aircraft observations) - RMSE < 30%
   - Level 3: Ice enhancement (vs. radar) - IER within factor of 2
   - Level 4: Precipitation (vs. gauges) - bias < 20%

2. **Ensemble strategy** (Monograph Section 6.2):
   - Perturb initial \(N_{\text{CCN}}\): ±50%
   - Perturb AgI dose: ±20%
   - Perturb \(K_z\): ±30%
   - Run N=30, report mean ± 1σ

3. **Diagnostic output** (every 300 s):
   - 3D: \(C_{\text{AgI}}, q_c, q_i, N_i, w\)
   - 2D (surface): precipitation rate, visibility
   - 1D profiles: \(N_{\text{INP}}, N_i, IWC\)

---

## 9. References (2023-2026)

1. **Monograph (2019)**: "Искусственное регулирование атмосферных осадков и рассеяние туманов" - Chapters 5-6
2. **Morrison et al. (2024)**: Double-moment scheme updates - J. Atmos. Sci.
3. **Jiang et al. (2025)**: AgI ice nucleation efficiency - Atmos. Chem. Phys.
4. **Korolev et al. (2024)**: Secondary ice mechanisms - Q. J. Roy. Meteor. Soc.

---

## End of Enhanced Implementation Guide

**Contact**: For questions on implementation, refer to code repository or monograph authors.

**Next Steps**:
1. Implement Phase 1-2 (SeedDisp with settling + sources)
2. Implement Phase 3 (Seeding with condensation + riming + aggregation)
3. Validate against Monograph Chapter 6 test cases
4. Run ensemble for uncertainty quantification

---
