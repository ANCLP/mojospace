# Particle Physics, Smoothing Kernels, and Numerical Integrators

## 1. Overview

Particle physics simulation engines compute physical trajectories by:
1. Identifying interacting particle pairs (via Spatial Hashing or Trees).
2. Evaluating inter-particle forces (gravity, electrostatics, Lennard-Jones, or SPH pressure/viscosity).
3. Advancing positions and velocities across discrete time intervals $\Delta t$ using numerical integrators.

This document presents the mathematical formulations, derivations, and academic references for the physics and kernel methods in `mojospace`.

---

## 2. Smoothed Particle Hydrodynamics (SPH) Smoothing Kernels

In SPH, continuous continuum fields $A(\mathbf{r})$ (density, pressure, velocity) are approximated by a discrete summation over neighboring particles within smoothing length $h$:

$$A(\mathbf{r}) \approx \sum_j \frac{m_j}{\rho_j} A_j W(\|\mathbf{r} - \mathbf{r}_j\|, h)$$

where $W(r, h)$ is an isotropic smoothing kernel.

### 2.1 The Standard Cubic Spline Kernel (Monaghan 1992)

Let $q = \frac{r}{h}$, where $r = \|\mathbf{r}\|$. The cubic spline kernel is defined as:

$$W(r, h) = \sigma_d \begin{cases} 1 - \frac{3}{2}q^2 + \frac{3}{4}q^3, & 0 \le q \le 1 \\ \frac{1}{4}(2 - q)^3, & 1 < q \le 2 \\ 0, & q > 2 \end{cases}$$

#### Normalization Constant $\sigma_d$:
To satisfy the partition of unity $\int_{\mathbb{R}^d} W(r, h) \, d\mathbf{r} = 1$:
- **In 2D ($d=2$)**:
  $$\sigma_2 = \frac{10}{7\pi h^2}$$
- **In 3D ($d=3$)**:
  $$\sigma_3 = \frac{1}{\pi h^3}$$

#### First Derivative / Kernel Gradient:
$$\nabla W(\mathbf{r}, h) = \frac{\partial W}{\partial r} \frac{\mathbf{r}}{r}$$

$$\frac{\partial W}{\partial r} = \frac{\sigma_d}{h} \begin{cases} -3q + \frac{9}{4}q^2, & 0 \le q \le 1 \\ -\frac{3}{4}(2 - q)^2, & 1 < q \le 2 \\ 0, & q > 2 \end{cases}$$

---

### 2.2 SPH Spiky Pressure Kernel (Müller et al. 2003)

For pressure force computation in incompressible or weakly compressible fluids, the cubic spline's gradient approaches zero as $r \to 0$, causing particle clustering (the "tensile instability"). Müller et al. introduced the **Spiky kernel**, whose gradient does not vanish at $r \to 0$:

$$W_{\text{spiky}}(r, h) = \frac{15}{\pi h^6} (h - r)^3, \quad 0 \le r \le h$$

$$\nabla W_{\text{spiky}}(\mathbf{r}, h) = - \frac{45}{\pi h^6} (h - r)^2 \frac{\mathbf{r}}{r}, \quad 0 \le r \le h$$

---

## 3. Molecular & Electrostatic Potentials

### 3.1 Lennard-Jones (12-6) Potential

The Lennard-Jones potential models intermolecular interaction between neutral atoms:

$$V_{\text{LJ}}(r) = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 \right]$$

where:
- $\epsilon$ is the depth of the potential well.
- $\sigma$ is the distance at which the potential is zero ($r_{\text{min}} = 2^{1/6}\sigma$ is the equilibrium distance).

#### Force Vector:
$$\mathbf{F}_{\text{LJ}}(\mathbf{r}) = - \nabla V_{\text{LJ}}(r) = \frac{24\epsilon}{r^2} \left[ 2\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 \right] \mathbf{r}$$

For computational efficiency with cutoff $r_{\text{cut}} = 2.5\sigma$:
- Evaluate using $s = (\sigma^2 / r^2)^3$.
- Repulsive term $\propto s^2$, attractive dispersion term $\propto s$.

---

### 3.2 Coulomb Electrostatic Force (with Softening)

For charged particles with charges $q_1, q_2$:

$$\mathbf{F}_{\text{Coulomb}} = \frac{1}{4\pi\varepsilon_0} \frac{q_1 q_2}{(r^2 + \epsilon_e^2)^{3/2}} (\mathbf{r}_2 - \mathbf{r}_1)$$

where $\epsilon_e$ is the electrostatic softening length.

---

## 4. Symplectic Numerical Time Integrators

In Hamiltonian systems (such as orbital mechanics and molecular dynamics), non-symplectic methods (like standard Forward Euler or Runge-Kutta 4) exhibit artificial energy drift over long integration times.

### 4.1 Velocity Verlet Integrator (Second-Order Symplectic)

The Velocity Verlet algorithm is time-reversible, preserves phase space volume, and exhibits bounded energy error:

1. **Half-Step Velocity Update**:
   $$\mathbf{v}\left(t + \frac{\Delta t}{2}\right) = \mathbf{v}(t) + \frac{1}{2}\mathbf{a}(t)\Delta t$$

2. **Full-Step Position Update**:
   $$\mathbf{r}(t + \Delta t) = \mathbf{r}(t) + \mathbf{v}\left(t + \frac{\Delta t}{2}\right)\Delta t = \mathbf{r}(t) + \mathbf{v}(t)\Delta t + \frac{1}{2}\mathbf{a}(t)\Delta t^2$$

3. **Force & Acceleration Evaluation**:
   Compute new accelerations $\mathbf{a}(t + \Delta t) = \frac{\mathbf{F}(\mathbf{r}(t + \Delta t))}{m}$ using Spatial Hashing or Barnes-Hut.

4. **Second Half-Step Velocity Update**:
   $$\mathbf{v}(t + \Delta t) = \mathbf{v}\left(t + \frac{\Delta t}{2}\right) + \frac{1}{2}\mathbf{a}(t + \Delta t)\Delta t$$

---

### 4.2 Symplectic Euler Integrator (First-Order Symplectic)

A lightweight alternative for real-time physics engines:
$$\mathbf{v}(t + \Delta t) = \mathbf{v}(t) + \mathbf{a}(t)\Delta t$$
$$\mathbf{r}(t + \Delta t) = \mathbf{r}(t) + \mathbf{v}(t + \Delta t)\Delta t$$

---

## 5. Published Sources & Citations

1. **Monaghan, J. J.** (1992). *Smoothed Particle Hydrodynamics*. **Annual Review of Astronomy and Astrophysics**, 30(1), 543–574.
   - DOI: [10.1146/annurev.aa.30.090192.002551](https://doi.org/10.1146/annurev.aa.30.090192.002551)
2. **Müller, M., Charypar, D., & Gross, M.** (2003). *Particle-based fluid simulation for interactive applications*. In **Proceedings of the 2003 ACM SIGGRAPH/Eurographics symposium on Computer animation**, pp. 154–159.
   - DOI: [10.5555/846276.846298](https://doi.org/10.5555/846276.846298)
3. **Verlet, L.** (1967). *Computer "Experiments" on Classical Fluids. I. Thermodynamical Properties of Lennard-Jones Molecules*. **Physical Review**, 159(1), 98–103.
   - DOI: [10.1103/PhysRev.159.98](https://doi.org/10.1103/PhysRev.159.98)
4. **Hairer, E., Lubich, C., & Wanner, G.** (2006). *Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations*. Springer Series in Computational Mathematics, Vol. 31, Springer Berlin Heidelberg.
   - DOI: [10.1007/3-540-30666-8](https://doi.org/10.1007/3-540-30666-8)
