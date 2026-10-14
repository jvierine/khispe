//! MPI slab-decomposed 3-D incompressible Boussinesq KHI pilot.
//!
//! This is a deliberately compact first 3-D extension of the validated 2-D
//! morphology model.  It advances velocity with a projection method, a neutral
//! buoyancy perturbation, and a conservative ion-layer tracer.  HDF5 frames
//! include velocity, vorticity, and four materially advected periodic Clebsch
//! labels that exactly represent the *base-shear* vorticity at t=0 as the sum
//! of two Clebsch pairs.  They are diagnostics, not yet a vortex-particle
//! integrator.

use anyhow::{bail, ensure, Result};
use mpi::{collective::SystemOperation, topology::SimpleCommunicator, traits::*};
use serde::{Deserialize, Serialize};
use std::{f64::consts::PI, path::Path, time::Instant};

#[derive(Clone, Serialize, Deserialize)]
#[serde(default)]
struct Config {
    nx: usize,
    ny: usize,
    nz: usize,
    lx_m: f64,
    ly_m: f64,
    lz_m: f64,
    altitude_m: f64,
    shear_m: f64,
    wind_m_s: f64,
    advection_m_s: f64,
    n2_s2: f64,
    viscosity_m2_s: f64,
    diffusivity_m2_s: f64,
    dt_s: f64,
    end_s: f64,
    output_s: f64,
    seed_m_s: f64,
    seed_wavelength_m: f64,
    seed_pairing_m_s: f64,
    seed_pairing_phase_rad: f64,
    spanwise_seed_m_s: f64,
    spanwise_wavelength_m: f64,
    spanwise_phase_rad: f64,
    meridional_shear_center_m: f64,
    meridional_wind_fraction: f64,
    layer_sigma_m: f64,
    peak_m3: f64,
    poisson_rtol: f64,
    ion_mode: String,
    b_t: f64,
    dip_deg: f64,
    nu_in_s: f64,
    ion_mass_amu: f64,
    electric_east_v_m: f64,
    electric_north_v_m: f64,
    electric_up_v_m: f64,
    layer_descent_m_s: f64,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            nx: 96,
            ny: 32,
            nz: 64,
            lx_m: 86e3,
            ly_m: 32e3,
            lz_m: 6e3,
            altitude_m: 116e3,
            shear_m: 1060.,
            wind_m_s: 70.67,
            advection_m_s: -100.,
            n2_s2: 4e-4,
            viscosity_m2_s: 30.,
            diffusivity_m2_s: 5.,
            dt_s: 0.25,
            end_s: 120.,
            output_s: 30.,
            seed_m_s: 2.5,
            seed_wavelength_m: 21.5e3,
            seed_pairing_m_s: 1.5,
            seed_pairing_phase_rad: 0.7,
            spanwise_seed_m_s: 0.35,
            spanwise_wavelength_m: 32e3,
            spanwise_phase_rad: 0.35,
            meridional_shear_center_m: 0.,
            meridional_wind_fraction: 0.,
            layer_sigma_m: 255.,
            peak_m3: 1e11,
            poisson_rtol: 1e-8,
            ion_mode: "tracer".into(),
            b_t: 35e-6,
            dip_deg: 20.,
            nu_in_s: 300.,
            ion_mass_amu: 55.845,
            electric_east_v_m: 0.,
            electric_north_v_m: 0.,
            electric_up_v_m: 0.,
            layer_descent_m_s: 0.,
        }
    }
}

struct Grid<'a> {
    world: &'a SimpleCommunicator,
    nx: usize,
    ny: usize,
    nz: usize,
    start: usize,
    dx: f64,
    dy: f64,
    dz: f64,
}
impl Grid<'_> {
    fn len(&self) -> usize {
        self.nx * self.ny * self.nz
    }
    fn k(&self, i: usize, j: usize, z: usize) -> usize {
        (i * self.ny + j) * self.nz + z
    }
    fn sum(&self, x: f64) -> f64 {
        let mut y = 0.;
        self.world
            .all_reduce_into(&x, &mut y, SystemOperation::sum());
        y
    }
    fn max(&self, x: f64) -> f64 {
        let mut y = 0.;
        self.world
            .all_reduce_into(&x, &mut y, SystemOperation::max());
        y
    }
    fn mean0(&self, a: &mut [f64]) {
        let m = self.sum(a.iter().sum()) / (self.len() * self.world.size() as usize) as f64;
        for x in a {
            *x -= m;
        }
    }
    fn halo(&self, a: &[f64]) -> Vec<f64> {
        let slab = self.ny * self.nz;
        let mut h = vec![0.; (self.nx + 2) * slab];
        h[slab..(self.nx + 1) * slab].copy_from_slice(a);
        let r = self.world.rank();
        let n = self.world.size();
        let left = self.world.process_at_rank((r + n - 1) % n);
        let right = self.world.process_at_rank((r + 1) % n);
        mpi::point_to_point::send_receive_into(
            &a[(self.nx - 1) * slab..],
            &right,
            &mut h[..slab],
            &left,
        );
        mpi::point_to_point::send_receive_into(
            &a[..slab],
            &left,
            &mut h[(self.nx + 1) * slab..],
            &right,
        );
        h
    }
    fn lap(&self, a: &[f64]) -> Vec<f64> {
        let h = self.halo(a);
        let slab = self.ny * self.nz;
        let mut r = vec![0.; self.len()];
        for i in 0..self.nx {
            for j in 0..self.ny {
                for z in 0..self.nz {
                    let k = self.k(i, j, z);
                    let v = a[k];
                    let jm = (j + self.ny - 1) % self.ny;
                    let jp = (j + 1) % self.ny;
                    let dn = if z == 0 { v } else { a[k - 1] };
                    let up = if z + 1 == self.nz { v } else { a[k + 1] };
                    r[k] = (h[i * slab + j * self.nz + z] - 2. * v
                        + h[(i + 2) * slab + j * self.nz + z])
                        / self.dx.powi(2)
                        + (a[self.k(i, jm, z)] - 2. * v + a[self.k(i, jp, z)]) / self.dy.powi(2)
                        + (dn - 2. * v + up) / self.dz.powi(2);
                }
            }
        }
        r
    }
    fn dot(&self, a: &[f64], b: &[f64]) -> f64 {
        self.sum(a.iter().zip(b).map(|(x, y)| x * y).sum())
    }
    fn poisson_neumann(&self, rhs: &[f64], tol: f64) -> Result<Vec<f64>> {
        let mut b = rhs.to_vec();
        self.mean0(&mut b);
        let mut x = vec![0.; self.len()];
        let mut r = b.clone();
        let mut p = r.clone();
        let mut rr = self.dot(&r, &r);
        let threshold = tol.powi(2) * rr.max(1e-30);
        for _ in 0..5000 {
            if rr <= threshold {
                self.mean0(&mut x);
                return Ok(x);
            }
            let (gx, gy, gz) = self.pressure_gradient(&p);
            let ap: Vec<_> = self.div(&gx, &gy, &gz).into_iter().map(|q| -q).collect();
            let pap = self.dot(&p, &ap);
            ensure!(pap > 0., "3-D pressure CG lost definiteness");
            let a = rr / pap;
            for q in 0..x.len() {
                x[q] += a * p[q];
                r[q] -= a * ap[q];
            }
            self.mean0(&mut r);
            let next = self.dot(&r, &r);
            let beta = next / rr;
            for q in 0..p.len() {
                p[q] = r[q] + beta * p[q];
            }
            rr = next;
        }
        bail!("3-D pressure Poisson did not converge")
    }
    fn div(&self, u: &[f64], v: &[f64], w: &[f64]) -> Vec<f64> {
        let uh = self.halo(u);
        let slab = self.ny * self.nz;
        let mut d = vec![0.; self.len()];
        for i in 0..self.nx {
            for j in 0..self.ny {
                for z in 0..self.nz {
                    let k = self.k(i, j, z);
                    let jm = (j + self.ny - 1) % self.ny;
                    let wu = if z + 1 == self.nz { 0. } else { w[k + 1] };
                    // Face-to-cell divergence. Horizontal boundaries are
                    // periodic; w[0] is the lower wall face and the upper wall
                    // face outside the array is identically zero.
                    d[k] = (u[k] - uh[i * slab + j * self.nz + z]) / self.dx
                        + (v[k] - v[self.k(i, jm, z)]) / self.dy
                        + (wu - w[k]) / self.dz;
                }
            }
        }
        d
    }
    // G = -D^T for the face-to-cell divergence above. Thus -DG is the
    // nearest-neighbour Neumann Laplacian: there are no checkerboard pressure
    // null modes and the projection uses one compatible operator pair.
    fn pressure_gradient(&self, p: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let ph = self.halo(p);
        let slab = self.ny * self.nz;
        let mut gx = vec![0.; self.len()];
        let mut gy = gx.clone();
        let mut gz = gx.clone();
        for i in 0..self.nx {
            for j in 0..self.ny {
                for z in 0..self.nz {
                    let k = self.k(i, j, z);
                    let jp = self.k(i, (j + 1) % self.ny, z);
                    gx[k] = (ph[(i + 2) * slab + j * self.nz + z] - p[k]) / self.dx;
                    gy[k] = (p[jp] - p[k]) / self.dy;
                    if z > 0 {
                        gz[k] = (p[k] - p[k - 1]) / self.dz;
                    }
                }
            }
        }
        (gx, gy, gz)
    }
    fn project(
        &self,
        u: &mut [f64],
        v: &mut [f64],
        w: &mut [f64],
        dt: f64,
        tol: f64,
    ) -> Result<()> {
        // The lower wall face is explicit; the upper wall face is implicit.
        for i in 0..self.nx {
            for j in 0..self.ny {
                w[self.k(i, j, 0)] = 0.;
            }
        }
        let rhs: Vec<_> = self.div(u, v, w).into_iter().map(|q| -q / dt).collect();
        let p = self.poisson_neumann(&rhs, tol)?;
        let (gx, gy, gz) = self.pressure_gradient(&p);
        for k in 0..self.len() {
            u[k] -= dt * gx[k];
            v[k] -= dt * gy[k];
            w[k] -= dt * gz[k];
        }
        Ok(())
    }
    fn advect(&self, a: &[f64], u: &[f64], v: &[f64], w: &[f64]) -> Vec<f64> {
        let ah = self.halo(a);
        let uh = self.halo(u);
        let slab = self.ny * self.nz;
        let mut sx = vec![0.; self.len()];
        let mut sy = sx.clone();
        let mut sz = sx.clone();
        for i in 0..self.nx {
            for j in 0..self.ny {
                for z in 0..self.nz {
                    let k = self.k(i, j, z);
                    let jm = (j + self.ny - 1) % self.ny;
                    let jp = (j + 1) % self.ny;
                    let hc = i * slab + j * self.nz + z;
                    sx[k] = minmod(a[k] - ah[hc], ah[hc + 2 * slab] - a[k]);
                    sy[k] = minmod(a[k] - a[self.k(i, jm, z)], a[self.k(i, jp, z)] - a[k]);
                    if z > 0 && z + 1 < self.nz {
                        sz[k] = minmod(a[k] - a[k - 1], a[k + 1] - a[k]);
                    }
                }
            }
        }
        let sh = self.halo(&sx);
        let mut r = vec![0.; self.len()];
        for i in 0..self.nx {
            for j in 0..self.ny {
                for z in 0..self.nz {
                    let k = self.k(i, j, z);
                    let jm = self.k(i, (j + self.ny - 1) % self.ny, z);
                    let jp = self.k(i, (j + 1) % self.ny, z);
                    let hc = i * slab + j * self.nz + z;
                    let fl = flux(
                        ah[hc] + 0.5 * sh[hc],
                        a[k] - 0.5 * sx[k],
                        0.5 * (uh[hc] + u[k]),
                    );
                    let fr = flux(
                        a[k] + 0.5 * sx[k],
                        ah[hc + 2 * slab] - 0.5 * sh[hc + 2 * slab],
                        0.5 * (u[k] + uh[hc + 2 * slab]),
                    );
                    let fy0 = flux(
                        a[jm] + 0.5 * sy[jm],
                        a[k] - 0.5 * sy[k],
                        0.5 * (v[jm] + v[k]),
                    );
                    let fy1 = flux(
                        a[k] + 0.5 * sy[k],
                        a[jp] - 0.5 * sy[jp],
                        0.5 * (v[k] + v[jp]),
                    );
                    let fz0 = if z == 0 {
                        0.
                    } else {
                        flux(
                            a[k - 1] + 0.5 * sz[k - 1],
                            a[k] - 0.5 * sz[k],
                            0.5 * (w[k - 1] + w[k]),
                        )
                    };
                    let fz1 = if z + 1 == self.nz {
                        0.
                    } else {
                        flux(
                            a[k] + 0.5 * sz[k],
                            a[k + 1] - 0.5 * sz[k + 1],
                            0.5 * (w[k] + w[k + 1]),
                        )
                    };
                    r[k] = -(fr - fl) / self.dx - (fy1 - fy0) / self.dy - (fz1 - fz0) / self.dz;
                }
            }
        }
        r
    }
    fn gather(&self, a: &[f64]) -> Vec<f64> {
        if self.world.rank() == 0 {
            let mut q = vec![0.; a.len() * self.world.size() as usize];
            self.world.process_at_rank(0).gather_into_root(a, &mut q);
            q
        } else {
            self.world.process_at_rank(0).gather_into(a);
            Vec::new()
        }
    }
}
fn minmod(a: f64, b: f64) -> f64 {
    if a * b <= 0. {
        0.
    } else {
        a.signum() * a.abs().min(b.abs())
    }
}
fn flux(a: f64, b: f64, v: f64) -> f64 {
    0.5 * (v * (a + b) - v.abs() * (b - a))
}

fn ion_velocity(u: f64, v: f64, w: f64, c: &Config) -> [f64; 3] {
    if c.ion_mode == "tracer" {
        return [u, v, w];
    }
    let mu = 1.602176634e-19 / (c.ion_mass_amu * 1.66053906660e-27 * c.nu_in_s);
    let beta = mu * c.b_t;
    let by = c.dip_deg.to_radians().cos();
    let bz = -c.dip_deg.to_radians().sin();
    let f = [
        u + mu * c.electric_east_v_m,
        v + mu * c.electric_north_v_m,
        w + mu * c.electric_up_v_m,
    ];
    let dot = f[1] * by + f[2] * bz;
    let den = 1. + beta * beta;
    [
        (f[0] + beta * (f[1] * bz - f[2] * by)) / den,
        (f[1] - beta * f[0] * bz + beta * beta * by * dot) / den,
        (f[2] + beta * f[0] * by + beta * beta * bz * dot) / den,
    ]
}
#[derive(Clone)]
struct State {
    u: Vec<f64>,
    v: Vec<f64>,
    w: Vec<f64>,
    b: Vec<f64>,
    n: Vec<f64>,
    a1: Vec<f64>,
    b1: Vec<f64>,
    a2: Vec<f64>,
    b2: Vec<f64>,
}
fn rhs(g: &Grid, s: &State, c: &Config) -> State {
    let lu = g.lap(&s.u);
    let lv = g.lap(&s.v);
    let lw = g.lap(&s.w);
    let lb = g.lap(&s.b);
    let ln = g.lap(&s.n);
    let velocities: Vec<_> = (0..g.len())
        .map(|k| ion_velocity(s.u[k], s.v[k], s.w[k], c))
        .collect();
    let vi: Vec<_> = velocities.iter().map(|v| v[0]).collect();
    let vj: Vec<_> = velocities.iter().map(|v| v[1]).collect();
    let vk: Vec<_> = velocities.iter().map(|v| v[2]).collect();
    let mut r = State {
        u: g.advect(&s.u, &s.u, &s.v, &s.w),
        v: g.advect(&s.v, &s.u, &s.v, &s.w),
        w: g.advect(&s.w, &s.u, &s.v, &s.w),
        b: g.advect(&s.b, &s.u, &s.v, &s.w),
        n: g.advect(&s.n, &vi, &vj, &vk),
        a1: g.advect(&s.a1, &s.u, &s.v, &s.w),
        b1: g.advect(&s.b1, &s.u, &s.v, &s.w),
        a2: g.advect(&s.a2, &s.u, &s.v, &s.w),
        b2: g.advect(&s.b2, &s.u, &s.v, &s.w),
    };
    for k in 0..g.len() {
        r.u[k] += c.viscosity_m2_s * lu[k];
        r.v[k] += c.viscosity_m2_s * lv[k];
        r.w[k] += s.b[k] + c.viscosity_m2_s * lw[k];
        r.b[k] += -c.n2_s2 * s.w[k] + c.diffusivity_m2_s * lb[k];
        r.n[k] += c.diffusivity_m2_s * ln[k];
    }
    r
}
fn add(a: &State, r: &State, dt: f64) -> State {
    State {
        u: a.u.iter().zip(&r.u).map(|(x, y)| x + dt * y).collect(),
        v: a.v.iter().zip(&r.v).map(|(x, y)| x + dt * y).collect(),
        w: a.w.iter().zip(&r.w).map(|(x, y)| x + dt * y).collect(),
        b: a.b.iter().zip(&r.b).map(|(x, y)| x + dt * y).collect(),
        n: a.n.iter().zip(&r.n).map(|(x, y)| x + dt * y).collect(),
        a1: a.a1.iter().zip(&r.a1).map(|(x, y)| x + dt * y).collect(),
        b1: a.b1.iter().zip(&r.b1).map(|(x, y)| x + dt * y).collect(),
        a2: a.a2.iter().zip(&r.a2).map(|(x, y)| x + dt * y).collect(),
        b2: a.b2.iter().zip(&r.b2).map(|(x, y)| x + dt * y).collect(),
    }
}
fn avg(a: &State, b: &State) -> State {
    State {
        u: a.u.iter().zip(&b.u).map(|(x, y)| 0.5 * (x + y)).collect(),
        v: a.v.iter().zip(&b.v).map(|(x, y)| 0.5 * (x + y)).collect(),
        w: a.w.iter().zip(&b.w).map(|(x, y)| 0.5 * (x + y)).collect(),
        b: a.b.iter().zip(&b.b).map(|(x, y)| 0.5 * (x + y)).collect(),
        n: a.n.iter().zip(&b.n).map(|(x, y)| 0.5 * (x + y)).collect(),
        a1: a.a1.iter().zip(&b.a1).map(|(x, y)| 0.5 * (x + y)).collect(),
        b1: a.b1.iter().zip(&b.b1).map(|(x, y)| 0.5 * (x + y)).collect(),
        a2: a.a2.iter().zip(&b.a2).map(|(x, y)| 0.5 * (x + y)).collect(),
        b2: a.b2.iter().zip(&b.b2).map(|(x, y)| 0.5 * (x + y)).collect(),
    }
}
fn project(g: &Grid, s: &mut State, dt: f64, tol: f64) -> Result<()> {
    g.project(&mut s.u, &mut s.v, &mut s.w, dt, tol)
}

fn vorticity(g: &Grid, s: &State) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let vh = g.halo(&s.v);
    let wh = g.halo(&s.w);
    let slab = g.ny * g.nz;
    let mut ox = vec![0.; g.len()];
    let mut oy = ox.clone();
    let mut oz = ox.clone();
    for i in 0..g.nx {
        for j in 0..g.ny {
            for z in 0..g.nz {
                let k = g.k(i, j, z);
                let jm = (j + g.ny - 1) % g.ny;
                let jp = (j + 1) % g.ny;
                let zd = if z == 0 { z } else { z - 1 };
                let zu = if z + 1 == g.nz { z } else { z + 1 };
                let den = if z == 0 || z + 1 == g.nz {
                    g.dz
                } else {
                    2. * g.dz
                };
                ox[k] = (s.w[g.k(i, jp, z)] - s.w[g.k(i, jm, z)]) / (2. * g.dy)
                    - (s.v[g.k(i, j, zu)] - s.v[g.k(i, j, zd)]) / den;
                oy[k] = (s.u[g.k(i, j, zu)] - s.u[g.k(i, j, zd)]) / den
                    - (wh[(i + 2) * slab + j * g.nz + z] - wh[i * slab + j * g.nz + z])
                        / (2. * g.dx);
                oz[k] = (vh[(i + 2) * slab + j * g.nz + z] - vh[i * slab + j * g.nz + z])
                    / (2. * g.dx)
                    - (s.u[g.k(i, jp, z)] - s.u[g.k(i, jm, z)]) / (2. * g.dy);
            }
        }
    }
    (ox, oy, oz)
}
fn write_frame(g: &Grid, s: &State, c: &Config, path: &Path, t: f64) -> Result<()> {
    let (ox, oy, oz) = vorticity(g, s);
    let velocities: Vec<_> = (0..g.len())
        .map(|k| ion_velocity(s.u[k], s.v[k], s.w[k], c))
        .collect();
    let vi: Vec<_> = velocities.iter().map(|v| v[0]).collect();
    let vj: Vec<_> = velocities.iter().map(|v| v[1]).collect();
    let vk: Vec<_> = velocities.iter().map(|v| v[2]).collect();
    let local = [
        ("ion_u_m_s", &vi),
        ("ion_v_m_s", &vj),
        ("ion_w_m_s", &vk),
        ("u_m_s", &s.u),
        ("v_m_s", &s.v),
        ("w_m_s", &s.w),
        ("buoyancy_m_s2", &s.b),
        ("ion_density_m-3", &s.n),
        ("omega_x_s-1", &ox),
        ("omega_y_s-1", &oy),
        ("omega_z_s-1", &oz),
        ("clebsch_alpha1", &s.a1),
        ("clebsch_beta1_m2_s-1", &s.b1),
        ("clebsch_alpha2", &s.a2),
        ("clebsch_beta2_m2_s-1", &s.b2),
    ];
    let all: Vec<_> = local.iter().map(|(_, q)| g.gather(q)).collect();
    if g.world.rank() == 0 {
        let f = hdf5::File::create(path)?;
        f.new_attr::<f64>().create("time_s")?.write_scalar(&t)?;
        f.new_attr::<i32>()
            .create("mpi_ranks")?
            .write_scalar(&g.world.size())?;
        f.new_attr_builder().with_data("3-D incompressible Boussinesq pilot; Clebsch labels represent base vorticity at t=0").create("status")?;
        f.new_dataset_builder()
            .with_data(serde_json::to_string_pretty(c)?.as_bytes())
            .create("config_json")?;
        let x: Vec<_> = (0..c.nx).map(|i| (i as f64 + 0.5) * g.dx).collect();
        let y: Vec<_> = (0..c.ny)
            .map(|j| (j as f64 + 0.5) * g.dy - c.ly_m / 2.)
            .collect();
        let z: Vec<_> = (0..c.nz)
            .map(|k| c.altitude_m + c.layer_descent_m_s * t + (k as f64 + 0.5) * g.dz - c.lz_m / 2.)
            .collect();
        f.new_dataset_builder().with_data(&x).create("x_m")?;
        f.new_dataset_builder().with_data(&y).create("y_m")?;
        f.new_dataset_builder().with_data(&z).create("altitude_m")?;
        for ((name, _), q) in local.iter().zip(&all) {
            f.new_dataset::<f64>()
                .shape((c.nx, c.ny, c.nz))
                .deflate(1)
                .create(*name)?
                .write_raw(q)?;
        }
    }
    Ok(())
}
fn run(world: &SimpleCommunicator) -> Result<()> {
    let wall_started = Instant::now();
    let args: Vec<_> = std::env::args().collect();
    let c: Config = if args.len() > 1 {
        serde_json::from_slice(&std::fs::read(&args[1])?)?
    } else {
        Config::default()
    };
    let dir = Path::new(
        args.get(2)
            .map(String::as_str)
            .unwrap_or("output/khi3d-pilot"),
    );
    std::fs::create_dir_all(dir)?;
    ensure!(
        c.nx % world.size() as usize == 0 && c.ny >= 8 && c.nz >= 8,
        "nx must divide MPI ranks; ny,nz >= 8"
    );
    let g = Grid {
        world,
        nx: c.nx / world.size() as usize,
        ny: c.ny,
        nz: c.nz,
        start: world.rank() as usize * c.nx / world.size() as usize,
        dx: c.lx_m / c.nx as f64,
        dy: c.ly_m / c.ny as f64,
        dz: c.lz_m / c.nz as f64,
    };
    let mut s = State {
        u: vec![0.; g.len()],
        v: vec![0.; g.len()],
        w: vec![0.; g.len()],
        b: vec![0.; g.len()],
        n: vec![0.; g.len()],
        a1: vec![0.; g.len()],
        b1: vec![0.; g.len()],
        a2: vec![0.; g.len()],
        b2: vec![0.; g.len()],
    };
    let kx = 2. * PI / c.seed_wavelength_m;
    let kp = kx / 2.;
    let ky = 2. * PI / c.spanwise_wavelength_m;
    for i in 0..g.nx {
        for j in 0..g.ny {
            for q in 0..g.nz {
                let k = g.k(i, j, q);
                let x = (g.start + i) as f64 * g.dx;
                let y = (j as f64 + 0.5) * g.dy - c.ly_m / 2.;
                let z = (q as f64 + 0.5) * g.dz - c.lz_m / 2.;
                // A weak periodic cross-stream variation makes the shear layer
                // genuinely three-dimensional while remaining compatible with
                // the periodic y boundary.  The displaced coordinate also
                // initializes the ion sheet on the local shear center.
                let meridional_phase = ky * y;
                let shear_center = c.meridional_shear_center_m * meridional_phase.sin();
                let local_wind =
                    c.wind_m_s * (1. + c.meridional_wind_fraction * meridional_phase.cos());
                let zr = z - shear_center;
                let wall = (PI * (q as f64 + 0.5) / c.nz as f64).sin();
                let e = (-0.5 * (zr / (2. * c.shear_m)).powi(2)).exp() * wall;
                let phase = kx * x + ky * y + c.spanwise_phase_rad;
                let psi = e
                    * (c.seed_m_s / kx * (kx * x).sin()
                        + c.seed_pairing_m_s / kp * (kp * x + c.seed_pairing_phase_rad).sin()
                        + c.spanwise_seed_m_s / kx * phase.sin());
                let ep = e * (-zr / (4. * c.shear_m.powi(2)))
                    + (-0.5 * (zr / (2. * c.shear_m)).powi(2)).exp() * PI / c.lz_m
                        * (PI * (q as f64 + 0.5) / c.nz as f64).cos();
                s.u[k] = c.advection_m_s
                    + local_wind * (zr / c.shear_m).tanh()
                    + ep * (c.seed_m_s / kx * (kx * x).sin()
                        + c.seed_pairing_m_s / kp * (kp * x + c.seed_pairing_phase_rad).sin()
                        + c.spanwise_seed_m_s / kx * phase.sin());
                s.v[k] = 0.;
                s.w[k] = -e
                    * (c.seed_m_s * (kx * x).cos()
                        + c.seed_pairing_m_s * (kp * x + c.seed_pairing_phase_rad).cos()
                        + c.spanwise_seed_m_s * phase.cos());
                s.n[k] = c.peak_m3 * (-0.5 * (zr / c.layer_sigma_m).powi(2)).exp();
                let base = c.advection_m_s + local_wind * (zr / c.shear_m).tanh();
                s.a1[k] = (2. * PI * x / c.lx_m).sin();
                s.b1[k] = -base * (2. * PI * x / c.lx_m).cos() * c.lx_m / (2. * PI);
                s.a2[k] = (2. * PI * x / c.lx_m).cos();
                s.b2[k] = base * (2. * PI * x / c.lx_m).sin() * c.lx_m / (2. * PI);
                let _ = psi;
            }
        }
    }
    project(&g, &mut s, 1., c.poisson_rtol)?;
    let mass0 = g.sum(s.n.iter().sum());
    let mut t = 0.;
    let mut frame = 0;
    let mut next = 0.;
    if world.rank() == 0 {
        println!(
            "3-D KHI pilot {}x{}x{}, {} ranks",
            c.nx,
            c.ny,
            c.nz,
            world.size()
        );
    }
    while t < c.end_s + 1e-9 {
        if t >= next - 1e-9 {
            write_frame(&g, &s, &c, &dir.join(format!("frame_{frame:05}.h5")), t)?;
            let divergence = g.div(&s.u, &s.v, &s.w);
            let divergence_l2 = g.dot(&divergence, &divergence).sqrt();
            let mass = g.sum(s.n.iter().sum());
            if world.rank() == 0 {
                println!(
                    "frame {frame:05}: t={t:.1} s, wall={:.1} s, ||div u||_2={divergence_l2:.3e}, ion mass drift={:.3e}",
                    wall_started.elapsed().as_secs_f64(),
                    mass / mass0 - 1.
                );
            }
            frame += 1;
            next += c.output_s;
        }
        if t >= c.end_s - 1e-9 {
            break;
        }
        let umax =
            s.u.iter()
                .zip(&s.v)
                .zip(&s.w)
                .fold(0_f64, |m, ((u, v), w)| {
                    m.max(u.abs() / g.dx + v.abs() / g.dy + w.abs() / g.dz)
                });
        let ionmax = (0..g.len()).fold(0_f64, |m, k| {
            let v = ion_velocity(s.u[k], s.v[k], s.w[k], &c);
            m.max(v[0].abs() / g.dx + v[1].abs() / g.dy + v[2].abs() / g.dz)
        });
        let rate = g.max(umax.max(ionmax))
            + 2. * c.viscosity_m2_s * (g.dx.powi(-2) + g.dy.powi(-2) + g.dz.powi(-2));
        let dt = c.dt_s.min(0.3 / rate).min(next - t).min(c.end_s - t);
        let r1 = rhs(&g, &s, &c);
        let mut q = add(&s, &r1, dt);
        project(&g, &mut q, dt, c.poisson_rtol)?;
        let r2 = rhs(&g, &q, &c);
        let mut q2 = add(&q, &r2, dt);
        project(&g, &mut q2, dt, c.poisson_rtol)?;
        s = avg(&s, &q2);
        project(&g, &mut s, dt, c.poisson_rtol)?;
        t += dt;
        let bad =
            s.u.iter()
                .chain(&s.v)
                .chain(&s.w)
                .chain(&s.b)
                .chain(&s.n)
                .any(|v| !v.is_finite())
                || s.n.iter().any(|&v| v < -1e-8 * c.peak_m3);
        ensure!(
            g.max(if bad { 1. } else { 0. }) == 0.,
            "3-D nonfinite state or negative ion density"
        );
    }
    let divmax = g.max(
        g.div(&s.u, &s.v, &s.w)
            .into_iter()
            .fold(0_f64, |m, x| m.max(x.abs())),
    );
    let mass_error = (g.sum(s.n.iter().sum()) - mass0) / mass0;
    if world.rank() == 0 {
        std::fs::write(
            dir.join("diagnostics.json"),
            serde_json::to_string_pretty(
                &serde_json::json!({"time_s":t,"frames":frame,"wall_seconds":wall_started.elapsed().as_secs_f64(),"max_divergence_s-1":divmax,"mass_relative_error":mass_error,"status":"3-D pilot, not an observational fit"}),
            )?,
        )?;
        println!("done t={t:0.1}s max|div u|={divmax:0.3e}");
    }
    Ok(())
}
fn main() {
    let universe = mpi::initialize().expect("MPI initialize");
    let world = universe.world();
    if let Err(e) = run(&world) {
        eprintln!("rank {}: {e:#}", world.rank());
        world.abort(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn projection_removes_discrete_divergence() {
        let universe = mpi::initialize().expect("MPI initialize");
        let world = universe.world();
        let g = Grid {
            world: &world,
            nx: 16,
            ny: 4,
            nz: 12,
            start: 0,
            dx: 2.,
            dy: 3.,
            dz: 1.,
        };
        let mut u = vec![0.; g.len()];
        let mut v = u.clone();
        let mut w = u.clone();
        for i in 0..g.nx {
            for j in 0..g.ny {
                for z in 0..g.nz {
                    let k = g.k(i, j, z);
                    u[k] = (2. * PI * i as f64 / 4.).sin() * (1. + z as f64 / 12.);
                    v[k] = (2. * PI * j as f64 / 4.).cos();
                    if z > 0 && z + 1 < g.nz {
                        w[k] = (z as f64).sin();
                    }
                }
            }
        }
        let before = g.dot(&g.div(&u, &v, &w), &g.div(&u, &v, &w)).sqrt();
        g.project(&mut u, &mut v, &mut w, 0.1, 1e-10).unwrap();
        let after = g.dot(&g.div(&u, &v, &w), &g.div(&u, &v, &w)).sqrt();
        assert!(
            after / before < 1e-8,
            "residual fraction {}",
            after / before
        );
        for i in 0..g.nx {
            for j in 0..g.ny {
                assert_eq!(w[g.k(i, j, 0)], 0.);
            }
        }
    }
    #[test]
    fn magnetic_reference_flow() {
        let mut c = Config::default();
        c.ion_mode = "collisional".into();
        let dip = c.dip_deg.to_radians();
        c.electric_north_v_m = c.advection_m_s * (-c.b_t * dip.sin());
        c.electric_up_v_m = -c.advection_m_s * c.b_t * dip.cos();
        let vi = ion_velocity(c.advection_m_s, 0., 0., &c);
        assert!((vi[0] - c.advection_m_s).abs() < 1e-10);
        assert!(vi[1].abs() < 1e-10 && vi[2].abs() < 1e-10);
        assert!(ion_velocity(c.advection_m_s + 20., 0., 0., &c)[2] > 0.);
    }
    #[test]
    fn momentum_balance() {
        let mut c = Config::default();
        c.ion_mode = "collisional".into();
        let u = [50., 3., 2.];
        let v = ion_velocity(u[0], u[1], u[2], &c);
        let omega = 1.602176634e-19 * c.b_t / (c.ion_mass_amu * 1.66053906660e-27);
        let by = c.dip_deg.to_radians().cos();
        let bz = -c.dip_deg.to_radians().sin();
        let cross = [v[1] * bz - v[2] * by, -v[0] * bz, v[0] * by];
        for k in 0..3 {
            assert!((c.nu_in_s * (v[k] - u[k]) - omega * cross[k]).abs() < 1e-9);
        }
    }
    #[test]
    fn limiter() {
        assert_eq!(minmod(2., 3.), 2.);
        assert_eq!(minmod(-2., 3.), 0.);
    }
}
