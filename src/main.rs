use anyhow::{Result, bail, ensure};
use mpi::{collective::SystemOperation, topology::SimpleCommunicator, traits::*};
use serde::{Deserialize, Serialize};
use std::{f64::consts::PI, path::Path};

#[derive(Clone, Serialize, Deserialize)]
#[serde(default)]
struct Config {
    nx: usize,
    nz: usize,
    lx_m: f64,
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
    seed_domain_m_s: f64,
    seed_domain_phase_rad: f64,
    ion_mode: String,
    b_t: f64,
    dip_deg: f64,
    nu_in_s: f64,
    ion_mass_amu: f64,
    layer_sigma_m: f64,
    peak_m3: f64,
    poisson_rtol: f64,
    validate_only: bool,
    restart_hdf5: Option<String>,
    // Prescribed phase motion of the entire layer in a translating vertical
    // coordinate. This is a kinematic pilot, not neutral vertical wind or
    // self-consistent electrodynamic ion convergence.
    layer_descent_m_s: f64,
    // Phenomenological metallic-ion trapping rate. In the translating layer
    // frame v_i,z = v_response,z - kappa*z. No chemical density loss.
    ion_convergence_s: f64,
    electric_east_v_m: f64,
    electric_north_v_m: f64,
    electric_up_v_m: f64,
    seed_packet_center_m: f64,
    seed_packet_half_width_m: f64,
    seed_packet_taper_m: f64,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            nx: 128,
            nz: 128,
            lx_m: 12000.,
            lz_m: 6000.,
            altitude_m: 115000.,
            shear_m: 350.,
            wind_m_s: 50.,
            advection_m_s: -30.,
            n2_s2: 4e-4,
            viscosity_m2_s: 30.,
            diffusivity_m2_s: 5.,
            dt_s: 0.25,
            end_s: 600.,
            output_s: 10.,
            seed_m_s: 0.5,
            seed_wavelength_m: 12000.,
            seed_pairing_m_s: 0.,
            seed_pairing_phase_rad: 0.,
            seed_domain_m_s: 0.,
            seed_domain_phase_rad: 0.,
            ion_mode: "collisional".into(),
            b_t: 3.5e-5,
            dip_deg: 20.,
            nu_in_s: 300.,
            ion_mass_amu: 55.845,
            layer_sigma_m: 300.,
            peak_m3: 1e11,
            poisson_rtol: 1e-9,
            validate_only: false,
            restart_hdf5: None,
            layer_descent_m_s: 0.,
            ion_convergence_s: 0.,
            electric_east_v_m: 0.,
            electric_north_v_m: 0.,
            electric_up_v_m: 0.,
            seed_packet_center_m: 0.,
            seed_packet_half_width_m: 0.,
            seed_packet_taper_m: 0.,
        }
    }
}
struct Grid<'a> {
    world: &'a SimpleCommunicator,
    nx: usize,
    nz: usize,
    start: usize,
    dx: f64,
    dz: f64,
}
impl Grid<'_> {
    fn len(&self) -> usize {
        self.nx * self.nz
    }
    fn halo(&self, a: &[f64]) -> Vec<f64> {
        let mut h = vec![0.; (self.nx + 2) * self.nz];
        h[self.nz..(self.nx + 1) * self.nz].copy_from_slice(a);
        let rank = self.world.rank();
        let size = self.world.size();
        let left = self.world.process_at_rank((rank + size - 1) % size);
        let right = self.world.process_at_rank((rank + 1) % size);
        mpi::point_to_point::send_receive_into(
            &a[(self.nx - 1) * self.nz..],
            &right,
            &mut h[..self.nz],
            &left,
        );
        mpi::point_to_point::send_receive_into(
            &a[..self.nz],
            &left,
            &mut h[(self.nx + 1) * self.nz..],
            &right,
        );
        h
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
    fn dot(&self, a: &[f64], b: &[f64]) -> f64 {
        self.sum(a.iter().zip(b).map(|(a, b)| a * b).sum())
    }
    fn lap(&self, a: &[f64], neumann: bool) -> Vec<f64> {
        let h = self.halo(a);
        let mut r = vec![0.; self.len()];
        for i in 0..self.nx {
            for j in 0..self.nz {
                let k = i * self.nz + j;
                let v = a[k];
                let down = if j == 0 {
                    if neumann { v } else { 0. }
                } else {
                    a[k - 1]
                };
                let up = if j + 1 == self.nz {
                    if neumann { v } else { 0. }
                } else {
                    a[k + 1]
                };
                r[k] = (h[i * self.nz + j] - 2. * v + h[(i + 2) * self.nz + j]) / self.dx.powi(2)
                    + (down - 2. * v + up) / self.dz.powi(2);
            }
        }
        r
    }
    // SPD -Laplacian, periodic x and homogeneous Dirichlet streamfunction at z walls.
    fn poisson(&self, omega: &[f64], guess: &[f64], tol: f64) -> Result<Vec<f64>> {
        let mut x = guess.to_vec();
        let ax = self.lap(&x, false);
        let mut r: Vec<_> = omega.iter().zip(ax).map(|(b, a)| b + a).collect();
        let mut p = r.clone();
        let mut rr = self.dot(&r, &r);
        let threshold = tol.powi(2) * self.dot(omega, omega).max(1e-30);
        for _ in 0..4000 {
            if rr <= threshold {
                return Ok(x);
            }
            let ap: Vec<_> = self.lap(&p, false).into_iter().map(|a| -a).collect();
            let pap = self.dot(&p, &ap);
            ensure!(pap > 0., "Poisson CG lost positive definiteness");
            let alpha = rr / pap;
            for k in 0..x.len() {
                x[k] += alpha * p[k];
                r[k] -= alpha * ap[k];
            }
            let next = self.dot(&r, &r);
            let beta = next / rr;
            for k in 0..x.len() {
                p[k] = r[k] + beta * p[k];
            }
            rr = next;
        }
        bail!("Poisson CG did not converge")
    }
    fn velocity(&self, psi: &[f64], c: &Config) -> (Vec<f64>, Vec<f64>) {
        let h = self.halo(psi);
        let mut u = vec![0.; self.len()];
        let mut w = u.clone();
        for i in 0..self.nx {
            for j in 0..self.nz {
                let k = i * self.nz + j;
                let z = (j as f64 + 1.) * self.dz - c.lz_m / 2.;
                let down = if j == 0 { 0. } else { psi[k - 1] };
                let up = if j + 1 == self.nz { 0. } else { psi[k + 1] };
                u[k] = c.advection_m_s
                    + c.wind_m_s * (z / c.shear_m).tanh()
                    + (up - down) / (2. * self.dz);
                w[k] = -(h[(i + 2) * self.nz + j] - h[i * self.nz + j]) / (2. * self.dx);
            }
        }
        (u, w)
    }
    // Conservative second-order MUSCL/Rusanov scalar fluxes. Reflecting/no-flux
    // z walls; two-stage SSP RK2 applies the same flux operator to ion density.
    fn transport(&self, a: &[f64], u: &[f64], w: &[f64]) -> Vec<f64> {
        let h = self.halo(a);
        let uh = self.halo(u);
        let mut sx = vec![0.; self.len()];
        let mut sz = sx.clone();
        for i in 0..self.nx {
            for j in 0..self.nz {
                let k = i * self.nz + j;
                sx[k] = minmod(a[k] - h[i * self.nz + j], h[(i + 2) * self.nz + j] - a[k]);
                if j > 0 && j + 1 < self.nz {
                    sz[k] = minmod(a[k] - a[k - 1], a[k + 1] - a[k]);
                }
            }
        }
        let sh = self.halo(&sx);
        let mut fx = vec![0.; (self.nx + 1) * self.nz];
        let mut fz = vec![0.; self.nx * (self.nz + 1)];
        for i in 0..=self.nx {
            for j in 0..self.nz {
                let l = i * self.nz + j;
                let r = (i + 1) * self.nz + j;
                fx[l] = flux(
                    h[l] + 0.5 * sh[l],
                    h[r] - 0.5 * sh[r],
                    0.5 * (uh[l] + uh[r]),
                );
            }
        }
        for i in 0..self.nx {
            for j in 1..self.nz {
                let k = i * self.nz + j;
                fz[i * (self.nz + 1) + j] = flux(
                    a[k - 1] + 0.5 * sz[k - 1],
                    a[k] - 0.5 * sz[k],
                    0.5 * (w[k - 1] + w[k]),
                );
            }
        }
        let mut rhs = vec![0.; self.len()];
        for i in 0..self.nx {
            for j in 0..self.nz {
                let k = i * self.nz + j;
                rhs[k] = -(fx[(i + 1) * self.nz + j] - fx[k]) / self.dx
                    - (fz[i * (self.nz + 1) + j + 1] - fz[i * (self.nz + 1) + j]) / self.dz;
            }
        }
        rhs
    }
    fn gather(&self, a: &[f64]) -> Vec<f64> {
        if self.world.rank() == 0 {
            let mut all = vec![0.; a.len() * self.world.size() as usize];
            self.world.process_at_rank(0).gather_into_root(a, &mut all);
            all
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
fn packet_window(x: f64, c: &Config) -> f64 {
    if c.seed_packet_half_width_m <= 0. {
        return 1.;
    }
    let d = (x - c.seed_packet_center_m).abs() - c.seed_packet_half_width_m;
    if d <= 0. {
        1.
    } else if d >= c.seed_packet_taper_m || c.seed_packet_taper_m <= 0. {
        0.
    } else {
        0.5 * (1. + (PI * d / c.seed_packet_taper_m).cos())
    }
}
fn ion_velocity(u: f64, w: f64, c: &Config) -> (f64, f64) {
    if c.ion_mode == "tracer" {
        return (u, w);
    }
    // nu(v-u)=(q/m)E+Omega v x b. The y component participates
    // algebraically even though this initial spatial domain is an x-z slice.
    // x=east, y=north, z=up; B=(0, B cos I, -B sin I).
    let mobility = 1.602176634e-19 / (c.ion_mass_amu * 1.66053906660e-27 * c.nu_in_s);
    let beta = mobility * c.b_t;
    let by = c.dip_deg.to_radians().cos();
    let bz = -c.dip_deg.to_radians().sin();
    let fx = u + mobility * c.electric_east_v_m;
    let fy = mobility * c.electric_north_v_m;
    let fz = w + mobility * c.electric_up_v_m;
    let dot = fy * by + fz * bz;
    let cross_x = fy * bz - fz * by;
    let cross_z = fx * by;
    (
        (fx + beta * cross_x) / (1. + beta * beta),
        (fz + beta * cross_z + beta * beta * bz * dot) / (1. + beta * beta),
    )
}
struct State {
    omega: Vec<f64>,
    b: Vec<f64>,
    n: Vec<f64>,
    psi: Vec<f64>,
}

fn read_restart(g: &Grid, c: &Config, path: &Path) -> Result<(State, f64)> {
    let file = hdf5::File::open(path)?;
    let time = file.attr("time_s")?.read_scalar::<f64>()?;
    let read_local = |name: &str| -> Result<Vec<f64>> {
        let dataset = file.dataset(name)?;
        ensure!(
            dataset.shape() == [c.nx, c.nz],
            "restart dataset {name} has shape {:?}, expected [{}, {}]",
            dataset.shape(),
            c.nx,
            c.nz
        );
        let global = dataset.read_raw::<f64>()?;
        let first = g.start * g.nz;
        let last = (g.start + g.nx) * g.nz;
        Ok(global[first..last].to_vec())
    };
    let omega = read_local("omega_s-1")?;
    let b = read_local("buoyancy_m_s2")?;
    let n = read_local("ion_density_m-3")?;
    ensure!(
        omega.iter().chain(&b).chain(&n).all(|x| x.is_finite()),
        "restart contains nonfinite values"
    );
    ensure!(
        n.iter().all(|x| *x >= -1e-8 * c.peak_m3),
        "restart contains negative density"
    );
    let psi = g.poisson(&omega, &vec![0.; g.len()], c.poisson_rtol)?;
    Ok((State { omega, b, n, psi }, time))
}

fn next_frame_number(dir: &Path) -> usize {
    std::fs::read_dir(dir)
        .into_iter()
        .flatten()
        .filter_map(|entry| entry.ok())
        .filter_map(|entry| {
            let name = entry.file_name();
            let name = name.to_str()?;
            name.strip_prefix("frame_")?
                .strip_suffix(".h5")?
                .parse::<usize>()
                .ok()
        })
        .max()
        .map_or(0, |n| n + 1)
}
fn rhs(g: &Grid, s: &State, c: &Config) -> Result<(State, Vec<f64>, Vec<f64>, f64)> {
    let psi = g.poisson(&s.omega, &s.psi, c.poisson_rtol)?;
    let (u, w) = g.velocity(&psi, c);
    let mut total = s.omega.clone();
    for i in 0..g.nx {
        for j in 0..g.nz {
            let z = (j as f64 + 1.) * g.dz - c.lz_m / 2.;
            total[i * g.nz + j] -= c.wind_m_s / c.shear_m / (z / c.shear_m).cosh().powi(2);
        }
    }
    let mut ro = g.transport(&total, &u, &w);
    let mut rb = g.transport(&s.b, &u, &w);
    let mut vi = vec![0.; g.len()];
    let mut wi = vi.clone();
    for k in 0..g.len() {
        (vi[k], wi[k]) = ion_velocity(u[k], w[k], c);
        let z = ((k % g.nz) as f64 + 1.) * g.dz - c.lz_m / 2.;
        wi[k] -= c.ion_convergence_s * z;
    }
    let mut rn = g.transport(&s.n, &vi, &wi);
    let lo = g.lap(&total, true);
    let lb = g.lap(&s.b, true);
    let ln = g.lap(&s.n, true);
    let bh = g.halo(&s.b);
    let mut rate: f64 = 0.;
    for i in 0..g.nx {
        for j in 0..g.nz {
            let k = i * g.nz + j;
            // curl(b e_z) in our omega=dx(w)-dz(u) convention is +dx(b).
            ro[k] += (bh[(i + 2) * g.nz + j] - bh[i * g.nz + j]) / (2. * g.dx)
                + c.viscosity_m2_s * lo[k];
            rb[k] += -c.n2_s2 * w[k] + c.diffusivity_m2_s * lb[k];
            rn[k] += c.diffusivity_m2_s * ln[k];
            rate = rate
                .max(u[k].abs() / g.dx + w[k].abs() / g.dz)
                .max(vi[k].abs() / g.dx + wi[k].abs() / g.dz);
        }
    }
    rate += 2. * c.viscosity_m2_s.max(c.diffusivity_m2_s) * (1. / g.dx.powi(2) + 1. / g.dz.powi(2));
    Ok((
        State {
            omega: ro,
            b: rb,
            n: rn,
            psi: psi.clone(),
        },
        u,
        w,
        g.max(rate),
    ))
}
fn write_frame(g: &Grid, s: &State, c: &Config, path: &Path, t: f64) -> Result<()> {
    let (_, u, w, _) = rhs(g, s, c)?;
    let fields = [
        ("omega_s-1", &s.omega),
        ("buoyancy_m_s2", &s.b),
        ("ion_density_m-3", &s.n),
        ("u_m_s", &u),
        ("w_m_s", &w),
    ];
    let all: Vec<_> = fields.iter().map(|(_, a)| g.gather(a)).collect();
    let mut error = None;
    if g.world.rank() == 0 {
        let result = (|| -> Result<()> {
            let file = hdf5::File::create(path)?;
            file.new_attr::<f64>().create("time_s")?.write_scalar(&t)?;
            file.new_attr::<f64>()
                .create("layer_center_altitude_m")?
                .write_scalar(&(c.altitude_m + c.layer_descent_m_s * t))?;
            file.new_attr::<i32>()
                .create("mpi_ranks")?
                .write_scalar(&g.world.size())?;
            file.new_dataset_builder()
                .with_data(serde_json::to_string_pretty(c)?.as_bytes())
                .create("config_json")?;
            let x: Vec<_> = (0..c.nx).map(|i| (i as f64 + 0.5) * g.dx).collect();
            let z: Vec<_> = (0..c.nz)
                .map(|j| {
                    c.altitude_m + c.layer_descent_m_s * t + (j as f64 + 1.) * g.dz - c.lz_m / 2.
                })
                .collect();
            file.new_dataset_builder().with_data(&x).create("x_m")?;
            file.new_dataset_builder()
                .with_data(&z)
                .create("altitude_m")?;
            for ((name, _), a) in fields.iter().zip(&all) {
                file.new_dataset::<f64>()
                    .shape((c.nx, c.nz))
                    .create(*name)?
                    .write_raw(a)?;
            }
            Ok(())
        })();
        if let Err(e) = result {
            error = Some(e)
        }
    }
    let mut failed = i32::from(error.is_some());
    g.world.process_at_rank(0).broadcast_into(&mut failed);
    if failed != 0 {
        bail!("HDF5 output failed: {:?}", error)
    }
    Ok(())
}
fn run(world: &SimpleCommunicator) -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("--compare") {
        ensure!(args.len() == 4, "--compare frame_a.h5 frame_b.h5");
        let a = hdf5::File::open(&args[2])?;
        let b = hdf5::File::open(&args[3])?;
        for name in [
            "omega_s-1",
            "buoyancy_m_s2",
            "ion_density_m-3",
            "u_m_s",
            "w_m_s",
        ] {
            let ad = a.dataset(name)?;
            let bd = b.dataset(name)?;
            ensure!(ad.shape() == bd.shape(), "shape mismatch");
            let x = ad.read_raw::<f64>()?;
            let y = bd.read_raw::<f64>()?;
            let scale = x.iter().chain(&y).fold(1e-15_f64, |v, x| v.max(x.abs()));
            let diff = x
                .iter()
                .zip(&y)
                .fold(0_f64, |v, (x, y)| v.max((x - y).abs()))
                / scale;
            ensure!(diff < 1e-7, "rank comparison failed for {name}: {diff}");
            if world.rank() == 0 {
                println!("{name}: max relative difference={diff:.3e}");
            }
        }
        return Ok(());
    }
    let c: Config = if args.len() > 1 {
        serde_json::from_slice(&std::fs::read(&args[1])?)?
    } else {
        Config::default()
    };
    ensure!(
        c.nx >= 8 && c.nz >= 8 && c.nx % world.size() as usize == 0,
        "nx must be >=8 and divisible by MPI ranks"
    );
    ensure!(
        c.nx / world.size() as usize >= 2,
        "at least two x columns per rank required"
    );
    ensure!(
        c.lx_m > 0.
            && c.lz_m > 0.
            && c.shear_m > 0.
            && c.dt_s > 0.
            && c.end_s >= 0.
            && c.output_s > 0.
            && c.nu_in_s > 0.
            && c.ion_mass_amu > 0.
            && c.poisson_rtol > 0.
            && c.layer_sigma_m > 0.
            && c.seed_wavelength_m > 0.
            && c.peak_m3 > 0.,
        "invalid positive parameter"
    );
    ensure!(
        c.n2_s2 >= 0.
            && c.ion_convergence_s >= 0.
            && c.viscosity_m2_s >= 0.
            && c.diffusivity_m2_s >= 0.
            && ["tracer", "collisional"].contains(&c.ion_mode.as_str()),
        "invalid model option"
    );
    let dir = Path::new(args.get(2).map(String::as_str).unwrap_or("output/pilot"));
    std::fs::create_dir_all(dir)?;
    let g = Grid {
        world,
        nx: c.nx / world.size() as usize,
        nz: c.nz,
        start: world.rank() as usize * c.nx / world.size() as usize,
        dx: c.lx_m / c.nx as f64,
        dz: c.lz_m / (c.nz + 1) as f64,
    };
    if c.validate_only {
        let exact: Vec<f64> = (0..g.len())
            .map(|k| {
                let i = k / g.nz;
                let j = k % g.nz;
                (2. * PI * (g.start + i) as f64 / c.nx as f64).sin()
                    * (PI * (j + 1) as f64 / (c.nz + 1) as f64).sin()
            })
            .collect();
        let discrete = g
            .lap(&exact, false)
            .into_iter()
            .map(|x| -x)
            .collect::<Vec<_>>();
        let solved = g.poisson(&discrete, &vec![0.; g.len()], c.poisson_rtol)?;
        let err = g.max(
            exact
                .iter()
                .zip(&solved)
                .fold(0_f64, |v, (a, b)| v.max((a - b).abs())),
        );
        ensure!(err < 1e-7, "manufactured discrete Poisson error {err}");
        let lambda = (2. * PI / c.lx_m).powi(2) + (PI / c.lz_m).powi(2);
        let continuum = exact.iter().map(|x| x * lambda).collect::<Vec<_>>();
        let recovered = g.poisson(&continuum, &vec![0.; g.len()], c.poisson_rtol)?;
        let truncation = g.max(
            exact
                .iter()
                .zip(&recovered)
                .fold(0_f64, |v, (a, b)| v.max((a - b).abs())),
        );
        let transported = g.transport(&exact, &vec![30.; g.len()], &vec![0.; g.len()]);
        ensure!(
            g.sum(transported.iter().sum()).abs() < 1e-10,
            "flux conservation failed"
        );
        if world.rank() == 0 {
            let report = serde_json::json!({"nx":c.nx,"nz":c.nz,"ranks":world.size(),"discrete_poisson_max_error":err,"continuum_poisson_max_error":truncation});
            std::fs::write(
                dir.join("validation.json"),
                serde_json::to_string_pretty(&report)?,
            )?;
            println!("{report}");
        }
        return Ok(());
    }
    let mut psi = vec![0.; g.len()];
    let mut n = psi.clone();
    for i in 0..g.nx {
        for j in 0..g.nz {
            let x = (g.start + i) as f64 * g.dx;
            let z = (j as f64 + 1.) * g.dz - c.lz_m / 2.;
            let k = i * g.nz + j;
            let envelope = (-0.5 * (z / (2. * c.shear_m)).powi(2)).exp()
                * (PI * (j as f64 + 1.) / (c.nz + 1) as f64).sin();
            let pairing_wavelength = 2. * c.seed_wavelength_m;
            psi[k] = envelope
                * packet_window(x, &c)
                * (c.seed_m_s * c.seed_wavelength_m / (2. * PI)
                    * (2. * PI * x / c.seed_wavelength_m).sin()
                    + c.seed_pairing_m_s * pairing_wavelength / (2. * PI)
                        * (2. * PI * x / pairing_wavelength + c.seed_pairing_phase_rad).sin()
                    + c.seed_domain_m_s * c.lx_m / (2. * PI)
                        * (2. * PI * x / c.lx_m + c.seed_domain_phase_rad).sin());
            n[k] = c.peak_m3 * (-0.5 * (z / c.layer_sigma_m).powi(2)).exp();
        }
    }
    let omega = g.lap(&psi, false).into_iter().map(|x| -x).collect();
    let mut s = State {
        omega,
        b: vec![0.; g.len()],
        n,
        psi,
    };
    let mut time = 0.;
    let mut frame = 0;
    let mut next = 0.;
    if let Some(restart) = c.restart_hdf5.as_deref() {
        (s, time) = read_restart(&g, &c, Path::new(restart))?;
        ensure!(
            c.end_s > time,
            "continuation end_s={} must exceed restart time_s={time}",
            c.end_s
        );
        frame = next_frame_number(dir);
        next = time + c.output_s;
        if world.rank() == 0 {
            println!("restarting from {restart} at t={time:.3} s; next frame={frame}");
        }
    }
    let mass0 = g.sum(s.n.iter().sum());
    let start = std::time::Instant::now();
    if world.rank() == 0 {
        println!(
            "2D x-z Boussinesq KHI pilot: {}x{}, {} MPI ranks; Ri_center={:.4}; ion mode={}",
            c.nx,
            c.nz,
            world.size(),
            c.n2_s2 / (c.wind_m_s / c.shear_m).powi(2),
            c.ion_mode
        );
    }
    loop {
        if time >= next - 1e-9 || time >= c.end_s - 1e-9 {
            write_frame(&g, &s, &c, &dir.join(format!("frame_{frame:05}.h5")), time)?;
            frame += 1;
            next += c.output_s;
        }
        if time >= c.end_s - 1e-9 {
            break;
        }
        let (a, _, _, rate) = rhs(&g, &s, &c)?;
        let dt = c
            .dt_s
            .min(0.35 / rate.max(1e-30))
            .min(0.2 / c.n2_s2.sqrt().max(1e-30))
            .min(c.end_s - time)
            .min(next - time);
        ensure!(dt > 1e-10, "timestep collapsed");
        let stage = State {
            omega: s
                .omega
                .iter()
                .zip(&a.omega)
                .map(|(x, r)| x + dt * r)
                .collect(),
            b: s.b.iter().zip(&a.b).map(|(x, r)| x + dt * r).collect(),
            n: s.n.iter().zip(&a.n).map(|(x, r)| x + dt * r).collect(),
            psi: a.psi,
        };
        let (b, _, _, rate2) = rhs(&g, &stage, &c)?;
        ensure!(dt * rate2 <= 0.5, "stage CFL exceeded; reduce dt_s");
        for k in 0..g.len() {
            s.omega[k] = 0.5 * (s.omega[k] + stage.omega[k] + dt * b.omega[k]);
            s.b[k] = 0.5 * (s.b[k] + stage.b[k] + dt * b.b[k]);
            s.n[k] = 0.5 * (s.n[k] + stage.n[k] + dt * b.n[k]);
        }
        s.psi = b.psi;
        time += dt;
        let bad = s
            .omega
            .iter()
            .chain(&s.b)
            .chain(&s.n)
            .any(|x| !x.is_finite())
            || s.n.iter().any(|x| *x < -1e-8 * c.peak_m3);
        ensure!(
            g.max(if bad { 1. } else { 0. }) == 0.,
            "nonfinite state or negative density"
        );
    }
    let mass = g.sum(s.n.iter().sum());
    let (_, _, w, _) = rhs(&g, &s, &c)?;
    let rms = (g.sum(w.iter().map(|x| x * x).sum()) / (c.nx * c.nz) as f64).sqrt();
    if world.rank() == 0 {
        let report = serde_json::json!({"time_s":time,"ranks":world.size(),"mass_relative_error":(mass-mass0)/mass0,"vertical_wind_rms_m_s":rms,"wall_seconds":start.elapsed().as_secs_f64(),"frames":frame,"status":"pilot_not_observationally_validated"});
        std::fs::write(
            dir.join("diagnostics.json"),
            serde_json::to_string_pretty(&report)?,
        )?;
        println!("{report}");
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
    fn limiter() {
        assert_eq!(minmod(1., -2.), 0.);
        assert_eq!(minmod(-3., -2.), -2.);
    }
    #[test]
    fn upwind_flux() {
        assert_eq!(flux(2., 9., 3.), 6.);
        assert_eq!(flux(2., 9., -3.), -27.);
    }
    #[test]
    fn ion_force_balance() {
        let c = Config::default();
        let (vx, vz) = ion_velocity(50., 2., &c);
        let beta = 1.602176634e-19 * c.b_t / (c.ion_mass_amu * 1.66053906660e-27 * c.nu_in_s);
        let by = c.dip_deg.to_radians().cos();
        let bz = -c.dip_deg.to_radians().sin();
        let vy = (-beta * 50. * bz + beta * beta * by * bz * 2.) / (1. + beta * beta);
        assert!((vx - 50. - beta * (vy * bz - vz * by)).abs() < 1e-12);
        assert!((vz - 2. - beta * vx * by).abs() < 1e-12);
    }
    #[test]
    fn magnetic_convergence_and_reference_flow() {
        let mut c = Config::default();
        let by = c.b_t * c.dip_deg.to_radians().cos();
        let bz = -c.b_t * c.dip_deg.to_radians().sin();
        c.electric_north_v_m = c.advection_m_s * bz;
        c.electric_up_v_m = -c.advection_m_s * by;
        let (vx, wz) = ion_velocity(c.advection_m_s, 0., &c);
        assert!((vx - c.advection_m_s).abs() < 1e-12);
        assert!(wz.abs() < 1e-12);
        // Negative dU/dz: ions rise below the zero crossing and sink above it.
        assert!(ion_velocity(c.advection_m_s + 20., 0., &c).1 > 0.);
        assert!(ion_velocity(c.advection_m_s - 20., 0., &c).1 < 0.);
        c.b_t = 0.;
        c.electric_north_v_m = 0.;
        c.electric_up_v_m = 0.;
        assert!(ion_velocity(c.advection_m_s - 20., 0., &c).1.abs() < 1e-12);
    }
    #[test]
    fn localized_packet_has_no_boundary_seed() {
        let mut c = Config::default();
        c.lx_m = 172000.;
        c.seed_packet_center_m = 86000.;
        c.seed_packet_half_width_m = 64500.;
        c.seed_packet_taper_m = 10750.;
        assert_eq!(packet_window(0., &c), 0.);
        assert_eq!(packet_window(c.lx_m, &c), 0.);
        assert_eq!(packet_window(c.seed_packet_center_m, &c), 1.);
        let midpoint =
            c.seed_packet_center_m + c.seed_packet_half_width_m + 0.5 * c.seed_packet_taper_m;
        assert!((packet_window(midpoint, &c) - 0.5).abs() < 1e-12);
    }
}
