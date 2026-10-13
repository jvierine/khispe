use anyhow::{Result, ensure};
use mpi::traits::*;
use serde::Serialize;
use std::{f64::consts::PI, path::Path};

const N2: f64 = 4.0e-4;
const PACK: usize = 24;

#[derive(Clone, Copy, Debug, Default, Serialize)]
struct Wave {
    lambda_x_m: f64,
    theta_rad: f64,
    lambda_z_m: f64,
    mean_u_m_s: f64,
    omega_intrinsic_s: f64,
    omega_observed_s: f64,
    kx_m: f64,
}

#[derive(Clone, Debug, Serialize)]
struct Fit {
    train_rmse_km: f64,
    heldout_rmse_km: f64,
    bic_train: f64,
    coefficients: Vec<f64>,
    waves: Vec<Wave>,
}

#[derive(Clone, Copy)]
struct Obs {
    x: f64,
    t: f64,
    z: f64,
    panel: usize,
}

fn splitmix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}
fn uniform(index: u64, lane: u64) -> f64 {
    let x = splitmix(index ^ lane.wrapping_mul(0xd2b74407b1ce6e93));
    ((x >> 11) as f64) * (1.0 / ((1u64 << 53) as f64))
}
fn sample_wave(index: u64, offset: u64) -> Wave {
    let lerp = |a: f64, b: f64, u: f64| a + (b - a) * u;
    let lambda_x_m = lerp(12e3, 50e3, uniform(index, offset));
    let theta_rad = lerp(0., 70f64.to_radians(), uniform(index, offset + 1));
    let lambda_z_m =
        (3e3f64.ln() + (30e3f64.ln() - 3e3f64.ln()) * uniform(index, offset + 2)).exp();
    let mean_u_m_s = lerp(-180., 40., uniform(index, offset + 3));
    let kx_m = 2. * PI / lambda_x_m;
    let kh = kx_m / theta_rad.cos();
    let m = 2. * PI / lambda_z_m;
    let omega_intrinsic_s = N2.sqrt() * kh / (kh * kh + m * m).sqrt();
    let omega_observed_s = omega_intrinsic_s + kx_m * mean_u_m_s;
    Wave {
        lambda_x_m,
        theta_rad,
        lambda_z_m,
        mean_u_m_s,
        omega_intrinsic_s,
        omega_observed_s,
        kx_m,
    }
}

fn basis(o: Obs, waves: &[Wave]) -> Vec<f64> {
    let mut b = vec![1., o.t / 600., o.x / 40_000.];
    for w in waves {
        let phase = w.kx_m * o.x - w.omega_observed_s * o.t;
        b.push(phase.cos());
        b.push(phase.sin());
    }
    b
}
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for k in 0..n {
        let p = (k..n).max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))?;
        if a[p][k].abs() < 1e-12 {
            return None;
        }
        a.swap(k, p);
        b.swap(k, p);
        let pivot = a[k][k];
        for j in k..n {
            a[k][j] /= pivot;
        }
        b[k] /= pivot;
        for i in 0..n {
            if i != k {
                let f = a[i][k];
                for j in k..n {
                    a[i][j] -= f * a[k][j];
                }
                b[i] -= f * b[k];
            }
        }
    }
    Some(b)
}
fn fit_linear(obs: &[Obs], waves: &[Wave]) -> Option<(Vec<f64>, f64)> {
    let p = 3 + 2 * waves.len();
    let mut ata = vec![vec![0.; p]; p];
    let mut aty = vec![0.; p];
    for &o in obs {
        let v = basis(o, waves);
        for i in 0..p {
            aty[i] += v[i] * o.z;
            for j in 0..p {
                ata[i][j] += v[i] * v[j];
            }
        }
    }
    let q = solve(ata, aty)?;
    let sse = obs
        .iter()
        .map(|&o| {
            let v = basis(o, waves);
            let e = v.iter().zip(&q).map(|(a, b)| a * b).sum::<f64>() - o.z;
            e * e
        })
        .sum();
    Some((q, sse))
}
fn rmse(obs: &[Obs], waves: &[Wave], q: &[f64]) -> f64 {
    (obs.iter()
        .map(|&o| {
            let v = basis(o, waves);
            let e = v.iter().zip(q).map(|(a, b)| a * b).sum::<f64>() - o.z;
            e * e
        })
        .sum::<f64>()
        / obs.len() as f64)
        .sqrt()
}
fn search(obs: &[Obs], nw: usize, count: u64, rank: i32, size: i32) -> Fit {
    let mut best_sse = f64::INFINITY;
    let mut best_q = Vec::new();
    let mut best_w = Vec::new();
    let mut i = rank as u64;
    while i < count {
        let waves: Vec<_> = (0..nw).map(|j| sample_wave(i, 1 + 4 * j as u64)).collect();
        if let Some((q, sse)) = fit_linear(obs, &waves) {
            if sse < best_sse {
                best_sse = sse;
                best_q = q;
                best_w = waves;
            }
        }
        i += size as u64;
    }
    Fit {
        train_rmse_km: (best_sse / obs.len() as f64).sqrt(),
        heldout_rmse_km: f64::NAN,
        bic_train: f64::NAN,
        coefficients: best_q,
        waves: best_w,
    }
}
fn pack(f: &Fit) -> [f64; PACK] {
    let mut p = [f64::NAN; PACK];
    p[0] = f.train_rmse_km;
    p[1] = f.waves.len() as f64;
    for (j, w) in f.waves.iter().enumerate() {
        let k = 2 + 7 * j;
        p[k] = w.lambda_x_m;
        p[k + 1] = w.theta_rad;
        p[k + 2] = w.lambda_z_m;
        p[k + 3] = w.mean_u_m_s;
        p[k + 4] = w.omega_intrinsic_s;
        p[k + 5] = w.omega_observed_s;
        p[k + 6] = w.kx_m;
    }
    let s = 16;
    for (i, &v) in f.coefficients.iter().enumerate() {
        p[s + i] = v;
    }
    p
}
fn unpack(p: &[f64; PACK]) -> Fit {
    let nw = p[1] as usize;
    let mut waves = Vec::new();
    for j in 0..nw {
        let k = 2 + 7 * j;
        waves.push(Wave {
            lambda_x_m: p[k],
            theta_rad: p[k + 1],
            lambda_z_m: p[k + 2],
            mean_u_m_s: p[k + 3],
            omega_intrinsic_s: p[k + 4],
            omega_observed_s: p[k + 5],
            kx_m: p[k + 6],
        });
    }
    Fit {
        train_rmse_km: p[0],
        heldout_rmse_km: f64::NAN,
        bic_train: f64::NAN,
        coefficients: p[16..16 + 3 + 2 * nw].to_vec(),
        waves,
    }
}

fn finalize(mut f: Fit, train: &[Obs], test: &[Obs]) -> Fit {
    f.train_rmse_km = rmse(train, &f.waves, &f.coefficients);
    f.heldout_rmse_km = rmse(test, &f.waves, &f.coefficients);
    let dof = 3 + 6 * f.waves.len();
    let sse = f.train_rmse_km.powi(2) * train.len() as f64;
    f.bic_train = train.len() as f64 * (sse / train.len() as f64).ln()
        + dof as f64 * (train.len() as f64).ln();
    f
}

fn model_fields(f: &Fit, x: &[f32], times: &[i32]) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = x.len() * times.len();
    let mut eta = Vec::with_capacity(n);
    let mut u = Vec::with_capacity(n);
    let mut wv = Vec::with_capacity(n);
    let mut buoy = Vec::with_capacity(n);
    let t0 = times[0] as f64;
    for &ti in times {
        for &xi in x {
            let o = Obs {
                x: xi as f64 * 1000.,
                t: ti as f64 - t0,
                z: 0.,
                panel: 0,
            };
            let v = basis(o, &f.waves);
            eta.push(v.iter().zip(&f.coefficients).map(|(a, b)| a * b).sum());
            let (mut uu, mut ww, mut bb) = (0., 0., 0.);
            for (j, wa) in f.waves.iter().enumerate() {
                let ph = wa.kx_m * o.x - wa.omega_observed_s * o.t;
                let c = f.coefficients[3 + 2 * j];
                let s = f.coefficients[4 + 2 * j];
                let e = c * ph.cos() + s * ph.sin();
                let wi = 1000. * wa.omega_intrinsic_s * (c * ph.sin() - s * ph.cos());
                let kh = wa.kx_m / wa.theta_rad.cos();
                let m = 2. * PI / wa.lambda_z_m;
                uu += -(m / kh) * wi * wa.theta_rad.cos();
                ww += wi;
                bb += -N2 * 1000. * e;
            }
            u.push(uu);
            wv.push(ww);
            buoy.push(bb);
        }
    }
    (eta, u, wv, buoy)
}

fn run(world: &mpi::topology::SimpleCommunicator) -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let input = args
        .get(1)
        .map(String::as_str)
        .unwrap_or("data/derived/figure3_digitized.h5");
    let output = args
        .get(2)
        .map(String::as_str)
        .unwrap_or("output/figure3-fit/fit.h5");
    let h = hdf5::File::open(input)?;
    let ds = h.dataset("figure3/log10_ne_m3")?;
    ensure!(ds.shape() == [15, 69, 140], "unexpected digitized shape");
    let d = ds.read_raw::<f32>()?;
    let x = h.dataset("figure3/distance_km")?.read_raw::<f32>()?;
    let z = h.dataset("figure3/altitude_km")?.read_raw::<f32>()?;
    let times = h.dataset("figure3/time_seconds_lt")?.read_raw::<i32>()?;
    let mut obs = Vec::new();
    for p in 0..15 {
        for ix in 0..140 {
            if x[ix] < -25. || x[ix] > 48. {
                continue;
            }
            let (mut sw, mut swz) = (0., 0.);
            for iz in 0..69 {
                let v = d[p * 69 * 140 + iz * 140 + ix] as f64;
                if v.is_finite() {
                    let wt = (v - 10.95).max(0.).powi(2);
                    sw += wt;
                    swz += wt * z[iz] as f64;
                }
            }
            if sw > 0. {
                obs.push(Obs {
                    x: x[ix] as f64 * 1000.,
                    t: (times[p] - times[0]) as f64,
                    z: swz / sw,
                    panel: p,
                });
            }
        }
    }
    let train: Vec<_> = obs.iter().copied().filter(|o| o.panel < 6).collect();
    let test: Vec<_> = obs.iter().copied().filter(|o| o.panel >= 6).collect();
    let base_q = fit_linear(&train, &[]).unwrap().0;
    let base = finalize(
        Fit {
            train_rmse_km: 0.,
            heldout_rmse_km: 0.,
            bic_train: 0.,
            coefficients: base_q,
            waves: vec![],
        },
        &train,
        &test,
    );
    let local1 = search(&train, 1, 100_000, world.rank(), world.size());
    let local2 = search(&train, 2, 100_000, world.rank(), world.size());
    let mut best1 = local1.clone();
    let mut best2 = local2.clone();
    if world.rank() == 0 {
        for r in 1..world.size() {
            let mut a = [0.; PACK];
            world.process_at_rank(r).receive_into(&mut a);
            let f = unpack(&a);
            if f.train_rmse_km < best1.train_rmse_km {
                best1 = f;
            }
            let mut b = [0.; PACK];
            world.process_at_rank(r).receive_into(&mut b);
            let f = unpack(&b);
            if f.train_rmse_km < best2.train_rmse_km {
                best2 = f;
            }
        }
    } else {
        world.process_at_rank(0).send(&pack(&local1));
        world.process_at_rank(0).send(&pack(&local2));
    }
    if world.rank() == 0 {
        best1 = finalize(best1, &train, &test);
        best2 = finalize(best2, &train, &test);
        let fits = [base, best1, best2];
        let out = Path::new(output);
        std::fs::create_dir_all(out.parent().unwrap())?;
        let file = hdf5::File::create(out)?;
        file.new_dataset_builder()
            .with_data(&x)
            .create("distance_km")?;
        file.new_dataset_builder()
            .with_data(&times)
            .create("time_seconds_lt")?;
        let names = ["background", "one_gw", "two_gw"];
        for (name, f) in names.iter().zip(&fits) {
            let g = file.create_group(name)?;
            g.new_attr::<f64>()
                .create("train_rmse_km")?
                .write_scalar(&f.train_rmse_km)?;
            g.new_attr::<f64>()
                .create("heldout_rmse_km")?
                .write_scalar(&f.heldout_rmse_km)?;
            g.new_attr::<f64>()
                .create("bic_train")?
                .write_scalar(&f.bic_train)?;
            g.new_dataset_builder()
                .with_data(&f.coefficients)
                .create("coefficients")?;
            let wave_json = serde_json::to_vec_pretty(&f.waves)?;
            g.new_dataset_builder()
                .with_data(&wave_json)
                .create("wave_parameters_json")?;
            let (eta, u, w, b) = model_fields(f, &x, &times);
            for (n, a) in [
                ("ridge_altitude_km", eta),
                ("u_prime_m_s", u),
                ("w_prime_m_s", w),
                ("buoyancy_m_s2", b),
            ] {
                g.new_dataset::<f64>()
                    .shape((15, 140))
                    .create(n)?
                    .write_raw(&a)?;
            }
        }
        let report = serde_json::json!({"status":"figure_derived_not_raw_data","article_doi":"10.1029/2025GL116932","training_panels":"a-f","heldout_panels":"g-o","mpi_ranks":world.size(),"fits":fits});
        std::fs::write(
            out.with_extension("json"),
            serde_json::to_string_pretty(&report)?,
        )?;
        println!("{}", serde_json::to_string_pretty(&report)?);
    }
    Ok(())
}
fn main() {
    let u = mpi::initialize().expect("MPI");
    let w = u.world();
    if let Err(e) = run(&w) {
        eprintln!("rank {}: {e:#}", w.rank());
        w.abort(1);
    }
}
