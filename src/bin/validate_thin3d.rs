use anyhow::{ensure, Result};

fn main() -> Result<()> {
    let path = std::env::args()
        .nth(1)
        .expect("usage: validate_thin3d FRAME.h5");
    let file = hdf5::File::open(path)?;

    for name in [
        "u_m_s",
        "v_m_s",
        "w_m_s",
        "buoyancy_m_s2",
        "ion_density_m-3",
        "omega_y_s-1",
    ] {
        let dataset = file.dataset(name)?;
        let shape = dataset.shape();
        ensure!(shape.len() == 3, "{name} is not three-dimensional");
        let (nx, ny, nz) = (shape[0], shape[1], shape[2]);
        let values = dataset.read_raw::<f64>()?;
        let mut maximum_spread = 0.0_f64;
        let mut maximum_magnitude = 0.0_f64;
        for i in 0..nx {
            for k in 0..nz {
                let first = values[(i * ny) * nz + k];
                maximum_magnitude = maximum_magnitude.max(first.abs());
                for j in 1..ny {
                    let value = values[(i * ny + j) * nz + k];
                    maximum_magnitude = maximum_magnitude.max(value.abs());
                    maximum_spread = maximum_spread.max((value - first).abs());
                }
            }
        }
        let relative = maximum_spread / maximum_magnitude.max(f64::MIN_POSITIVE);
        println!(
            "{name}: max|q|={maximum_magnitude:.9e}, max y spread={maximum_spread:.9e}, relative={relative:.3e}"
        );
    }

    let u = file.dataset("u_m_s")?.read_raw::<f64>()?;
    let v = file.dataset("v_m_s")?.read_raw::<f64>()?;
    let w = file.dataset("w_m_s")?.read_raw::<f64>()?;
    let shape = file.dataset("u_m_s")?.shape();
    let (nx, ny, nz) = (shape[0], shape[1], shape[2]);
    let x = file.dataset("x_m")?.read_raw::<f64>()?;
    let y = file.dataset("y_m")?.read_raw::<f64>()?;
    let z = file.dataset("altitude_m")?.read_raw::<f64>()?;
    let (dx, dy, dz) = (x[1] - x[0], y[1] - y[0], z[1] - z[0]);
    let index = |i: usize, j: usize, k: usize| (i * ny + j) * nz + k;
    let mut divergence_squared = 0.0;
    let mut lower_wall_speed = 0.0_f64;
    for i in 0..nx {
        let im = (i + nx - 1) % nx;
        for j in 0..ny {
            let jm = (j + ny - 1) % ny;
            lower_wall_speed = lower_wall_speed.max(w[index(i, j, 0)].abs());
            for k in 0..nz {
                let q = index(i, j, k);
                let wu = if k + 1 == nz {
                    0.0
                } else {
                    w[index(i, j, k + 1)]
                };
                let divergence = (u[q] - u[index(im, j, k)]) / dx
                    + (v[q] - v[index(i, jm, k)]) / dy
                    + (wu - w[q]) / dz;
                divergence_squared += divergence * divergence;
            }
        }
    }
    println!(
        "projection: rms divergence={:.3e} s^-1, lower-wall max|w|={lower_wall_speed:.3e} m/s",
        (divergence_squared / (nx * ny * nz) as f64).sqrt()
    );
    Ok(())
}
