use anyhow::Result;
use rustplotlib::*;
use serde_json::Value;
use std::{fs, path::PathBuf};
fn num(c: &Value, k: &str) -> f64 {
    c[k].as_f64().unwrap()
}
fn ion_up(u: f64, w: f64, c: &Value) -> f64 {
    let mu = 1.602176634e-19 / (num(c, "ion_mass_amu") * 1.66053906660e-27 * num(c, "nu_in_s"));
    let beta = mu * num(c, "b_t");
    let dip = num(c, "dip_deg").to_radians();
    let by = dip.cos();
    let bz = -dip.sin();
    let fx = u + mu * num(c, "electric_east_v_m");
    let fy = mu * num(c, "electric_north_v_m");
    let fz = w + mu * num(c, "electric_up_v_m");
    (fz + beta * fx * by + beta * beta * bz * (fy * by + fz * bz)) / (1. + beta * beta)
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let h = hdf5::File::open(&args[1])?;
    let out = PathBuf::from(&args[2]);
    fs::create_dir_all(out.join("ne"))?;
    fs::create_dir_all(out.join("ion_up"))?;
    fs::create_dir_all(out.join("vorticity"))?;
    let mut keys = h.member_names()?;
    keys.sort();
    let g = h.group(&keys[0])?;
    let c: Value = serde_json::from_slice(&g.dataset("config_json")?.read_raw::<u8>()?)?;
    let shape = g.dataset("u_m_s")?.shape();
    let (nx, nz) = (shape[0], shape[1]);
    let lx = num(&c, "lx_m");
    let dx = lx / nx as f64;
    let z = g.dataset("altitude_m")?.read_raw::<f64>()?;
    let dz = z[1] - z[0];
    let rel: Vec<_> = z
        .iter()
        .map(|x| (x - num(&c, "altitude_m")) / 1000.)
        .collect();
    let xedges: Vec<_> = (0..=nx)
        .map(|i| (i as f64 * dx - lx / 2.) / 1000.)
        .collect();
    let zedges: Vec<_> = (0..=nz)
        .map(|j| rel[0] + (j as f64 - 0.5) * dz / 1000.)
        .collect();
    let mut peak: f64 = 0.;
    let mut vort_peak: f64 = 0.;
    for key in &keys {
        let g = h.group(key)?;
        let u = g.dataset("u_m_s")?.read_raw::<f32>()?;
        let w = g.dataset("w_m_s")?.read_raw::<f32>()?;
        let omega = curl(&u, &w, nx, nz, dx, dz);
        vort_peak = omega.iter().fold(vort_peak, |m, v| m.max(v.abs()));
        for (i, &uu) in u.iter().enumerate() {
            peak = peak.max(ion_up(uu as f64, w[i] as f64, &c).abs());
        }
    }
    let limit = (peak / 5.).ceil().max(1.) * 5.;
    let vort_limit = (vort_peak / 0.01).ceil().max(1.) * 0.01;
    println!("Vorticity fixed limits +/- {vort_limit} s^-1");
    println!("Ion vertical-drift fixed limits +/- {limit} m/s; max {peak}");
    let mut tex = Latex::default();
    tex.font_pt = 10.;
    tex.typeset(keys.iter().map(|key| {
        let t = h
            .group(key)
            .unwrap()
            .attr("time_s")
            .unwrap()
            .read_scalar::<f64>()
            .unwrap();
        format!("$t={:.1}\\,\\mathrm{{min}}$", t / 60.)
    }))?;
    for (k, key) in keys.iter().enumerate() {
        let g = h.group(key)?;
        let t = g.attr("time_s")?.read_scalar::<f64>()?;
        let q = g.dataset("log10_ne_quantized")?.read_raw::<u16>()?;
        let u = g.dataset("u_m_s")?.read_raw::<f32>()?;
        let w = g.dataset("w_m_s")?.read_raw::<f32>()?;
        let omega = curl(&u, &w, nx, nz, dx, dz);
        let mut density = vec![0.; nx * nz];
        let mut vz = density.clone();
        let mut vort = density.clone();
        let shift = num(&c, "advection_m_s") * t / dx;
        for j in 0..nz {
            for i in 0..nx {
                let p = i as f64 + shift;
                let a = p.floor() as i64;
                let f = p - p.floor();
                let aa = a.rem_euclid(nx as i64) as usize * nz + j;
                let bb = (a + 1).rem_euclid(nx as i64) as usize * nz + j;
                vort[j * nx + i] = (1. - f) * omega[aa] + f * omega[bb];
                density[j * nx + i] =
                    10f64.powf(((1. - f) * q[aa] as f64 + f * q[bb] as f64) / 4096.);
                vz[j * nx + i] = (1. - f) * ion_up(u[aa] as f64, w[aa] as f64, &c)
                    + f * ion_up(u[bb] as f64, w[bb] as f64, &c);
            }
        }
        for (mode, values) in [("ne", density), ("ion_up", vz), ("vorticity", vort)] {
            if args.get(3).is_some_and(|selected| selected != mode) {
                continue;
            }
            let mut fig = Figure::new(IntendedUse::ArticleWide);
            fig.height_mm = 105.;
            fig.raster_dpi = 180.;
            fig.tight_export = false;
            let (left_width, sc, label, ticks) = if mode == "ne" {
                (
                    151.,
                    fig.scale(Scale::new(
                        [1e10, 2e12],
                        Normalization::Log,
                        ColorMap::ember(),
                    )?),
                    r"$n_e\,[\mathrm{m}^{-3}]$",
                    vec![1e10, 1e11, 2e12],
                )
            } else if mode == "vorticity" {
                (
                    151.,
                    fig.scale(Scale::new(
                        [-vort_limit, vort_limit],
                        Normalization::Linear,
                        ColorMap::diverging(),
                    )?),
                    r"$\Omega\,[\mathrm{s}^{-1}]$",
                    vec![-vort_limit, 0., vort_limit],
                )
            } else {
                (
                    116.,
                    fig.scale(Scale::new(
                        [-limit, limit],
                        Normalization::Linear,
                        ColorMap::diverging(),
                    )?),
                    r"$v_{i,z}\,[\mathrm{m\,s}^{-1}]$",
                    vec![-limit, 0., limit],
                )
            };
            let mut ax = Axes::new(
                [-lx / 2000., lx / 2000.],
                [zedges[0], *zedges.last().unwrap()],
                r"$x-U_0t\,[\mathrm{km}]$",
                r"$z-z_c(t)\,[\mathrm{km}]$",
            );
            ax.xticks = vec![-80., -40., 0., 40., 80.];
            ax.yticks = vec![-2., -1., 0., 1., 2.];
            let mut p = Panel::new(Rect::new(19., 10., left_width, 67.), ax);
            let mut field = Field::new(xedges.clone(), zedges.clone(), values, sc)?;
            field.interpolate = true;
            p.push(Layer::Pcolormesh(field));
            fig.panels.push(p);
            let mut bar =
                Colorbar::new(Rect::new(25., 89., left_width - 12., 2.), sc, label, ticks);
            bar.horizontal = true;
            fig.colorbars.push(bar);
            fig.labels.push(Label::new(
                format!("$t={:.1}\\,\\mathrm{{min}}$", t / 60.),
                90.,
                4.,
            ));
            if mode == "ion_up" {
                let mut ax = Axes::new(
                    [-180., -20.],
                    [zedges[0], *zedges.last().unwrap()],
                    r"$U\,[\mathrm{m\,s}^{-1}]$",
                    "",
                );
                ax.xticks = vec![-175., -100., -25.];
                ax.yticks = vec![];
                let mut panel = Panel::new(Rect::new(145., 10., 29., 67.), ax);
                panel.push(Layer::Line(Line::new(
                    rel.iter()
                        .map(|&zz| {
                            [
                                num(&c, "advection_m_s")
                                    + num(&c, "wind_m_s")
                                        * (zz * 1000. / num(&c, "shear_m")).tanh(),
                                zz,
                            ]
                        })
                        .collect(),
                    "#252525",
                )));
                fig.panels.push(panel);
                fig.labels.push(Label::new("Background wind", 159., 4.));
            }
            fig.save(out.join(mode).join(format!("{k:04}.png")))?;
            if [0, 360, 720].contains(&k) {
                fig.save(out.join(format!("{mode}_{k:04}.pdf")))?;
            }
        }
        if k % 60 == 0 {
            println!("Rendered {k}/{}", keys.len());
        }
    }
    fs::write(
        out.join("provenance.json"),
        serde_json::to_vec_pretty(
            &serde_json::json!({"renderer":"rustplotlib","script":"plots/rustplotlib/src/main.rs","source_hdf5":args[1],"frames":keys.len(),"ion_velocity":"collisional Lorentz-neutral-drag solution, model frame, positive upward; no artificial phase-descent velocity added","ion_velocity_limits_m_s":[-limit,limit],"vorticity":"d_x w - d_z u; centered finite differences of saved float32 velocities, periodic x, one-sided vertical endpoints","vorticity_limits_s-1":[-vort_limit,vort_limit],"camera":"mean zonal drift comoving","mask":false,"config":c}),
        )?,
    )?;
    Ok(())
}

fn curl(u: &[f32], w: &[f32], nx: usize, nz: usize, dx: f64, dz: f64) -> Vec<f64> {
    let mut out = vec![0.; nx * nz];
    for i in 0..nx {
        for j in 0..nz {
            let im = (i + nx - 1) % nx;
            let ip = (i + 1) % nx;
            let jm = j.saturating_sub(1);
            let jp = (j + 1).min(nz - 1);
            out[i * nz + j] = (w[ip * nz + j] as f64 - w[im * nz + j] as f64) / (2. * dx)
                - (u[i * nz + jp] as f64 - u[i * nz + jm] as f64) / ((jp - jm) as f64 * dz);
        }
    }
    out
}
