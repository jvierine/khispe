use anyhow::{Result, ensure};
use font8x8::{BASIC_FONTS, UnicodeFonts};
use std::{fs, path::Path};
const WIDTH: usize = 960;
const HEIGHT: usize = 440;
fn text(p: &mut [u8], x: usize, y: usize, s: &str) {
    for (i, c) in s.chars().enumerate() {
        if let Some(g) = BASIC_FONTS.get(c) {
            for (dy, row) in g.iter().enumerate() {
                for dx in 0..8 {
                    if row & (1 << dx) != 0 {
                        for sy in 0..2 {
                            for sx in 0..2 {
                                let xx = x + i * 16 + dx * 2 + sx;
                                let yy = y + dy * 2 + sy;
                                if xx < WIDTH && yy < HEIGHT {
                                    p[yy * WIDTH + xx] = 255;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
fn palette() -> Vec<u8> {
    let mut p = vec![0; 256 * 3];
    for k in 0..128 {
        let s = k as f64 / 127.;
        let colors = if s < 0.5 {
            let t = s * 2.;
            [30. + 225. * t, 80. + 175. * t, 220. + 35. * t]
        } else {
            let t = (s - 0.5) * 2.;
            [255., 255. - 215. * t, 255. - 220. * t]
        };
        for j in 0..3 {
            p[k * 3 + j] = colors[j] as u8;
        }
    }
    for k in 128..254 {
        let t = (k - 128) as f64 / 125.;
        let anchors = [
            [16., 18., 60.],
            [110., 35., 150.],
            [225., 70., 90.],
            [255., 230., 70.],
        ];
        let i = (t * 3.).floor().min(2.) as usize;
        let a = t * 3. - i as f64;
        for j in 0..3 {
            p[k * 3 + j] = (anchors[i][j] * (1. - a) + anchors[i + 1][j] * a) as u8;
        }
    }
    p[254 * 3..255 * 3].copy_from_slice(&[8, 20, 30]);
    p[255 * 3..].copy_from_slice(&[230, 240, 245]);
    p
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let dir = Path::new(args.get(1).map(String::as_str).unwrap_or("output/pilot"));
    let mut paths: Vec<_> = fs::read_dir(dir)?
        .filter_map(|x| x.ok().map(|x| x.path()))
        .filter(|p| p.extension().and_then(|x| x.to_str()) == Some("h5"))
        .collect();
    paths.sort();
    ensure!(!paths.is_empty(), "no frames");
    let mut frames = Vec::new();
    let mut maxwind: f64 = 0.;
    let mut nz = 0;
    let mut nx = 0;
    let mut altitude = Vec::new();
    let mut x = Vec::new();
    for path in &paths {
        let f = hdf5::File::open(path)?;
        let shape = f.dataset("w_m_s")?.shape();
        nx = shape[0];
        nz = shape[1];
        let w = f.dataset("w_m_s")?.read_raw::<f64>()?;
        let n = f.dataset("ion_density_m-3")?.read_raw::<f64>()?;
        maxwind = w.iter().fold(maxwind, |a, b| a.max(b.abs()));
        let t = f.attr("time_s")?.read_scalar::<f64>()?;
        altitude = f.dataset("altitude_m")?.read_raw::<f64>()?;
        x = f.dataset("x_m")?.read_raw::<f64>()?;
        frames.push((t, w, n));
    }
    maxwind = maxwind.max(0.01);
    let output = dir.join("neutral-ion-pilot.gif");
    let mut encoder = gif::Encoder::new(
        fs::File::create(&output)?,
        WIDTH as u16,
        HEIGHT as u16,
        &palette(),
    )?;
    encoder.set_repeat(gif::Repeat::Infinite)?;
    for (t, w, n) in &frames {
        let mut p = vec![254; WIDTH * HEIGHT];
        text(
            &mut p,
            24,
            14,
            &format!("SYISR KHI PILOT | t={t:5.0} s | 2D Rust/MPI"),
        );
        text(&mut p, 56, 52, "NEUTRAL VERTICAL WIND");
        text(&mut p, 540, 52, "ION DENSITY: E=0 PILOT");
        for panel in 0..2 {
            let left = if panel == 0 { 56 } else { 540 };
            for py in 0..256 {
                for px in 0..380 {
                    let i = px * nx / 380;
                    let j = nz - 1 - py * nz / 256;
                    let k = i * nz + j;
                    let color = if panel == 0 {
                        ((w[k] / maxwind * 0.5 + 0.5).clamp(0., 1.) * 127.) as u8
                    } else {
                        128 + (((n[k].max(1e7).log10() - 7.) / 4.).clamp(0., 1.) * 125.) as u8
                    };
                    p[(80 + py) * WIDTH + left + px] = color;
                }
            }
            text(
                &mut p,
                left,
                342,
                &format!("x: {:.1}-{:.1} km", x[0] / 1000., x[nx - 1] / 1000.),
            );
            for px in 0..380 {
                let v = if panel == 0 {
                    (px as f64 / 379. * 127.) as u8
                } else {
                    128 + (px as f64 / 379. * 125.) as u8
                };
                for py in 0..8 {
                    p[(368 + py) * WIDTH + left + px] = v;
                }
            }
        }
        text(
            &mut p,
            56,
            384,
            &format!("-{maxwind:.1} to +{maxwind:.1} m/s"),
        );
        text(&mut p, 540, 384, "1e7 to 1e11 ions/m^3");
        text(
            &mut p,
            24,
            412,
            &format!(
                "z={:.1}-{:.1} km | NOT A FIT TO SYISR DATA",
                altitude[0] / 1000.,
                altitude[nz - 1] / 1000.
            ),
        );
        let frame = gif::Frame {
            width: WIDTH as u16,
            height: HEIGHT as u16,
            buffer: p.into(),
            delay: 10,
            ..Default::default()
        };
        encoder.write_frame(&frame)?;
    }
    let rms = |w: &[f64]| (w.iter().map(|x| x * x).sum::<f64>() / w.len() as f64).sqrt();
    let first = rms(&frames[0].1);
    let last = rms(&frames.last().unwrap().1);
    let report = serde_json::json!({"initial_vertical_wind_rms_m_s":first,"final_vertical_wind_rms_m_s":last,"rms_growth_factor":last/first,"fixed_colour_limit_m_s":maxwind,"frames":frames.len(),"status":"exploratory_2D_pilot_not_observation_fit","renderer":"src/bin/render-pilot.rs"});
    fs::write(
        dir.join("pilot-report.json"),
        serde_json::to_string_pretty(&report)?,
    )?;
    println!("{report}\n{}", output.display());
    Ok(())
}
