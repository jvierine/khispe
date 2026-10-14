use anyhow::{ensure, Result};
use rustplotlib::*;

fn edges_km(centers: &[f64]) -> Vec<f64> {
    let spacing = centers[1] - centers[0];
    (0..=centers.len())
        .map(|i| (centers[0] + (i as f64 - 0.5) * spacing) / 1_000.0)
        .collect()
}

fn density_cut(file: &hdf5::File) -> Result<Vec<f64>> {
    let dataset = file.dataset("ion_density_m-3")?;
    let shape = dataset.shape();
    let density = dataset.read_raw::<f64>()?;
    match shape.as_slice() {
        [nx, nz] => {
            let mut cut = Vec::with_capacity(nx * nz);
            for k in 0..*nz {
                for i in 0..*nx {
                    cut.push(density[i * nz + k]);
                }
            }
            Ok(cut)
        }
        [nx, ny, nz] => {
            let mut cut = Vec::with_capacity(nx * nz);
            let j = ny / 2;
            for k in 0..*nz {
                for i in 0..*nx {
                    cut.push(density[(i * ny + j) * nz + k]);
                }
            }
            Ok(cut)
        }
        _ => anyhow::bail!("ion density must be two- or three-dimensional"),
    }
}

fn main() -> Result<()> {
    let arguments: Vec<_> = std::env::args().collect();
    ensure!(
        arguments.len() == 4,
        "usage: thin_validation REFERENCE_2D.h5 THIN_3D.h5 OUTPUT.png"
    );
    let files = [
        hdf5::File::open(&arguments[1])?,
        hdf5::File::open(&arguments[2])?,
    ];
    let times = [
        files[0].attr("time_s")?.read_scalar::<f64>()?,
        files[1].attr("time_s")?.read_scalar::<f64>()?,
    ];
    ensure!((times[0] - times[1]).abs() < 1e-9, "frame times differ");

    let mut figure = Figure::new(IntendedUse::ArticleWide);
    figure.height_mm = 100.0;
    figure.raster_dpi = 300.0;
    let scale = figure.scale(Scale::new(
        [1e10, 2e12],
        Normalization::Log,
        ColorMap::ember(),
    )?);
    for (column, (file, title)) in files
        .iter()
        .zip(["a) 2-D reference", "b) Thin 3-D central cut"])
        .enumerate()
    {
        let x = file.dataset("x_m")?.read_raw::<f64>()?;
        let z = file.dataset("altitude_m")?.read_raw::<f64>()?;
        let x_edges = edges_km(&x);
        let z_edges = edges_km(&z);
        let left = 17.0 + column as f64 * 87.0;
        let axes = Axes::new(
            [x_edges[0], *x_edges.last().unwrap()],
            [115.0, 121.0],
            "Zonal distance (km)",
            if column == 0 { "Altitude (km)" } else { "" },
        );
        let mut panel = Panel::new(Rect::new(left, 13.0, 79.0, 48.0), axes);
        let mut field = Field::new(x_edges, z_edges, density_cut(file)?, scale)?;
        field.interpolate = true;
        panel.push(Layer::Pcolormesh(field));
        figure.panels.push(panel);
        figure.labels.push(Label::new(title, left + 39.5, 8.0));
    }
    let mut colorbar = Colorbar::new(
        Rect::new(34.0, 80.0, 126.0, 2.0),
        scale,
        r"$n_e\;[\mathrm{m}^{-3}]$",
        vec![1e10, 1e11, 2e12],
    );
    colorbar.horizontal = true;
    figure.colorbars.push(colorbar);
    figure.labels.push(Label::new(
        format!(r"$t={:.1}\,\mathrm{{min}}$", times[0] / 60.0),
        97.0,
        3.0,
    ));
    figure.save(&arguments[3])?;
    Ok(())
}
