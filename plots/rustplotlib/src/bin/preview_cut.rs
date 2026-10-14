use anyhow::{ensure, Result};
use rustplotlib::*;

fn cell_edges(centers: &[f64]) -> Vec<f64> {
    let spacing = centers[1] - centers[0];
    (0..=centers.len())
        .map(|i| (centers[0] + (i as f64 - 0.5) * spacing) / 1_000.0)
        .collect()
}

fn main() -> Result<()> {
    let arguments: Vec<_> = std::env::args().collect();
    ensure!(
        arguments.len() == 3,
        "usage: preview_cut FRAME.h5 OUTPUT.png"
    );
    let input = hdf5::File::open(&arguments[1])?;
    let x = input.dataset("x_m")?.read_raw::<f64>()?;
    let z = input.dataset("altitude_m")?.read_raw::<f64>()?;
    let dataset = input.dataset("ion_density_m-3")?;
    let shape = dataset.shape();
    let density = dataset.read_raw::<f64>()?;
    let mut values = Vec::with_capacity(x.len() * z.len());
    match shape.as_slice() {
        [nx, nz] => {
            ensure!((*nx, *nz) == (x.len(), z.len()), "inconsistent 2-D shape");
            for k in 0..z.len() {
                for i in 0..x.len() {
                    values.push(density[i * z.len() + k]);
                }
            }
        }
        [nx, ny, nz] => {
            ensure!((*nx, *nz) == (x.len(), z.len()), "inconsistent 3-D shape");
            let j = ny / 2;
            for k in 0..z.len() {
                for i in 0..x.len() {
                    values.push(density[(i * ny + j) * z.len() + k]);
                }
            }
        }
        _ => anyhow::bail!("ion density must be a 2-D field or a 3-D volume"),
    }

    let x_edges = cell_edges(&x);
    let z_edges = cell_edges(&z);
    let mut figure = Figure::new(IntendedUse::ArticleWide);
    figure.height_mm = 105.0;
    figure.raster_dpi = 300.0;
    let scale = figure.scale(Scale::new(
        [1e10, 2e12],
        Normalization::Log,
        ColorMap::ember(),
    )?);
    let axes = Axes::new(
        [x_edges[0], *x_edges.last().unwrap()],
        [z_edges[0], *z_edges.last().unwrap()],
        "Zonal distance (km)",
        "Altitude (km)",
    );
    let mut panel = Panel::new(Rect::new(20.0, 12.0, 145.0, 62.0), axes);
    let mut field = Field::new(x_edges, z_edges, values, scale)?;
    field.interpolate = true;
    panel.push(Layer::Pcolormesh(field));
    figure.panels.push(panel);
    let mut colorbar = Colorbar::new(
        Rect::new(25.0, 88.0, 135.0, 2.0),
        scale,
        r"$n_e\;[\mathrm{m}^{-3}]$",
        vec![1e10, 1e11, 2e12],
    );
    colorbar.horizontal = true;
    figure.colorbars.push(colorbar);
    let time = input.attr("time_s")?.read_scalar::<f64>()?;
    figure.labels.push(Label::new(
        format!(r"$t={:.1}\,\mathrm{{min}}$", time / 60.0),
        92.0,
        5.0,
    ));
    figure.save(&arguments[2])?;
    Ok(())
}
