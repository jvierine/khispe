use anyhow::Result;
use rustplotlib::*;
use std::f64::consts::PI;
fn arrow(p: &mut Panel, a: [f64; 2], b: [f64; 2], color: &str) {
    let mut l = Line::new(vec![a, b], color);
    l.width_pt = 1.;
    p.push(Layer::Line(l));
    let dx = b[0] - a[0];
    let dz = b[1] - a[1];
    let n = (dx * dx + dz * dz * 100.).sqrt();
    let ux = dx / n;
    let uz = dz / n;
    for sign in [-1., 1.] {
        let end = [
            b[0] - 2. * ux + sign * 0.8 * uz * 10.,
            b[1] - 2. * uz - sign * 0.08 * ux,
        ];
        let mut l = Line::new(vec![b, end], color);
        l.width_pt = 1.;
        p.push(Layer::Line(l));
    }
}
fn main() -> Result<()> {
    let h = hdf5::File::open("data/ri008-memo/plot_data.h5")?;
    let out = "../syisr-khi-memo-overleaf/figures";
    let x = h.dataset("x_edges_km")?.read_raw::<f64>()?;
    let z = h.dataset("z_edges_km")?.read_raw::<f64>()?;
    let mut f = Figure::new(IntendedUse::ArticleWide);
    f.height_mm = 181.;
    let ne = f.scale(Scale::new(
        [1e10, 2e12],
        Normalization::Log,
        ColorMap::ember(),
    )?);
    let vort = f.scale(Scale::new(
        [-0.08, 0.08],
        Normalization::Linear,
        ColorMap::diverging(),
    )?);
    f.labels.push(Label::new("Electron density", 56., 4.));
    f.labels.push(Label::new("Neutral vorticity", 142., 4.));
    for (row, key) in ["120", "240", "480"].iter().enumerate() {
        let g = h.group(key)?;
        let t = g.attr("time_s")?.read_scalar::<f64>()?;
        for (col, name, sc) in [(0, "density", ne), (1, "vorticity", vort)] {
            let xx = 19. + col as f64 * 86.;
            let yy = 12. + row as f64 * 49.;
            let mut ax = Axes::new(
                [-45., 5.],
                [-3., 3.],
                if row == 2 {
                    r"$x-U_ct\,[\mathrm{km}]$"
                } else {
                    ""
                },
                if col == 0 {
                    r"$\zeta\,[\mathrm{km}]$"
                } else {
                    ""
                },
            );
            ax.xticks = vec![-40., -20., 0.];
            ax.yticks = vec![-2., 0., 2.];
            let mut p = Panel::new(Rect::new(xx, yy, 70., 33.), ax);
            let mut field = Field::new(
                x.clone(),
                z.clone(),
                g.dataset(name)?.read_raw::<f64>()?,
                sc,
            )?;
            field.interpolate = true;
            p.push(Layer::Pcolormesh(field));
            if row < 2 {
                let targets = if row == 0 {
                    [-36.79, -15.96]
                } else {
                    [-35.78, -16.63]
                };
                arrow(
                    &mut p,
                    [targets[0] - 8., 2.2],
                    [targets[0], 0.],
                    if col == 0 { "#ffffff" } else { "#151515" },
                );
                arrow(
                    &mut p,
                    [targets[1] + 8., 2.2],
                    [targets[1], 0.],
                    if col == 0 { "#ffffff" } else { "#151515" },
                );
            } else {
                arrow(
                    &mut p,
                    [-39., 2.2],
                    [-29.73, 0.],
                    if col == 0 { "#ffffff" } else { "#151515" },
                );
            }
            f.labels.push(Label::new(
                format!("$t={:.0}\\,\\mathrm{{min}}$", t / 60.),
                xx + 35.,
                yy - 4.,
            ));
            f.panels.push(p);
        }
    }
    for (xx, sc, label, ticks) in [
        (19., ne, r"$n_e\,[\mathrm{m}^{-3}]$", vec![1e10, 1e11, 2e12]),
        (
            105.,
            vort,
            r"$\Omega\,[\mathrm{s}^{-1}]$",
            vec![-0.08, 0., 0.08],
        ),
    ] {
        let mut bar = Colorbar::new(Rect::new(xx, 160., 70., 2.), sc, label, ticks);
        bar.horizontal = true;
        f.colorbars.push(bar);
    }
    f.save(format!("{out}/ri008_pairing.pdf"))?;
    f.save(format!("{out}/ri008_pairing.png"))?;
    let mut f = Figure::new(IntendedUse::ArticleWide);
    f.height_mm = 83.;
    let mut ax = Axes::new(
        [-86., 86.],
        [-4.5, 4.5],
        r"$x\,[\mathrm{km}]$",
        r"$w'(x,0,0)\,[\mathrm{m\,s}^{-1}]$",
    );
    ax.xticks = vec![-80., -40., 0., 40., 80.];
    ax.yticks = vec![-4., 0., 4.];
    let mut p = Panel::new(Rect::new(19., 12., 88., 52.), ax);
    for (mode, color) in [(0, "#245774"), (1, "#b64442"), (2, "#181818")] {
        p.push(Layer::Line(Line::new(
            (0..1001)
                .map(|j| {
                    let x = -86. + 172. * j as f64 / 1000.;
                    let primary = -2.5 * (2. * PI * x / 21.5).cos();
                    let pairing = -1.5 * (2. * PI * x / 43. + 0.7).cos();
                    [
                        x,
                        match mode {
                            0 => primary,
                            1 => pairing,
                            _ => primary + pairing,
                        },
                    ]
                })
                .collect(),
            color,
        )));
    }
    f.panels.push(p);
    f.labels
        .push(Label::new("Primary, subharmonic, sum", 63., 5.));
    let mut ax = Axes::new(
        [-180., -20.],
        [-3., 3.],
        r"$U(\zeta)\,[\mathrm{m\,s}^{-1}]$",
        r"$\zeta\,[\mathrm{km}]$",
    );
    ax.xticks = vec![-175., -100., -25.];
    ax.yticks = vec![-2., 0., 2.];
    let mut p = Panel::new(Rect::new(133., 12., 40., 52.), ax);
    p.push(Layer::Line(Line::new(
        (0..401)
            .map(|j| {
                let z = -3. + 6. * j as f64 / 400.;
                [-100. - 74.95331880577403 * (z / 1.06).tanh(), z]
            })
            .collect(),
        "#181818",
    )));
    f.panels.push(p);
    f.labels.push(Label::new("Background shear", 153., 5.));
    f.save(format!("{out}/ri008_initial.pdf"))?;
    f.save(format!("{out}/ri008_initial.png"))?;
    Ok(())
}
