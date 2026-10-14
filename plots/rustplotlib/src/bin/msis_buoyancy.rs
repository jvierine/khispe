use anyhow::Result;
use rustplotlib::*;
fn main() -> Result<()> {
    let base = std::env::args()
        .nth(1)
        .unwrap_or("../syisr-simulation-memos-overleaf".into());
    let h = hdf5::File::open(format!("{base}/data/msis_buoyancy.h5"))?;
    let z = h.dataset("altitude_km")?.read_raw::<f64>()?;
    let n = h.dataset("N_s-1")?.read_raw::<f64>()?;
    let mut f = Figure::new(IntendedUse::ArticleWide);
    f.height_mm = 105.;
    let mut ax = Axes::new(
        [0.01, 0.035],
        [70., 150.],
        r"$N\,[\mathrm{s}^{-1}]$",
        r"Height $z\,[\mathrm{km}]$",
    );
    ax.xticks = vec![0.01, 0.02, 0.03];
    ax.yticks = vec![70., 90., 110., 130., 150.];
    ax.secondary_x = Some(SecondaryAxis::new(
        "Buoyancy period [min]",
        vec![3., 4., 5., 7., 10.],
        |period| std::f64::consts::TAU / (60. * period),
    )?);
    let mut p = Panel::new(Rect::new(18., 23., 88., 68.), ax);
    let mut l = Line::new(vec![[0.02, 70.], [0.02, 150.]], "#777777");
    l.dash = Some("4 3".into());
    p.push(Layer::Line(l));
    let mut l = Line::new(
        n.iter().zip(z.iter()).map(|(&n, &z)| [n, z]).collect(),
        "#245b78",
    );
    l.width_pt = 1.3;
    p.push(Layer::Line(l));
    f.panels.push(p);
    let hv = hdf5::File::open(format!("{base}/data/viscosity_reference.h5"))?;
    let zv = hv.dataset("altitude_km")?.read_raw::<f64>()?;
    let nu = hv.dataset("kinematic_viscosity_m2_s")?.read_raw::<f64>()?;
    let mut ax = Axes::new(
        [-1., 4.5],
        [70., 150.],
        r"$\log_{10}(\nu/[\mathrm{m}^2\,\mathrm{s}^{-1}])$",
        "",
    );
    ax.xticks = vec![0., 1., 2., 3., 4.];
    ax.yticks = vec![];
    let mut p = Panel::new(Rect::new(118., 23., 88. / 1.618033988749895, 68.), ax);
    let mut fixed = Line::new(
        vec![[30_f64.log10(), 70.], [30_f64.log10(), 150.]],
        "#777777",
    );
    fixed.dash = Some("4 3".into());
    p.push(Layer::Line(fixed));
    let mut line = Line::new(
        nu.iter()
            .zip(zv.iter())
            .map(|(&v, &z)| [v.log10(), z])
            .collect(),
        "#a24b3d",
    );
    line.width_pt = 1.3;
    p.push(Layer::Line(line));
    f.panels.push(p);
    f.labels.push(Label::new("Kinematic viscosity", 145., 13.));
    f.save(format!("{base}/figures/msis_buoyancy.pdf"))?;
    f.save(format!("{base}/figures/msis_buoyancy.png"))?;
    Ok(())
}
