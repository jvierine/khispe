use anyhow::Result;
use rustplotlib::*;
fn main()->Result<()> {
 let base=std::env::args().nth(1).unwrap_or("../syisr-simulation-memos-overleaf".into());
 let h=hdf5::File::open(format!("{base}/data/msis_buoyancy.h5"))?;
 let z=h.dataset("altitude_km")?.read_raw::<f64>()?;
 let n=h.dataset("N_s-1")?.read_raw::<f64>()?;
 let mut f=Figure::new(IntendedUse::ArticleColumn);
 f.height_mm=90.;
 let mut ax=Axes::new([70.,150.],[0.01,0.035],r"Altitude $z\,[\mathrm{km}]$",r"$N\,[\mathrm{s}^{-1}]$");
 ax.xticks=vec![70.,90.,110.,130.,150.];ax.yticks=vec![0.01,0.02,0.03];
 let mut p=Panel::new(Rect::new(16.,11.,65.,51.),ax);
 let mut l=Line::new(vec![[70.,0.02],[150.,0.02]],"#777777");l.dash=Some("4 3".into());p.push(Layer::Line(l));
 let mut l=Line::new(z.iter().zip(n.iter()).map(|(&a,&b)|[a,b]).collect(),"#245b78");l.width_pt=1.3;p.push(Layer::Line(l));
 f.panels.push(p);
 f.labels.push(Label::new("NRLMSIS 2.1 over SYISR",48.5,4.));
 f.labels.push(Label::new(r"$N=0.02\,\mathrm{s}^{-1}$ (simulation)",48.,83.));
 f.save(format!("{base}/figures/msis_buoyancy.pdf"))?;
 f.save(format!("{base}/figures/msis_buoyancy.png"))?;
 Ok(())
}
