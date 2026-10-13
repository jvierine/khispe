use anyhow::Result;
fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let radar = a.iter().any(|s| s == "--radar");
    let visualization = radar || a.iter().any(|s| s == "--visualization");
    let mut paths: Vec<_> = std::fs::read_dir(&a[1])?
        .filter_map(|v| v.ok().map(|v| v.path()))
        .filter(|v| {
            v.file_name()
                .unwrap()
                .to_string_lossy()
                .starts_with("frame_")
                && v.extension().is_some_and(|v| v == "h5")
        })
        .collect();
    paths.sort();
    let out = hdf5::File::create(&a[2])?;
    for path in &paths {
        let src = hdf5::File::open(path)?;
        let g = out.create_group(&path.file_stem().unwrap().to_string_lossy())?;
        let t = src.attr("time_s")?.read_scalar::<f64>()?;
        g.new_attr::<f64>().create("time_s")?.write_scalar(&t)?;
        for name in [
            "ion_density_m-3",
            "omega_s-1",
            "u_m_s",
            "w_m_s",
            "altitude_m",
            "x_m",
        ] {
            if !radar && (name == "u_m_s" || name == "w_m_s") {
                continue;
            }
            if visualization && name == "omega_s-1" {
                continue;
            }
            let d = src.dataset(name)?;
            let v = d.read_raw::<f64>()?;
            if visualization && name == "ion_density_m-3" {
                let q: Vec<u16> = v
                    .iter()
                    .map(|x| (x.max(1.).log10() * 4096.).round().clamp(0., 65535.) as u16)
                    .collect();
                g.new_dataset::<u16>()
                    .shape(d.shape())
                    .deflate(3)
                    .create("log10_ne_quantized")?
                    .write_raw(&q)?;
                g.new_attr::<f64>()
                    .create("log10_density_quantum")?
                    .write_scalar(&(1. / 4096.))?;
            } else if name == "ion_density_m-3"
                || name == "omega_s-1"
                || name == "u_m_s"
                || name == "w_m_s"
            {
                let q: Vec<f32> = v.iter().map(|x| *x as f32).collect();
                g.new_dataset::<f32>()
                    .shape(d.shape())
                    .deflate(3)
                    .create(name)?
                    .write_raw(&q)?;
            } else {
                g.new_dataset::<f64>()
                    .shape(d.shape())
                    .deflate(3)
                    .create(name)?
                    .write_raw(&v)?;
            }
        }
        let cfg = src.dataset("config_json")?.read_raw::<u8>()?;
        g.new_dataset_builder()
            .with_data(&cfg)
            .create("config_json")?;
    }
    println!("Packed {} physical density frames", paths.len());
    Ok(())
}
