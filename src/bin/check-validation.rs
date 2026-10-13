use anyhow::{Result, ensure};
use serde_json::{Value, json};
fn read(path: &str) -> Result<Value> {
    Ok(serde_json::from_slice(&std::fs::read(path)?)?)
}
fn main() -> Result<()> {
    let a = read("output/validation/poisson-1/validation.json")?;
    let b = read("output/validation/poisson-fine/validation.json")?;
    let ratio = a["continuum_poisson_max_error"].as_f64().unwrap()
        / b["continuum_poisson_max_error"].as_f64().unwrap();
    ensure!(
        (3.5..4.5).contains(&ratio),
        "Poisson refinement ratio {ratio}"
    );
    let mut mass: f64 = 0.;
    for ranks in [1, 2, 4] {
        let d = read(&format!("output/validation/ranks-{ranks}/diagnostics.json"))?;
        mass = mass.max(d["mass_relative_error"].as_f64().unwrap().abs());
    }
    ensure!(mass < 1e-12, "ion mass conservation {mass}");
    let initial = hdf5::File::open("output/validation/control/frame_00000.h5")?;
    let final_state = hdf5::File::open("output/validation/control/frame_00002.h5")?;
    let wind = final_state
        .dataset("w_m_s")?
        .read_raw::<f64>()?
        .iter()
        .fold(0_f64, |v, x| v.max(x.abs()));
    ensure!(
        wind < 1e-12,
        "seed-free control grew vertical velocity {wind}"
    );
    let n0 = initial.dataset("ion_density_m-3")?.read_raw::<f64>()?;
    let n1 = final_state.dataset("ion_density_m-3")?.read_raw::<f64>()?;
    let tracer_error = n0
        .iter()
        .zip(&n1)
        .fold(0_f64, |v, (a, b)| v.max((a - b).abs()))
        / 1e11;
    ensure!(
        tracer_error < 1e-12,
        "uniform horizontal advection changed horizontally uniform tracer"
    );
    let report = json!({"poisson_refinement_error_ratio":ratio,"maximum_ion_mass_relative_error":mass,"seed_free_max_vertical_wind_m_s":wind,"tracer_control_relative_error":tracer_error,"status":"numerical_checks_only_not_SYISR_validation"});
    std::fs::write(
        "output/validation/report.json",
        serde_json::to_string_pretty(&report)?,
    )?;
    println!("{report}");
    Ok(())
}
