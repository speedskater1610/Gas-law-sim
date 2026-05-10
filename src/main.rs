use rand::prelude::*;
use rand_distr::{Normal, Distribution};
use std::fmt;
use std::collections::VecDeque;

// CONSTANTS

const R: f64 = 0.082057; // L*atm/(mol*K)

// GAS SPECIES  (Van der Waals constants + molar mass)

#[derive(Debug, Clone, Copy)]
pub struct GasSpecies {
    pub name: &'static str,
    pub a: f64,
    pub b: f64,
    pub molar_mass: f64, // kg/mol
}

impl fmt::Display for GasSpecies {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} (a={:.3} L^2atm/mol^2, b={:.4} L/mol)", self.name, self.a, self.b)
    }
}

pub const GAS_IDEAL: GasSpecies = GasSpecies { name: "Ideal", a: 0.000, b: 0.0000, molar_mass: 0.028 };
pub const GAS_CO2:   GasSpecies = GasSpecies { name: "CO2",   a: 3.640, b: 0.0427, molar_mass: 0.044 };
pub const GAS_N2:    GasSpecies = GasSpecies { name: "N2",    a: 1.390, b: 0.0391, molar_mass: 0.028 };
pub const GAS_O2:    GasSpecies = GasSpecies { name: "O2",    a: 1.360, b: 0.0318, molar_mass: 0.032 };
pub const GAS_CL2:   GasSpecies = GasSpecies { name: "Cl2",   a: 6.580, b: 0.0562, molar_mass: 0.071 };
pub const GAS_HE:    GasSpecies = GasSpecies { name: "He",    a: 0.034, b: 0.0238, molar_mass: 0.004 };
pub const GAS_NH3:   GasSpecies = GasSpecies { name: "NH3",   a: 4.170, b: 0.0371, molar_mass: 0.017 };

pub const ALL_GASES: [GasSpecies; 7] = [GAS_IDEAL, GAS_CO2, GAS_N2, GAS_O2, GAS_CL2, GAS_HE, GAS_NH3];

// ideal gass law
pub fn ideal_pressure(n: f64, v: f64, t: f64) -> f64 {
    (n * R * t) / v
}

pub fn ideal_volume(n: f64, p: f64, t: f64) -> f64 {
    (n * R * t) / p
}

// van der waals equation   (P + a n^(2) / V^(2))(V - nb) = nRT
//   -> P = nRT/(V-nb) - a n^(2) / V^(2)

pub fn vdw_pressure(n: f64, v: f64, t: f64, gas: &GasSpecies) -> Option<f64> {
    let vm = v / n;
    if vm <= gas.b { return None; }
    Some((R * t) / (vm - gas.b) - gas.a / (vm * vm))
}

pub fn vdw_volume(n: f64, p: f64, t: f64, gas: &GasSpecies) -> f64 {
    let mut v = ideal_volume(n, p, t);
    for _ in 0..200 {
        let pv = match vdw_pressure(n, v, t, gas) { Some(x) => x, None => break };
        let err = pv - p;
        if err.abs() < 1e-10 { break; }
        let vm = v / n;
        let dp_dvm = -(R * t) / (vm - gas.b).powi(2) + 2.0 * gas.a / vm.powi(3);
        v -= err / (dp_dvm / n);
        v = v.max(n * gas.b * 1.001 + 1e-9);
    }
    v
}

pub fn compressibility_factor(p: f64, v: f64, n: f64, t: f64) -> f64 {
    (p * v) / (n * R * t)
}
// PARTICLE (3D pose + heading for MCL)
#[derive(Debug, Clone)]
pub struct Particle {
    pub pos: [f64; 3],     // position in box [0,1]^3
    pub vel: [f64; 3],     // velocity
    pub heading: [f64; 2], // azimuth, elevation (in radians)
    pub weight: f64,
}

impl Particle {
    pub fn random(rng: &mut impl Rng, speed_sigma: f64) -> Self {
        let normal = Normal::new(0.0, speed_sigma).unwrap();
        Particle {
            pos: [rng.gen(), rng.gen(), rng.gen()],
            vel: [normal.sample(rng), normal.sample(rng), normal.sample(rng)],
            heading: [
                rng.gen::<f64>() * 2.0 * std::f64::consts::PI,
                rng.gen::<f64>() * std::f64::consts::PI - std::f64::consts::FRAC_PI_2,
            ],
            weight: 1.0,
        }
    }

    pub fn speed(&self) -> f64 {
        self.vel.iter().map(|v| v * v).sum::<f64>().sqrt()
    }
}

/// MCL - monte carlo Localization
///
/// From paper §Methods:
///   "generates a large set of random particle positions and updates their
///    likelihood based on observed particle interactions and system constraints.
///    Over many iterations, particles that better match expected behavior are
///    weighted more heavily."

pub struct MCL {
    pub cloud: Vec<Particle>,
    pub n_particles: usize,
    pub target_speed: f64,
    pub noise_sigma: f64,
}

impl MCL {
    pub fn new(n_particles: usize, target_speed: f64, rng: &mut impl Rng) -> Self {
        let cloud = (0..n_particles).map(|_| Particle::random(rng, target_speed)).collect();
        MCL { cloud, n_particles, target_speed, noise_sigma: 0.01 }
    }

    /// Motion model: thermal random walk with elastic wall collisions
    pub fn predict(&mut self, rng: &mut impl Rng, dt: f64) {
        let vel_noise = Normal::new(0.0, self.noise_sigma).unwrap();
        let pos_noise = Normal::new(0.0, self.noise_sigma * dt * 0.1).unwrap();
        for p in &mut self.cloud {
            for v in &mut p.vel { *v += vel_noise.sample(rng); }
            let spd = p.speed();
            if spd > self.target_speed * 3.0 && spd > 1e-12 {
                let s = self.target_speed * 3.0 / spd;
                for v in &mut p.vel { *v *= s; }
            }
            for i in 0..3 {
                p.pos[i] += p.vel[i] * dt + pos_noise.sample(rng);
                if p.pos[i] < 0.0 { p.pos[i] = -p.pos[i]; p.vel[i] = p.vel[i].abs(); }
                if p.pos[i] > 1.0 { p.pos[i] = 2.0 - p.pos[i]; p.vel[i] = -p.vel[i].abs(); }
                p.pos[i] = p.pos[i].clamp(0.0, 1.0);
            }
            let [vx, vy, vz] = p.vel;
            p.heading[0] = vy.atan2(vx);
            p.heading[1] = (vz / (p.speed() + 1e-12)).asin();
        }
    }

    /// Observation model: weight by maxwell boltzmann speed likelihood
    /// and by pressure ratio consistency
    pub fn update_weights(&mut self, observed_speed: f64, pressure_ratio: f64) {
        let sigma_v = self.target_speed * 0.3;
        let sigma_p = 0.2;
        let mut w_sum = 0.0;
        for p in &mut self.cloud {
            let spd_err = (p.speed() - observed_speed) / sigma_v;
            let w_speed = (-0.5 * spd_err * spd_err).exp();
            let p_err = (p.weight.max(1e-300).ln() - pressure_ratio.max(1e-300).ln()) / sigma_p;
            let w_pressure = (-0.5 * p_err * p_err).exp();
            p.weight = (w_speed * w_pressure).max(1e-300);
            w_sum += p.weight;
        }
        if w_sum > 0.0 {
            for p in &mut self.cloud { p.weight /= w_sum; }
        } else {
            let w = 1.0 / self.n_particles as f64;
            for p in &mut self.cloud { p.weight = w; }
        }
    }

    /// Systematic resampling (triggered when N_eff < N/2)
    pub fn resample(&mut self, rng: &mut impl Rng) -> bool {
        let n_eff = self.effective_sample_size();
        if n_eff >= (self.n_particles as f64) / 2.0 { return false; }

        let mut cum = vec![0.0f64];
        for p in &self.cloud { cum.push(cum.last().unwrap() + p.weight); }

        let step = 1.0 / self.n_particles as f64;
        let start: f64 = rng.gen::<f64>() * step;
        let vel_noise = Normal::new(0.0, self.noise_sigma * 0.5).unwrap();
        let mut new_cloud = Vec::with_capacity(self.n_particles);
        let mut j = 0usize;

        for i in 0..self.n_particles {
            let target = start + i as f64 * step;
            while j < self.cloud.len() - 1 && cum[j + 1] < target { j += 1; }
            let mut clone = self.cloud[j].clone();
            for v in &mut clone.vel { *v += vel_noise.sample(rng); }
            clone.weight = step;
            new_cloud.push(clone);
        }
        self.cloud = new_cloud;
        true
    }

    pub fn effective_sample_size(&self) -> f64 {
        let sq: f64 = self.cloud.iter().map(|p| p.weight * p.weight).sum();
        if sq < 1e-300 { self.n_particles as f64 } else { 1.0 / sq }
    }

    pub fn estimated_speed(&self) -> f64 {
        self.cloud.iter().map(|p| p.weight * p.speed()).sum()
    }

    pub fn rms_speed(&self) -> f64 {
        let ms: f64 = self.cloud.iter().map(|p| p.weight * p.speed().powi(2)).sum();
        ms.sqrt()
    }
}

/// PID controller
///
/// From paper §Methods:
///   "proportional term (P) corrects deviations based on current errors
///    (calculated with VdW equation), integral (I) accounts for accumulated
///    errors, derivative (D) predicts future errors."
/// Tuning via Ziegler-Nichols (cited: Passaro 2021 PMID 34577364)

pub struct PidController {
    pub kp: f64,
    pub ki: f64,
    pub kd: f64,
    integral: f64,
    pub prev_error: f64,
    history: VecDeque<f64>,
}

impl PidController {
    pub fn new(kp: f64, ki: f64, kd: f64) -> Self {
        PidController { kp, ki, kd, integral: 0.0, prev_error: 0.0, history: VecDeque::new() }
    }

    // zeigler nichols PID:  Kp=0.6Ku, Ki=1.2Ku/Tu, Kd=0.075Ku*Tu
    // This PID was hand tuned until I got desired results
    pub fn ziegler_nichols(ku: f64, tu: f64) -> Self {
        Self::new(0.6 * ku, 1.2 * ku / tu, 0.075 * ku * tu)
    }

    pub fn compute(&mut self, setpoint: f64, measured: f64, dt: f64) -> f64 {
        let error = setpoint - measured;
        self.integral = (self.integral + error * dt).clamp(-10.0, 10.0);
        let derivative = (error - self.prev_error) / dt.max(1e-12);
        self.prev_error = error;
        self.history.push_back(error);
        if self.history.len() > 100 { self.history.pop_front(); }
        self.kp * error + self.ki * self.integral + self.kd * derivative
    }
}

// Extended kalman filter
//
// State:       x = [P (atm), T (K), rho (mol/L), v_rms (m/s)]
// f(x):        nonlinear (VdW based state transition)
// H = I4:      direct state observation

fn mat4_mul(a: &[f64; 16], b: &[f64; 16]) -> [f64; 16] {
    let mut c = [0.0f64; 16];
    for i in 0..4 { for k in 0..4 { for j in 0..4 {
        c[i*4+j] += a[i*4+k] * b[k*4+j];
    }}}
    c
}

fn mat4_add(a: &[f64; 16], b: &[f64; 16]) -> [f64; 16] {
    let mut c = [0.0f64; 16];
    for i in 0..16 { c[i] = a[i] + b[i]; }
    c
}

fn mat4_sub(a: &[f64; 16], b: &[f64; 16]) -> [f64; 16] {
    let mut c = [0.0f64; 16];
    for i in 0..16 { c[i] = a[i] - b[i]; }
    c
}

fn mat4_transpose(a: &[f64; 16]) -> [f64; 16] {
    let mut b = [0.0f64; 16];
    for i in 0..4 { for j in 0..4 { b[j*4+i] = a[i*4+j]; }}
    b
}

fn mat4_inv_approx(a: &[f64; 16]) -> [f64; 16] {
    // Diagonal approximation so valid when off diagonals are small
    let mut b = [0.0f64; 16];
    for i in 0..4 {
        let d = a[i*4+i];
        b[i*4+i] = if d.abs() > 1e-15 { 1.0 / d } else { 0.0 };
    }
    b
}

fn identity4() -> [f64; 16] {
    let mut m = [0.0f64; 16];
    for i in 0..4 { m[i*4+i] = 1.0; }
    m
}

#[derive(Debug, Clone)]
pub struct Ekf {
    pub x: [f64; 4],         // state estimate
    pub cov: [f64; 16],      // covariance P
    pub q: [f64; 16],        // process noise Q
    pub r_noise: [f64; 16],  // measurement noise R
    pub innovation: f64,
    pub kalman_gain: [f64; 4],
    pub nees: f64,
}

impl Ekf {
    pub fn new(p0: f64, t0: f64, rho0: f64, v0: f64) -> Self {
        let mut cov = [0.0f64; 16];
        cov[0] = 0.01; cov[5] = 100.0; cov[10] = 1e-4; cov[15] = 1e6;

        let mut q = [0.0f64; 16];
        q[0] = 1e-4; q[5] = 1.0; q[10] = 1e-6; q[15] = 1e3;

        let mut r = [0.0f64; 16];
        r[0] = 5e-3; r[5] = 4.0; r[10] = 1e-5; r[15] = 5e2;

        Ekf { x: [p0, t0, rho0, v0], cov, q, r_noise: r,
              innovation: 0.0, kalman_gain: [0.0; 4], nees: 0.0 }
    }

    /// Linearised state transition jacobian F = df/dx around state estimate
    // When Writing this for robotics the first time I was so confused about the math
    fn jacobian_f(&self, gas: &GasSpecies, n: f64, v: f64) -> [f64; 16] {
        let t = self.x[1];
        let vm = v / n;
        let dp_dt = if vm > gas.b { R / (vm - gas.b) } else { 0.0 };
        let dv_dt = if t > 0.0 { self.x[3] / (2.0 * t) } else { 0.0 };

        let mut f = identity4();
        f[0*4+0] = 0.9;             // P is auto regressive
        f[0*4+1] = dp_dt * 0.1;     // dP/dT coupling
        f[1*4+1] = 1.0;             // T held constant between steps
        f[2*4+2] = 1.0;             // rho = n/V
        f[3*4+1] = dv_dt * 0.05;    // dv/dT
        f[3*4+3] = 0.999;
        f
    }

    /// Nonlinear state transition f(x)
    fn state_transition(&self, gas: &GasSpecies, n: f64, v: f64) -> [f64; 4] {
        let [p, t, _rho, _v] = self.x;
        let p_vdw = vdw_pressure(n, v, t, gas).unwrap_or(p);
        let p_new = 0.9 * p + 0.1 * p_vdw;
        let t_new = t;
        let rho_new = n / v;
        // maxwell bboltzmann RMS speed v = sqrt(3 RT/M),  R in J/(mol*K) = 8.314
        // Zig comptime would have been nice for this
        let v_rms = (3.0 * 8.314 * t / gas.molar_mass).sqrt();
        [p_new, t_new, rho_new, v_rms]
    }

    /// EKF predict:  x^(-) = f(x),   P^(-) = F P F^(T) + Q
    pub fn predict(&mut self, gas: &GasSpecies, n: f64, v: f64) {
        let x_pred = self.state_transition(gas, n, v);
        let f  = self.jacobian_f(gas, n, v);
        let ft = mat4_transpose(&f);
        let fp   = mat4_mul(&f, &self.cov);
        let fpft = mat4_mul(&fp, &ft);
        self.cov = mat4_add(&fpft, &self.q);
        self.x = x_pred;
    }

    /// EKF update (H = I):  K = P S^(-1),  x^(+) = x^(-) + K ν,  P^(+) = (I-KH) P
    pub fn update(&mut self, z: [f64; 4]) {
        // innovation covariance S = P + R
        let s = mat4_add(&self.cov, &self.r_noise);
        let s_inv = mat4_inv_approx(&s);

        // kalman gain K = P S^(-1)
        let k = mat4_mul(&self.cov, &s_inv);
        for i in 0..4 { self.kalman_gain[i] = k[i*4+i]; }

        // innovation
        let innov: [f64; 4] = [z[0]-self.x[0], z[1]-self.x[1], z[2]-self.x[2], z[3]-self.x[3]];
        self.innovation = innov[0];

        // state update
        for i in 0..4 { self.x[i] += self.kalman_gain[i] * innov[i]; }

        // covariance update  P^(+) = (I - K) P
        let i_k = mat4_sub(&identity4(), &k);
        self.cov = mat4_mul(&i_k, &self.cov);

        // NEES = v^(t) S^(-1) ν
        self.nees = (0..4).map(|i| {
            let si_row: f64 = (0..4).map(|j| s_inv[i*4+j] * innov[j]).sum();
            innov[i] * si_row
        }).sum();
    }
}

// scenario & data point
#[derive(Debug, Clone)]
pub struct Scenario {
    pub name: String,
    pub gas: GasSpecies,
    pub temperature: f64,
    pub n_moles: f64,
    pub v_start: f64,
    pub v_end: f64,
    pub v_steps: usize,
}

#[derive(Debug, Clone)]
pub struct DataPoint {
    pub volume: f64,
    pub p_ideal: f64,
    pub p_vdw: f64,
    pub p_ekf: f64,
    pub z_ideal: f64,
    pub z_real: f64,
    pub z_ekf: f64,
    pub ekf_nees: f64,
    pub ekf_gain_p: f64,
    pub innovation: f64,
    pub mcl_n_eff: f64,
    pub pid_sigma: f64,
    pub resampled: bool,
}

// full Sim Engine

pub struct SimEngine {
    pub scenario: Scenario,
    pub rng: ThreadRng,
    pub mcl: MCL,
    pub pid: PidController,
    pub ekf: Ekf,
}

impl SimEngine {
    pub fn new(scenario: Scenario) -> Self {
        let mut rng = thread_rng();
        let gas = &scenario.gas;
        let t = scenario.temperature;
        let n = scenario.n_moles;
        let v = scenario.v_start;

        // maxwell boltzmann RMS speed (m/s)
        let v_rms = (3.0 * 8.314 * t / gas.molar_mass).sqrt();
        // Scale to sim units for MCL cloud (arbitrary, consistent)
        let target_sim = (v_rms * 0.0001).max(0.001);

        let mcl = MCL::new(300, target_sim, &mut rng);

        // ziegler nichols PID (Ku, Tu estimated from VdW nonlinearity)
        let pid = PidController::ziegler_nichols(2.0, 10.0);

        let p0   = vdw_pressure(n, v, t, gas).unwrap_or(ideal_pressure(n, v, t));
        let rho0 = n / v;
        let ekf  = Ekf::new(p0, t, rho0, v_rms);

        SimEngine { scenario, rng, mcl, pid, ekf }
    }

    /// Full P-V sweep: for each volume step, run INNER_STEPS of MCL & PID,
    /// then one EKF predict/update cycle.
    pub fn run(&mut self) -> Vec<DataPoint> {
        let sc = self.scenario.clone();
        let gas = sc.gas;
        const INNER: usize = 20;
        let dt = 0.05;

        // Log spaced volumes from v_start -> v_end
        let volumes: Vec<f64> = (0..=sc.v_steps).map(|i| {
            let frac = i as f64 / sc.v_steps as f64;
            ((1.0 - frac) * sc.v_start.ln() + frac * sc.v_end.ln()).exp()
        }).collect();

        let mut out = Vec::with_capacity(volumes.len());

        for &v in &volumes {
            let n = sc.n_moles;
            let t = sc.temperature;

            let p_ideal = ideal_pressure(n, v, t);
            let p_vdw   = vdw_pressure(n, v, t, &gas).unwrap_or(p_ideal);

            // Update MCL target speed for current T
            let v_rms_real = (3.0 * 8.314 * t / gas.molar_mass).sqrt();
            self.mcl.target_speed = (v_rms_real * 0.0001).max(0.001);

            let mut resampled = false;
            let mut n_eff = 0.0;

            for _ in 0..INNER {
                self.mcl.predict(&mut self.rng, dt);

                // Noisy observations
                let noise_p: f64 = Normal::new(0.0, 0.02).unwrap().sample(&mut self.rng);
                let p_ratio = (p_vdw / p_ideal.max(1e-10)) + noise_p;
                let obs_spd = self.mcl.target_speed * (1.0 + noise_p * 0.1);

                self.mcl.update_weights(obs_spd, p_ratio);
                if self.mcl.resample(&mut self.rng) { resampled = true; }
                n_eff = self.mcl.effective_sample_size();

                // PID: drive N_eff -> N / 2 by adjusting process noise sigma
                let pid_out = self.pid.compute(self.mcl.n_particles as f64 / 2.0, n_eff, dt);
                self.mcl.noise_sigma = (self.mcl.noise_sigma + pid_out * 0.001).clamp(0.001, 0.2);
            }

            // EKF step
            self.ekf.predict(&gas, n, v);

            let nm_v: f64 = Normal::new(0.0, 50.0).unwrap().sample(&mut self.rng);
            let z_meas = [
                p_vdw  + Normal::new(0.0, (p_vdw.abs() * 0.02).max(0.001)).unwrap().sample(&mut self.rng),
                t      + Normal::new(0.0, 2.0).unwrap().sample(&mut self.rng),
                n / v  + Normal::new(0.0, 0.0001).unwrap().sample(&mut self.rng),
                v_rms_real + nm_v,
            ];
            self.ekf.update(z_meas);

            let p_ekf = self.ekf.x[0];

            out.push(DataPoint {
                volume:     v,
                p_ideal,
                p_vdw,
                p_ekf,
                z_ideal: compressibility_factor(p_ideal, v, n, t),
                z_real:  compressibility_factor(p_vdw,   v, n, t),
                z_ekf:   compressibility_factor(p_ekf,   v, n, t),
                ekf_nees:   self.ekf.nees,
                ekf_gain_p: self.ekf.kalman_gain[0],
                innovation: self.ekf.innovation,
                mcl_n_eff:  n_eff,
                pid_sigma:  self.mcl.noise_sigma,
                resampled,
            });
        }
        out
    }
}

// ASCII Plot helpers

fn make_plot(width: usize, height: usize) -> Vec<Vec<char>> {
    vec![vec![' '; width]; height]
}

fn plot_series(grid: &mut Vec<Vec<char>>,
               xs: &[f64], ys: &[f64], ch: char,
               x_min: f64, x_max: f64, y_min: f64, y_max: f64) {
    let w = grid[0].len();
    let h = grid.len();
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        if !y.is_finite() { continue; }
        let col = ((x - x_min) / (x_max - x_min + 1e-12) * (w - 1) as f64) as usize;
        let row_f = (y - y_min) / (y_max - y_min + 1e-12) * (h - 1) as f64;
        let row = (h - 1).saturating_sub(row_f as usize);
        if col < w && row < h { grid[row][col] = ch; }
    }
}

fn draw_axes(grid: &mut Vec<Vec<char>>) {
    let h = grid.len();
    let w = grid[0].len();
    for r in 0..h { grid[r][0] = '|'; }
    for c in 0..w { grid[h-1][c] = '-'; }
    grid[h-1][0] = '+';
}

fn print_grid(grid: &Vec<Vec<char>>) {
    for row in grid {
        println!("    {}", row.iter().collect::<String>());
    }
}

fn hr(ch: char) {
    println!("{}", ch.to_string().repeat(72)); 
}

fn banner(s: &str) {
    println!();
    hr('=');
    let pad = (72usize.saturating_sub(s.len())) / 2;
    println!("{}{}", " ".repeat(pad), s);
    hr('=');
}

fn section(s: &str) {
    println!();
    println!("  -- {} {}", s, "-".repeat(60usize.saturating_sub(s.len() + 6)));
}

fn print_pv(data: &[DataPoint], title: &str) {
    section(&format!("P-V Diagram: {}", title));
    let vols: Vec<f64> = data.iter().map(|d| d.volume).collect();
    let pi:   Vec<f64> = data.iter().map(|d| d.p_ideal).collect();
    let pv:   Vec<f64> = data.iter().map(|d| d.p_vdw).collect();
    let pe:   Vec<f64> = data.iter().map(|d| d.p_ekf).collect();

    let v_min = vols.iter().cloned().fold(f64::INFINITY, f64::min);
    let v_max = vols.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let all_p: Vec<f64> = pi.iter().chain(pv.iter()).chain(pe.iter())
        .cloned().filter(|x| x.is_finite()).collect();
    let p_min = all_p.iter().cloned().fold(f64::INFINITY, f64::min).max(0.0);
    let p_max = all_p.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    let (w, h) = (62, 18);
    let mut g = make_plot(w, h);
    draw_axes(&mut g);
    // draw Z=1 reference ( horizontal at p=nRT/V so shown as dotted )
    plot_series(&mut g, &vols, &pi, '.', v_min, v_max, p_min, p_max);
    plot_series(&mut g, &vols, &pv, '#', v_min, v_max, p_min, p_max);
    plot_series(&mut g, &vols, &pe, 'x', v_min, v_max, p_min, p_max);

    println!("    P (atm)");
    println!("    {:>8.3} |", p_max);
    print_grid(&g);
    println!("    {:>8.3} +{}", p_min, "-".repeat(w));
    println!("    {:>10.4}{:>52.4}  V(L)", v_max, v_min);
    println!("    . = Ideal (PV=nRT)   # = VdW real gas   x = EKF estimate");
}

fn print_z(data: &[DataPoint], title: &str) {
    section(&format!("Compressibility Z: {}", title));
    let vols: Vec<f64> = data.iter().map(|d| d.volume).collect();
    let zr:   Vec<f64> = data.iter().map(|d| d.z_real.clamp(-1.0, 6.0)).collect();
    let ze:   Vec<f64> = data.iter().map(|d| d.z_ekf.clamp(-1.0, 6.0)).collect();

    let v_min = vols.iter().cloned().fold(f64::INFINITY, f64::min);
    let v_max = vols.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let z_min = 0.0f64;
    let z_max = zr.iter().chain(ze.iter()).cloned()
        .filter(|x| x.is_finite()).fold(1.5f64, f64::max);

    let (w, h) = (62, 14);
    let mut g = make_plot(w, h);
    draw_axes(&mut g);
    // Z=1 dashed line
    let z1_ys: Vec<f64> = vols.iter().map(|_| 1.0).collect();
    plot_series(&mut g, &vols, &z1_ys, '-', v_min, v_max, z_min, z_max);
    plot_series(&mut g, &vols, &zr, '#', v_min, v_max, z_min, z_max);
    plot_series(&mut g, &vols, &ze, 'x', v_min, v_max, z_min, z_max);

    println!("    Z");
    println!("    {:>6.3} |", z_max);
    print_grid(&g);
    println!("    {:>6.3} +{}", z_min, "-".repeat(w));
    println!("    {:>8.4}{:>54.4}  V(L)", v_max, v_min);
    println!("    - = Z=1 (ideal)   # = VdW real gas   x = EKF estimate");
}

fn print_table(data: &[DataPoint]) {
    section("Numerical Results (sampled rows)");
    println!();
    println!("  {:>7}  {:>9}  {:>9}  {:>9}  {:>7}  {:>7}  {:>7}  {:>7}  {:>7}  Resamp",
             "V(L)", "P_ideal", "P_VdW", "P_EKF", "Z_ideal", "Z_VdW", "Z_EKF", "NEES", "N_eff");
    println!("  {}", "-".repeat(87));
    let step = (data.len() / 14).max(1);
    for (i, d) in data.iter().enumerate() {
        if i % step == 0 || i == data.len() - 1 {
            println!("  {:>7.4}  {:>9.3}  {:>9.3}  {:>9.3}  {:>7.4}  {:>7.4}  {:>7.4}  {:>7.3}  {:>7.1}  {}",
                d.volume, d.p_ideal, d.p_vdw, d.p_ekf,
                d.z_ideal.min(9.9999), d.z_real.min(9.9999), d.z_ekf.min(9.9999),
                d.ekf_nees, d.mcl_n_eff,
                if d.resampled { "YES *" } else { "no" });
        }
    }
    println!("  * = MCL systematic resample triggered");
}

fn print_summary(data: &[DataPoint], sc: &Scenario) {
    section("Summary");
    let n = data.len();
    if n == 0 { return; }

    let avg_z: f64 = data.iter().filter(|d| d.z_real.is_finite())
        .map(|d| d.z_real).sum::<f64>() / n as f64;
    let max_dp: f64 = data.iter().map(|d|
        ((d.p_vdw - d.p_ideal) / d.p_ideal.abs().max(1e-12) * 100.0).abs()
    ).fold(0.0f64, f64::max);
    let first = &data[0];
    let last  = &data[n-1];

    println!("  Gas:          {}", sc.gas);
    println!("  Temperature:  {} K    n = {} mol    V: {:.3} -> {:.3} L",
             sc.temperature, sc.n_moles, sc.v_start, sc.v_end);
    println!();
    println!("  Mean Z (real gas):       {:.5}", avg_z);
    println!("  Z deviation at min V:    {:.5}  ({:.2}% from Z=1)",
             last.z_real - 1.0, (last.z_real - 1.0).abs() * 100.0);
    println!("  Max |delta_P| vs ideal:  {:.2}%", max_dp);
    println!();
    println!("  HIGH volume:   Z_real = {:.5}   Z_EKF = {:.5}", first.z_real, first.z_ekf);
    println!("  LOW  volume:   Z_real = {:.5}   Z_EKF = {:.5}", last.z_real,  last.z_ekf);
    println!();
    println!("  EKF final state: P={:.4} atm  T={:.1} K  rho={:.5} mol/L  v_rms={:.1} m/s",
             last.p_ekf, sc.temperature, sc.n_moles/last.volume, last.p_ekf);
    println!("  EKF NEES={:.4}  |innovation|={:.4}  K_P={:.4}",
             last.ekf_nees, last.innovation.abs(), last.ekf_gain_p);
}

fn main() {
    banner("GAS LAW SIMULATOR");
    println!();
    println!("  Physics:     Ideal (PV=nRT)  vs  Van der Waals real gas");
    println!("  Algorithms:  MCL particle filter, Ziegler-Nichols PID,");
    println!("               Extended Kalman Filter (EKF)");
    println!("  State:       x = [P, T, rho, v_rms]^T  (4-D)");
    println!("  Gases:       CO2, N2, O2, Cl2, He, NH3, Ideal");

    // Low T, High P - CO2
    {
        banner("SCENARIO 1: Low T / High P  [T=50K, n=5mol, V: 5->0.1 L, CO2]");
        println!("  Hypothesis: Z << 1 at low V; large P deviation from ideal.");
        println!("  Paper result: ~18% deviation in P; Z well below 1 at low V.");
        let sc = Scenario {
            name: "Low-T High-P CO2".into(), gas: GAS_CO2,
            temperature: 50.0, n_moles: 5.0,
            v_start: 5.0, v_end: 0.1, v_steps: 60,
        };
        let mut eng = SimEngine::new(sc.clone());
        let data = eng.run();
        print_pv(&data, &sc.name);
        print_z(&data, &sc.name);
        print_table(&data);
        print_summary(&data, &sc);
    }

    // High T, Low P - CO2
    {
        banner("SCENARIO 2: High T / Low P  [T=1000K, n=0.4mol, V: 20->0.5 L, CO2]");
        println!("  Hypothesis: Z ~= 1 throughout; deviation < 2%.");
        println!("  Paper result: 0.6% Z deviation; lines nearly coincide.");
        let sc = Scenario {
            name: "High-T Low-P CO2".into(), gas: GAS_CO2,
            temperature: 1000.0, n_moles: 0.4,
            v_start: 20.0, v_end: 0.5, v_steps: 60,
        };
        let mut eng = SimEngine::new(sc.clone());
        let data = eng.run();
        print_pv(&data, &sc.name);
        print_z(&data, &sc.name);
        print_table(&data);
        print_summary(&data, &sc);
    }

    // differnet gas comparison 
    {
        banner("SCENARIO 3: Cross-Gas Comparison  [T=200K, n=2mol, V=0.5L]");
        println!("  Ranked by Z deviation; Cl2/NH3 largest, He near-ideal.");
        println!();
        println!("  {:>8}  {:>7}  {:>7}  {:>10}  {:>10}  {:>8}",
                 "Gas", "a", "b", "P_ideal", "P_VdW", "Z_real");
        println!("  {}", "-".repeat(58));
        let t = 200.0; let n = 2.0; let v = 0.5;
        let mut rows: Vec<(GasSpecies, f64, f64, f64)> = ALL_GASES.iter().map(|gas| {
            let pi = ideal_pressure(n, v, t);
            let pv = vdw_pressure(n, v, t, gas).unwrap_or(pi);
            let z  = compressibility_factor(pv, v, n, t);
            (*gas, pi, pv, z)
        }).collect();
        rows.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap());
        for (gas, pi, pv, z) in &rows {
            println!("  {:>8}  {:>7.3}  {:>7.4}  {:>10.3}  {:>10.3}  {:>8.5}",
                     gas.name, gas.a, gas.b, pi, pv, z);
        }
    }

    // EKF convergence / N2 moderate conditions
    {
        banner("SCENARIO 4: EKF Convergence Detail  [N2, T=400K, n=1mol, V: 2->0.2L]");
        let sc = Scenario {
            name: "N2 EKF Tracking".into(), gas: GAS_N2,
            temperature: 400.0, n_moles: 1.0,
            v_start: 2.0, v_end: 0.2, v_steps: 40,
        };
        let mut eng = SimEngine::new(sc.clone());
        let data = eng.run();
        section("EKF State Tracking vs Volume");
        println!();
        println!("  {:>7}  {:>9}  {:>9}  {:>9}  {:>8}  {:>8}  {:>8}  {:>7}",
                 "V(L)", "P_VdW", "P_EKF", "Innov", "NEES", "K_P", "N_eff", "PID_sig");
        println!("  {}", "-".repeat(74));
        let step = (data.len() / 12).max(1);
        for (i, d) in data.iter().enumerate() {
            if i % step == 0 || i == data.len() - 1 {
                println!("  {:>7.4}  {:>9.4}  {:>9.4}  {:>9.4}  {:>8.4}  {:>8.5}  {:>8.1}  {:>7.4}",
                    d.volume, d.p_vdw, d.p_ekf, d.innovation,
                    d.ekf_nees, d.ekf_gain_p, d.mcl_n_eff, d.pid_sigma);
            }
        }
        print_summary(&data, &sc);
    }
}
