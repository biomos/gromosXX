//! Logging and debugging utilities for GROMOS-RS
//!
//! Provides structured logging with different verbosity levels,
//! timing information, and debug output for troubleshooting.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

/// Global log level
static LOG_LEVEL: AtomicUsize = AtomicUsize::new(LogLevel::Info as usize);

/// Timing marker for performance measurement
static mut START_TIME: Option<Instant> = None;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(usize)]
pub enum LogLevel {
    Debug = 0,
    Info = 1,
    Warn = 2,
    Error = 3,
    Silent = 4,
}

impl LogLevel {
    pub fn from_verbosity(verbose: usize) -> Self {
        match verbose {
            0 => LogLevel::Info,
            1 => LogLevel::Debug,
            _ => LogLevel::Debug,
        }
    }
}

/// Set the global log level
pub fn set_log_level(level: LogLevel) {
    LOG_LEVEL.store(level as usize, Ordering::Relaxed);
}

/// Get the current log level
pub fn get_log_level() -> LogLevel {
    match LOG_LEVEL.load(Ordering::Relaxed) {
        0 => LogLevel::Debug,
        1 => LogLevel::Info,
        2 => LogLevel::Warn,
        3 => LogLevel::Error,
        _ => LogLevel::Silent,
    }
}

/// Start the global timer
pub fn start_timer() {
    unsafe {
        START_TIME = Some(Instant::now());
    }
}

/// Get elapsed time since timer start
pub fn elapsed_time() -> f64 {
    unsafe {
        START_TIME.map_or(0.0, |start| start.elapsed().as_secs_f64())
    }
}

/// Format elapsed time as string
pub fn format_time() -> String {
    let elapsed = elapsed_time();
    format!("[{:7.3}s]", elapsed)
}

/// Log a debug message
#[macro_export]
macro_rules! log_debug {
    ($($arg:tt)*) => {
        if $crate::logging::get_log_level() <= $crate::logging::LogLevel::Debug {
            eprintln!("{} [DEBUG] {}", $crate::logging::format_time(), format!($($arg)*));
        }
    };
}

/// Log an info message
#[macro_export]
macro_rules! log_info {
    ($($arg:tt)*) => {
        if $crate::logging::get_log_level() <= $crate::logging::LogLevel::Info {
            eprintln!("{} [INFO]  {}", $crate::logging::format_time(), format!($($arg)*));
        }
    };
}

/// Log a warning message
#[macro_export]
macro_rules! log_warn {
    ($($arg:tt)*) => {
        if $crate::logging::get_log_level() <= $crate::logging::LogLevel::Warn {
            eprintln!("{} [WARN]  {}", $crate::logging::format_time(), format!($($arg)*));
        }
    };
}

/// Log an error message
#[macro_export]
macro_rules! log_error {
    ($($arg:tt)*) => {
        if $crate::logging::get_log_level() <= $crate::logging::LogLevel::Error {
            eprintln!("{} [ERROR] {}", $crate::logging::format_time(), format!($($arg)*));
        }
    };
}

/// Format time duration
pub fn format_duration(seconds: f64) -> String {
    if seconds < 1.0 {
        format!("{:.1} ms", seconds * 1000.0)
    } else if seconds < 60.0 {
        format!("{:.2} s", seconds)
    } else if seconds < 3600.0 {
        format!("{:.1} min", seconds / 60.0)
    } else {
        format!("{:.2} h", seconds / 3600.0)
    }
}

/// Performance timer for code sections
pub struct Timer {
    name: String,
    start: Instant,
}

impl Timer {
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            start: Instant::now(),
        }
    }

    pub fn elapsed(&self) -> f64 {
        self.start.elapsed().as_secs_f64()
    }

    pub fn stop(&self) {
        let elapsed = self.elapsed();
        log_debug!("  {} took {}", self.name, format_duration(elapsed));
    }
}

impl Drop for Timer {
    fn drop(&mut self) {
        if get_log_level() <= LogLevel::Debug {
            self.stop();
        }
    }
}

/// Energy diagnostics
pub struct EnergyDiagnostics {
    pub step: usize,
    pub time: f64,
    pub kinetic: f64,
    pub potential: f64,
    pub bond: f64,
    pub angle: f64,
    pub dihedral: f64,
    pub improper: f64,
    pub lj: f64,
    pub coulomb: f64,
    pub temperature: f64,
}

impl EnergyDiagnostics {
    pub fn print(&self, verbose: bool) {
        if verbose {
            eprintln!("\n╔══════════════════════════════════════════════════════════════╗");
            eprintln!("║  Energy Diagnostics - Step {:6}  Time {:8.3} ps        ║",
                self.step, self.time);
            eprintln!("╠══════════════════════════════════════════════════════════════╣");
            eprintln!("║  Kinetic Energy:    {:12.4} kJ/mol                    ║", self.kinetic);
            eprintln!("║  Potential Energy:  {:12.4} kJ/mol                    ║", self.potential);
            eprintln!("║  Total Energy:      {:12.4} kJ/mol                    ║",
                self.kinetic + self.potential);
            eprintln!("╠══════════════════════════════════════════════════════════════╣");
            eprintln!("║  Bond:              {:12.4} kJ/mol                    ║", self.bond);
            eprintln!("║  Angle:             {:12.4} kJ/mol                    ║", self.angle);
            eprintln!("║  Dihedral:          {:12.4} kJ/mol                    ║", self.dihedral);
            eprintln!("║  Improper:          {:12.4} kJ/mol                    ║", self.improper);
            eprintln!("║  Lennard-Jones:     {:12.4} kJ/mol                    ║", self.lj);
            eprintln!("║  Coulomb:           {:12.4} kJ/mol                    ║", self.coulomb);
            eprintln!("╠══════════════════════════════════════════════════════════════╣");
            eprintln!("║  Temperature:       {:12.2} K                         ║", self.temperature);
            eprintln!("╚══════════════════════════════════════════════════════════════╝");
        } else {
            log_info!("Step {:6}  E_tot: {:10.4}  E_kin: {:10.4}  E_pot: {:10.4}  T: {:6.1} K",
                self.step, self.kinetic + self.potential, self.kinetic, self.potential, self.temperature);
        }
    }
}

/// Force diagnostics
pub struct ForceDiagnostics {
    pub max_force: f32,
    pub max_force_atom: usize,
    pub avg_force: f32,
    pub num_atoms: usize,
}

impl ForceDiagnostics {
    pub fn print(&self) {
        log_debug!("  Force statistics:");
        log_debug!("    Maximum force: {:.4} kJ/(mol·nm) on atom {}", self.max_force, self.max_force_atom + 1);
        log_debug!("    Average force: {:.4} kJ/(mol·nm)", self.avg_force);
        log_debug!("    Atoms: {}", self.num_atoms);
    }
}

/// Trajectory frame diagnostics
pub struct FrameDiagnostics {
    pub frame: usize,
    pub time: f64,
    pub num_atoms: usize,
    pub box_volume: f32,
    pub center_of_mass: (f32, f32, f32),
}

impl FrameDiagnostics {
    pub fn print(&self) {
        log_debug!("Frame {}: time={:.3} ps, atoms={}, volume={:.3} nm³, COM=({:.3}, {:.3}, {:.3})",
            self.frame, self.time, self.num_atoms, self.box_volume,
            self.center_of_mass.0, self.center_of_mass.1, self.center_of_mass.2);
    }
}

/// Progress bar for long calculations
pub struct ProgressBar {
    total: usize,
    current: usize,
    width: usize,
    last_print: Instant,
}

impl ProgressBar {
    pub fn new(total: usize) -> Self {
        Self {
            total,
            current: 0,
            width: 50,
            last_print: Instant::now(),
        }
    }

    pub fn update(&mut self, current: usize) {
        self.current = current;

        // Only print every 0.5 seconds to avoid spam
        if self.last_print.elapsed().as_secs_f64() > 0.5 || current == self.total {
            self.print();
            self.last_print = Instant::now();
        }
    }

    pub fn print(&self) {
        if get_log_level() > LogLevel::Info {
            return;
        }

        let progress = self.current as f64 / self.total as f64;
        let filled = (progress * self.width as f64) as usize;
        let bar: String = (0..self.width)
            .map(|i| if i < filled { '█' } else { '░' })
            .collect();

        eprint!("\r  [{}] {:>3.0}% ({}/{})",
            bar, progress * 100.0, self.current, self.total);

        if self.current == self.total {
            eprintln!();
        }
    }
}

impl Drop for ProgressBar {
    fn drop(&mut self) {
        if self.current == self.total && get_log_level() <= LogLevel::Info {
            eprintln!();
        }
    }
}
