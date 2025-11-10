//! Generate GROMOS simulation scripts and input files
//!
//! This is a Rust port of GROMOS++ mk_script functionality,
//! focusing on IMD file generation and basic validation.
//!
//! Usage:
//!   mk_script --topo system.top --coord system.cnf --output sim.imd [options]

use std::path::PathBuf;
use std::fs::File;
use std::io::Write;
use clap::Parser;

/// Generate GROMOS simulation scripts and input files (GROMOS++ mk_script port)
#[derive(Parser, Debug)]
#[command(name = "mk_script")]
#[command(about = "Generate GROMOS simulation scripts and input files", long_about = None)]
struct Args {
    /// Topology file (.top)
    #[arg(short = 't', long)]
    topo: PathBuf,

    /// Coordinate file (.cnf)
    #[arg(short = 'c', long)]
    coord: PathBuf,

    /// Output IMD file
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Simulation title
    #[arg(long, default_value = "GROMOS simulation")]
    title: String,

    /// Number of MD steps
    #[arg(long, default_value_t = 10000)]
    nstlim: usize,

    /// Time step (ps)
    #[arg(long, default_value_t = 0.002)]
    dt: f64,

    /// Initial temperature (K)
    #[arg(long, default_value_t = 300.0)]
    temp: f64,

    /// Temperature coupling time (ps)
    #[arg(long, default_value_t = 0.1)]
    tau_t: f64,

    /// Pressure coupling (0=no, 1=isotropic, 2=anisotropic)
    #[arg(long, default_value_t = 0)]
    pcouple: i32,

    /// Reference pressure (bar)
    #[arg(long, default_value_t = 1.01325)]
    pres: f64,

    /// Pressure coupling time (ps)
    #[arg(long, default_value_t = 0.5)]
    tau_p: f64,

    /// SHAKE constraints (1=none, 2=H-bonds, 3=all)
    #[arg(long, default_value_t = 2)]
    ntc: i32,

    /// Long-range electrostatics (0=cutoff, 1=RF, 2=PME)
    #[arg(long, default_value_t = 1)]
    nlrele: i32,

    /// Reaction field permittivity (0.0 = infinity)
    #[arg(long, default_value_t = 0.0)]
    epsrf: f64,

    /// Short-range cutoff (nm)
    #[arg(long, default_value_t = 0.8)]
    rcutp: f64,

    /// Long-range cutoff (nm)
    #[arg(long, default_value_t = 1.4)]
    rcutl: f64,

    /// Trajectory output frequency (steps)
    #[arg(long, default_value_t = 100)]
    ntwx: usize,

    /// Energy output frequency (steps)
    #[arg(long, default_value_t = 10)]
    ntwe: usize,

    /// Random seed
    #[arg(long, default_value_t = 12345)]
    seed: i64,

    /// Initial velocities (0=read, 1=generate)
    #[arg(long, default_value_t = 1)]
    ntivel: i32,
}

fn main() {
    let args = Args::parse();

    // Validate input files exist
    if !args.topo.exists() {
        eprintln!("Error: Topology file not found: {}", args.topo.display());
        std::process::exit(1);
    }

    if !args.coord.exists() {
        eprintln!("Error: Coordinate file not found: {}", args.coord.display());
        std::process::exit(1);
    }

    println!("Generating IMD file:");
    println!("  Topology:    {}", args.topo.display());
    println!("  Coordinates: {}", args.coord.display());
    println!("  Output:      {}", args.output.display());
    println!();

    // Generate IMD file
    match generate_imd(&args) {
        Ok(()) => {
            println!("✓ Successfully generated: {}", args.output.display());
            println!();
            print_summary(&args);
        }
        Err(e) => {
            eprintln!("Error generating IMD file: {}", e);
            std::process::exit(1);
        }
    }
}

fn generate_imd(args: &Args) -> std::io::Result<()> {
    let mut file = File::create(&args.output)?;

    // TITLE block
    writeln!(file, "TITLE")?;
    writeln!(file, "  {}", args.title)?;
    writeln!(file, "END")?;
    writeln!(file)?;

    // SYSTEM block (would parse from topology in full implementation)
    writeln!(file, "SYSTEM")?;
    writeln!(file, "#  NPM    NSM")?;
    writeln!(file, "     1      0")?;
    writeln!(file, "END")?;
    writeln!(file)?;

    // STEP block
    writeln!(file, "STEP")?;
    writeln!(file, "#  NSTLIM      T        DT")?;
    writeln!(file, "  {:>7}    0.0    {:>8.5}", args.nstlim, args.dt)?;
    writeln!(file, "END")?;
    writeln!(file)?;

    // BOUNDCOND block
    writeln!(file, "BOUNDCOND")?;
    writeln!(file, "#   NTB  NDFMIN")?;
    writeln!(file, "      1       0")?;
    writeln!(file, "END")?;
    writeln!(file)?;

    // MULTIBATH block (temperature coupling)
    writeln!(file, "MULTIBATH")?;
    writeln!(file, "#  ALGORITHM")?;
    writeln!(file, "          1")?;
    writeln!(file, "#  NUM")?;
    writeln!(file, "      1")?;
    writeln!(file, "#  TEMP0         TAU")?;
    writeln!(file, "  {:>6.2}      {:>6.4}", args.temp, args.tau_t)?;
    writeln!(file, "#  DOFSET: num, last_atom_index")?;
    writeln!(file, "      1     0")?;
    writeln!(file, "END")?;
    writeln!(file)?;

    // PRESSURESCALE block (if requested)
    if args.pcouple > 0 {
        writeln!(file, "PRESSURESCALE")?;
        writeln!(file, "#  COUPLE  SCALE  COMP       TAUP   VIRIAL")?;
        writeln!(file, "        {}      1     1   {:>8.5}        2", args.pcouple, args.tau_p)?;
        writeln!(file, "#  SEMIANISOTROPIC COUPLINGS(X, Y, Z)")?;
        writeln!(file, "      1     1     1")?;
        writeln!(file, "#  PRES0(1...3,1...3)")?;
        writeln!(file, "  {:>10.5}    0.00000    0.00000", args.pres)?;
        writeln!(file, "    0.00000  {:>10.5}    0.00000", args.pres)?;
        writeln!(file, "    0.00000    0.00000  {:>10.5}", args.pres)?;
        writeln!(file, "END")?;
        writeln!(file)?;
    }

    // FORCE block
    writeln!(file, "FORCE")?;
    writeln!(file, "# NTF array")?;
    writeln!(file, "# 1: bonds, 2: angles, 3: imp. dihedrals, 4: dihedrals,")?;
    writeln!(file, "# 5: electrostatic, 6: Lennard-Jones")?;
    writeln!(file, "#  bonds  angles  impdih  dih     ele     LJ")?;
    if args.ntc >= 2 {
        writeln!(file, "      0      1       1      1       1      1")?;
    } else {
        writeln!(file, "      1      1       1      1       1      1")?;
    }
    writeln!(file, "# NRE(1..2)")?;
    writeln!(file, "#  NRE: last_atom_index(1), last_atom_index(2)")?;
    writeln!(file, "    0       0")?;
    writeln!(file, "END")?;
    writeln!(file)?;

    // CONSTRAINT block
    writeln!(file, "CONSTRAINT")?;
    writeln!(file, "#  NTC  NTCP  NTCP0(1)  NTCS  NTCS0(1)")?;
    writeln!(file, "    {}     0         0      1         0", args.ntc)?;
    writeln!(file, "END")?;
    writeln!(file)?;

    // PAIRLIST block
    writeln!(file, "PAIRLIST")?;
    writeln!(file, "#  ALGORITHM  NSNB  RCUTP   RCUTL     SIZE  TYPE")?;
    writeln!(file, "          0     5  {:>5.2}   {:>5.2}    {:>5.2}     0",
             args.rcutp, args.rcutl, args.rcutp / 2.0)?;
    writeln!(file, "END")?;
    writeln!(file)?;

    // NONBONDED block
    writeln!(file, "NONBONDED")?;
    writeln!(file, "#  NLRELE  APPAK    RCRF   EPSRF  NSLFEXCL")?;
    writeln!(file, "       {}    0.0   {:>5.2}  {:>6.1}         1",
             args.nlrele, args.rcutl, args.epsrf)?;
    writeln!(file, "END")?;
    writeln!(file)?;

    // INITIALISE block
    writeln!(file, "INITIALISE")?;
    writeln!(file, "#  NTIVEL  NTISHK  NTINHT  NTINHB  NTISHI     NTIRTC  NTICOM")?;
    writeln!(file, "        {}       0       0       0       0          0       0", args.ntivel)?;
    writeln!(file, "#  NTIR    NTIG      IG     TEMPI")?;
    writeln!(file, "      0       0  {:>6}  {:>6.2}", args.seed, args.temp)?;
    writeln!(file, "END")?;
    writeln!(file)?;

    // WRITETRAJ block
    writeln!(file, "WRITETRAJ")?;
    writeln!(file, "#  NTWX  NTWE  NTWV  NTWF  NTWE_SPEC")?;
    writeln!(file, "  {:>5} {:>5}     0     0          0", args.ntwx, args.ntwe)?;
    writeln!(file, "END")?;
    writeln!(file)?;

    // PRINTOUT block
    writeln!(file, "PRINTOUT")?;
    writeln!(file, "#  NTPR")?;
    writeln!(file, "  {:>5}", args.ntwe)?;
    writeln!(file, "END")?;
    writeln!(file)?;

    Ok(())
}

fn print_summary(args: &Args) {
    println!("Simulation Parameters:");
    println!("  Steps:       {}", args.nstlim);
    println!("  Time step:   {} ps", args.dt);
    println!("  Total time:  {:.2} ps", args.nstlim as f64 * args.dt);
    println!();
    println!("Temperature:");
    println!("  Target:      {} K", args.temp);
    println!("  Coupling:    {} ps", args.tau_t);
    println!();
    if args.pcouple > 0 {
        println!("Pressure:");
        println!("  Target:      {} bar", args.pres);
        println!("  Coupling:    {} ps", args.tau_p);
        println!("  Type:        {}", match args.pcouple {
            1 => "Isotropic",
            2 => "Anisotropic",
            _ => "Unknown"
        });
        println!();
    }
    println!("Constraints:");
    println!("  Type:        {}", match args.ntc {
        1 => "None",
        2 => "H-bonds",
        3 => "All bonds",
        _ => "Unknown"
    });
    println!();
    println!("Electrostatics:");
    println!("  Method:      {}", match args.nlrele {
        0 => "Cutoff",
        1 => "Reaction Field",
        2 => "PME",
        _ => "Unknown"
    });
    if args.nlrele == 1 {
        println!("  RF epsilon:  {}", if args.epsrf == 0.0 {
            "infinity (conducting)".to_string()
        } else {
            args.epsrf.to_string()
        });
    }
    println!("  Cutoff:      {} nm", args.rcutl);
    println!();
    println!("Output:");
    println!("  Trajectory:  every {} steps", args.ntwx);
    println!("  Energy:      every {} steps", args.ntwe);
}
