/// Replica Exchange Molecular Dynamics (REMD)
///
/// This module implements replica exchange methods for enhanced sampling:
/// - Temperature Replica Exchange (T-REMD)
/// - 2D Temperature-Lambda Replica Exchange (T-Lambda-REPEX)
///
/// Based on the Metropolis-Hastings criterion with Boltzmann-weighted
/// energy differences.

use crate::replica::{Replica, ReplicaId, ReplicaInfo};
use crate::topology::Topology;
use crate::configuration::Configuration;
use crate::integrator::Integrator;
use crate::algorithm::thermostats::BerendsenThermostatParameters;
use crate::algorithm::constraints::ShakeParameters;
use rand::Rng;
use rayon::prelude::*;
use std::sync::Arc;

/// Exchange type for replica exchange
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExchangeType {
    /// Temperature exchange only
    Temperature,
    /// Lambda exchange only (for FEP)
    Lambda,
    /// 2D temperature and lambda exchange
    TemperatureLambda,
}

/// Exchange scheme for determining partners
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExchangeScheme {
    /// Sequential exchange (0↔1, 2↔3, ...)
    Sequential,
    /// Odd-even exchange alternating (0↔1, 2↔3 then 1↔2, 3↔4)
    OddEven,
    /// Random pairs
    Random,
}

/// Statistics for replica exchange
#[derive(Debug, Clone, Default)]
pub struct ExchangeStatistics {
    /// Total number of exchange attempts
    pub attempts: usize,

    /// Number of successful exchanges
    pub accepted: usize,

    /// Acceptance rate
    pub acceptance_rate: f64,

    /// Per-replica pair acceptance counts [i][j] = exchanges between i and j
    pub pair_acceptances: Vec<Vec<usize>>,

    /// Per-replica pair attempt counts
    pub pair_attempts: Vec<Vec<usize>>,
}

impl ExchangeStatistics {
    /// Create new statistics for n replicas
    pub fn new(n_replicas: usize) -> Self {
        ExchangeStatistics {
            attempts: 0,
            accepted: 0,
            acceptance_rate: 0.0,
            pair_acceptances: vec![vec![0; n_replicas]; n_replicas],
            pair_attempts: vec![vec![0; n_replicas]; n_replicas],
        }
    }

    /// Update statistics after exchange attempt
    pub fn record_attempt(&mut self, id1: ReplicaId, id2: ReplicaId, accepted: bool) {
        self.attempts += 1;
        self.pair_attempts[id1][id2] += 1;
        self.pair_attempts[id2][id1] += 1;

        if accepted {
            self.accepted += 1;
            self.pair_acceptances[id1][id2] += 1;
            self.pair_acceptances[id2][id1] += 1;
        }

        self.acceptance_rate = self.accepted as f64 / self.attempts as f64;
    }

    /// Get acceptance rate for a specific pair
    pub fn pair_acceptance_rate(&self, id1: ReplicaId, id2: ReplicaId) -> f64 {
        if self.pair_attempts[id1][id2] == 0 {
            0.0
        } else {
            self.pair_acceptances[id1][id2] as f64 / self.pair_attempts[id1][id2] as f64
        }
    }
}

/// Replica Exchange Controller
///
/// Manages multiple replicas and coordinates exchange attempts
pub struct ReplicaController {
    /// All replicas
    pub replicas: Vec<Replica>,

    /// Shared topology
    topology: Arc<Topology>,

    /// Exchange type
    exchange_type: ExchangeType,

    /// Exchange scheme
    exchange_scheme: ExchangeScheme,

    /// Number of MD steps between exchange attempts
    exchange_interval: usize,

    /// Number of equilibration steps (no exchanges)
    equilibration_steps: usize,

    /// Current step counter
    current_step: usize,

    /// Whether to rescale velocities after temperature exchange
    rescale_velocities: bool,

    /// Exchange statistics
    pub statistics: ExchangeStatistics,

    /// Random number generator
    rng: rand::rngs::ThreadRng,

    /// Odd/even alternation state
    odd_even_state: bool,
}

impl ReplicaController {
    /// Create a new replica exchange controller
    pub fn new(
        replicas: Vec<Replica>,
        topology: Arc<Topology>,
        exchange_type: ExchangeType,
        exchange_scheme: ExchangeScheme,
        exchange_interval: usize,
        equilibration_steps: usize,
        rescale_velocities: bool,
    ) -> Self {
        let n_replicas = replicas.len();

        ReplicaController {
            replicas,
            topology,
            exchange_type,
            exchange_scheme,
            exchange_interval,
            equilibration_steps,
            current_step: 0,
            rescale_velocities,
            statistics: ExchangeStatistics::new(n_replicas),
            rng: rand::thread_rng(),
            odd_even_state: false,
        }
    }

    /// Get number of replicas
    pub fn n_replicas(&self) -> usize {
        self.replicas.len()
    }

    /// Run MD for all replicas in parallel
    pub fn run_md_parallel<I: Integrator + Send + Sync + Clone>(
        &mut self,
        integrator: &I,
        n_steps: usize,
        thermostat_params: Option<&BerendsenThermostatParameters>,
        shake_params: Option<&ShakeParameters>,
    ) {
        // Create thermostat parameters for each replica
        let thermostat_params_vec: Vec<Option<BerendsenThermostatParameters>> = self
            .replicas
            .iter()
            .map(|replica| {
                thermostat_params.map(|params| BerendsenThermostatParameters {
                    target_temperature: replica.info.temperature,
                    coupling_time: params.coupling_time,
                })
            })
            .collect();

        // Run MD for all replicas in parallel
        self.replicas
            .par_iter_mut()
            .zip(thermostat_params_vec.par_iter())
            .for_each(|(replica, thermo_params)| {
                let mut local_integrator = integrator.clone();
                replica.run_md_steps(
                    &self.topology,
                    &mut local_integrator,
                    n_steps,
                    thermo_params.as_ref(),
                    shake_params,
                );
            });

        self.current_step += n_steps;
    }

    /// Determine exchange pairs based on scheme
    fn determine_exchange_pairs(&mut self) -> Vec<(ReplicaId, ReplicaId)> {
        let n = self.n_replicas();
        let mut pairs = Vec::new();

        match self.exchange_scheme {
            ExchangeScheme::Sequential => {
                // Try all sequential pairs: (0,1), (2,3), (4,5), ...
                for i in (0..n - 1).step_by(2) {
                    pairs.push((i, i + 1));
                }
            }
            ExchangeScheme::OddEven => {
                // Odd-even scheme: alternate between even and odd starting pairs
                if self.odd_even_state {
                    // Even: (0,1), (2,3), (4,5), ...
                    for i in (0..n - 1).step_by(2) {
                        pairs.push((i, i + 1));
                    }
                } else {
                    // Odd: (1,2), (3,4), (5,6), ...
                    for i in (1..n - 1).step_by(2) {
                        pairs.push((i, i + 1));
                    }
                }
                self.odd_even_state = !self.odd_even_state;
            }
            ExchangeScheme::Random => {
                // Generate random non-overlapping pairs
                let mut indices: Vec<usize> = (0..n).collect();

                // Shuffle indices
                use rand::seq::SliceRandom;
                indices.shuffle(&mut self.rng);

                // Create pairs from shuffled indices
                for i in (0..n - 1).step_by(2) {
                    pairs.push((indices[i], indices[i + 1]));
                }
            }
        }

        pairs
    }

    /// Calculate exchange probability using Metropolis criterion
    ///
    /// For temperature exchange:
    ///   delta = (β1 - β2) * (E2 - E1)
    ///   where β = 1/(kB*T)
    ///
    /// For 2D temperature-lambda exchange:
    ///   delta = β1 * (E12 - E11) - β2 * (E22 - E21)
    ///   where Eij = energy of config i with parameters j
    ///
    /// Accept if exp(-delta) > random(0,1)
    pub fn calculate_exchange_probability(
        &self,
        replica1: &Replica,
        replica2: &Replica,
    ) -> (f64, f64) {
        let beta1 = replica1.info.beta();
        let beta2 = replica2.info.beta();

        let e1 = replica1.potential_energy();
        let e2 = replica2.potential_energy();

        let delta = match self.exchange_type {
            ExchangeType::Temperature => {
                // Temperature exchange only
                // delta = (β1 - β2) * (E2 - E1)
                (beta1 - beta2) * (e2 - e1)
            }
            ExchangeType::Lambda => {
                // Lambda exchange only (same temperature)
                // For FEP, need to recalculate energies with swapped lambdas
                // For now, simplified version
                // TODO: Implement proper lambda-dependent energy calculation
                (beta1 - beta2) * (e2 - e1)
            }
            ExchangeType::TemperatureLambda => {
                // 2D exchange: need to calculate E12 and E21
                // E11 = e1 (config 1, params 1)
                // E22 = e2 (config 2, params 2)
                // E12 = energy of config 1 with params 2
                // E21 = energy of config 2 with params 1
                //
                // For now, simplified without recalculating cross-energies
                // TODO: Implement proper cross-energy calculation
                beta1 * (e2 - e1) - beta2 * (e2 - e1)
            }
        };

        // Calculate acceptance probability
        let probability = if delta < 0.0 {
            1.0 // Always accept if energy difference favors exchange
        } else {
            (-delta).exp()
        };

        (probability, delta)
    }

    /// Attempt exchange between two replicas
    pub fn attempt_exchange(&mut self, id1: ReplicaId, id2: ReplicaId) -> bool {
        // Calculate exchange probability
        let (probability, _delta) = {
            let replica1 = &self.replicas[id1];
            let replica2 = &self.replicas[id2];
            self.calculate_exchange_probability(replica1, replica2)
        };

        // Generate random number
        let random = self.rng.gen::<f64>();

        // Accept exchange?
        let accepted = random < probability;

        // Update replica info
        self.replicas[id1].info.partner_id = Some(id2);
        self.replicas[id2].info.partner_id = Some(id1);
        self.replicas[id1].info.exchange_probability = probability;
        self.replicas[id2].info.exchange_probability = probability;
        self.replicas[id1].info.exchange_accepted = accepted;
        self.replicas[id2].info.exchange_accepted = accepted;

        if accepted {
            self.perform_exchange(id1, id2);
        }

        // Record statistics
        self.statistics.record_attempt(id1, id2, accepted);

        accepted
    }

    /// Perform the actual exchange
    fn perform_exchange(&mut self, id1: ReplicaId, id2: ReplicaId) {
        // Exchange configurations
        {
            let (left, right) = if id1 < id2 {
                let (left_slice, right_slice) = self.replicas.split_at_mut(id2);
                (&mut left_slice[id1], &mut right_slice[0])
            } else {
                let (left_slice, right_slice) = self.replicas.split_at_mut(id1);
                (&mut right_slice[0], &mut left_slice[id2])
            };

            left.exchange_configuration(right);
        }

        // Rescale velocities if temperature exchange
        if self.rescale_velocities && self.exchange_type != ExchangeType::Lambda {
            let temp1 = self.replicas[id1].info.temperature;
            let temp2 = self.replicas[id2].info.temperature;

            self.replicas[id1].rescale_velocities(&self.topology, temp1);
            self.replicas[id2].rescale_velocities(&self.topology, temp2);
        }
    }

    /// Attempt exchanges for all determined pairs
    pub fn attempt_all_exchanges(&mut self) {
        // Only attempt exchanges after equilibration
        if self.current_step < self.equilibration_steps {
            return;
        }

        // Determine pairs
        let pairs = self.determine_exchange_pairs();

        // Attempt exchanges for all pairs
        for (id1, id2) in pairs {
            self.attempt_exchange(id1, id2);
        }
    }

    /// Run replica exchange simulation
    pub fn run<I: Integrator + Send + Sync + Clone>(
        &mut self,
        integrator: &I,
        total_steps: usize,
        thermostat_params: Option<&BerendsenThermostatParameters>,
        shake_params: Option<&ShakeParameters>,
    ) {
        let mut steps_completed = 0;

        while steps_completed < total_steps {
            // Run MD for exchange_interval steps
            let steps_to_run = std::cmp::min(
                self.exchange_interval,
                total_steps - steps_completed,
            );

            self.run_md_parallel(integrator, steps_to_run, thermostat_params, shake_params);

            steps_completed += steps_to_run;

            // Attempt exchanges
            if steps_completed % self.exchange_interval == 0 {
                self.attempt_all_exchanges();
            }
        }
    }

    /// Get replica by ID
    pub fn get_replica(&self, id: ReplicaId) -> Option<&Replica> {
        self.replicas.get(id)
    }

    /// Get mutable replica by ID
    pub fn get_replica_mut(&mut self, id: ReplicaId) -> Option<&mut Replica> {
        self.replicas.get_mut(id)
    }

    /// Print exchange statistics
    pub fn print_statistics(&self) {
        println!("=== Replica Exchange Statistics ===");
        println!("Total attempts: {}", self.statistics.attempts);
        println!("Accepted: {}", self.statistics.accepted);
        println!(
            "Overall acceptance rate: {:.2}%",
            self.statistics.acceptance_rate * 100.0
        );
        println!("\nPer-pair acceptance rates:");

        for i in 0..self.n_replicas() {
            for j in (i + 1)..self.n_replicas() {
                let rate = self.statistics.pair_acceptance_rate(i, j);
                if self.statistics.pair_attempts[i][j] > 0 {
                    println!(
                        "  Replicas {} ↔ {}: {:.2}% ({}/{})",
                        i,
                        j,
                        rate * 100.0,
                        self.statistics.pair_acceptances[i][j],
                        self.statistics.pair_attempts[i][j]
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exchange_statistics() {
        let mut stats = ExchangeStatistics::new(4);

        stats.record_attempt(0, 1, true);
        stats.record_attempt(0, 1, false);
        stats.record_attempt(2, 3, true);

        assert_eq!(stats.attempts, 3);
        assert_eq!(stats.accepted, 2);
        assert!((stats.acceptance_rate - 2.0 / 3.0).abs() < 1e-10);

        let rate_01 = stats.pair_acceptance_rate(0, 1);
        assert!((rate_01 - 0.5).abs() < 1e-10);

        let rate_23 = stats.pair_acceptance_rate(2, 3);
        assert!((rate_23 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_metropolis_criterion() {
        // Test that probability is 1.0 when delta < 0
        let delta_negative = -5.0;
        let prob = if delta_negative < 0.0 {
            1.0
        } else {
            (-delta_negative).exp()
        };
        assert_eq!(prob, 1.0);

        // Test that probability is exp(-delta) when delta > 0
        let delta_positive = 2.0;
        let prob = if delta_positive < 0.0 {
            1.0
        } else {
            (-delta_positive).exp()
        };
        assert!((prob - (-2.0_f64).exp()).abs() < 1e-10);
    }
}
