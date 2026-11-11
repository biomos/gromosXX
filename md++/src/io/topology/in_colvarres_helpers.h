/**
 * @file in_colvarres_helpers.h
 * Internal helpers shared by the collective-variable input parsers.
 * Not part of the installed reader interface.
 */
#ifndef INCLUDED_IN_COLVARRES_HELPERS_H
#define INCLUDED_IN_COLVARRES_HELPERS_H

#include "../../stdheader.h"
#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"

namespace io {
namespace colvar_input {

bool line_is_data(const std::string &line);

std::string colvar_type_from_int(int type);

void validate_rah(int rah, const std::string &type,
                  const std::string &blockname);

std::vector<std::string> tokens_from_line(const std::string &line);

bool decode_method(int raw,
                   int &averaging,
                   int &force_scale,
                   bool &constraining,
                   const std::string &type,
                   const std::string &blockname);

const simulation::Parameter::colvar_bias_spec *indexed_colvar_spec(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    const std::string &type,
    unsigned int index);

bool colvar_index_exists(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    unsigned int index);

const simulation::Parameter::colvar_bias_spec *nth_colvar_spec(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    const std::string &type,
    size_t index);

const simulation::Parameter::colvar_bias_spec *nth_colvar_spec(
    const simulation::Simulation &sim,
    const std::string &type,
    size_t index);

unsigned int count_colvar_specs(const simulation::Simulation &sim,
                                const std::string &type);

unsigned int count_colvar_specs(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    const std::string &type);

void make_coordnum_virtual_atoms(
    const std::vector<std::vector<int> > &atomgroup,
    const std::vector<int> &vatypes,
    std::vector<std::vector<util::Virtual_Atom> > &atoms);

void read_distance_difference_geometries(const std::vector<std::string> &buffer,
    const std::string &blockname,
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    bool perturbed, topology::Topology &topo, std::set<std::string> &block_read);

bool matching_wall_settings(const simulation::Parameter::colvar_bias_spec &a,
                            const simulation::Parameter::colvar_bias_spec &b);

void pair_angular_walls(std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    const std::vector<simulation::Parameter::colvar_bias_spec> &ordinary_specs,
    const topology::Topology &topo, bool perturbed);

void read_colvar_branches(const std::map<std::string, std::vector<std::string> > &blocks,
    simulation::Simulation &sim, std::set<std::string> &read, std::ostream &os, bool quiet);

void read_eds_colvar_blocks(const std::map<std::string, std::vector<std::string> > &blocks,
    topology::Topology &topo, simulation::Simulation &sim, std::set<std::string> &block_read);

} // namespace colvar_input
} // namespace io

#endif
