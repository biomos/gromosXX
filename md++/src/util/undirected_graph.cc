/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include "undirected_graph.h"
#include <limits>

namespace util {
	const UndirectedGraph::index_t UndirectedGraph::unexplored_id =
		std::numeric_limits<UndirectedGraph::index_t>::max();

	void UndirectedGraph::resize(const std::size_t& N) {
		m_vertices.resize(N);
	}

	bool UndirectedGraph::adjacent(const index_t& i, const index_t& j) const {
		for (edges_cont_t::const_iterator it = m_vertices[i].begin();
				it != m_vertices[i].end(); ++it) {
			if (*it == j)
				return true;
		}
		return false;
	}

	void UndirectedGraph::add_edge(const index_t& i, const index_t& j) {
		m_vertices[i].push_back(j);
		m_vertices[j].push_back(i);
	}

	void UndirectedGraph::remove_edge(const index_t& i, const index_t& j) {
		m_vertices[i].remove(j);
		m_vertices[j].remove(i);
	}

	UndirectedGraph::component_cont_t UndirectedGraph::connected_components(index_t& num_components) const {
		// set default value to maximum to test if undiscovered and component index starts from 0
		const vertices_cont_t::size_type N = number_of_vertices();
		component_cont_t components(N, unexplored_id);
		index_t component_index = 0;

		for (index_t i = 0; i < N; ++i) {
			if (components[i] == unexplored_id) {
				//start recurcsion
				explore_vertex(i, component_index, components);
				++component_index;
			}
		}
		num_components = component_index + 1;
		return components;
	}

	void UndirectedGraph::explore_vertex(const index_t& i, const index_t& component_index,
			component_cont_t& components) const{
		if (components[i] == unexplored_id) {
			components[i] = component_index;
			for (edges_cont_t::const_iterator it = m_vertices[i].begin();
					it != m_vertices[i].end(); ++it) {
				explore_vertex(*it, component_index, components);
			}
		}
	}
}
