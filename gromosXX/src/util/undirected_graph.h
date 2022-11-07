#include <list>
#include <vector>

namespace util {
	/**
	* @class UndirectedGraph
	* @ingroup util
	* 
	* Basic undirected graph implentation with "connected components"
	* capability
	*/
	class UndirectedGraph {
	public:
		/**
		 * typedefs
		 */
		typedef unsigned int index_t; // used to store incices
		typedef std::list<index_t> edges_cont_t; // connected vertices for one vertex
		typedef std::vector<edges_cont_t> vertices_cont_t; // storage of all vertices
		typedef std::vector<index_t> component_cont_t; // component index container

		/**
		 * Constructor
		 * @param[in] N number of vertices
		 */
		UndirectedGraph(const std::size_t& N) : m_vertices(N) {}

		/**
		 * resize graph to new number of vertices
		 * @param[in] N new number of vertices
		 */
		void resize(const std::size_t& N);

		/**
		 * checks if two vertices are connected by an edje
		 * @param[in] i index of vertex one
		 * @param[in] j index of vertex two
		 * @return true if adjacent, false otherwise
		 */
		bool adjacent(const index_t& i, const index_t& j) const;

		/**
		 * add edge between two vertices
		 * @param[in] i index of vertex one
		 * @param[in] j index of vertex two
		 */
		void add_edge(const index_t& i, const index_t& j);

		/**
		 * remove edge between two vertices
		 * @param[in] i index of vertex one
		 * @param[in] j index of vertex two
		 */
		void remove_edge(const index_t& i, const index_t& j);

		/**
		 * determines groups of vertices connected by edges
		 * @param[out] num_components number of detected components
		 * @return vector with group index for each vertex
		 */
		component_cont_t connected_components(index_t& num_components) const;

		/**
		 * determines groups of vertices connected by edges
		 * @return vector with group index for each vertex
		 */
		inline component_cont_t connected_components() const {
			index_t num_components = 0;
			return connected_components(num_components);
		}

		/**
		 * @return number of vertices
		 */
		inline vertices_cont_t::size_type number_of_vertices() const { return m_vertices.size(); }

	private:
		/**
		 * explores connected vertices recursively
		 * @param[in] i vertex index to explore
		 * @param[in] component_index index of component being explored
		 * @param[out] components vector, in which component_index is written to
		 */
		void explore_vertex(const index_t& i, const index_t& component_index,
				component_cont_t& components) const;

		/**
		 * vector of lists of connected vertices
		 */
		vertices_cont_t m_vertices;

		/**
		 * index used to indicate unexplored vertices
		 */
		static const index_t unexplored_id;
	};
}
