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

/**
 * @file algorithm_sequence.h
 * the Algorithm class.
 */

#ifndef INCLUDED_ALGORITHM_ALGORITHM_H
#define INCLUDED_ALGORITHM_ALGORITHM_H

namespace algorithm
{
  /**
   * @class Algorithm_Sequence
   * contains the specific algorithms.
   * An almost clean implementation of the Strategy pattern.
   */
  class Algorithm_Sequence : public std::vector<IAlgorithm *>
  {
  public:
    /**
     * Constructor
     */
    Algorithm_Sequence(bool clean = true);

    /**
     * Destructor
     */
    ~Algorithm_Sequence();

    /**
     * init
     */
    int init(topology::Topology &topo, 
	     configuration::Configuration &conf,
	     simulation::Simulation &sim,
	     std::ostream & os = std::cout,
	     bool quiet = false);

    /**
     * calculate all interactions.
     */
    int run(topology::Topology &topo, 
	    configuration::Configuration &conf,
	    simulation::Simulation &sim);

    /**
     * print timing information
     */
    int print_timing(std::ostream & os);

    /**
     * @brief algorithm accessor - templated version
     *  - casts directly
     *  - performs runtime casting check (dynamic_cast) in debug mode
     *  - NO check (static_cast) in release mode!
     * 
     * @tparam T type of the returned pointer
     * @tparam FailIfNotFound bool fail if not found (otherwise return nullptr)
     * @param name name of the algorithm
     * @param file file this method was called from (optional)
     * @param line line this method was called from (optional)
     * @return T* the pointer to the algorithm
     * @throws runtime_error if pointer not found (DEBUG mode: also if pointer cannot be cast)
     */
    template <typename T, bool FailIfNotFound = false>
    T * algorithm(std::string name, const char* file = __FILE__, int line = __LINE__) {
      using ValueType = std::remove_pointer_t<typename Algorithm_Sequence::value_type>;
      static_assert(std::is_base_of<Algorithm, T>::value, 
        "Template parameter T must derive from Algorithm");

      static_assert(std::is_base_of<Algorithm, ValueType>::value, 
        "Algorithm_Sequence::value_type must derive from Algorithm");
      
      for(Algorithm_Sequence::iterator 
      it = begin(), to = end();
          it != to;
          ++it){
        
        if ((*it)->name == name) {
          T* alg_p = nullptr;
    #ifndef NDEBUG
          alg_p = dynamic_cast<T*>(*it);
          if (!alg_p) {
            std::ostringstream oss;
            oss << "Invalid dynamic_cast to requested algorithm type at "
                << file << ":" << line << "\n"
                << "Failed to cast algorithm named '" << name << "' (internal error)";
            throw std::runtime_error(oss.str());
          }
    #else
          // No check in release mode!
          alg_p = static_cast<T*>(*it);
    #endif
          return alg_p;
        }
      }
      // Not found
      if constexpr (FailIfNotFound) {
        std::ostringstream oss;
        oss << "Algorithm named '" << name << "' not found at "
            << file << ":" << line << " (internal error)";
        throw std::runtime_error(oss.str());
      } else {
        return nullptr;
      }
    }
    
    /**
     * print algorithm sequence.
     * this is a nice debugging function.
     */
    void printSequence();
    
  protected:
    bool clean;

  };
  
} // algorithm

#ifdef NDEBUG
  #define GET_ALGORITHM(seq, T, Strict, name) seq.template algorithm<T,Strict>(name)
#else
  #define GET_ALGORITHM(seq, T, Strict, name) seq.template algorithm<T,Strict>(name, __FILE__, __LINE__)
#endif

#endif