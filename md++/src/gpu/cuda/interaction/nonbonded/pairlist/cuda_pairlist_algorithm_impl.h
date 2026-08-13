#pragma once

namespace interaction {
  /**
   * @class CUDA_Pairlist_Algorithm_Impl
   * GPU-native implementation backing CUDA_Pairlist_Algorithm.
   *
   * Not a Backend-templated pImpl: CUDA_Pairlist_Algorithm is only ever
   * compiled/included under USE_CUDA (its single call site,
   * create_nonbonded.cc, guards the #include the same way), so there is
   * no CPU counterpart to keep in sync with. Per PLAN.md D5, there is
   * exactly one GPU-native pairlist and it is never dispatched through
   * make_algorithm's Backend-template mechanism -- the earlier
   * Backend-templated shape bought no genericity (the <gpuBackend>
   * specialization shared nothing with the primary template body) and
   * was dropped for this plain class + the standard .h/.cc/.cu split.
   */
  class CUDA_Pairlist_Algorithm_Impl {
    public:
      CUDA_Pairlist_Algorithm_Impl();

      int init(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim,
        std::ostream &os = std::cout,
        bool quiet = false);

      void set_cutoff(double const cutoff_short, double const cutoff_long);

      void prepare_cog(configuration::Configuration & conf,
                  topology::Topology & topo,
                  simulation::Simulation & sim);

      /**
       * put the chargegroups into the box
       */
      template<math::boundary_enum b>
      void _prepare_cog(configuration::Configuration & conf,
                        topology::Topology & topo);

      /**
       * order atoms / chargegroups based on their grid cell?
       */
      void reorder(configuration::Configuration & conf,
                  topology::Topology & topo,
                  simulation::Simulation & sim);

    protected:
      /**
       * squared shortrange cutoff.
       */
      double m_cutoff_short_2;
      /**
       * squared longrange cutoff.
       */
      double m_cutoff_long_2;
      /**
       * longrange cutoff.
       */
      double m_cutoff_long;
      /**
       * shortrange cutoff.
       */
      double m_cutoff_short;

    private:
      /**
       * chargegroup center of geometry array.
       */
      math::CuVArray m_cg_cog;
      /**
       * chargegroup cell indices array.
       */
      gpu::cuvector<ushort4> m_cg_cells;
      /**
       * flat copy of m_cg_cells[i].w (Morton cell index), one entry per
       * chargegroup -- feeds a Thrust sort-by-key in the block-build step
       * (TILE_PAIRLIST_DESIGN.md §3 step 2), which needs a plain array,
       * not one field of a ushort4.
       */
      gpu::cuvector<unsigned> m_cg_sort_key;
  };
}
