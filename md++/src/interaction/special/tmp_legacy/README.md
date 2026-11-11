# Archived geometric restraint interactions

This directory retains the independent distance, angle, dihedral, perturbed and
EDS-distance implementations for reference. They remain in the build sources,
but are no longer registered in `create_special.cc`. Old geometric restraint
activation modes are rejected with migration errors; use COLVARRES/@colvarres.
The legacy EDS distance route is also retired; its replacement uses EDSRESTRAINTS in @colvarres for single-s EDS FORM 1. The archived implementation is used for state-energy/force equivalence checks, not registered for MD execution.

Angle and dihedral mode 3 constraints retain their separate input readers and
constraint algorithms. Distance-field input retains @distrest. Neither relies
on registering these archived restraint interactions.

New geometry and common bias functionality belongs in `../colvar/`. DFUNCT
remains a separate interaction with its own energy accounting. Installed legacy
headers retain the `interaction/special/tmp_legacy/` path. The original relocation
preserved non-include function bodies; shared virtual-atom corrections are
separate behavioral changes.

See `doc/colvar_compact_output.md` for the current output and migration policy.
