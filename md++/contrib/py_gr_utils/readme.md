## simulation tests
### structure
`ene_ana.py` is a script that reads *tre files using numpy and is used to generate hard-coded values for a short simulation (e.g. 20 steps) that is used by the CI to assert if these values are reproduced upon a given commit.
Additionally, it can be used as a stand-alone program for energy analyses.

`helper_fnc` stores generalized tests (classes) that can be used to create specific simulation tests

`tests_short` store simulation tests that are always performed upon commit

`tests_extensive` is a set of more extensive tests that are done on request only

input files needed for the tests are stored and maintained in a separate repository: gromos_test_files
