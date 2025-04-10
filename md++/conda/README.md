## To create the conda packages:
1. Have a working conda, with conda-build installed
2. Set correct version in VERSION and conda/meta.yaml
3. `cd` to md++ root
4. `conda-build .`
5. `anaconda upload -u OUR_CHANNEL_NAME PATH_TO_MD_BUILD.conda`