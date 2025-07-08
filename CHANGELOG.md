# Changelog


Versions are of the form MAJOR.MINOR.PATCH and this changelog tries to conform
to [Common Changelog](https://common-changelog.org)


## Unreleased

### Fixed

- Fix pipeline for running only per sample analyses. ([#9](https://github.com/metagenlab/zshoman/pull/9)) (Niklaus Johner)
- Fix pipeline for samples lacking any eukaryotic contigs. ([#18](https://github.com/metagenlab/zshoman/pull/18)) (Niklaus Johner)

### Changed

- Avoid waiting for all samples before starting gene counts and functional annotations. ([#23](https://github.com/metagenlab/zshoman/pull/23)) (Niklaus Johner)
- Enable nf-boost cleanup and make the default error strategy to ignore. ([#21](https://github.com/metagenlab/zshoman/pull/21)) (Niklaus Johner)
- Publish log files to output directory. ([#9](https://github.com/metagenlab/zshoman/pull/9)) (Niklaus Johner)
- Avoid rerunning processes leading to output already present in the output directory. ([#10](https://github.com/metagenlab/zshoman/pull/10), [#13](https://github.com/metagenlab/zshoman/pull/13)) (Niklaus Johner)
- Optimize resources requested for METAEUK process. ([#11](https://github.com/metagenlab/zshoman/pull/11)) (Niklaus Johner)
- Retry spades on error code 250. ([#17](https://github.com/metagenlab/zshoman/pull/17)) (Niklaus Johner)

### Added

- Add script to copy input files from a source to the input folder. ([#30](https://github.com/metagenlab/zshoman/pull/30)) (Niklaus Johner)
- Add support for 5 lanes or experiments for a given sample. ([#29](https://github.com/metagenlab/zshoman/pull/29)) (Niklaus Johner)
- Add support for samples with multiple single-end experiments. ([#28](https://github.com/metagenlab/zshoman/pull/28)) (Niklaus Johner)
- Add post-processing script calculate abundances of complete KEGG modules. ([#24](https://github.com/metagenlab/zshoman/pull/24)) (Niklaus Johner)
- Add script to merge mOTUS outputs into a single table. ([#20](https://github.com/metagenlab/zshoman/pull/20)) (Niklaus Johner)
- Add script to filter out already processed samples from input. ([#19](https://github.com/metagenlab/zshoman/pull/19)) (Niklaus Johner)
- Add support for per sample runs in annotations post-processing script. ([#24](https://github.com/metagenlab/zshoman/pull/24)) (Niklaus Johner)
- Add option *resume_from_output* to avoid using output directory to skip processes. ([#16](https://github.com/metagenlab/zshoman/pull/16)) (Niklaus Johner)
- Add support for multi-lane samples. ([#15](https://github.com/metagenlab/zshoman/pull/15)) (Niklaus Johner)
- Add various scripts for cleaning up the work and output directories. ([#14](https://github.com/metagenlab/zshoman/pull/14), [#19](https://github.com/metagenlab/zshoman/pull/19), [#22](https://github.com/metagenlab/zshoman/pull/22)) (Niklaus Johner)
- Add post-processing script to check the quality of the results. ([#11](https://github.com/metagenlab/zshoman/pull/11)) (Niklaus Johner)
- Add log file to filter_scaffolds process. ([#9](https://github.com/metagenlab/zshoman/pull/9)) (Niklaus Johner)
- Add script to collect (copy and rename) output files to new location. ([#8](https://github.com/metagenlab/zshoman/pull/8)) (Niklaus Johner)
- Add nf-boost plugin. ([#9](https://github.com/metagenlab/zshoman/pull/9)) (Niklaus Johner)


## 0.9.0 - 2024-10-31

*Initial release*
