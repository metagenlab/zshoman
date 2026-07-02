# Changelog


Versions are of the form MAJOR.MINOR.PATCH and this changelog tries to conform
to [Common Changelog](https://common-changelog.org)

## Unreleased

### Fixed

- Fix loading genes from output folder. ([#55](https://github.com/metagenlab/zshoman/pull/55)) (Niklaus Johner)

### Changed

- Improve pipeline disk space usage. ([#55](https://github.com/metagenlab/zshoman/pull/55)) (Niklaus Johner)
- Disable read error correction on retry for SPADES. ([#55](https://github.com/metagenlab/zshoman/pull/55)) (Niklaus Johner)
- Speed up gene catalog indexing by increasing batch size and requested resources. ([#55](https://github.com/metagenlab/zshoman/pull/55)) (Niklaus Johner)
- Add more checkpoints to avoid running unnecessary processes when restarting from output. ([#55](https://github.com/metagenlab/zshoman/pull/55)) (Niklaus Johner)
- Parallelize script to copy downloaded files. ([#50](https://github.com/metagenlab/zshoman/pull/50)) (Niklaus Johner)
- Do not skip read error correction during assembly. ([#51](https://github.com/metagenlab/zshoman/pull/51)) (Niklaus Johner)
- Use metaspades instead of spades for assembly. ([#49](https://github.com/metagenlab/zshoman/pull/49)) (Niklaus Johner)
- Support runs using the gene catalog in clean-up and pre- and post-processing scripts. ([#47](https://github.com/metagenlab/zshoman/pull/47)) (Niklaus Johner)
- Adapt clean-up and pre- and post-processing scripts to new samplesheets. ([#47](https://github.com/metagenlab/zshoman/pull/47)) (Niklaus Johner)
- Change name of output files from bbmerge to simplify downstream usage. ([#45](https://github.com/metagenlab/zshoman/pull/45)) (Niklaus Johner)
- Avoid I/O bottleneck when restarting from output to calculate the gene catalog. ([#43](https://github.com/metagenlab/zshoman/pull/43)) (Niklaus Johner)
- Avoid redoing contig classification and gene calling when already done. ([#41](https://github.com/metagenlab/zshoman/pull/41)) (Niklaus Johner)
- Support samples with more than 5 lanes/experiments. ([#40](https://github.com/metagenlab/zshoman/pull/40)) (Farid Chaabane, Niklaus Johner)
- Improve gene profiling by using coverage instead of number of reads. ([#38](https://github.com/metagenlab/zshoman/pull/38)) (Niklaus Johner)
- Use MMSEQS2 LINCLUST instead of CDHIT to create the gene catalog. ([#37](https://github.com/metagenlab/zshoman/pull/37)) (Niklaus Johner)
- Set threshold for per sample gene profiling to 99% identity instead of 95% ([#37](https://github.com/metagenlab/zshoman/pull/37)) (Niklaus Johner)

### Added

- Add MAG reconstruction workflow. ([#54](https://github.com/metagenlab/zshoman/pull/54)) (Sedreh Nassirnia)
- Add documentation of the output and its analysis. ([#44](https://github.com/metagenlab/zshoman/pull/44)) (Niklaus Johner)
- Add support for phanta tables in post-processing script to merge outputs. ([#41](https://github.com/metagenlab/zshoman/pull/41)) (Niklaus Johner)

## 1.0.0 - 2025-12-12

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
