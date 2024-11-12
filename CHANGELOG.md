# Changelog


Versions are of the form MAJOR.MINOR.PATCH and this changelog tries to conform
to [Common Changelog](https://common-changelog.org)


## Unreleased

### Fixed

- Fix pipeline for running only per sample analyses. ([#9](https://github.com/metagenlab/zshoman/pull/9)) (Niklaus Johner)

### Changed

- Publish log files to output directory. ([#9](https://github.com/metagenlab/zshoman/pull/9)) (Niklaus Johner)
- Avoid rerunning processes leading to output already present in the output directory. ([#10](https://github.com/metagenlab/zshoman/pull/10)) (Niklaus Johner)

### Added

- Add log file to filter_scaffolds process. ([#9](https://github.com/metagenlab/zshoman/pull/9)) (Niklaus Johner)
- Add script to collect (copy and rename) output files to new location. ([#8](https://github.com/metagenlab/zshoman/pull/8)) (Niklaus Johner)
- Add nf-boost plugin and enable automatic clean-up. ([#9](https://github.com/metagenlab/zshoman/pull/9)) (Niklaus Johner)


## 0.9.0 - 2024-10-31

*Initial release*
