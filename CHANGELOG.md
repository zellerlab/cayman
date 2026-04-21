# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [v0.11.1] - 2026-04-21
This update allows the user to supply an alignment file as input instead of reads directly
[v0.11.1]: https://github.com/zellerlab/cayman/compare/v0.11.0...v0.11.1

### Added
- `--sam` optional input flag added to pass SAM alignment file
- `--bam` optional input flag added to pass BAM alignment file
- `--import_readcounts` <int> and `-unmarked_orphans` flags added to correct counting from alignment file
- `--keep_alignment_file` <path> added which will save the SAM alignment file in case reads were passed as input

### Breaking
- `--cpus_for_alignment` flag got renamed `--threads` in the profile subparser
- `bwa_index` is no longer a positional argument but now is only required with read input. pass it after the `--bwa_index` flag


## [v0.11.0] - 2026-03-20
This is a purely performance driven update. IO remains exactly the same - just 4-5x faster

### Added
- API entry points to main function + handle exit codes
- API entry points to cayman annotate
- Unit-tests for CLI of cayman annotate
- Ability to read HMMs from either a single file (e.g. h3m), dir or list of file paths
- Added HMMLoader, CazyHit, ThresholdTable, Sequences and CazyResultsTable classes
- Let polars handle column wise dataframe operations (ultimately entrirely replacing pandas)

### Fixed
- Fix compatibility with pyHMMER [v0.12](https://github.com/althonos/pyhmmer/releases/tag/v0.12.0)
- Refactored cayman annotate functions
- Multiprocessing now uses threads (pyHMMER default) instead of processes
- Sequences now load as binary [DigitalSequenceBlocks](https://pyhmmer.readthedocs.io/en/stable/api/easel/block.html)
- HMM loading happens in parallel to scanning with lazy HMMLoader class
- HMM hits which are discarded due to lacking cutoffs are now never scanned

## [v0.10.2] - 2025-12-05
[v0.10.2]: https://github.com/zellerlab/cayman/compare/v0.10.1...v0.10.2

### Added
- Protein start and end coordinates added to annotation output

### Fixed
- Clarified README
- improved descriptions for CLI options in README
- added conda install to README

## [v0.10.1] - 2024-12-17
[v0.10.1]: https://github.com/zellerlab/cayman/compare/v0.9.7...v0.10.1

### Added
- Cayman annotator for annotating protein fasta files with CAZymes
- Argument parser now has subparsers for either protein annotation or read mapping

### Fixed
- Container definitions updated
- build is now in gitignore

## [v0.9.7] - 2024-08-07
[v0.9.7]: 

Initial release.

