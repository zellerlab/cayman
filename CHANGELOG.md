# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [v0.12.0] - 2026-03-21
This update fundamentally changes how proteins are annotated.
Instead of annotating with all hmm-folds for a cazyme family,
only one fold is randomly (seeded) selected for annotation.
This removes the need for per-residue merging of majority-vote annotations.
[v0.12.0]: https://github.com/zellerlab/cayman/compare/v0.11.0...v0.12.0

### Breaking
- Start and end coordinates will now direcly reflect hmm boundaries
- p-values now actually reflect the p-values calculated by HMMER.
- overlaps between annotation hits are now resolved by selecting the more significat hit
- API for `CazyResultsTable` methods `apply_thresholds` and `disentangle_domains` operate on self and return an instance of self

### Added
- There is now an optional `--seed` argument in the cli for `annotate_proteome`
- Changing the seed will now change the results since this will cause different hmms per fold to be selected!
- API for HMMLoader now has optional `blacklist` parameter
- Tests for HMMLoader

### Fixed
- Oddly split annotations like (from unittests):
    AVQ15184.1 GH13 120->144; GH13_11 145->493; GH13 494->508
    Now become:
    AVQ15184.1 GH13_11, 145->496
    This is because overlaps are now treated differently:
        1. There are no more overlaps between folds which need to be reconciled
        2. Overlaps between hmm families are resolved by taking the lowest p-value hit
    Previously results were aggregated per residue,
    creating un-biological start-end coordinates as shown here.
- No more doube geometric averaging of p-values over folds per residue and over residues in a sequence
- Much faster annotation runtime -> Selecting 1 hmm per fold slashes runtime by factor 5
- Some annotations which were being surpressed due to p-value avering (I think) are now being found


## [v0.11.0] - 2026-03-20
This is a purely performance driven update. IO remains exactly the same - just 4-5x faster
[v0.11.0]: https://github.com/zellerlab/cayman/compare/v0.10.1...v0.11.0

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

