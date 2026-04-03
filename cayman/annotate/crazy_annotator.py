from __future__ import annotations
from dataclasses import dataclass

import os
import logging
import itertools
import collections
import warnings
import hashlib
from typing import IO, List, Iterable, Optional, Iterator, Literal
from pathlib import Path

import polars
import pyhmmer
import numpy as np
from pyhmmer.hmmer import hmmsearch


alphabet = pyhmmer.easel.Alphabet.amino()

logger = logging.getLogger(__name__)


def is_parquet_file(path: Path) -> bool:
    # Check if the file exists
    if not path.is_file():
        raise FileNotFoundError(f"{path} does not exist.")

    with open(path, "rb") as file:
        peek = file.peek(4)

    return peek.startswith(b"PAR1")


class Sequences:
    def __init__(self, sequences: pyhmmer.easel.SequenceBlock):
        self.sequences = sequences

    @ classmethod
    def read_sequences_from_file(self, path, digital = True) -> Sequences:
        """
        Read sequences from file and store them in the attribute sequences
        :param path: path to the sequence file
        :return: None
        """
        with pyhmmer.easel.SequenceFile(
            path,
            digital=digital,
            alphabet=alphabet,
        ) as f:
            return Sequences(f.read_block())


class HMMLoader(Iterable[pyhmmer.plan7.HMM]):
    def __init__(
        self,
        hmms: Iterable[pyhmmer.plan7.HMM],
        seed: int = 42,
        blacklist: Optional[List[str]] = None,
    ):  
        self.seed = seed
        if blacklist is not None:
            # We need to materialize hmms since we iterate twice by necessity
            hmms = list(hmms)
            hmms.sort(key=lambda x: x.name)
            self.hmms = (
                hmm
                for hmm in hmms
                if hmm.name in self.select_hmms(
                    hmms=hmms,
                    blacklist=blacklist,
                    seed=self.seed,
                )
            )
            if not self.hmms:
                raise ValueError(
                    "No hmms passed the selection! Are you sure the cutoff file was correct?"
                )
        else:
            self.hmms = hmms

        if self.seed != 42 and blacklist is None:
            warnings.warn(
                "You changed the seed from the default (42)"
                "but did not pass a blacklist."
                "This will bypass the hmm-fold selection!"
            )

    def __iter__(self) -> Iterator[pyhmmer.plan7.HMM]:
        return iter(self.hmms)

    @staticmethod
    def select_hmms(hmms, blacklist: List[str], seed: int = 42) -> List[str]:
        tbl_data = collections.defaultdict(list)
        for hmm in hmms:
            splits = hmm.name.split(".")[0].split("__")
            tbl_data["hmm"].append(hmm.name)
            tbl_data["family"].append(splits[1])
            tbl_data["fold"].append(int(splits[2][-1]))

        df = (
            polars.from_dict(
                data=tbl_data,
                schema_overrides={
                    "hmm": polars.Utf8,
                    "family": polars.Utf8,
                    "fold": polars.UInt8,
                }
            )
            # For hmms without thresholds but with other family members with trehsolds
            # We need to select one of the hmm with treshold
            # In other words, select hmms not in the blacklist.
            .filter(~polars.col("hmm").is_in(blacklist))
        )
        # we want to group by family and then randomly select one fold from each family
        # to do this I assign a seeded random index to each row
        # and then select the lowest index within each family
        # (group_by maintains order within group)
        return (
            df
            .with_columns(polars.int_range(0, polars.len()).shuffle(seed).alias("rand"))
            .sort("rand")
            .group_by("family")
            .agg(
                polars.first("fold", "hmm"),
            )
            ["hmm"].to_list()
        )

    @staticmethod
    def _read_hmm_from_file(path: Path|str) -> Iterator[pyhmmer.plan7.HMM]:
        """
        Read HMM from file
        :param path: path to the HMM file
        :return: None
        """
        with pyhmmer.plan7.HMMFile(path) as hmm_file:
            for hmm in hmm_file:
                yield hmm

    @staticmethod
    def md5_of_file(path: Path) -> str:
        hasher = hashlib.md5()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                hasher.update(chunk)
        return hasher.hexdigest()

    @classmethod
    def read_hmms(
        self,
        file_with_paths: Optional[Path]=None,
        hmmdb_path: Optional[Path]=None,
        blacklist: Optional[List[str]] = None,
        seed: int = 42,
    ) -> HMMLoader:
        """
        Read HMMs from a either a file with paths to HMM files
        Or from a directory containing HMM files
        Or from a single HMM file
        :param file_with_paths: path to the file with paths
        :param hmmdb_path: path to directory of HMMs or single HMM file
        :return: 'HMMLoader'
        """

        if file_with_paths is None and hmmdb_path is not None:
            if hmmdb_path.is_file():
                return HMMLoader(
                    self._read_hmm_from_file(hmmdb_path),
                    blacklist=blacklist,
                    seed=seed,
                )
            elif hmmdb_path.is_dir():
                return HMMLoader(
                    itertools.chain.from_iterable(
                        self._read_hmm_from_file(os.path.join(hmmdb_path, f))
                        for f in os.listdir(hmmdb_path)
                    ),
                    blacklist=blacklist,
                    seed=seed,
                )
            else:
                raise ValueError(
                    f"Encounterd path which isnt a file or a dir at {str(hmmdb_path)}"
                )

        elif hmmdb_path is None and file_with_paths is not None:
            with open(file_with_paths, 'r') as _in:
                return HMMLoader(
                    itertools.chain.from_iterable(
                        self._read_hmm_from_file(f.strip())
                        for f in _in
                    ),
                    blacklist=blacklist,
                    seed=seed,
                )

        else:
            raise ValueError(
                "Pass either a file with paths or a path to HMMLoader.read_hmms()"
            )

    def write_to_h3m_file(
        self,
        output: Path,
    ):
        """Write the hmms to a single h3m binary file"""
        output = Path(output).with_suffix(f".seed{self.seed}_selected.h3m")
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "wb") as f:
            for hmm in self.hmms:
                hmm.write(f, binary=True)


# TODO in the setup.py download the zenodo repo
# runa static vonverter script wihcih concats and writes them to a binarz h3m file
# and saves them someplace
# the converted hm3 file should live in the annotated folder
# but shouldnt be commited which means add it to the gitignore
# This makes the cli flag hmmdb optional. only if set will the loader take what the user asks

@dataclass(frozen=True)
class CazyHit:
    hmm_name: str
    seq_name: str
    start: int
    end: int
    pvalue: float


class CazyAnnotator:

    def __init__(
        self,
        hmms: Iterable[pyhmmer.plan7.HMM],
    ):
        self.hmms = hmms
        self.background = pyhmmer.plan7.Background(alphabet)

    @staticmethod
    def _annotate(
        hmms: Iterable[pyhmmer.plan7.HMM],
        sequences: pyhmmer.easel.SequenceBlock,
        threads: int = 1,
        blacklist: Optional[List[str]] = None, # List of hmm names to skip
    ) -> Iterator[CazyHit]:

        if blacklist is not None:
            hmms = (
                hmm
                for hmm in hmms
                if hmm.name not in blacklist
            )

        # TODO add a progress bar?
        hmm_hits = hmmsearch(hmms, sequences)
        for hits in hmm_hits:
            for hit in hits.reported:
                for domain in hit.domains.reported:
                    yield CazyHit(
                        hmm_name=hits.query.name,
                        seq_name=hit.name,
                        start=domain.env_from,
                        end=domain.env_to,
                        pvalue=domain.pvalue,
                        # domain.i_evalue, # Do not collect unless you set Z
                        # domain.c_evalue, # Do not collect unless you set DomZ
                    )
        # hmm names are only needed to match to the thresholds
        # so once we integrate the thresholds, we dont need the names anymore

    def annotate(
            self,
            sequences: pyhmmer.easel.SequenceBlock,
            threads: int = 1,
            blacklist: Optional[List[str]] = None, # List of hmm names to skip
        ) -> CazyResultsTable:
        return CazyResultsTable.create_from_TopHits(
            self._annotate(
                hmms=self.hmms,
                sequences=sequences,
                threads=threads,
                blacklist=blacklist,
            ),
        )

    def curate_annotations(
        self,
        thresholds: ThresholdTable,
        cazy_results: CazyResultsTable,
    ) -> polars.DataFrame:

        annotations_filtered = (
            cazy_results
            .apply_thresholds(thresholds)
            .disentangle_domains()
            .table
            .rename(
                {"start": "domain_start", "end": "domain_end"}
            )
            # Transform to nucleotide coordinates
            .with_columns(
                start = polars.col("domain_start")*3 - 2,
                end = polars.col("domain_end")*3 - 2,
                annotLength = (polars.col("domain_end")-polars.col("domain_start"))*3
            )
            # bring the columns in order
            .select(
                [
                    'sequenceID',
                    'start',
                    'end',
                    'pvalue',
                    'family',
                    'annotLength',
                    'domain_start',
                    'domain_end'
                ]
            )
            .sort(by=["sequenceID", "domain_start"])
        )
        return annotations_filtered


class CazyResultsTable:
    schema = polars.Schema(
        {
            "hmm": polars.String,
            "sequenceID": polars.String,
            "start": polars.UInt32,
            "end": polars.UInt32,
            "pvalue": polars.Float64,
            "family": polars.String,
            "familyType": polars.String,
            "fold": polars.UInt8,
        }
    )

    def __init__(self, table: polars.DataFrame):
        self.table = table.select(
            [
                "hmm",
                "sequenceID",
                "start",
                "end",
                "pvalue",
                "family",
                "familyType",
                "fold",
            ]
        )
        self._validate()

    def _validate(self):
        if isinstance(self.table, polars.DataFrame):
            if self.table.schema != self.schema:
                raise ValueError(f"Got wrong schema. Expected {self.schema}")
        else:
            raise TypeError(f"Expected polars DataFrame but got {type(self.table)}")

    @classmethod
    def load(cls, file: str | Path) -> CazyResultsTable:
        if is_parquet_file(Path(file)):
            tbl = polars.read_parquet(file)
        else:
            tbl = polars.read_csv(file, separator="\t", schema_overrides=cls.schema)

        return CazyResultsTable(tbl)

    def dump(
        self,
        file: str | Path | IO[bytes],
        format: Literal["tsv", "parquet", "csv"] = "parquet",
    ):
        if format == "parquet":
            self.table.write_parquet(file)
        elif format == "tsv":
            self.table.write_csv(file, separator="\t")
        elif format == "csv":
            self.table.write_csv(file)
        else:
            raise NotImplementedError

    @classmethod
    def create_from_TopHits(
        cls,
        hits: Iterable[CazyHit],
    ) -> CazyResultsTable:
        tbl_data = collections.defaultdict(list)
        for hit in hits:
            tbl_data["hmm"].append(hit.hmm_name)
            tbl_data["sequenceID"].append(hit.seq_name)
            tbl_data["start"].append(hit.start)
            tbl_data["end"].append(hit.end)
            tbl_data["pvalue"].append(hit.pvalue)

        if not tbl_data:
            raise ValueError("Did not get any hits!")

        df = (
            polars.from_dict(
                data=tbl_data,
                schema_overrides=cls.schema
            )
            # Process hmm name into family, familyType and fold columns
            .with_columns(
                polars.col("hmm").str.extract(r"^(.*?)__", 1).alias("family"),
                polars.col("hmm").str
                .extract(r"fold(\d+)", 1) # capture digits after "fold"
                .cast(polars.UInt8)
                .alias("fold")
            )
            .with_columns(
                polars.col("family").str.replace_all(r"[_\d]", "").alias("familyType")
            )
        )
        return CazyResultsTable(df)

    def apply_thresholds(self, thresholds: ThresholdTable) -> CazyResultsTable:
        self.table = (
            self.table
            .join(
                other=thresholds.table,
                on=["family", "fold", "familyType"],
                how="left",
            )
            # NOTE i think this code is necessary
            # but unfortunately my current tests dont actually show the difference
            # because for all those where we apply the median treshold
            # they fail the median threshold which is the same as failing agaist null
            .join(
                thresholds.familyType_median_cutoffs,
                on="familyType",
                how="left",
            )
            # NOTE
            # where the cutoff is none, replace it with the median cutoff
            # of the enitre Cazyme familyType (e.g. CBM or GH)
            .with_columns(
                polars.coalesce(
                    polars.col("cutoff"),
                    polars.col("median_cutoff")
                ).alias("cutoff")
            )
            # filter by fold and family specific pvalue treshold
            .filter(polars.col("pvalue")<polars.col("cutoff"))
            .select(self.schema.names())
        )
        return self

    @staticmethod
    def _disentangle_group(
        gene_group_df: polars.DataFrame,
        start_col_name: str,
        end_col_name: str,
        threshold_column: str,
    ) -> polars.DataFrame:

        df = gene_group_df.sort(threshold_column)
        starts = df[start_col_name].to_numpy()
        ends = df[end_col_name].to_numpy()
        n = len(starts)
        keep = np.ones(n, dtype=bool)

        for i in range(n):
            if not keep[i]:
                continue
            # mask out everything overlapping domain i
            keep[i+1:] &= (ends[i+1:] < starts[i]) | (starts[i+1:] > ends[i])

        return df.filter(keep)

    def disentangle_domains(
        self,
    ) -> CazyResultsTable:

        def wrapper(group_df: polars.DataFrame) -> polars.DataFrame:
            return self._disentangle_group(
                group_df,
                start_col_name="start",
                end_col_name="end",
                threshold_column="pvalue",
            )

        self.table = self.table.group_by("sequenceID").map_groups(wrapper)
        return self


class ThresholdTable:
    schema = polars.Schema(
        {
            # "Unnamed: 0": polars.UInt32,
            "cutoff": polars.Float64,
            # "F1": polars.Float64,
            # "predictedPos": polars.UInt32,
            # "Ps": polars.UInt32,
            "family": polars.String,
            "fold": polars.UInt8,
        }
    )
    # NOTE that the cutoffs v3 table includes useless columns
    # that we never ever use anyway. I just skip them to keep the code clean

    def __init__(self, table: polars.DataFrame):
        self.table = table.select(
            [
                # "Unnamed: 0",
                "cutoff",
                # "F1",
                # "predictedPos",
                # "Ps",
                "family",
                "fold",
            ]
        )
        self._validate()
        self._preprocess() # TODO this is all unnecessary and can be saved in the tbl

    def _validate(self):
        if isinstance(self.table, polars.DataFrame):
            if self.table.schema != self.schema:
                raise ValueError(f"Got wrong schema. Expected {self.schema}")
        else:
            raise TypeError(f"Expected polars DataFrame but got {type(self.table)}")

    @classmethod
    def load(cls, file: str | Path) -> ThresholdTable:
        if is_parquet_file(Path(file)):
            tbl = polars.read_parquet(file)
        else:
            tbl = polars.read_csv(file, separator=",", schema_overrides=cls.schema)

        return ThresholdTable(tbl)

    @property
    def hmms_which_will_be_skipped(self) -> List[str]:
        """
        If only some folds of a CAZy family have thresholds,
        hits to those HMMs without thresholds are discarded.

        This function identifies the HMMs missing from the treshold table (which
        therefore lack a threshold). Computationally, it is most efficient to simply
        never annotate with these hmms in the first place

        Returns list of str of hmm names
        """
        # make a df of all the families with cutoffs
        # and all the folds that should be there
        expected = (
            self.table
            .select(["family"])
            .unique()
            .with_columns(
                fold=polars.lit([1, 2, 3, 4, 5])
            )
            .explode("fold")
        )

        # do an anti join of the expected family-fold cobminations with existing ones
        # to get a df of all the misssing family-fold combinations
        # for which we lack threshodls
        missing = (
            expected
            .join(
                other=self.table.select(["family", "fold"]),
                on=["family", "fold"],
                how="anti",
            )
            .with_columns(
                hmm_name=polars.col("family")
                + "__" + polars.col("family") +
                "__fold" + polars.col("fold").cast(polars.String) + ".mafft"
            )
            ["hmm_name"].to_list()
        )
        return missing

    @property
    def familyType_median_cutoffs(self) -> polars.DataFrame:
        return (
            self.table
            .group_by("familyType")
            .agg(
                polars.col("cutoff").median().alias("median_cutoff"),
            )
        )

    def _preprocess(self):
        counts = (
            self.table
            .group_by("family")
            .agg(
                polars.col("fold").len().alias("num_folds"),
            )
        )
        self.table = (
            self.table
            .with_columns(
                polars.col("family")
                # we keep only the letters
                .str.replace_all(r"[_\d]", "") # remove underscores and digits
                .alias("familyType")
            )
            .join(
                other=counts,
                how="left",
                on="family",
            )
        )