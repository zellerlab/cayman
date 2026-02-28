from __future__ import annotations

import os
import logging
import itertools
from typing import List, Iterable, Optional, Iterator, Tuple
from pathlib import Path

import pandas as pd
import pyhmmer
from pyhmmer.hmmer import hmmsearch

# from tqdm import tqdm


alphabet = pyhmmer.easel.Alphabet.amino()

logger = logging.getLogger(__name__)

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


class HMM_Loader:
    def __init__(self, hmms: Iterable[pyhmmer.plan7.HMM]):
        self.hmms = hmms

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

    @classmethod
    def read_hmms(
        self,
        file_with_paths: Optional[Path]=None,
        hmmdb_path: Optional[Path]=None,
    ) -> HMM_Loader:
        """
        Read HMMs from a either a file with paths to HMM files
        Or from a directory containing HMM files
        Or from a single HMM file
        :param file_with_paths: path to the file with paths
        :param hmmdb_path: path to directory of HMMs or single HMM file
        :return: 'HMM_Loader'
        """

        if file_with_paths is None and hmmdb_path is not None:
            if hmmdb_path.is_file():
                return HMM_Loader(self._read_hmm_from_file(hmmdb_path))
            elif hmmdb_path.is_dir():
                return HMM_Loader(
                    itertools.chain.from_iterable(
                        self._read_hmm_from_file(os.path.join(hmmdb_path, f))
                        for f in os.listdir(hmmdb_path)
                    )
                )
            else:
                ValueError(
                    f"Encounterd path which isnt a file or a dir at {str(hmmdb_path)}"
                )

        elif hmmdb_path is None and file_with_paths is not None:
            with open(file_with_paths, 'r') as _in:
                return HMM_Loader(
                    itertools.chain.from_iterable(
                        self._read_hmm_from_file(f.strip())
                        for f in _in
                    )
                )

        else:
            raise ValueError(
                "Pass either a file with paths or a path to HMM_Loader.read_hmms()"
            )


# TODO in the setup.py download the zenodo repo
# runa static vonverter script wihcih concats and writes them to a binarz h3m file
# and saves them someplace
# the converted hm3 file should live in the annotated folder
# but shouldnt be commited which means add it to the gitignore
# This makes the cli flag hmmdb optional. only if set will the loader take what the user asks

class CazyAnnotator:

    def __init__(
        self,
        hmms: Iterable[pyhmmer.plan7.HMM],
    ):
        self.hmms = hmms
        self.annotations_by_family_and_fold = []
        self.annotations_by_family_and_fold_filtered =  []
        self.annotations_filtered = None
        self.background = pyhmmer.plan7.Background(alphabet)

    @staticmethod
    def _annotate(
        hmms: Iterable[pyhmmer.plan7.HMM],
        sequences: pyhmmer.easel.SequenceBlock,
        threads: int = 1,
    ) -> Tuple[List[pd.DataFrame], List[str]]:

        hmm_hits = hmmsearch(hmms, sequences)
        dfs = []
        hmm_names = []
        for hits in hmm_hits:
            hmm_names.append(hits.query.name)
            res: List[Tuple[str, int, int, int|float]] = []
            for hit in hits.reported:
                name = hit.name
                for domain in hit.domains.reported:
                    res.append(
                        (
                            name,
                            domain.env_from,
                            domain.env_to,
                            domain.pvalue,
                            # domain.i_evalue, # Do not collect unless you set Z
                            # domain.c_evalue, # Do not collect unless you set DomZ
                        )
                    )
            if len(res)==0:
                df = pd.DataFrame(
                    data=None,
                    columns=["moduleID", "start", "end", "pvalue"],  # ty:ignore[invalid-argument-type]
                )
            # TODO why are we calling drop_duplicates? how couuld that happen unless
            # we have duplicate sequences
            df = pd.DataFrame(
                res,
                columns=["moduleID", "start", "end", "pvalue"],  # ty:ignore[invalid-argument-type]
            ).drop_duplicates()

            dfs.append(df)

        # TODO currently i return a list of dfs and a list of hmm names
        # better to have it return one big df with an extra column for hmm names
        # hmm names are only needed to match to the thresholds
        # so once we integrate the thresholds, we dont need the names anymore

        # TODO Also this returned df should then be its own class on which the fns:
        # curate_annotations, resolve_overlapping_annotations and merge_annots operate

        return dfs, hmm_names

    def annotate(
            self,
            sequences: pyhmmer.easel.SequenceBlock,
            threads: int = 1,
        ):
        self.annotations_by_family_and_fold = self._annotate(
            hmms=self.hmms,
            sequences=sequences,
            threads=threads,
        )

    def curate_annotations(self, precomputed_hmm_cutoffs):
        # TODO: Clean this code up..
        from re import sub
        
        median_cutoffs = pd.read_csv(precomputed_hmm_cutoffs)
        median_cutoffs['familyType'] = [sub(r"\d", "", x.replace("_", "")) for x in median_cutoffs['family']]
        median_cutoffs = median_cutoffs.groupby('familyType').median(numeric_only=True)['cutoff']
        cutoffs_all = pd.read_csv(precomputed_hmm_cutoffs)
        family_fold_results_filtered = []
        for index, (family_fold_results, hmm_name) in enumerate(
            zip(*self.annotations_by_family_and_fold, strict=True)
        ):
            family = hmm_name.split("__")[0]
            familyType = sub(r"\d", "", family).replace('_', '')
            fold = hmm_name.split("fold")[1].split('.mafft')[0]
            #print(family, fold)
            cutoffs_this = cutoffs_all[cutoffs_all['family'] == family]
            cutoffs_this = cutoffs_this[cutoffs_this['fold'] == int(fold)]       

            # Exclude empty annotations...
            if family_fold_results.shape[0] == 0:
                continue     

            if cutoffs_this.shape[0] != 1:
                #print("Cant find optimal cutoff for {}, presumably because CV didn't happen, presumably because small family.".format(hmm_name))
                if any([True if x == family else False for x in cutoffs_all['family']]):
                    # If we have other folds of this family, it's okay
                    #print("Ignoring this fold")
                    continue
                else:
                    # if not, simply take median p-value cutoff
                    cutoff = median_cutoffs[familyType]
                    #print("Taking median p-value cutoff {}".format(cutoff))
            else:
                cutoff = list(cutoffs_this['cutoff'])[0]

            #print(GMGC_annots.head())
            family_fold_results_filtered = family_fold_results[family_fold_results['pvalue'] < cutoff].copy(deep=True)
            family_fold_results_filtered['family'] = family
            family_fold_results_filtered['fold'] = fold
            self.annotations_by_family_and_fold_filtered.append(family_fold_results_filtered)
        self.annotations_filtered = pd.concat(self.annotations_by_family_and_fold_filtered)
        del self.annotations_by_family_and_fold
        del self.annotations_by_family_and_fold_filtered
        #self.annotations_filtered = self.annotations_filtered.iloc[:, 1:]
        #GMGC_annots_all.drop('pvalue', axis = 1, inplace = True)

        self.annotations_filtered.rename(columns={'moduleID': 'sequenceID'}, inplace=True)
        
        fold_counts = cutoffs_all.groupby('family').count()['fold'].to_frame().rename({'fold' : 'num_folds'}, axis = 'columns')
        fold_counts['num_folds'] = [int(x)for x in fold_counts['num_folds']]
        annotations_with_fold_counts = pd.merge(self.annotations_filtered, fold_counts, on = 'family', how = 'left')
        annotations_with_fold_counts = annotations_with_fold_counts.groupby(['sequenceID', "family"])

        annotations_with_fold_counts_series = []
        names = []
        for name, group in annotations_with_fold_counts:
            names.append(name)
            annotations_with_fold_counts_series.append(group)

        annotations_with_fold_counts_series = pd.Series(annotations_with_fold_counts_series)
        
        logger.info("Merging annotations...")
        self.annotations_filtered = pd.concat(list(annotations_with_fold_counts_series.apply(CazyAnnotator.merge_annots)))

        tmp2 = self.annotations_filtered.groupby("sequenceID")
        aSeries = []
        names = []
        for name, group in tmp2:
            names.append(name)
            aSeries.append(group)
        aSeries = pd.Series(aSeries)
        logger.info("Resolving overlapping annotations...")
        self.annotations_filtered = pd.concat(list(aSeries.apply(CazyAnnotator.resolve_overlapping_annotations)))

        # Write start and end coordinates out as integers and not floats
        self.annotations_filtered['start'] = [int(x) for x in self.annotations_filtered['start']]
        self.annotations_filtered['end'] = [int(x) for x in self.annotations_filtered['end']]
        self.annotations_filtered['annotLength'] = self.annotations_filtered['end'] - self.annotations_filtered['start']
        self.annotations_filtered = self.annotations_filtered[[True if x >= 10 else False for x in list(self.annotations_filtered['annotLength'])]]

        self.annotations_filtered['start_protein'] = self.annotations_filtered['start']
        self.annotations_filtered['end_protein'] = self.annotations_filtered['end']

        ## TRANSFORM INTO NUCLEOTIDE COORDINATES!
        self.annotations_filtered['start'] = [(x * 3) - 2 for x in self.annotations_filtered['start']]
        self.annotations_filtered['end'] = [x * 3 for x in self.annotations_filtered['end']]
        self.annotations_filtered['annotLength'] = [x * 3 for x in self.annotations_filtered['annotLength']]
        #self.annotations_filtered.to_csv(base_dir + "/" + gc_name + "_all_v3_FINAL.csv", index = False)        

    @staticmethod
    def resolve_overlapping_annotations(df):
        from collections import defaultdict  
        if df.shape[0] == 1:
            return(df)

        if all(df['pvalue'] == 0):
            df['pvalue'] = 1E-500
        else:
            smallestNonZero = min([x for x in df['pvalue'] if x > 0])
            df['pvalue'] = [x if x > 0 else smallestNonZero for x in df['pvalue']]

        #numberFolds = len(list(set(list(df['fold']))))
        df['start'] = [int(x) for x in df['start']]
        df['end'] = [int(x) for x in df['end']]
        ranges = [list(range(x,y)) for x, y in zip(df['start'], df['end'])]
        cc = defaultdict(int)
        # Get geomtric mean p-value over folds within a residue
        p_vals = list(df['pvalue'])
        families = list(df['family'])
        p_vals_dict = defaultdict(dict)
        for i_r in range(len(ranges)):
            for i in ranges[i_r]:
                ## cc is like a depth of coverage dict.
                cc[i] += 1
                p_vals_dict[i][i_r] = p_vals[i_r]
        #print(p_vals_dict)

        ranges = ranges
        final_ranges = dict()
        # loop over residues, memorizing the 'range to take' by taking that range that correspond to the smallest p-value
        # Save it in final_ranges
        for residue in sorted(cc):
            #coverage = cc[residue]
            p_values = p_vals_dict[residue]
            #print(p_values)
            smallest_p_value = None
            range_i_to_take = None
            #print("I have {} p-values".format(len(p_values)))
            for range_index, p_value in p_values.items():
                #print(range_index, p_value)
                if smallest_p_value == None or p_value < smallest_p_value:
                    #print("Updating smallest_p_value")
                    smallest_p_value = p_value
                    range_i_to_take = range_index
            final_ranges[residue] = range_i_to_take

        # Loop over chosen ranges
        oldRangeIndex = None
        rangesOut = []
        currentRange = []
        for residue in sorted(cc):
            currentRangeIndex = final_ranges[residue]
            #print(oldRangeIndex, currentRangeIndex)
            if oldRangeIndex == None or (oldRangeIndex) == currentRangeIndex:
                #print("Appending to current range")
                currentRange.append([residue, currentRangeIndex])
            else:
                if len(currentRange) > 1:
                    rangesOut.append(currentRange)
                currentRange = []
            oldRangeIndex = currentRangeIndex

        # Append final range...
        if len(currentRange) > 1:
            rangesOut.append(currentRange)

        ranges = []
        # There's a one-off problem with the index, fixing it here
        for i, r in enumerate(rangesOut):
            if i >= 1:
                ranges.append([r[0][0] -1, r[-1][0]])
            else:
                ranges.append([r[0][0], r[-1][0]])
        #ranges = [[x[0][0], x[-1][0]] for x in rangesOut]
        pvals_fin = [p_vals[x[0][1]] for x in rangesOut]
        families = [families[x[0][1]] for x in rangesOut]

        out = pd.DataFrame({"sequenceID" : [list(df['sequenceID'])[0]] * len(ranges),
                        "start" : [x[0] for x in ranges],
                        "end" : [x[1] for x in ranges],
                        'pvalue' : pvals_fin,
                        'family' : families})
        return(out)

    @staticmethod
    def merge_annots(df):
        from collections import defaultdict
        #from scipy.stats.mstats import gmean as gmean_old # This import triggers the subnormals warning - let's not get into the details 
        from statistics import geometric_mean
        import numpy as np
    
        # Emulate scipy.stats.mstats.gmean
        def gmean(iterable):
            if all(x == 0 for x in iterable):
                return 0
            else:
                return(geometric_mean(iterable))

        numFolds = list(df['num_folds'])[0]
        if (np.isnan(numFolds)):
            numFolds = 1
        foldCutoff = round(numFolds/2)

        if all(df['pvalue'] == 0):
            df['pvalue'] = 1E-500
        else:
            smallestNonZero = min([x for x in df['pvalue'] if x > 0])
            df['pvalue'] = [x if x > 0 else smallestNonZero for x in df['pvalue']]
        #numberFolds = len(list(set(list(df['fold']))))
        ranges = [list(range(x,y)) for x, y in zip(df['start'], df['end'])]
        cc = defaultdict(int)
        # Get geomtric mean p-value over folds within a residue
        p_vals = list(df['pvalue'])
        p_vals_dict = defaultdict(list)
        for i_r in range(len(ranges)):
            for i in ranges[i_r]:
                cc[i] += 1
                p_vals_dict[i].append(p_vals[i_r])
        for index, entry in p_vals_dict.items():
            p_vals_dict[index] = gmean(entry)
        # Merge ranges
        ranges = []
        currentRange = []
        oldKey = None
        for key in sorted(cc):
            value = cc[key]
            if value > foldCutoff:
                if oldKey == None or (oldKey+1) == key :
                    currentRange.append(key)
                else:
                    #print('appending range')
                    if len(currentRange) > 1:
                        ranges.append(currentRange)
                    currentRange = []
                oldKey = key
        if len(currentRange) > 1:
            ranges.append(currentRange)
        ranges = [range(x[0], x[-1]) for x in ranges]
        # Get get geomtric mean of p-values within a range
        pvals_fin = []
        tmp = []
        for r in ranges:
            for i in r:
                tmp.append(p_vals_dict[i])
            pvals_fin.append(gmean(tmp))
        ranges = [list(x) for x in ranges]
        #return(ranges)
        ranges = [[x[0], x[-1]] for x in ranges]
        out = pd.DataFrame({"sequenceID" : [list(df['sequenceID'])[0]] * len(ranges),
                        "start" : [x[0] for x in ranges],
                        "end" : [x[1] for x in ranges],
                        'pvalue' : pvals_fin,
                            # CAREFUL: This only makes sense if you group by family, otherwise this will introduce bugs!
                        'family' : [list(df['family'])[0]] * len(ranges)})
        return(out)
