# pylint: disable=C0103,C0301,C0116

""" module docstring """
import logging
import pathlib
import errno

# pylint: disable=W0611
from gqlib.db.db_import import SmallDatabaseImporter
from gqlib.profilers import RegionQuantifier
from gqlib.runners.alignment_runner import BwaMemRunner
from gqlib.ui.validation import check_bwa_index, check_input_reads

from .annotate.crazy_annotator import CazyAnnotator, HMMLoader, Sequences, ThresholdTable


logger = logging.getLogger(__name__)

def run_profile(args):

    # ------------------------ validate inputs and args --------------------------------
    if args.input_type == "fastq":
        input_data = check_input_reads(
            fwd_reads=args.reads1,
            rev_reads=args.reads2,
            single_reads=args.singles,
            orphan_reads=args.orphans,
        )

        if not check_bwa_index(args.bwa_index):
            raise ValueError(f"{args.bwa_index} is not a valid bwa index.")

    else:
        input_file = args.bam if args.input_type == "bam" else args.sam
        if not pathlib.Path(input_file).is_file():
            raise ValueError(f"{input_file} does not not exist!")

    if not pathlib.Path(args.annotation_db).is_file():
        raise ValueError(
            f"{args.annotation_db} is not a valid annotation database"
        )

    if pathlib.Path(args.out_prefix).parent != pathlib.Path('.'):
        pathlib.Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)

    db_format = None
    with open(args.annotation_db, "r") as db_in:
        print(args.annotation_db)
        for line in db_in:
            line = line.strip()
            if line and line[0] != "#":
                if line.find(",") != -1:
                    db_format = "hmmer"
                    break
                if line.find("\t") != -1:
                    db_format = "bed"
                    break
    if db_format is None:
        logger.error("Cannot determine database format in %s.", args.annotation_db)
        raise ValueError(f"Cannot determine database format in {args.annotation_db}.")
    
    logger.info("Identified database format as `%s`.", db_format)

    # --------------------------- Load Database and Quanifier --------------------------

    # for the domain mode counting, we need the SmallDatabaseImporter with "cazy"
    db_importer = SmallDatabaseImporter(
        logger, args.annotation_db, single_category="cazy", db_format=db_format,
    )
    logger.info("Finished loading database.")

    profiler = RegionQuantifier(
        db=db_importer,
        out_prefix=args.out_prefix,
        ambig_mode="1overN",
        reference_type="domain",
    )

    # ------------------------------ Dealing with FastQ input --------------------------
    if args.input_type == "fastq":

        # NOTE could also allow use of minimap2 here
        aln_runner = BwaMemRunner(
            args.threads,
            args.bwa_index,
            sample_id=pathlib.Path(args.out_prefix).name,
        )

        for input_type, *reads in input_data:

            logger.info("Running %s alignment: %s", input_type, ",".join(reads))
            proc, call = aln_runner.run(
                reads,
                single_end_reads=input_type == "single",
                alignment_file=args.keep_alignment_file,
            )

            try:
                profiler.count_alignments(
                    aln_stream=proc.stdout,
                    aln_format="sam",
                    min_identity=args.min_identity,
                    min_seqlen=args.min_seqlen,
                )

            except Exception as err:
                if isinstance(err, ValueError) and str(err).strip() == "file does not contain alignment data":
                    # pylint: disable=W1203
                    logger.error("Failed to align. This could have different reasons:")
                    logger.error(f"* Is `{args.aligner}` installed and on the path? Type `bwa mem` and see what happens.")
                    logger.error("* Syntax errors or missing files. Please try running the aligner call below manually to troubleshoot the problem.")
                    logger.error("* Alignment stream was interrupted, perhaps due to a memory issue.")
                    
                    logger.error("Aligner call was:")
                    logger.error("%s", call)
                    
                    return 1

                logger.error("Encountered problems digesting the alignment stream:")
                logger.error("%s", err)
                logger.error("Aligner call was:")
                logger.error("%s", call)
                logger.error("Shutting down.")
                
                return 1

    # -------------------- dealing with SAM or BAM input -------------------------------
    else:
        input_file = args.bam if args.input_type == "bam" else args.sam

        logger.info("Profiling %s alignment:", input_file)

        profiler.count_alignments(
            aln_stream=input_file,
            aln_format=args.input_type,
            min_identity=args.min_identity,
            min_seqlen=args.min_seqlen,
            external_readcounts=args.import_readcounts,
            unmarked_orphans=args.unmarked_orphans,
        )

    # ---------------------- Finalize the output ---------------------------------------

    profiler.finalise(
        restrict_reports=("raw", "rpkm",),
        report_category=False,
        report_unannotated=False,
        dump_counters=args.debug,
    )

    return 0


def run_proteome_annotation(args):

    if args.cutoffs is None:
        args.cutoffs = pathlib.Path(args.hmmdb, "cutoffs.csv")

    logger.info("Reading HMMs...")
    annotator = CazyAnnotator(
        hmms=HMMLoader.read_hmms(
            hmmdb_path=pathlib.Path(args.hmmdb),
            # file_with_paths=pathlib.Path(args.file_with_hmm_paths), # TODO
        )
    )
    logger.info("Reading sequences...")
    seqs = Sequences.read_sequences_from_file(path=args.proteins)
    
    logger.info("Loading Tresholds...")
    tresholds = ThresholdTable.load(args.cutoffs)

    logger.info("Annotating sequences (can take a few minutes; be patient)")
    results = annotator.annotate(
        sequences=seqs.sequences,
        threads=args.threads,
        blacklist=tresholds.hmms_which_will_be_skipped,
    )
    logger.info("Filtering and merging annotations over folds")
    annotations_filtered = annotator.curate_annotations(
        thresholds=tresholds,
        cazy_results=results,
    )
    annotations_filtered.to_csv(args.output_file, index=False)

    return 0
