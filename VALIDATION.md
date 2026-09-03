# Validation notes

The repository follows the same practical reproducibility ideas used in `mice_transplant_2025`, while keeping the measles workflow much smaller.

Before Slurm submission, `scripts/validate_repository.py` checks:

- all required source files;
- clean JSON notebooks with unique cell IDs;
- Python syntax of every script;
- exactly 96 new study-table records and 16 legacy UDP records;
- full patient identifiers, including the distinction between `СА1_013` and `СА1_017`;
- sequencing-index and `.1/.2/.3` timepoint consistency;
- the new healthy-control mapping;
- legacy UDP0097-UDP0105 case and UDP0106-UDP0112 control rules;
- UDP0108 exclusion from biological inference;
- pinned Python/MiXCR/matplotlib-venn environment;
- required preprocessing and longitudinal-analysis logic.

During the real run, preprocessing additionally checks the filesystem. Every new included sample must have four complete R1/R2 lane pairs L001-L004. Every included legacy sample must be found. A full run stops on missing input instead of producing a partial result set.

After notebook execution, the runner verifies the required result tables, figures 01-10 and patient aaV Venn outputs.

The final remaining external prerequisites are the real Aldan-3 FASTQ files, the correct old `pe150_res` location and a valid MiXCR license for the user account.
