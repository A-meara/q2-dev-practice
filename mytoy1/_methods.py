# ----------------------------------------------------------------------------
# Copyright (c) 2025, AMeara.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from skbio.alignment import global_pairwise_align_nucleotide, TabularMSA
from skbio import DNA

from q2_types.feature_data import DNAIterator


def duplicate_table(table: pd.DataFrame) -> pd.DataFrame:
    return table

#type hints using : and -> after the argument
#uses workaround to get first dna sequence from multiple seqs
def nw_align(seq1: DNA,
             seq2: DNA,
             gap_open_penalty: float = 5,
             gap_extend_penalty: float = 2,
             match_score: float = 1,
             mismatch_score: float = -2) -> TabularMSA:

    msa, _, _ = global_pairwise_align_nucleotide(
        seq1=seq1, seq2=seq2, gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty, match_score=match_score,
        mismatch_score=mismatch_score
    )

    return msa

# want sequence file as input
#
def seqcount(sequences: DNAIterator) -> int:

    # outputstr = "number of sequences: %d" % count
    # click.secho(outputstr,fg='green', bg='black')
    # click.echo("hello")
    seqs = list(sequences)  # consumes the iterator
    return len(seqs)