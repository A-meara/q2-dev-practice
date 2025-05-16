# ----------------------------------------------------------------------------
# Copyright (c) 2025, AMeara.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from skbio.alignment import global_pairwise_align_nucleotide, TabularMSA

from q2_types.feature_data import DNAIterator, DNAFASTAFormat

import os
import click
import q2cli.util




def duplicate_table(table: pd.DataFrame) -> pd.DataFrame:
    return table

#type hints using : and -> after the argument
#uses workaround to get first dna sequence from multiple seqs
def nw_align(seq1: DNAIterator,
             seq2: DNAIterator,
             gap_open_penalty: float = 5,
             gap_extend_penalty: float = 2,
             match_score: float = 1,
             mismatch_score: float = -2) -> TabularMSA:
    seq1 = next(iter(seq1))
    seq2 = next(iter(seq2))

    msa, _, _ = global_pairwise_align_nucleotide(
        seq1=seq1, seq2=seq2, gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty, match_score=match_score,
        mismatch_score=mismatch_score
    )

    return msa

# want sequence file as input ->
# def seqcount(sequences: DNAFASTAFormat ) -> int:
#
#     count = 0
#     with open('sequences') as my_fasta:
#         for line in my_fasta:
#             if line.startswith('>'):
#                 count += 1
#     outputstr = "number of sequences: %d" % count
#     click.secho(outputstr,fg='green', bg='black')
#
#     return count