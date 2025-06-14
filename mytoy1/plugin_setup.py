# ----------------------------------------------------------------------------
# Copyright (c) 2025, AMeara.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib
from qiime2.plugin import Citations, Plugin, Float, Range
from q2_types.feature_table import FeatureTable, Frequency
from mytoy1 import __version__
from mytoy1._methods import nw_align, duplicate_table, seqcount
from q2_types.feature_data import FeatureData, AlignedSequence, Sequence

from mytoy1 import (
    SingleDNASequence, SingleRecordDNAFASTAFormat,
    SingleRecordDNAFASTADirectoryFormat,
    TotalSeqCount, MySequenceCountFormat,
    SingleRecordSeqCountDirectoryFormat
)

citations = Citations.load("citations.bib", package="mytoy1")

plugin = Plugin(
    name="mytoy1",
    version=__version__,
    website="https://example.com",
    package="mytoy1",
    description="https://github.com/A-meara",
    short_description="first draft",
    # The plugin-level citation of 'Caporaso-Bolyen-2024' is provided as
    # an example. You can replace this with citations to other references
    # in citations.bib.
    citations=[citations['Caporaso-Bolyen-2024']]
)

plugin.methods.register_function(
    function=duplicate_table,
    inputs={'table': FeatureTable[Frequency]},
    parameters={},
    outputs=[('new_table', FeatureTable[Frequency])],
    input_descriptions={'table': 'The feature table to be duplicated.'},
    parameter_descriptions={},
    output_descriptions={'new_table': 'The duplicated feature table.'},
    name='Duplicate table',
    description=("Create a copy of a feature table with a new uuid. "
                 "This is for demonstration purposes only. 🧐"),
    citations=[]
)

plugin.methods.register_function(
    function=nw_align,
    inputs={'seq1': SingleDNASequence,
            'seq2': SingleDNASequence},
    parameters={
        'gap_open_penalty': Float % Range(0, None, inclusive_start=False),
        'gap_extend_penalty': Float % Range(0, None, inclusive_start=False),
        'match_score': Float % Range(0, None, inclusive_start=False),
        'mismatch_score': Float % Range(None, 0, inclusive_end=True)},
    outputs={'aligned_sequences': FeatureData[AlignedSequence]},
    input_descriptions={'seq1': 'The first sequence to align.',
                        'seq2': 'The second sequence to align.'},
    parameter_descriptions={
        'gap_open_penalty': ('The penalty incurred for opening a new gap. By '
                             'convention this is a positive number.'),
        'gap_extend_penalty': ('The penalty incurred for extending an existing '
                               'gap. By convention this is a positive number.'),
        'match_score': ('The score for matching characters at an alignment '
                        'position. By convention, this is a positive number.'),
        'mismatch_score': ('The score for mismatching characters at an '
                           'alignment position. By convention, this is a '
                           'negative number.')},
    output_descriptions={
        'aligned_sequences': 'The pairwise aligned sequences.'
    },
    name='Pairwise global sequence alignment.',
    description=("Align two DNA sequences using Needleman-Wunsch (NW). "
                 "This is a Python implementation of NW, so it is very slow! "
                 "This action is for demonstration purposes only. 🐌"),
    citations=[citations['Needleman1970']]
)

# Register semantic types
plugin.register_semantic_types(SingleDNASequence,TotalSeqCount)

# Register formats
plugin.register_formats(SingleRecordDNAFASTAFormat,
                        SingleRecordDNAFASTADirectoryFormat,
                        MySequenceCountFormat,
                        SingleRecordSeqCountDirectoryFormat
                        )

# Define and register new ArtifactClass
plugin.register_artifact_class(SingleDNASequence,
                               SingleRecordDNAFASTADirectoryFormat,
                               description="A single DNA sequence.")

plugin.register_artifact_class(TotalSeqCount,
                               SingleRecordSeqCountDirectoryFormat,
                               description="Total DNA sequence count.")

plugin.methods.register_function(
    function=seqcount,
    inputs={'sequences': FeatureData[Sequence]},
    parameters={},
    outputs={'seq_count': TotalSeqCount},
    input_descriptions={'sequences': 'File containing the sequences to be counted'},
    parameter_descriptions={},
    #output_descriptions={
    #    'aligned_sequences': 'The pairwise aligned sequences.'
    #},
    name='Sequence Counter',
    description="Count my sequences!!"
    #citations=[citations['Needleman1970']]
)

importlib.import_module('mytoy1._transformers')

