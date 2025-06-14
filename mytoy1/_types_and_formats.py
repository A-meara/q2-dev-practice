# ----------------------------------------------------------------------------
# Copyright (c) 2025, AMeara.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from skbio import DNA
from skbio.io import UnrecognizedFormatError

from qiime2.plugin import SemanticType, TextFileFormat, ValidationError, model

SingleDNASequence = SemanticType("SingleDNASequence")

TotalSeqCount = SemanticType("TotalSeqCount")


class SingleRecordDNAFASTAFormat(TextFileFormat):

    def _confirm_single_record(self):
        with self.open() as fh:
            try:
                # DNA.read(..., validate = False) disables checking to ensure
                # that all characters in the sequence are IUPAC DNA characters
                # by scikit-bio.
                # This will be validated independently by _validate_, with user
                # control over how much of the sequence is read during
                # validation, to manage the runtime of validation.
                _ = DNA.read(fh, seq_num=1, validate=False)
            except UnrecognizedFormatError:
                raise ValidationError(
                    "At least one sequence record must be present, but none "
                    "were found."
                )

            try:
                _ = DNA.read(fh, seq_num=2, validate=False)
            except ValueError:
                # if there is no second record, a ValueError should be raised
                # when we try to access the second record
                pass
            else:
                raise ValidationError(
                    "At most one sequence record must be present, but more "
                    "than one record was found."
                )

    def _confirm_acgt_only(self, n_chars):
        with self.open() as fh:
            seq = DNA.read(fh, seq_num=1, validate=False)
            validation_seq = seq[:n_chars]
            validation_seq_len = len(validation_seq)
            non_definite_chars_count = \
                validation_seq_len - validation_seq.definites().sum()
        if non_definite_chars_count > 0:
            raise ValidationError(
                f"{non_definite_chars_count} non-ACGT characters detected "
                f"during validation of {validation_seq_len} positions."
            )

    def _validate_(self, level):
        validation_level_to_n_chars = {'min': 50, 'max': None}
        self._confirm_single_record()
        self._confirm_acgt_only(validation_level_to_n_chars[level])

SingleRecordDNAFASTADirectoryFormat = model.SingleFileDirectoryFormat(
    'SingleRecordDNAFASTADirectoryFormat', 'sequence.fasta',
    SingleRecordDNAFASTAFormat)


class MySequenceCountFormat(TextFileFormat):
    # make sure only number in the format
    # havent included higher level, as it is just a toy
    def validate(self, level):
        with self.open() as fh:
            content = fh.read().strip()
            if not content.startswith("Sequences contained: "):
                raise ValidationError(
                    f"Expected line to start with 'Sequences contained: ', got: {content}"
                )
            try:
                count_str = content.split(":", 1)[1].strip()
                int(count_str)
            except (IndexError, ValueError):
                raise ValidationError(
                    f"Could not parse integer from line: {content}"
                )

SingleRecordSeqCountDirectoryFormat = model.SingleFileDirectoryFormat(
    'SingleRecordSeqCountDirectoryFormat', 'seqcount.txt',
    MySequenceCountFormat)

