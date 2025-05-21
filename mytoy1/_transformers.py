from skbio import DNA

from mytoy1 import SingleRecordDNAFASTAFormat, MySequenceCountFormat

from .plugin_setup import plugin


# Define and register transformers
@plugin.register_transformer
def _1(ff: SingleRecordDNAFASTAFormat) -> DNA:
    # by default, DNA.read will read the first sequence in the file
    with ff.open() as fh:
        return DNA.read(fh)

@plugin.register_transformer
def _2(myint: int) -> MySequenceCountFormat:
    ff = MySequenceCountFormat()
    with ff.open() as fh:
        fh.write(f"Sequences contained: {myint}\n")
    return ff