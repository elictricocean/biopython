# Copyright 2010 by Kyle Ellrott.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.Application import _Option, AbstractCommandline, _Switch, _Argument
import os


class HmmScanCommandline(AbstractCommandline):
    """
    >>> from Bio.Align.Applications import HmmScanCommandline
    >>> hmmer_bin="/usr/bin/hmmscan"
    >>> hmmPath="/dev/null"
    >>> seqPath="/dev/null"
    >>> hmmscan_cmdline = HmmScanCommandline(hmmer_bin, hmm=hmmPath, input=seqPath)
    >>> print hmmscan_cmdline
    /usr/bin/hmmscan /dev/null /dev/null    
    """
    def __init__(self, cmd="hmmscan", **kwargs):
        assert cmd is not None
        self.parameters = [
           _Switch(["--cut_ga", "cut_ga"],
                "Gathering Cutoff"),
           _Switch(["--cut_nc", "cut_nc"],
                "Noise Cutoff"),
           _Switch(["--cut_tc", "cut_tc"],
                "Trusted Cutoff"),
           _Switch(["-h", "help"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description;  ignore other arguments."),
           _Switch(["--acc", "accession"],
                    "prefer accessions over names in output"),
           _Option(["--cpu", "cpu"],
                    "number of parallel CPU workers to use for multithreads"),
           _Option(["-o", "out"],
                    "Output File", filename=True, equate=False ),
           _Argument(["hmm"],
                      "HMM Library",
                      checker_function=os.path.exists,
                      filename=True,
                      is_required=True),
            _Argument(["input"],
                      "FASTA Query file",
                      checker_function=os.path.exists,
                      filename=True,
                      is_required=True),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class HmmAlignCommandline(AbstractCommandline):
    """
    """
    def __init__(self, cmd="hmmalign", **kwargs):
        assert cmd is not None
        self.parameters = [
           _Switch(["-h", "help"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description;  ignore other arguments."),
           _Switch(["--acc", "accession"],
                    "prefer accessions over names in output"),
           _Option(["--cpu", "cpu"],
                    "number of parallel CPU workers to use for multithreads"),
           _Argument(["hmm"],
                      "HMM Library",
                      checker_function=os.path.exists,
                      filename=True,
                      is_required=True),
           _Argument(["input"],
                      "FASTA Query file",
                      checker_function=os.path.exists,
                      filename=True,
                      is_required=True),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class HmmPressCommandline(AbstractCommandline):
    """
    """
    def __init__(self, cmd="hmmpress", **kwargs):
        assert cmd is not None
        self.parameters = [
           _Switch(["-h", "help"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description;  ignore other arguments."),
           _Switch(["-f", "force"],
                    "force: overwrite any previous pressed files"),
           _Argument(["hmm"],
                      "HMM Library",
                      checker_function=os.path.exists,
                      filename=True,
                      is_required=True),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
