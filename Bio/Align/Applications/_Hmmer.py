
from Bio.Application import _Option, AbstractCommandline, _Switch, _Argument
import os


class HmmScanCommandline(AbstractCommandline):
    """
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
           _Argument(["hmm"],
                      "HMM Library",
                      checker_function=os.path.exists,
                      #types=["file"],
                      is_required=True),
            _Argument(["input"],
                      "FASTA Query file",
                      checker_function=os.path.exists,
                      #types=["file"],
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
                      #types=["file"],
                      is_required=True),
           _Argument(["input"],
                      "FASTA Query file",
                      checker_function=os.path.exists,
                      #types=["file"],
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
                      #types=["file"],
                      is_required=True),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
