# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Unittests for Bio.Align.Applications interface for MAFFT

This code is part of the Biopython distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.
"""

import sys
import os
import unittest
import subprocess
from Bio import MissingExternalDependencyError
from Bio.HMMER.Applications import HmmScanCommandline, HmmPressCommandline

hmmscan_exe=None
if sys.platform=="win32":
    raise MissingExternalDependencyError("Testing with HMMER3 not implemented on Windows yet")
else:
    import commands
    output = commands.getoutput("hmmscan -h")
    if "not found" not in output and "hmmscan" in output:
        hmmscan_exe = "hmmscan"
        
if not hmmscan_exe:
    raise MissingExternalDependencyError(\
        "Install HMMER3 if you want to use the Bio.HMMER.Applications wrapper.")


class HmmScanApplication(unittest.TestCase):

    def setUp(self):
        self.infile1 = "Fasta/f002"
        self.hmmfile = "HMMER/PF07980.hmm"
        cmdline = HmmPressCommandline(hmm=self.hmmfile, force=True)
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        child = subprocess.Popen(str(cmdline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        stdoutdata, stderrdata = child.communicate()
        
    def tearDown(self):
        if os.path.isfile("HMMER/PF07980.hmm.h3f"):
            os.remove("HMMER/PF07980.hmm.h3f")
        if os.path.isfile("HMMER/PF07980.hmm.h3i"):
            os.remove("HMMER/PF07980.hmm.h3i")
        if os.path.isfile("HMMER/PF07980.hmm.h3m"):
            os.remove("HMMER/PF07980.hmm.h3m")
        if os.path.isfile("HMMER/PF07980.hmm.h3p"):
            os.remove("HMMER/PF07980.hmm.h3p")

    def test_hmmscan(self):
        cmdline = HmmScanCommandline(hmm=self.hmmfile, input=self.infile1)
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        child = subprocess.Popen(str(cmdline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        stdoutdata, stderrdata = child.communicate()    

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
