
import unittest
import Bio.HMMER

class TestCluster(unittest.TestCase):
	def test_outparse(self):		
		handle = open( "HMMER/103l_A.hmmscan" )
		results = Bio.HMMER.parse( handle )
		resArray = []
		for res in results:
			resArray.append( res )
		assert len(resArray)==1
		assert resArray[0][0].name=="Phage_lysozyme"
		assert resArray[0][0].evalue==3.1e-23
		assert (str(resArray[0][0])).rstrip()=="103l_A     24    150     24    151 N/A        Phage_lysozyme       N/A     1   109    -1     81.8  3.100e-23   1"

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
