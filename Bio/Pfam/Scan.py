

import re
import os
import subprocess
import string
from Bio import HMMER
import Bio.SeqIO
import StringIO


"""

This module requires HMMER3 and Pfam24.  It is based on the code found at 
ftp://ftp.sanger.ac.uk/pub/rdf/PfamScanBeta/


Pfam24:
ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam24.0/Pfam-A.hmm.gz
ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam24.0/Pfam-A.hmm.dat.gz

DB Prep:
# gunzip Pfam-A.hmm.gz
# gunzip Pfam-A.hmm.dat.gz
# hmmpress Pfam-A.hmm

HMMER 3:
ftp://selab.janelia.org/pub/software/hmmer3/hmmer-3.0b2.tar.gz

The PATH must contain 'hmmscan' in order to work

Module Usage:

>>> import Bio.Pfam.Scan.PfamScan
>>> ps = Bio.Pfam.Scan.PfamScan(dir=baseDir, file=fastaFile)
>>> ps.search()
>>> print ps

Or:

>>> ps = Bio.Pfam.Scan.PfamScan(dir=sys.argv[2])
>>> ps.setComment( False )
>>> handle = open( sys.argv[1] )
>>> seqDB = Bio.SeqIO.parse( handle, "fasta" )
>>> for seq in seqDB:
>>>     ps.scan(seq)
>>>     outStr = str(ps)
>>>     print outStr

Or:

>>> fileIn = gzip.GzipFile( sys.argv[1] )
>>> hmmReader = Bio.HMMER.HMMResultsIO()
>>> results = hmmReader.parseMultiHMMER3( fileIn )
>>> fileIn.close()
>>> pfam = Bio.Pfam.Scan.PfamScan( dir=sys.argv[2] )
>>> pfam.setResults( results )
>>> pfam.removeClanDup()
>>> print pfam.ascii_out( )

>>> pfam.calcSignificant()
>>> for res in results:
>>>     for unit in res:
>>>         if unit.sig:
>>>             print res.seqName, pfam.unitAcc( unit ), pfam.unitClan( unit ), unit.name, unit.evalue, unit.bits, unit.seqFrom, unit.seqTo

	
"""

def get_fasta(title, seq):
	"""A quick and dirty FASTA generator"""
	out_str = ">%s\n" % title
	for i in (range(0, len(seq)+1, 60)):
		out_str += "%s\n" % seq[i:i+60]
	return out_str


class MissingDBFile(Exception):
	def __init__(self, fileName):
		self.fileName = fileName
	def __str__(self):
		return "Missing DB file: %s" % (self.fileName)


class PfamScan:
	"""PfamScan Class.  Needs a baseDir where the Pfam_HMM, Pfam-A.scan.dat and active_site.dat files are located"""

	re_stockline = re.compile(r'^#=GF (..)\s+(.*)$')
	re_stockend  = re.compile(r'^//')
	re_ga = re.compile(r'^([^\s]*);')

	def __init__(self, dir=None, db="Pfam-A.hmm", file=None):
		"""
Setup Pfam Scan Parser
dir  : Name of directory where Pfam related files are kept
file : Path to fasta file to be scanned
db   : Database name.  By default Pfam-A.hmm.  
"""
		self._dir = dir
		self.file = file
		self._hmmlib = [ db ]
		self._as = False		
		self._read_ = {}
		self.comment = True
		
		for ext in [ "h3f", "h3i", "h3m", "h3p", "dat" ]:
			extPath = os.path.join( self._dir, "%s.%s" % (db, ext) )
			if ( not os.path.exists( extPath ) ) :
				raise MissingDBFile( extPath )

		self._read_pfam_data()

	def setComment(self, comment):
		self.comment = comment

	def search(self, file=None):
		if file is None:
			handle = open(self.file)
		else:
			handle = open(file)			
		seqdb = Bio.SeqIO.parse(handle, "fasta")
		for seq in seqdb:
			self.scan(seq)
		handle.close()

	def reset(self):
		self.results = None

	def scan(self, seq):
		"""Run hmmscan against a sequence"""

		#if self.pfamInfo is None:
		#	self.LoadStockInfo( "%s/Pfam-A.scan.dat" % self.baseDir )

		tmpPath = "/tmp/%s.%d" % ( "pfamScanHmmer", os.getpid() ) 
		tmpFile = open( tmpPath, "w" )
		tmpFile.write( get_fasta(seq.id, seq.seq) )
		tmpFile.close()
		dbPath = "%s/Pfam-A.hmm" % (self._dir)
		
		pipe = subprocess.Popen( ["hmmscan", "--notextw", dbPath, tmpPath] , stdout=subprocess.PIPE).stdout
		self.hmmParser = HMMER.HMMResultsIO()
		self.setResults( self.hmmParser.parseHMMER3( pipe ) )
		pipe.close()
		os.unlink( tmpPath )		
		return self.results

	def setResults(self, results):
		self.hmmParser = HMMER.HMMResultsIO()
		if not isinstance( results, list ):
			self.results = [ results ]
		else:
			self.results = results

	def calcSig(self):
		self.calcSignificant()
	
	def calcSignificant(self):
		self.loadPfam()
		for result in self.results:
			for unit in result.units:
					unit.sig = False
					if ( float(result.seqs[ unit.name ].bits) >= float(self._seqGA[ unit.name ]) ):
						if ( float(unit.bits) >= float(self._domGA[ unit.name ]) ):
							unit.sig = True

	def unitAcc(self, unit):
		self.loadPfam()
		return self._accmap[ unit.name ]

	def getPfamClan( self, acc ):
		self.loadPfam()
		try:
			return self._clanmap[ acc ]
		except KeyError:
			return "No_clan"
			
	def unitClan(self, unit):
		self.loadPfam()
		try:
			return self._clanmap[ unit.name ]
		except KeyError:
			return "No_clan"

	def __str__(self):
		self.calcSignificant()
		out = []
		if self.comment:
			out.append( "# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n#" ) 
			if self._as:
				out.append( "# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan> <predicted_active_site_residues>" )
			else:
				out.append( "# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>" )
		for res in self.results:
			outStr = self.hmmParser.write_ascii_out( res, self )
			if len( outStr ):
				out.append( outStr )
		return "%s" % ( "\n".join(out) )
		
	def removeClanDup(self):
		self.loadPfam()
		newresults = []
		for res in self.results:
			newresults.append( res.remove_overlaps_by_clan( self._clanmap, self._nested ) )
		self.results = newresults

	def ascii_out( self, **kwargs ):
		self.calcSignificant()
		out = []
		for res in self.results:
			outStr= self.hmmParser.write_ascii_out( res, self, **kwargs )
			if len( outStr ):
				out.append( outStr )
		return "\n".join(out)

	def loadPfam(self):
		if sorted(self._read_.keys()) != self._hmmlib:
			self._read_pfam_data() 

	def _read_pfam_data(self):			
		self._accmap    = {};
		self._nested    = {};
		self._clanmap   = {};
		self._desc      = {};
		self._seqGA     = {};
		self._domGA     = {};
		self._type      = {};
		self._model_len = {};
		id = None
		for hmmlib in self._hmmlib:
			scandat = self._dir + '/' + hmmlib + '.dat'
			SCANDAT = open( scandat )
			for line in SCANDAT:
				res = re.search( r'^\#=GF ID\s+(\S+)', line )
				if res:
					id = res.group(1)
					self._nested[ id ] = {}
					continue
				res = re.search( r'^\#=GF\s+AC\s+(\S+)', line )
				if res:
					self._accmap[id] = res.group(1)
					continue
				res = re.search( r'^\#=GF\s+DE\s+(.+)', line )
				if res:
					self._desc[id] = res.groups(1)
					continue
				res = re.search( r'^\#=GF\s+GA\s+(\S+)\;\s+(\S+)\;', line )
				if res:
					self._seqGA[ id ] = res.group(1)
					self._domGA[ id ] = res.group(2)
					continue
				res = re.search(r'^\#=GF\s+TP\s+(\S+)', line )
				if res:
					self._type[ id ] = res.group(1)
					continue
				res = re.search( r'^\#=GF\s+ML\s+(\d+)', line )
				if res:
					self._model_len[ id ] = res.group(1)
					continue
				res = re.search( r'^\#=GF\s+NE\s+(\S+)', line )
				if res:
					self._nested[ id ][ res.group(1) ] = 1
					try:
						self._nested[ res.group(1) ][ id ] = 1
					except KeyError:
						self._nested[ res.group(1) ] = { id:1 }						
					continue
				res = re.search( r'^\#=GF\s+CL\s+(\S+)', line )
				if res:
					self._clanmap[ id ] = res.group(1)
			SCANDAT.close()
			# set a flag to show that we've read the data files already
			self._read_[ hmmlib ] = 1

	
