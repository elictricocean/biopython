
"""
Based on code from PfamScan from: ftp://ftp.sanger.ac.uk/pub/rdf/PfamScanBeta/

Copyright (c) 2009: Genome Research Ltd.

Authors: Jaina Mistry (jm14@sanger.ac.uk), John Tate (jt6@sanger.ac.uk)

This is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
or see the on-line version at http://www.gnu.org/copyleft/gpl.txt


"""

import re
import string
import os
import subprocess
import Bio.Seq
import Bio.SeqRecord


def parseHMMER3( fh ):
	io = HMMResultsIO()
	return io.parseHMMER3( fh )

def parseMultiHMMER3(fh):
	io = HMMResultsIO()
	return io.parseMultiHMMER3( fh )
		
class FormatError(Exception):
	def __init__(self, message):
		self.message = message
	def __str__(self):
		return str(self.message)


class EOF(Exception):
	pass


class MissingDBFile(Exception):
	def __init__(self, fileName):
		self.fileName = fileName
	def __str__(self):
		return "Missing DB file: %s" % (self.fileName)


instanceNum = 0

class HMMER3Scan:
	"""HMMER3 Scan Program Control Class."""

	re_stockline = re.compile(r'^#=GF (..)\s+(.*)$')
	re_stockend  = re.compile(r'^//')
	re_ga = re.compile(r'^([^\s]*);')

	def __init__(self, db, score=None):
		"""
Setup HMMERPfam Scan Parser
db   : Path to HMMER database
"""
		global instanceNum
		self._hmmlib = [ db ]
		self._as = False		
		self._read_ = {}
		self.score = score
		self.cutoff = None
		self.comment = True
		self.params = []
		self.results = []
		self.instanceNum = instanceNum
		instanceNum += 1
		for ext in [ "h3f", "h3i", "h3m", "h3p" ]:
			extPath = "%s.%s" % (db, ext)
			if ( not os.path.exists( extPath ) ) :
				raise MissingDBFile( extPath )

	def cut_ga(self):
		"""use profile's GA gathering cutoffs to set all thresholding"""
		self.score = None
		self.cutoff = "--cut_ga"

	def cut_nc(self):
		"""use profile's NC noise cutoffs to set all thresholding"""
		self.score = None
		self.cutoff = "--cut_nc"

	def cut_tc(self):
		"""use profile's TC trusted cutoffs to set all thresholding"""
		self.score = None
		self.cutoff = "--cut_tc"

	def use_acc(self):
		"""prefer accessions over names in alignment units"""
		self.params.extend( "--acc" )

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
		self.results = []

	def scan(self, seqIn):
		self.reset()		
		scanSet = []
		if isinstance(seqIn, Bio.SeqRecord.SeqRecord):
			scanSet = [ seqIn ]
		elif isinstance(seqIn, str):
			scanSet = [ Bio.SeqRecord( Bio.Seq.Seq(seqIn), id="query" ) ]
		elif isinstance(seqIn, Bio.Seq.Seq):
			scanSet = [ Bio.SeqRecord( seqIn, id="query" ) ]
			
		else:
			scanSet = seqIn	
			
		"""Run hmmscan against a sequence"""
		tmpPath = "/tmp/%s.%d.%d" % ( "pfamScanHmmer", os.getpid(), self.instanceNum ) 
		tmpFile = open( tmpPath, "w" )
		for seq in scanSet:
			tmpFile.write( seq.format("fasta") )
		tmpFile.close()
			
		for db in self._hmmlib:
			cmdline = ["hmmscan", "--notextw" ]
			if self.score is not None:
				cmdline.extend( ['-T', "%e" % (self.score) ] )
			if self.cutoff is not None:
				cmdline.append( self.cutoff )
			cmdline.extend( [ db, tmpPath] )
			#cmdline.extend( self.params )
			pipe = subprocess.Popen( cmdline, stdout=subprocess.PIPE).stdout
			for result in parseMultiHMMER3( pipe ):
				yield result
			pipe.close()
		os.unlink( tmpPath )

class HMMIO:
	"""This is a dummy class right now.  Inserted to follow the structure of the Perl Bio.Pfam.HMM"""
	def __init__(self):
		pass
		
	def readHMM(self, hmm):
		pass
	
	def writeHMM(self):
		pass


class HMMMatch:
	"""This is a dummy class right now.  Inserted to follow the structure of the Perl Bio.Pfam.HMM"""
	def __init__(self):
		self.evalue = None
		self.bits = None	
		self.name = None
		self.bias = None

class HMM:
	def __init__(self):
		self.hmmVersion = None
		self.hmmName = None
		self.hmmAcc = None 
		self.hmmAlpha = None
		self.hmmMsvStats = None
		self.hmmViterbiStats = None
		self.hmmForwardStats = None
		self.version = None
		self.name = None
		self.accession = None
		self.description = None
		self.length = None
		self.alpha = None
		self.rf = None
		self.cs = None
		self.map = None
		self.date = None
		self.buildLine = None
		self.searchMethod = None
		self.nSeq = None
		self.msvStats = None
		self.viterbiStats = None
		self.forwardStats = None
		self.effn = None
		self.cksum = None
		self.seqGA = None
		self.domGA = None
		self.seqTC = None
		self.domTC = None
		self.seqNC = None
		self.domNC = None
		self.emissionLines = None
		self.mapPos = None
		self.compLines = None

class HMMResults:
	def __init__(self, **kwargs):
		self.seedName = None
		self.seqName = None
		self.program = None
		self.units = []		
		self.seqs = {}
		for a in kwargs:
			if hasattr( self, a ):
				setattr( self, a, kwargs[a] )
			else:
				raise NameError( "arg %s not recognized" % (a) )

	#def __str__(self):
	#	out = []
	#	for a in self.units:
	#		out.append( str( a ) )
	#	return "\n".join( out )

	def __iter__(self):
		self.curUnit = 0
		return self
		
	def next(self):
		self.curUnit += 1
		try:
			return self.units[ self.curUnit - 1 ]
		except IndexError:
			raise StopIteration

	def addHMMSeq(self, hmmSeq):
		self.seqs[ hmmSeq.name ] = hmmSeq	

	def eachHMMSeq(self):
		pass

	def addHMMUnit(self, hmmunit):
		self.units.append( hmmunit )
	
	def domainBitsCutoffFromEvalue(self):
		pass

	def lowestTrue(self):
		pass

	def highestNoise(self):
		pass

	def applyEdits(self):
		pass

	def remove_overlaps_by_clan(self, clanmap, nested):
		new = HMMResults()
		new.seqName = self.seqName
		for unit in sorted( self.units, key=lambda a: float(a.evalue) ):
			#check if it overlaps before adding
			o = False
			for u in new.units:
				if( clanmap.has_key( unit.name ) and clanmap.has_key( u.name ) and (clanmap[ unit.name ] == clanmap[ u.name ] ) ):
					if( self.overlap( unit, u ) ):
						if nested[ unit.name ].has_key( u.name ):
							continue
						else:
							o=True
							break
			if not o:
				if not new.seqs.has_key( self.seqs[unit.name].name ):
					new.addHMMSeq( HMMSequence( 
						name       = self.seqs[unit.name].name,
						desc       = self.seqs[unit.name].desc,
						bits       = self.seqs[unit.name].bits,
						evalue     = self.seqs[unit.name].evalue,
						numberHits = self.seqs[unit.name].numberHits
					) )
				new.addHMMUnit(unit)
		return new
		
		
	def overlap(self, unit1, unit2):
		[ u1, u2 ] = sorted( [ unit1, unit2 ], key=lambda a: int(a.seqFrom) ) 
		if( int(u2.seqFrom) <= int(u1.seqTo) ):
			return True
		return False


class HMMResultsIO:
	"""HMMResult parsing class"""
	re_query = re.compile(r'^Query:\s*([^\s]*).*\[L=(.*)\]')
	re_targetStart = re.compile(r'^>>\s*([^\s]*)')
	re_fileEnd = re.compile(r'^Internal pipeline statistics summary:')
	re_space = re.compile(r'\s+')

	re_stockline = re.compile(r'^#=GF (..)\s+(.*)$')
	re_stockend  = re.compile(r'^//')
	re_ga = re.compile(r'^([^\s]*);')

	def __init__(self):
		pass
	
	def parseHMMER3( self, fh ):
		hmmRes = HMMResults()
		self._readHeader( fh,  hmmRes )
		self._readSeqHits( fh, hmmRes )
		self._readUnitHits( fh, hmmRes )
		self._readFooter( fh, hmmRes)
		return hmmRes

	def parseMultiHMMER3(self, fh):
		hmmResAll = []
		program = None
		try:
			while(1):
				hmmRes = HMMResults()
				self._readHeader( fh, hmmRes )
				if(hmmRes.program):
					program = hmmRes.program
				else:
					hmmRes.program = program
				self._readSeqHits( fh, hmmRes )
				self._readUnitHits( fh, hmmRes )
				#hmmResAll.append( hmmRes )
				self._readFooter(fh, hmmRes)
				yield hmmRes
		except EOF:
			pass
		#return hmmResAll

	def parseSplitHMMER3(self, fh):
		hmmRes = HMMResults()
		for filename in files:
			fh = open( filename )
			self._readHeader( fh, hmmRes )
			self._readSeqHits( fh, hmmRes )
			self._readUnitHits( fh, hmmRes )
			self._readFooter(fh, hmmRes)
		return hmmRes

	
	def convertHMMSearch(self):
		pass

	def writePFAMOUT(self):
		pass

	def parsePFAMOUT(self):
		pass


	def _readHeader( self, hs, hmmRes ):
		while ( 1 ):
			line = hs.readline()
			if not line:
				raise EOF
			if ( line.startswith( "Scores for complete" ) ):
				return				
			res = re.search(r'^# query HMM file:\s+(\S+)', line)
			if ( res ):
				hmmRes.hmmName = res.group(1)
				continue				
			res = re.search(r'^# target sequence database:\s+(\S+)', line)
			if ( res ):
				hmmRes.seqDB = res.group(1)
				continue				
			res = re.search( r'^output directed to file:\s+(\S+)', line)
			if res:
				hmmRes.thisFile = res.group(1)
				continue
			res = re.search( r'^Query:\s+(\S+)\s+\[M\=(\d+)\]', line) 
			if res:
				hmmRes.seedName = res.group(1)
				hmmRes.hmmLength = res.group(2)
				continue
			res = re.search( r'^Query:\s+(\S+)\s+\[L\=(\d+)\]', line )
			if res:
				hmmRes.seqName = res.group(1)
				hmmRes.seqLength = res.group(2)
				continue				
			res = re.search(r'^sequence E-value threshold: <= (\d+)', line)
			if res:
				hmmRes.evalueThr = res.group(1)
				continue
			res = re.search(r'^# Random generator seed:      (\d+)', line)
			if res:
				hmmRes.randSeedNum = res.group(1)
				continue
			res = re.search(r'^Description:\s+(.*)', line)
			if res:
				hmmRes.description = res.group(1)
				continue
			res = re.search(r'^# (phmmer|hmmsearch|hmmscan)', line)
			if res:
				hmmRes.program = res.group(1)
				continue
			#raise FormatError( "Failed to parse %s in sequence section\n" % (line) )

	def _readSeqHits( self, hs, hmmRes ):
		while (1):
			line = hs.readline()
			if line is None:
				raise EOF
#Match a line like this
# E-value  score  bias    E-value  score  bias    exp  N  Sequence Description
#    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------
#      4e-83  285.8  10.0    5.3e-83  285.5   7.0    1.1  1  Q14SN3.1 Q14SN3_9HEPC Polyprotein (Fragment).
			res = re.search(r'^Domain and alignment annotation for each [sequence|model]', line )
			if res:
				return
			res = re.search(r'^Domain annotation for each model ', line)
			if res:
				return
				
			res = re.search( r'^\s+(E-value|---)', line )
			if res:
				continue
			res = re.search(r'^$', line )
			if res:
				continue
			res = re.search(r'No hits detected that satisfy reporting thresholds', line)
			if res:
				continue
				
			#Assume that we have a sequence match
			sMatch = re.split( r'\s+', line )
			if not (len( sMatch ) >= 10 ):
				raise FormatError( "Expected at least 10 pieces of data: %s" % (line) )
			desc = "-"
			if ( len(sMatch) >= 11 ):
				desc = " ".join( sMatch[ 10: ] )
			hmmRes.addHMMSeq(
				HMMSequence(
						evalue     = sMatch[1],
						bits       = sMatch[2],
						bias       = sMatch[3],
						exp        = sMatch[7],
						numberHits = sMatch[8],
						name       = sMatch[9],
						desc       = desc         
				)
			)
			
	def _readUnitHits( self, hs, hmmRes ):
#Parse the domain hits section
#>> P37935.1  MAAY4_SCHCO Mating-type protein A-alpha Y4.
#     # bit score    bias    E-value ind Evalue hmm from   hmm to    ali from   ali to    env from   env to    ali-acc
#   --- --------- ------- ---------- ---------- -------- --------    -------- --------    -------- --------    -------
#     1     244.0     0.5    9.5e-76    1.7e-70        1      146 [.        1      145 [.        1      146 [.    0.99
#
#  Alignments for each domain:
#  == domain 1    score: 244.0 bits;  conditional E-value: 9.5e-76
#      SEED   1 medrlallkaisasakdlvalaasrGaksipspvkttavkfdplptPdldalrtrlkeaklPakaiksalsayekaCarWrsdleeafdktaksvsPanlhllealrirlyteqvekWlvqvlevaerWkaemekqrahiaatmgp 146
#               m+++la+l++isa+akd++ala+srGa+++ +p++tt+++fd+l++P+ld++rtrl+ea+lP+kaik++lsaye+aCarW++dleeafd+ta+s+sP+n+++l++lr+rly+eqv+kWl++vl+v+erWkaemekqrahi+atmgp
#  P37935.1   1 MAELLACLQSISAHAKDMMALARSRGATGS-RPTPTTLPHFDELLPPNLDFVRTRLQEARLPPKAIKGTLSAYESACARWKHDLEEAFDRTAHSISPHNFQRLAQLRTRLYVEQVQKWLYEVLQVPERWKAEMEKQRAHINATMGP 145
#               899***************************.******************************************************************************************************************8 PP

		while (1):
			line = hs.readline()
			if not line:
				raise EOF
			res = re.search( r'^Internal' , line )
			if res:
				return
		
			res = re.search( r'\>\>\s+(\S+)', line )
			if res:
				seqId = res.group(1)
				self._readUnitData( seqId, hs, hmmRes )
				if hmmRes.eof:
					break
				continue
			

	def _readUnitData(self, id, hs, hmmRes ):
		hmmName = hmmRes.seedName
		seqName = hmmRes.seqName
# bit score    bias    E-value ind Evalue hmm from   hmm to    ali from   ali to    env from   env to    ali-acc
#   --- --------- ------- ---------- ---------- -------- --------    -------- --------    -------- --------    -------
#     1     244.0     0.5    9.5e-76    1.7e-70        1      146 [.        1      145 [.        1      146 [.    0.99
#
#  Alignments for each domain:

		units = []
		align   = 1;
		recurse = 0;
		eof = 0;
		nextSeqId = None
		while 1:
			line = hs.readline()
			if not line:
				raise EOF
			
			res = re.search(r'^[(\/\/|Internal)]', line )
			if res:
				align   = 0
				recurse = 0
				eof = 1
				continue
			res = re.search( r'^\>\>\s+(\S+)', line)
			if res:
				nextSeqId = res.group(1)
				align     = 0
				recurse   = 1
				break
			res = re.search( r'\[No individual domains that satisfy reporting thresholds', line )
			if res:
				break
			res = re.search( r'^\s+Alignments for each domain:', line)
			if res:
				align   = 1
				recurse = 0
				break
			res = re.search(r'^\s+(#\s+score|---)', line)
			if res:
				continue
			res = re.search( r'^$', line )
			if res:
				continue
			res = re.search( r'^\s+\d+\s+', line )
			if res:
				dMatch = re.split( r'\s+', string.rstrip(line) )
				if len(dMatch) != 17:
					raise FormatError( "Expected 16 elements of datam, got %d: %s " % (len(dMatch), line) )
				units.append(
					HMMUnit(
					    seqName   = seqName,
						name      = id,
						domain    = dMatch[1],
						bits      = float(dMatch[3]),
						bias      = dMatch[4],
						domEvalue = float(dMatch[5]),
						evalue    = dMatch[6],
						hmmFrom   = int(dMatch[7]),
						hmmTo     = int(dMatch[8]),
						seqFrom   = int(dMatch[10]),
						seqTo     = int(dMatch[11]),
						envFrom   = int(dMatch[13]),
						envTo     = int(dMatch[14]),
						aliAcc    = dMatch[16]
					)
				)
				continue
			raise FormatError( "Did not parse line: %s" % (line) );

#  == domain 1    score: 244.0 bits;  conditional E-value: 9.5e-76
#      SEED   1 medrlallkaisasakdlvalaasrGaksipspvkttavkfdplptPdldalrtrlkeaklPakaiksalsayekaCarWrsdleeafdktaksvsPanlhllealrirlyteqvekWlvqvlevaerWkaemekqrahiaatmgp 146
#               m+++la+l++isa+akd++ala+srGa+++ +p++tt+++fd+l++P+ld++rtrl+ea+lP+kaik++lsaye+aCarW++dleeafd+ta+s+sP+n+++l++lr+rly+eqv+kWl++vl+v+erWkaemekqrahi+atmgp
#  P37935.1   1 MAELLACLQSISAHAKDMMALARSRGATGS-RPTPTTLPHFDELLPPNLDFVRTRLQEARLPPKAIKGTLSAYESACARWKHDLEEAFDRTAHSISPHNFQRLAQLRTRLYVEQVQKWLYEVLQVPERWKAEMEKQRAHINATMGP 145
#               899***************************.******************************************************************************************************************8 PP
#
# OR....
#
#  == domain 1    score: 27.6 bits;  conditional E-value: 7.4e-10
#   PF00018  17 LsfkkGdvitvleksee.eWwkaelkdg.keGlvPsnYvep 55 
#               L++++Gd+++++++++e++Ww++++++++++G++P+n+v+p
#  P15498.4 617 LRLNPGDIVELTKAEAEqNWWEGRNTSTnEIGWFPCNRVKP 657
#               7899**********9999*******************9987 PP

		if (align):
			pattern1 = None
			pattern2 = None
			if ( hmmName and hmmRes.program == 'hmmsearch'):
				pattern1 = re.compile(r'^\s+%s\s+\d+\s+(\S+)\s+\d+' % (hmmName) )
				pattern2 = re.compile(r'^\s+%s\s+\d+\s+(\S+)\s+\d+' % (id) )
			if ( seqName and hmmRes.program == 'hmmscan'):
				tmpSeqName = seqName
				tmpSeqName = re.sub(r'\|', r'\\|', tmpSeqName)
				pattern1 = re.compile(r'^\s+%s\s+\d+\s+(\S+)\s+\d+' % (id))
				pattern2 = re.compile(r'^\s+%s\s+\d+\s+(\S+)\s+\d+' % (tmpSeqName) )
			if ( seqName and hmmRes.program == 'phmmer'):
				seqName = re.sub(r'\|', r'\\|', seqName)
				id = re.sub( r'\|', r'\\|', id)
				pattern1 = re.compile( r'^\s+%s\s+\d+\s+(\S+)\s+\d+' % (seqName))
				pattern2 = re.compile( r'^\s+%s\s+\d+\s+(\S+)\s+\d+/' % (id))
			
			recurse = 0
			matchNo = None
			hmmlen = 0
			while 1:
				line = hs.readline()
				if not line:
					raise EOF
				#print string.strip(line)
				res = pattern1.search( line )
				if res:
					try:
						units[ matchNo - 1 ].hmmalign['hmm'] += res.group( 1 )
					except KeyError:
						units[ matchNo - 1 ].hmmalign['hmm'] = res.group( 1 )
					hmmlen = len( res.group(1) )
					continue
				res = pattern2.search( line )
				if res:
					try:
						units[ matchNo - 1 ].hmmalign['seq'] += res.group( 1 )
					except KeyError:
						units[ matchNo - 1 ].hmmalign['seq'] = res.group( 1 )
					continue
				res = re.search(r'^\s+([0-9\*\.]+)\s+PP$', line) 
				if res:
					pp = res.group(1)
					try:
						units[ matchNo - 1 ].hmmalign['pp'] += pp
					except KeyError:
						units[ matchNo - 1 ].hmmalign['pp'] = pp						
					continue
				res = re.search( r'^\s+==\s+domain\s+(\d+)', line )
				if res:
					matchNo = int( res.group(1) )
					#print "Match", matchNo
					continue
				res = re.search('^\s+(.*)\s+$', line )
				if res:
					# $1 is *not* the match - this fails if there are prepended
					# or appended spaces
					# $units[ $matchNo - 1 ]->hmmalign->{match} .= $1;
					# Let's get a right substring based on the HMM length
					m1 = line[-hmmlen:]
					try:
						units[ matchNo - 1 ].hmmalign['match'] += m1
					except KeyError:
						units[ matchNo - 1 ].hmmalign['match'] = m1						
					continue
				res = re.search( r'^$', line )
				if res:
					continue
				res = re.search( '^[(\/\/|Internal)]', line )
				if res:
					align   = 0
					recurse = 0
					eof = 1
					break
				res = re.search( r'^\>\>\s+(\S+)', line)
				if res:
					nextSeqId = res.group(1)
					recurse   = 1
					break
				raise FormatError( "Did not parse |%s| in units" % (line) )

		hmmRes.eof = eof
		for u in units:
			hmmRes.addHMMUnit(u)

		if (recurse and nextSeqId):
			self._readUnitData( nextSeqId, hs, hmmRes )
       
	def parseHMMER2(self):
		pass

	def parseHMMER1(self):
		pass

	def writeScoresFile(self):
		pass

	def _readAlign(self):
		pass

	def _readFooter(self, fh, hmmRes ):
		while 1:
			line = fh.readline()
			if not line:
				raise EOF
			res = re.search(r'\/\/', line )
			if res:
				break

	def write_ascii_out(self, HMMResults, scanData, ga=False, e_seq=None, e_dom=None, b_seq=None, b_dom=None):
		scanData._max_seqname = 20 # unless($scanData->{_max_seqname} or $scanData->{_max_seqname} < 1);
		if (e_seq or e_dom):
			if not e_seq:
				e_seq = e_dom
			if not e_dom:
				e_dom = "10"
		elif(b_seq or b_dom):
			if not b_seq:
				b_seq = b_dom
			if not b_dom:
				b_dom = "0" 
		else:
			ga = True
		out_array = []
		#for unit in ( sorted( HMMResults.units, lambda x,y : x.seqFrom < y.seqFrom) ):  
		for unit in HMMResults.units:
			if( re.search(r'Pfam\-B', unit.name) ):
				if not ( float(HMMResults.seqs[ unit.name ].evalue) <= 0.001 and float(unit.evalue <= 0.001) ):
					continue
				out_array.append(  ("%-" + str(scanData._max_seqname) + "s %6d %6d %6d %6d %-10s %-16s %7s %5d %5d %5d %8s %9s %3s %-8s\n") % \
					( 
						HMMResults.seqName,
						unit.seqFrom,
						unit.seqTo,
						unit.envFrom,
						unit.envTo,
						scanData._accmap[ unit.name ],
						unit.name,
						"Pfam-B",
						unit.hmmFrom,
						unit.hmmTo,
						scanData._model_len[ unit.name ],
						unit.bits,
						unit.evalue,
						"NA",
						"NA"
					)
					)
			else:
				#Filter results based on thresholds
				if(ga):
					if ( not unit.sig ):
						continue
				if (e_seq):
					if not (float(HMMResults.seqs[ unit.name ].evalue) <= float(e_seq) and float(unit.evalue) <= float(e_dom) ):
						continue
				if(b_seq):
					if not (float(HMMResults.seqs[ unit.name ].bits) >= b_seq and float(unit.bits) >= float(b_dom)):
						continue
	    
				#print unit
				clan = "No_clan" 
				try:
					clan = scanData._clanmap[ unit.name ]
				except KeyError:
					pass
				out_array.append( ("%-" + str(scanData._max_seqname) + "s %6d %6d %6d %6d %-10s %-16s %7s %5d %5d %5d %8s %9s %3d %-8s ") % (
					HMMResults.seqName,
					int(unit.seqFrom),
					int(unit.seqTo),
					int(unit.envFrom),
					int(unit.envTo),
					scanData._accmap[ unit.name ],
					unit.name,
					scanData._type[ unit.name ],
					int(unit.hmmFrom),
					int(unit.hmmTo),
					int(scanData._model_len[ unit.name ]),
					unit.bits,
					unit.evalue,
					1, #int(unit.sig), 
					clan )
					)
				#if (unit->{'act_site'}):
				#	outStr += predicted_active_site[@{$unit->{'act_site'}}]";
				
				#if($scanData->{_align}){
				#    print $fh sprintf( "%-10s %s\n", "#HMM",   $unit->hmmalign->{hmm} );
				#    print $fh sprintf( "%-10s %s\n", "#MATCH", $unit->hmmalign->{match} );
				#    print $fh sprintf( "%-10s %s\n", "#PP",   $unit->hmmalign->{pp});
				#    print $fh sprintf( "%-10s %s\n", "#SEQ",   $unit->hmmalign->{seq});

		return "\n".join( out_array )
		
class HMMSequence(HMMMatch):
	def __init__(self, **kwargs):
		HMMMatch.__init__(self)
		self.sumEvalue = None
		self.H2mode = None
		self.sumScore = None
		self.desc = None
		self.numberHits = None
		self.exp = None
		self.hmmUnits = []		
		for a in kwargs:
			if hasattr( self, a ):
				setattr( self, a, kwargs[a] )
			else:
				raise TypeError( "__init__() got an unexpected keyword argument '%s'" % (a) )
				

	def addHMMUnit ( self, hmmUnit ):
		self.hmmUtils.append( hmmUnit )

class HMMUnit(HMMMatch):
	"""
	Describes one HMMER alignment
	includes:
	self.Domain : Is Domain?
	self.name : Name of HMM hit
	self.proteinCoos = 
	self.seqEvalue = sequence Evalue (float)
	self.domain = 
	self.seqFrom = 
	self.seqTo = 
	self.domEvalue = 
	self.hmmalign = {}
	self.hmmFrom = 
	self.hmmTo = 
	self.envFrom =
	self.envTo =
	self.coreFrom =
	self.coreTo = 
	self.aliAcc = 
	self.sig = 
	"""
	def __init__(self, **kwargs):
		HMMMatch.__init__(self)
		self.seqName = None
		self.Domain = None
		self.proteinCoos = None
		self.seqEvalue = None
		self.domain = None
		self.seqFrom = None
		self.seqTo = None
		self.domEvalue = None
		self.hmmalign = {}
		self.hmmFrom = None
		self.hmmTo = None
		self.envFrom = None
		self.envTo = None
		self.coreFrom = None
		self.coreTo = None
		self.aliAcc = None
		self.sig = None
		for a in kwargs:
			if hasattr( self, a ):
				setattr( self, a, kwargs[a] )
			else:
				raise TypeError( "__init__() got an unexpected keyword argument '%s'" % (a) )
	
	def __str__(self):
		return ("%s %6d %6d %6d %6d %-10s %-16s %7s %5d %5d %5d %8s %9s %3d %-8s ") % (
					self.seqName,
					int(self.seqFrom),
					int(self.seqTo),
					int(self.envFrom),
					int(self.envTo),
					getattr(self, "acc", "N/A"),
					self.name,
					getattr(self, "type", "N/A"),
					int(self.hmmFrom),
					int(self.hmmTo),
					int(getattr(self, "model_len", "-1")), 
					self.bits,
					self.evalue,
					1, #int(unit.sig), 
					"" )

