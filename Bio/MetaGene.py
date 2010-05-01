
"""MetaGene is a parser written to extra Genes from metagenomic 
sequence samples based on predictions from MetaGeneAnnotator.
The MetaGeneAnnotator can be found at http://metagene.cb.k.u-tokyo.ac.jp/metagene/

Usage:
	>>> from Bio.MetaGene import MetaGeneParser
	>>> mgHandle = open( mgOutPath )  #Path to the metagene output file
	>>> seqHandle = open( fastaPath ) #path to the fasta file that MetaGeneAnnotator was run on
	>>> seq = Bio.SeqIO.parse( seqHandle, "fasta" )
	>>> myMG = Bio.MetaGene.MetaGeneParser()
	>>> mgList = myMG.parse( seq, mgHandle )

"""

import sys
import os
import re
import string
import Bio.SeqRecord
from Bio.Alphabet import IUPAC

re_comment = re.compile(r'^# ')
re_entry   = re.compile(r'^[0-9]+\s+[0-9]+')
re_space   = re.compile(r'\s+')

class MetaGeneData:
	"""Holds metagene prediction data, such as:
start
stop
strand
frame
state ( does the gene fall off the end of the DNA fragment )
geneScore
model ( self, bacterial, archaea, phage )
"""
	modelCode = { 's' : 'self', 'b' : 'bacterial', 'a' : 'archaea', 'p' : 'phage' }

	def __init__(self):
		pass

	def parseLine(self, line):	
		tmp = re_space.split( line )
		self.start  = int(tmp[1])-1
		self.stop   = int(tmp[2])
		self.strand = tmp[3]
		self.frame  = int(tmp[4])
		self.state   = tmp[5]
		self.geneScore = float(tmp[6])
		self.model     = tmp[7]

	def hasStart(self):
		"""Returns True if predicted gene has start codon contained on the 
		sequence fragment, False if the start of the sequence falls off of the fragment."""
		if self.state[0] == '1':
			return True
		return False

	def hasStop(self):
		"""Returns True if predicted gene has stop codon contained on the 
		sequence fragment, False if the end of the sequence falls off of the fragment."""
		if self.state[1] == '1':
			return True
		return False

	def getModel(self):
		"""Which model was used to predict the gene: self, bacterial, archaea, phage """
		try:
			return MetaGeneData.modelCode[ self.model ]
		except KeyError:
			return None

	def __str__(self):
		return "%s:%s%s" % ( self.start, self.stop, self.strand) 

class MetaGene(Bio.SeqRecord.SeqRecord):
	"""Extension of the Bio.SeqRecord.SeqRecord object.  
	Contains the the amino acid translation of the MetaGeneAnnotator prediction.  
	Additional MetaGene prediction information is stored in the self.metaGene member"""
	def __init__(self, parentDNA, mgData):
		newDNA = None
		if (mgData.strand == '+'):
			newDNA = parentDNA.seq[ mgData.start+mgData.frame:mgData.stop ]
		if (mgData.strand == '-'):
			newDNA = parentDNA.seq[ mgData.start:mgData.stop-mgData.frame ].reverse_complement()
		assert newDNA
		#print newDNA
		Bio.SeqRecord.SeqRecord.__init__( self, Bio.Seq.Seq( str(newDNA) ), id="%s|%s" % (parentDNA.id, str(mgData) ),
			description="Metagene %s=%s" % (mgData.getModel(), mgData.geneScore)  )
		self.metaGene = mgData

class MetaGeneRunner:
	"""Run MetaGeneAnnotator as a subprocess and translate the predicted genes into amino acid sequences"""
	def __init__(self, metagene_path = "mga_linux_ia32" ):
		"""metagene_path : The path to the MetaGeneAnnotator binary"""
		self.metagene_path = metagene_path

	def run(self, nucleotide_file ):
		(child_stdin, child_stdout) = os.popen2( "%s %s" % (self.metagene_path, nucleotide_file ) )
		child_stdin.close()
		self.parse( child_stdout )
		child_stdout.close()
		my_seq = seq.fasta( nucleotide_file )
		self.set_fasta(my_seq.read())	

class MetaGeneParser:	
	def parse( self, seqRecord, mgHandle ):
		"""Creates generator that returns the amino acid translations of MetaGeneAnnotator predictions.
		seqRecord : A SeqRecord generator created by parsing the same fasta file that was sent to MetaGeneAnnotator
		mgHandle  : A file handle to the output of MetaGeneAnnotator"""
		state = 0
		cur_i = 0
		cur_entry = -1
		curSeq = seqRecord.next()
		for line in mgHandle:
			if ( re_comment.search( line ) ):
				if cur_i%3 == 0:
					tmp = re_space.split( line )
					curID = tmp[1]
					while curID != curSeq.id:
						curSeq = seqRecord.next()
						assert curSeq is not None
					assert( curID == curSeq.id )
				cur_i += 1
			else:
				mgData = MetaGeneData( )
				mgData.parseLine( line )
				gene = MetaGene( curSeq, mgData )	
				yield gene



