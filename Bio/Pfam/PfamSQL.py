

import re

re_quoteSplit = re.compile(r'\'\s*\'')
re_quoteMiddle = re.compile(r'^\s*\'(.*)\'\s*$')
re_newlinecode = re.compile(r'\\[Nn]')

class GOEntry:
	def __init__(self, **kwargs):
		self.auto_pfamA = None
		self.go_id = None
		self.term = None
		self.category = None
		for a in kwargs:
			if hasattr( self, a ):
				setattr(self, a, kwargs[a])
			else:
				raise ValueError("Unknown arg %s" % (a))

	def id(self):
		return self.go_id.split(":")[1]

	def __str__(self):
		return "%s=%s" % (self.go_id, self.term)

class GOTable:
	"""Parse the SQL flatfile found at 
ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam24.0/database_files/gene_ontology.txt.gz
"""
	def __init__(self):
		self.autoPfamHash = {}
	
	def read(self, handle):
		for line in handle:
			tmp = re_quoteSplit.split( re_quoteMiddle.sub( r'\1', line) )
			#col order: auto_pfamA, go_id, term, category
			entry = GOEntry( auto_pfamA=tmp[0], go_id=tmp[1], term=tmp[2], category=tmp[3] )
			try:
				self.autoPfamHash[ tmp[0] ].append( entry )
			except KeyError:
				self.autoPfamHash[ tmp[0] ] = [ entry ] 

	def entry(self, id):
		try:
			if isinstance( id, PfamAEntry):
				return self.autoPfamHash[ id.auto_pfamA ]
			else:
				return self.autoPfamHash[ id ]
		except KeyError:
			return []
			
class PfamAEntry:
	def __init__(self, **kwargs):
		self.auto_pfamA = None
		self.pfamA_acc  = None
		self.type = None
		self.comment = None
		for a in kwargs:
			if hasattr( self, a ):
				setattr(self, a, kwargs[a])
			else:
				raise ValueError("Unknown arg %s" % (a))		

	def __str__(self):
		return "%s %s %s" % ( self.pfamA_acc, self.type, self.comment )


class PfamATable:
	"""Parse the SQL flat file found at:
ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam24.0/database_files/pfamA.txt.gz
"""
	def __init__(self):
		self.autoPfamHash = {}
		self.accPfamHash = {}
		
	def read(self, handle):
		for line in handle:
			tmp = re_quoteSplit.split( re_quoteMiddle.sub( r'\1', re_newlinecode.sub("",line) ) )
			#col order:
			#0 auto_pfamA
			#1 pfamA_acc
			#2 pfamA_id
			#3 previous_id
			#4 description
			#5 author
			#6 deposited_by
			#7 seed_source
			#8 type
			#9 comment
			#10 sequence_GA
			#11 domain_GA
			#12 sequence_TC
			#13 domain_TC
			#14 sequence_NC
			#15 domain_NC
			#16 buildMethod
			#17 model_length
			#18 searchMethod
			#19 msv_lambda
			#20 msv_mu
			#21 viterbi_lambda
			#22 viterbi_mu
			#23 forward_lambda
			#24 forward_tau
			#25 num_seed
			#26 num_full
			#27 updated
			#28 created
			#29 version
			#30 number_archs
			#31 number_species
			#32 number_structures
			#33 number_ncbi
			#34 number_meta
			#35 average_length
			#36 percentage_id
			#37 average_coverage
			#38 change_status
			#39 seed_consensus
			#40 full_consensus

			entry = PfamAEntry( auto_pfamA=tmp[0], pfamA_acc=tmp[1], type=tmp[8], comment=tmp[9] )
			self.autoPfamHash[ tmp[0] ] = entry
			self.accPfamHash[ tmp[1] ] = entry


	def keys(self):
		return self.accPfamHash.keys()
		
	def entry(self, str):
		return self.accPfamHash[ str.split(".")[0] ]

if __name__ == '__main__':
	import sys
	dir = sys.argv[1]
	
		
	pfamATable = PfamATable()
	file = open( "%s/pfamA.txt" % (dir))
	pfamATable.read( file )
	file.close()

	goTable = GOTable()
	file = open( "%s/gene_ontology.txt" % (dir) )
	goTable.read( file )
	file.close()

	for a in pfamATable.keys():
		entry = pfamATable.entry(a)
		print entry.pfamA_acc, "GO:", ",".join( [ b.id() for b in goTable.entry( entry ) ] )

