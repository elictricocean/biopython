
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Interfaces import AlignmentIterator
from Bio.Align import MultipleSeqAlignment
import re

class EOF(Exception):
    pass
    
class HMMER3Iterator(AlignmentIterator):
    def next(self):
        ids = None
        seqs = None
        
        try:
            if len( self.outList ) > 0:
                return self.outList.pop()
        except AttributeError:
            self.outList = []
        
        while 1:
            try:
                headInfo = {}
                _readHeader( self.handle, headInfo )
                seqHits = []
                _readSeqHits( self.handle, seqHits )
                unitHits = []
                _readUnitHits( self.handle, headInfo, unitHits )
                #print headInfo
                #print seqHits            
                #print unitHits
                records = [] #Alignment obj will put them all in a list anyway
                for i in range(len(seqHits)):
                    #print seqHits[i]
                    #print unitHits[i]
                    seq = SeqRecord(Seq(unitHits[i]['seq'], self.alphabet),
                                   id = headInfo['seqName'], description = headInfo['description'],
                                   annotations = {})
                    hmm = SeqRecord(Seq(unitHits[i]['hmm'], self.alphabet),
                                   id = seqHits[i]['name'],
                                   annotations = {})
                    alignment = MultipleSeqAlignment( [seq, hmm], self.alphabet)
                    self.outList.append( alignment )
                return self.outList.pop()
            except EOF:
                raise StopIteration


def _readHeader( hs, hmmRes ):
    while ( 1 ):
        line = hs.readline()
        if not line:
            raise EOF
        if ( line.startswith( "Scores for complete" ) ):
            return                
        res = re.search(r'^# query HMM file:\s+(\S+)', line)
        if ( res ):
            hmmRes['hmmName'] = res.group(1)
            continue                
        res = re.search(r'^# target sequence database:\s+(\S+)', line)
        if ( res ):
            hmmRes['seqDB'] = res.group(1)
            continue                
        res = re.search( r'^output directed to file:\s+(\S+)', line)
        if res:
            hmmRes['thisFile'] = res.group(1)
            continue
        res = re.search( r'^Query:\s+(\S+)\s+\[M\=(\d+)\]', line) 
        if res:
            hmmRes['seedName'] = res.group(1)
            hmmRes['hmmLength'] = res.group(2)
            continue
        res = re.search( r'^Query:\s+(\S+)\s+\[L\=(\d+)\]', line )
        if res:
            hmmRes['seqName'] = res.group(1)
            hmmRes['seqLength'] = res.group(2)
            continue                
        res = re.search(r'^sequence E-value threshold: <= (\d+)', line)
        if res:
            hmmRes['evalueThr'] = float(res.group(1))
            continue
        res = re.search(r'^# Random generator seed:      (\d+)', line)
        if res:
            hmmRes['randSeedNum'] = res.group(1)
            continue
        res = re.search(r'^Description:\s+(.*)', line)
        if res:
            hmmRes['description'] = res.group(1)
            continue
        res = re.search(r'^# (phmmer|hmmsearch|hmmscan)', line)
        if res:
            hmmRes['program'] = res.group(1)
            continue
        #raise FormatError( "Failed to parse %s in sequence section\n" % (line) )


def _readSeqHits( hs, hmmRes ):
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
        hmmRes.append(
            {
                    'evalue'     : float(sMatch[1]),
                    'bits'       : float(sMatch[2]),
                    'bias'       : float(sMatch[3]),
                    'exp'        : float(sMatch[7]),
                    'numberHits' : int(sMatch[8]),
                    'name'       : sMatch[9],
                    'desc'       : desc         
            }
        )

def _readUnitHits( hs, headInfo, hmmRes ):
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
            _readUnitData( headInfo, seqId, hs, hmmRes )
            if headInfo['eof']:
                break
            continue
            

def _readUnitData(headInfo, id, hs, hmmRes ):
    #hmmName = hmmRes.seedName
    #seqName = hmmRes.seqName
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
            dMatch = re.split( r'\s+', line.rstrip() )
            if len(dMatch) != 17:
                raise FormatError( "Expected 16 elements of datam, got %d: %s " % (len(dMatch), line) )
            units.append(
                {
                    #'seqName'   : seqName,
                    'name'      : id,
                    'domain'    : dMatch[1],
                    'bits'      : float(dMatch[3]),
                    'bias'      : dMatch[4],
                    'domEvalue' : float(dMatch[5]),
                    'evalue'    : float(dMatch[6]),
                    'hmmFrom'   : int(dMatch[7]),
                    'hmmTo'     : int(dMatch[8]),
                    'seqFrom'   : int(dMatch[10]),
                    'seqTo'     : int(dMatch[11]),
                    'envFrom'   : int(dMatch[13]),
                    'envTo'     : int(dMatch[14]),
                    'aliAcc'    : dMatch[16]
                }
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
        if ( headInfo.has_key('hmmName') and headInfo['program'] == 'hmmsearch'):
            pattern1 = re.compile(r'^\s+%s\s+\d+\s+(\S+)\s+\d+' % (hmmName) )
            pattern2 = re.compile(r'^\s+%s\s+\d+\s+(\S+)\s+\d+' % (id) )
        if ( headInfo.has_key('seqName') and headInfo['program'] == 'hmmscan'):
            tmpSeqName = headInfo['seqName']
            tmpSeqName = re.sub(r'\|', r'\\|', tmpSeqName)
            pattern1 = re.compile(r'^\s+%s\s+\d+\s+(\S+)\s+\d+' % (id))
            pattern2 = re.compile(r'^\s+%s\s+\d+\s+(\S+)\s+\d+' % (tmpSeqName) )
        if ( headInfo.has_key('seqName') and headInfo['program'] == 'phmmer'):
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
                    units[ matchNo - 1 ]['hmm'] += res.group( 1 )
                except KeyError:
                    units[ matchNo - 1 ]['hmm'] = res.group( 1 )
                hmmlen = len( res.group(1) )
                continue
            res = pattern2.search( line )
            if res:
                try:
                    units[ matchNo - 1 ]['seq'] += res.group( 1 )
                except KeyError:
                    units[ matchNo - 1 ]['seq'] = res.group( 1 )
                continue
            res = re.search(r'^\s+([0-9\*\.]+)\s+PP$', line) 
            if res:
                pp = res.group(1)
                try:
                    units[ matchNo - 1 ]['pp'] += pp
                except KeyError:
                    units[ matchNo - 1 ]['pp'] = pp                        
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
                    units[ matchNo - 1 ]['match'] += m1
                except KeyError:
                    units[ matchNo - 1 ]['match'] = m1                        
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

    headInfo['eof'] = eof
    for u in units:
        hmmRes.append(u)

    if (recurse and nextSeqId):
        self._readUnitData( nextSeqId, hs, hmmRes )

   