# Copyright 2010 by Kyle Ellrott.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parser to interpret the HMMER 3 output file

This parser returns a Biopython object representation of the
data given by HMMER 3.0. 
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Interfaces import AlignmentIterator
from Bio.Align import MultipleSeqAlignment
import re
import string

class EOF(Exception):
    pass


class ParseError(Exception):
    def __init__(self, line):
        self.line = line
    def __str__(self):
        return "UNKNOWN LINE:%s" % (self.line)


reEnd = re.compile(r'^//$')
reProgram = re.compile(r'\# (hmmscan|hmmsearch) ::')
reQueryDef = re.compile(r'^Query:\s+(\S+)\s+\[[ML]=(\d+)\]')#(r'^Query:\s+(\w+)\s+\[[ML]=(266)\]')
reDescDef  = re.compile('^Description:\s*(.*)')
reScores = re.compile(r'Scores for complete')
reScoreFormat = re.compile( r'^\s+(E-value|---)' )
reHitScoreFormat = re.compile( r'^\s+(#\s+score|---)' )
reSpace = re.compile(r'\s+')
reDomains = re.compile(r'^Domain.*annotation.*:')
reAlignDesc = re.compile( r'^\>\>\s+(\S+)' )
reAlignHead = re.compile( r'^\s+Alignments for each domain:' )
reDomainAlignHead = re.compile( r'\s+== domain (\d+)' )
reTailStart = re.compile(r'^(Internal pipeline statistics summary|-------------------------------------)')
reTailInfoLine = re.compile(r'^(.*):\s+(\d+)')
reNoHits = re.compile(r'^\s+\[No (hits|targets) detected')
class HMMER3Iterator(AlignmentIterator):
    def next(self):
        if not hasattr( self, 'state' ):
            self.state = { 'eof' : False }
        ids = None
        seqs = None
        
        try:
            if len( self.outList ) > 0:
                return self.outList.pop()
        except AttributeError:
            self.outList = []
        
        while 1:
            try:
                self._parse( self.handle )
                self._prepResults()
                if len( self.outList ) > 0:
                    return self.outList.pop()               
            except EOF:
                self._prepResults()
                if len( self.outList ) > 0:
                    return self.outList.pop()
                raise StopIteration
        
    def _prepResults(self):
        #print "Query", self.state['queryName'], self.state['queryLen']
        #print "HitINFO:", self.hitInfo
        #print "HitRecord", self.hitRecord
        #print "HitAlign:", self.alignMap        
        if ( len(self.hitInfo) == 0 and not self.state['eof'] ):
            alignment = MultipleSeqAlignment( [], self.alphabet)
            self.outList.append( alignment )
        
        for hit in self.hitInfo:
            for domain in self.hitRecord[ hit ]:
                queryStr  = "".join( self.alignMap[ hit ][ domain ][ 'query' ] )
                targetStr = "".join( self.alignMap[ hit ][ domain ][ 'target' ] )
                
                query = SeqRecord(Seq(queryStr, self.alphabet),
                    id = self.state['queryName'], description = self.state.get( 'desc', "" ),
                    annotations = {})
                    
                target = SeqRecord(Seq(targetStr, self.alphabet),
                    id = hit,
                    annotations = {})
                alignment = HMMERAlign( [query,target], self.alphabet)
                alignment._annotations = self.hitRecord[ hit ][ domain ]
                alignment._annotations[ 'seqName' ] = self.state['queryName']
                alignment._annotations[ 'hmmName' ] = hit
                
                
                self.outList.append( alignment )
    
    def _parse( self, handle ):
        self.state[ 'fileLoc' ] = "start"
        self.hitInfo = {}
        self.hitRecord = {}
        self.alignMap = {}
        self.tailInfo = {}
        for line in handle:
            if reEnd.search( line ):
                return
            
            res = reProgram.search( line ) 
            if res:
                self.state['program'] = res.group(1)
                self.state['eof'] = False
                continue
            res = reQueryDef.search( line )
            if res:
                self.state['queryName'] = res.group(1)
                self.state['queryLen'] = int(res.group(2))
                continue
            res = reDescDef.search( line )
            if res:
                self.state['desc'] = res.group(1)
                continue
            if (line.startswith('#')):
                continue
            if (line=="\n"):
                continue
            res = reScores.search(line)
            if res:
                self.state['fileLoc'] = "scores"
                continue
            res = reScoreFormat.search( line )
            if res:
                continue
            res = reDomains.search( line )
            if res:
                continue
            
            res = reAlignDesc.search(line)
            if res:
                self.state[ 'curAlign' ] = res.group(1)
                self.alignMap[ self.state[ 'curAlign' ] ] = {}
                continue
            
            res = reHitScoreFormat.search( line )
            if res:
                self.state['fileLoc'] = 'hitScore'
                continue
            
            res = reAlignHead.search( line )
            if res:
                #self.state['fileLoc'] = 'hitAlign'
                continue
            
            res = reTailStart.search( line )
            if res:
                self.state['fileLoc'] = 'tail'
                continue
            res = reNoHits.search( line )
            if res:
                continue
            
            if self.state[ 'fileLoc' ] == "scores":
                tmp = reSpace.split( line )
                if ( len(tmp) > 10 ):
                    desc = " ".join( tmp[ 10: ] )
                    self.hitInfo[ tmp[9] ] = {
                        'name'       : tmp[9],
                        'evalue'     : float(tmp[1]),
                        'bits'       : float(tmp[2]),
                        'bias'       : float(tmp[3]),
                        'exp'        : float(tmp[7]),
                        'numberHits' : int(tmp[8])
                    }
                    self.hitRecord[ tmp[9] ] = {}
                    continue
            if self.state[ 'fileLoc' ] == 'hitScore':
                tmp = reSpace.split( line.rstrip() )
                if len(tmp) == 17:
                    self.hitRecord[ self.state[ 'curAlign' ] ][ tmp[1] ] = {
                        #seqName   = seqName,
                        #name      = id,
                        'bits'      : float(tmp[3]),
                        'bias'      : tmp[4],
                        'domEvalue' : float(tmp[5]),
                        'evalue'    : float(tmp[6]),
                        'hmmFrom'   : int(tmp[7]),
                        'hmmTo'     : int(tmp[8]),
                        'seqFrom'   : int(tmp[10]),
                        'seqTo'     : int(tmp[11]),
                        'envFrom'   : int(tmp[13]),
                        'envTo'     : int(tmp[14]),
                        'aliAcc'    : tmp[16]                    
                    }
                    continue
            
            res = reDomainAlignHead.search( line )
            if res:
                self.state[ 'fileLoc' ] = 'domainAlign_1'
                self.state[ 'curDomain' ] = res.group(1)
                self.alignMap[ self.state[ 'curAlign' ] ][ self.state[ 'curDomain' ] ] = { 'query' : [], 'target' : [], 'qRange' : [], 'tRange' : [] }
                continue

            if self.state[ 'fileLoc' ].startswith('domainAlign_'):
                if line.rstrip().endswith(' CS'):
                    continue
                if line.rstrip().endswith(' PP'):
                    continue
                if line.rstrip().endswith(' RF'):
                    continue
            
            if self.state[ 'fileLoc' ] == 'domainAlign_1':
                tmp = reSpace.split( line )
                range = ( int(tmp[2]), int(tmp[4]) )
                self.alignMap[ self.state[ 'curAlign' ] ][ self.state[ 'curDomain' ] ][ 'query' ].append( tmp[3] )
                self.alignMap[ self.state[ 'curAlign' ] ][ self.state[ 'curDomain' ] ][ 'qRange' ].append( range )
                self.state[ 'fileLoc' ] = 'domainAlign_2'
                continue
            
            if self.state[ 'fileLoc' ] == 'domainAlign_2':
                self.state[ 'fileLoc' ] = 'domainAlign_3'
                continue
            
            if self.state[ 'fileLoc' ] == 'domainAlign_3':
                tmp = reSpace.split( line )
                range = ( int(tmp[2]), int(tmp[4]) )
                self.alignMap[ self.state[ 'curAlign' ] ][ self.state[ 'curDomain' ] ][ 'target' ].append( tmp[3] )
                self.alignMap[ self.state[ 'curAlign' ] ][ self.state[ 'curDomain' ] ][ 'tRange' ].append( range )
                self.state[ 'fileLoc' ] = 'domainAlign_1'
                continue
            
            if self.state[ 'fileLoc' ] == 'tail':
                res = reTailInfoLine.search(line)
                if res:
                    self.tailInfo[ res.group(1) ] = res.group(2)
                    continue
            raise ParseError( line )
        self.state[ 'eof' ] = True
        raise EOF
    
    
    

class HMMERAlign(MultipleSeqAlignment):    
    def __str__( self ):
        outStr = "%s %6d %6d %6d %6d %-16s %5d %5d %5d %8.1f %10.3e" % (
                    self._annotations['seqName'],
                    int(self._annotations['seqFrom']),
                    int(self._annotations['seqTo']),
                    int(self._annotations['envFrom']),
                    int(self._annotations['envTo']),
                    self._annotations['hmmName'],
                    int(self._annotations['hmmFrom']),
                    int(self._annotations['hmmTo']),
                    int(getattr(self, "model_len", "-1")), 
                    self._annotations['bits'],
                    self._annotations['evalue']
                    )
        return outStr

    
