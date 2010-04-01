
import re

reEnd = re.compile(r'^//')
reGF  = re.compile(r'^\#=GF (..)\s+(.*)$')

def parse(handle, format, **kwargs):
    if format in ("stockholmInfo"):
        return parseStockholmInfo( handle, **kwargs )


class StockholmEntry:
    def __init__(self, dataSet):
        self.id = dataSet[ 'ID' ][0]
        self.comment = " ".join( dataSet.get( 'CC', [] ) )
        self.accession = dataSet.get( 'AC', [None] )[0]
        self.definition = dataSet.get( 'DE', [None] )[0]
        self.previous  = dataSet.get('PI', [None] )[0]

def parseStockholmInfo( handle, **kwargs ):
    out = {}
    curSet = {}
    for line in handle:
        if line.startswith( "#=GF " ):
            res = reGF.search( line )
            [type, data] = res.groups()
            try:
                curSet[ type ].append( data )
            except KeyError:
                curSet[ type ] = [ data ]
        elif line.startswith('//'):
            yield StockholmEntry( curSet )
            curSet = {}

