# Copyright 2000 by Jeffrey Chang, Brad Chapman.  All rights reserved.
# Copyright 2006-2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to work with GenBank formatted files.

Rather than using Bio.GenBank, you are now encouraged to use Bio.SeqIO with
the "genbank" or "embl" format names to parse GenBank or EMBL files into
SeqRecord and SeqFeature objects (see the Biopython tutorial for details).

Also, rather than using Bio.GenBank to search or download files from the NCBI,
you are now encouraged to use Bio.Entrez instead (again, see the Biopython
tutorial for details).

Currently the ONLY reason to use Bio.GenBank directly is for the RecordParser
which turns a GenBank file into GenBank-specific Record objects.  This is a
much closer representation to the raw file contents that the SeqRecord
alternative from the FeatureParser (used in Bio.SeqIO).

Classes:
Iterator              Iterate through a file of GenBank entries
ErrorFeatureParser    Catch errors caused during parsing.
FeatureParser         Parse GenBank data in SeqRecord and SeqFeature objects.
RecordParser          Parse GenBank data into a Record object.

Exceptions:
ParserFailureError    Exception indicating a failure in the parser (ie.
                      scanner or consumer)
LocationParserError   Exception indiciating a problem with the spark based
                      location parser.


17-MAR-2009: added wgs, wgs_scafld for GenBank whole genome shotgun master records.
These are GenBank files that summarize the content of a project, and provide lists of
scaffold and contig files in the project. These will be in annotations['wgs'] and
annotations['wgs_scafld']. These GenBank files do not have sequences. See
http://groups.google.com/group/bionet.molbio.genbank/browse_thread/thread/51fb88bf39e7dc36

http://is.gd/nNgk
for more details of this format, and an example.
Added by Ying Huang & Iddo Friedberg
"""
import cStringIO

# other Biopython stuff
from Bio import SeqFeature
from Bio.ParserSupport import AbstractConsumer
from Bio import Entrez

# other Bio.GenBank stuff
import LocationParser
from utils import FeatureValueCleaner
from Scanner import GenBankScanner

#Constants used to parse GenBank header lines
GENBANK_INDENT = 12
GENBANK_SPACER = " " * GENBANK_INDENT

#Constants for parsing GenBank feature lines
FEATURE_KEY_INDENT = 5
FEATURE_QUALIFIER_INDENT = 21
FEATURE_KEY_SPACER = " " * FEATURE_KEY_INDENT
FEATURE_QUALIFIER_SPACER = " " * FEATURE_QUALIFIER_INDENT

class Iterator:
    """Iterator interface to move over a file of GenBank entries one at a time.
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle with GenBank entries to iterate through.
        o parser - An optional parser to pass the entries through before
        returning them. If None, then the raw entry will be returned.
        """
        self.handle = handle
        self._parser = parser

    def next(self):
        """Return the next GenBank record from the handle.

        Will return None if we ran out of records.
        """
        if self._parser is None:
            lines = []
            while True:
                line = self.handle.readline()
                if not line : return None #Premature end of file?
                lines.append(line)
                if line.rstrip() == "//" : break
            return "".join(lines)
        try:
            return self._parser.parse(self.handle)
        except StopIteration:
            return None

    def __iter__(self):
        return iter(self.next, None)

class ParserFailureError(Exception):
    """Failure caused by some kind of problem in the parser.
    """
    pass

class LocationParserError(Exception):
    """Could not Properly parse out a location from a GenBank file.
    """
    pass
                                                          
class FeatureParser:
    """Parse GenBank files into Seq + Feature objects.
    """
    def __init__(self, debug_level = 0, use_fuzziness = 1, 
                 feature_cleaner = FeatureValueCleaner()):
        """Initialize a GenBank parser and Feature consumer.

        Arguments:
        o debug_level - An optional argument that species the amount of
        debugging information the parser should spit out. By default we have
        no debugging info (the fastest way to do things), but if you want
        you can set this as high as two and see exactly where a parse fails.
        o use_fuzziness - Specify whether or not to use fuzzy representations.
        The default is 1 (use fuzziness).
        o feature_cleaner - A class which will be used to clean out the
        values of features. This class must implement the function 
        clean_value. GenBank.utils has a "standard" cleaner class, which
        is used by default.
        """
        self._scanner = GenBankScanner(debug_level)
        self.use_fuzziness = use_fuzziness
        self._cleaner = feature_cleaner

    def parse(self, handle):
        """Parse the specified handle.
        """
        self._consumer = _FeatureConsumer(self.use_fuzziness, 
                                          self._cleaner)
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class RecordParser:
    """Parse GenBank files into Record objects
    """
    def __init__(self, debug_level = 0):
        """Initialize the parser.

        Arguments:
        o debug_level - An optional argument that species the amount of
        debugging information the parser should spit out. By default we have
        no debugging info (the fastest way to do things), but if you want
        you can set this as high as two and see exactly where a parse fails.
        """
        self._scanner = GenBankScanner(debug_level)

    def parse(self, handle):
        """Parse the specified handle into a GenBank record.
        """
        self._consumer = _RecordConsumer()
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class _BaseGenBankConsumer(AbstractConsumer):
    """Abstract GenBank consumer providing useful general functions.

    This just helps to eliminate some duplication in things that most
    GenBank consumers want to do.
    """
    # Special keys in GenBank records that we should remove spaces from
    # For instance, \translation keys have values which are proteins and
    # should have spaces and newlines removed from them. This class
    # attribute gives us more control over specific formatting problems.
    remove_space_keys = ["translation"]

    def __init__(self):
        pass

    def _split_keywords(self, keyword_string):
        """Split a string of keywords into a nice clean list.
        """
        # process the keywords into a python list
        if keyword_string == "" or keyword_string == ".":
            keywords = ""
        elif keyword_string[-1] == '.':
            keywords = keyword_string[:-1]
        else:
            keywords = keyword_string
        keyword_list = keywords.split(';')
        clean_keyword_list = [x.strip() for x in keyword_list]
        return clean_keyword_list

    def _split_accessions(self, accession_string):
        """Split a string of accession numbers into a list.
        """
        # first replace all line feeds with spaces
        # Also, EMBL style accessions are split with ';'
        accession = accession_string.replace("\n", " ").replace(";"," ")

        return [x.strip() for x in accession.split() if x.strip()]

    def _split_taxonomy(self, taxonomy_string):
        """Split a string with taxonomy info into a list.
        """
        if not taxonomy_string or taxonomy_string==".":
            #Missing data, no taxonomy
            return []
        
        if taxonomy_string[-1] == '.':
            tax_info = taxonomy_string[:-1]
        else:
            tax_info = taxonomy_string
        tax_list = tax_info.split(';')
        new_tax_list = []
        for tax_item in tax_list:
            new_items = tax_item.split("\n")
            new_tax_list.extend(new_items)
        while '' in new_tax_list:
            new_tax_list.remove('')
        clean_tax_list = [x.strip() for x in new_tax_list]

        return clean_tax_list

    def _clean_location(self, location_string):
        """Clean whitespace out of a location string.

        The location parser isn't a fan of whitespace, so we clean it out
        before feeding it into the parser.
        """
        #Originally this imported string.whitespace and did a replace
        #via a loop.  It's simpler to just split on whitespace and rejoin
        #the string - and this avoids importing string too.  See Bug 2684.
        return ''.join(location_string.split())

    def _remove_newlines(self, text):
        """Remove any newlines in the passed text, returning the new string.
        """
        # get rid of newlines in the qualifier value
        newlines = ["\n", "\r"]
        for ws in newlines:
            text = text.replace(ws, "")

        return text

    def _normalize_spaces(self, text):
        """Replace multiple spaces in the passed text with single spaces.
        """
        # get rid of excessive spaces
        text_parts = text.split(" ")
        text_parts = filter(None, text_parts)
        return ' '.join(text_parts)

    def _remove_spaces(self, text):
        """Remove all spaces from the passed text.
        """
        return text.replace(" ", "")

    def _convert_to_python_numbers(self, start, end):
        """Convert a start and end range to python notation.

        In GenBank, starts and ends are defined in "biological" coordinates,
        where 1 is the first base and [i, j] means to include both i and j.

        In python, 0 is the first base and [i, j] means to include i, but
        not j. 

        So, to convert "biological" to python coordinates, we need to 
        subtract 1 from the start, and leave the end and things should
        be converted happily.
        """
        new_start = start - 1
        new_end = end

        return new_start, new_end

class _FeatureConsumer(_BaseGenBankConsumer):
    """Create a SeqRecord object with Features to return.

    Attributes:
    o use_fuzziness - specify whether or not to parse with fuzziness in
    feature locations.
    o feature_cleaner - a class that will be used to provide specialized
    cleaning-up of feature values.
    """
    def __init__(self, use_fuzziness, feature_cleaner = None):
        from Bio.SeqRecord import SeqRecord
        _BaseGenBankConsumer.__init__(self)
        self.data = SeqRecord(None, id = None)
        self.data.id = None
        self.data.description = ""

        self._use_fuzziness = use_fuzziness
        self._feature_cleaner = feature_cleaner

        self._seq_type = ''
        self._seq_data = []
        self._cur_reference = None
        self._cur_feature = None
        self._cur_qualifier_key = None
        self._cur_qualifier_value = None
        self._expected_size = None

    def locus(self, locus_name):
        """Set the locus name is set as the name of the Sequence.
        """
        self.data.name = locus_name

    def size(self, content):
        """Record the sequence length."""
        self._expected_size = int(content)

    def residue_type(self, type):
        """Record the sequence type so we can choose an appropriate alphabet.
        """
        self._seq_type = type

    def data_file_division(self, division):
        self.data.annotations['data_file_division'] = division

    def date(self, submit_date):
        self.data.annotations['date'] = submit_date 

    def definition(self, definition):
        """Set the definition as the description of the sequence.
        """
        if self.data.description:
            #Append to any existing description
            #e.g. EMBL files with two DE lines.
            self.data.description += " " + definition
        else:
            self.data.description = definition

    def accession(self, acc_num):
        """Set the accession number as the id of the sequence.

        If we have multiple accession numbers, the first one passed is
        used.
        """
        new_acc_nums = self._split_accessions(acc_num)

        #Also record them ALL in the annotations
        try:
            #On the off chance there was more than one accession line:
            for acc in new_acc_nums:
                #Prevent repeat entries
                if acc not in self.data.annotations['accessions']:
                    self.data.annotations['accessions'].append(acc)
        except KeyError:
            self.data.annotations['accessions'] = new_acc_nums

        # if we haven't set the id information yet, add the first acc num
        if self.data.id is None:
            if len(new_acc_nums) > 0:
                #self.data.id = new_acc_nums[0]
                #Use the FIRST accession as the ID, not the first on this line!
                self.data.id = self.data.annotations['accessions'][0]

    def wgs(self, content):
        self.data.annotations['wgs'] = content.split('-')

    def add_wgs_scafld(self, content):
        self.data.annotations.setdefault('wgs_scafld',[]).append(content.split('-'))

    def nid(self, content):
        self.data.annotations['nid'] = content

    def pid(self, content):
        self.data.annotations['pid'] = content

    def version(self, version_id):
        #Want to use the versioned accession as the record.id
        #This comes from the VERSION line in GenBank files, or the
        #obsolete SV line in EMBL.  For the new EMBL files we need
        #both the version suffix from the ID line and the accession
        #from the AC line.
        if version_id.count(".")==1 and version_id.split(".")[1].isdigit():
            self.accession(version_id.split(".")[0])
            self.version_suffix(version_id.split(".")[1])
        else:
            #For backwards compatibility...
            self.data.id = version_id

    def project(self, content):
        """Handle the information from the PROJECT line as a list of projects.

        e.g.
        PROJECT     GenomeProject:28471

        or:
        PROJECT     GenomeProject:13543  GenomeProject:99999

        This is stored as dbxrefs in the SeqRecord to be consistent with the
        projected switch of this line to DBLINK in future GenBank versions.
        Note the NCBI plan to replace "GenomeProject:28471" with the shorter
        "Project:28471" as part of this transition.
        """
        content = content.replace("GenomeProject:", "Project:")
        self.data.dbxrefs.extend([p for p in content.split() if p])

    def dblink(self, content):
        """Store DBLINK cross references as dbxrefs in our record object.

        This line type is expected to replace the PROJECT line in 2009. e.g.

        During transition:
        
        PROJECT     GenomeProject:28471
        DBLINK      Project:28471
                    Trace Assembly Archive:123456

        Once the project line is dropped:

        DBLINK      Project:28471
                    Trace Assembly Archive:123456

        Note GenomeProject -> Project.

        We'll have to see some real examples to be sure, but based on the
        above example we can expect one reference per line.
        """
        #During the transition period with both PROJECT and DBLINK lines,
        #we don't want to add the same cross reference twice.
        if content.strip() not in self.data.dbxrefs:
            self.data.dbxrefs.append(content.strip())

    def version_suffix(self, version):
        """Set the version to overwrite the id.

        Since the verison provides the same information as the accession
        number, plus some extra info, we set this as the id if we have
        a version.
        """
        #e.g. GenBank line:
        #VERSION     U49845.1  GI:1293613
        #or the obsolete EMBL line:
        #SV   U49845.1
        #Scanner calls consumer.version("U49845.1")
        #which then calls consumer.version_suffix(1)
        #
        #e.g. EMBL new line:
        #ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.
        #Scanner calls consumer.version_suffix(1)
        assert version.isdigit()
        self.data.annotations['sequence_version'] = int(version)

    def db_source(self, content):
        self.data.annotations['db_source'] = content.rstrip()

    def gi(self, content):
        self.data.annotations['gi'] = content

    def keywords(self, content):
        self.data.annotations['keywords'] = self._split_keywords(content)

    def segment(self, content):
        self.data.annotations['segment'] = content

    def source(self, content):
        #Note that some software (e.g. VectorNTI) may produce an empty
        #source (rather than using a dot/period as might be expected).
        if content == "":
            source_info = ""
        elif content[-1] == '.':
            source_info = content[:-1]
        else:
            source_info = content
        self.data.annotations['source'] = source_info

    def organism(self, content):
        self.data.annotations['organism'] = content

    def taxonomy(self, content):
        """Records (another line of) the taxonomy lineage.
        """
        lineage = self._split_taxonomy(content)
        try:
            self.data.annotations['taxonomy'].extend(lineage)
        except KeyError:
            self.data.annotations['taxonomy'] = lineage
        
    def reference_num(self, content):
        """Signal the beginning of a new reference object.
        """
        # if we have a current reference that hasn't been added to
        # the list of references, add it.
        if self._cur_reference is not None:
            self.data.annotations['references'].append(self._cur_reference)
        else:
            self.data.annotations['references'] = []

        self._cur_reference = SeqFeature.Reference()

    def reference_bases(self, content):
        """Attempt to determine the sequence region the reference entails.

        Possible types of information we may have to deal with:
        
        (bases 1 to 86436)
        (sites)
        (bases 1 to 105654; 110423 to 111122)
        1  (residues 1 to 182)
        """
        # first remove the parentheses or other junk
        ref_base_info = content[1:-1]

        all_locations = []
        # parse if we've got 'bases' and 'to'
        if ref_base_info.find('bases') != -1 and \
            ref_base_info.find('to') != -1:
            # get rid of the beginning 'bases'
            ref_base_info = ref_base_info[5:]
            locations = self._split_reference_locations(ref_base_info)
            all_locations.extend(locations)
        elif (ref_base_info.find("residues") >= 0 and
              ref_base_info.find("to") >= 0):
            residues_start = ref_base_info.find("residues")
            # get only the information after "residues"
            ref_base_info = ref_base_info[(residues_start + len("residues ")):]
            locations = self._split_reference_locations(ref_base_info)
            all_locations.extend(locations)

        # make sure if we are not finding information then we have
        # the string 'sites' or the string 'bases'
        elif (ref_base_info == 'sites' or
              ref_base_info.strip() == 'bases'):
            pass
        # otherwise raise an error
        else:
            raise ValueError("Could not parse base info %s in record %s" %
                             (ref_base_info, self.data.id))

        self._cur_reference.location = all_locations

    def _split_reference_locations(self, location_string):
        """Get reference locations out of a string of reference information
        
        The passed string should be of the form:

            1 to 20; 20 to 100

        This splits the information out and returns a list of location objects
        based on the reference locations.
        """
        # split possibly multiple locations using the ';'
        all_base_info = location_string.split(';')

        new_locations = []
        for base_info in all_base_info:
            start, end = base_info.split('to')
            new_start, new_end = \
              self._convert_to_python_numbers(int(start.strip()),
                                              int(end.strip()))
            this_location = SeqFeature.FeatureLocation(new_start, new_end)
            new_locations.append(this_location)
        return new_locations

    def authors(self, content):
        if self._cur_reference.authors:
            self._cur_reference.authors += ' ' + content
        else:
            self._cur_reference.authors = content

    def consrtm(self, content):
        if self._cur_reference.consrtm:
            self._cur_reference.consrtm += ' ' + content
        else:
            self._cur_reference.consrtm = content

    def title(self, content):
        if self._cur_reference is None:
            import warnings
            warnings.warn("GenBank TITLE line without REFERENCE line.")
        elif self._cur_reference.title:
            self._cur_reference.title += ' ' + content
        else:
            self._cur_reference.title = content

    def journal(self, content):
        if self._cur_reference.journal:
            self._cur_reference.journal += ' ' + content
        else:
            self._cur_reference.journal = content

    def medline_id(self, content):
        self._cur_reference.medline_id = content

    def pubmed_id(self, content):
        self._cur_reference.pubmed_id = content

    def remark(self, content):
        """Deal with a reference comment."""
        if self._cur_reference.comment:
            self._cur_reference.comment += ' ' + content
        else:
            self._cur_reference.comment = content

    def comment(self, content):
        try:
            self.data.annotations['comment'] += "\n" + "\n".join(content)
        except KeyError:
            self.data.annotations['comment'] = "\n".join(content)

    def features_line(self, content):
        """Get ready for the feature table when we reach the FEATURE line.
        """
        self.start_feature_table()

    def start_feature_table(self):
        """Indicate we've got to the start of the feature table.
        """
        # make sure we've added on our last reference object
        if self._cur_reference is not None:
            self.data.annotations['references'].append(self._cur_reference)
            self._cur_reference = None

    def _add_feature(self):
        """Utility function to add a feature to the SeqRecord.

        This does all of the appropriate checking to make sure we haven't
        left any info behind, and that we are only adding info if it
        exists.
        """
        if self._cur_feature:
            # if we have a left over qualifier, add it to the qualifiers
            # on the current feature
            self._add_qualifier()

            self._cur_qualifier_key = ''
            self._cur_qualifier_value = ''
            self.data.features.append(self._cur_feature)
            
    def feature_key(self, content):
        # if we already have a feature, add it on
        self._add_feature()

        # start a new feature
        self._cur_feature = SeqFeature.SeqFeature()
        self._cur_feature.type = content

        # assume positive strand to start with if we have DNA or cDNA
        # (labelled as mRNA). The complement in the location will 
        # change this later if something is on the reverse strand
        if self._seq_type.find("DNA") >= 0 or self._seq_type.find("mRNA") >= 0:
            self._cur_feature.strand = 1

    def location(self, content):
        """Parse out location information from the location string.

        This uses a comprehensive but slow spark based parser to do the
        parsing, and then translates the results of the parse into appropriate
        Location objects.
        """
        # --- first preprocess the location for the spark parser
        
        # we need to clean up newlines and other whitespace inside
        # the location before feeding it to the parser.
        # locations should have no whitespace whatsoever based on the
        # grammer
        location_line = self._clean_location(content)

        # Older records have junk like replace(266,"c") in the
        # location line. Newer records just replace this with
        # the number 266 and have the information in a more reasonable
        # place. So we'll just grab out the number and feed this to the
        # parser. We shouldn't really be losing any info this way.
        if location_line.find('replace') != -1:
            comma_pos = location_line.find(',')
            location_line = location_line[8:comma_pos]

        # feed everything into the scanner and parser
        try:
            parse_info = \
                       LocationParser.parse(LocationParser.scan(location_line))
        # spark raises SystemExit errors when parsing fails
        except SystemExit:
            raise LocationParserError(location_line)

        # print "parse_info:", repr(parse_info)
        
        # add the parser information the current feature
        self._set_location_info(parse_info, self._cur_feature)

    def _set_function(self, function, cur_feature):
        """Set the location information based on a function.

        This handles all of the location functions like 'join', 'complement'
        and 'order'.

        Arguments:
        o function - A LocationParser.Function object specifying the
        function we are acting on.
        o cur_feature - The feature to add information to.
        """
        assert isinstance(function, LocationParser.Function), \
               "Expected a Function object, got %s" % function
        
        if function.name == "complement":
            # mark the current feature as being on the opposite strand
            cur_feature.strand = -1
            # recursively deal with whatever is left inside the complement
            for inner_info in function.args:
                self._set_location_info(inner_info, cur_feature)
        # deal with functions that have multipe internal segments that
        # are connected somehow.
        # join and order are current documented functions.
        # one-of is something I ran across in old files. Treating it
        # as a sub sequence feature seems appropriate to me.
        # bond is some piece of junk I found in RefSeq files. I have
        # no idea how to interpret it, so I jam it in here
        elif (function.name == "join" or function.name == "order" or
              function.name == "one-of" or function.name == "bond"):
            self._set_ordering_info(function, cur_feature)
        elif (function.name == "gap"):
            assert len(function.args) == 1, \
              "Unexpected number of arguments in gap %s" % function.args
            # make the cur information location a gap object
            position = self._get_position(function.args[0].local_location)
            cur_feature.location = SeqFeature.PositionGap(position)
        else:
            raise ValueError("Unexpected function name: %s" % function.name)

    def _set_ordering_info(self, function, cur_feature):
        """Parse a join or order and all of the information in it.

        This deals with functions that order a bunch of locations,
        specifically 'join' and 'order'. The inner locations are
        added as subfeatures of the top level feature
        """
        # for each inner element, create a sub SeqFeature within the
        # current feature, then get the information for this feature
        cur_feature.location_operator = function.name
        for inner_element in function.args:
            new_sub_feature = SeqFeature.SeqFeature()
            # inherit the type from the parent
            new_sub_feature.type = cur_feature.type 
            # add the join or order info to the location_operator
            new_sub_feature.location_operator = function.name
            # inherit references and strand from the parent feature
            new_sub_feature.ref = cur_feature.ref
            new_sub_feature.ref_db = cur_feature.ref_db
            new_sub_feature.strand = cur_feature.strand

            # set the information for the inner element
            self._set_location_info(inner_element, new_sub_feature)

            # now add the feature to the sub_features
            cur_feature.sub_features.append(new_sub_feature)

        # set the location of the top -- this should be a combination of
        # the start position of the first sub_feature and the end position
        # of the last sub_feature

        # these positions are already converted to python coordinates 
        # (when the sub_features were added) so they don't need to
        # be converted again
        feature_start = cur_feature.sub_features[0].location.start
        feature_end = cur_feature.sub_features[-1].location.end
        cur_feature.location = SeqFeature.FeatureLocation(feature_start,
                                                          feature_end)
        # Historically a join on the reverse strand has been represented
        # in Biopython with both the parent SeqFeature and its children
        # (the exons for a CDS) all given a strand of -1.  Likewise, for
        # a join feature on the forward strand they all have strand +1.
        # However, we must also consider evil mixed strand examples like
        # this, join(complement(69611..69724),139856..140087,140625..140650)
        strands = set(sf.strand for sf in cur_feature.sub_features)
        if len(strands)==1:
            cur_feature.strand = cur_feature.sub_features[0].strand
        else:
            cur_feature.strand = None # i.e. mixed strands

    def _set_location_info(self, parse_info, cur_feature):
        """Set the location information for a feature from the parse info.

        Arguments:
        o parse_info - The classes generated by the LocationParser.
        o cur_feature - The feature to add the information to.
        """
        # base case -- we are out of information
        if parse_info is None:
            return
        # parse a location -- this is another base_case -- we assume
        # we have no information after a single location
        elif isinstance(parse_info, LocationParser.AbsoluteLocation):
            self._set_location(parse_info, cur_feature)
            return
        # parse any of the functions (join, complement, etc)
        elif isinstance(parse_info, LocationParser.Function):
            self._set_function(parse_info, cur_feature)
        # otherwise we are stuck and should raise an error
        else:
            raise ValueError("Could not parse location info: %s"
                             % parse_info)

    def _set_location(self, location, cur_feature):
        """Set the location information for a feature.

        Arguments:
        o location - An AbsoluteLocation object specifying the info
        about the location.
        o cur_feature - The feature to add the information to.
        """
        # check to see if we have a cross reference to another accession
        # ie. U05344.1:514..741
        if location.path is not None:
            cur_feature.ref = location.path.accession
            cur_feature.ref_db = location.path.database
        # now get the actual location information
        cur_feature.location = self._get_location(location.local_location)

    def _get_location(self, range_info):
        """Return a (possibly fuzzy) location from a Range object.

        Arguments:
        o range_info - A location range (ie. something like 67..100). This
        may also be a single position (ie 27).

        This returns a FeatureLocation object.
        If parser.use_fuzziness is set at one, the positions for the
        end points will possibly be fuzzy.
        """
        if isinstance(range_info, LocationParser.Between):
            if not (range_info.low.val+1 == range_info.high.val \
            or range_info.low.val==self._expected_size \
            and range_info.high.val==1):
                raise ValueError(range_info)
            #A between location like "67^68" (one based counting) is a
            #special case (note it has zero length). In python slice
            #notation this is 67:67, a zero length slice.  See Bug 2622
            #Further more, on a circular genome of length N you can have
            #a location N^1 meaning the junction at the origin. See Bug 3098.
            pos = self._get_position(range_info.low)
            return SeqFeature.FeatureLocation(pos, pos)
            #NOTE - We can imagine between locations like "2^4", but this
            #is just "3".  Similarly, "2^5" is just "3..4"
        # check if we just have a single base
        elif not(isinstance(range_info, LocationParser.Range)):
            #A single base like "785" becomes [784:785] in python
            s_pos = self._get_position(range_info)
            # move the single position back one to be consistent with how
            # python indexes numbers (starting at 0)
            s_pos.position = s_pos.position  - 1
            e_pos = self._get_position(range_info)
            return SeqFeature.FeatureLocation(s_pos, e_pos)
        # otherwise we need to get both sides of the range
        else:
            # get *Position objects for the start and end
            start_pos = self._get_position(range_info.low)
            end_pos = self._get_position(range_info.high)

            start_pos.position, end_pos.position = \
              self._convert_to_python_numbers(start_pos.position,
                                              end_pos.position)
            #If the start location is a one-of position, we also need to
            #adjust their positions to use python counting.
            if isinstance(start_pos, SeqFeature.OneOfPosition):
                for p in start_pos.position_choices:
                    p.position -= 1
                
            return SeqFeature.FeatureLocation(start_pos, end_pos)

    def _get_position(self, position):
        """Return a (possibly fuzzy) position for a single coordinate.

        Arguments:
        o position - This is a LocationParser.* object that specifies
        a single coordinate. We will examine the object to determine
        the fuzziness of the position.

        This is used with _get_location to parse out a location of any
        end_point of arbitrary fuzziness.
        """
        # case 1 -- just a normal number
        if (isinstance(position, LocationParser.Integer)):
            final_pos = SeqFeature.ExactPosition(position.val) 
        # case 2 -- we've got a > sign
        elif isinstance(position, LocationParser.LowBound):
            final_pos = SeqFeature.AfterPosition(position.base.val)
        # case 3 -- we've got a < sign
        elif isinstance(position, LocationParser.HighBound):
            final_pos = SeqFeature.BeforePosition(position.base.val)
        # case 4 -- we've got 100^101
        # Is the extension is zero in this example?
        elif isinstance(position, LocationParser.Between):
            #NOTE - We don't *expect* this code to get called!
            #We only except between locations like 3^4 (consecutive)
            #which are handled in _get_location.  We don't expect
            #non consecutive variants like "2^5" as this is just "3..4".
            #Similarly there is no reason to expect composite locations
            #like "(3^4)..6" which should just be "4..6".
            final_pos = SeqFeature.BetweenPosition(position.low.val,
                                 position.high.val-position.low.val)
        # case 5 -- we've got (100.101)
        elif isinstance(position, LocationParser.TwoBound):
            final_pos = SeqFeature.WithinPosition(position.low.val,
                                position.high.val-position.low.val)
        # case 6 -- we've got a one-of(100, 110) location
        elif isinstance(position, LocationParser.Function) and \
                        position.name == "one-of":
            # first convert all of the arguments to positions
            position_choices = []
            for arg in position.args:
                # we only handle AbsoluteLocations with no path
                # right now. Not sure if other cases will pop up
                assert isinstance(arg, LocationParser.AbsoluteLocation), \
                  "Unhandled Location type %r" % arg
                assert arg.path is None, "Unhandled path in location"
                position = self._get_position(arg.local_location)
                position_choices.append(position)
            final_pos = SeqFeature.OneOfPosition(position_choices)
        # if it is none of these cases we've got a problem!
        else:
            raise ValueError("Unexpected LocationParser object %r" %
                             position)

        # if we are using fuzziness return what we've got
        if self._use_fuzziness:
            return final_pos
        # otherwise return an ExactPosition equivalent
        else:
            return SeqFeature.ExactPosition(final_pos.location)

    def _add_qualifier(self):
        """Add a qualifier to the current feature without loss of info.

        If there are multiple qualifier keys with the same name we
        would lose some info in the dictionary, so we append a unique
        number to the end of the name in case of conflicts.
        """
        # if we've got a key from before, add it to the dictionary of
        # qualifiers
        if self._cur_qualifier_key:
            key = self._cur_qualifier_key
            value = "".join(self._cur_qualifier_value)
            if self._feature_cleaner is not None:
                value = self._feature_cleaner.clean_value(key, value)
            # if the qualifier name exists, append the value
            if key in self._cur_feature.qualifiers:
                self._cur_feature.qualifiers[key].append(value)
            # otherwise start a new list of the key with its values
            else:
                self._cur_feature.qualifiers[key] = [value]

    def feature_qualifier_name(self, content_list):
        """When we get a qualifier key, use it as a dictionary key.
        
        We receive a list of keys, since you can have valueless keys such as
        /pseudo which would be passed in with the next key (since no other
        tags separate them in the file)
        """
        for content in content_list:
            # add a qualifier if we've got one
            self._add_qualifier()

            # remove the / and = from the qualifier if they're present
            qual_key = content.replace('/', '')
            qual_key = qual_key.replace('=', '')
            qual_key = qual_key.strip()
            
            self._cur_qualifier_key = qual_key
            self._cur_qualifier_value = []
        
    def feature_qualifier_description(self, content):
        # get rid of the quotes surrounding the qualifier if we've got 'em
        qual_value = content.replace('"', '')
        
        self._cur_qualifier_value.append(qual_value)

    def contig_location(self, content):
        """Deal with CONTIG information."""
        #Historically this was stored as a SeqFeature object, but it was
        #stored under record.annotations["contig"] and not under
        #record.features with the other SeqFeature objects.
        #
        #The CONTIG location line can include additional tokens like
        #Gap(), Gap(100) or Gap(unk100) which are not used in the feature
        #location lines, so storing it using SeqFeature based location
        #objects is difficult.
        #
        #We now store this a string, which means for BioSQL we are now in
        #much better agreement with how BioPerl records the CONTIG line
        #in the database.
        #
        #NOTE - This code assumes the scanner will return all the CONTIG
        #lines already combined into one long string!
        self.data.annotations["contig"] = content

    def origin_name(self, content):
        pass

    def base_count(self, content):
        pass

    def base_number(self, content):
        pass

    def sequence(self, content):
        """Add up sequence information as we get it.

        To try and make things speedier, this puts all of the strings
        into a list of strings, and then uses string.join later to put
        them together. Supposedly, this is a big time savings
        """
        new_seq = content.replace(' ', '')
        new_seq = new_seq.upper()

        self._seq_data.append(new_seq)

    def record_end(self, content):
        """Clean up when we've finished the record.
        """
        from Bio import Alphabet
        from Bio.Alphabet import IUPAC
        from Bio.Seq import Seq, UnknownSeq

        #Try and append the version number to the accession for the full id
        if self.data.id is None:
            assert 'accessions' not in self.data.annotations, \
                   self.data.annotations['accessions']
            self.data.id = self.data.name #Good fall back?
        elif self.data.id.count('.') == 0:
            try:
                self.data.id+='.%i' % self.data.annotations['sequence_version']
            except KeyError:
                pass
        
        # add the last feature in the table which hasn't been added yet
        self._add_feature()

        # add the sequence information
        # first, determine the alphabet
        # we default to an generic alphabet if we don't have a
        # seq type or have strange sequence information.
        seq_alphabet = Alphabet.generic_alphabet

        # now set the sequence
        sequence = "".join(self._seq_data)

        if self._expected_size is not None \
        and len(sequence) != 0 \
        and self._expected_size != len(sequence):
            import warnings
            warnings.warn("Expected sequence length %i, found %i (%s)." \
                          % (self._expected_size, len(sequence), self.data.id))

        if self._seq_type:
            # mRNA is really also DNA, since it is actually cDNA
            if self._seq_type.find('DNA') != -1 or \
               self._seq_type.find('mRNA') != -1:
                seq_alphabet = IUPAC.ambiguous_dna
            # are there ever really RNA sequences in GenBank?
            elif self._seq_type.find('RNA') != -1:
                #Even for data which was from RNA, the sequence string
                #is usually given as DNA (T not U).  Bug 2408
                if "T" in sequence and "U" not in sequence:
                    seq_alphabet = IUPAC.ambiguous_dna
                else:
                    seq_alphabet = IUPAC.ambiguous_rna
            elif self._seq_type.find('PROTEIN') != -1:
                seq_alphabet = IUPAC.protein  # or extended protein?
            # work around ugly GenBank records which have circular or
            # linear but no indication of sequence type
            elif self._seq_type in ["circular", "linear"]:
                pass
            # we have a bug if we get here
            else:
                raise ValueError("Could not determine alphabet for seq_type %s"
                                 % self._seq_type)

        if not sequence and self.__expected_size:
            self.data.seq = UnknownSeq(self._expected_size, seq_alphabet)
        else:
            self.data.seq = Seq(sequence, seq_alphabet)

class _RecordConsumer(_BaseGenBankConsumer):
    """Create a GenBank Record object from scanner generated information.
    """
    def __init__(self):
        _BaseGenBankConsumer.__init__(self)
        import Record
        self.data = Record.Record()

        self._seq_data = []
        self._cur_reference = None
        self._cur_feature = None
        self._cur_qualifier = None
        
    def wgs(self, content):
        self.data.wgs = content.split('-')

    def add_wgs_scafld(self, content):
        self.data.wgs_scafld.append(content.split('-'))

    def locus(self, content):
        self.data.locus = content

    def size(self, content):
        self.data.size = content

    def residue_type(self, content):
        self.data.residue_type = content

    def data_file_division(self, content):
        self.data.data_file_division = content

    def date(self, content):
        self.data.date = content

    def definition(self, content):
        self.data.definition = content

    def accession(self, content):
        for acc in self._split_accessions(content):
            if acc not in self.data.accession:
                self.data.accession.append(acc)

    def nid(self, content):
        self.data.nid = content

    def pid(self, content):
        self.data.pid = content

    def version(self, content):
        self.data.version = content

    def db_source(self, content):
        self.data.db_source = content.rstrip()

    def gi(self, content):
        self.data.gi = content

    def keywords(self, content):
        self.data.keywords = self._split_keywords(content)

    def project(self, content):
        self.data.projects.extend([p for p in content.split() if p])

    def dblink(self, content):
        self.data.dblinks.append(content)

    def segment(self, content):
        self.data.segment = content

    def source(self, content):
        self.data.source = content

    def organism(self, content):
        self.data.organism = content

    def taxonomy(self, content):
        self.data.taxonomy = self._split_taxonomy(content)

    def reference_num(self, content):
        """Grab the reference number and signal the start of a new reference.
        """
        # check if we have a reference to add
        if self._cur_reference is not None:
            self.data.references.append(self._cur_reference)

        import Record
        self._cur_reference = Record.Reference()
        self._cur_reference.number = content

    def reference_bases(self, content):
        self._cur_reference.bases = content

    def authors(self, content):
        self._cur_reference.authors = content

    def consrtm(self, content):
        self._cur_reference.consrtm = content

    def title(self, content):
        if self._cur_reference is None:
            import warnings
            warnings.warn("GenBank TITLE line without REFERENCE line.")
            return
        self._cur_reference.title = content

    def journal(self, content):
        self._cur_reference.journal = content

    def medline_id(self, content):
        self._cur_reference.medline_id = content
        
    def pubmed_id(self, content):
        self._cur_reference.pubmed_id = content

    def remark(self, content):
        self._cur_reference.remark = content
        
    def comment(self, content):
        self.data.comment += "\n".join(content)

    def primary_ref_line(self,content):
        """Data for the PRIMARY line"""
        self.data.primary.append(content)

    def primary(self,content):
        pass
    
    def features_line(self, content):
        """Get ready for the feature table when we reach the FEATURE line.
        """
        self.start_feature_table()

    def start_feature_table(self):
        """Signal the start of the feature table.
        """
        # we need to add on the last reference
        if self._cur_reference is not None:
            self.data.references.append(self._cur_reference)

    def feature_key(self, content):
        """Grab the key of the feature and signal the start of a new feature.
        """
        # first add on feature information if we've got any
        self._add_feature()

        import Record
        self._cur_feature = Record.Feature()
        self._cur_feature.key = content

    def _add_feature(self):
        """Utility function to add a feature to the Record.

        This does all of the appropriate checking to make sure we haven't
        left any info behind, and that we are only adding info if it
        exists.
        """
        if self._cur_feature is not None:
            # if we have a left over qualifier, add it to the qualifiers
            # on the current feature
            if self._cur_qualifier is not None:
                self._cur_feature.qualifiers.append(self._cur_qualifier)

            self._cur_qualifier = None
            self.data.features.append(self._cur_feature)

    def location(self, content):
        self._cur_feature.location = self._clean_location(content)

    def feature_qualifier_name(self, content_list):
        """Deal with qualifier names
        
        We receive a list of keys, since you can have valueless keys such as
        /pseudo which would be passed in with the next key (since no other
        tags separate them in the file)
        """
        import Record
        for content in content_list:
            # the record parser keeps the /s -- add them if we don't have 'em
            if content.find("/") != 0:
                content = "/%s" % content
            # add on a qualifier if we've got one
            if self._cur_qualifier is not None:
                self._cur_feature.qualifiers.append(self._cur_qualifier)

            self._cur_qualifier = Record.Qualifier()
            self._cur_qualifier.key = content

    def feature_qualifier_description(self, content):
        # if we have info then the qualifier key should have a ='s
        if self._cur_qualifier.key.find("=") == -1:
            self._cur_qualifier.key = "%s=" % self._cur_qualifier.key
        cur_content = self._remove_newlines(content)
        # remove all spaces from the value if it is a type where spaces
        # are not important
        for remove_space_key in self.__class__.remove_space_keys:
            if self._cur_qualifier.key.find(remove_space_key) >= 0:
                cur_content = self._remove_spaces(cur_content)
        self._cur_qualifier.value = self._normalize_spaces(cur_content)

    def base_count(self, content):
        self.data.base_counts = content

    def origin_name(self, content):
        self.data.origin = content

    def contig_location(self, content):
        """Signal that we have contig information to add to the record.
        """
        self.data.contig = self._clean_location(content) 

    def sequence(self, content):
        """Add sequence information to a list of sequence strings.

        This removes spaces in the data and uppercases the sequence, and
        then adds it to a list of sequences. Later on we'll join this
        list together to make the final sequence. This is faster than
        adding on the new string every time.
        """
        new_seq = content.replace(' ', '')
        self._seq_data.append(new_seq.upper())

    def record_end(self, content):
        """Signal the end of the record and do any necessary clean-up.
        """
        # add together all of the sequence parts to create the
        # final sequence string
        self.data.sequence = "".join(self._seq_data)
        # add on the last feature
        self._add_feature()


