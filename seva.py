import sys
sys.path
sys.path.append('/Users/bbartley/Dev/git/libSBOL/release/wrapper/Mac_64_2/sbol')

from sbol import *

SO_FIVE_PRIME_STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE = SO + "0001975"
SO_THREE_PRIME_STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE = SO + "0001976"
SO_FIVE_PRIME_STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE = SO + "0001975"
SO_THREE_PRIME_STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE = SO + "0001976"
SO_ENGINEERED_GENE = SO + "0000280"
SO_ORIGIN = SO + "0000296"
SO_RESTRICTION_SITE = SO + "0001687"

SEVA_VECTOR_DOC = Document()  # Used globally by Cargo and SEVAVector classes

class Cargo(ComponentDefinition, PythonicInterface):               # Inheriting from PythonicInterface converts the C++ style API to a more Pythonic API
    def __init__(self, id = 'example'):
        ComponentDefinition.__init__(self, id)
        SEVA_VECTOR_DOC.addComponentDefinition(self)
        self.default_sequence = SEVA_VECTOR_DOC.sequences.create(id + '_seq')
        self.default_sequence.elements = 'TTAATTAAAGCGGATAACAATTTCACACAGGAGGCCGCCTAGGCCGCGGCCGCGCGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGCGGCCGCGTCGTGACTGGGAAAACCCTGGCGACTAGT'
        self.sequence = self.default_sequence
        self.name = id
        sa = self.sequenceAnnotations.create('PacI')
        r = sa.locations.createRange('r')
        r.start = 1
        r.end = 8
        sa = self.sequenceAnnotations.create('R24')
        r = sa.locations.createRange('r')
        r.start = 9
        r.end = 32
        sa = self.sequenceAnnotations.create('SfiI')
        r = sa.locations.createRange('r')
        r.start = 33
        r.end = 45
        sa = self.sequenceAnnotations.create('AvrII')
        r = sa.locations.createRange('r')
        r.start = 38
        r.end = 43
        sa = self.sequenceAnnotations.create('NotI')
        r = sa.locations.createRange('r')
        r.start = 46
        r.end = 53
        sa = self.sequenceAnnotations.create('EcoRI')
        r = sa.locations.createRange('r')
        r.start = 56
        r.end = 61
        sa = self.sequenceAnnotations.create('SacI')
        r = sa.locations.createRange('r')
        r.start = 62
        r.end = 67
        sa = self.sequenceAnnotations.create('KpnI')
        r = sa.locations.createRange('r')
        r.start = 68
        r.end = 73
        sa = self.sequenceAnnotations.create('SmaI')
        r = sa.locations.createRange('r')
        r.start = 72
        r.end = 77
        sa = self.sequenceAnnotations.create('BamHI')
        r = sa.locations.createRange('r')
        r.start = 77
        r.end = 82
        sa = self.sequenceAnnotations.create('XbaI')
        r = sa.locations.createRange('r')
        r.start = 83
        r.end = 88
        sa = self.sequenceAnnotations.create('SalI')
        r = sa.locations.createRange('r')
        r.start = 89
        r.end = 94
        sa = self.sequenceAnnotations.create('PstI')
        r = sa.locations.createRange('r')
        r.start = 95
        r.end = 100
        sa = self.sequenceAnnotations.create('SphI')
        r = sa.locations.createRange('r')
        r.start = 101
        r.end = 106
        sa = self.sequenceAnnotations.create('HindIII')
        r = sa.locations.createRange('r')
        r.start = 107
        r.end = 112
        sa = self.sequenceAnnotations.create('NotI_2')
        r = sa.locations.createRange('r')
        r.start = 113
        r.end = 120
        sa = self.sequenceAnnotations.create('F24')
        r = sa.locations.createRange('r')
        r.start = 121
        r.end = 144
        sa = self.sequenceAnnotations.create('SpeI')
        r = sa.locations.createRange('r')
        r.start = 145
        r.end = 150

    def insert(self, insert_component_definition, upstream_restriction_site, downstream_restriction_site):
        if not insert_component_definition.doc:
           raise Exception("Insertion failed. The specified ComponentDefinition is not associated with a Document")
        insert = insert_component_definition.copy(SEVA_VECTOR_DOC)
        if not insert.sequence:
           raise Exception("Insertion failed. The specified ComponentDefinition is not associated with a Sequence")
        sequence_id = insert.sequence
        sequence = insert_component_definition.doc.getSequence(sequence_id).copy(SEVA_VECTOR_DOC)

        sa_upstream = None
        sa_downstream = None
        for sa in self.sequenceAnnotations:
            if sa.displayId == upstream_restriction_site:
                sa_upstream = sa
            if sa.displayId == downstream_restriction_site:
                sa_downstream = sa

        # Validate arguments
        if not (sa_upstream and sa_downstream): raise Exception('Invalid restriction sites specified for insertion')
        if sa_upstream.locations['r'].start >= sa_downstream.locations['r'].start:
            raise Exception('Insertion failed. %s is not downstream of %s' %(sa_downstream.displayId, sa_upstream.displayId))

        # Partition cargo region into upstream, deleted, and downstream regions relative to insertion site
        downstream = []
        upstream = []
        deleted = []
        for sa in self.sequenceAnnotations:
            # if sa.locations['r'].start <= sa_upstream.locations['r'].start:
            #     upstream.append(sa)
            # elif sa.locations['r'].start >= sa_downstream.locations['r'].start:
            #     downstream.append(sa)
            if sa.precedes(sa_upstream):
                upstream.append(sa)
            elif sa.follows(sa_downstream):
                downstream.append(sa)
            else:
                deleted.append(sa)

        # Recalculate nucleotide sequence
        upstream_nucleotides = self.default_sequence.elements[:sa_upstream.locations['r'].end]
        downstream_nucleotides = self.default_sequence.elements[(sa_downstream.locations['r'].start-1):]
        self.default_sequence.elements = upstream_nucleotides + sequence.elements + downstream_nucleotides

        # Delete region between restriction sites
        for sa in deleted:
            self.sequenceAnnotations.remove(sa.identity)
        deletion_size = sa_downstream.locations['r'].start - sa_upstream.locations['r'].end

        insert = self.sequenceAnnotations.create(insert_component_definition.displayId)
        r = insert.locations.createRange('r')
        r.start = sa_upstream.locations['r'].end + 1
        r.end = sa_upstream.locations['r'].end + len(sequence.elements)
        for sa in downstream:
            sa.locations['r'].start = sa.locations['r'].start + len(sequence.elements) - deletion_size + 1
            sa.locations['r'].end = sa.locations['r'].end + len(sequence.elements) - deletion_size + 1


class SEVAVector(ComponentDefinition, PythonicInterface):               # Inheriting from PythonicInterface converts the C++ style API to a more Pythonic API

    def __init__(self, id = 'example', marker_id = 1, origin_id = 1, cargo = None):
        ComponentDefinition.__init__(self, id)

        if type(marker_id) != int : raise Exception()
        if type(origin_id) != int : raise Exception()
        if type(cargo) != Cargo : raise Exception()
        if not cargo : raise Exception()
        self.name = id

        SEVA_selection_markers = {
                            1 : 'Ap',
                            2 : 'Km',
                            3 : 'Cm',
                            4 : 'Sm',
                            5 : 'Tc',
                            6 : 'Gm'
                            }
        SEVA_origins = {    
                            1 : 'R6K',
                            2 : 'RK2',
                            3 : 'pBBR1',
                            4 : 'pRO1600_ColE1',
                            5 : 'RSF1010',
                            6 : 'p15A',
                            7 : 'pSC101',
                            8 : 'pUC',
                            9 : 'pBBR322_ROP'
                            }

        SEVA_VECTOR_DOC.addComponentDefinition(self)

        # Get current configuration options, so we can restore them later
        homespace = getHomespace()
        sbol_typed_uris = Config.getOption('sbol_typed_uris')
        sbol_compliant_uris = Config.getOption('sbol_compliant_uris')

        # Set configuration options for SEVA parts collections
        setHomespace('http://seva.cnb.csic.es')
        Config.setOption('sbol_typed_uris', False)
        Config.setOption('sbol_compliant_uris', True)

        SEVA_PARTS = Document('template.xml')
        SEVA_PARTS.append('collections.xml')
        selection_markers = {}
        origins = {}
        restriction_sites = {}
        for cd in SEVA_PARTS.componentDefinitions:
            if len(cd.roles):
                if cd.roles[0] == SO_GENE:
                    selection_markers[cd.name] = cd
                if cd.roles[0] == SO_ORIGIN:
                    origins[cd.name] = cd
                if cd.roles[0] == SO_RESTRICTION_SITE:
                    restriction_sites[cd.name] = cd


        selection_marker = SEVA_PARTS.componentDefinitions[ 'cd/' + SEVA_selection_markers [marker_id] ]
        origin = SEVA_PARTS.componentDefinitions[ 'cd/' + SEVA_origins [ origin_id ]]
        t1 = SEVA_PARTS.componentDefinitions['cd/T1']
        pacI = SEVA_PARTS.componentDefinitions['cd/PacI']
        speI = SEVA_PARTS.componentDefinitions['cd/SpeI']
        t0 = SEVA_PARTS.componentDefinitions['cd/T0']
        sanDI = SEVA_PARTS.componentDefinitions['cd/SanDI']
        swaI = SEVA_PARTS.componentDefinitions['cd/SwaI']
        pshAI = SEVA_PARTS.componentDefinitions['cd/PshAI']
        oriT =SEVA_PARTS.componentDefinitions['cd/OriT']
        fseI = SEVA_PARTS.componentDefinitions['cd/FseI']
        ascI = SEVA_PARTS.componentDefinitions['cd/AscI']

        selection_marker_seq = SEVA_PARTS.sequences['seq/seq' + SEVA_selection_markers [marker_id]]
        origin_seq = SEVA_PARTS.sequences['seq/seq' + origin.name]
        t1_seq = SEVA_PARTS.sequences['seq/seqT1']
        pacI_seq = SEVA_PARTS.sequences['seq/seqPacI']
        speI_seq = SEVA_PARTS.sequences['seq/seqSpeI']
        t0_seq = SEVA_PARTS.sequences['seq/seqT0']
        sanDI_seq = SEVA_PARTS.sequences['seq/seqSanDI']
        swaI_seq = SEVA_PARTS.sequences['seq/seqSwaI']
        pshAI_seq = SEVA_PARTS.sequences['seq/seqPshAI']
        oriT_seq =SEVA_PARTS.sequences['seq/seqOriT']
        fseI_seq = SEVA_PARTS.sequences['seq/seqFseI']
        ascI_seq = SEVA_PARTS.sequences['seq/seqAscI']

        Config.setOption('sbol_compliant_uris', False)
        spacer = SEVA_VECTOR_DOC.componentDefinitions.create('http://seva.cnb.csic.es/cd/spacer') 
        spacer.displayId = 'spacer'       
        spacer.name = 'spacer'
        spacer_seq = Sequence('http://seva.cnb.csic.es/seq/spacer', 'CAATAATTACG')
        spacer_seq.displayId = 'spacer'
        spacer.sequence = spacer_seq

        # Restore configuration options
        setHomespace(homespace)
        Config.setOption('sbol_typed_uris', sbol_typed_uris)
        Config.setOption('sbol_compliant_uris', sbol_compliant_uris)

        # Copy parts into user's namespace
        selection_marker = selection_marker.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        origin = origin.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        t1 = t1.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        pacI = pacI.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        speI = speI.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        t0 = t0.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        sanDI = sanDI.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        swaI = swaI.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        pshAI = pshAI.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        oriT = oriT.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        fseI = fseI.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        ascI = ascI.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')

        selection_marker_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        origin_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        t1_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        pacI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        speI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        t0_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        sanDI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        swaI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        pshAI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        oriT_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        fseI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        ascI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')

        self.assemblePrimaryStructure([cargo, t0, sanDI, spacer, swaI, selection_marker, pshAI, oriT, fseI, origin, ascI, t1])
        compiled_sequence = self.compile() 
        response = SEVA_VECTOR_DOC.write(id + '.xml')
        print(response)

        # Clear global cache
        cd_ids = []
        for cd in SEVA_VECTOR_DOC.componentDefinitions:
            cd_ids.append(cd.identity)
        for cd_id in cd_ids:
            SEVA_VECTOR_DOC.componentDefinitions.remove(cd_id)
        seq_ids = []
        for seq in SEVA_VECTOR_DOC.sequences:
            seq_ids.append(seq.identity)
        for seq_id in seq_ids:
            SEVA_VECTOR_DOC.sequences.remove(seq_id)




