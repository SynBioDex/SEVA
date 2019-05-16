import sys
sys.path
sys.path.append('/Users/bbartley/Dev/git/libSBOL/release/wrapper/Mac_64_2')

from sbol import *

SO_FIVE_PRIME_STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE = SO + "0001975"
SO_THREE_PRIME_STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE = SO + "0001976"
SO_FIVE_PRIME_STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE = SO + "0001975"
SO_THREE_PRIME_STICKY_END_RESTRICTION_ENZYME_CLEAVAGE_SITE = SO + "0001976"
SO_ENGINEERED_GENE = SO + "0000280"
SO_ORIGIN = SO + "0000296"
SO_RESTRICTION_SITE = SO + "0001687"
SO_SCAR = SO + "0001953"

SEVA_VECTOR_DOC = Document()  # Used globally by Cargo and SEVAVector classes

SEVA_SELECTION_MARKERS = {'1': 'Ap', 
                          '2': 'Km',
                          '3': 'Cm', 
                          '4': 'Sm', 
                          '5': 'Tc', 
                          '6': 'Gm', 
                          '7': 'Tmp', 
                          '8': 'Apr'
                         }
SEVA_ORIGINS = {'1': 'R6K',
                '2': 'RK2',
                '2a': 'RK2_int_phiC31',
                '2b': 'RK2_SCP2',
                '2c': 'pIJ101_RK2',
                '3': 'pBBR1',
                '3a': 'pBBR1_int_phiC31',
                '3b': 'pBBR1_SCP2',
                '4': 'pRO1600_ColE1',
                '5': 'RSF1010',
                '6': 'p15A',
                '7': 'pSC101',
                '8': 'pUC',
                '8a': 'pUC_int_phiC31',
                '8b': 'pUC_SCP2',
                '8c': 'pIJ101_pUC',
                '9': 'pBBR322_ROP'
                }
SEVA_CARGOES = {'1': 'Mcs_Default',
                '2': 'lacZ_pUC19',
                '2S': 'lacZ_pUC19_I_SceI',
                '3': 'lacZ_pUC18',
                '3X': 'lacZalpha_pUC18_IIs_sites_up',
                '3Y': 'lacZalpha_pUC18_IIs_sites_down',
                '4': 'lacIq_Ptrc',
                '4R': 'lacIq_Ptac',
                '4F': 'T5_lacO',
                '4E': 'T7_lacO',
                '5': 'lacZ',
                '5T': 'lacZ_translational',
                '6': 'luxCDABE',
                '7': 'GFP',
                '7Y': 'YFP',
                '7C': 'CFP',
                '7D': 'DsRed2',
                '7R': 'mCherry',
                '7F': 'EcFbFP',
                '7M': 'msfGFP',
                '7V': 'GFP_LVA',
                '8': 'xylS_Pm',
                '8S': 'xylS_Pm_RBS',
                '9': 'alkS_PalkB',
                '10': 'araC_pBAD',
                '11': 'ChnR_PchnB',
                '12': 'CprK1_PDB3',
                '13': 'PEM7',
                '14': 'cI857_PL',
                '15': 'cos_T7_ccdB_T3',
                'alpha': 'Hok_sok'}

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

        # insert = insert_component_definition.copy(SEVA_VECTOR_DOC)
        # sequence = insert_component_definition.sequence.copy(SEVA_VECTOR_DOC) 
        # insert.sequence = sequence

        # if not insert.sequence:
        #    raise Exception("Insertion failed. ComponentDefinition is not associated with a Sequence")

        if not insert_component_definition.sequence:
           raise Exception("Insertion failed. Insert ComponentDefinition is not associated with a Sequence")
        sequence = insert_component_definition.sequence

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

    def __init__(self, id = 'example', marker_id = '1', origin_id = '1', cargo_id = '1'):
        ComponentDefinition.__init__(self, id)

        if not marker_id in SEVA_SELECTION_MARKERS.keys() : raise Exception()
        if not origin_id in SEVA_ORIGINS.keys() : raise Exception()
        if not type(cargo_id) == Cargo and not cargo_id in SEVA_CARGOES.keys(): raise Exception()
        self.name = id

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
        cargoes = {}
        for cd in SEVA_PARTS.componentDefinitions:
            if len(cd.roles):
                if cd.roles[0] == SO_GENE:
                    selection_markers[cd.displayId] = cd
                if cd.roles[0] == SO_ORIGIN:
                    origins[cd.displayId] = cd
                if cd.roles[0] == SO_RESTRICTION_SITE:
                    restriction_sites[cd.displayId] = cd
                if cd.roles[0] == SO_ENGINEERED_GENE:
                    cargoes[cd.displayId] = cd

        if not type(cargo_id) == Cargo:
            c = Cargo(SEVA_CARGOES[cargo_id])
            c.insert(cargoes[SEVA_CARGOES[cargo_id]], 'PstI', 'SpeI')
        else:
            c = cargo_id

        selection_marker = SEVA_PARTS.componentDefinitions[ 'cd/' + SEVA_SELECTION_MARKERS[marker_id] ]
        origin = SEVA_PARTS.componentDefinitions[ 'cd/' + SEVA_ORIGINS[ origin_id ]]
        t1 = SEVA_PARTS.componentDefinitions['cd/T1']
        # pacI = SEVA_PARTS.componentDefinitions['cd/PacI']
        # speI = SEVA_PARTS.componentDefinitions['cd/SpeI']
        t0 = SEVA_PARTS.componentDefinitions['cd/T0']
        sanDI = SEVA_PARTS.componentDefinitions['cd/SanDI']
        swaI = SEVA_PARTS.componentDefinitions['cd/SwaI']
        # pshAI = SEVA_PARTS.componentDefinitions['cd/PshAI']

        pshAI = None
        if marker_id == '1':
            pshAI = SEVA_PARTS.componentDefinitions['cd/PshAI_v1']
        elif marker_id == '2':
            pshAI = SEVA_PARTS.componentDefinitions['cd/PshAI_v2']
        elif marker_id == '3':
            pshAI = SEVA_PARTS.componentDefinitions['cd/PshAI_v3']
        elif marker_id == '4':
            pshAI = SEVA_PARTS.componentDefinitions['cd/PshAI_v4']
        elif marker_id == '5':
            pshAI = SEVA_PARTS.componentDefinitions['cd/PshAI_v5']
        elif marker_id == '6':
            pshAI = SEVA_PARTS.componentDefinitions['cd/PshAI_v6']
        elif marker_id == '7':
            pshAI = SEVA_PARTS.componentDefinitions['cd/PshAI_v7']
        elif marker_id == '8':
            pshAI = SEVA_PARTS.componentDefinitions['cd/PshAI_v8']

        oriT =SEVA_PARTS.componentDefinitions['cd/OriT']
        fseI = SEVA_PARTS.componentDefinitions['cd/FseI']
        ascI = SEVA_PARTS.componentDefinitions['cd/AscI']

        selection_marker_seq = SEVA_PARTS.sequences['seq/seq' + SEVA_SELECTION_MARKERS[marker_id]]
        origin_seq = SEVA_PARTS.sequences['seq/seq' + origin.displayId]
        t1_seq = SEVA_PARTS.sequences['seq/seqT1']
        # pacI_seq = SEVA_PARTS.sequences['seq/seqPacI']
        speI_seq = SEVA_PARTS.sequences['seq/seqSpeI']
        t0_seq = SEVA_PARTS.sequences['seq/seqT0']
        sanDI_seq = SEVA_PARTS.sequences['seq/seqSanDI']
        swaI_seq = SEVA_PARTS.sequences['seq/seqSwaI']
        # pshAI_seq = SEVA_PARTS.sequences['seq/seqPshAI']
        oriT_seq =SEVA_PARTS.sequences['seq/seqOriT']
        fseI_seq = SEVA_PARTS.sequences['seq/seqFseI']
        ascI_seq = SEVA_PARTS.sequences['seq/seqAscI']

        pshAI_seq = None
        if marker_id == '1':
            pshAI_seq = SEVA_PARTS.sequences['seq/seqPshAI_v1']
        elif marker_id == '2':
            pshAI_seq = SEVA_PARTS.sequences['seq/seqPshAI_v2']
        elif marker_id == '3':
            pshAI_seq = SEVA_PARTS.sequences['seq/seqPshAI_v3']
        elif marker_id == '4':
            pshAI_seq = SEVA_PARTS.sequences['seq/seqPshAI_v4']
        elif marker_id == '5':
            pshAI_seq = SEVA_PARTS.sequences['seq/seqPshAI_v5']
        elif marker_id == '6':
            pshAI_seq = SEVA_PARTS.sequences['seq/seqPshAI_v6']
        elif marker_id == '7':
            pshAI_seq = SEVA_PARTS.sequences['seq/seqPshAI_v7']
        elif marker_id == '8':
            pshAI_seq = SEVA_PARTS.sequences['seq/seqPshAI_v8']

        Config.setOption('sbol_compliant_uris', False)
        scar1 = SEVA_VECTOR_DOC.componentDefinitions.create('http://seva.cnb.csic.es/cd/scar1') 
        scar1.displayId = 'scar1'       
        scar1.name = 'scar1'
        scar1.roles = [SO_SCAR]
        if marker_id in ['4', '5', '7']: 
            scar1_seq = Sequence('http://seva.cnb.csic.es/seq/scar1', 'CAATAATTACGATTTACGT')
        else:
            scar1_seq = Sequence('http://seva.cnb.csic.es/seq/scar1', 'CAATAATTACG')
        scar1_seq.displayId = 'scar1_seq'
        scar1.sequence = scar1_seq

        # Conditional scar
        if marker_id in ['2', '6', '8']: 
            scar2 = SEVA_VECTOR_DOC.componentDefinitions.create('http://seva.cnb.csic.es/cd/scar2') 
            scar2.displayId = 'scar2'       
            scar2.name = 'scar2'
            scar2.roles = [SO_SCAR]
            scar2_seq = Sequence('http://seva.cnb.csic.es/seq/scar2', 'CGCGCGTTGTC')
            scar2_seq.displayId = 'scar2_seq'
            scar2.sequence = scar2_seq

        # Restore configuration options
        setHomespace(homespace)
        Config.setOption('sbol_typed_uris', sbol_typed_uris)
        Config.setOption('sbol_compliant_uris', sbol_compliant_uris)

        # Copy parts into user's namespace
        selection_marker = selection_marker.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        origin = origin.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        t1 = t1.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        # pacI = pacI.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        # speI = speI.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
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
        # pacI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        # speI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        t0_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        sanDI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        swaI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        pshAI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        oriT_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        fseI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')
        ascI_seq.copy(SEVA_VECTOR_DOC, 'http://seva.cnb.csic.es')

        # Scar after pshAI is conditional on antibiotic marker
        primary_structure = [c, t0, sanDI, scar1, swaI, selection_marker, pshAI, oriT, fseI, origin, ascI, t1]
        if marker_id in ['2', '6', '8']: 
            primary_structure = [c, t0, sanDI, scar1, swaI, selection_marker, pshAI, scar2, oriT, fseI, origin, ascI, t1]

        self.assemblePrimaryStructure(primary_structure)
        compiled_sequence = self.compile() 

        response = SEVA_VECTOR_DOC.write(id + '.xml')
        print(response)
        response = SEVA_VECTOR_DOC.convert('GenBank', id + '.gb')
        print(response)
        # response = SEVA_VECTOR_DOC.convert('FASTA', id + '.fasta')
        # print(response)

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





