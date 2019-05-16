from sbol import *

SO_ORIGIN = SO + "0000296"

doc = Document('collections.xml')
selection_markers = []
origins = []
for cd in doc.componentDefinitions:
	if SO_GENE in cd.roles:
		selection_markers.append(cd)
	if SO_ORIGIN in cd.roles:
		origins.append(cd)
selection_markers_doc = Document()
for cd in selection_markers:
	new_cd = cd.copy(selection_markers_doc)
	seq = doc.getSequence(cd.sequence)
	new_seq = seq.copy(selection_markers_doc)
	new_cd.wasDerivedFrom = None
	new_seq.wasDerivedFrom = None
selection_markers_doc.displayId = "selection_markers"
selection_markers_doc.name = "selection markers"
selection_markers_doc.description = "Selection markers for SEVA vectors"
selection_markers_doc.version = '1'
print selection_markers_doc.write('SEVA_selection_markers.xml')

origins_doc = Document()
for cd in origins:
	new_cd = cd.copy(origins_doc)
	seq = doc.getSequence(cd.sequence)
	new_seq = seq.copy(origins_doc)
	new_cd.wasDerivedFrom = None
	new_seq.wasDerivedFrom = None

origins_doc.displayId = "origins_of_replication"
origins_doc.name = "origins of replication"
origins_doc.description = "Origins of replication for SEVA vectors"
origins_doc.version = '1'
print origins_doc.write('SEVA_origins.xml')

seva_subparts_doc = Document('template.xml')
seva_subparts_doc.displayId = "seva_template"
seva_subparts_doc.name = "SEVA template"
seva_subparts_doc.description = "Subparts for SEVA vector template"
seva_subparts_doc.version = '1'

sevahub = PartShop('http://sevahub.es')
sevahub.login('bartleyba@sbolstandard.org')
response = sevahub.submit(origins_doc, 'http://sevahub.es/user/bartleyba/%s/%s_collection/1' %(origins_doc.displayId, origins_doc.displayId), 1)
print response
response = sevahub.submit(selection_markers_doc, 'http://sevahub.es/user/bartleyba/%s/%s_collection/1' %(selection_markers_doc.displayId, selection_markers_doc.displayId), 1)
print response
response = sevahub.submit(seva_subparts_doc, 'http://sevahub.es/user/bartleyba/%s/%s_collection/1' %(seva_subparts_doc.displayId, seva_subparts_doc.displayId), 1)
print response
