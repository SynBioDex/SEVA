import sys
sys.path
sys.path.append('/Users/bbartley/Dev/git/libSBOL/release/wrapper/Mac_64_2')

from seva import *

### Create default SEVA vectors
Config.setOption('validate', True)
Config.setOption('sbol_typed_uris', False)
setHomespace('http://seva.cnb.csic.es')
for origin_id in SEVA_ORIGINS.keys():
    for marker_id in SEVA_SELECTION_MARKERS.keys():
        for cargo_id in SEVA_CARGOES.keys():
            vector_name = 'pSEVA_%s_%s_%s' %(marker_id, origin_id, cargo_id)
            print vector_name
            v = SEVAVector(vector_name, marker_id, origin_id, cargo_id)
            v.name = vector_name
