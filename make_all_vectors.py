import sys
sys.path
sys.path.append('/Users/bbartley/Dev/git/libSBOL/release/wrapper/Mac_64_2/sbol')

from seva import *

### Create default SEVA vectors
Config.setOption('validate', True)
# Config.setOption('sbol_typed_uris', False)
setHomespace('http://seva.cnb.csic.es')
for origin_id in range(1, 10):
    for marker_id in range(1, 7):
        vector_name = 'pSEVA%d%d1' %(marker_id, origin_id)
    	print vector_name
        c = Cargo('cargo')
        v = SEVAVector(vector_name, marker_id, origin_id, c)
        v.name = vector_name

