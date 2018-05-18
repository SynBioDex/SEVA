from seva import *

### Create default SEVA vectors
Config.setOption('validate', True)
setHomespace('http://seva.cnb.csic.es')
for origin_id in range(1, 10):
    for marker_id in range(1, 7):
        vector_name = 'pSEVA%d%d1' %(marker_id, origin_id)
    	print vector_name
        c = Cargo('cargo')
        v = SEVAVector(vector_name, marker_id, origin_id, c)
        v.name = vector_name

