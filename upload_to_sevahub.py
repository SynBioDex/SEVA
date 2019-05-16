from sbol import *

setHomespace('http://seva.cnb.csic.es/')

# Submit vectors one at a time
sevahub = PartShop('http://sevahub.es')
sevahub.login('bartleyba@sbolstandard.org')
for origin_id in range(1, 10):
    for marker_id in range(1, 7):
        vector_name = 'pSEVA%d%d1' %(marker_id, origin_id)
        doc = Document(vector_name + '.xml')
        doc.displayId = vector_name
        doc.name = vector_name
        doc.description = vector_name
        doc.version = '1'
        response = sevahub.submit(doc, 'http://sevahub.es/public/%s/%s_collection/1' %(vector_name, vector_name), 1);
        print(response)


# Integrate all vectors
# doc = Document()
# for origin_id in range(1, 10):
#     for marker_id in range(1, 7):
#     	print(vector_name)
#         vector_name = 'pSEVA%d%d1' %(marker_id, origin_id)
#         doc.append(vector_name + '.xml')

# doc.displayId = "pSEVA"
# doc.name = "pSEVA"
# doc.description = "SEVA vectors with default carge"
# doc.version = '1'
# doc.write('pSEVA.xml')
#response = sevahub.submit(doc, 'http://sevahub.es/user/bartleyba/%s/%s_collection/1' %(vector_name, vector_name), 1);
# print(response)
