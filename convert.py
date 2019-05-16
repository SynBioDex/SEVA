import sys
sys.path
sys.path.append('/Users/bbartley/Dev/git/libSBOL/release/wrapper/Mac_64_2')

from sbol import *

#doc = Document('pSEVA111.xml')
doc = Document()
cd = doc.componentDefinitions.create('cd')
cd.sequence = Sequence('s')
cd.sequence.elements = 'ACTG'

print(doc.convert('GenBank', 'test.gb'))
print(doc.convert('FASTA', 'test.fasta'))
'''
request = { 'options': {'language' : 'GenBank',
                        'test_equality': False,
                        'check_uri_compliance': False,
                        'check_completeness': False,
                        'check_best_practices': False,
                        'fail_on_first_error': False,
                        'provide_detailed_stack_trace': False,
                        'subset_uri': '',
                        'uri_prefix': '',
                        'version': '',
                        'insert_type': False,
                        'main_file_name': 'main file',
                        'diff_file_name': 'comparison file',
                                },
            'return_file': True,
            'main_file': file
          }
'''
