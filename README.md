## SEVA in SBOL

The script seva.py constructs an SBOL representation of plasmids that conform to the Standard European Vector Architecture (SEVA).  Custom cargos can be ``cloned'' into the SEVA multiple cloning site.

## DEPENDENCIES

The seva.py script depends on **pySBOL**. For installation options, go to [Installation Page](https://pysbol2.readthedocs.io/en/latest/installation.html).

## EXAMPLE

```
# Create a custom insert in the Cargo region
doc = Document()
cd = doc.componentDefinitions.create('insert')
seq = doc.sequences.create('insert')
seq.elements = 'actg'
cd.sequence = seq
c = Cargo('my_cargo')
c.insert(cd, 'PstI', 'SpeI')  # insert the cargo between the PstI and SpeI restriction sites

# Write the SEVA vector as SBOL
marker_id = 1  # According to SEVA standard, the selection marker can be 1 - 6
origin_id = 1  # The origin marker can be 1 - 9
v = SEVAVector('pSEVA11my_cargo', marker_id, origin_id, c)  # Writes the file pSEVAmy_cargo.xml
```