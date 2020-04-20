# PyDesc cookbook


## Configuration


## Load structures

To quickly load structure and be able to read and edit it, simply run:

```python
from pydesc.structure import StructureLoader

structure_loader = StructureLoader()
models_list = structure_loader.load_structures('1no5')
structure = models_list[0]  # if there is more than one model in pdb file
                            # they are stored in resulting list

for mer in structure:
    print(mer)          # residues, nucleotides or others (ions or ligands)
    for atom in mer:    # or mer.atoms
        print(atom)
    for pseudoatom in mer.pseudoatoms:
        print(pseudoatom)

```


### There is more

To gain more control over process of loading structure, note that there are details
that could be changed. Lets take a look on the example (explained below):

```python
from pydesc import dbhandler
from pydesc.mers import factories as mer_factories
from pydesc.structure import StructureLoader


db_handler = dbhandler.MetaHandler(mode=(3, 2, 1))
file_parser = dbhandler.MetaParser()
mer_factory = mer_factories.BioPythonMerFactory()

structure_loader = StructureLoader(
    handler=db_handler,
    parser=file_parser,
    mer_factory=mer_factory,
    )

models_db = structure_loader.load_structures(code="1no5")  # download from db
models_hd = structure_loader.load_structures(path='tests/data/dbtest/PDB/1nkl')
                                                           # download from hard drive

```

`db_handler = dbhandler.MetaHandler(mode=(3, 2, 1))`
Handlers are objects responsible for handling structure files.
There is at least one for each data base.
Meta handler contains them all.
Learn more [here](#structure-file-handlers).

`file_parser = dbhandler.MetaParser()`
Pydesc uses external file parsers, e.g. from BioPython and MDTraj.
MetaParser consists of PDBParser and mmCIFParser and executes both on given input.
Learn more [here](#file-parsers)

`mer_factory = mer_factories.BioPythonMerFactory()`
Mer factory is object responsible for creation of structures building blocks.
Its settings are important whenever one needs to change way of how and which mers are
 created.
Learn more [here](#mers).
 
`models_db = structure_loader.load_structures(code="1no5")`
That line loads structure using file handler.
Argument `code` is passed to handler, so it can be any string comprehensive for such 
 object, e.g. "cath://1no5" for MetaHandler or "bio://1no5/1" for BioUnitHandler.
Learn more [here](#structure-file-handlers).

`models_hd = structure_loader.load_structures(path='tests/data/dbtest/PDB/1nkl')`
This command loads file form hard drive. Handler is surpassed, file content is given
 directly to parser.

#### Structure file handlers

Handlers are objects responsible for handling structure files.
There is at least one for each data base.
Meta handler contains them all.
Handlers store downloaded files in local cache.
Therefor users can set their behaviour by setting modes and their priorities:
* in mode 3 handler reads local cache
* in mode 2 handler copies file from local db to local cache (overwrite)
* in mode 1 handler downloads file from remote db to local cache (overwrite)
Handler tries to get files in different modes 

#### File parsers

PyDesc does not provide any structure file parser itself.
Instead it uses BioPython parsers to read `.pdb` and `.cif` files and MDTraj parsers
to read different trajectories.
`StructureLoader` can accept:
* BioPython's PDBParser
* BioPython's MMCIFParser
* PyDesc's MetaParser
* any object with method and signature such as: `get_structure(code: str, file_path: 
str)`

### Mers

### Substructures

### Descriptors

### Number converters

## Trajectories


## Selections

## Contact maps

### Contact Criteria

## Structure comparison

### Overfit

### Compdesc

### FitDesc

## Alignments

## DescMOL
