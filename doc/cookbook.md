# PyDesc cookbook

<!--

Add sections using template below.
After adding section -- run md_toc and paste result as TOC.

examples that are to be tested place in "```python <code> "
examples that should not be tested by pytest -- in " ```python <code>" (note space at
 start).

==== TEMPLATE

    ## Template

    Description

    ### Configuration

    There is no related configuration.

    ### API

    TBD

    ### Simple usage

    TBD

    ### There is more 

    TBD

==== TEMPLATE END
-->

  - [Configuration](#configuration)
  - [Loading structures](#loading-structures)
    - [Configuration](#configuration-1)
    - [API](#api)
    - [Simple usage](#simple-usage)
    - [Database handlers](#database-handlers)
    - [File parsers](#file-parsers)
    - [Number converters](#number-converters)
  - [AtomSet - structure building block](#atomset---structure-building-block)
    - [Configuration](#configuration-2)
    - [API](#api-1)
    - [Simple usage (AtomSet trivia)](#simple-usage-atomset-trivia)
    - [AtomSet factories](#atomset-factories)
    - [Full-atom representation](#full-atom-representation)
    - [Martini](#martini)
    - [Backbone trace (P-trace, CA-trace)](#backbone-trace-p-trace-ca-trace)
    - [User defined representation](#user-defined-representation)
  - [Substructures](#substructures)
    - [Configuration](#configuration-3)
    - [API](#api-2)
    - [Simple usage](#simple-usage-1)
    - [Descriptors](#descriptors)
  - [Trajectories](#trajectories)
    - [Configuration](#configuration-4)
    - [API](#api-3)
    - [Simple usage](#simple-usage-2)
  - [Selections](#selections)
  - [Contact maps](#contact-maps)
    - [Configuration](#configuration-5)
    - [API](#api-4)
    - [Simple usage](#simple-usage-3)
    - [Changing contact criteria](#changing-contact-criteria)
    - [Contact maps for substructures](#contact-maps-for-substructures)
    - [Contacts between two substructures](#contacts-between-two-substructures)
  - [Frequency maps](#frequency-maps)
    - [Configuration](#configuration-6)
    - [API](#api-5)
    - [Simple usage](#simple-usage-4)
  - [Contact Criteria](#contact-criteria)
    - [Pre-defined criteria](#pre-defined-criteria)
    - [Customizing criteria](#customizing-criteria)
  - [Structure comparison](#structure-comparison)
    - [Overfit](#overfit)
    - [Compdesc](#compdesc)
    - [FitDesc](#fitdesc)
  - [Alignments](#alignments)
    - [Configuration](#configuration-7)
    - [API](#api-6)
    - [Loading alignments](#loading-alignments)
    - [Alignment objects](#alignment-objects)
    - [Editing alignments](#editing-alignments)
    - [Alignment comparison](#alignment-comparison)
    - [Saving alignment to a file](#saving-alignment-to-a-file)
    - [Advanced alignment loading](#advanced-alignment-loading)
  - [Integration with PyMOL](#integration-with-pymol)
    - [Integrating PyMOL with PyDesc](#integrating-pymol-with-pydesc)
    - [Configuration](#configuration-8)
    - [API](#api-7)
  - [Geometry](#geometry)

## Configuration

TBD

## Loading structures

Basis for any study one can perform using PyDesc is loading structure.

### Configuration

```python
from pydesc.config import ConfigManager
import pydesc.chemistry
import pydesc.dbhandler

ConfigManager.dbhandler.cachedir
ConfigManager.dbhandler.pdb_handler
ConfigManager.chemistry.solvent
```

* `cachedir` -- path to directory where downloaded structure files are stored.
* `pdb_handler` -- path to local copy of PDB.
* `solvent` -- list of strings (particle names) to be skipped while loading file.

### API

Api for structure loading is available as module `pydesc.api.structure` and contains two
 functions:
```python
from pydesc.api.structure import get_structures
from pydesc.api.structure import get_structures_from_file

structures = get_structures("1no5")

re_loaded_structures = get_structures_from_file("./biodb/pdb/no/1no5.pdb")
```
`get_structures` takes structure id (for PDB, CATH, SCOP or BioUnit) as string and
 returns list of structure, if it was possible to get to the file.
 It first tries to load a local copy from cache dir, then it tries to download it (and
 store in cache dir).
 It returns a list with as many structures, as many models there were in the file.
 It is often just one, but for NMR or trajectories its longer.
 
To learn more about trajectories see also [this section](#trajectories).
 
`get_structure_from_file` also loads list of structure, but from single file, 
 anywhere on the local machine.
 It can be file in PDB or mmCIF format.

Both methods passes `common_converter` flag to `load_structures` method, which was 
described [here](#number-converters).

### Simple usage

To load locally stored structure and be able to read and edit it, simply run:

```python
from pydesc.structure import StructureLoader

path = "tests/data/test_structures/prots_only/2BLL.pdb"
structure_loader = StructureLoader()
with open(path) as file_handler:
    structures = structure_loader.load_structures([file_handler])
structure = structures[0]
assert len(structure)
for atom_set in structure:
    assert atom_set
    print(atom_set)
    for atom in atom_set.atoms:
        assert atom
        print(atom)
```

During initialization of `StructureLoader` additional arguments can be passed to 
 alter the way structures are loaded.
Note that `load_structures` method takes list of file handlers, as it can read from 
 many files at the time (in order to get common mer numbering or when molecules from 
 bio-unit db are loaded).

To get files from a remote database, an additional database handler is needed:

```python
from pydesc.structure import StructureLoader
from pydesc.dbhandler import MetaHandler

handler = MetaHandler()
structure_loader = StructureLoader()
with handler.open("1no5") as files:
    structures = structure_loader.load_structures(files)
assert len(structures)
```

Handler provides method `open`, which takes reference code (in this case -- to PDB), and
 returns config manager, that then used in with-statement returns a list of file 
 handlers to local cache dir, where it downloaded PDB or mmCIF files.

### Database handlers

There are different handlers in PyDesc, facilitating access to different resources:
* `PDBHandler` -- PDB files in PDB database;
* `MMCIFHandler` -- mmCIF files in PDB database;
* `SCOPHandler` -- PDB files in SCOP database;
* `BioUnitHandler` -- bundles of PDB files in PDBBioUnit database;
* `MetaHandler` -- containing all the above.

Each of them has `get_file` method, which accepts reference code to database resource.
This method returns list of open file handlers.
As this approach requires user to deal with each open file handler, there is also 
 `open` method, returning a context manager.
The latter passes all arguments to `get_file`, but returns context manager, ready 
 to use in with-statement.
It will close all open file handlers when leaving with-statement code.

During initialization one can set handler operation mode, which is a list of integers.
Different numbers denote certain behavior:
* in mode 3 handler reads a local cache;
* in mode 2 handler copies file from a local db to local cache (overwrite);
* in mode 1 handler downloads file from a remote db to local cache (overwrite).

Also order of numbers matter, e.g. "[3, 1]" mean that handler will first try to access
 a local copy of a file, then download it, while "[1, 3]" will result it attempt to
 download and overwrite a local copy, but read it anyway if remote access failed.
 "[1]" would mean handler will always try to access a remote database, even if local 
 copy was available.
Mode "2" works only for PDB database (so `PDBHandler` and `MMCIFHandler`), but requires 
 a synchronised copy od PDB and setting `ConfigManager.dbhandler.pdb_handler` path.

```python
from pydesc.dbhandler import MMCIFHandler
from pydesc.dbhandler import MetaHandler
from pydesc.dbhandler import PDBHandler

meta = MetaHandler(mode=[3, 1])
pdb = PDBHandler(mode=[1])
mmcif = MMCIFHandler(mode=[3])
with meta.open("unit://3pi2") as files:
    assert type(files) is list
    assert len(files) > 1
```
In example above we initialize three handlers: one meta handler and two different 
 PDB handlers.
PDB-file handler can only access the remote database.
mmCIF-file handler can only access a local cache.
`MetaHandler` will access a local copy first, but on failure it will try to download 
 one.

In this example we also download files from PDBBioUnit database to show how 
 `MetaHandler` works.
It is contains all other handlers and if one simply passes a code, it will search all 
 databases.
It is possible to force it to search one particular database, as in the example above,
 using search pattern `database-name://code` or in case of PDBBioUnit: 
 `unit://code/unit`.
 Read the documentation to that class for more details.

As much as `MetaHandler` is convenient, for sake of speed it is better to use 
 handler related to the particular databases if it is clear where the data should be 
 taken from.

### File parsers

PyDesc does not provide any structure file parser itself.
Instead, it uses BioPython parsers to read `.pdb` and `.cif` files and MDTraj parsers
 to read different trajectories.
`StructureLoader` can accept:
* BioPython's `PDBParser`
* BioPython's `MMCIFParser`
* PyDesc's `MetaParser`
* any object with method and signature such as:
  `get_structure(code: str, file_path: str) -> Model`
  where `Model` is a BioPython's class.

By default, `StructureLoader` uses `MetaParser`, which can deal with both PDB and mmCIF
 files.
That means `MetaParser` tries both parsers when loading a file, which does not 
 matter when dealing with single file, but might matter in case of large amount of them.
Both parsers contained in `MetaParser` are set to not log any warnings.

```python
from pydesc.parsers import MetaParser
from pydesc.structure import StructureLoader
from Bio.PDB import MMCIFParser

file_parser = MetaParser(QUIET=False)
structure_loader = StructureLoader(parser=file_parser)

mmcif_parser = MMCIFParser(QUIET=False)
mmcif_loader = StructureLoader(parser=mmcif_parser)
```

### Number converters

Each structure comes with object responsible for translating PDB ids to atoms set ids.
 Such an object is called converter.
 By default, each loaded structure comes with a separate converter.
 When loading large amount of structures that share the same numbering, a converter 
 can be shared as well, e.g. for NMR structures or handful of trajectory frames 
 stored as different pdb files.
 To do so just set flag `common_converter` to True.

```python
from pydesc.dbhandler import PDBHandler
from pydesc.structure import StructureLoader

with PDBHandler().open("1DUF") as files:
    structures = StructureLoader().load_structures(files, common_converter=True)
```

## AtomSet - structure building block

`AtomSet` is a base class for objects representing mer of biopolymers and other
 chemical entities that occur during structure or trajectory analysis.
 In Pydesc each structure comprises some instances of `AtomSet` subclasses.
 There are two important ones:
 * `Mer` -- representing parts of biopolymers like residues or nucleotides.
 * `Ligand` -- all other particles.
 Both of them are abstract as well, which means there are other subclasses inheriting
 from them.

![](images/atom_set.png)

PyDesc provides some classes representing sets of atom in full-atom representation.
 Those will be described as an example.
 Other representations require separate implementation.
 Some basic coarse-grain representation are also implemented.
 Browse through `pydesc.chemistry` module to learn more.
 This section covers common features of `AtomSet` subclasses and how to implement own
  subclasses.

### Configuration

TBD

### API

Right now api related to this section is limited to single convenience function for
 full-atom residues:
```python
from pydesc.api.full_atom import calculate_residues_angles
from pydesc.api.structure import get_structures

structure = get_structures("2bll")[0]
mer1 = structure[0]
mer2 = structure[1]
print(mer1.dynamic_features)
print(mer1.angles)
print(mer1.dynamic_features)
calculate_residues_angles(structure)
print(mer2.dynamic_features)
```
Angles are stored in dictionary `dynamic_features`, but calculation is lazy
, performed only when one tries to get `angles` attribute.
It is suboptimal to call the same algorithm on each particular residue, much faster
 method (vectorized) is invoked by function `calculate_residues_angles`, that takes
 structure as argument.
 It sets angles for all full-atom residues in given structure, even if attribute
 `angles` aws not accessed (like in case of `mer2` and all subsequent residues).

### Simple usage (AtomSet trivia)

`AtomSet` instances are iterable.
 Basic iterator iterates over atoms, returns `Atom` instances.
 They store data typical for `.pdb` file and inherit methods from `Coord` class
 , which is basically vectors:
```python
from pydesc.api.structure import get_structures

structure = get_structures("1no5")[0]

first_residue = structure[0]

for atom in first_residue:
    print(atom.element, atom.vector, atom.occupancy, atom.b_factor)

atom1, atom2, *_ = tuple(first_residue)

vector_a1_a2 = atom1 - atom2
distance_a1_a2 = vector_a1_a2.calculate_length()
direction = vector_a1_a2.get_unit_vector()
```
Another way to get atoms from `AtomSet` instance is by its name:
```python
from pydesc.api.structure import get_structures

structure = get_structures("1no5")[0]

first_residue = structure[0]
ca_atom = first_residue.atoms["CA"]
cb_atom = first_residue.atoms["CB"]
```
Beside atoms, `AtomSet` instances can also store pseudoatoms.
Any type will have geometrical center for sure.
```python
from pydesc.api.structure import get_structures

structure = get_structures("1no5")[0]

residue = structure[0]
ion = structure[100]
compound = structure[105]

for atoms in (residue, ion, compound):
    print(atoms.gc)
```
Pseudoatoms are calculated dynamically and cached when accessed first time.
 Each `AtomSet`, beside pseudoatoms available in `pseudoatoms` attribute, stores also
 `dynamic_features`, which are also calculated dynamically and represents all sort of
 features that can be calculated for set of mers, but are not pseudoatoms, like base
 planes or torsion angles.
 Bear that in mind when working with trajectories (changing frames resets cache!).

Distinguish between backbone and side chain atoms for mers leads us to additional
 iterators in `Mer` subclasses:
```python
from pydesc.api.structure import get_structures

structure = get_structures("1no5")[0]
residue = structure[0]

if residue.is_chainable():
    print("Backbone")
    for atom in residue.iter_bb_atoms():
        print(atom)
    print("Side chain:")
    for atom in residue.iter_nbb_atoms():
        print(atom)
    print(residue.next_mer)
    print(residue.prev_mer)
```
`Mer` subclasses also return `True` when one calls `is_chainable` on them, just to
 make it easier.

If PyDesc detects chemical bound between two mers, it sets their attributes `next_mer
` and `prev_mer` (otherwise they are set to `None`).

### AtomSet factories

While parsing a file, instances of `AtomSet` subclasses are produced by factories.
 By changing factory one might influence the way mers or ligands are created, while
 loading structure.

Structure loader parses file, then passes every item that has its own id in pdb or
 cif file to factory.
 By default, factory tries to produce `Residues` or `Nucleotides` out of what it gets
 , then -- `MonoatomicIon` or `Compound`.
 PDB id "2BLL" refers to structure that comprises only residues (no ions, no other
 particles).
 Knowing that one could want to change default behaviour and set factory to produce
 Residues only:

```python
from pydesc.chemistry import factories as as_factories
from pydesc.chemistry.full_atom import Residue
from pydesc.dbhandler import MetaHandler
from pydesc.structure import StructureLoader

classes = [Residue]

as_factory = as_factories.BioPythonAtomSetFactory(classes)

structure_loader = StructureLoader(
    atom_set_factory=as_factory,
)

with MetaHandler().open("2bll") as files:
    models_db = structure_loader.load_structures(files)
```

In this case factory only tries to produce `Residue` instances.
 Structure loader tries to get some kind of representation for every set of atoms
 that has id in file and raises `ValueError` if it is not possible.
 That could be the case if loaded file contained incomplete residues or something
 other than residue.
 Limiting possible types of mers speeds up structure loading, but optimization is not
 the main reason to do so.

Creating default factory is equivalent of: `BioPythonAtomSetFactory([Residue,
 Nucleotide, MonoatomicIon, Compound])`.
 All classes come from `pydesc.chemistry.full_atom` module.
 That setting allows for successful load of almost any structure, but produces
 incorrectly classified objects, e.g. `Compound` instances out of residues that lack
 backbone atoms.
 
That leads us to actual use case for this feature: representations other than full
-atom, for example CaBS, MARTINI or C-alpha- or P-trace require change of that
 setting, otherwise structure loading gives misleading results.
 For example nucleic acids represented as P-trace are loaded as chains full of
 `MonoatomicIon` instances.

### Full-atom representation

By default, PyDesc expects loaded files to store structures in full-atom representation.
 In this representation residues and nucleotides provide some features.

Residues in pydesc store pseudoatom `cbx`, which stands for "extended CB" and is
 calculated by extending vector CA->CB by 1.
 It is useful during structure comparison, as it grants extra sensitivity for
 comparison methods for secondary structure type.
 PyDesc also calculates torsion angles (as dynamic feature).
 Note that calculating them individually for each mer is not efficient, so to get all
 of them rather then values for single residues of interest use
 `calculate_residues_angles_vectorized` from `pydesc.chemistry.full_atom` module or
 even more convenient `pydesc.api.full_atom.calculate_residues_angles` function.
 Another interesting pseudoatom is moving average CA.
 Size of frame is configurable through `ConfigManager.chemistry.residue.moving_average`:
 Residues also store geometrical center of side chain, called `rc`.

```python
from pydesc.api.structure import get_structures

structure = get_structures("2bll")[0]
residue = structure[0]

print(residue.cbx)
print(residue.angles)
print(residue.backbone_average)
```

Nucleotides come with `rc`, `prc`, `nx` and `ring_center` pseudoatoms.
 `rc` is geometrical center of side chain (all atoms except backbone atoms; for
  historical reasons we preserved the name, which stands for "residue center").
 `prc` is center of ring for pyrimidines or center of five-member ring for purines.
 `nx` is point one gets when extending vector along glycosidic bond by 1.4.
 `ring_center` is center of only ring for pyrimidines or center of six-member ring
 for purines.
 Note that for pyrimidines `ring_center` and `prc` is the same pseudoatom, and is
 not calculated twice if accessed from different attributes.
 An additional feature calculated for nucleotides is `ring_plane`, which returns `Plane`
 instance.
 Plane stores vector perpendicular to itself.
 Having two planes enables calculation of dihedral angles or bisection planes.
 See more information about this class [here](#geometry).
```python
from pydesc.api.structure import get_structures

structure = get_structures("1kis")[0]
pyrimidine = structure[6]
purine = structure[7]

print(purine.rc)
print(purine.prc)
print(purine.nx)
print(purine.ring_center)
print(purine.prc == purine.ring_center)
print(pyrimidine.prc == pyrimidine.ring_center)

print(purine.ring_plane)
```

### Martini

For now, it is possible to load protein structures in Martini representation:
```python
from pydesc.chemistry.factories import BioPythonAtomSetFactory
from pydesc.chemistry.martini import MartiniResidue
from pydesc.structure import StructureLoader

factory = BioPythonAtomSetFactory(classes=[MartiniResidue])
loader = StructureLoader(atom_set_factory=factory)
path = "tests/data/test_structures/martini/gpcr_d.pdb"
with open(path) as file_handler:
    stc = loader.load_structures([file_handler])[0]

for residue in stc:
    print(residue, residue.atoms, residue.last_sc)
```
Residues in this representation consist of `BB` pseudoatom representing backbone, and
 some number of `SC<no>` pseudoatoms representing sets of 4 side chain atoms.
 Property `last_cs` provides access to `SC` pseudoatom of the highest index. Glycine and
 alanine has no `SC` atom, so for them this property returns `BB` pseudoatom.

#### Martini configuration

```python
from pydesc.config import ConfigManager
import pydesc.chemistry.martini

ConfigManager.chemistry.martiniresidue.bb_bond_threshold
ConfigManager.chemistry.martiniresidue.backbone_atoms
ConfigManager.chemistry.martiniresidue.indicators
```
* `bb_bond_threshold` -- max distance between BB atoms for mers to be considered
  consecutive
* `backbone_atoms` -- "BB"
* `indicators` -- "BB" + "last_sc"

### Backbone trace (P-trace, CA-trace)

Representation in which only CA for residues or P for nucleotides have very limited 
 usage for descriptor-based functionality, but can be useful in analysing results.
 It is compact in terms of memory usage, therefore great for storage purposes.

Typically, results of multiple alignment are stored in this representation to save 
 memory.
 Original structure is then superimposed onto that CA-trace representation.

CA-trace provides pseudoatom `mpp`, which is calculated as sum of vectors anchored in 
 adjacent CAs, pointing at CA of chosen residue, scaled to have length of 2.53Å.
 It is meant as fast to calculate, preliminary prediction of cbx position.

```python
from pydesc.chemistry.bbtrace import CATrace
from pydesc.chemistry.factories import BioPythonAtomSetFactory
from pydesc.dbhandler import MetaHandler
from pydesc.structure import StructureLoader

factory = BioPythonAtomSetFactory(classes=[CATrace])
loader = StructureLoader(atom_set_factory=factory)
with MetaHandler().open("1KAN") as fh:
    stc, = loader.load_structures(fh)

for ca_residue in stc:
    assert ca_residue.atoms["CA"]
    assert ca_residue.mpp
```

Similarly, for P-trace:

```python
from pydesc.chemistry.bbtrace import PTrace
from pydesc.chemistry.factories import BioPythonAtomSetFactory
from pydesc.chemistry.full_atom import Ligand
from pydesc.chemistry.full_atom import MonoatomicIon
from pydesc.dbhandler import MetaHandler
from pydesc.selection import AtomSetExactType
from pydesc.structure import StructureLoader

factory = BioPythonAtomSetFactory(classes=[PTrace, MonoatomicIon, Ligand])
loader = StructureLoader(atom_set_factory=factory)
with MetaHandler().open("1AGN") as fh:
    stc, = loader.load_structures(fh)

for p_nucleotide in AtomSetExactType(PTrace).create_structure(stc):
    assert p_nucleotide.atoms["P"]
```
There is no additional pseudoatoms for P-trace.

In the example above `MonoatomicIon` and `Ligand` are also included, as structure 
 "1AGN" does contain more than just P-trace.
 This is also a reason why selection is used to filter out ions in assertion code block.

#### BB-trace configuration

```python
from pydesc.config import ConfigManager
import pydesc.chemistry.bbtrace

ConfigManager.chemistry.catrace.bb_bond_threshold
ConfigManager.chemistry.catrace.backbone_atoms
ConfigManager.chemistry.catrace.indicators
ConfigManager.chemistry.catrace.mpp_length
ConfigManager.chemistry.ptrace.bb_bond_threshold
ConfigManager.chemistry.ptrace.backbone_atoms
ConfigManager.chemistry.ptrace.indicators
```
* `bb_bond_threshold` -- max distance between CA or P atoms for mers to be considered
 consecutive
* `backbone_atoms` -- "CA" for "catrace", "P" for "ptrace"
* `indicators` -- backbone atom + `mpp` for "catrace"
* `mpp_length` -- by default set to 2.53 (which is average `cbx` distance from "CA" 
  atom); determines length of vector between "CA" and `mpp` pseudoatom when the 
  latter is created.

### User defined representation

PyDesc does not cover a lot of possible representation of biopolymers, so users
 interested in its functionality and operation on structures stored in representation
 not on available list might be interested in implementing their own subclass.

Let us take a closer look at incremental implementation of residues in Martini
 representation as an example:
```python
from pydesc.chemistry.base import Mer
from pydesc.config import ConfigManager

ConfigManager.chemistry.new_branch("newmartiniresidue")

class NewMartiniResidue(Mer):
    pass
```
Basically that is all that is needed for `AtomSetFactory` to use this class, but it
 is not very interesting (with user defined subclass one cas pass it to factory the
 same way as in examples [here](#atomset-factories)).
 While using it one would notice that those mers do not even form valid chains, they
 do not store information about succession and precession in polymer. This
 functionality is obvious and common for mers, therefore it depends on settings
 rather than attributes or methods:
```python
from pydesc.chemistry.base import Mer
from pydesc.config import ConfigManager

ConfigManager.chemistry.new_branch("newmartiniresidue")
ConfigManager.chemistry.newmartiniresidue.set_default("backbone_atoms", ("BB",))
ConfigManager.chemistry.newmartiniresidue.set_default("bb_bond_threshold", 5.0)

class NewMartiniResidue(Mer):
    pass
```
To learn more about settings, look [here](#configuration).

It is important to name branch storing settings exactly as class in lowercase and
 place it in `chemistry` branch of config manager.
 Settings valid for a class are also valid for all its subclasses (so there would
 be no need for separate settings for `class MartiniSpecialRes(MartiniResidue): pass`
 unless, of course, user wants to have them different).

Setting `backbone_atoms` is a sequence of names of atoms that belong to backbone.
 Order is important, first and last names on that list are taken into account when
 checking if two mers form a chemical bond.
 Threshold for that bond is simple distance stored in `bb_bond_threshold` setting.
 In case that was not enough, it is possible to overwrite method `has_bond` that
 takes single argument: instance of the same class and returns True or False.
 On that basis PyDesc decides if two subsequent (in a file) mers are part of the same
 polymer chain and sets attributes `next_mer` and `prev_mer`.
 With `MartiniResidue` defined above loaded mers would have those attributes set.

There are more settings possibly relevant to chainable atom sets, such as `indicators`,
 `crucial_distances` etc. Learn more about them in the section "configuration" of this
 chapter.

Different representation have different features that might be interesting for users,
 but may require calculation and are not read from file directly.
 PyDesc distinguishes between pseudoatoms (points) and other features (so called
 "dynamic features").
 While implementing new class, it is possible to add algorithms that should be
 performed to get them.

If one stores structures in Martini representation, all points such as `BB` or `SC1`
 will be loaded as atoms, not pseudoatoms.
 For this reason there is point in adding more pseudoatoms to this representation.
 It might be interesting to have quick access to pseudoatom that is most distant from
 backbone.
 That is pseudoatom, but the one read from file, so we implement this as regular
 python property:
```python
from pydesc.chemistry.base import Mer
from pydesc.config import ConfigManager

ConfigManager.chemistry.new_branch("newmartiniresidue")
ConfigManager.chemistry.newmartiniresidue.set_default("backbone_atoms", ("BB",))
ConfigManager.chemistry.newmartiniresidue.set_default("bb_bond_threshold", 5.0)


class NewMartiniResidue(Mer):

    @property
    def last_cs(self):
        self_len = len(self.atoms)
        try:
            return self.atoms[f"SC{self_len - 1}"]
        except KeyError:
            return self.atoms["BB"]
```

To show how to implement pseudoatoms in terms of PyDesc, lets implement geometrical
 center of side chain atoms:
```python
from pydesc.chemistry.base import Mer
from pydesc.chemistry.base import Pseudoatom
from pydesc.chemistry.base import register_pseudoatom
from pydesc.config import ConfigManager

ConfigManager.chemistry.new_branch("newmartiniresidue")
ConfigManager.chemistry.newmartiniresidue.set_default("backbone_atoms", ("BB",))
ConfigManager.chemistry.newmartiniresidue.set_default("bb_bond_threshold", 5.0)

# Note decorator register_pseudoatom instead of property

class NewMartiniResidue(Mer):
        
    @register_pseudoatom
    def scc(self):
        sc_atoms = [name for name in self.atoms if 'SC' in name]
        if not len(sc_atoms):
            return self.atoms["BB"]
        vec = self.atoms[sc_atoms[0]]
        for atom in sc_atoms[1:]:
            vec += self.atoms[atom]
        register_pseudoatom()
        return Pseudoatom(numpy_vec=vec.vector, name="SCC")
```
In this case we use decorator `register_pseudoatom` to indicate that result of
 calculation should be stored in `pseudoatoms` dictionary of instance.
 Each calculation for each instance is performed once, unless the cache is reset.
 That happens when user calls `reset_dynamic_cache` method or when structure is used
 as topology for a trajectory and user changes frame.

In case of other features that requires calculation, but are not represented by
 points (like angles, planes etc.) -- use decorator `register_dynamic_feature` to
 decorate method that returns any object.
 Result will be stored in dictionary similar to `pseudoatoms`, but won't be mixed with
 them.
 Method `reset_dynamic_cache` clears both dictionaries.

## Substructures

Analysis of biomolecule structure often requires picking up only subset of mers present
 in a biopolymer.
PyDesc enables multiple ways of picking up such subsets, depending on their purpose:
 they can be created by slicing structures (described "Simple usage" in this section),
  based on selections (see
 [selections](#selections)) or alignments (see [alignments](#alignments)) or
 calculated based on contact maps (like descriptors, see 
 [this subsection](#descriptors)).

There are many types of substructures in PyDesc.
Some of them are used internally by library, some may be useful to users:
* chains -- subset of mers marked with common chain name.
* segment -- any subset of subsequent chainable mers.
* element -- special subclass of segments of specific length, typically 5.
  Elements are crucial part of local descriptors.
  They are meant for internal usage.
* contacts -- pair of elements making up descriptors.
  They, too, are meant for internal usage.
* local descriptor ("descriptors" for short) -- pivotal type of substructures used by 
  structural alignment methods implemented in PyDesc library.
  Descriptors can be interesting objects to study, but they are implemented for sake of
  aforementioned methods.
* partial structures -- any other subset of mers coming from single structure.

Chains, segments and partial structures will probably be most commonly encountered while
 working with structures in PyDesc.

### Configuration

There is no related configuration.

### API

TBD

### Simple usage

Getting chain is different from getting other structures, but easy, so that is where we
 will begin:

#### Getting chains

```python
from pydesc.api.structure import get_structures

structure, = get_structures("1no5")
for chain in structure.chains:
    print(chain.chain_name, len(chain))
chain_a = structure.get_chain("A")

assert chain_a.chain_name == "A"
assert len(chain_a) == 106
```

Chains contain all atom sets that were marked with particular chain name.
 That includes ligands.
 To get only chainable atom sets, learn about [selections](#selections).
 Names of chains are case-sensitive and can be longer than one letter.

#### Ranges and lists of indices

Method of getting substructure used internally by PyDesc, but also available for users
 is via structure (or substructure) indexing.
 To utilise that, one needs to know mers index numbers.
 Those are assigned to mers in context of the whole structure, thus it differs from
 indexing of python data structures like lists, tuples and arrays.

```python
from pydesc.api.structure import get_structures

structure, = get_structures("1no5")     # two chains: A (106 mers) and B (108 mers)
segment = structure[12: 23]
same_segment = structure[[i for i in range(12, 24)]]
partial_structure = structure[90: 120]
chain_b = structure.get_chain("B")
segment_b = chain_b[106: 207]

single_slice = structure[1:1]
assert segment[-1].ind == 23            # different from python lists
assert type(segment).__name__ == "Segment"
assert type(partial_structure).__name__ == "PartialStructure"
assert type(segment_b) == type(segment)
assert type(single_slice) == type(segment)
```

Indexing `structure[12: 23]` returns a segment, because all mers in that range are 
 residues (so they are chainable) and they are subsequent in chemical chain.
That is not true for range 90-120, which contains residues, ions and complex ligand  
 from chain "A" and residues from chain "B".

This method is used internally by PyDesc.
 It is fast, but might be inconvenient for users, because:
 * It requires mers' indices, not PDB ids.
 * It expects mers' indices IN STRUCTURE, not their index in given substructure.
    Note that `chain_b[106]` returns the first mer in that substructure.

It also includes the "stop" index mer into returned structure.
 It also means that slice starting and ending on the same mer index returns segment
 instead of (maybe) expected single mer in `structure[1:1]`.
 That is, given atom set of index `1` is a mer, and not ligand, in which case partial
 structure would be returned.
 That behaviour is helpful during automatic segment creation, when it assured that 
 slicing returns iterable object containing mers and ligands.

If so, how can one get atom set of known PDB id or of known index in substructure
(e.g. first mer in chain)?
Let us start with the latter.

 To do this, simply convert a substructure of choice to list or tuple:
```python
from pydesc.api.structure import get_structures

structure, = get_structures("1no5")
chain_b = structure.get_chain("B")
first_chain_b_mer = tuple(chain_b)[0]

assert first_chain_b_mer.ind == 106

for atom_set in chain_b:
    print(atom_set)
```
Substructure, even partial structures, are iterable and preserve order of atom sets
 from original structure.
 What is more, they are sets -- single mer will not occur more than once in any 
 substructure.
 Bear that in mind, as it might lead to unexpected behaviour when using iterable as 
 index:

```python
from pydesc.api.structure import get_structures

structure, = get_structures("1no5")
inds = [1, 2, 3, 1]
part = structure[inds]
assert len(part) != len(inds)

single_mer_part = structure[[1]]
assert type(part) == type(single_mer_part)
```
Index `1` occurred twice in given list, so it was reduced.

Similarly to slice referring to only one atom set, iterable with only one element
 returns segment or partial structure, not a single mer.

#### Indexing with PDB ids

That kind of indexing requires parsing PDB ids and using number converter, so it takes 
 more time than using bare indexing.
Each structure and its substructure has PDB indexer attached to it as `pdb_ids` 
 attribute.
It enables all the types of indexing enabled by regular indexing method:
* single id
* list of ids
* slice

And more, because it also accepts wild cards.
In example below all those modes are utilised to get some substructures that we 
 have seen in previous examples and more.
```python

from pydesc.api.structure import get_structures

structure, = get_structures("1no5")

mer_A22 = structure.pdb_ids["A:22"]             # particular mer
chain_b = structure.pdb_ids["B:"]               # chain wild card
mers_no_22 = structure.pdb_ids["22"]            # mer id wild card
segment = structure.pdb_ids["A:18": "A:23"]     # range
chain_a_residues = structure.pdb_ids[: "A:105"] # range with implicit start
partial = structure.pdb_ids[["A:12", "A:14", "B:45"]]

assert len(partial) == 3
```

Using wild cards in slices and lists of PDB ids is not served.

~~~
Note that wild card for chain ends with ":".
Chain names can be numbers too, so this is needed to distinguish between mer and 
chain wild cards.
~~~

### Descriptors

TBD

## Trajectories

PyDesc copes with trajectories using MDTraj library.
Trajectories supposed to be a dynamic version of structures, so it should be possible 
to do with them whatever is possible with structures.

### Configuration

There is no related configuration.

### API

To make it easier to deal with trajectories, API enables some useful functions. Below
 is an example usage of `api.trajectory.freeze_frame` and `api.trajectory.from_frames`.

```python
from pydesc import api
from pydesc.api import trajectory 

from pydesc.structure import StructureLoader
from pydesc.structure import TrajectoryLoader

structure_loader = StructureLoader()
trajectory_loader = TrajectoryLoader()

topology_path = 'tests/data/test_trajectories/topologies/mdm2.pdb'
with open(topology_path) as fh:
    topology = structure_loader.load_structures([fh])[0]  # note [0]

trajectory_path = 'tests/data/test_trajectories/xtc/mdm2_5frames.xtc'
dynamic_structure = trajectory_loader.load_trajectory(trajectory_path, topology)

frames = []
for i in range(5):
    dynamic_structure.set_frame(i)
    frame = api.trajectory.freeze_frame(dynamic_structure)  # convert frame into structure
    frames.append(frame)

new_trajectory = api.trajectory.from_frames(frames)         # convert frames to trajectory
```

`api.trajectory.freeze_frames` turns the current frame into static structure. Frozen 
frames no longer change coords when `set_frame` is called.

`api.trajectory.from_frames` turns list of static structures (of the same structure) 
into a trajectory, so it is a reverse operation. It is meant to deal with NMR 
structures and trajectories stored in pdb files.

### Simple usage

Basically, to work with a trajectory, one needs topology and trajectory file.
Topology is to be loaded with `StructureLoader`, trajectory -- with 
`TrajectoryLoader`, both coming from `structure` submodule.

```python
from pydesc.structure import StructureLoader
from pydesc.structure import TrajectoryLoader

structure_loader = StructureLoader()
trajectory_loader = TrajectoryLoader()

topology_path = 'tests/data/test_trajectories/topologies/mdm2.pdb'
with open(topology_path) as fh:
    topology = structure_loader.load_structures([fh])[0]  # note [0]

trajectory_path = 'tests/data/test_trajectories/xtc/mdm2_5frames.xtc'
trajectory = trajectory_loader.load_trajectory(trajectory_path, topology)
last_frame = trajectory.get_n_frames()  # 5 in this case

residue0 = trajectory[0]            # GLU
segment = trajectory[20:30]         # LYS20-PHE30
position1 = residue0.ca.vector      # 
for residue in segment:
    print(residue.ca.vector)


trajectory.set_frame(3)
current_frame = trajectory.get_frame()

position2 = residue0.ca.vector      # all values shifted
for residue in segment:             # as frame shifted
    print(residue.ca.vector)
```

Trajectory loader takes two arguments: path to trajectory file and structure object
 returned by structure loader.
Trajectory object has all the features structure has, plus methods `set_frame`, 
`get_frame` and `get_n_frames`.

From trajectory one can derive all substructures available in PyDesc like segments or
 contacts, even descriptors to see how they change with time (although the same set 
 of residues in different frames might not be able to form a descriptor any more).
 It is important to understand that changing trajectory frame will make implicit 
 changes to derivatives as well.

To avoid that, it is possible to freeze a frame. See how in "API" section.
<!--
TODO: Add link
-->

## Selections

TBD

## Contact maps

Contact maps stores data about contacts between mers or ligands (we will refer to
 them as "mers" to keep it shorter or "atom sets" as common base class for those
  objects is called `AtomSet`).
Each pair of atom sets therefore have certain contact value, depending on definition
 of contact.
Hence, contact map is a matrix of contacts values associated with loaded structure
 and contact criterion.

### Configuration

There is no related configuration.
See configuration section for Contact Criteria.
<!--
TODO: Add link
-->

### API

API for contact maps is stored in two submodules:
* `pydesc.api.cmaps`
* `pydesc.api.criteria`

See [this section](#pre-defined-criteria) to learn more about `criteria`.

Module `cmaps` helps perform quick calculation of contact maps:

```python
from pydesc.api.cmaps import calculate_contact_map
from pydesc.api.structure import get_structures

structure, = get_structures("1no5")
contact_map = calculate_contact_map(structure)
```

If no criterion was passed, default one is used. To change that simply add argument:

```python
from pydesc.api.cmaps import calculate_contact_map
from pydesc.api.criteria import get_default_protein_criterion
from pydesc.api.structure import get_structures

structure, = get_structures("1no5")
criterion = get_default_protein_criterion()
contact_map = calculate_contact_map(structure, criterion)
```

### Simple usage

To calculate contact map one needs structure and contact criterion:

```python
from io import StringIO

from pydesc.api.structure import get_structures
from pydesc.contacts import ContactMapCalculator
from pydesc.contacts.geometrical import PointsDistanceCriterion

structure, = get_structures("2bll")

criterion = PointsDistanceCriterion("rc", 7.0, 0.5)
cm_calculator = ContactMapCalculator(structure, criterion)
contact_map = cm_calculator.calculate_contact_map()

file_like = StringIO()
with file_like as file_h:
    contact_map.to_string(file_h)

for ids, contact_value in contact_map:
    residue1, residue2 = ids

    if contact_value == 0:
        value_string = " not"
    elif contact_value == 1:
        value_string = " probably"
    else:
        value_string = " "
    print(f"Residue {residue1} and {residue2} are{value_string} in contact.")

cv_12_34 = contact_map.get_contact_value(12, 34)
print(f"Contact value for residues 12 and 34 is: {cv_12_34}.")

m12_contacts = contact_map.get_atom_set_contacts(12)
print("Residue 12 is in contact with:")
for ind, value in m12_contacts:
    print(f"-residue {ind} (value: {value})")
```

Structure "2bll" contains only residues.

`criterion = PointsDistanceCriterion("rc", 7.0, 0.5)` creates simple contact criterion
  based on distance between geometrical centers of side chains. More information in 
  [this section](#contact-criteria). Right now its only important to know that to 
  calculate contact map users need to provide contact criterion, which is an instance of 
  class implementing certain interface.
`cm_calculator = ContactMapCalculator(structure, criterion)` prepares contact map 
  calculation for given structure according to given contact criterion.
`contact_map = cm_calculator.calculate_contact_map()` performs calculation and returns
  ContactMap instance. This object stores calculated data.

`file_like = StringIO()` creates a stream that mimics file handler. Along with
 `with file_like as file_h: contact_map.to_string(file_h)` it could be replaced with
 `with open(<path>) as fh: contact_map.to_string(fh)` in order to save results to a
  file.
 `to_string` method converts pydesc inds to PDB ids before writing them to file-like
  object.

`for ids, contact_value in contact_map: ...` utilises fact that contact maps are
 iterable.
 Iterator returns tuples storing two elements: atom set ids and contact value.
 `ids` is nested tuple with two integers. Contact value is integer 1 or 2.
 Contacts of value 0 are skipped by iterator.

`cv_12_34 = contact_map.get_contact_value(12, 34)` returns value of contact between two
 sets of atoms.
 When having access to AtomSet objects -- their ids are stored as `ind` attribute.
 See [this section](#atomset---structure-building-block) for more information.
 
`m12_contacts = contact_map.get_atom_set_contacts(12)` returns list tuples storing
 ind of residue in contact and value of that contact. Zero-valued contacts are skipped.

### Changing contact criteria

As different types of research might require different types of criteria -- it is
 possible to use other criteria, for example those defined in PyDesc instead of
  default one (see [this section](#pre-defined-criteria) for more information), or
   even combine them to get new, more complex, criteria (see [this section
   ](#combining-existing-criteria-into-complex-ones) for more information).
 It is also possible to change thresholds or other settings (see [this section
 ](#customizing-targets-and-parameters) for more information).
 In case of customized coarse-grained representation or complicated criteria it might
  also make sense to implement new classes (see [this section](#extending-base
  -classes) for more information).

Whatever is the origin of criteria other than default, using them is as simple as
 passing them to `ContactMapCalculator` during initialization:

```python
from pydesc.api.criteria import get_default_protein_criterion
from pydesc.api.structure import get_structures
from pydesc.contacts import ContactMapCalculator

structure, = get_structures("1no5")

criterion = get_default_protein_criterion()
cm_calculator = ContactMapCalculator(structure, contact_criterion=criterion)
contact_map = cm_calculator.calculate_contact_map()
```
In this case we use pre-defined criterion that was designed by our team.
 It calculates contacts between residues only and is efficient in terms of time of
  computation.

### Contact maps for substructures

It is possible to calculate contact map for a substructure alone.
 To do so substructure should be passed as first argument instead of whole structure
  to `ContactMapCalculator` initialization:
```python
from pydesc.api.criteria import get_rc_distance_criterion
from pydesc.api.structure import get_structures
from pydesc.contacts import ContactMapCalculator

structure, = get_structures("1no5")
chain_A = structure.chains[0]

criterion = get_rc_distance_criterion()
cm_calculator = ContactMapCalculator(chain_A, criterion)
chain_A_contact_map = cm_calculator.calculate_contact_map()
```
In this case we use chain, but that could be any substructure: segment, descriptor or
 custom subset created on base of selection. Visit other sections to get more
 information about creating [descriptors](#descriptors), substructures directly from
 other structures ([here](#substructures)) or from [selection](#selections).

### Contacts between two substructures

Instead of intra-structure contacts, one might want to calculate contacts between two
 substructures, e.g. between two different chains. PyDesc enables that as long as two
 substructures comes from the same structure (basically: from the same loaded file).
 To calculate such contact maps one needs to define two selections that can produce
 desired substructures. Learn more about selections from [this section](#selections).

Let us calculate contacts between peptide (chain F) and T-Cell receptor (TRC, chains D
 and E):

```python
from pydesc.api.criteria import get_default_protein_criterion
from pydesc.api.structure import get_structures
from pydesc.contacts import ContactMapCalculator
from pydesc.selection import ChainSelection

structure, = get_structures("4h26")
criterion = get_default_protein_criterion()

peptide = ChainSelection("F")
trc = ChainSelection("D") + ChainSelection("E")

cm_calculator = ContactMapCalculator(structure, criterion, selections=(peptide, trc))
inter_contact_map = cm_calculator.calculate_contact_map()

for i in inter_contact_map:
    print(i)
```

## Frequency maps

Given sequence of contact maps one could be interested in calculating frequency
 maps, i.e. map storing frequencies of contacts instead of their values.

### Configuration

There is no related configuration.

### API

Single function related to frequency maps is stored in `pydesc.api.cmaps` submodule.
It turns sequence of contact maps into frequency map:
```python
from pydesc.api.cmaps import calculate_contact_map
from pydesc.api.cmaps import create_frequency_map_from_contact_maps
from pydesc.api.structure import get_structures

cmaps = []
for structure in get_structures("2ljp"):  # 20 NMR structure
    cmap = calculate_contact_map(structure)
    cmaps.append(cmap)

fmap = create_frequency_map_from_contact_maps(cmaps)
for ids, frequency in fmap:
    print(f"Residues {ids[0]} and {ids[1]} contact frequency: {frequency}.")
```

### Simple usage

Frequency maps are mostly useful when working with trajectories.
 Similarly to contact maps, there is class performing calculation:
```python
from pydesc.contacts.maps import FrequencyMapCalculator
from pydesc.api.trajectory import from_frames
from pydesc.api.structure import get_structures
from pydesc.api.criteria import get_default_protein_criterion

structures = get_structures("2ljp")
trajectory = from_frames(structures)
criterion = get_default_protein_criterion()
fm_calculator = FrequencyMapCalculator(trajectory, criterion)
fmap = fm_calculator.calculate_frequency_map()

frequencies3 = fmap.get_contacts_frequencies(3)
frequency_3_118 = fmap.get_contact_frequency(3, 118)
occurs3 = fmap.get_contacts_occurs(3)
occurs_3_118 = fmap.get_contact_occurs(3, 118)

for inds, frequency in fmap:
    print(f"{inds[0]} - {inds[1]} frequency: {frequency}")
```
Frequency maps, as contact maps, are iterable. Iterator returns tuple of ids 
 and frequency. However there are methods specific to that kind of maps:
* `get_contacts_frequencies` -- takes single id of atom set and returns list of
 results in format (<id of atom set in contact>, <frequency>)
* `get_contact_frequency` -- returns frequency of contact between two atom sets.
* `get_contacts_occurs` -- similar to `get_contacts_frequencies` in terms of
 arguments, but returns number of frames on which contacts were present.
* `get_contact_occurs` -- takes two ids and returns number of frames with contact.

Contacts of value 1 (uncertain) are counted as 0.5, therefore occurs are float.
Since this might be confusing -- consider using binary criteria (returning only 
contact values 0 and 2) with frequency maps. See next section for more information.

## Contact Criteria

In PyDesc contact criteria use three-value logic with possible values being:
* 0 -- no contact
* 1 -- possible contact
* 2 -- sure contact

If it does not make sense to have uncertain value for the criterion of usage -- value of
 1 is dropped.

PyDesc implements several useful criteria.
 They are all accessible through module `pydesc.contacts.criteria`.
 Users can use default criteria or create customized ones.

### Pre-defined criteria

Pre-defined criteria can be obtained by calling functions stored in `pydesc.contacts
.criteria` module.

Right now there are three worth mentioning:
* geometrical center of atom set distance criterion
* geometrical center of side chain distance criterion
* default protein criterion

```python
from pydesc.api.criteria import get_default_protein_criterion
from pydesc.api.criteria import get_gc_distance_criterion
from pydesc.api.criteria import get_rc_distance_criterion

gc_criterion = get_gc_distance_criterion()
rc_criterion = get_rc_distance_criterion()
protein_criterion = get_default_protein_criterion()
```
Gc criterion works with any kind of atom set (even the one designed by user).
 It is fast, but probably works only as heuristics for further, more precise and
 adjusted to type of research or calculation.

Rc criterion works for full-atom mers only.
 It is fast, but not very accurate.
 Default thresholds are:
 * distance up to 7.0Å gives contact value 2
 * distance between 7.0Å and 8.0Å gives 1
 Those values can be changed by passing additional arguments to function creating
  this criterion:
 `threshold` and `margin`. Actual thresholds are calculated as `threshold - margin` and 
 `threshold + margin`.

Default protein criterion is example of complex criterion.
 According to that criterion residues are in contact when:
 * their CAs are less than 5.5Å (value 2) or 6.5Å (value 1) away
 OR
 * their CBXs are less than 6.0Å (value 2) or 7.0Å (value 1) away, and they are
  pointing at each other.
 To make sure that CBXs are pointing at each other we calculate difference between
  two distances: CAs anc CBXs. This distance is expected to be less than 0.75Å with
   margin 0.05Å.

### Customizing criteria

There are different ways of introducing customized criteria to calculation of contact
 maps with PyDesc:
* using existing classes with modified parameters
* setting custom targets (like residues, nucleotides etc.)
* combining existing or customized criteria to generate new ones
* extending existing base classes

#### Customizing targets and parameters

One of the most robust criteria in terms of performance is the one taking into
 account distance between 
 geometrical center of residues. By default, PyDesc uses criterion that was defined
  like this:
```python
from pydesc.contacts.geometrical import PointsDistanceCriterion

criterion = PointsDistanceCriterion("gc", 8.5, 0.5)
```
Usage of this criterion during calculation of contact map will give value 2 to each
 pair of atom sets such that their geometrical centers ("gc" is the name of
  pseudoatom representing geometrical center of set of atoms; more about it
  [here](#atomset---structure-building-block)) are closer than 8.5Å - 0.5Å = 8.0Å.
 Distances in range 8.0Å up to 9.0Å (7.5 + 0.5) will result in contact value of 1
 , which means uncertain contact. Everything above 9.0Å will be marked as 0.

Obvious way of customizing criterion is passing different arguments to its
 initializer, e.g.
```python
from pydesc.contacts.geometrical import PointsDistanceCriterion
from pydesc.chemistry.full_atom import Residue
from pydesc.selection import AtomSetSubclass

residues_selection = AtomSetSubclass(Residue)

criterion = PointsDistanceCriterion("cbx", 4.5, 0.25)
criterion.set_selection(residues_selection)
```
In this case we created criterion similar to one of PyDesc defaults, but this time it
 calculates distances between extended carbon β (cbx), which will give:
 - 2 for residues with cbxs less than 4.25Å away
 - 1 for residues with cbxs 4.25 to 4.75 away
 - 0 for residues with cbxs more than 4.75 away

Since that criterion only makes sens to residues, we added an extra line there -- 
 `criterion.set_selection(residues_selection)`, which tells how the structure will be 
 narrowed down during calculation of this particular criterion.
 In this case -- only residues will be used.

That features give another possibility of criterion customization.
 It allows to specify for which types of atom sets calculation of contacts will be
  performed.
 For example, it enables excluding ligands from calculation for certain criteria by
 passing this selection:
 ```python
from pydesc.chemistry.base import Mer
from pydesc.selection import AtomSetSubclass
from pydesc.api.criteria import get_gc_distance_criterion

selection = AtomSetSubclass(Mer)
criterion = get_gc_distance_criterion()
criterion.set_selection(selection)
```
In this case distance criterion based on distance of geometrical center, by default
 taking into account mer and all types of ligands (as all have `gc` pseudoatom
 ) -- works only for mers (i.e. residues and nucleotides).

For simple criteria the same can be achieved (preferred way) by calculating contact
 map for structure lacking ligands (see [that part](#contact-maps-for-substructures)).
 Setting selections is necessary in case of complex criteria, e.g. when
 overall criterion is alternative of, lets say, aromatic and aliphatic residues-specific
 criteria.
 See [next section](#combining-existing-criteria-into-complex-ones) for more
 information.

Take a closer look at built-in documentation for different submodules of `pydesc
.contacts` to learn more about already implemented classes of contact criteria.

#### Combining existing criteria into complex ones

PyDesc provides complex criteria that enable logical entanglement of criteria.
 Three complex criteria are available:
 * conjunction (and)
 * alternative (or)
 * exclusive disjunction (xor)

For mixed structures it might make sense to use two alternative criteria, i.e. one
 for residues and other one for nucleotides:
```python
from pydesc.contacts.base import ContactsAlternative
from pydesc.contacts.geometrical import PointsDistanceCriterion
from pydesc.selection import AtomSetSubclass
from pydesc.chemistry.full_atom import Nucleotide, Residue

nucleotide_selection = AtomSetSubclass(Nucleotide)
residue_selection = AtomSetSubclass(Residue)

nucleotide_criterion = PointsDistanceCriterion('ring_center', 7.0, 1.0)
nucleotide_criterion.set_selection(nucleotide_selection)

residue_criterion = PointsDistanceCriterion('cbx', 5.0, 1.0)
residue_criterion.set_selection(residue_selection)

final_criterion = ContactsAlternative(residue_criterion, nucleotide_criterion)
```
It is also possible to set selection for `final_criterion`, but in case of complex
 criteria that method sets selections for all sub criteria, so use that with caution
 , as that is not always wanted (as in the case above).

#### Extending base classes

If none of criteria defined in PyDesc is enough, and it is impossible to achieve what
 is needed by combining or customizing settings of existing criteria or result is 
 suboptimal for some reason, it might be a good idea to implement new criterion class.
 To achieve that task it might be necessary to get familiar with some details of
 implementation of AtomSet and contact criteria.
 While this section covers contact criteria, more information about AtomSet is 
 [here](#atomset---structure-building-block).

Easiest way of implementing own criterion class is to extend class `pydesc.contacts
.base.ContactCriterion`, which already provides some methods. In particular, it
 provides methods `calculate_contacts` and `calculate_inter_contacts`. They probably
 should not be overwritten. What they do is narrowing structure down to attached
 selections, then they both call `_fill_contact_matrix` in way specific to type of
 contacts to be calculated (intra- or inter-structural, respectively).
 To extend this subclass one could overwrite `_fill_contact_matrix` (exclusive) or
 `_calculate_contact`. `ContactCriterion` implements former in rather suboptimal
 way: it iterates over two nested loops over AtomSet instances from both selections
 , then calls not implemented method `_calculate_contact`.
 Overwriting `_fill_contact_matrix` is good choice if there is space for
  vectorization, while latter is used when criterion is only easy to calculate when two
 sets of atoms are accessed directly.

Below are two examples:

```python
from pydesc.contacts.base import ContactCriterion
from pydesc.chemistry.full_atom import Residue
from pydesc.selection import AtomSetSubclass

class MyCriterion(ContactCriterion):

    def __init__(self, threshold, margin):
        self.threshold = threshold
        self.margin = margin
        super().__init__()

    def __str__(self):
        return "my criterion"

    def _calculate_contact(self, atom_set1, atom_set2):
        diff = (atom_set1.ca - atom_set2.ca).calculate_length()
        if diff < self.threshold:
            return 2
        elif diff > self.threshold + self.margin:
            return 1
        return 0

criterion = MyCriterion(6.5, 1.0)
criterion.set_selection(AtomSetSubclass(Residue))
```
Signature for this method is `(atom_set1, atom_set2)`, both objects will be AtomSet
 subclass instances.
 It should return single integer equal to 0, 1 or 2. In this example implemented
  criterion has the same functionality as `PointsDistanceCriterion` in PyDesc
  , however interpretation of threshold is different, which might suit some users more.
 
In case of criteria that could be vectorised it is recommended to do so.
 Let us take a look at criterion below.
 `PointsDistanceCriterion` already covers that scope, but it is introduced as example
  easy to understand.
 Such a criterion calculates distance between alpha carbons and has hard coded thresholds:
```python
import numpy
from scipy.spatial.distance import cdist

from pydesc.contacts.base import ContactCriterion
from pydesc.chemistry.full_atom import Residue
from pydesc.selection import AtomSetSubclass

class MyCriterion(ContactCriterion):

    def __init__(self):
        super().__init__()

    def _fill_contact_matrix(self, atom_sets1, atom_sets2, matrix):
        structure_ids1 = [atom_set.ind for atom_set in atom_sets1]
        structure_ids2 = [atom_set.ind for atom_set in atom_sets2]
        points1 = numpy.array([residue.CA.vector for residue in atom_sets1])
        points2 = numpy.array([residue.CA.vector for residue in atom_sets2])
        dist_mtx = cdist(points1, points2)

        possible_contacts = dist_mtx <= 6.0
        points1_indexes, points2_indexes = numpy.where(possible_contacts)
        inds1, inds2 = structure_ids1[points1_indexes], structure_ids2[points2_indexes]
        matrix[inds1, inds2] = 1
        sure_contacts = dist_mtx <= 4.5
        points1_indexes, points2_indexes = numpy.where(sure_contacts)
        inds1, inds2 = structure_ids1[points1_indexes], structure_ids2[points2_indexes]
        matrix[inds1, inds2] = 2

        return matrix

criterion = MyCriterion()
criterion.set_selection(AtomSetSubclass(Residue))
```
Arguments `atom_sets1` and `atom_sets2` are sequences of residues, it should be easy to
 turn them into a vector or array of coordinates or some other features.
 `matrix` is sparse, square matrix that should be changed and returned.
 Each row and column corresponds with the same atom set's ind coming from original
  structure.
 In other words, for each set of atoms, its ind (attribute `ind`) is both index of
  column and row that should show all its contact values.

In this case first block of code in fast, vectorized calculation of distances between
 CA atoms, while second sets appropriate contact values in a matrix representing
  contact map.
 This particular class is useless as it is much more convenient to use
  `PointsDistanceCriterion` like this:
```python
from pydesc.contacts.geometrical import PointsDistanceCriterion
from pydesc.chemistry.full_atom import Residue
from pydesc.selection import AtomSetSubclass

criterion = PointsDistanceCriterion("CA", 5.25, 0.75)
criterion.set_selection(AtomSetSubclass(Residue))
```
It was introduced to help understand the concept.

## Structure comparison

### Overfit

TBD

### Compdesc

TBD

### FitDesc

TBD

## Alignments

Alignment is a mapping between residues or nucleotides of two or more biomolecules.
It could be represented as graph, where vertices are residues or nucleotides (mers)
 and where edges indicate that mer from one structure corresponds with mer from other
 structure(s).
In proper alignment there should be as many subsets of vertices as many structures
 are aligned.
 Edges would be allowed only between different subsets of vertices and every pair of
 structures would be a bijection (with edges representing relation).
 In multiple alignment (alignment of more than two structures) that relation should be
 reflexive, meaning that if mer `a` from structure `A` is aligned with mer `b` from
 structure `B` and mer `b` from structure `B` was aligned with mer `c` from structure 
 `C`, that would lead us to conclusion that mer `a` from `A` is aligned with `c` from
 `C`.
 That allows for two kind of improper alignments: ones that align mers within 
 structures, and those that align single mer with more than one mer from other 
 structure.
 First kind in not served by PyDesc at all.
 Second kind might be dealt with using PyDesc.
 Intention was to enable editing in order to fix such improper alignments.
 Bear in mind that there are methods that will result in unexpected behaviour while 
 called for improper alignments though.
In PyDesc we assume that an edge between two mers from different structures indicates
 structural similarity.
It is possible to use PyDesc to work with aligned sequences as well, but note that this
 might result in unexpected behaviour of some methods.

PyDesc is able to read, write and calculate alignments.
Formats available for reading:
* FASTA
* CSV/TSV
* PAL v1.0
* PAL v2.0

Formats available for writing:
* FASTA
* CSV/TSV
* PAL v2.0

### Configuration

```python
from pydesc.config import ConfigManager
import pydesc.chemistry

ConfigManager.chemistry.mer.code
ConfigManager.chemistry.mer.additional_code
```

This configuration is going to be needed only, when saving files to CSV or FASTA files.
* `code` -- dict mapping (up-to-)3-letter names of mers to 1-letter code.
* `additional_code` -- same dict storing extended mapping.

Those are default dicts.
 Custom representation might use different dicts for this mapping.
 Usually they are stored in `ConfigManager.chemistry.<representation-name>.code`.

### API

There are four convenience functions available in api module.

Two of them are related to loading alignments:
* `get_loader` -- returns loader bonded to given alignment path.
 To learn more about using loaders see [simple](#loading-alignments) and 
 [advanced](#advanced-alignment-loading) examples.
* `load_alignment` -- load alignment from given file, pair or multiple, using
 appropriate loader. Choice of loader depends on file extension.
 For unknown extension this function tries to use FASTA loader.

```python
from pydesc.api.alignment import get_loader
from pydesc.api.alignment import load_alignment
from pydesc.api.structure import get_structures_from_file

alignment_path = 'tests/data/test_alignments/pal/sars_pair.pal'
loader = get_loader(alignment_path)

pdb_path = "tests/data/test_alignments/structures/sars_pair.pdb"
aligned_structures = get_structures_from_file(pdb_path)

alignment = load_alignment(alignment_path, aligned_structures)

same_alignment = loader.load_alignment(aligned_structures)

assert alignment == same_alignment
```
~~~
Note that in the example above two aligned structures are stored in a single pdb file.
~~~
Two other functions are related to processing alignments.
Both take alignment as an argument and return dictionary, where structures aligned in
 given alignment are keys:
* `get_selections` -- stores corresponding selections as values (see 
 [this section](#selections) to learn more).
* `get_partial_structures` -- stores partial structures as values (see
 [this section](#substructures) to learn more).
 Note that this function will return partial structures even if aligned part is a
 segment, chain or whole structure.
 A returned object is treated just as container for atom sets with an interface of a 
 structure.

```python
from pydesc.api.alignment import get_partial_structures
from pydesc.api.alignment import get_selections
from pydesc.api.alignment import load_alignment
from pydesc.api.structure import get_structures_from_file

pdb_path = "tests/data/test_alignments/structures/sars_pair.pdb"
structures = get_structures_from_file(pdb_path)
alignment_path = 'tests/data/test_alignments/pal/sars_pair.pal'
alignment = load_alignment(alignment_path, structures)

for structure in structures:
    print(structure.name, len(structure))

selection_dct = get_selections(alignment)
for k, v in selection_dct.items():
    print(k, v)

substructures_dct = get_partial_structures(alignment)
for k, v in substructures_dct.items():
    print(k, v)
```

### Loading alignments

Loading alignments is done by alignment loader stored in `pydesc.alignments.loaders`
 submodule.
There is one loader for each format: FASTA, CSV/TSV and PAL, and each accepts slightly 
 different arguments, but all need path to a file as first positional argument:

```python
from pydesc.alignment.loaders import CSVLoader
from pydesc.alignment.loaders import FASTALoader
from pydesc.alignment.loaders import PALLoader
from pydesc.api.structure import get_structures_from_file

pdb_path = "tests/data/test_alignments/structures/sars_pair.pdb"
structures = get_structures_from_file(pdb_path)

csv_path = "tests/data/test_alignments/csv/sars_pair.csv"
with open(csv_path) as file_handler:
    csv_loader = CSVLoader(file_handler)
csv_alignment = csv_loader.load_alignment(structures)
print(csv_alignment)

fasta_path = "tests/data/test_alignments/fasta/sars_pair.fasta"
with open(fasta_path) as file_handler:
    fasta_loader = FASTALoader(file_handler)
fasta_alignment = fasta_loader.load_alignment(structures)
print(fasta_alignment)

pal_path = "tests/data/test_alignments/pal/sars_pair.pal"
with open(pal_path) as file_handler:
    pal_loader = PALLoader(file_handler)
pal_alignment = pal_loader.load_alignment(structures)
print(pal_alignment)
```
Example above shows how to load the same alignment of structures of some protein 
 present in SARS1 and SARS2.
Each file stores the same alignment in a different format.

CSVLoader accepts additional argument, `delimiter` set to "\t" by default, which
 makes it TSV loader in fact.
 Passing other delimiter enables reading alignments separated with other character,
 e.g. comma.
 Note that setting delimiter to ":" might cause unexpected behaviour, as this format
 uses colons to separate fields in PDB ids.

FASTALoader accepts `miss_match_characters`, which by default is set to ".-".
 That argument tells loader which characters denote missing mer.
 It could be any sequence of characters.

PALLoader accepts no additional arguments.

#### Sample files

Sample TSV file looks like this:
```csv
MOL1	MOL2
B:382:E	E:335:L
B:383:C	E:336:C
B:384:D	E:337:P
B:385:F	E:338:F
B:386:S	E:339:G
```
Note that this is just 5 first lines.

Sample FASTA alignment looks like this:
```fasta
>MOL1 [B:382-485, B:496-502, B:516-524, B:558-585]
ECDFSPLLSGTPPQVYNFKRLVFTNCNYNLTKLLSLFSVNDFTCSQISPAAIASNCYSSLILD
YFSYPLSMKSDLSVSSAGPISQFNYKQSFSNPTCLILATVPKYSYINKQLVNANQYSGSTVAM
TEQLQMGFGITVQYGTDTNSVC
>MOL2 [E:335-343, E:345-439, E:450-456, E:459-467, E:492-499, E:506-525]
LCPFGEVFNTRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYAD
SFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNYLYRLFSNLKPFERDLQSYGF
QPQPYRVVVLSFELLHAPATVC
```
Ranges in "[]" might be omitted if whole structure was aligned, otherwise loader 
 will be confused while mapping sequence to mers.

Sample PAL file looks like this:
```pal
v:2.0
2
MOL1
MOL2
>MOL1 MOL2
@3
382 -- 390 <--> 335 -- 343	(0,1)
391 -- 485 <--> 345 -- 439	(0,1)
496 -- 502 <--> 450 -- 456	(0,1)
```
This format is used by PyDesc, DEDAL and DAMA.
It was meant to be easy to read by both machines and humans and present all pair 
 alignments in multiple alignment.
If top line "v:x.x" is missing, that indicates version 1.0.
Second line tells how many structures were aligned and is followed by structure
 labels (one per line).
This format consists of sections, each showing single pair alignment.
 Each section starts with line containing ">" and labels of alignment structures.
 Following line tells how many aligned ranges there are in pair alignment ("@3" in
 sample above).
 Each range shows PDB id of segment start and end, possibly with a chain and insertion
 code:
 `A:382 -- A:390i <--> 335 -- 343`
 or in version 1.0:
 `A382 -- A390i <--> 335 -- 343`

Version 2.0 allows for chain names longer than one character (at cost of accepting 
 colon in chain name).

### Alignment objects

Loaders method `load_alignment` returns either `PairAlignment` or `MultipleAlignment`.
 They both implement common interface, although Multiple alignment extends it slightly.
 Main difference is details of implementation of some methods.

Most importantly, alignments store PyDesc inds of mer objects in arrays.
 Each alignment has `structures` attribute, so it is always possible to tell alignment 
 of which structures an object represents.

There are two iterators: one that iterates over rows, and other iterating over columns.
 Each column corresponds with single structure (line in CSV files).
 Each row stores inds of aligned mers from subsequent structures.

It is also possible to crop alignment using standard notation with "[]".
 In that case new alignment object is created.
 Method allows for single index, slice or sequence of indices.
 Arrays are indexed from 0 and indices are relative to first aligned row, not mers
 in any structure, however such approach is also possible (see below).
```python
from pydesc.alignment.loaders import PALLoader
from pydesc.api.structure import get_structures_from_file

pdb_path = "tests/data/test_alignments/structures/sars_pair.pdb"
structures = get_structures_from_file(pdb_path)

pal_path = "tests/data/test_alignments/pal/sars_pair.pal"
with open(pal_path) as file_handler:
    loader = PALLoader(file_handler)
alignment = loader.load_alignment(structures)

left_stc, right_stc = alignment.structures   # this attribute should be identical to 
                                             # "structures" variable
for index, row in enumerate(alignment.iter_rows()):
    left_ind, right_ind = row
    left_mer_name = left_stc[left_ind].name
    right_mer_name = right_stc[right_ind].name
    print(f"{left_mer_name} aligned with {right_mer_name}")

for structure, aligned_inds in zip(alignment.structures, alignment.iter_columns()):
    structure_id = id(structure)
    print(f"Inds aligned for structure {structure_id}:")
    print(",".join([str(i) for i in aligned_inds]))

first_5_rows = alignment[0:5]
also_those_rows = alignment[[0, 1, 2, 3, 4]]
for row1, row2 in zip(first_5_rows.iter_rows(), also_those_rows.iter_rows()):
    assert (row1 == row2).all()     # check if all inds are the same
assert len(first_5_rows) == len(also_those_rows) == 5
```

Given two alignment it is possible to determine what common structures they align with 
 method `get_common_structures`, which returns set.
 See first example [here](#additional-methods) to see it in action.

### Editing alignments
~~~
Note that all methods presented here return new alignment, none works in place.
~~~
In most cases array index-based cropping of alignments is not very useful.
 It is much more common to identify rows in alignments with mers from certain structures
 they relate to.
```python
from pydesc.alignment.loaders import PALLoader
from pydesc.api.structure import get_structures_from_file

pdb_path = "tests/data/test_alignments/structures/sars_pair.pdb"
structures = get_structures_from_file(pdb_path)
sars1, sars2 = structures

pal_path = "tests/data/test_alignments/pal/sars_pair.pal"
with open(pal_path) as file_handler:
    loader = PALLoader(file_handler)
alignment = loader.load_alignment(structures)

helix = sars1[0:9]
helix_inds = [mer.ind for mer in helix]     # they are 0-9 obviously
helix_alignment = alignment.extract_aligned_with(sars1, helix_inds)
assert len(helix_alignment) == 10

strand = sars1[18:23]
strand_inds = [mer.ind for mer in strand]
strand_alignment = alignment.extract_aligned_with(sars1, strand_inds)
assert len(strand_alignment) == 6

partial_alignment = strand_alignment.sum_rows(helix_alignment)
assert len(partial_alignment) == 16

sorted_partial_alignment = partial_alignment.sort()
for row in sorted_partial_alignment.iter_rows():
    print(row)
```
In this example we assume we know location of interesting helix and strand in one of
 aligned proteins.
 We pass indices of mer making up those substructures and are ready to extract rows 
 aligning only residues of interest using `extract_aligned_with` method.
 It takes two arguments:
 * structure for which we know indices of mers we want to extract
 * sequence of indices

It could be done in one step with `helix_inds + strand_inds` passed as second argument.
 In this example however we also wanted to present method `sum_rows`.
 This method takes other alignment and extracts from it columns that contain structures
 present in first alignment, then glues both alignments together, creating a new one.
 That means that:
 * second alignment align more structures than first one, but it will be narrowed
  down in the process;
 * order of structures in second alignment does not matter, columns will be ordered
  accordingly to first alignment.

Since that (and many other) operations can mess up order of rows, method `sort` is
 provided to deal with that problem.
 It orders rows in ascending order of mer indices for longest aligned structure, then 
 tries to fit remaining rows (with gaps for longest structure) in the alignment so
 that they succeed rows aligning mer for other structures with lower mer indices.
 This approach fails for heavily gaped alignment or ones with cyclic permutations
  -- those might need some manual adjustments.

Working with (initially) pair alignment, which are generally well-formed usually
 presents no difficulties.
Let us take a look at another example, involving multiple alignments.
 In case of those one can be interested in narrowing down not only rows, but also
 columns (removing structures from alignment).
 Editing those alignments can lead to creation of pair alignment with gaps or multiple
 alignments with "empty" columns (containing only gaps).
 There are methods to deal with both of those problems.

```python
from pydesc.alignment.loaders import PALLoader
from pydesc.api.structure import get_structures_from_file

pdb_path = "tests/data/test_alignments/structures/kinases5_dama.pdb"
structures = get_structures_from_file(pdb_path)     # 1LUF, 2SRC, 1BO1, 1CDK, 1IA9

pal_path = "tests/data/test_alignments/pal/kinases5_dama.pal"
with open(pal_path) as file_handler:
    loader = PALLoader(file_handler)
ory_alignment = loader.load_alignment(structures)

alignment = ory_alignment.close()
print(f"Original alignment length={len(ory_alignment)}")
print(f"Closed alignment length={len(alignment)}")

structure_1cdk = structures[3]

inds_range = list(range(191, 210))
three_rows = alignment.extract_aligned_with(structure_1cdk, inds_range)
    # 1BO1 and 1IA9 lack any aligned mers here

for row in three_rows.iter_rows():
    print(row)

fixed_three_rows = three_rows.drop_empty_structures()
assert len(fixed_three_rows.structures) == 3

pair_alignment = alignment.limit_to_structures(structures[3:])
fixed_pair_alignment = pair_alignment.prune()
assert len(fixed_pair_alignment) < len(pair_alignment)
```
Alignment read from PAL format are stored as set of pair alignment glued together.
 Even if they are proper alignments (look [here](#alignments) to learn about criteria),
 at this stage there might be several rows aligning the same mer from single structure 
 with different mers in other structures.
 That alone does not mean the alignment is not proper one -- not unless one is unable
 to close it.
 Method `close` returns new alignment.
 If in original one there were rows that might have been merged into one row, in new
 one they will appear as one row.
 For instance in example above there were rows:
```csv
1LUF	2SRC	1CDK
184	347	---
---	347	199
184	---	199
```
which were turned into single row:
```csv
1LUF	2SRC	1CDK
184	347	199
```
...and so on for all rows.
In most cases, right after reading PAL file with multiple alignment -- one wants to
 close that alignment (not always though, and it is impossible to perform opposite
 operation, therefore it was left for a programmer to decide).
 FASTA and TSV/CSV files usually store already closed alignments, but if it is
 otherwise, the very same method should be used to close them.
 Note that this method calculates transitions to add rows if they are missing.
 For example if there was no third row aligning mer no. 184 from `1LUF` with mer no.
 199 from `1CDK`, it would be added on-the-flight.

In example above there are 5 structures aligned.
 It is possible to pick up a set of rows that contains only gaps for certain
 structures, as shown above, where a part of alignment gives "empty" columns for two
 structures.
 Pruning such columns is done by calling `drop_empty_structures` method.

To select only some structures (columns), one could use `limit_to_structures` method, 
 that accepts sequence of structures as argument.
 Structures passed as argument have to be subset of actually aligned structures.
 This method returns pair or multiple alignment, but it is impossible to select only
 single column.
 To get indices aligned in specific single structure, try `api` for alignments 
 ([this section](#api-6)) or lear how to iterate columns [here](#simple-usage-5)

Since multiple alignments contain gaps, it is possible create from them alignments
 that contain rows aligning only one mer or no mers at all.
 That might be desired, otherwise `prune` method drops such rows.
 In example above it was called on pair alignment, but it works for multiple
 alignments as well (as the only method presented in this example).

There is also method `set_columns_order` not shown in example above.
 It accepts single argument -- list of all structures aligned in object for which method
 was called.
 It then returns a new alignment with reordered columns.

#### Additional methods

`PairAlignment` implements some additional methods.
 One of them is `transit`.
 It calculates alignment between structures `B` and `C`, when `A-B` and `B-C` alignments
 were given.
 Here is an example utilising parts of multiple alignment known from paragraphs above:

```python
from pydesc.api.structure import get_structures_from_file
from pydesc.api.alignment import load_alignment

pdb_path = "tests/data/test_alignments/structures/kinases5_dama.pdb"
structures = get_structures_from_file(pdb_path)     # 1LUF, 2SRC, 1BO1, 1CDK, 1IA9

pal_path = "tests/data/test_alignments/pal/kinases5_dama.pal"
ory_alignment = load_alignment(pal_path, structures).close()
pair_alignment1 = ory_alignment.limit_to_structures(structures[:2])
pair_alignment2 = ory_alignment.limit_to_structures(structures[1:3])

common_structure = pair_alignment1.get_common_structures(pair_alignment2)
assert len(common_structure) == 1
common_structure = max(common_structure)
assert common_structure is structures[1]

transition_alignment = pair_alignment1.transit(pair_alignment2)
two_structures = transition_alignment.structures
pair_alignment3 = ory_alignment.limit_to_structures(two_structures).prune()

assert set(transition_alignment.structures) == set(pair_alignment3.structures)
assert transition_alignment.is_consistent_with(pair_alignment3)
assert pair_alignment3.is_consistent_with(transition_alignment)
print(transition_alignment)
print(pair_alignment3)
assert transition_alignment.is_internally_consistent()
```

Calculated transition alignment is 120 long, which is different from alignment of the
 same structures read from file.
 We will take a look at differences in the next chapter, but this example shows how to
 calculate new alignment from two already existing ones.
 It also shows that this new alignment is consistent with actual pair alignment read 
 from a file, thanks to `is_consistent_with` method specific to pair alignments.
 
This method returns boolean values.
 It always returns ture if given objects align different structures.
 Otherwise, it might take a while to calculate value.
 In such case, algorithm checks if there is any mer aligned with different mers in both 
 alignments and returns false if it finds one.
 For internally inconsistent pair alignments (those having single mer from one 
 structure aligned with several mers from other) -- only one (last) alignment of
 inconsistent mer is taken into account.
 That might lead to unexpected behaviour, so it is a good idea to check internal 
 consistency in case o surprising results.

To do so simply call `is_internally_consistent`, as shown in example above.
 It might be needed to manually remove some rows from alignment.
 Alignment are immutable, so that means picking subsets of rows (which results in new
 (sub)alignments and merging them with `sum_rows`.

### Alignment comparison

Analysing alignments requires understanding alignment objects and there are no methods 
 for this specific purpose.
 We just need to use all we know about alignments to extract information we want.
 In previous example we had an alignment of structures 1LUF, 2SRC and 1BO1.
 First we extracted two pair alignment from it: 1LUF-2SRC and 2SRC-1BO1.
 Then we calculated transition to gain 1LUF-1BO1 and noticed that it was different from
 pair alignment 1LUF-1BO1 extracted from original multiple alignment (different, yet 
 consistent).
 It is easy to understand where those differences come from: not all rows from alignment
 1LUF-1BO1 could be inferred from transition.
 Let us make sure it was the case.
 First part of the example is the same chain of operations:

```python
from pydesc.api.structure import get_structures_from_file
from pydesc.api.alignment import load_alignment
from pydesc.alignment.base import DASH

pdb_path = "tests/data/test_alignments/structures/kinases5_dama.pdb"
structures = get_structures_from_file(pdb_path)     # 1LUF, 2SRC, 1BO1, 1CDK, 1IA9

pal_path = "tests/data/test_alignments/pal/kinases5_dama.pal"
ory_alignment = load_alignment(pal_path, structures)
ory_alignment = ory_alignment.limit_to_structures(structures[:3]).close()
    # we just need 1LUF, 2SRC, 1BO1
pair_alignment1 = ory_alignment.limit_to_structures(structures[:2])
pair_alignment2 = ory_alignment.limit_to_structures(structures[1:3])
transition_alignment = pair_alignment1.transit(pair_alignment2)
two_structures = transition_alignment.structures
pair_alignment3 = ory_alignment.limit_to_structures(two_structures).prune()

    # extract different rows
excluded_rows = []
for row in pair_alignment3.iter_rows():
    if row not in transition_alignment.inds:
        excluded_rows.append(row)
excluded_mers_stc0 = [row[0] for row in excluded_rows]
    # get ids of residues from 1LUF aligned in excluded rows
stc0 = pair_alignment3.structures[0]    # 1LUF
    # get only excluded rows from multiple alignment
    # where middle structures is 2SRC
difference = ory_alignment.extract_aligned_with(stc0, excluded_mers_stc0)
for row in difference.iter_rows():
    assert row[1] is DASH
    print(row)
```

Now it is clear that each of the rows not included in transition alignment contains
 alignment of residues from 1LUF and 1BO1, which were not aligned with any residue
 in 2SRC.

Alignments do not have `__eq__`, `__lt__` etc. methods defined (with one exception),
 because it would be ambiguous if answer concerns space of structures or rows.
 They do have `__len__` though, which tells how many rows there are in an object.

To compare structures -- take a look at `structures` attribute.
 Their order corresponds with order of columns in `inds` attribute.

To compare rows one needs to narrow down structure space (with `limit_to_structures`)
 and order alignments in the same way (with `set_columns_order`), then possibly performs
 operations on `inds` arrays.

Exception mentioned above regards `PairAlignment` class, which implements `__eq__`.
 Equality between two pair alignments occurs when two objects are the same, OR 
 they are different, but align the same structures (object equality) and all rows are
 the same.
 `__eq__` methods deals fine with alignment even if order of columns is different.

### Saving alignment to a file

To save an alignment to a file one needs alignment, saver and file handler, as in the
 example below:

```python
from io import StringIO

from pydesc.alignment.savers import CSVSaver
from pydesc.alignment.savers import FASTASaver
from pydesc.alignment.savers import PALSaver
from pydesc.api.alignment import load_alignment
from pydesc.api.structure import get_structures_from_file
from pydesc.chemistry import ConfigManager

# in this case we load alignment, but it could be calculated or crafted as well...
pdb_path = "tests/data/test_alignments/structures/kinases5_dama.pdb"
structures = get_structures_from_file(pdb_path)
pal_path = "tests/data/test_alignments/pal/kinases5_dama.pal"
alignment = load_alignment(pal_path, structures).close()

# ...here we could do some editing.
# We assume 'alignment' is an object to be saved.

mocked_pal_file = StringIO()
pal_saver = PALSaver()
pal_saver.save(mocked_pal_file, alignment)

residue_3_to_1_map = ConfigManager.chemistry.mer.code

mocked_csv_file = StringIO()
csv_saver = CSVSaver(delimiter=",", sequence_dict=residue_3_to_1_map)
csv_saver.save(mocked_csv_file, alignment)

structure_iterator = enumerate(alignment.structures)
better_names = {structure: "stc-%i" % i for i, structure in structure_iterator}

mocked_fasta_file = StringIO()
fasta_saver = FASTASaver(wrap_at=60, sequence_dict=residue_3_to_1_map)
fasta_saver.save(mocked_fasta_file, alignment, names=better_names)

for mocked_file in mocked_pal_file, mocked_csv_file, mocked_fasta_file:
    mocked_file.seek(0)
    file_content = mocked_file.read()
    print(file_content)
    assert len(len(file_content))
```

Assuming one has an alignment object, next thing is to get and create saver object.
 There are different savers for different formats.
 None of them requires any arguments, but they might be useful.
 They are described below.
 Each saver has method `save`, which takes file-like object and alignment to be saved.
 File-like object can be file handler, so in the example above we could replace:
 ```
mocked_pal_file = StringIO()
pal_saver = PALSaver()
pal_saver.save(mocked_pal_file, alignment)
 ```
with
 ```
with open("test.pal", "w") as file_handler:
    PALSaver().save(file_handler, alignment)
 ```
.

Each saver has `save` method and each of them accepts named key argument `names`,
 as shown on fasta saver example.
 That means this argument is optional, but if passed, it has to be explicitly named
 (with `names=...`).
 This argument is expected to be a dict mapping structures (from given alignment) to
 their names.
 When it is not provided -- names from structure objects are taken (`name` attribute).

`PALSaver` does not take additional arguments during initialization.
 For now, it always saves pal version 2.0.
~~~
 PAL files are bonded to their structure files.
 Trying load pal alignment with different structures (with shifted mer indices) will
 lead to errors.
~~~

`CSVSaver` takes optional arguments:
 * `delimiter` -- string to put between mers in rows (default is tabulation);
 * `dash` -- string to put instead of dash (default is "-");
 * `sequence_dict` -- dict mapping 3-letter names of mers to 1-letter code.
If `sequence_dict` is not given, names are translated to 1-letter code basing on their
 type.
 That approach might lead to unexpected errors when given structures are flawed (e.g.
 contain incomplete mers that were loaded as ions or ligands).
 Instead of fixing the structures, which might not be the most important thing at the
 time, it is possible to simply pass a universal mapping for all atom sets present in
 the structures.
 That ugly workaround might save some time and only exists for formats that require
 sequence (CSV and FASTA).

`FASTASaver` takes:
 * `wrap_at` -- integer or `None`, number of characters to fit in single line.
   If (default) `None` is given, whole sequence will be placed in single line.
 * `dash` -- string to put instead of dash (default is "-");
 * `sequence_dict` -- dict mapping 3-letter names of mers to 1-letter code.
    See description of corresponding argument in `CSVSaver`.

~~~
FASTA anc CSV files are bonded to biomolecule sequence,
therefore might require additional sequence mapping or carefully loaded structures.
~~~

### Advanced alignment loading

When dealing with multiple alignment, especially in multiple processes, on clusters or
 very elaborated scripts, it might be necessary to have a bit more control over 
 alignment loading.
 With loaders in PyDesc it is possible to read some metadata from alignment before 
 creating alignment object.
 It is also possible to load only part of alignment.

Loader is bonded to a file, and (in a current version) it reads the file when it is 
 created.
 Right at start it tries to read metadata, so it will raise an error if it could not do
 this.
 Kind of error and exact type of file structure required at this stage depends on 
 file format, thus -- loader class.
 E.g. `PALLoader` requires number of structures in 1st or 2nd line.

After loader object was created, user can call `read_metadata` method.
 It returns a dictionary.
 Its structure really depends on a file format, but it will contain 'labels' key.
 Value for this key will be list of strings -- names of structures aligned in the 
 alignment file.
 When calling `load_alignment` and passing list of structures, they will be mapped on
 that list of labels (preserving order).
 Other metadata depends on file format, e.g. FASTA stores ranges as metadata.

In example below we have single PDB file with 5 structures and some alignment of three
 of them.
 Unfortunately, alignment file calls structures "MOL1", "MOL2" and so on.
 Assuming that we work with larger number of alignment files,
 each aligning from 2 to 5 of those structures and refer to them as "MOLx",
 where x could be integer from range 1-5, we could read metadata before loading
 alignment to pass correct structures.

```python
from pydesc.api.structure import get_structures_from_file
from pydesc.alignment.loaders import PALLoader

pdb_path = "tests/data/test_alignments/structures/kinases5_dama.pdb"
structures = get_structures_from_file(pdb_path)     
names = "1LUF", "2SRC", "1BO1", "1CDK", "1IA9"
structures_dct = {}
for ind, (name, structure) in enumerate(zip(names, structures), 1):
    structure.name = name
    mol_no = "MOL%i" % ind
    structures_dct[mol_no] = structure
    
pal_path = "tests/data/test_alignments/pal/kinases42.pal"
with open(pal_path) as file_handler:
    loader = PALLoader(file_handler)
metadata = loader.read_metadata()
alignment_structures = []
for label in metadata['labels']:
    alignment_structures.append(structures_dct[label])

alignment = loader.load_alignment(alignment_structures)

    # reading only pair alignment from the same file
local_name_map = {
    "MOL4": structures[3],
    "MOL5": structures[4],
}
pair_alignment = loader.load_alignment_mapping(local_name_map)
```

Below comment in example above, we use method `load_alignment_mapping` and pass the 
 mapping, to load only part of alignment.
 In that case we read pair alignment from a file that contains alignment of three 
 structures.
 In combination with `read_metadata`, preparing a local map should be possible.
 It allows loading only alignments of interest, without need to loading whole 
 structures just to load a multiple alignment file to discard them later.

## Integration with PyMOL

PyDesc was meant to be used in scripts and programs, but it can also help in EDA.
 Sometimes opportunity to take a look at structures we work with can greatly help in
 understanding the problem, whether it is development of some tool or research.

Knowing this we introduced some functions that utilise broadly used visualisation
 software for biomolecules -- PyMOL. Integration is not full, there is no PyMOL plugin
 associated with PyDesc, rather a bunch of useful scripts. Reason for that being our
 understanding of how PyDesc could be used, which is: as scripts, probably run on 
 clusters without graphical environment. Keeping that in mind we imagine that reason
 why one would need to run PyDesc with PyMOL is either to take a look at structures or 
 contact maps in the middle of some kind of script or at the very end of research, to 
 visualise results. Either way -- not using PyMOL to perform any heavy computation.

Basic functionality covered by provided functions covers:
* structure visualisation
* contact map visualisation
* creating selection representing different PyDesc objects
* (future) trajectory visualisation
* (future) visualisation of PyDesc alignments as superposition:
    * of structures
    * of descriptors

### Integrating PyMOL with PyDesc

PyMOL uses its ovn copy od Python, so just having PyDesc installed in Python is not
 enough.
There are two ways of integrating PyMOL with PyDesc:
1. install PyDesc in PyMOL's copy of Python
1. build PyMOL from the source and include installed PyDesc library

#### Install PyDesc

TBD

#### Build own PyMOL from source

1. clone repo `git clone https://github.com/schrodinger/pymol-open-source/blob/master`
1. `cd pymol-open-source`
1. install all pre-requirements
1. run `PREFIX_PATH=<site-packages> python setup.py install --prefix=<install-dir>`,
 where `<site-packages>` is path to `site-packages` dir with PyDesc installed of your
  local Python 3.6+ and `<install-dir>` is directory to install pymol to.
1. to run pymol execute `<install-dir>/bin/pymol`

Trick for `freetype`, which sometimes is not recognised by C++ compiler:
- locate you freetype directory (`<freetype-dir>`)
- set variable `CPATH=<freetype-dir>` in the same line where `python setup.py` is called

### Configuration

There is no related configuration.

### API

All functions related to PyMOL are provided in `api` module.
This module requires PyMOL to be installed, otherwise it is impossible (well, very
 hard) to import it.

It provides subsequent functions:
* draw_structures
* draw_contact
* draw_contact_maps
* select

#### Visualising structures

 ```python
from pydesc.api.pymol import draw_structures
from pydesc.structure import StructureLoader

loader = StructureLoader()
structures = loader.load_structures("2ljp") # 20 NMR structure

draw_structures(structures)
 ```

For more detailed description of structure loading look [here](#loading-structures).

Line `draw_structures(structures)` creates single object in PyMOL and puts 20
 structures loaded by PyDesc as separate states of that object.
Every structure that we want to visualise in PyMOL is added to registry with its PyMOL 
 object name and state. Changing any of it in PyMOL will have unpredicted consequences.
 
`draw_structures(structures, split_states=True)` would load them as 20 separate objects.

Note that is order to load single structure one still needs to pass list (with single 
 item) as first argument.

#### Visualising trajectories

To draw trajectories use different function:
 ```python
from pydesc.api.pymol import draw_trajectory
from pydesc.api.structure import get_structures
from pydesc.api.trajectory import from_frames

structures = get_structures("2ljp")
trajectory = from_frames(structures)

draw_trajectory(trajectory)
 ```
In this case, if visualization was the only point, it would make no sense to
 turn sequence of structures into a trajectory just to visualise it, but a) there
 are other reasons to keep structures as trajectories (e.g. it is more compact)
 and b) some data is already stored as trajectory.

It is more costly to work with trajectories and their visualization, so it 
 is recommended to load small batches of frames at once or use pymol in 
 command-line mode (with `-c` option, which runs scripts for batch processing)
 if possible.

#### Visualising contact maps

To see how to draw contact maps consider following code:
 ```python
from pydesc.api.criteria import get_default_protein_criterion
from pydesc.api.pymol import draw_contact_map
from pydesc.api.pymol import draw_structures
from pydesc.contacts import ContactMapCalculator
from pydesc.structure import StructureLoader

loader = StructureLoader()
structures = loader.load_structures("2ljp") # 20 NMR structures

draw_structures(structures)

criterion = get_default_protein_criterion()
for structure in structures:
    cm_calculator = ContactMapCalculator(structure, criterion)
    contact_map = cm_calculator.calculate_contact_map()
    draw_contact_map(contact_map, structure)
 ```

Note that `draw_structures(structures)` line is required in order to have
 contact maps drawn. To learn more about this function read
 [this section](#visualising-structures).

For description of contact map calculation look [here](#contact-maps).

Line `draw_contact_maps(contact_maps)` will add two contact maps: red one showing all 
 certain (according to given contact criterion) contacts (contact value 2); and orange 
 one for uncertain contacts (contact value 1).
 Since maps were calculated for structures stored in different states of single PyMOL 
 object -- to see then one needs to change visualised state (in PyMOL press ">" and "<"
 in right bottom panel to change states).
If structures were loaded with `split_states` argument set to `True
`, `draw_contact_maps` would simply draw two separate maps for each structure.

Note that to draw single contact map one still needs to pass single-item list as first
 argument.

Function `draw_contact_maps` also takes optional arguments:
* `split_contacts` -- by default set to False, otherwise it creates a different object 
 for each contact in contact map rather than storing all contacts as two contact maps.
* `structures` -- list of structures to draw contact maps on. By default, each map is 
 drawn on structure it was calculated for. These arguments allow drawing contact map 
 calculated for one conformation on another conformation. It only makes sense to draw 
 contact maps calculated for single conformation on another conformation of the same 
 molecule, so e.g. different NMR frames, different x-ray crystals of the same protein 
 or different frames of molecular dynamics. See example below.
* `point` -- name of atom or pseudoatom to be pointed at by two ends of dashed line 
 representing contact. By default, it is `gc`, which is geometrical center of set of
 atoms (or its side chain). Function will skip contacts for mers that lack chosen
 (pseudo)atom (`gc` is always present). That option makes sense only for contact
 maps for very specific cases, e.g. for nucleic acids it makes sense to chang it to
 (`ring_center`).

Let's take a closer look at `2ljp` NMR structure.
 ```python
from pydesc.api.criteria import get_default_protein_criterion
from pydesc.api.pymol import draw_contact_map
from pydesc.api.pymol import draw_structures
from pydesc.contacts import ContactMapCalculator
from pydesc.structure import StructureLoader

loader = StructureLoader()
structures = loader.load_structures("2ljp") 

draw_structures(structures, split_states=True)

criterion = get_default_protein_criterion()
cm_calculator = ContactMapCalculator(structures[0], criterion)
contact_map = cm_calculator.calculate_contact_map()

draw_contact_map(contact_map, structures=[structures[1]])
 ```
In that case we see how contacts calculated for frame 1 would "stretch" while mobile 
 N-terminus changes conformation, floating away from middle of protein.

Additionally, it is possible to draw single contacts with `draw_contact` function. It 
 takes three positional arguments:
 * structure or substructure
 * two atom set indicators
For example, to previous example one could add:
 ```
draw_contact(structure[0], 12, 34)
segment = structure[0][12:34]
draw_contact(segment, 13, 33)
 ```
to add single contacts, even though they are not present in contact maps.

Similarly to `draw_contact_maps`, `draw_contact` takes optional arguments:
 * `point` -- doing the same for both functions.
 * `contact_name` -- to change default name to the one passed as this argument.
 * `gap` -- float, length of gap between dashes of drawn line.

#### Visualising frequency maps

In this example we turn NMR structure into a trajectory, calculate frequency 
 map and then visualize it:
 ```python
from pydesc.api.criteria import get_default_protein_criterion
from pydesc.api.pymol import draw_frequency_map
from pydesc.api.pymol import draw_trajectory
from pydesc.api.structure import get_structures
from pydesc.api.trajectory import from_frames
from pydesc.contacts.maps import FrequencyMapCalculator

structures = get_structures("2ljp")
trajectory = from_frames(structures)
criterion = get_default_protein_criterion()
fm_calculator = FrequencyMapCalculator(trajectory, criterion)
fmap = fm_calculator.calculate_frequency_map()

draw_trajectory(trajectory)
draw_frequency_map(fmap, trajectory)
```

Method `draw_frequency_map` draws contacts to every frame, so it might take a while.
 Contacts are always split. Frequencies are represented both as gaps in dashed lines
 (solid lines for frequency 1.0) and color (red for greater frequency, through green
 to blue for lowest frequencies).

Similarly to `draw_contact_maps`, `draw_frequency_map` allows setting `point` argument.
 Additionally, it takes `frames` which is expected to be a list of frames over which 
 contacts are to be sketched. As adding so many contacts to pymol is costly -- this 
 provides opportunity to save some time, especially if there is only one frame to be
 turned into a picture.
 In that case remember that `frames` is still expected to be a list of one element,
 e.g. `frames=[0]`.

#### Visualising pseudoatoms

It is possible to calculate some pseudoatoms for mers and entities.
 They are often important for contact maps and descriptors.
 To visualise them one could use `draw_pseudoatoms` function:
 ```python
from pydesc.api.pymol import draw_pseudoatoms
from pydesc.api.pymol import draw_structures
from pydesc.structure import StructureLoader

loader = StructureLoader()
structures = loader.load_structures("2ljp") # 20 NMR structure

draw_structures(structures)
for structure in structures:
    draw_pseudoatoms(structure, "gc")
 ```
Executing this example should display additional pseudoatoms near each residue.
 They will all be accessible in PyMOL as `2ljp_gc`, unless additional argument 
 `split_states=True` is passed. With that argument set to `True`, each pseudoatom 
 has its own name.
 It is also possible to display pseudoatoms as "bonded" with some other (pseudo)atoms 
 of their set of atoms. To do so simply add `anchor_name=<(psuedo)atom name>`, e.g.:
 `draw_pseudoatoms(structure, 'cbx', anchor_name='CA')`
 Putting that line inside last loop should display `cbx` pseudoatoms "bonded" with 
 alpha carbons.
 Word "bonded" is abused here, as those pseudoatoms are in fact bonded with another
 pseudoatoms created at exact position of actual anchor.

#### Selecting PyDesc substructures in PyMOL

Term "selection" means something different in PyDesc and in PyMOL.
Different types of substructures, e.g. descriptors, that one can create in PyDesc can 
 be turned into selections in PyMOL, if their parent structure was already drawn.
 To do so, one needs function `select`, which takes two arguments:
 * (sub)structure (PyDesc object), e.g. segment or descriptor.
 * name of selection to be created (its optional, by default name "sele" is passed).

Selections can be later visualised in different ways.

Note that PyMOL selections have no sense of object states, which might cause some 
 unexpected behaviour for structures loaded as different states.

## Geometry

TBD

## END
