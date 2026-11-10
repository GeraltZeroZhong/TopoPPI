# Interface representations and physical descriptors

`topoppi_representation` provides native-mesh, UV and direct-table encodings,
surface-property queries, and a named physicochemical descriptor dictionary.
It ships with TopoPPI 3.0 alongside `topoppi_chemistry`.

```bash
python -m pip install 'topoppi[representation]'
```

The representation extra installs compression and descriptor dependencies.
Importing TopoPPI or `topoppi_representation` leaves these dependencies unloaded.
Optional analysis scripts use the additional dependencies listed in
[`scripts/analysis/README.md`](../scripts/analysis/README.md).

## Compute descriptors from a structure

```bash
topoppi representation descriptors \
  --structure complex.pdb --chains-a A,B --chains-b C \
  --output descriptors --residues --contacts
```

The command accepts PDB and mmCIF structures and writes `descriptors.csv`, with
one row and 477 named columns, plus `metadata.json` with the structure, partner
groups and calculation settings. `--residues` adds `residues.csv` with
per-residue SASA, chemical classes, interface regions and contact counts;
`--contacts` adds `contacts.csv` with residue pairs within 5.5 Å. Residue-table
side labels follow sorted partner chain groups, recorded in the metadata.
Choose a new output directory for each export to keep its descriptor and optional
tables together.

`--output-dir` is an alias for `--output`. Add `--json` for descriptor values and
saved-file locations on standard output; progress messages go to standard error.
The SASA defaults are 960 points per atom and a 1.4 Å probe radius. Set
`--n-points` and `--probe-radius` to change them. The same command is available
as `python -m topoppi_representation descriptors`.

From Python:

```python
import topoppi
from topoppi_representation import physical_descriptors, residue_physics

structure = topoppi.load_structure("complex.pdb")
interface = topoppi.get_interface(structure, "A,B", "C")
values = physical_descriptors(interface, n_points=960)
residues, contacts = residue_physics(interface, n_points=960)
```

`physical_descriptors` returns the 477 named fields of
`physicochemical_477_v1`. `residue_physics` returns tables of all residues in
both partners and their cross-partner contacts. Both accept a Core `Interface`
or `GeometryRepresentation` and reuse its selected atoms. Native mesh generation
is optional for these atom-based calculations.

The dictionary includes solvent-accessible surface area (SASA) for isolated
partners and the bound complex, buried surface area, relative solvent
accessibility (RSA), interface regions, residue composition and contacts.
SASA uses Biopython's Shrake–Rupley implementation with a default probe radius
of 1.4 Å.
RSA uses the Tien et al. theoretical maxima in Biopython's Wilke table. The
contact definition is any selected cross-partner heavy-atom pair at most 5.5 Å
apart. Buried area sums both partners' accessibility losses. Chain groups are
sorted; sum and absolute-difference summaries are partner-symmetric.

The dictionary's polar, apolar, charged and unknown groups are defined in
`descriptors.physical.CLASSES`. Affinity and Hotspot models use their
model-specific chemical classifications.

## Spatial descriptors for substitution queries

Prepare one decoded residue-aware UV surface and reuse it across substitutions:

```python
from topoppi_representation import prepare_spatial_descriptors

geometry = topoppi.compute_geometry(interface)
cache, encoded, report = prepare_spatial_descriptors(geometry)
global_values, global_names, _ = cache.describe()
values, names, coverage = cache.describe([("G", "A", 42, "", "A")])
```

Each substitution tuple contains the wild-type amino acid, author chain ID,
residue number, insertion code and replacement amino acid. Match these values
to your structure. `describe()` returns 60 global descriptors; a substitution
query appends 61 fields, for 121 in total. Local queries use Gaussian kernels
with length scales of 4 and 12 Å along shortest graph paths, truncated at three
times the length scale. `coverage` records uncovered sites. Local summaries use
covered sites and return NaN when all queried sites lack coverage.

Preparation requires residue-aware OptCuts. Pass `optcuts_bin=` to select an
executable, or use the standard Core installation. The default complete encoding
budget is one quarter of the uncompressed native source's byte size, including
the identities and SASA values of all partner residues. `encoded.blob` contains
the selected representation; `report` records its size and settings.
`SpatialDescriptorCache(decode_uv(blob))`
can query a saved encoding. The Hotspot plugin uses this same preparation and
query API for its bundled mutation-derived model.

## Reuse native geometry and UV layouts

```python
from topoppi_representation import from_geometry, save_source, native_surface

geometry = topoppi.compute_geometry(interface)
source = from_geometry(geometry)
save_source(source, "interface_source.npz")
surface = native_surface(source)
```

For UV coding, pass existing layouts explicitly:

```python
source = from_geometry(
    geometry,
    layouts={(0, 0): {"ordinary": first_patch_corner_uv}},
)
```

Each key is `(partner_side, patch_index)`. Each array has shape
`(n_faces, 3, 2)` in that patch's face-corner order. Supply the selected layout
for every patch used by the encoder. Core atoms, geometry and UV coordinates
remain unchanged. The adapter stores native coordinates as float32 and records
the maximum rounding error.

Uniformly weighted residue-aware layouts can be derived from a complete source with
`topoppi_representation.layouts.augment_uniform_layout(source,
optcuts_bin=...)`. This uses native geometry, unit residue weights, a residue
fragmentation weight of 20 and a distortion bound of 4.1. Other stored layouts
are preserved. The returned metadata reports the result for each patch and
the inputs used to calculate the layout.

`ComplexSource` contains residue identities, amino-acid identities and native
patches. Sample and family identifiers are optional metadata fields.
`Surface` carries XYZ coordinates (Å), physical adjacency, partner and patch
membership, and a sparse cell-by-residue area matrix (Å²). Additional residue
fields must be supplied explicitly. The source NPZ and companion JSON describe
these inputs; their file sizes differ from the uncompressed source byte count.

## Encode and decode

```python
from topoppi_representation import encode_uv, decode_uv

encoded = encode_uv(source, layout="ordinary", target_cells=1024,
                    quant_step=0.001)
decoded = decode_uv(encoded.blob)
print(encoded.nbytes)
```

| Representation | Encode/decode functions | Decoded object |
| --- | --- | --- |
| Native mesh | `encode_native`, `decode_native` | Surface, including decoded triangle faces |
| UV cells | `encode_uv`, `decode_uv` | Surface with cell positions, areas and adjacency |
| Direct physical table | `encode_direct`, `decode_direct` | Named scalar values and explicit missing values |

Native encoding uses optional quadric error metric (QEM) simplification,
area-conserving residue transport and Draco geometry compression. UV encoding
uses area-conserving cells, a Haar wavelet transform, quantization and Zstandard.
Separate slots retain disconnected pieces occupying the same pixel. Layout names
distinguish ordinary (`ordinary`), contact-weighted residue-aware (`aware`) and
uniformly weighted residue-aware (`aware_uniform`) coordinates.

`target_cells` is an approximate occupied-cell budget; disconnected slots and
boundaries can increase the realized count. `Encoded.nbytes` counts the complete
bitstream, including headers, connectivity, identities and carried residue fields.
The decoder accepts those bytes alone. `Encoded.metadata` contains diagnostics.
`canonical_source_bytes(source)` returns the uncompressed native source payload
used for source-relative size comparisons.

## Query physical properties

`descriptors.surface` provides `hydrophobic_organization`,
`proximity_composition`, `physical_targets` and `local_queries`. Hydrophobic
organization measures the largest connected hydrophobic component as a fraction
of each partner's hydrophobic surface area. Here, mass denotes transported
surface area in Å². Chemical proximity uses an area-product Gaussian kernel with
a 2 Å length scale and a 6 Å cutoff. Local hydropathy uses physical graph distance
from a residue's surface support, with caller-selected length scales.

`reference_quadrature` and `source.surface_quadrature` integrate native triangle
dual regions while preserving carried residue area. Query results include
applicability and surface-coverage information. Transported source area and area
recomputed from quantized coordinates are separate quantities.

`descriptors.physical.recover_from_sasa(table)` reconstructs the 430 dictionary
fields determined by partner side, amino-acid identity and the three SASA columns.
`recover_from_carrier(surface)` uses the same calculation when those columns are
explicitly carried in the decoded representation. The remaining 47 dictionary
fields require atomic contact information.

## Optional readouts

The analysis scripts compare a chosen representation with its source and fit
physical-target recovery models or graph readouts with specified depth and width.
They accept explicit inputs and data partitions. Their model classes live in the optional
`topoppi_representation.analysis` modules, so saved models retain stable import
paths. Basic codecs and descriptors run independently of PyTorch and scikit-learn.

See the [analysis script guide](../scripts/analysis/README.md) for input schemas,
commands, outputs and dependency requirements.
