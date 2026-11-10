# Chemistry plugin

Chemistry is installed by default with TopoPPI. Its `topoppi_chemistry` code
package and entry points are bundled in the `topoppi` distribution, together
with dependencies on ProLIF, MDAnalysis and RDKit. Chemistry is installed and
removed with Core. Affinity and Hotspot remain optional.
Existing `topoppi[chemistry]` installation commands still work.

ProLIF is the supplied chemical interaction detection backend. Backend discovery
uses metadata; ProLIF, MDAnalysis and RDKit load when detection runs.

In the desktop app, open **Plugins → Manage Plugins...**, select **Chemical
interactions** and choose **Open**. The next map's backend setting becomes
visible in Advanced. The default selection is ProLIF.

## Use selected Core atoms

```python
import topoppi

structure = topoppi.load_structure("complex.cif")
interface = topoppi.get_interface(structure, "A", "B")
result = topoppi.compute_interactions(interface, backend="prolif")
# Equivalent: topoppi.run("chemistry.detect", interface, backend="prolif")
print(result.engine, result.interactions)
```

Both `Interface` and `GeometryRepresentation` supply the selected atoms and
coordinates. Detection reuses these objects. ProLIF operates on separate,
hydrogenated RDKit copies while preserving shared atoms and surface geometry.
A restored native geometry remains usable after source removal. The interaction format and ProLIF adapter require **one author
chain per partner**. Multi-character mmCIF chain identifiers and uniquely resolved PDB insertion codes are supported.

```bash
topoppi plugins info chemistry
topoppi chemistry backends
topoppi chemistry detect complex.cif --chain-a A --chain-b B --output interactions.json
topoppi complex.cif -A A -B B --interactions interactions.json -o interface.png
```

Direct detection writes canonical records and object provenance. The mapping
pipeline's generated sidecars also carry the schema-4 source checksum, chain
direction and software version. Formal benchmark runs require these additional
provenance fields.

The mapping default remains `--interaction-source prolif`. `--interactions`
accepts an existing interaction document; the older `--prolif` spelling is an
alias. The desktop selector discovers installed backend names. Selecting
`geometric` uses Core distance-defined contact partners. Chemistry also provides
a distance-based chemical-label heuristic for diagnostic display, separate from
backend interaction detection.

## Records and rendering

`topoppi.InteractionResult` contains `engine`, `chain_a`, `chain_b`,
`interactions`, and optional `metadata`. Each interaction has `res_a_seq`,
`res_b_seq` and `interaction`; residue tokens preserve insertion codes. Multiple
records for a residue pair retain their multiplicity; partner degree counts
distinct partner residues. The adapter preserves the caller's partner direction.

Canonical types are `HydrogenBond`, `Ionic`, `PiStacking`, `PiCation`,
`Hydrophobic`, `HalogenBond`, `MetalCoordination`, `VdWContact`, `PolarContact`
and `Other`. Core reads records and renders annotations without calling the
backend. Legacy ProLIF names are decoded as file-format aliases. Saved atlases
retain the recorded backend and annotations, including the old ProLIF atlas
keys for compatibility.

## Adding a backend

A backend package registers a callable in its own `pyproject.toml`:

```toml
[project.entry-points."topoppi.chemistry.backends"]
example = "example_chemistry:detect"
```

```python
from topoppi import InteractionResult

def detect(interface, *, context):
    context.check_cancelled()
    first, second = interface.partners
    # Perceive interactions from copies of the selected atoms and coordinates.
    # Return records using author residue identities and canonical type names.
    records = perceive(first, second)
    return InteractionResult(
        engine="example", chain_a=first.chains[0], chain_b=second.chains[0],
        interactions=records,
    )
```

`perceive` above represents the backend's own implementation. It must support
the documented single-chain format and report unsupported partner groups.
Operate on copies of the supplied atoms and arrays, preserving native surfaces
and independence from the original structure file. Use the supplied
`ExecutionContext` for logging and cancellation.
The returned engine must match the registered name. Give each backend a unique
name; duplicate providers for one name produce an error.

After installation, restart Python (or refresh plugin discovery) and select the
registered name through `backend=`, `--backend`, `--interaction-source` or the
desktop selector. Core discovers the backend through its registered entry point.
Formal TopoPPI benchmarks use the prescribed ProLIF records and provenance
requirements.
