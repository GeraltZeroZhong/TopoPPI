"""Encode one supplied interface and measure physical recovery, size and CPU cost."""

from __future__ import annotations

import argparse
import json
import time
import tracemalloc
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Representation source NPZ with its JSON description")
    parser.add_argument("--codec", required=True, choices=["native", "uv", "direct"])
    parser.add_argument("--layout", choices=["ordinary", "aware", "aware_uniform"], default="ordinary")
    parser.add_argument("--cells", type=int)
    parser.add_argument("--quant-step", type=float, default=0.001)
    parser.add_argument("--fraction", type=float, default=1.0)
    parser.add_argument("--bits", type=int, default=18)
    parser.add_argument("--subdivision", type=int, default=1)
    parser.add_argument("--output", required=True, type=Path, help="New output directory")
    args = parser.parse_args()
    import topoppi_representation as representation
    from topoppi_representation.descriptors.surface import compare_targets, physical_targets
    from topoppi_representation.source import surface_quadrature

    source = representation.load_source(args.source)
    native = representation.native_surface(source)
    reference = physical_targets(native, representation.reference_quadrature(source, args.subdivision))
    tracemalloc.start()
    wall, cpu = time.perf_counter(), time.process_time()
    if args.codec == "native":
        code = representation.encode_native(source, args.fraction, args.bits)
    elif args.codec == "uv":
        code = representation.encode_uv(source, args.layout, args.cells, args.quant_step)
    else:
        values = {f"O_side{side}": reference["O_by_side"][str(side)]["L_H"] for side in (0, 1)}
        values.update({f"X_{name}": value for name, value in zip(reference["X_states"], reference["X"], strict=False)})
        code = representation.encode_direct(values, units={name: "dimensionless" for name in values})
    encoding_wall, encoding_cpu = time.perf_counter() - wall, time.process_time() - cpu
    wall, cpu = time.perf_counter(), time.process_time()
    decoded = getattr(representation, f"decode_{args.codec}")(code.blob)
    decoding_wall, decoding_cpu = time.perf_counter() - wall, time.process_time() - cpu
    wall, cpu = time.perf_counter(), time.process_time()
    if args.codec == "direct":
        comparison = {
            "stored_queries": decoded,
            "maximum_stored_query_error": max(abs(decoded[name] - value) for name, value in values.items()),
        }
    else:
        proximity = surface_quadrature(decoded, subdivision=args.subdivision) if decoded.faces is not None else decoded
        comparison = compare_targets(reference, physical_targets(decoded, proximity))
    query_wall, query_cpu = time.perf_counter() - wall, time.process_time() - cpu
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    result = {
        "codec": args.codec,
        "parameters": {k: v for k, v in vars(args).items() if k not in {"source", "output"}},
        "encoded_bytes": code.nbytes,
        "canonical_source_bytes": len(representation.canonical_source_bytes(source)),
        "encode_wall_seconds": encoding_wall,
        "encode_cpu_seconds": encoding_cpu,
        "decode_wall_seconds": decoding_wall,
        "decode_cpu_seconds": decoding_cpu,
        "query_wall_seconds": query_wall,
        "query_cpu_seconds": query_cpu,
        "python_traced_peak_bytes": peak,
        "comparison": comparison,
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "encoded.bin").write_bytes(code.blob)
    (args.output / "comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
