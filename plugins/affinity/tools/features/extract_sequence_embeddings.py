#!/usr/bin/env python3
"""Prepare structure-matched sequence mappings and extract frozen ESMC/ESM-2 embeddings."""
import argparse
from pathlib import Path

from topoppi_affinity.sequence_embeddings import extract_embeddings, prepare_sequences

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--manifest', action='append', type=Path)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--checkpoint', type=Path, help='Local checkpoint; required for ESM2')
parser.add_argument('--device', default='cuda')
parser.add_argument('--model', choices=['esmc_600m', 'esm2_650m'], default='esmc_600m')
parser.add_argument('--prepare-only', action='store_true')
parser.add_argument('--limit', type=int)
args = parser.parse_args()
if args.manifest:
    print('Prepared unique sequences:', prepare_sequences(args.manifest, args.output), flush=True)
if not args.prepare_only:
    if args.model == 'esm2_650m' and args.checkpoint is None:
        parser.error('--checkpoint is required for ESM2 embedding')
    print(extract_embeddings(args.output, args.checkpoint, device=args.device, limit=args.limit, model_name=args.model))
