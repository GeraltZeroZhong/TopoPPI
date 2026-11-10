"""TopoPPI Affinity: prepare resources, predict and review interface-linked results."""
from __future__ import annotations

import argparse
import contextlib
import importlib.util
import json
import signal
import sys
import threading
import time
import zipfile
from pathlib import Path


def _log(message):
    print(message, file=sys.stderr, flush=True)


class _DownloadProgress:
    def __init__(self):
        self.last = {}

    def __call__(self, event):
        resource = event.get('resource', 'Resource')
        now = time.monotonic()
        if event.get('event') != 'downloading' or resource not in self.last or now - self.last[resource] >= 1:
            self.last[resource] = now
            current, total = event.get('current', 0), event.get('total')
            suffix = f' ({current / total:.0%})' if total else ''
            _log(event.get('message', resource) + suffix)


def analysis_summary(document):
    report = document.report
    columns = [name for name in ('variant_id', 'mutations', 'predicted', 'target', 'status',
        'predicted_kd_M', 'predicted_ddg_kcal_mol', 'anchored_mutant_kd_M', 'components_used', 'error')
        if name in document.predictions]
    return {'analysis': str(document.path), 'prediction_status': report['prediction_status'],
            'layout_status': report['layout_status'], 'layout_attempt_status': report.get('layout_attempt_status'),
            'model': report.get('model', {}).get('scope'),
            'device': report.get('device'), 'partners': report.get('interface', {}).get('groups'),
            'warnings': report.get('warnings', []), 'scope': report.get('scope_assignment'),
            'errors': {key: report[key] for key in ('interface_error', 'prediction_error') if key in report},
            'predictions': document.predictions[columns].to_dict('records')}


def _print_summary(summary, as_json=False):
    if as_json:
        print(json.dumps(summary, indent=2, allow_nan=False))
        return
    print(f"Prediction: {summary['prediction_status']} | Layout: {summary['layout_status']}")
    if summary.get('layout_attempt_status'):
        print(f"Latest atlas attempt: {summary['layout_attempt_status']}")
    print(f"Results: {summary['analysis']}")
    if summary.get('model'):
        print(summary['model'])
    if summary.get('scope'):
        print(f"Applicability ({summary['scope'].get('context', 'unspecified')}): {summary['scope'].get('user_basis') or 'User-selected model scope'}")
    for warning in summary.get('warnings', []):
        print(f'Notice: {warning}')
    for message in summary.get('errors', {}).values():
        print(f'Error: {message}')
    print('Variant\tSubstitutions\tPrediction\tUnits\tStatus')
    for row in summary.get('predictions', []):
        value = row.get('predicted', '')
        units = {'pkd': 'pKd', 'ddg': 'kcal/mol'}.get(row.get('target'), 'unavailable')
        print(f"{row['variant_id']}\t{row['mutations'] or 'whole complex'}\t{value}\t{units}\t{row['status']}")
        if row.get('predicted_kd_M') != '' and 'predicted_kd_M' in row:
            print(f"  Kd: {row['predicted_kd_M']} M")
        if row.get('anchored_mutant_kd_M') != '' and 'anchored_mutant_kd_M' in row:
            print(f"  Measured-WT-anchored mutant Kd: {row['anchored_mutant_kd_M']} M")
        if row.get('error'):
            print(f"  {row['error']}")



def _print_resource_status(status, *, tasks=None):
    print(f"Model cache: {status['cache_dir']}")
    for task in tasks or status['models']:
        model = status['models'][task]
        if model['ready']:
            print(f"{model['name']}: ready")
        else:
            print(f"{model['name']}: preparation needed")
            print("  Required: " + ", ".join(model['missing_resources']))
            print(f"  Prepare: topoppi-affinity models download --task {task}")
    if not status['esmc']['installed'] and any(task in {'known', 'unseen'} for task in (tasks or status['models'])):
        print(f"Optional neural models: shared ESMC download {status['esmc']['size_bytes'] / 1e9:.2f} GB.")
    if any(model['ready'] for model in status['models'].values()):
        print("Open predictions: topoppi-affinity gui")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    from topoppi_affinity import __version__

    prog = 'topoppi affinity' if Path(sys.argv[0]).name == 'topoppi affinity' else 'topoppi-affinity'
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument('--version', action='version', version=f'TopoPPI Affinity {__version__}')
    commands = parser.add_subparsers(dest='command', required=True)
    review = commands.add_parser('review', help='Open linked 3D, contact graph and atlas review without model downloads')
    review.add_argument('analysis', type=Path, help='Analysis directory, analysis.json or saved review.json')
    review.add_argument('--view', choices=['structure', 'contacts', 'atlas'], default='atlas')
    analyze = commands.add_parser('analyze', help='Predict from a structure; optionally generate an atlas')
    analyze.add_argument('--structure', required=True, type=Path)
    analyze.add_argument('--chains-a', required=True, help='Complete partner group, e.g. A,B')
    analyze.add_argument('--chains-b', required=True)
    analyze.add_argument('--task', choices=['absolute', 'mutation', 'known', 'unseen'], help='Scientific task; defaults to mutation when substitutions are supplied, otherwise absolute')
    analyze.add_argument('--model', type=Path, help='Advanced: explicit compatible model file')
    analyze.add_argument('--output', '--output-dir', type=Path, help='New directory; defaults to a unique folder under topoppi-results/')
    analyze.add_argument('--retry', action='store_true', help='Create a new sibling of --output and reuse an unchanged saved interface')
    analyze.add_argument('--check-only', action='store_true', help='Validate structure, variants, model and environment without expensive geometry or prediction')
    analyze.add_argument('--json', action='store_true', help='Print a JSON result summary; all progress goes to stderr')
    changes = analyze.add_mutually_exclusive_group()
    changes.add_argument('--mutations', action='append', default=[], help='One variant per flag; comma joins substitutions, e.g. A:G42V,1:A7F')
    changes.add_argument('--variants', type=Path, help='CSV variant_id,mutations or structured variant_id,chain,position,wt,mt[,insertion]')
    analyze.add_argument('--source', default='user', help='Structure source/numbering namespace')
    analyze.add_argument('--complex-id')
    analyze.add_argument('--mutation-context', choices=['known', 'unseen'])
    analyze.add_argument('--scope-note', default='', help='Basis for selecting known/unseen applicability; no automatic family inference')
    analyze.add_argument('--temperature', type=float, default=298.15)
    analyze.add_argument('--wt-kd', type=float, help='Optional measured WT Kd in M, e.g. 1e-8 for 10 nM')
    analyze.add_argument('--sequence-dir', type=Path, help='Advanced: existing structural sequence registry')
    analyze.add_argument('--cache-dir', type=Path, help='Override the shared model and embedding cache for this run')
    analyze.add_argument('--device', default='auto', help='auto (available CUDA or CPU), cpu, cuda or cuda:N')
    priors = analyze.add_mutually_exclusive_group()
    priors.add_argument('--prior-csv', type=Path, help='Advanced: explicit model prior with matching valid variant IDs')
    priors.add_argument('--prior-repository', type=Path, help='Advanced: author source; installed priors are found automatically')
    analyze.add_argument('--prior-checkpoint', type=Path)
    analyze.add_argument('--layout', action='store_true', help='Generate both atlases after saving predictions')
    analyze.add_argument('--optcuts-bin', default='OptCuts_bin')
    analyze.add_argument('--layout-timeout', type=float, default=600.)
    predict = commands.add_parser('predict', help='Advanced: prediction from named feature NPZs')
    predict.add_argument('--model', required=True, type=Path)
    for name in ['features', 'plm-features', 'prior']:
        predict.add_argument('--' + name, type=Path)
    predict.add_argument('--descriptor-features', dest='legacy_features', metavar='PATH', type=Path,
                         help='Conventional interface-descriptor NPZ')
    predict.add_argument('--temperature', type=float, default=298.15)
    predict.add_argument('--wt-kd', type=float)
    predict.add_argument('--output', required=True, type=Path)
    inspect = commands.add_parser('inspect', help='Show scores, units, warnings and branch status from saved results')
    inspect.add_argument('analysis', type=Path)
    inspect.add_argument('--json', action='store_true', help='Machine-readable result summary')
    inspect.add_argument('--full', action='store_true', help='Full provenance JSON, including training metadata')
    render = commands.add_parser('render', help='Export native 3D or saved atlas; no model needed')
    render.add_argument('analysis', type=Path)
    render.add_argument('--output', required=True, type=Path)
    render.add_argument('--side', type=int, choices=[0, 1], default=0, help='0: Partner A, 1: Partner B')
    render.add_argument('--view', choices=['surface', 'atlas'], default='surface')
    render.add_argument('--variant')
    render.add_argument('--descriptor', choices=['contact_degree', 'footprint_area_A2', 'surface_area_A2'])
    layout = commands.add_parser('layout', help='Generate/retry atlases from a saved interface without repeating prediction')
    layout.add_argument('analysis', type=Path)
    layout.add_argument('--optcuts-bin', default='OptCuts_bin')
    layout.add_argument('--layout-timeout', type=float, default=600.)
    layout.add_argument('--json', action='store_true')
    gui = commands.add_parser('gui', help='Open optional analysis/review panel')
    gui.add_argument('analysis', nargs='?', type=Path)
    models = commands.add_parser('models', help='Prepare or inspect cached release models and official upstream resources')
    resource_commands = models.add_subparsers(dest='resource_command', required=True)
    download = resource_commands.add_parser('download', help='Prepare bundled models or download explicitly selected neural resources')
    selection = download.add_mutually_exclusive_group()
    selection.add_argument('--all', action='store_true', help='Prepare both bundled tasks: affinity and mutation effects')
    selection.add_argument('--task', action='append', choices=['absolute', 'mutation', 'known', 'unseen'],
                           help='Task to prepare (default: absolute); repeat for multiple tasks')
    download.add_argument('--cache-dir', type=Path, help='Choose and remember where large models are stored')
    download.add_argument('--no-esmc', action='store_true', help='Prepare estimators/priors only for feature-only inference')
    download.add_argument('--catalog', type=Path, help='Advanced: release catalog supplied with local or mirrored archives')
    status = resource_commands.add_parser('status', help='List cache paths and model readiness without downloading')
    status.add_argument('--cache-dir', type=Path)
    for resource_parser in (download, status):
        resource_parser.add_argument('--json', action='store_true', help='Print complete resource information as JSON')
    doctor = commands.add_parser('doctor', help='Check installed versions, device and resources; no downloads')
    doctor.add_argument('--cache-dir', type=Path)
    example = commands.add_parser('example', help='Export the installed 1BVK example; optionally run absolute affinity')
    example.add_argument('--output', type=Path, required=True)
    example.add_argument('--run', action='store_true')
    example.add_argument('--device', default='auto')
    example.add_argument('--cache-dir', type=Path)
    args = parser.parse_args(argv)
    values = vars(args).copy()
    command = values.pop('command')
    cancel = threading.Event()
    previous_handler = None
    if threading.current_thread() is threading.main_thread() and command not in {'gui', 'review'}:
        previous_handler = signal.getsignal(signal.SIGINT)
        def interrupt(_signum, _frame):
            if cancel.is_set():
                raise KeyboardInterrupt
            cancel.set()
            _log('Cancelling after the current operation; press Ctrl+C again to stop immediately')
        signal.signal(signal.SIGINT, interrupt)
    try:
        if command == 'analyze':
            from topoppi_affinity.workflow import analyze as run
            from topoppi_affinity.workflow import load_analysis, validate_inputs
            as_json, check_only = values.pop('json'), values.pop('check_only')
            with contextlib.redirect_stdout(sys.stderr):
                if check_only:
                    checked = validate_inputs(**values)
                else:
                    report = run(**values, cancel_event=cancel, log=_log)
            if check_only:
                if cancel.is_set():
                    return 130
                print(json.dumps({'status': 'valid', 'output': str(checked['output']), 'device': checked['device'],
                    'model': checked['model_info']['scope'], 'valid_variants': len(checked['requests']) - len(checked['row_errors']),
                    'invalid_variants': {checked['requests'].loc[i, 'variant_id']: error for i, error in checked['row_errors'].items()}}, indent=2))
                return 0 if not checked['row_errors'] else 2
            _print_summary(analysis_summary(load_analysis(report['output_directory'])), as_json)
            if (cancel.is_set() or report['prediction_status'] == 'cancelled' or report['layout_status'] == 'cancelled'
                    or report.get('layout_attempt_status') == 'cancelled'):
                return 130
            return 0 if report['prediction_status'] == 'ok' else 2
        if command == 'predict':
            from topoppi_affinity.inference import predict_features
            output = values.pop('output')
            with contextlib.redirect_stdout(sys.stderr):
                result = predict_features(**values)
            output.parent.mkdir(parents=True, exist_ok=True)
            result.to_csv(output, index=False)
            print(result.to_string(index=False))
            print(f'Results: {output}')
        elif command in {'inspect', 'render', 'layout'}:
            from topoppi_affinity.workflow import layout_existing, load_analysis
            if command == 'layout':
                with contextlib.redirect_stdout(sys.stderr):
                    report = layout_existing(args.analysis, optcuts_bin=args.optcuts_bin, layout_timeout=args.layout_timeout,
                                             cancel_event=cancel, log=_log)
                _print_summary(analysis_summary(load_analysis(args.analysis)), args.json)
                attempt_status = report.get('layout_attempt_status', report['layout_status'])
                if cancel.is_set() or attempt_status == 'cancelled':
                    return 130
                return 0 if attempt_status == 'ok' else 2
            document = load_analysis(args.analysis)
            if command == 'inspect':
                if args.full:
                    print(json.dumps(document.report, indent=2))
                else:
                    _print_summary(analysis_summary(document), args.json)
            else:
                import matplotlib.pyplot as plt
                if args.variant:
                    row, visible, hidden = document.selection(args.variant, args.side)
                    if row.status not in {'ok', 'base_only'}:
                        raise ValueError(f"Variant {args.variant!r} has no completed prediction: {row.error}")
                    if not visible:
                        other_sides = sorted({site['side'] for site in hidden if site['coverage'] == 'footprint'})
                        suggestion = f"; use --side {other_sides[0]} to view its footprint" if other_sides else ''
                        _log(f"Displaying partner {args.side}: selected variant {args.variant!r} has no retained footprint here{suggestion}")
                patches, viz, style = document.view(args.side, args.view, args.variant, args.descriptor)
                if not patches:
                    raise ValueError('No retained interface on this side')
                args.output.parent.mkdir(parents=True, exist_ok=True)
                figure = viz.plot_patches(patches, output_file=str(args.output), show=False, style_config=style)
                plt.close(figure)
                print(args.output)
        elif command == 'models':
            from topoppi_affinity.resources import download_models, resource_status
            selected_tasks = None
            if args.resource_command == 'download':
                selected_tasks = ('absolute', 'mutation') if args.all else tuple(dict.fromkeys(args.task or ['absolute']))
                with contextlib.redirect_stdout(sys.stderr):
                    result = download_models(tasks=selected_tasks, include_esmc=not args.no_esmc,
                        cache_dir=args.cache_dir, progress=_DownloadProgress(), cancel_event=cancel, catalog_path=args.catalog)
            else:
                result = resource_status(cache_dir=args.cache_dir)
            if args.json:
                print(json.dumps(result, indent=2))
            else:
                _print_resource_status(result, tasks=selected_tasks)
        elif command == 'doctor':
            from topoppi import __version__ as core_version

            from topoppi_affinity.resources import resource_status
            from topoppi_affinity.runtime import resolve_device
            dependencies = {name: importlib.util.find_spec(name) is not None for name in ('torch', 'esm', 'sklearn', 'pandas')}
            result = {'core_version': core_version, 'plugin_version': __version__, 'dependencies': dependencies,
                      'resources': resource_status(cache_dir=args.cache_dir)}
            try:
                result['device'] = resolve_device('auto', neural=False)
            except (ValueError, ImportError, RuntimeError) as exc:
                result['device_error'] = str(exc)
            print(json.dumps(result, indent=2))
            return 0 if all(dependencies[name] for name in ('sklearn', 'pandas')) and 'device_error' not in result else 2
        elif command == 'example':
            from topoppi_affinity.resources import export_example
            configuration = export_example(args.output, reuse=args.run)
            if args.run:
                from topoppi_affinity.workflow import analyze as run
                from topoppi_affinity.workflow import load_analysis
                with contextlib.redirect_stdout(sys.stderr):
                    report = run(configuration['structure'], configuration['chains_a'], configuration['chains_b'],
                        task='absolute', output=args.output / 'results', device=args.device, cache_dir=args.cache_dir,
                        source='software_example', cancel_event=cancel, log=_log, retry=True)
                _print_summary(analysis_summary(load_analysis(report['output_directory'])))
                if cancel.is_set() or report['prediction_status'] == 'cancelled' or report['layout_status'] == 'cancelled':
                    return 130
                return 0 if report['prediction_status'] == 'ok' else 2
            print(json.dumps(configuration, indent=2))
        elif command == 'review':
            from topoppi_affinity.atlas_review.gui import main as review_main
            review_main([str(args.analysis), '--view', args.view])
        elif command == 'gui':
            from topoppi_affinity.gui import main as gui_main
            gui_main(args.analysis)
        return 0
    except (InterruptedError, KeyboardInterrupt):
        _log('Cancelled. Completed outputs and cached resources are retained; retry with the same resource destination or analyze --retry.')
        return 130
    except zipfile.BadZipFile as exc:
        parser.exit(2, f'Invalid or incomplete NPZ archive: {exc}. Restore the complete saved file and retry.\n')
    except (OSError, ValueError, KeyError, ImportError, RuntimeError) as exc:
        parser.exit(2, f'{type(exc).__name__}: {exc}\n')
    finally:
        if previous_handler is not None:
            signal.signal(signal.SIGINT, previous_handler)


if __name__ == '__main__':
    raise SystemExit(main())
