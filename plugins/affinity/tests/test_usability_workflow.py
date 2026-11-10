"""User-facing execution contracts: validate early, resume, cancel and stream scores."""
import hashlib
import json
import threading
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pytest
from sklearn.dummy import DummyRegressor

from topoppi_affinity import cli, sequence_embeddings, workflow
from topoppi_affinity.features import parse_mutations
from topoppi_affinity.models import AffinityRegressor

FIXTURE = Path(__file__).parent / 'fixtures/tiny_complex.pdb'


@pytest.fixture
def model(tmp_path):
    bundle = AffinityRegressor(DummyRegressor(strategy='constant', constant=1.5).fit([[1]], [1.5]),
        ('contact.residue_pairs',), 'ddg', 'test constant', {'scope': 'software test only'})
    path = tmp_path / 'model.joblib'
    joblib.dump(bundle, path)
    return path


@pytest.mark.parametrize('changes,match', [({'temperature': -1}, 'Temperature'),
    ({'wt_kd': 0}, 'WT Kd'), ({'mutations': ['WA1F']}, 'WT identity'),
    ({'mutations': []}, 'substitution'), ({'chains_b': 'A'}, 'overlap')])
def test_invalid_requests_fail_before_geometry_or_output(tmp_path, model, monkeypatch, changes, match):
    def unexpected(*args, **kwargs):
        raise AssertionError('Expensive geometry must not run for invalid input')
    monkeypatch.setattr(workflow, 'prepare_interface', unexpected)
    args = dict(structure=FIXTURE, chains_a='A', chains_b='B', model=model,
                output=tmp_path / 'analysis', mutations=['GA1W'])
    args.update(changes)
    with pytest.raises(ValueError, match=match):
        workflow.analyze(**args)
    assert not (tmp_path / 'analysis').exists()


@pytest.mark.parametrize('existing', ['predictions.csv', 'requests.csv', 'interface.npz', 'features/descriptors.npz'])
def test_existing_result_files_are_preserved_before_an_analysis_report_exists(tmp_path, model, monkeypatch, existing):
    output = tmp_path / 'analysis'
    original = output / existing
    original.parent.mkdir(parents=True)
    original.write_bytes(b'existing result to preserve')
    monkeypatch.setattr(workflow, 'prepare_interface', lambda *a, **k: pytest.fail('Output conflict reached geometry'))
    args = dict(structure=FIXTURE, chains_a='A', chains_b='B', model=model,
                output=output, mutations=['GA1W'])
    with pytest.raises(ValueError, match='Output already contains result files'):
        workflow.analyze(**args)
    assert original.read_bytes() == b'existing result to preserve'
    assert not (output / 'analysis.json').exists()
    retry = workflow.validate_inputs(**args, retry=True)
    assert retry['output'] == output.with_name('analysis-retry')


def test_output_directory_can_retain_unrelated_notes(tmp_path, model):
    output = tmp_path / 'analysis'
    output.mkdir()
    notes = output / 'notes.txt'
    notes.write_text('Experiment notes')
    result = workflow.analyze(FIXTURE, 'A', 'B', model, output, mutations=['GA1W'])
    assert result['prediction_status'] == 'ok'
    assert notes.read_text() == 'Experiment notes'


@pytest.mark.parametrize('artifact', ['interface', 'features'])
def test_cli_reports_truncated_npz_without_traceback(tmp_path, model, capsys, artifact):
    output = tmp_path / 'analysis'
    workflow.analyze(FIXTURE, 'A', 'B', model, output, mutations=['GA1W'])
    path = output / ('interface.npz' if artifact == 'interface' else 'features/descriptors.npz')
    content = path.read_bytes()
    path.write_bytes(content[:len(content) // 2])
    argv = ['inspect', str(output), '--json'] if artifact == 'interface' else [
        'predict', '--model', str(model), '--features', str(path), '--output', str(tmp_path / 'scores.csv')]
    with pytest.raises(SystemExit) as stopped:
        cli.main(argv)
    assert stopped.value.code == 2
    captured = capsys.readouterr()
    assert 'Invalid or incomplete NPZ archive' in captured.err
    assert 'Traceback' not in captured.err
    assert not (tmp_path / 'scores.csv').exists()


@pytest.mark.parametrize('table,message', [
    ('experiment,value\ncontrol,7.2\n', 'Saved predictions table is missing columns'),
    ('variant_id,mutations,predicted,status,error,sites\n', 'Saved predictions table has no rows'),
])
def test_saved_prediction_table_must_contain_result_columns_and_rows(tmp_path, model, table, message):
    output = tmp_path / 'analysis'
    workflow.analyze(FIXTURE, 'A', 'B', model, output, mutations=['GA1W'])
    (output / 'predictions.csv').write_text(table)
    with pytest.raises(ValueError, match=message):
        workflow.load_analysis(output)


def test_cancelled_result_can_reopen_a_header_only_prediction_table(tmp_path, model):
    output = tmp_path / 'analysis'
    report = workflow.analyze(FIXTURE, 'A', 'B', model, output, mutations=['GA1W'])
    report['prediction_status'] = 'cancelled'
    workflow._write_report(report, output)
    (output / 'predictions.csv').write_text('variant_id,mutations,predicted,status,error,sites\n')
    document = workflow.load_analysis(output)
    assert document.predictions.empty
    assert document.report['prediction_status'] == 'cancelled'


def test_cli_render_rejects_failed_variant_with_row_error(tmp_path, model, capsys):
    output = tmp_path / 'analysis'
    workflow.analyze(FIXTURE, 'A', 'B', model, output, mutations=['GA1W', 'WA1F'])
    image = tmp_path / 'failed.png'
    with pytest.raises(SystemExit) as stopped:
        cli.main(['render', str(output), '--variant', 'variant_002', '--output', str(image)])
    assert stopped.value.code == 2
    assert 'WT identity mismatch' in capsys.readouterr().err
    assert not image.exists()


def test_cli_render_explains_variant_on_other_partner(tmp_path, model, capsys):
    output = tmp_path / 'analysis'
    workflow.analyze(FIXTURE, 'A', 'B', model, output, mutations=['GA1W'])
    image = tmp_path / 'other.png'
    assert cli.main(['render', str(output), '--variant', 'variant_001', '--side', '1', '--output', str(image)]) == 0
    assert 'use --side 0 to view its footprint' in capsys.readouterr().err
    assert image.stat().st_size > 1000


def test_numeric_chains_structured_csv_and_leading_zero_ids(tmp_path, model):
    pdb = tmp_path / 'numbered.pdb'
    lines = FIXTURE.read_text().splitlines(keepends=True)
    pdb.write_text(''.join(line[:21] + {'A': '1', 'B': '2'}[line[21]] + line[22:]
        if line.startswith(('ATOM  ', 'HETATM')) and line[21] in {'A', 'B'} else line for line in lines))
    variants = pd.DataFrame([{'variant_id': '001', 'chain': '1', 'position': '1', 'wt': 'G', 'mt': 'W'},
                             {'variant_id': '001', 'chain': '2', 'position': '1', 'wt': 'G', 'mt': 'A'}])
    report = workflow.analyze(pdb, '1', '2', model, tmp_path / 'analysis', variants=variants)
    doc = workflow.load_analysis(report['output_directory'])
    assert doc.predictions.variant_id.tolist() == ['001']
    assert set(doc.residues.chain) == {'1', '2'}
    assert doc.selection('001', 0)[1] == ['1:GLY:1']
    assert doc.selection('001', 1)[1] == ['2:GLY:1']
    assert parse_mutations('AA:G-3AV')[0]['insertion'] == 'A'
    with pytest.raises(ValueError, match='Unknown variant ID'):
        doc.selection('typo')


def test_prediction_event_precedes_layout_and_retry_reuses_geometry(tmp_path, model, monkeypatch):
    events = []
    def receive(event):
        if event['stage'] == 'prediction_ready':
            doc = workflow.load_analysis(event['analysis_json'])
            assert float(doc.predictions.predicted.iloc[0]) == 1.5
        events.append(event)
    def fail_layout(*args, **kwargs):
        assert any(e['stage'] == 'prediction_ready' for e in events)
        raise RuntimeError('solver unavailable')
    monkeypatch.setattr('topoppi.pipeline.layout_interface', fail_layout)
    path = tmp_path / 'analysis'
    original = workflow.analyze(FIXTURE, 'A', 'B', model, path, mutations=['GA1W'], layout=True, event_callback=receive)
    assert original['prediction_status'] == 'ok' and original['layout_status'] == 'error'
    def unexpected(*args, **kwargs):
        raise AssertionError('Unchanged native geometry should be reused')
    monkeypatch.setattr(workflow, 'prepare_interface', unexpected)
    second = workflow.analyze(FIXTURE, 'A', 'B', model, path, mutations=['GA1F'], retry=True)
    assert second['output_directory'] != str(path)
    assert second['reused_interface'] == str(path / 'interface.npz')
    assert second['analysis_id'] != original['analysis_id']
    before = (path / 'predictions.csv').read_bytes()
    workflow.layout_existing(path)
    assert (path / 'predictions.csv').read_bytes() == before


def test_cancel_during_preparation_is_not_an_error(tmp_path, model, monkeypatch):
    def cancelled(*args, **kwargs):
        raise InterruptedError('Interface preparation cancelled')
    monkeypatch.setattr(workflow, 'prepare_interface', cancelled)
    report = workflow.analyze(FIXTURE, 'A', 'B', model, tmp_path / 'analysis', mutations=['GA1W'])
    assert report['prediction_status'] == report['interface_status'] == 'cancelled'
    assert workflow.load_analysis(report['output_directory']).interface is None


def test_shared_embedding_cache_skips_model_load_and_cancel_is_checked(tmp_path, monkeypatch):
    sequences = tmp_path / 'sequences'
    sequences.mkdir()
    pd.DataFrame([{'sequence_id': 'sequence_000007', 'sequence': 'AG', 'length': 2}]).to_csv(sequences / 'sequences.csv', index=False)
    shared = tmp_path / 'shared'
    shared.mkdir()
    path = shared / (hashlib.sha256(b'AG').hexdigest() + '.npz')
    np.savez_compressed(path, embedding=np.ones((2, 1152)), sequence='AG', model_name='esmc_600m')
    def unexpected(*args, **kwargs):
        raise AssertionError('Cached embeddings must not load/download a model')
    monkeypatch.setattr('topoppi_affinity.resources.esmc_checkpoint', unexpected)
    events = []
    result = sequence_embeddings.extract_embeddings(sequences, shared_cache=shared, progress=events.append)
    assert result['newly_embedded'] == 0
    assert (sequences / 'embeddings/sequence_000007.npz').exists()
    assert events[-1]['current'] == events[-1]['total'] == 1
    cancel = threading.Event()
    cancel.set()
    with pytest.raises(InterruptedError):
        sequence_embeddings.extract_embeddings(sequences, shared_cache=shared, cancel_event=cancel)


@pytest.mark.parametrize('output_flag', ['--output', '--output-dir'])
def test_json_cli_remains_clean_and_inspect_contains_scores(tmp_path, model, monkeypatch, capsys, output_flag):
    from topoppi_affinity import inference
    real = inference.predict_features
    def chatty(*args, **kwargs):
        print('Library progress which must go to stderr')
        return real(*args, **kwargs)
    monkeypatch.setattr(inference, 'predict_features', chatty)
    path = tmp_path / 'analysis'
    args = ['analyze', '--structure', str(FIXTURE), '--chains-a', 'A', '--chains-b', 'B',
            '--model', str(model), '--mutations', 'GA1W', output_flag, str(path), '--json']
    assert cli.main(args) == 0
    output = capsys.readouterr()
    result = json.loads(output.out)
    assert result['predictions'][0]['predicted'] == 1.5
    assert 'Library progress' in output.err
    assert cli.main(['inspect', str(path), '--json']) == 0
    assert json.loads(capsys.readouterr().out)['predictions'][0]['target'] == 'ddg'
    assert cli.main(['inspect', str(path)]) == 0
    assert 'kcal/mol' in capsys.readouterr().out


def test_example_export_then_run_can_reuse_original_inputs(tmp_path, monkeypatch, capsys):
    path = tmp_path / 'example'
    assert cli.main(['example', '--output', str(path)]) == 0
    capsys.readouterr()
    called = {}
    def run(*args, **kwargs):
        called.update(kwargs)
        raise RuntimeError('test stops before actual prediction')
    monkeypatch.setattr(workflow, 'analyze', run)
    with pytest.raises(SystemExit):
        cli.main(['example', '--output', str(path), '--run'])
    assert called['retry'] is True and called['task'] == 'absolute'
    assert 'already exist' not in capsys.readouterr().err
