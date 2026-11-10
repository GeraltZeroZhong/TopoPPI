"""Command-line entry point for benchmark and sensitivity studies."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import fields
from pathlib import Path
from typing import Mapping, Optional, Sequence

from topoppi import __version__
from topoppi.benchmarking.evidence_bundle import (
    read_json_object,
    validate_benchmark_evidence_bundle,
)
from topoppi.benchmarking.runner import BenchmarkRunner
from topoppi.benchmarking.sensitivity import SensitivityBenchmarkRunner, write_sensitivity_plan
from topoppi.config import (
    BenchmarkConfig,
    OptCutsConfig,
    ParameterizationConfig,
    SurfaceConfig,
    TopologyConfig,
    benchmark_config_from_dict,
)
from topoppi.errors import TopoPPIError
from topoppi.json_utils import dump_json, dump_json_atomic


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="topoppi-benchmark",
        description="Check inputs, run benchmarks, and compare TopoPPI parameter settings.",
        epilog="Start with: topoppi-benchmark preflight benchmark.json",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    preflight = subparsers.add_parser(
        "preflight",
        help="Check benchmark inputs, configuration, and planned resources",
        description="Check benchmark inputs, configuration, and planned resources before a run.",
    )
    preflight.add_argument("config", help="Benchmark JSON configuration")
    preflight_output = preflight.add_mutually_exclusive_group()
    preflight_output.add_argument("--json", dest="json_output", action="store_true", help="Print full JSON")
    preflight_output.add_argument("--output-json", help="Write the full preflight report to this path")

    run = subparsers.add_parser("run", help="Run a benchmark configuration",
                               description="Run the configured structures and write benchmark reports.")
    run.add_argument("config", help="Benchmark JSON configuration")
    run.add_argument("--json", dest="json_output", action="store_true", help="Print full result JSON")
    run.add_argument(
        "--confirm-formal-benchmark",
        action="store_true",
        help="Start a benchmark configured with formal_mode=true",
    )

    plan = subparsers.add_parser(
        "plan-sensitivity",
        help="Write a sensitivity plan and its scenario configurations",
        description="Generate scenario configurations for a parameter sensitivity study.",
    )
    plan.add_argument("config", help="Baseline benchmark JSON configuration")
    plan.add_argument("axes", help="JSON mapping of sensitivity axes to numeric values")
    plan.add_argument("--plan-root", required=True, help="Directory for the plan and scenario configs")
    plan.add_argument("--design", choices=("one_factor", "factorial"), default="one_factor",
                      help="Vary one axis at a time or include all parameter combinations (default: one_factor)")
    plan.add_argument("--json", dest="json_output", action="store_true", help="Print full result JSON")

    sensitivity_preflight = subparsers.add_parser(
        "preflight-sensitivity",
        help="Check inputs and resources for each planned scenario",
        description="Check the configurations and inputs for every scenario in a sensitivity plan.",
    )
    sensitivity_preflight.add_argument("plan", help="sensitivity_plan.json")
    sensitivity_preflight_output = sensitivity_preflight.add_mutually_exclusive_group()
    sensitivity_preflight_output.add_argument("--json", dest="json_output", action="store_true", help="Print full JSON")
    sensitivity_preflight_output.add_argument("--output-json", help="Write the full preflight report to this path")

    sensitivity_run = subparsers.add_parser(
        "run-sensitivity",
        help="Run every scenario in a saved sensitivity plan",
        description="Run a sensitivity plan and collect the scenario results.",
    )
    sensitivity_run.add_argument("plan", help="sensitivity_plan.json")
    sensitivity_run.add_argument("--json", dest="json_output", action="store_true", help="Print full result JSON")
    sensitivity_run.add_argument(
        "--confirm-formal-benchmark",
        action="store_true",
        help="Start scenarios configured with formal_mode=true",
    )

    verify = subparsers.add_parser("verify", help="Check a completed benchmark report and its artifacts",
                                  description="Check that a benchmark report and its recorded artifacts agree.")
    verify.add_argument("report", help="Path to benchmark_report.json")
    verify.add_argument("--json", dest="json_output", action="store_true", help="Print full result JSON")
    return parser


def _read_json(path: str) -> object:
    with Path(path).open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def _resolve_path(value: object, base: Path) -> str:
    path = Path(str(value)).expanduser()
    return str(path if path.is_absolute() else (base / path).resolve())


def load_benchmark_config(path: str) -> BenchmarkConfig:
    config_path = Path(path).resolve()
    payload = _read_json(str(config_path))
    if not isinstance(payload, Mapping):
        raise ValueError("Benchmark config must be a JSON object.")
    data = dict(payload)
    allowed_fields = {item.name for item in fields(BenchmarkConfig)}
    unknown_fields = sorted(set(data) - allowed_fields)
    if unknown_fields:
        label = "field" if len(unknown_fields) == 1 else "fields"
        raise ValueError(f"Unknown benchmark config {label}: {', '.join(unknown_fields)}.")
    section_types = {
        "surface": SurfaceConfig,
        "topology": TopologyConfig,
        "parameterization": ParameterizationConfig,
        "optcuts": OptCutsConfig,
    }
    for section_name, section_type in section_types.items():
        if section_name not in data:
            continue
        section = data[section_name]
        if not isinstance(section, Mapping):
            raise ValueError(f"Benchmark config field '{section_name}' must be a JSON object.")
        unknown_section_fields = sorted(set(section) - {item.name for item in fields(section_type)})
        if unknown_section_fields:
            names = ", ".join(f"{section_name}.{name}" for name in unknown_section_fields)
            label = "field" if len(unknown_section_fields) == 1 else "fields"
            raise ValueError(f"Unknown benchmark config {label}: {names}.")
        data[section_name] = dict(section)
    base = config_path.parent
    for name in ("input_folder", "output_root", "manifest_path", "coordinate_audit_path"):
        if data.get(name):
            data[name] = _resolve_path(data[name], base)
    optcuts = dict(data.get("optcuts", {}))
    binary = str(optcuts.get("optcuts_bin") or "").strip()
    if binary:
        binary_path = Path(binary).expanduser()
        relative_candidate = base / binary_path
        if binary_path.is_absolute() or binary_path.parent != Path(".") or relative_candidate.exists():
            optcuts["optcuts_bin"] = _resolve_path(binary, base)
    data["optcuts"] = optcuts
    try:
        return benchmark_config_from_dict(data)
    except TypeError as exc:
        raise ValueError(f"Invalid benchmark config: {exc}") from None


def _emit(
    payload: object,
    *,
    json_output: bool = False,
    output_path: str | None = None,
    summary: Sequence[str] = (),
) -> None:
    if output_path:
        target = Path(output_path).resolve()
        dump_json_atomic(payload, target)
        print(f"JSON report: {target}")
        return
    if json_output:
        dump_json(payload, sys.stdout)
        sys.stdout.write("\n")
        return
    for line in summary:
        print(line)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command in {"preflight", "preflight-sensitivity"} and args.output_json:
            source = args.config if args.command == "preflight" else args.plan
            if Path(args.output_json).expanduser().resolve() == Path(source).expanduser().resolve():
                raise ValueError("Choose a preflight report path separate from the input configuration or plan.")
        if args.command == "preflight":
            output = BenchmarkRunner(load_benchmark_config(args.config)).preflight()
            state = "ready" if output["ready"] else "blocked"
            summary = [
                f"Preflight: {state}",
                (
                    f"Structures: {output['accepted_job_count']} accepted, "
                    f"{output['resumed_structure_count']} resumed, "
                    f"{output['remaining_job_count']} remaining"
                ),
                f"Worker processes: {output['planned_worker_process_count']} planned",
            ]
            summary.extend(f"Blocker: {blocker}" for blocker in output["blockers"])
            _emit(
                output,
                json_output=args.json_output,
                output_path=args.output_json,
                summary=summary,
            )
            return 0 if output["ready"] else 2
        if args.command == "run":
            config = load_benchmark_config(args.config)
            if config.formal_mode and not args.confirm_formal_benchmark:
                parser.error("formal_mode=true requires --confirm-formal-benchmark")
            runner = BenchmarkRunner(config)
            preflight = runner.preflight()
            output = runner.run()
            report_path = str(Path(config.output_root) / config.report_filename)
            result_summary = output.get("summary", {})
            attempted = int(result_summary.get("attempted_structure_count", 0))
            failed = int(result_summary.get("failed_structure_count", 0))
            summary = [
                "Benchmark: complete",
                f"Report: {report_path}",
                (f"Structures: {attempted} attempted, {preflight['resumed_structure_count']} resumed, {failed} failed"),
            ]
            if config.execution_profile.strip().lower() == "operational_optcuts":
                usable = int(result_summary.get("operational_scientifically_usable_structure_count", 0))
                summary.append(f"Scientifically usable: {usable}/{attempted}")
            else:
                complete = int(result_summary.get("complete_comparison_structure_count", 0))
                summary.append(f"Complete comparisons: {complete}/{attempted}")
            result = {
                "status": "complete",
                "report": report_path,
                "resumed_structure_count": preflight["resumed_structure_count"],
                "summary": result_summary,
            }
            _emit(
                result,
                json_output=args.json_output,
                summary=summary,
            )
            return 0
        if args.command == "plan-sensitivity":
            config = load_benchmark_config(args.config)
            axes = _read_json(args.axes)
            if not isinstance(axes, Mapping):
                raise ValueError("Sensitivity axes file must contain a JSON object.")
            plan_path = write_sensitivity_plan(
                config,
                axes,
                args.plan_root,
                design=args.design,
            )
            result = {"status": "planned_not_run", "plan": plan_path}
            _emit(
                result,
                json_output=args.json_output,
                summary=(f"Sensitivity plan: {plan_path}",),
            )
            return 0
        if args.command == "preflight-sensitivity":
            output = SensitivityBenchmarkRunner(args.plan).preflight()
            state = "ready" if output["ready"] else "blocked"
            summary = [
                f"Sensitivity preflight: {state}",
                f"Scenarios: {output['scenario_count']}",
                f"Worker processes: {output['planned_worker_process_count']} planned",
            ]
            summary.extend(
                f"Blocker [{scenario['scenario_id']}]: {blocker}"
                for scenario in output["scenarios"]
                for blocker in scenario["blockers"]
            )
            _emit(
                output,
                json_output=args.json_output,
                output_path=args.output_json,
                summary=summary,
            )
            return 0 if output["ready"] else 2
        if args.command == "run-sensitivity":
            output = SensitivityBenchmarkRunner(args.plan).run(confirm_formal=args.confirm_formal_benchmark)
            records = output.get("scenarios", [])
            failed = sum(1 for record in records if isinstance(record, Mapping) and record.get("status") == "failed")
            results_path = Path(args.plan).resolve().parent / "sensitivity_results.json"
            _emit(
                output,
                json_output=args.json_output,
                summary=(
                    "Sensitivity study: complete",
                    f"Results: {results_path}",
                    f"Scenarios: {output.get('scenario_count', 0)} total, {failed} failed",
                ),
            )
            return 0
        if args.command == "verify":
            report_path = Path(args.report).expanduser().resolve()
            report = read_json_object(report_path, "Benchmark report")
            report_sha256 = validate_benchmark_evidence_bundle(report_path, report)
            result = {
                "status": "verified",
                "report": str(report_path),
                "report_sha256": report_sha256,
                "schema_version": report.get("schema_version"),
                "topoppi_version": report.get("topoppi_version"),
            }
            _emit(
                result,
                json_output=args.json_output,
                summary=(
                    "Evidence bundle: verified",
                    f"Report: {report_path}",
                    f"SHA-256: {report_sha256}",
                ),
            )
            return 0
    except (TopoPPIError, ValueError, OSError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"topoppi-benchmark: {exc}", file=sys.stderr)
        return 2
    return 2


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ["build_parser", "load_benchmark_config", "main"]
