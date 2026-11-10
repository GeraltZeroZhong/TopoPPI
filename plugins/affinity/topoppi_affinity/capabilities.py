"""Affinity analysis of an existing Core geometry representation."""
from __future__ import annotations


def predict(prepared, *, model=None, output=None, context=None, **options):
    """Predict using existing atoms and geometry; preserve the supplied object.

    Path-based CLI callers continue to use ``workflow.analyze``. This capability
    accepts Core geometry directly and does not reopen its original structure.
    Fitted models and feature transformations remain owned by this plugin.
    """
    from topoppi.execution import ExecutionContext

    from topoppi_affinity.workflow import analyze

    context = context or ExecutionContext()
    context.check_cancelled()
    values = dict(context.config_for("affinity"))
    values.update(options)
    configured_model = values.pop("model", None)
    configured_output = values.pop("output", context.output_dir)
    model = model if model is not None else configured_model
    output = output if output is not None else configured_output
    values.setdefault("source", prepared.source.get("source", "user"))
    values.setdefault("complex_id", prepared.source.get("complex_id"))
    values.setdefault("log", context.log)
    values.setdefault("cancel_event", context.cancel_event)
    values.setdefault("event_callback", context.progress)
    return analyze(prepared.source["path"], *[",".join(partner.chains) for partner in prepared.partners],
                   model=model, output=output, prepared=prepared, **values)
