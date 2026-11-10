"""Deprecated command-line aliases; no task implementation belongs here."""
import warnings


def extension_arguments(arguments):
    name, setup = arguments.plugin, arguments.plugin_setup
    if arguments.affinity:
        warnings.warn("--affinity is deprecated; use --plugin affinity", DeprecationWarning, stacklevel=2)
        name = name or "affinity"
    if arguments.affinity_setup:
        warnings.warn("--affinity-setup is deprecated; use --plugin-setup", DeprecationWarning, stacklevel=2)
        setup = setup or arguments.affinity_setup
    return name, setup
