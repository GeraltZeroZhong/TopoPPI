"""Minimal plugin: summarize a supplied Core geometry object."""
from topoppi.plugins import PluginMetadata


def describe(geometry, *, context):
    context.check_cancelled()
    return {"partners": [list(p.chains) for p in geometry.partners],
            "atoms": [len(p.atoms) for p in geometry.partners],
            "patches": [len(p.patches) for p in geometry.partners],
            "contact_pairs": sum(len(partners) for partners in geometry.contacts().values())}


def main(argv=None):
    import argparse
    import json

    import topoppi

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structure")
    parser.add_argument("--chains-a", required=True)
    parser.add_argument("--chains-b", required=True)
    args = parser.parse_args(argv)
    interface = topoppi.get_interface(topoppi.load_structure(args.structure), args.chains_a, args.chains_b)
    print(json.dumps(topoppi.run("stats.describe", topoppi.compute_geometry(interface))))
    return 0


class StatsPlugin:
    metadata = PluginMetadata("stats", "0.1.0", "Shared interface geometry counts.", ("stats.describe",))

    def register(self, registry):
        registry.register_capability("stats.describe", "topoppi_example_stats:describe")
        registry.register_command("stats", "topoppi_example_stats:main")


plugin = StatsPlugin()
