"""Helpers for reading TRAPI node/edge bindings across schema versions.

The structure of a binding value (the value mapped to a single query-graph
node/edge key inside a result's ``node_bindings`` or an analysis's
``edge_bindings``) changed in TRAPI 1.5.0:

* TRAPI <= 1.4.x: a list of objects, each carrying a single ``id`` --
  e.g. ``[{"id": "MONDO:0005148"}]``
* TRAPI >= 1.5.0: a single object carrying an ``ids`` list --
  e.g. ``{"ids": ["MONDO:0005148"]}``

These helpers normalize both shapes so callers always get a flat list of
Knowledge Graph identifiers regardless of the TRAPI version of the message.
"""


def binding_ids(binding):
    """Return the list of Knowledge Graph identifiers for a single binding value.

    ``binding`` is the value mapped to one query-graph node/edge key, i.e. an
    entry of ``result["node_bindings"]`` or ``analysis["edge_bindings"]``.
    Tolerant of both the pre-1.5 (list of ``{"id": ...}``) and 1.5+
    (``{"ids": [...]}``) formats.
    """
    if binding is None:
        return []
    # TRAPI >= 1.5.0: {"ids": [...]}
    if isinstance(binding, dict):
        return list(binding.get("ids", []))
    # TRAPI <= 1.4.x: [{"id": ...}, ...]
    ids = []
    for entry in binding:
        if isinstance(entry, dict) and "id" in entry:
            ids.append(entry["id"])
    return ids
