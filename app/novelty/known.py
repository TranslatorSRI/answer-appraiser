def find_known_results(message):
    inferring_sources = [
        "infores:aragorn",
        "infores:arax",
        "infores:biothings-explorer",
        "infores:improving-agent",
        "infores:robokop",
        "infores:unsecret-agent",
    ]
    known_result_ids = []
    unknown_result_ids = []
    results = message["results"]
    knowledge_graph = message["knowledge_graph"]
    for idres, result in enumerate(results):
        for analysis in result.get("analyses") or []:
            for eb in analysis["edge_bindings"].values():
                for element in eb:
                    edge_id = element["id"]
                    edge = knowledge_graph["edges"][edge_id]
                    for source in edge["sources"]:
                        if source["resource_role"] == "primary_knowledge_source":
                            if source["resource_id"] not in inferring_sources:
                                known_result_ids.append(idres)
                                break
                            else:
                                unknown_result_ids.append(idres)
    return known_result_ids, unknown_result_ids
