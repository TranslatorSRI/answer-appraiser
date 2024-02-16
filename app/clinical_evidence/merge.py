"""Script for combining multiple KGX node and edge files."""

import glob
import json
import jsonlines


if __name__ == "__main__":
    node_files = glob.glob("./kgx/*_nodes.jsonl")
    edge_files = glob.glob("./kgx/*_edges.jsonl")

    node_ids = set()
    nodes = {}

    for node_file in node_files:
        print(f"merging {node_file}")
        with jsonlines.open(node_file) as reader:
            for node in reader:
                if node["id"] not in node_ids:
                    nodes[node["id"]] = node

    with open("nodes_merged.json", "w") as f:
        json.dump(nodes, f, indent=2)

    supported_categories = ["biolink:Drug", "biolink:ChemicalEntity", "biolink:Disease"]

    edges = {}

    for edge_file in edge_files:
        print(f"merging {edge_file}")
        with jsonlines.open(edge_file) as reader:
            for edge in reader:
                # get subject node
                if type(edge["subject"]) == dict:
                    # from ehr
                    edge["subject"] = edge["subject"]["id"]
                subject_node = nodes[edge["subject"]]

                # check to see if subject is drug/chemical or disease
                if any(
                    category in subject_node["category"]
                    for category in supported_categories
                ):
                    # only now check if object is also drug/chemical or disease
                    if type(edge["object"]) == dict:
                        # from ehr
                        edge["object"] = edge["object"]["id"]
                    object_node = nodes[edge["object"]]

                    if any(
                        category in object_node["category"]
                        for category in supported_categories
                    ):
                        # save edge
                        # TODO: we will want to save the full edge so we can put it in an aux graph,
                        # which will require TRAPIfying all the edges
                        # For now, just save the bits we need to calculate the score
                        save_edge = {
                            "subject": edge["subject"],
                            "object": edge["object"],
                            "log_odds_ratio": 0,
                            "log_odds_ratio_95_ci": [0, 0],
                            "total_sample_size": 0,
                        }
                        if "association" in edge:
                            # from ehr
                            save_edge["predicate"] = (
                                f"biolink:{edge['association']['predicate']}"
                            )
                            for attribute in edge["association"]["edge_attributes"]:
                                if (
                                    attribute["attribute_type_id"]
                                    == "biolink:supporting_data_source"
                                ):
                                    save_edge["supporting_data_source"] = attribute[
                                        "value_url"
                                    ]
                                if (
                                    attribute["attribute_type_id"]
                                    == "biolink:has_supporting_study_result"
                                ):
                                    for sub_attribute in attribute["attributes"]:
                                        if (
                                            sub_attribute["attribute_type_id"]
                                            == "biolink:log_odds_ratio"
                                        ):
                                            save_edge["log_odds_ratio"] = sub_attribute[
                                                "value"
                                            ]
                                        if (
                                            sub_attribute["attribute_type_id"]
                                            == "biolink:log_odds_ratio_95_ci"
                                            and type(sub_attribute["value"]) != str
                                        ):
                                            save_edge["log_odds_ratio_95_ci"] = (
                                                sub_attribute["value"]
                                            )
                                        if (
                                            sub_attribute["attribute_type_id"]
                                            == "biolink:total_sample_size"
                                        ):
                                            save_edge["total_sample_size"] = (
                                                sub_attribute["value"]
                                            )

                        elif "biolink:supporting_data_source" in edge:
                            # from icees
                            save_edge["supporting_data_source"] = edge[
                                "biolink:supporting_data_source"
                            ]
                            save_edge["predicate"] = edge["predicate"]
                            save_edge["log_odds_ratio"] = edge.get("log_odds_ratio", 0)
                            save_edge["log_odds_ratio_95_ci"] = edge.get(
                                "log_odds_ratio_95_ci", [0, 0]
                            )
                            save_edge["total_sample_size"] = int(
                                edge.get("total_sample_size", 0)
                            )

                        elif "log_odds_analysis_result" in edge:
                            # from cohd
                            save_edge["supporting_data_source"] = edge[
                                "supporting_data_source"
                            ]
                            save_edge["predicate"] = edge["predicate"]
                            save_edge["log_odds_ratio"] = edge[
                                "log_odds_analysis_result"
                            ].get("log_odds_ratio", 0)
                            save_edge["log_odds_ratio_95_ci"] = edge[
                                "log_odds_analysis_result"
                            ].get("log_odds_ratio_95_ci", [0, 0])
                            save_edge["total_sample_size"] = edge[
                                "log_odds_analysis_result"
                            ].get("total_sample_size", 0)

                        if save_edge["log_odds_ratio"] > 10:
                            save_edge["log_odds_ratio"] = 10
                            save_edge["log_odds_ratio_95_ci"] = [10, 10]

                        # save edge
                        if f"{save_edge['subject']}_{save_edge['object']}" in edges:
                            edges[
                                f"{save_edge['subject']}_{save_edge['object']}"
                            ].append(save_edge)
                        else:
                            edges[f"{save_edge['subject']}_{save_edge['object']}"] = [
                                save_edge
                            ]

    print("Writing output edges...")
    with open("edges_merged.json", "w") as f:
        json.dump(edges, f, indent=2)
    print("Merge Complete!")
