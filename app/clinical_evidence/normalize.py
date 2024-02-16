"""Simple script for normalizing kgx node and edge files."""

import httpx
import jsonlines

normalizer = "https://nodenormalization-sri.renci.org/get_normalized_nodes"
node_file_path = ""
edge_file_path = ""
output_node_file_path = "./output_normalized_nodes.jsonl"
output_edge_file_path = "./output_normalized_edges.jsonl"

node_map = {}
nodes = []

if __name__ == "__main__":
    with jsonlines.open(node_file_path) as reader:
        for node in reader:
            try:
                with httpx.Client() as client:
                    response = client.post(normalizer, json={"curies": [node["id"]]})
                    response.raise_for_status()
                    normalized_node = response.json()
                    new_node = {}
                    if normalized_node[node["id"]] is None:
                        if "category" in node and type(node["category"]) == str:
                            node["category"] = [node["category"]]
                        new_node = node
                    else:
                        new_node["id"] = normalized_node[node["id"]]["id"]["identifier"]
                        new_node["name"] = normalized_node[node["id"]]["id"].get(
                            "label", ""
                        )
                        new_node["category"] = normalized_node[node["id"]]["type"]
                    nodes.append(new_node)
                    node_map[node["id"]] = new_node["id"]
            except:
                print(node["id"])

    with open(output_node_file_path, "w") as f:
        writer = jsonlines.Writer(f)
        writer.write_all(nodes)
        writer.close()

    edges = []

    with jsonlines.open(edge_file_path) as reader:
        for edge in reader:
            edge["subject"]["id"] = node_map[edge["subject"]["id"]]
            edge["object"]["id"] = node_map[edge["object"]["id"]]
            edges.append(edge)

    with open(output_edge_file_path, "w") as f:
        writer = jsonlines.Writer(f)
        writer.write_all(edges)
        writer.close()
