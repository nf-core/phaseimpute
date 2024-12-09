#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 25/05/2024
@author: LouisLeNezet
Main script to the metro maps
"""

import argparse
import pandas as pd
from lxml import html, etree

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str)

args = parser.parse_args()
file_path = args.file

with open(file_path) as f:
    page = f.read()

tree = html.fromstring(page)
pre = tree.find_class("mermaid")

mermaid = pre[0].text_content().split("\n")
mermaid = [x.strip() for x in mermaid if x.strip() != ""]

# Extract relationships from mermaid
rels = [x for x in mermaid if "-->" in x]
rels = pd.DataFrame(rels, columns=["Rel"])
rels[["From", "To"]] = rels["Rel"].str.split(" --> ", expand=True)

# Extract graph structure
graph = [x for x in mermaid if "-->" not in x]
graph = pd.DataFrame(graph, columns=["Nodes"])
names = graph["Nodes"].str.findall("^v\\d+")
graph["Names"] = [x[0] if len(x)==1 else x for x in names]

# Extract nodes to delete
to_del = []
for pat in ['[" "]']:
    to_del = to_del + [x.split(pat)[0] for x in graph["Nodes"] if pat in x]

# Delete unwanted nodes
for x in to_del:
    all_to = rels.loc[(rels["From"] == x)]
    if len(all_to) == 0:
        rels = rels.loc[(rels["To"] != x)]
        graph = graph.loc[(graph["Names"] != x)]
    elif len(all_to) == 1:
        print(x, "link to one")
    else :
        print(x, "link to more than one")

new_graph="\n".join(graph["Nodes"])
new_rels = "\n".join(rels["Rel"])

new_mermaid="\n".join(["\n",new_graph, new_rels])

with open("new.html", "wb") as f:
    tree.find_class("mermaid")[0].text = new_mermaid
    f.write(etree.tostring(tree))
