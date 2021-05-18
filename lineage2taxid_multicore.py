#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import logging
# 3rd party
import networkx as nx

# multithreading
from multiprocessing import Pool
from functools import partial
from itertools import repeat

# from networkx.algorithms.dag import descendants
# from networkx.algorithms.lowest_common_ancestors import lowest_common_ancestor
# from networkx.algorithms.shortest_paths.unweighted import bidirectional_shortest_path

# Statics Variables
THREADS=48
# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

# argparse
desc = 'Map taxonomic lineages to taxids'
epi = """DESCRIPTION:
If you have a set of taxonomic lineages (eg., generated from GTDB-Tk)
and need then taxids for the taxonomic lineages, then this script is for you!

Example lineages: 
* d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Ruminococcus_E;s__Ruminococcus_E bromii_B
* d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola plebeius_A

Method:
* The nodes.dmp & names.dmp files are loaded as a graph
* For each lineage in the "table_file":
  * From the finest-to-coarsest taxonomic resolution:
    * Find node in graph with that taxonomic classification 
      * Note: captilization invariant
      * If only 1 node (so no duplicate naming), then return that taxid & rank
      * If not, go to next-finest taxonomic resolution
    * If no taxonomic classification can be mapped:
      * Use "NA" and warn the user

Notes:
* The input table file can be compressed via gzip or bzip2
* The input table is assumed to be tab-delimited
* The output is writted to STDOUT
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('table_file', metavar='table_file', type=str,
                    help='Tab-delim input table containing GTDB taxonomic lineages')
parser.add_argument('names_dmp', metavar='names_dmp',
                    help='NCBI names.dmp file. Only needed if providing NCBI taxids')
parser.add_argument('nodes_dmp', type=str, default=None,
                    help='NCBI nodes.dmp file. Only needed if providing NCBI taxids')
parser.add_argument('--lineage-column', type=str, default='classification',
                    help='Column name that contains the lineages')
parser.add_argument('--taxid-column', type=str, default='taxid',
                    help='Name of taxid column that will be appended to the input table')
parser.add_argument('--taxid-rank-column', type=str, default='taxid_rank',
                    help='Name of taxid-rank column that will be appended to the input table')
parser.add_argument('--version', action='version', version='0.0.1')

# functions
def _open(infile, mode='rb'):
    """
    Openning of input, regardless of compression
    """
    if infile.endswith('.bz2'):
        return bz2.open(infile, mode)
    elif infile.endswith('.gz'):
        return gzip.open(infile, mode)
    else:
        return open(infile)

def _decode(x):
    """
    Decoding input, if needed
    """
    try:
        x = x.decode('utf-8')
    except AttributeError:
        pass
    return x

def load_dmp(names_dmp_file, nodes_dmp_file):
    """
    Loading NCBI names/nodes dmp files as DAG
    Arguments:
      names_dmp_file : str, names.dmp file
      nodes_dmp_file : str, nodes.dmp file 
    Return:
      network.DiGraph object
    """
    regex = re.compile(r'\t\|\t')
    regex_sci_name = re.compile(r'scientific name\t\|')
    # nodes
    logging.info('Loading file: {}'.format(names_dmp_file))
    idx = {}    # {taxid : name}
    with open(names_dmp_file) as inF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            line = regex.split(line)
            if not regex_sci_name.match(line[-1]):
                continue
            # print(line)
            idx[int(line[0])] = line[1].lower()
    # names
    logging.info('Loading file: {}'.format(nodes_dmp_file))
    G = nx.DiGraph()
    G.add_node(1, rank = 'root', name = 'root',prt_id = 1, prt_name = 'root')
    with open(nodes_dmp_file) as inF:
        for line in inF:
            line = line.rstrip()
            if line == '':
                continue
            line = regex.split(line)
            taxid_child = int(line[0])
            taxid_parent = int(line[1])
            rank_child = line[2]
            name_child = idx[taxid_child].lower()
            name_parent = idx[taxid_parent].lower()
            # adding node
            G.add_node(taxid_child, rank = rank_child, name = name_child, prt_id = taxid_parent, prt_name = name_parent)
            # adding edge
            if taxid_parent == 1:
                G.add_edge(1, taxid_child)
            else:
                G.add_edge(taxid_parent, taxid_child)
    idx.clear()
    logging.info('  No. of nodes: {}'.format(G.number_of_nodes()))
    logging.info('  No. of edges: {}'.format(G.number_of_edges()))
    return G

def lineage2taxid(lineage, G):

    line = _decode(line).rstrip().split('\t')
    lineage = line[1]
    lineage = lineage.split(';')
    print(lineage)
    for cls in lineage[::-1]:
        print(cls)
        if '__' in cls:
            cls = cls.split('__')[1]
            print(cls)
        else:
            pass
        cls = cls.lower()
        print(cls)
        nodes = [x for x,y in G.nodes(data=True) if y['name'] == cls]
        print(nodes)
        if len(nodes) == 1:
            with open(file="gtdb_tax_final_formatted.tsv", mode="a+") as of:
                of.write('\t'.join(line + [str(nodes[0]), str(G.nodes[nodes[0]]['rank']),str( G.nodes[nodes[0]]['prt_id']),str(G.nodes[nodes[0]]['prt_name'])]))
            # status
            # if i > 0 and (i+1) % 100 == 0:
                # logging.info('  Records processed: {}'.format(i+1))
        else:
            msg = 'Could not find a taxid for lineage: {}'
            logging.warning(msg.format(';'.join(lineage)))
    # return ['NA', 'NA','NA','NA']

def parse_lineage_table(table_file, lineage_column, G,
                        taxid_column, taxid_rank_column):
    """
    Parsing lineage and finding taxid
    """
    logging.info('Parsing file: {}'.format(table_file))

    pool = Pool(THREADS)
    with _open(table_file) as inF:
        pool.starmap(lineage2taxid,zip(inF,repeat(G)))

def main(args):
    """
    Main interface
    """
    # loading dmp as DAG
    G = load_dmp(args.names_dmp, args.nodes_dmp)
    # lineage2taxid
    parse_lineage_table(args.table_file, args.lineage_column, G=G,
                        taxid_column = args.taxid_column,
                        taxid_rank_column = args.taxid_rank_column)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
