import pandas as pd
from torch_geometric.data import Data
import torch
from features import atom_to_feature_vector,bond_to_feature_vector
from rdkit import Chem
import numpy as np
import csv
import pickle
import os
from tqdm import tqdm
def smiles2graph(smiles_string,inchy_key):

    file = "./pickle_molecule_graph/" + inchy_key + ".pkl"
    if os.path.exists(file):
        with open(file, "rb") as f:
            try:
                graph = pickle.load(f)
            except:
                print(file)
                exit()
        return graph
    mol = Chem.MolFromSmiles(smiles_string)

    # atoms
    atom_features_list = []

    for atom in mol.GetAtoms():
        atom_features_list.append(atom_to_feature_vector(atom))
    x = np.array(atom_features_list, dtype = np.int64)


    # bonds
    num_bond_features = 3  # bond type, bond stereo, is_conjugated
    if len(mol.GetBonds()) > 0: # mol has bonds
        edges_list = []
        edge_features_list = []
        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()

            edge_feature = bond_to_feature_vector(bond)

            # add edges in both directions
            edges_list.append((i, j))
            edge_features_list.append(edge_feature)
            edges_list.append((j, i))
            edge_features_list.append(edge_feature)

        # data.edge_index: Graph connectivity in COO format with shape [2, num_edges]
        edge_index = np.array(edges_list, dtype = np.int64).T

        # data.edge_attr: Edge feature matrix with shape [num_edges, num_edge_features]
        edge_attr = np.array(edge_features_list, dtype = np.int64)

    else:   # mol has no bonds
        edge_index = np.empty((2, 0), dtype = np.int64)
        edge_attr = np.empty((0, num_bond_features), dtype = np.int64)

    x = convert_to_single_emb(torch.tensor(x))
    graph = dict()
    graph['edge_index'] = torch.tensor(edge_index)
    graph['edge_feat'] = torch.tensor(edge_attr)
    graph['node_feat'] = x
    graph['num_nodes'] = len(x)
    with open("./pickle_molecule_graph/"+inchy_key+".pkl","wb")as f:
        pickle.dump(graph, f)
    return graph

class MoleculeDataset:
    def __init__(self,molecule_graphs,molecule_labels,molecule_ID,names,smiles=None):

        super().__init__()

        self.molecule_graphs = molecule_graphs
        self.molecule_labels = molecule_labels
        self.molecule_ID=molecule_ID
        self.names=names
        self.smiles=smiles
        self.data=[]
        for i in range(len(molecule_graphs)):
            data=Data(x=molecule_graphs[i]["node_feat"],edge_index=molecule_graphs[i]["edge_index"],
                      edge_attr=molecule_graphs[i]["edge_feat"],y=molecule_labels[i])
            data.protein_pdbID=names[i]
            data.edge_num=molecule_graphs[i]["edge_index"].shape[1]
            if self.smiles:
                data.smiles=smiles[i]
            self.data.append(data)

    def __getitem__(self,idx):
        return self.data[idx]

    def __len__(self):
        return len(self.molecule_graphs)

def convert_to_single_emb(x, offset=512):
    feature_num = x.size(1) if len(x.size()) > 1 else 1
    feature_offset = 1 + \
                     torch.arange(0, feature_num * offset, offset, dtype=torch.long)
    x = x + feature_offset
    return x


def preprocess_item(item):
    edge_attr, edge_index, x = item['edge_feat'], item['edge_index'], item['node_feat']
    # print("x",x)
    # print("edge_attr",edge_attr)
    # print("edge_index",edge_index)
    N = x.size(0)
    x = convert_to_single_emb(x)

    # node adj matrix [N, N] bool
    adj = torch.zeros([N, N], dtype=torch.bool)
    adj[edge_index[0, :], edge_index[1, :]] = True

    # edge feature here
    if len(edge_attr.size()) == 1:
        edge_attr = edge_attr[:, None]
    attn_edge_type = torch.zeros([N, N, edge_attr.size(-1)], dtype=torch.long)
    attn_edge_type[edge_index[0, :], edge_index[1, :]
    ] = convert_to_single_emb(edge_attr) + 1

    shortest_path_result, path = algos.floyd_warshall(adj.numpy())
    max_dist = np.amax(shortest_path_result)
    edge_input = algos.gen_edge_input(max_dist, path, attn_edge_type.numpy())
    spatial_pos = torch.from_numpy((shortest_path_result)).long()
    attn_bias = torch.zeros(
        [N + 1, N + 1], dtype=torch.float)  # with graph token

    ppp={}
    # combine
    ppp['node_feat'] = x
    ppp['adj'] = adj
    ppp['attn_bias'] = attn_bias
    ppp['attn_edge_type'] = attn_edge_type
    ppp['spatial_pos'] = spatial_pos
    ppp['in_degree'] = adj.long().sum(dim=1).view(-1)
    ppp['out_degree'] = adj.long().sum(dim=0).view(-1)
    ppp['edge_input'] = torch.from_numpy(edge_input).long()

    return ppp



def train_molecules():

    np.random.seed(42)
    molecule_graphs={}
    molecule_labels={}
    with open("./scaffold_split/protein_train.csv",'r')as f:
        reader=csv.reader(f)
        next(reader)
        for m,row in tqdm(enumerate(reader)):
            if row[0] not in molecule_graphs.keys():
                molecule_graphs[row[0]]=[]
                molecule_labels[row[0]]=[]
            if float(row[5])==1:
                molecule_labels[row[0]].append(1)
            elif float(row[5])==0:
                molecule_labels[row[0]].append(0)
            graph = smiles2graph(row[1],row[2])

            molecule_graphs[row[0]].append(graph)

        train_molecule_graphs=[]
        train_molecule_labels=[]
        train_names=[]
        train_molecule_ID={}
        start_index=0
        for train_name in molecule_graphs.keys():
            train_molecule_ID[train_name]=[i for i in range(start_index,start_index+len(molecule_graphs[train_name]))]
            start_index=start_index+len(molecule_graphs[train_name])
            train_names.extend([train_name for _ in range(len(molecule_graphs[train_name]))])
            train_molecule_graphs.extend(molecule_graphs[train_name])
            train_molecule_labels.extend(molecule_labels[train_name])

        train_set=MoleculeDataset(train_molecule_graphs,train_molecule_labels,train_molecule_ID,train_names)
    f.close()

    molecule_graphs={}
    molecule_labels={}
    with open("./scaffold_split/protein_inductive.csv",'r')as f1:
        reader=csv.reader(f1)
        next(reader)
        for row in tqdm(reader):
            if row[0] not in molecule_graphs.keys():
                molecule_graphs[row[0]]=[]
                molecule_labels[row[0]]=[]
            if float(row[5]) == 1:
                molecule_labels[row[0]].append(1)
            elif float(row[5]) == 0:
                molecule_labels[row[0]].append(0)
            graph = smiles2graph(row[1],row[2])
            molecule_graphs[row[0]].append(graph)

        valid_molecule_ID={}
        valid_molecule_graphs=[]
        valid_molecule_labels=[]
        valid_names=[]
        start_index=0
        for valid_name in molecule_graphs.keys():
            valid_molecule_ID[valid_name]=[i for i in range(start_index,start_index+len(molecule_graphs[valid_name]))]
            start_index=start_index+len(molecule_graphs[valid_name])
            valid_names.extend([valid_name for _ in range(len(molecule_graphs[valid_name]))])
            valid_molecule_graphs.extend(molecule_graphs[valid_name])
            valid_molecule_labels.extend(molecule_labels[valid_name])

        val_set=MoleculeDataset(valid_molecule_graphs,valid_molecule_labels,valid_molecule_ID,valid_names)
    f1.close()
    return train_set,val_set

def test_molecules():
    molecule_graphs = {}
    molecule_labels = {}
    with open("./scaffold_split/protein_transductive.csv", 'r') as f1:
        reader = csv.reader(f1)
        next(reader)
        for row in tqdm(reader):
            if row[0] not in molecule_graphs.keys():
                molecule_graphs[row[0]] = []
                molecule_labels[row[0]] = []
            if float(row[4]) < 1000:
                molecule_labels[row[0]].append(1)
            elif float(row[4]) >= 100000:
                molecule_labels[row[0]].append(0)
            else:
                continue
            graph = smiles2graph(row[1],row[2])
            molecule_graphs[row[0]].append(graph)

    test_molecule_ID = {}
    test_molecule_graphs = []
    test_molecule_labels = []
    test_names = []
    start_index = 0
    for test_name in molecule_graphs.keys():
        test_molecule_ID[test_name] = [i for i in range(start_index, start_index + len(molecule_graphs[test_name]))]
        start_index = start_index + len(molecule_graphs[test_name])
        test_names.extend([test_name for _ in range(len(molecule_graphs[test_name]))])
        test_molecule_graphs.extend(molecule_graphs[test_name])
        test_molecule_labels.extend(molecule_labels[test_name])

    test_set = MoleculeDataset(test_molecule_graphs, test_molecule_labels, test_molecule_ID, test_names)
    return test_set


