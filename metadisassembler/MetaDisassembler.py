#!/usr/bin/env python
import os
import sys
import glob
import re
import math
import pickle
import datetime
import argparse
from copy import deepcopy

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from cairosvg import svg2png
from PIL import Image
import networkx as nx

from .Library import Library
from .Compound import Compound


class MetaDisassembler(Library):
    """This class uses Compound class and Library class"""

    def __init__(self):
        super().__init__()
        self.useChirality = True # consideration of stereo
        self.history = [] # split history
        self.query_specific_bul = Library() # query-specific BUL
        self.detected_pbu = []
        self.seed_list = [] # combinations after one split
        self.relation = [] # fragment to fragments relationship
        self.result = []
        self.n_img = 20 # maximum number of output images
        self.output_folder = "./output" # save folder
        self.discriminate_same_compounds = True
        self.match_bu_idx = {} # the index of a building unit having matched a fragment
        self.highlight_bonds_in_each_fragment = {}
        self.color_info = {}
        self.regenerate = False
        self.time_limit = 300 # set a time limit [s]
        self.show_stereo = True # show stereochemistry
        self.output_color = False # output color allocation information


    ### -*- Initial Settings -*- ###

    def input_query(self, query):
        if re.match("^C\d{5}$", query):
            self.name = query
            self._input_from_kegg(query)
        elif re.search(".mol$", query):
            name_tmp = query.replace(".mol", "")
            if re.search("/", name_tmp):
                name_tmp = name_tmp.split("/")[-1]
            self.name = name_tmp
            self._input_molfile(query)
        elif re.search("^InChI", query):
            self.name = query.split("/")[1]
            self._input_inchi(query)
        elif re.match("^C\d{8}$", query):
            self.name = query
            self._input_from_knapsack(query)
        else:
            try:
                self._input_smiles(query)
                self.name = self.inchis[0].split("/")[1]
            except:
                print("Invalid query :-(")
                return False

        if not len(self.inchis) or not self.inchis[0]:
            print("Invalid query :-(")
            return False
        else:
            return True

    def _initialize(self):
        self._setup_frag_limit()
        self._count_n_aromatic_bonds()
        self._generate_query_specific_bul()

        return True

    def _setup_frag_limit(self):
        if self.cpds[0].n_atoms <= 4:
            self.frag_limit = self.cpds[0].n_atoms
        elif self.cpds[0].n_atoms <= 15:
            self.frag_limit = math.ceil(self.cpds[0].n_atoms / 2)
        elif self.cpds[0].n_atoms <= 27:
            self.frag_limit = math.ceil(self.cpds[0].n_atoms / 3)
        elif self.cpds[0].n_atoms <= 53:
            self.frag_limit = math.ceil(self.cpds[0].n_atoms / 4)
        else:
            self.frag_limit = math.ceil(self.cpds[0].n_atoms / 6)

        return True

    def _count_n_aromatic_bonds(self):
        cnt = 0
        for b in self.cpds[0].mol.GetBonds():
            if b.GetBondType() == rdkit.Chem.rdchem.BondType.AROMATIC:
                cnt += 1
        if cnt >= 17:
            self.ambiguity = True
        else:
            self.ambiguity = False

        return True

    def _get_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("query",
                            type=str,
                            help="MDL_Molfile, SMILES, InChI, KEGG_COMPOUND_ID, KNApSAcK_COMPOUND_ID")
        parser.add_argument("-t",
                            "--time",
                            type=int,
                            default=300,
                            help="set a time limit [s] [default: 300]")
        parser.add_argument("--hide",
                            action="store_true",
                            default=False,
                            help="hide stereochemistry [default: False]")
        parser.add_argument("-c",
                            "--color",
                            action="store_true",
                            default=False,
                            help="output color allocation information [default: False]")

        args = parser.parse_args()

        if not args.query:
            return False

        self.query = args.query
        self.time_limit = args.time
        self.show_stereo = not(args.hide)
        self.output_color = args.color

        return True

    def _generate_query_specific_bul(self):
        def return_bu_mol(bu_idx):
            bu_mol = bul.cpds[bu_idx].mol
            if self.cpds[0].mol.HasSubstructMatch(bu_mol, useChirality=self.useChirality):
                return bu_mol
            else:
                return False

        def return_mcs_mol(bu_idx, ambiguous_rate=0.2):
            bu_mol = bul.cpds[bu_idx].mol
            mcs = rdFMCS.FindMCS([self.cpds[0].mol, bu_mol],
                                 bondCompare=rdFMCS.BondCompare.CompareAny)

            if len(mcs.smartsString):
                mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
                tmp_lib = Library()
                tmp_lib._input_rdkmol(mcs_mol)
                if bul.cpds[bu_idx].n_atoms == tmp_lib.cpds[0].n_atoms:
                    cnt = 0
                    for bond in tmp_lib.cpds[0].graph.edges(data=True):
                        if int(bond[2]["row"][2]) >= 5:
                            cnt += 1
                    if float(cnt / mcs_mol.GetNumBonds()) <= ambiguous_rate:
                        inchi = tmp_lib.inchis[0]
                        if inchi == "" or inchi not in query_specific_bul.inchis:
                            return mcs_mol
            else:
                return False

        with open(os.path.dirname(os.path.abspath(__file__)) +
                  "/data/bul.pickle", "br") as f:
            bul = pickle.load(f)

        query_specific_bul = Library()
        for i in range(len(bul.cpds)):
            bu_mol = return_bu_mol(i)
            if bu_mol:
                query_specific_bul._input_rdkmol(bu_mol, name=bul.names[i])
                continue

            if self.ambiguity:
                mcs_mol = return_mcs_mol(i)
                if mcs_mol:
                    query_specific_bul._input_rdkmol(mcs_mol, name=bul.names[i])
                    continue

        cpd_name = {}
        for i, cpd in enumerate(query_specific_bul.cpds):
            cpd_name[cpd] = query_specific_bul.names[i]

        query_specific_bul_sorted = Library()
        for i, cpd in enumerate(sorted(query_specific_bul.cpds, key=lambda x: (x.n_atoms, x.n_bonds), reverse=True)):
            query_specific_bul_sorted._input_rdkmol(cpd.mol, name=cpd_name[cpd])

        self.query_specific_bul = query_specific_bul_sorted

        return True

    ### -*- Combinatorial Optimization -*- ###

    def _detect_pbu(self):
        pbu = {
            "Shikimic acid": "InChI=1S/C7H10O5/c8-4-1-3(7(11)12)2-5(9)6(4)10/h1,4-6,8-10H,2H2,(H,11,12)/t4-,5-,6-/m1/s1",
            "Betalamic acid": "InChI=1S/C9H11NO4/c1-2-5-3-6(8(11)12)10-7(4-5)9(13)14/h2-3,7,10H,4H2,1H3,(H,11,12)(H,13,14)/b5-2-/t7-/m0/s1",
            "D-Glucose": "InChI=1S/C6H12O5/c7-1-4-6(10)5(9)3(8)2-11-4/h3-10H,1-2H2/t3-,4+,5+,6+/m0/s1",
            "D-Glucuronate": "InChI=1S/C6H10O6/c7-2-1-12-5(6(10)11)4(9)3(2)8/h2-5,7-9H,1H2,(H,10,11)/t2-,3+,4-,5-/m0/s1",
            "D-Mannose": "InChI=1S/C6H12O5/c7-1-4-6(10)5(9)3(8)2-11-4/h3-10H,1-2H2/t3-,4-,5-,6-/m1/s1",
            "D-Galactose": "InChI=1S/C6H12O5/c7-1-4-6(10)5(9)3(8)2-11-4/h3-10H,1-2H2/t3-,4+,5+,6-/m0/s1",
            "L-Rhamnose": "InChI=1S/C6H12O4/c1-3-5(8)6(9)4(7)2-10-3/h3-9H,2H2,1H3/t3-,4-,5-,6-/m0/s1",
            }

        for i, inchi in enumerate(self.query_specific_bul.inchis):
            if inchi in pbu.values():
                self.detected_pbu.append([inchi, self.query_specific_bul.cpds[i].n_atoms])

        if len(self.detected_pbu):
            self.detected_pbu = sorted(self.detected_pbu, key=lambda x: x[1], reverse=True)
            return True
        else:
            return False

    def _check_fragment_is_bu(self, target_id):
        for i, cpd in enumerate(self.query_specific_bul.cpds):
            if len(self.cpds[target_id].mol.GetSubstructMatch(cpd.mol, useChirality=self.useChirality)) == self.cpds[target_id].n_atoms:
                self._get_highlight_bonds(target_id, i)
                self.match_bu_idx[target_id] = i
                if target_id == 0:
                    self.result.append({"target": [0], "cut_bond_id": []})

                return True

        return False

    def cut_selected_bonds(self, i, bondids):
        new_cpd = deepcopy(self.cpds[i])
        for bond in list(new_cpd.graph.edges(data=True)):
            if bond[2]["index"] in bondids:
                new_cpd.cut_bond(bond[0], bond[1])

        comment = sys._getframe().f_code.co_name + ": " + ",".join([str(n) for n in bondids])

        if self._split_cpd(i, new_cpd, comment):
            return True
        else:
            return False

    def _split_cpd(self, i, new_cpd, comment=""):
        new_cpds = []
        new_inchis = []
        destinations = []
        for subg in nx.connected_component_subgraphs(new_cpd.graph):
            if not new_cpd._adjust_charges(new_cpd):
                continue
            partial_molblock = new_cpd._get_partial_molblock(subg)
            new_mol = Chem.MolFromMolBlock(partial_molblock)
            try:
                new_inchi = Chem.MolToInchi(new_mol)
            except:
                return False
            new_cpd2 = Compound()
            new_cpd2.mol = new_mol
            new_cpd2._input_molblock(partial_molblock)
            new_cpd2.inchi = new_inchi
            #new_cpd2.fit2d = False
            new_cpds.append(new_cpd2)
            new_inchis.append(new_inchi)
        for new_cpd2, new_inchi in zip(new_cpds, new_inchis):
            k = self.add_cpd(i, new_cpd2, new_inchi)
            destinations.append(k)

        self.history.append(
            {"source": i, "target": destinations, "comment": comment})

        return True


    def _divide_into_fragments(self, target_id, query_id):
        target_mol = self.cpds[target_id].mol
        query_mol = self.query_specific_bul.cpds[query_id].mol

        if not target_mol.HasSubstructMatch(query_mol, useChirality=self.useChirality):
            return False

        boundary_bonds = []
        scaffold_atom_ids = target_mol.GetSubstructMatch(query_mol, useChirality=self.useChirality)
        for atom_idx in scaffold_atom_ids:
            for bond in target_mol.GetAtomWithIdx(atom_idx).GetBonds():
                if bond.GetOtherAtomIdx(atom_idx) not in scaffold_atom_ids:
                    boundary_bonds.append(bond.GetIdx())

        if not len(boundary_bonds):
            return False

        if self.cut_selected_bonds(target_id, boundary_bonds):
            return True
        else:
            return False

    def _get_highlight_bonds(self, target_id, query_id):
        target_mol = self.cpds[target_id].mol
        bu_mol = self.query_specific_bul.cpds[query_id].mol
        atom_idx_target = target_mol.GetSubstructMatch(bu_mol)
        atom_idx_bu = bu_mol.GetSubstructMatch(bu_mol)

        atom_idx_relation = {idx:atom_idx_target[i] for i, idx in enumerate(atom_idx_bu)}
        match_bonds = []

        for edge in self.query_specific_bul.cpds[query_id].graph.edges():
            a1 = atom_idx_relation[edge[0]]
            a2 = atom_idx_relation[edge[1]]
            a1_orig = int(self.cpds[target_id].graph.nodes[a1]["row"][12]) - 1
            a2_orig = int(self.cpds[target_id].graph.nodes[a2]["row"][12]) - 1
            match_bonds.append(self.cpds[0].graph.edges[a1_orig, a2_orig]["index"])

        self.highlight_bonds_in_each_fragment[target_id] = match_bonds

        return True

    def _generate_relation(self):
        result = []
        for i in range(len(self.history)):
            source = self.history[i]["source"]
            target = self.history[i]["target"]

            bonds = self.cpds[source].graph.edges(data=True)
            cut_bond_id = list(map(int, self.history[i]["comment"].split(": ")[1].split(",")))
            atom_pair = []
            for bond_id in cut_bond_id:
                for x in bonds:
                    if x[2]["index"] == bond_id:
                        a = int(self.cpds[source].graph.nodes[x[0]]["row"][12])
                        b = int(self.cpds[source].graph.nodes[x[1]]["row"][12])
                        atom_pair.append([a, b])

            cut_bond_id_modified = []
            for x in self.cpds[0].graph.edges(data=True):
                for (j, k) in atom_pair:
                    if set([x[0], x[1]]) == set([j-1, k-1]):
                        cut_bond_id_modified.append(x[2]["index"])

            if len(target):
                result.append({"source":source, "target":sorted(target),
                               "cut_bond_id":sorted(cut_bond_id_modified)})

        tmp = []
        for x in result:
            if x not in tmp:
                tmp.append(x)
        result = tmp

        return result

    def _generate_seed(self):
        dt_start = datetime.datetime.today()
        self._detect_pbu()

        matched_bu = {0:{None}} # Matched BUs for each target
        if not len(self.detected_pbu):
            for i in range(len(self.query_specific_bul.cpds)):
                if i <= math.ceil(len(self.query_specific_bul.cpds) / 2) - 1:
                    if self._divide_into_fragments(0, i):
                        for k in self.history[-1]["target"]:
                            if k not in matched_bu.keys():
                                matched_bu[k] = {i}
                            else:
                                matched_bu[k].add(i)
            conds = [True]
            conds.extend([False for _ in range(len(self.cpds) - 1)])

        else:
            init = True
            for j, x in enumerate(self.detected_pbu):
                for i in range(len(self.query_specific_bul.cpds)):
                    if x[0] == self.query_specific_bul.inchis[i]:
                        all_checked = False
                        while not all_checked:
                            if init:
                                if self._divide_into_fragments(0, i):
                                    for k in self.history[-1]["target"]:
                                        matched_bu[k] = {0}
                                    conds = [True]
                                    conds.extend([False for _ in range(len(self.cpds) - 1)])
                                    init = False
                            else:
                                for l, cpd in enumerate(self.cpds):
                                    if not conds[l]:
                                        if self._check_fragment_is_bu(l):
                                            conds[l] = True

                                        else:
                                            len_pre = len(self.cpds)
                                            if self._divide_into_fragments(l, i):
                                                conds[l] = True
                                                conds.extend([False for _ in range(len(self.cpds) - len_pre)])
                                                for k in self.history[-1]["target"]:
                                                    if k not in matched_bu.keys():
                                                        matched_bu[k] = {0}
                                                    else:
                                                        matched_bu[k].add(0)

                            tmp = False
                            cpd2 = self.query_specific_bul.cpds[i]
                            for l, cpd in enumerate(self.cpds):
                                if not tmp and not conds[l]:
                                    if len(cpd.mol.GetSubstructMatch(cpd2.mol, useChirality=self.useChirality)) == cpd.n_atoms:
                                        tmp = True
                                        break

                            if not tmp:
                                all_checked = True

                        break

        limit_line = len(self.query_specific_bul.inchis) - 1
        for i, cpd in enumerate(self.query_specific_bul.cpds):
            if cpd.n_atoms == 1:
                limit_line = i
                break

        while set(conds) != {True}:
            for i, cpd in enumerate(self.cpds):
                dt_tmp = datetime.datetime.today()
                if (dt_tmp - dt_start).seconds >= self.time_limit:
                    print("Over the time limit.")
                    return False

                if not conds[i]:
                    if self._check_fragment_is_bu(i):
                        conds[i] = True
                    else:
                        cut_id = min(list(matched_bu[i]) + [limit_line])
                        for j in range(len(self.query_specific_bul.cpds) - cut_id):
                            len_pre = len(self.cpds)
                            if self._divide_into_fragments(i, cut_id + j):
                                conds[i] = True
                                conds.extend([False for _ in range(len(self.cpds) - len_pre)])
                                for k in self.history[-1]["target"]:
                                    if k not in matched_bu.keys():
                                        matched_bu[k] = {cut_id + j}
                                    else:
                                        matched_bu[k].add(cut_id + j)

        self.relation = self._generate_relation()
        seed_list = [{"target":x["target"],
                      "cut_bond_id":x["cut_bond_id"]} for x in self.relation if x["source"] == 0]

        return seed_list

    def _optimize_seed(self):
        dt_start = datetime.datetime.today()
        def cond_check(list1, list2):
            if set(list1) <= set(list2):
                return True
            else:
                return False

        def replace_combi(seed_list):
            while len(seed_list):
                dt_tmp = datetime.datetime.today()
                if (dt_tmp - dt_start).seconds >= self.time_limit:
                    print("Over the time limit.")
                    return False
                seed = seed_list[0]["target"]

                if len(seed) > self.frag_limit:
                    seed_list.pop(0)
                    continue

                if cond_check(seed, goal_list):
                    seed_list_optimized.append(seed_list[0])
                    seed_list.pop(0)

                    tmp = min([len(lst["target"]) for lst in seed_list_optimized]) + 1
                    if tmp < self.frag_limit:
                        self.frag_limit = tmp
                    continue

                cut_bond_id_curr = deepcopy(seed_list[0]["cut_bond_id"])
                for i, s in enumerate(seed):
                    if s not in goal_list:
                        target_tmp = deepcopy(seed[0:i] + seed[(i+1):])
                        seed_list.pop(0)
                        for x in self.relation:
                            if x["source"] == s:
                                target_tmp_tmp =  deepcopy(target_tmp)
                                target_tmp_tmp.extend(x["target"])
                                target_tmp_tmp = sorted(set(target_tmp_tmp))
                                cut_bond_id_tmp = x["cut_bond_id"]
                                if target_tmp_tmp not in [x["target"] for x in seed_list]:
                                    seed_list.append({"target":target_tmp_tmp,
                                                      "cut_bond_id":sorted(cut_bond_id_curr + cut_bond_id_tmp)})
                        break

        source_list = set([self.history[i]["source"] for i in range(len(self.history))])
        goal_list = [i for i in range(len(self.cpds)) if i not in source_list]
        seed_list_copy = deepcopy(self.seed_list)
        seed_list_optimized = []

        replace_combi(seed_list_copy)

        seed_list_optimized2 = []
        for x in seed_list_optimized:
            if len(x) <= self.frag_limit:
                seed_list_optimized2.append(x)

        for x in seed_list_optimized2:
            if x not in self.result:
                self.result.append(x)

        return True

    ## -*- Refine Result -*- ##

    def _remove_unnatural_carbon(self, valence=4):
        result_tmp = deepcopy(self.result)
        result_new = []
        for i, x in enumerate(result_tmp):
            cond = True
            for k in x["target"]:
                n_atoms = self.cpds[k].n_atoms
                if n_atoms == 1 and self.cpds[k].graph.nodes(data=True)[0]["symbol"] == "C":
                    j = int(self.cpds[k].graph.nodes(data=True)[0]["row"][12]) - 1
                    count = 0
                    for edge in self.cpds[0].graph.edges(j, data = True):
                        count += edge[2]["order"]
                    if count >= valence:
                        cond = False
            if cond:
                result_new.append(x)

        self.result = result_new

        return True

    def _sort_result(self):
        result_tmp = deepcopy(self.result)
        vector_dict = {}
        n_frag_max = self.cpds[0].n_atoms
        for i, x in enumerate(result_tmp):
            vec_tmp = [n_frag_max - len(x["target"])]
            n_atoms_list = [self.cpds[k].n_atoms for k in x["target"]]
            vec_tmp.extend(sorted(n_atoms_list, reverse=True))
            vec_tmp.extend([0] * (n_frag_max + 1 - len(vec_tmp)))
            vec_tmp[0], vec_tmp[1] = vec_tmp[1], vec_tmp[0]
            vector_dict[i] = vec_tmp

        self.result = [result_tmp[i] for i, x in sorted(vector_dict.items(), key=lambda x:x[1], reverse=True)]

        return True

    ### -*- Output Image-*- ###

    def _generate_svg(self, result_id, params):
        color_list = [(1.00,0.50,0.50), (0.50,1.00,0.50), (0.50,1.00,1.00),
                      (1.00,0.95,0.50), (1.00,0.60,1.00), (0.60,0.60,1.00),
                      (1.00,0.70,0.70), (0.70,1.00,0.70), (0.70,1.00,1.00),
                      (1.00,0.95,0.70), (1.00,0.85,1.00), (0.85,0.85,1.00),
                      (1.00,0.85,0.85), (0.85,1.00,0.85), (0.85,1.00,1.00),
                      (1.00,0.95,0.85), (1.00,0.95,1.00), (0.95,0.95,1.00),
                      (1.00,0.95,0.95), (0.95,1.00,0.95), (0.95,1.00,1.00),
                      (1.00,0.95,0.95)]

        def thicken_lines_and_letters(svg):
            out_svg = ""
            for line in svg.split("\n"):
                if "path" in line:
                    if "stroke-width:2px" in line:
                        out_svg += re.sub(":2px", ":5px", line) + "\n"
                    else:
                        out_svg += re.sub(":1px", ":2px", line) + "\n"
                elif "text" in line:
                    out_svg += re.sub("normal", "bold", line) + "\n"
                else:
                    out_svg += line + "\n"

            return out_svg

        def adjust_thickness_of_lines_and_letters(svg):
            out_svg = ""
            for line in svg.split("\n"):
                if "path" in line:
                    out_svg += re.sub(":1px", ":2px", line) + "\n"
                elif "text" in line:
                    out_svg += re.sub("normal", "bold", line) + "\n"
                else:
                    out_svg += line + "\n"

            return out_svg

        def blacken_lines_and_letters(svg):
            out_svg = ""
            for line in svg.split("\n"):
                if "text" in line or re.search("stroke-width:[2-5]px", line):
                    out_svg += re.sub("#[a-fA-F0-9]{6}", "#000000", line)+ "\n"
                else:
                    out_svg += line + "\n"

            return out_svg

        def adjust_highlight_size(svg):
            r = 0
            for x in svg.split("\n"):
                if "ellipse" in x:
                    r = float(x.split("rx='")[1].split("' ry")[0])
                    break

            if r > 0:
                out_svg = ""
                for line in svg.split("\n"):
                    line = re.sub(":16px", ":" + str(r * 1.8) + "px", line)
                    out_svg += line + "\n"
            else:
                out_svg = svg

            return out_svg

        def take_corner(svg):
            out_svg = ""
            for line in svg.split("\n"):
                line = re.sub("miter", "round", line)
                line = re.sub("butt", "round", line)
                out_svg += line + "\n"

            return out_svg

        if self.regenerate:
            pass
        else:
            frag = sorted(self.result[result_id]["target"], key=lambda x: self.cpds[x].n_atoms, reverse=True)
            atom_list = [[int(atom[1]["row"][12]) - 1 for atom in self.cpds[i].graph.nodes(data=True)] for i in frag]
            bond_list = [self.highlight_bonds_in_each_fragment[i] for i in frag]

            color_atoms = {}
            for i, atom_tmp in enumerate(atom_list):
                for k in atom_tmp:
                    color_atoms[k] = color_list[i]
            color_bonds = {}
            for i, bond_tmp in enumerate(bond_list):
                for k in bond_tmp:
                    color_bonds[k] = color_list[i]

            highlight_atoms = range(self.cpds[0].n_atoms)
            highlight_bonds = set([e for inner_list in bond_list for e in inner_list])

            self.color_info[result_id] = {"color_atoms": color_atoms,
                                          "color_bonds": color_bonds,
                                          "highlight_atoms": highlight_atoms,
                                          "highlight_bonds": highlight_bonds}

        mol = self.cpds[0].mol
        drawer = rdMolDraw2D.MolDraw2DSVG(params["fig_size"][0], params["fig_size"][1])
        mol = rdMolDraw2D.PrepareMolForDrawing(mol,
                                               addChiralHs=False,
                                               wedgeBonds=self.show_stereo)
        drawer.SetFontSize(params["fontsize"])
        option = drawer.drawOptions()
        option.additionalAtomLabelPadding = 0.1
        option.padding = params["padding"]

        drawer.DrawMolecule(
            mol,
            highlightAtoms=self.color_info[result_id]["highlight_atoms"],
            highlightBonds=self.color_info[result_id]["highlight_bonds"],
            highlightAtomColors=self.color_info[result_id]["color_atoms"],
            highlightBondColors=self.color_info[result_id]["color_bonds"])

        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        svg = thicken_lines_and_letters(svg)
        svg = blacken_lines_and_letters(svg)
        svg = adjust_highlight_size(svg)
        svg = take_corner(svg)

        return svg

    def _generate_img(self):
        def set_params():
            params = {"fig_size": (600, 600),
                      "padding": 0.2,
                      "fontsize": 0.8}

            if self.cpds[0].n_atoms <= 7:
                params["fig_size"] = (600, 300)
            elif self.cpds[0].n_atoms <= 18:
                if self.cpds[0].mol.GetRingInfo().IsAtomInRingOfSize(5, 6):
                    params["fontsize"] = 0.8
            elif self.cpds[0].n_atoms >= 40 and self.cpds[0].n_atoms < 80:
                params["fig_size"] = (1000, 1000)
                params["fontsize"] = 0.8
            elif self.cpds[0].n_atoms >= 80:
                params["fig_size"] = (2000, 2000)
                params["fontsize"] = 1
            else:
                pass

            return params

        def expand2square(pil_img, background_color):
            width, height = pil_img.size
            if width == height:
                return pil_img
            elif width > height:
                result = Image.new(pil_img.mode, (width, width), background_color)
                result.paste(pil_img, (0, (width - height) // 2))
                return result
            else:
                result = Image.new(pil_img.mode, (height, height), background_color)
                result.paste(pil_img, ((height - width) // 2, 0))

                return result

        def white2transparence(png):
            org = Image.open(png)
            width, height = org.size
            if width != height:
                org = expand2square(org, (255, 255, 255))
            trans = Image.new("RGBA", org.size, (0, 0, 0, 0))
            width, height = org.size
            for x in range(width):
                for y in range(height):
                    px = org.getpixel((x, y))
                    if px[0] == 255 and px[1] == 255 and px[2] == 255:
                        continue
                    trans.putpixel((x, y), px)
            # trans = trans.resize(600, 600)
            trans.save(png)

            return True

        if not os.path.exists(self.output_folder):
            os.mkdir(self.output_folder)
        if not os.path.exists(self.output_folder + "/" + self.name):
            os.mkdir(self.output_folder + "/" + self.name)

        params = set_params()

        for i in range(min(len(self.result), self.n_img)):
            svg = self._generate_svg(i, params)
            img_path = self.output_folder + "/" + self.name + "/" + str(i) + ".png"
            svg2png(bytestring=svg, write_to=img_path)
            white2transparence(img_path)

        return True

    def regenerate_img(self):

        color_path = self.output_folder + '/' + self.name + '/color.pickle'

        try:
            with open(color_path, 'br') as f:
                self.color_info = pickle.load(f)
        except:
            return False

        self.result = [[] for _ in range(len(self.color_info))]
        self.regenerate = True
        self._generate_img()

        return True

    ### -*- Result Output -*- ###

    def _output_txt(self):
        if not os.path.exists(self.output_folder):
            os.mkdir(self.output_folder)
        if not os.path.exists(self.output_folder + "/" + self.name):
            os.mkdir(self.output_folder + "/" + self.name)

        f = open(self.output_folder + "/" + self.name + "/result.txt", "w")
        f.write("DATE            : {0:%Y/%m/%d %H:%M}\n".format(datetime.datetime.now()))
        f.write("COMPOUND        : " + self.name + "\n")
        f.write("RUN TIME[s]     : " + str(self.runtime) + "\n")
        if len(self.result):
            f.write("MIN # FRAGMENTS : " + str(len(self.result[0]["target"])) + "\n")
        else:
             f.write("MIN # FRAGMENTS : -\n")
        f.write("\n")
        f.write("RESULTS (Total " + str(len(self.result)) +"):\n")

        if not len(self.result):
            f.write("No results\n")
            f.write("\n")
        else:
            for i, x in enumerate(self.result):
                f.write(">>> " + str(i) + "\n")
                f.write("# fragments : " + str(len(x["target"])) + "\n")
                f.write("Cut_bond_ID : " + str(x["cut_bond_id"]) + "\n")
                for j, k in enumerate(x["target"]):
                    f.write(str(j) + " : " + str(self.query_specific_bul.names[self.match_bu_idx[k]]) + "\n")
                f.write("\n")

        f.write("END\n")
        f.close()

        return True

    def _output_molfile(self):
        f =  open(self.output_folder + "/" + self.name + "/" + self.name + ".mol", "w")
        f.write(self.cpds[0].get_molblock())
        f.close()

        return True

    def _output_color_info(self):
        f = open(self.output_folder + "/" + self.name + "/color.pickle", "bw")
        pickle.dump(self.color_info, f)
        f.close()

        return True

    def output_matched_bu(self, result_id=0):
        combi = self.result[result_id]['target']
        bus = []
        for frag in combi:
            bu = []
            bu_idx = self.match_bu_idx[frag]
            cpd = self.query_specific_bul.cpds[bu_idx]
            bu_n_atom = cpd.n_atoms
            bu_id = self.query_specific_bul.names[bu_idx]
            bu = [bu_n_atom, cpd.mol, bu_id]
            bus.append(bu)

        bus = sorted(bus, key=lambda x: x[0], reverse=True)

        bu_info = []
        for bu in bus:
            bu_info.append({
                'mol': bu[1],
                'bu_id': bu[2]
                }
            )

        return bu_info

    ### -*- Disassemble -*- ###

    def disassemble(self):
        dt_start = datetime.datetime.now()

        self._initialize()

        if not self._check_fragment_is_bu(0):
            self.seed_list = self._generate_seed()

            if not self.seed_list:
                print("No results.")
                return False

            self._optimize_seed()

            if self.cpds[0].mol.HasSubstructMatch(Chem.MolFromSmiles("C=CO")):
                self._remove_unnatural_carbon(valence=3)

            self._sort_result()

        self._generate_img()
        self.runtime = round((datetime.datetime.now() - dt_start).total_seconds(), 2)
        self._output_txt()
        self._output_molfile()
        if self.output_color:
            self._output_color_info()

        print("Finished.")

        return True
