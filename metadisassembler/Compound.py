import os
import urllib.request
from copy import deepcopy

import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib.pyplot as plt

class Compound:
    """
    Daylight Molfile and SDF
        - 2010
            http://infochim.u-strasbg.fr/recherche/Download/Fragmentor/MDL_SDF.pdf
        - 2003
            http://www.daylight.com/meetings/mug05/Kappler/ctfile.pdf
    RDKit 2018.09.1
        http://www.rdkit.org/docs/index.html
    NetworkX 2.2
        https://networkx.github.io/documentation/networkx-2.2/_downloads/networkx_reference.pdf
    """

    def __init__(self):
        self.heads = []
        self.tails = []
        self.n_atoms = None
        self.n_bonds = None
        self.fit2d = False
        self.graph = nx.Graph()
        self.mol = None
        self.allsingle = False # turns all bonds to single bonds

    def _atom_block_line(self, a):
        line = "%10s%10s%10s" % (a[0], a[1], a[2])
        line += " " + a[3] + " " * (2 - len(a[3]))
        for k in range(4, 16):
            if k >= len(a):
                a.append("0")
            line += "%3d" % (int(a[k]))
        return line

    def _bond_block_line(self, a):
        line = ""
        line += "%3d%3d%3d%3d" % (int(a[0]), int(a[1]), int(a[2]), int(a[3]))
        return line

    def cut_bond(self, i, j):
        if i not in self.graph.adj.keys():
            return False
        elif j not in self.graph.adj[i].keys():
            return False
        self.graph.remove_edge(i, j)

        return True

    def get_molblock(self, shownum=False):
        if self.graph.order() == 0:
            return False
        out_str = ""
        for i, line in enumerate(self.heads):
            if i == 1:
                out_str += "                    2D\n"
            elif i == 3:
                out_str += "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n" % (
                    self.graph.order(), self.graph.size())
            else:
                out_str += line + "\n"
        if shownum:
            for node in self.graph.nodes(data=True):
                if "row" in node[1].keys():
                    new_node = deepcopy(node)
                    new_node[1]["row"][13] = new_node[1]["row"][12]
                    out_str += self._atom_block_line(new_node[1]["row"]) + "\n"
        else:
            out_str += "\n".join([self._atom_block_line(atom[1]["row"])
                                  for atom in self.graph.nodes(data=True)]) + "\n"
        out_str += "\n".join([self._bond_block_line(bond[2]["row"])
                              for bond in self.graph.edges(data=True)]) + "\n"
        out_str += "\n".join(self.tails) + "\n"
        return out_str

    def _get_partial_molblock(self, subg):
        out_str = ""
        replace = sorted(subg.nodes())
        for i, line in enumerate(self.heads):
            if i == 3:
                out_str += "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n" % (
                    subg.order(), subg.size())
            else:
                out_str += line + "\n"

        nodeline_array = []
        for node in subg.nodes(data=True):
            if "row" in node[1].keys():
                line = self._atom_block_line(node[1]["row"])
                line = line[:54] + "  0" + line[57:]
                # num = int(line.split()[12])
                num = int(line[57:60])
                nodeline_array.append([num, line])
        out_str += "\n".join([line for num, line in sorted(nodeline_array)]) + "\n"
        for edge in subg.edges(data=True):
            if "desc" not in edge[2]:
                e0 = int(edge[2]["row"][0]) - 1
                e1 = int(edge[2]["row"][1]) - 1
                if e0 > e1:
                    edge[2]["desc"] = True
                else:
                    edge[2]["desc"] = False

            else:
                if edge[2]["desc"]:
                    e1, e0 = sorted([edge[0], edge[1]])
                else:
                    e0, e1 = sorted([edge[0], edge[1]])

            edge[2]["row"][0] = str(replace.index(e0) + 1)
            edge[2]["row"][1] = str(replace.index(e1) + 1)

        if subg.size() > 0:
            out_str += "\n".join([self._bond_block_line(edge[2]["row"])
                                  for edge in subg.edges(data=True)]) + "\n"
        out_str += "\n".join(self.tails) + "\n"

        return out_str


    def _adjust_charges(self, new_cpd):
        """
        If a carbon in a new_cpd has 5 or more hands, the cpd is removed.
        Adjusts the charge of nitrogen and oxygen.
        """

        for i, node in enumerate(new_cpd.graph.nodes(data=True)):
            count = 0
            if node[1]["symbol"] == "C":
                for edge in new_cpd.graph.edges(i, data = True):
                    count += edge[2]["order"]

            if count >= 5:
                return False

        for i, node in enumerate(new_cpd.graph.nodes(data=True)):
            count = 0
            if node[1]["symbol"] == "N":
                for edge in new_cpd.graph.edges(i, data = True):
                    count += edge[2]["order"]
            if count >= 5:
                return False
            elif count == 4:
                new_cpd._cationize(i)
            elif count == 2 or count == 3:
                new_cpd._neutralize(i)

        for i, node in enumerate(new_cpd.graph.nodes(data=True)):
            count = 0
            if node[1]["symbol"] == "O":
                for edge in new_cpd.graph.edges(i, data = True):
                    count += edge[2]["order"]
            if count >= 4:
                return False
            elif count == 3:
                new_cpd._cationize(i)
            elif count == 2:
                new_cpd._neutralize(i)

        return True

    def _input_from_kegg(self, cid):
        """
        e.g., cid = "C00002"
        This method downloads a Molfile specified by cid (KEGG Compound ID) from KEGG database,
        save it as ./kegg/cid.mol, and calls self._input_molfile to generate a Compound object.
        If ./kegg/cid.mol already exists, this method does not download the same file again.
        """
        kegg_dir = "kegg"
        if not os.path.exists("./" + kegg_dir):
            os.mkdir("./" + kegg_dir)
        if not os.path.exists("%s/%s.mol" % (kegg_dir, cid)):
            url = "http://www.genome.jp/dbget-bin/www_bget?-f+m+%s" % (cid)
            urllib.request.urlretrieve(url, "%s/%s.mol" % (kegg_dir, cid))
        self._input_molfile("%s/%s.mol" % (kegg_dir, cid))
        self.fit2d = True

        return True

    def _input_from_knapsack(self, cid):
        """
        e.g., cid = "C00002657"
        This method downloads a Molfile specified by cid (KNApSAcK ID) from KNApSAcK 3D database,
        save it as ./knapsack/cid.mol, and calls self._input_molfile to generate a Compound object.
        If ./knapsack/cid.mol already exists, this method does not download the same file again.
        The molfile from KNApSAcK contains 3D coordinates, so the 2D coordinates are recalculated.
        """
        knapsack_dir = "knapsack"
        if not os.path.isdir("./" + knapsack_dir):
            os.mkdir("./" + knapsack_dir)
        if not os.path.isfile("./{}/{}.mol".format(knapsack_dir, cid)):
            url = "http://knapsack3d.sakura.ne.jp/mol3d/{}.3d.mol".format(cid)
            urllib.request.urlretrieve(url,
                                       "./{}/{}.mol".format(knapsack_dir, cid))
        mol = Chem.MolFromMolFile("./{}/{}.mol".format(knapsack_dir, cid))
        Chem.MolToMolFile(mol, "./{}/{}.mol".format(knapsack_dir, cid))

        self._input_molfile("./{}/{}.mol".format(knapsack_dir, cid))
        self._compute_2d_coords()

        return True

    def _input_inchi(self, inchi):
        """
        inputs an InChI string to generate an RDKit mol,
        which is passed to self._input_rdkmol() where a molblock is generated.
        """
        self._input_rdkmol(Chem.MolFromInchi(inchi))
        return True

    def _input_smiles(self, smiles):
        """
        inputs a SMILES string to generate an RDKit mol,
        which is passed to self._input_rdkmol() where a molblock is generated.
        """
        self._input_rdkmol(Chem.MolFromSmiles(smiles))
        return True

    def _input_molfile(self, molfile):
        """
        This method inputs the path to the Molfile (not the Molfile itself).
        The data in the Molfile is referred to as molblock.
        This method calls self._input_molblock to incorporate varias information into self.
        """
        with open(molfile, "r") as f:
            molblock = f.readlines()
        molblock[1] = "                    2D\n"

        if "H" in [x[31] for x in molblock[4:(4+int(molblock[3][:3]))]]:
            cpd_tmp = Compound()
            molblock = "".join(molblock)
            cpd_tmp._input_molblock(molblock)
            try:
                # for removing H atoms
                molblock = Chem.MolToMolBlock(Chem.MolFromMolBlock(molblock))
                self._input_molblock(molblock)
                # complement stereochemistry information
                for edge in cpd_tmp.graph.edges(data=True):
                    if edge[2]["row"][3] != "0":
                        self.graph.edges[edge[0], edge[1]]["row"][3] = edge[2]["row"][3]
            except:
                return False
        else:
            molblock = "".join(molblock)
            self._input_molblock(molblock)

        if self.allsingle:
            self.mol = Chem.MolFromMolBlock(self.get_molblock())
        else:
            self.mol = Chem.MolFromMolBlock(molblock)

        return True

    def _input_rdkmol(self, rdkmol):
        """
        inputs an RDKit mol and generates a molblock.
        """
        self.mol = rdkmol
        self._input_molblock(Chem.MolToMolBlock(rdkmol))

        if self.fit2d is False:
            self._compute_2d_coords()

        return True

    def _input_molblock(self, molblock):
        n_atoms = 0
        n_bonds = 0
        replace = []
        cnt_ignored_atom = 0
        for i, line in enumerate(molblock.split("\n")):
            if i < 4:
                if i == 3:
                    n_atoms = int(line[:3])
                    n_bonds = int(line[3:6])
                self.heads.append(line)
            elif i < (4 + n_atoms):
                a = []
                a.append(line[:10].lstrip())
                a.append(line[10:20].lstrip())
                a.append(line[20:30].lstrip())
                a.append(line[30:33].rstrip().lstrip())
                for t in range(33, 69, 3):
                    a.append(line[t:t+3].lstrip())

                if a[12] == "0":
                    a[12] = str(self.graph.order() + 1)
                if float(a[0]) != 0.0:
                    self.fit2d = True
                if a[3] == "H":
                    continue
                replace.append(i - 4)
                self.graph.add_node(
                    self.graph.order(), symbol=a[3], row=a)

            elif i < (4 + n_atoms + n_bonds):
                a = [line[0:3].strip(), line[3:6].strip(),
                     line[6:9].strip(), line[9:12].strip()]

                if a[2] == "2":
                    a[3] = "0"
                if self.allsingle:
                    a[2] = "1"

                bondidx = i - 4 - n_atoms - cnt_ignored_atom
                if int(a[0]) - 1 not in replace:
                    cnt_ignored_atom += 1
                    continue
                if int(a[1]) - 1 not in replace:
                    cnt_ignored_atom += 1
                    continue

                b = [s for s in a]
                b[0] = str(replace.index(int(a[0]) - 1) + 1)
                b[1] = str(replace.index(int(a[1]) - 1) + 1)

                if int(a[0]) > int(a[1]):
                    desc = True
                else:
                    desc = False

                self.graph.add_edge(replace.index(int(a[0]) - 1),
                                    replace.index(int(a[1]) - 1), order=int(a[2]), index=bondidx, row=b, desc=desc)
            else:
                if "CHG" in line.split():
                    continue
                self.tails.append(line)
        self.n_atoms = n_atoms
        self.n_bonds = n_bonds

        return True

    def _compute_2d_coords(self):
        """
        calculates the 2D coordinates, and set them to self by calling self._set_coordinates()
        This method is called if self.fit2d == False in
        _input_from_kegg, _input_from_knapsack, _input_inchi, _input_rdkmol, and _input_molfile.
        """
        AllChem.Compute2DCoords(self.mol)
        self._set_coordinates()
        self.fit2d = True

        return True

    def _neutralize(self, i):
        self.graph.nodes[i]["row"][5] = "0"
        return True

    def _cationize(self, i):
        self.graph.nodes[i]["row"][5] = "3"
        return True

    def _set_coordinates(self):
        """
        obtains molblock from self.mol, extracts 2D coordinates from the molblock, and sets them to self.graph.node.
        Before using this method:
        - 2D coordinates must be calculated by AllChem.Compute2DCoords(self.mol) or _compute_2d_coords(self)
        - self._input_molblock() must be called to prepare self.graph[node]["row"] in advance.
        """
        molblock = Chem.MolToMolBlock(self.mol)
        l_molblock = molblock.split("\n")

        for i in range(4, 4 + self.n_atoms):
            line = l_molblock[i].split()
            if i - 4 in self.graph.node.keys():
                self.graph.node[i - 4]["row"][0] = line[0]
                self.graph.node[i - 4]["row"][1] = line[1]
                self.graph.node[i - 4]["row"][2] = line[2]

        return True

    def _get_coordinates(self):
        """
        The "row" attribute of self.graph.nodes contain the row in the atom block in the molblock,
        where the 0th and 1st columns (0-indexed) contains X and Y coordinates.
        This method extracts these 2D coorinates from self.graph.nodes as the nexted list.
        """
        l_rows = [node[1]["row"] for node in self.graph.nodes(data=True)]
        l_coordinates = [list(map(float, row[0:2])) for row in l_rows]

        return l_coordinates

    def draw_cpd(self, imagefile="mol.png", shownum=False):
        if not self.fit2d:
            rdDepictor.Compute2DCoords(self.mol)
            self._set_coordinates()
            self.fit2d = True
        if shownum:
            mol = Chem.MolFromMolBlock(self.get_molblock(shownum=True))
            Draw.MolToFile(mol, imagefile)
        else:
            Draw.MolToFile(self.mol, imagefile)
        print(imagefile)

        return True

    def draw_cpd_with_labels(self):
        if not self.fit2d:
            rdDepictor.Compute2DCoords(self.mol)
            self._set_coordinates()
            self.fit2d = True
        pos = self._get_coordinates()
        node_color = self._get_node_colors()
        node_label = self._get_node_labels()
        edge_labels = {}
        for edge in self.graph.edges(data=True):
            edge_labels[(edge[0], edge[1])] = edge[2]["index"]
        nx.draw(self.graph, pos, node_color=node_color, alpha=0.4)
        nx.draw_networkx_labels(self.graph, pos, fontsize=6, labels=node_label)
        nx.draw_networkx_edge_labels(
            self.graph, pos, edge_labels=edge_labels, fontsize=3, font_color="b")

        plt.draw()

    def _get_node_colors(self):
        """
        obtains the list of the atomic element symbols in self.graph.nodes,
        sorts the listをsort，gives the 0-indexed integer IDs for the symbols without redundancy,
        and returns the correspondence list between the integer IDs and the symbol_list.
        """
        symbol_list = \
            [node[1]["symbol"] for node in self.graph.nodes(data=True)]

        count = 0
        d_index = dict()
        for symbol in sorted(symbol_list):
            if symbol in d_index.keys():
                continue
            else:
                d_index[symbol] = count
                count += 1

        l_index = [d_index[symbol] for symbol in symbol_list]

        return l_index

    def _get_node_labels(self, start=0):
        """
        returns dict to show the label of the node.
        must be a dict for the use of NetworkX.
        This method is bacially for solving the 0-index and 1-index problem in the molblock.
        """
        d_label = {}
        l_nodes = self.graph.nodes(data=True)
        for node in l_nodes:
            d_label[node[0]] = str(int(node[1]["row"][12]) + start - 1)

        return d_label
