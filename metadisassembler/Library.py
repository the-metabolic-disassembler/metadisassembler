import pandas as pd
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdDepictor

from .Compound import Compound

class Library:
    """This class uses the Compound class"""

    def __init__(self):

        self.names = []
        self.cpds = []
        self.inchis = []
        self.newcomers = []
        self.fps = []
        self.digraph = nx.DiGraph()
        self.generations = []
        self.discriminate_same_compounds = False

    def add_cpd(self, i, new_cpd, new_inchi):
        if not self.discriminate_same_compounds:
            if new_inchi in self.inchis:
                k = self.inchis.index(new_inchi)
                if k == i:
                    return False
                self.digraph.add_edge(i, k)
                return k
            else:
                self.digraph.add_node(len(self.cpds))
                self.digraph.add_edge(i, len(self.cpds))
                self.newcomers.append(len(self.cpds))
                self.inchis.append(new_inchi)
                self.generations.append(self.generations[i] + 1)
                self.cpds.append(new_cpd)
                self.names.append(len(self.cpds) - 1)
                return len(self.cpds) - 1

        else:
            if new_inchi in self.inchis:
                k = self.inchis.index(new_inchi)
                if k == i:
                    return False
                new_cpd_atoms = sorted([int(x[1]["row"][12]) for x in new_cpd.graph.nodes(data=True)])
                for k in range(len(self.inchis)):
                    if new_cpd_atoms == sorted([int(x[1]["row"][12]) for x in self.cpds[k].graph.nodes(data=True)]):
                        self.digraph.add_edge(i, k)
                        return k

            self.digraph.add_node(len(self.cpds))
            self.digraph.add_edge(i, len(self.cpds))
            self.newcomers.append(len(self.cpds))
            self.inchis.append(new_inchi)
            self.generations.append(self.generations[i] + 1)
            self.cpds.append(new_cpd)
            self.names.append(len(self.cpds) - 1)
            return len(self.cpds) - 1

    def _append_cpd(self, cpd, name):
        """
        appends cpd to Library.
        """
        if name is None:
            name = len(self.cpds)
        self.digraph.add_node(len(self.cpds))

        try:
            inchi = Chem.MolToInchi(cpd.mol)
        except:
            #print("error ", name)
            return False

        self.inchis.append(inchi)
        self.cpds.append(cpd)
        self.generations.append(0)
        self.names.append(name)

        return True

    def calc_fingerprints(self, fingerprint="MorganFingerprint"):
        """
        converts mols into fingerprints

        Fingerprint
            - PatternFingerprint
            - RDKFingerprint
            - MorganFingerprint
            - LayeredFingerprint
            - PatternFingerprint
            - MACCSKeysFingerprint

        http://rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html
        https://www.rdkit.org/UGM/2012/Landrum_RDKit_UGM.Fingerprints.Final.pptx.pdf
        """
        self.fps = []
        for cpd in self.cpds:
            if fingerprint == "PatternFingerprint":
                fp = Chem.PatternFingerprint(cpd.mol, fpSize=1024)
            elif fingerprint == "RDKFingerprint":
                fp = Chem.RDKFingerprint(cpd.mol)
            elif fingerprint == "MorganFingerprint":
                fp = AllChem.GetMorganFingerprintAsBitVect(cpd.mol, 2, nBits=1024)
            elif fingerprint == "LayeredFingerprint":
                fp = Chem.LayeredFingerprint(
                    cpd.mol,
                    layerFlags=Chem.LayeredFingerprint_substructLayers)
            elif fingerprint == "PatternFingerprint":
                fp = Chem.PatternFingerprint(cpd.mol, fpSize=1024)
            elif fingerprint == "MACCSKeys":
                fp = AllChem.GetMACCSKeysFingerprint(m)
            else:
                return False

            self.fps.append(fp)

        return True

    def info(self):
        df = pd.DataFrame(
            columns=["index", "#pa", "parents", "inchi", "#ch", "children"])
        for nodeidx in self.digraph.nodes():
            df2 = pd.DataFrame([[
                nodeidx,
                len(list(self.digraph.predecessors(nodeidx))),
                list(self.digraph.predecessors(nodeidx)),
                self.inchis[nodeidx],
                len(list(self.digraph.successors(nodeidx))),
                list(self.digraph.successors(nodeidx))
                ]],
                columns=df.columns
            )

            df = df.append(df2, ignore_index=True)

        return df

    def _input_from_kegg(self, cid, name=None):
        cpd = Compound()
        cpd._input_from_kegg(cid)
        if not name:
            name = cid
        if self._append_cpd(cpd, name=name):
            return True
        else:
            return False

    def _input_from_knapsack(self, cid, name=None):
        cpd = Compound()
        cpd._input_from_knapsack(cid)
        if not name:
            name = cid
        if self._append_cpd(cpd, name=name):
            return True
        else:
            return False

    def _input_inchi(self, inchi, name=None):
        cpd = Compound()
        cpd._input_inchi(inchi)
        if self._append_cpd(cpd, name):
            return True
        else:
            return False

    def _input_molfile(self, molfile, name=None):
        cpd = Compound()
        cpd._input_molfile(molfile)
        if not name:
            name = molfile.split("/")[-1].split(".")[0]
        if self._append_cpd(cpd, name=name):
            return True
        else:
            return False

    def _input_smiles(self, smiles, name=None):
        cpd = Compound()
        cpd._input_smiles(smiles)
        if self._append_cpd(cpd, name):
            return True
        else:
            return True

    def _input_rdkmol(self, mol, name=None):
        cpd = Compound()
        cpd._input_rdkmol(mol)
        self._append_cpd(cpd, name)

        return True
