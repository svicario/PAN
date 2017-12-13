import subprocess
import io,os
from os.path import splitext
import pandas as pd
from Bio import SeqIO, AlignIO

import json
import re
from io import StringIO
from Bio import Phylo
from pandas.io.json import json_normalize
import numpy as np
import sys



class Assigner:
    ## param for cutadapt
    trim_exec = "cutadapt"
    #primer of amplicon
    p5primo="ATTAGATACCCYGGTAGTCC"
    p3primo="ACGAGCTGACGACARCCATG"
    #quality score
    q = "30"
    e = "0.1"
    Zero = "1"
    # param for vsearch
    merge_exec = '/home/saverio/Scaricati/vsearch-2.6.0/bin/vsearch'
    nosingle_exec = cluster_exec = derep_exec = merge_exec
    id="0.98"
    #minimal size of cluster across samples, within sample is at least 2
    minsize="4"
    # param align
    align_exec="cmalign"
    RefAli = 'LTPs128_SSU_alignment_marzo2017.bacteria.mask.stk'
    CM = 'LTPs128_SSU_bacteria_mask.cm'
    # param build tree
    treebuild_exec="raxmlHPC-PTHREADS"
    treebuildfast_exec='/home/saverio/Documenti/Cluster/Mangel/epa/bin/epa-ng'
    dofasttree=True
    refTree='LTPs128_SSU_tree_bacteria.newick'
    #Name convention
    fastq_ext=[".fastq"]
    #SILVA DB ref
    silvaNCBITax = 'All_SILVA_fullLineage_5maggio.xml'
    silvataxid = 'LTPs128_SSU_taxId.csv'
    phyloTaxonomy = 'phyloT_generated_tree_1495789075_newick.txt'
    #Convention for PAN
    pattern = "size="

    def __init__(self, folder=""):
        assert isinstance(folder, str)
        self.folder=folder
        self.origin = os.path.dirname(__file__) + os.path.sep

    def run(self, all_in_one=True):
        if all_in_one:
            Report=self.all_in_one_folder(self.folder)
        else:
            Report=self.folders_navigator(self.folder)

        with open("Report.csv","w") as handle:
            handle.write(Report.to_csv())
        with open("Report.html","w") as handle:
            handle.write(Report.to_html())
        jplace=self.pre_pan(self.folder+"OverallCleanNoSingle.fa",fast=self.dofasttree)
        #jplace=self.folder+"RAxML_portableTree.OverallCleanNoSingle.jplace"
        self.pan(jplace)
        return Report

    @staticmethod
    def common_suffix(str1, str2):
        suffix = ""
        for s1, s2 in zip(str1, str2):
            if s1 != s2: break
            suffix += s1
        return suffix

    @staticmethod
    def info_cutadapt(Trim):
        out = ""
        flag = False
        for r in str(Trim.stdout, "utf-8").split("\n"):
            if flag:
                out += r.strip() + "\n"
            if r.find("=== Summary") == 0:
                flag = True
            elif r.find("Pairs") == 0:
                flag = False
        # print(out)
        Out = pd.read_csv(io.StringIO(out), header=None, sep=":")
        Out.iloc[:, 0] = Out.iloc[:, 0].str.replace("\(.+", "")
        Out.set_index(0, inplace=True)
        Out.iloc[:, 0] = Out.iloc[:, 0].str.replace("\(.+", "")
        Out.iloc[:, 0] = Out.iloc[:, 0].str.replace(",", "").astype("float")
        # Out.iloc[:,0]=Out.iloc[:,0].str.replace(",",".")
        return Out.iloc[:, 0]

    @staticmethod
    def info_merge(Merge, key="Merged", pos=0):
        fuffa = str(Merge.stderr, 'utf-8').split("\n")
        # print(fuffa)
        # print(key)
        temp = [x for x in fuffa if x.find(key) > -1][0]
        # print(temp)
        return float(temp.split()[pos])

    def folders_navigator(self,path="./", prefix=""):
        if path == "": path = "./"
        A=pd.DataFrame()
        for folder in os.listdir(path):
            if os.path.isdir(folder)&(folder.find(prefix)==0):
                files_list=[]
                for file in os.listdir(path+"/"+folder):
                    name, ext = splitext(file)
                    if ext in self.fastq_ext:
                        files_list.append(file)
                files_list.sort()
                a=self.process_fq(path+"/"+files_list[0],path+"/"+files_list[1])
                a["path"]=path
                A[a.name]=a
        Final = [x +"NoSingle.fa" for x in A.columns.tolist()]
        subprocess.run(" ".join(["cat"]+Final+[">", path+"/"+"OverallClean.fa"]), shell=True)
        A["Overall"]= self.cluster_filter(path+"/"+"OverallClean", A=pd.Series(), ext=".fa",table=True, misize=self.minsize)
        return A

    def all_in_one_folder(self, path="./"):
        if path=="": path="./"
        A = pd.DataFrame()
        fq=[x for x in os.listdir(path) if splitext(x)[1] in self.fastq_ext]
        Pairs={}
        for ffq in fq:
            base="_".join(ffq.split("_")[:-2])
            Pairs.setdefault(base,[]).append(ffq)
        print(Pairs)
        print("Start Trim,Merge and Cluster")
        for pair in Pairs:
            print(pair)
            pair=Pairs[pair]
            pair.sort()
            InfastaqFor, InfastaqRev=pair
            a=self.process_fq( InfastaqFor, InfastaqRev, multicore=10)
            a["path"] = path
            A[a.name] = a
        print("Accross pairs clustering")
        Final = [x + "NoSingleS.fa" for x in A.columns.tolist()]
        print(Final)
        subprocess.run(" ".join(["cat"] + Final + [">", path + "/" + "OverallClean.fasta"]), shell=True)
        A["Overall"] = self.cluster_filter(path + "/" + "OverallClean", A=pd.Series(), ext=".fasta", table=True,
                                           misize=self.minsize)
        return A



    def cluster_filter(self,base, A=None, ext=".fqD", table=False, misize="2"):
        """
        Clusterize fastq file with size annotations and filter oligo/singleton
        :param base: base of file name,
        :param A: pandas Series to add up new stat
        :param ext: extension of file to be processed
        :return: pandas Series with stats on process
        """
        argomenticlust = {"--cluster_size": base + ext, "--id": self.id, "--centroids": base + ".fa"}
        if table:
            argomenticlust["--otutabout"]=base+".tab"
        Clust = subprocess.run([self.cluster_exec] + list(sum(argomenticlust.items(), ())) + ["--sizeout", "--sizein"],
                               stderr=subprocess.PIPE)
        A["Clusters"] = self.info_merge(Clust, key="Clusters", pos=1)
        print("Cluster")
        argomentinosingle = {"--sortbysize": base + ".fa", "--relabel": "Read_", "--output": base + "NoSingle.fa",
                             "--minsize": misize}
        NoSingle = subprocess.run([self.nosingle_exec, "--sizeout"] + list(sum(argomentinosingle.items(), ())),
                                  stderr=subprocess.PIPE)
        A["NoSingle"] = len(list(SeqIO.parse(base + "NoSingle.fa", "fasta")))
        return A

    def process_fq(self,InfastaqFor, InfastaqRev, multicore=10):
        base1 = ".".join(InfastaqFor.split(".")[:-1])
        base2 = ".".join(InfastaqRev.split(".")[:-1])
        base = self.common_suffix(base1, base2)

        Cleanedfq2 = base2 + ".fqClean"
        Cleanedfq1 = base1 + ".fqClean"

        argomentiTrim = {"-g": self.p5primo, "-G": self.p3primo, "-e": self.e, "-O": self.Zero, "-o": Cleanedfq1, "-p": Cleanedfq2,
                         "-q": self.q}
        Trim = subprocess.run([self.trim_exec] + list(sum(argomentiTrim.items(), ())) + [InfastaqFor, InfastaqRev],
                              stdout=subprocess.PIPE)
        # print(Trim)
        A = self.info_cutadapt(Trim)
        print("Trimmed")

        argomentiMerge = {"--fastq_mergepairs": Cleanedfq1, "--fastq_maxee_rate": self.e  ,"--reverse": Cleanedfq2, "--threads": str(multicore),
                          "--fastqout": base + ".fq"}
        Merge = subprocess.run([self.merge_exec] + list(sum(argomentiMerge.items(), ())),
                               stderr=subprocess.PIPE)
        A["Merged"] = self.info_merge(Merge, key="Merged")
        print("Merged")
        argomentiderep = {"--derep_fulllength": base + ".fq", "--output": base + ".fqD"}
        Derep = subprocess.run([self.derep_exec, "--sizeout"] + list(sum(argomentiderep.items(), ())),
                               stderr=subprocess.PIPE)
        A["Unique"] = self.info_merge(Derep, key="unique")
        print("Derep")
        A = self.cluster_filter(base,A)
        S = SeqIO.parse(base + "NoSingle.fa", "fasta")
        handle = open(base + "NoSingleS.fa", "w")
        handle.close()
        handle = open(base + "NoSingleS.fa", "a")
        entry = True
        truebase = base.split("/")[-1]
        batchsize = 10
        batch = []
        while entry:
            try:
                entry = next(S)
            except StopIteration:
                entry = False
                break
            entry.description += "sample=" + truebase
            entry.id = entry.name = entry.description
            batch.append(entry)
            if len(batch) > batchsize:
                SeqIO.write(batch, handle, "fasta")
                batch = []
        SeqIO.write(batch, handle, "fasta")
        A.name = base
        return A

    def pre_pan(self,CleanedfqClustered, multicore=10, fast=False):
        multicore = str(multicore)
        base = ".".join(CleanedfqClustered.split(".")[:-1])
        Ali = base + ".stk"
        print("Align")
        A = subprocess.run([self.align_exec, "--mapali", self.origin+self.RefAli, "-o", Ali, self.origin+self.CM, CleanedfqClustered], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        print(A.stderr)
        #print(A.stdout)
        A = AlignIO.read(Ali, "stockholm")
        for r in A:
            r.id = r.id.replace(':', '_')
            r.id = r.id.replace(';', '_')
        AlignIO.write(A, base + '.afa', "fasta")
        print("Place on Tree")
        if fast:
            A=AlignIO.read(base+".afa", "fasta")
            Q=[x for x in A if x.id[:5]=="Read_"]
            R=[x for x in A if x.id[:5]!="Read_"]
            for r in R:
                r.description=r.name=""
            A._records=Q
            with open(base+"Q.afa","w") as handle:
                handle.write(format(A,"fasta"))
            A._records = R
            with open(base + "R.afa","w") as handle:
                handle.write(format(A, "fasta"))
            print([self.treebuildfast_exec, "-s", base + 'R.afa',
                                '-t', self.origin + self.refTree, "-q", base + 'Q.afa',"-w", os.getcwd()])
            B = subprocess.run([self.treebuildfast_exec, "-s", base + 'R.afa',
                                '-t', self.origin + self.refTree, "-q", base + 'Q.afa',"-w", os.getcwd()],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
            jplace="epa_result.jplace"
        else:
            print(" ".join([self.treebuild_exec, '-T', multicore, "-s", base + '.afa',
                            '-t', self.origin+self.refTree, '-f', 'v', '-m', 'GTRGAMMA', "-n", base]))
            B = subprocess.run([self.treebuild_exec, '-T', multicore, "-s", base + '.afa',
                            '-t', self.origin+self.refTree, '-f', 'v', '-m', 'GTRGAMMA', "-n", base], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
            print(B.stderr)
            jplace="RAxML_portableTree."+base+".jplace"
        return jplace

    def pan(self, jplace):
        """
        Phylogenetic assignment namer
        :param jplace: tree placement from RAxML /(pplacer should be fine too) of queries on SILVA bacteria tree
        :return: dict with two table one with assignment per sequence and one per name
        """
        print("Give a Name")
        #### read jplace
        outputName = jplace.split(".")[-2]
        print("read Jplace")
        with open(jplace) as json_file:
            data = json.load(json_file)

        pl = pd.io.json.json_normalize(data, 'placements')
        #jplace may have other columns
        pl=pl[["n","p"]]
        pl.rename(columns={"n": 'Query', "p": 'Placements'}, inplace=True)
        pl["Query"] = pl.Query.apply(lambda x: x[0])
        #NEW construction of pl2 faster and robust to change of jplace
        pl.set_index("Query", inplace=True)
        pl2 = pl.Placements.apply(lambda x: pd.DataFrame(x, columns=data["fields"])[["edge_num", "like_weight_ratio"]])
        for i, j in zip(pl2.index, pl2):
            j["Query"] = i
        pl2=pd.concat(pl2.tolist())

        # OLD construction of pl2
        # pl2 = pd.DataFrame()
        # for row in pl.index:
        #     col = pl.ix[row, 'Placements']
        #     N = len(col)
        #     pl2 = pl2.append(pd.DataFrame([[i, pl.ix[row, 'Query']] for i in col],
        #                                   index=[row] * N, columns=['Placements', 'Query']))
        #
        # # add 'edge_num' and 'like_weight_ratio' columns
        # pl2['edge_num'] = [i[0] for i in pl2['Placements']]
        # pl2['like_weight_ratio'] = [i[2] for i in pl2['Placements']]
        # pl2['Query'] = [i[0] for i in pl2['Query']]

        edge_num_list = pl2['edge_num'].tolist()

        #del pl2['Placements']
        pl2.set_index("edge_num", inplace=True)
        del pl2.index.name
        # list of edge_num to search in the ref_tree (json)
        ref_tree = data['tree']
        ref_tree = re.sub(r"[\}]", "]", ref_tree)
        ref_tree = re.sub(r"[\{]", "[", ref_tree)

        handle = StringIO(ref_tree)
        ref_tree = Phylo.read(handle, "newick")
        ref_tree.clade.comment = "-1"

        def get_parent(tree, child_clade):
            node_path = tree.get_path(child_clade)
            try:
                par = node_path[-2]
            except IndexError:
                par = ref_tree.clade
            return par

        parents = {}
        for i in set(edge_num_list):
            ins_node = next(ref_tree.find_elements(comment=str(i)))
            parent = get_parent(ref_tree, ins_node)
            parents.setdefault(parent, []).append(i)

        print("done")
        ### Get Silva full lineages for Bacteria
        print("Get Silva full lineages for Bacteria")
        import xml.etree.ElementTree as ET

        tree = ET.parse(self.origin+self.silvaNCBITax)
        root = tree.getroot()
        taxa = root.findall("./ROOT/taxon")

        taxonomy = {}
        # lineageToTree={}
        for i in taxa:
            tax_lin = []
            for j in i.findall("./lineage/taxon"):
                if "rank" in j.attrib:
                    tax_lin.append((j.attrib["scientificName"], j.attrib["rank"]))
            taxonomy[i.attrib["taxId"]] = tax_lin

        import time
        print("done")
        ### Silva lineages to pandas dataframe
        print("Silva lineages to pandas dataframe")

        DF = pd.DataFrame()
        frames = []
        counter = 0
        old = A = time.time()
        for t, tt in taxonomy.items():
            counter += 1
            if (counter / 1000) == (counter / 1000.0):
                new = time.time()
                # print new-A
                # print new-old
                old = new
                # print counter
                # print len(frames)
            DFtemp = pd.DataFrame(tt)
            DFtemp.set_index(1, inplace=True)
            DFtemp.columns = [t]
            DFtemp = DFtemp.transpose()
            frames.append(DFtemp)

        DF = pd.concat(frames)
        DF = DF[["superkingdom", 'phylum', 'subphylum', 'class', 'subclass', 'order', 'suborder', 'family', 'subfamily',
                 'tribe', 'genus', 'species group', 'species subgroup', 'species', 'subspecies']]

        silva = pd.read_table(self.origin + self.silvataxid, sep='\t', header=None)
        silva = silva.iloc[:, [0, 1]]
        silva.columns = ['Accession', 'TaxID']
        silva['TaxID'] = silva['TaxID'].astype(str)

        DF2 = DF.merge(silva, left_index=True, right_on='TaxID', copy=True, indicator=False, how='outer')
        DF2.set_index('Accession', inplace=True)
        del DF2['TaxID']
        print("done")
        print("Getting Parent of insertion")
        parentsName = {}
        for i in parents:
            refname = DF2.loc[i.get_terminals()[0].name, :].tolist()
            value = [set(DF2.loc[x.name, :].tolist()) for x in i.get_terminals()]
            CommonTaxonName = set.intersection(*value)
            CommonTaxonName.discard(np.nan)
            parentsName[i] = sorted(CommonTaxonName, key=lambda x: refname.index(x))[-1]
            parentsName[i] = parentsName[i], DF2.columns[refname.index(parentsName[i])]

        taxon = [[parentsName[x]] * len(e) for x, e in parents.items()]
        taxon2, rank = zip(*sum(taxon, []))
        print("Building stats of assignments")
        assign_df = pd.DataFrame({"edge_num": sum(parents.values(), []), "Assignment": taxon2, "Rank": rank})
        assign_df.set_index("edge_num", inplace=True)

        output = pl2.merge(assign_df, left_index=True, right_index=True, copy=True, indicator=False)

        grouped = output.groupby(['Query', "Assignment", "Rank"])['like_weight_ratio'].sum()

        assignment = grouped.to_frame()
        assignment = assignment.reset_index()
        assignment.set_index('Query', inplace=True, drop=False)

        ByQueryAssignment = assignment.groupby(["Query"])

        phyloT_taxonomy = Phylo.read(self.origin+self.phyloTaxonomy, 'newick')

        ##aggregating likelihood values for parent nodes by NCBI taxonomy tree
        print("aggregating likelihood values for parent nodes by NCBI taxonomy tree")
        #self.ByQueryAssignment=ByQueryAssignment
        Bits = []
        for Query, ByTaxon in ByQueryAssignment:
            #to avoid
            ByTaxon=ByTaxon.copy()
            if ByTaxon.shape[0] != 1:
                #print("mult")
                TestN = [
                    next(phyloT_taxonomy.find_elements(name=re.sub("\)", "", re.sub("\(", "", re.sub("[ ]", "_", i)))))
                    for i in ByTaxon.Assignment]
                AggregatedValue = [sum(
                    [ByTaxon.like_weight_ratio[ByTaxon.Assignment == i.name.replace('_', ' ').replace('MAC', "(MAC)")]
                     for i in TestN if x.get_path(i) is not None]).tolist()[0] for x in TestN]
                #AggregatedValue = [i.values for i in AggregatedScore]
            else:
                #print("single")
                AggregatedValue = ByTaxon.like_weight_ratio
            ByTaxon.loc[:, "AggregLikeWeightRatio"] = AggregatedValue
            Bits.append(ByTaxon)

        AssignmentF = pd.concat(Bits)

        AssignmentF['Likelihood_Weight_Ratio'] = AssignmentF['AggregLikeWeightRatio']  # .str.get(0).astype(float)
        AssignmentF.Likelihood_Weight_Ratio.fillna(AssignmentF['AggregLikeWeightRatio'], inplace=True)
        print("clean assignment stat")
        def ReadNumber(x, pattern="size="):
            match = re.search(pattern + "[0-9]+", x)
            if match:
                res = int(match.group().replace(pattern, ""))
            else:
                res = 1
            return res

        AssignmentF['nread'] = AssignmentF.Query.apply(ReadNumber, self.pattern)
        N=AssignmentF.groupby("Query").nread.first().sum()
        del AssignmentF['like_weight_ratio']
        del AssignmentF.index.name
        del AssignmentF['AggregLikeWeightRatio']
        del AssignmentF['Query']

        report = AssignmentF.iloc[:, :]
        report.sort_values(by=['nread'], ascending=False, inplace=True)

        with open(outputName + '.html', 'w') as report_file:
            report_file.write(report.to_html(justify='left'))

        with open(outputName + '.csv', 'w') as report_file:
            report_file.write(report.to_csv())



        dummy=AssignmentF.groupby(["Assignment", "Rank"])
        relnread=dummy.nread.transform(lambda x: x / x.sum())
        AssignmentF["mean read assignement"] = AssignmentF.Likelihood_Weight_Ratio * relnread.values
        AssignmentF['read_%'] = (100 * AssignmentF.nread / N).round(2)
        AssignmentF["UniqueRead"] = 1
        del AssignmentF['Likelihood_Weight_Ratio']
        Taxon_report = AssignmentF.groupby(["Assignment", "Rank"]).sum()
        Taxon_report.sort_values(by=['nread'], ascending=False, inplace=True)
        #Taxon_report

        with open("Taxon_" + outputName + '.html', 'w') as report_file:
            report_file.write(Taxon_report.to_html(justify='left'))

        with open("Taxon_" + outputName + '.csv', 'w')as report_file:
            report_file.write(Taxon_report.to_csv())


if __name__ == "__main__":
    mode = sys.argv[1]
    print(mode)
    print(__file__)
    As = Assigner()
    if mode == "FULL":
        As.run(all_in_one=True)
    elif mode == "DENOISE":
        Report = As.all_in_one_folder(As.folder)
    elif mode == "NAME":
        jplace = sys.argv[2]
        As.pan(jplace)
    else:
        d = """
        Input files are supposed to stay in the folder where the library is called
        the library assume that the system has:
        cutadapt
        vsearch 
        raxml
        infernal
        the path to get those program could be set at the start of defintion of the class in the library as all other details of the different program.
        further the python library should include:
        pandas
        biopython 
        This library work as standalone in in three modes:

        1)Denoise+Align to ref+place on SILVA tree+ Name read. input are file *.fastq pairends with illumina naming convention
        python PanPipe.py FULL
        2) Denoise only (trim+demultiplicate+cluster+filter singleton+group across sample)
        python PanPipe.py DENOISE
        3) Give a name based on NCBI taxonomy and a given jplace file based on the SILVA tree distributed with library.
        python PanPipe.py NAME jplacefilename

        """
        print(d)
