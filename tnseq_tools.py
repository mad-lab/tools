import sys
import math
import numpy


class Gene:
    """

        This is a longer explanation, which may include math with latex syntax
        Then, you need to provide optional subsection in this order (just to be
        consistent and have a uniform documentation. Nothing prevent you to
        switch the order):

          - parameters using ``:param <name>: <description>``
          - type of the parameters ``:type <name>: <description>``
          - returns using ``:returns: <description>``
          - examples (doctest)
          - seealso using ``.. seealso:: text``
          - notes using ``.. note:: text``
          - warning using ``.. warning:: text``
          - todo ``.. todo:: text``


        :param orf: ORF ID. Must be unique.
        :param name: Human readable name of the ORF.
        :param reads: Read-counts data for the ORF.
        :param position: Position of TA sites for the ORF
        :param start: Start coordinate for the ORF.
        :param end: End coordinate for the ORF.
        :param strand: Strand for the ORF.
        :type orf: string
        :type name: string
        :type reads: list
        :type start: int
        :type end: int
        :type strand: string
        :returns: Gene object with a information about the ORF.
        :rtype: Gene object

        :Example:

        >>> import tnseq_tools
        >>> G = tnseq_tools.Gene("Rv0001", "dnaA", [[0,0,0,0,1,0,32]], start=1, end=1500, strand="+")
        >>> print G
        Rv0001 (dnaA): 2,7,4

        .. warning:: orf must be unique.
        .. seealso:: :class:`Genes`
    """

    def __init__(self, orf, name, reads, position, start=0, end=0, strand=""):
        """ """

        self.orf = orf
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.reads = reads
        self.position = position
        self.tosses = self.tossify()
        self.k = self.tosses.count("1")
        self.n = len(self.tosses)
        self.r = maxrun(self.tosses)
        self.s = self.getspan()
        self.t = self.getlength()

    def __getitem__(self, i):
        """Return read-counts at index i.

        :param i: Index position
        :type i: int
        :returns: Read-counts at index i
        :rtype: Gene object
        """
        return self.reads[:, i]

    def __unicode__(self):
        """Return a string representation of the object"""
        return "%s (%s): %d, %d, %d" % (self.orf, self.name, self.k, self.n, self.r)

    def theta(self):
        """ Return the insertion density ("theta") for the gene.

        :returns: Insertion density (i.e. Num of sites with insertions/ Total sites)
        :rtype: float
        """
        return float(self.k)/self.n

    def phi(self):
        """ Return the non-insertion density ("phi") for the gene.

        :returns: Non-insertion density (i.e. 1 - Num of sites with insertions/ Total sites)
        :rtype: float
        """
        return 1.0 - self.theta()

    def total_reads(self):
        """ Return the total reads for the gene.

        :returns: Insertion density (i.e. Num of sites with insertions/ Total sites)
        :rtype: float
        """
        return numpy.sum(self.reads, 1)

    def runs(self):
        """ Return list of all the runs of consecutive non-insertions.

        :returns: List of all the runs of consecutive non-insertions
        :rtype: List
        """
        combined_reads = numpy.sum(self.reads, 0)
        runs = []
        current_r = 0
        for read in combined_reads:
            if read > 0: # If ending a run of zeros
                #if current_r > 0: # If we were in a run, add to list
                runs.append(current_r)
                current_r = 0
            else:
                current_r += 1
        # If we ended in a run, add it
        if current_r > 0:
            runs.append(current_r)
        return(runs)

    def calculate_span(self):
       
        runs = self.runs()
        if len(self.raw_data) > 0:
            runs = self.runs(include_pos=True) or [(0,0)]
            maxrun = max(runs)
            return self.get_TA_coord(maxrun[1] + maxrun[0]-1) + 2 - self.get_TA_coord(maxrun[1]) 
        else:
            return -1

    def calculate_length(self, raw_data):
        if len(self.raw_data) > 0:
            return(self.raw_data[-1][0] + 2 - self.raw_data[0][0])
        else:
            return -1



class Genes:

    def __getitem__(self, i):
        if isinstance(i, int):
            return(self.genes[i])
        if isinstance(i, basestring):
            for j in range(len(self.genes)):
                if self.genes[j].orf == i:
                    return(self.genes[j])


    def __contains__(self, item):
        for gene in self.genes:
            if gene.orf == item: return(True)
        return(False)


    def __len__(self):
        return( len(self.genes) )


    def __str__(self):
        return("Genes Object (N=%d)" % len(self.genes))


    def add(self, gene):
        self.genes.append(gene)

    
    def __init__(self, P="", M=None, L=None, NC = False, BS=False, MID = False, filetype="trash", protTable="", gff=""):
        self.path = P
        self.min_read = M
        self.lanes = L
        self.nc = NC
        self.bs = BS
        self.mid = MID
        self.filetype = filetype
        self.protTable = protTable
        self.gffFile = gff

        hash = None; orf2name = None;
        if self.gffFile:
            hash = hash_gff_genes(self.gffFile)
            orf2name = get_gff_names(self.gffFile)
        elif self.protTable:
            hash = hash_prot_genes(self.protTable)
            orf2name = get_prot_names(self.protTable)

        if filetype.lower() == "trash":
            self.genes = self.fromTrashFile()
        elif filetype.lower() == "wig":
            assert hash and orf2name, "Error: filetype %s requires a GFF file or protein table."
            self.genes = self.fromWigFile(hash, orf2name)
        elif filetype.lower() == "igv":
            assert hash and orf2name, "Error: filetype %s requires a GFF file or protein table."
            self.genes = self.fromIGVFile(hash, orf2name)
        else:
            self.genes = []
            raise Exception("Error: Unknown file type: %s" % filetype) 


    #@classmethod
    def fromTrashFile(self):
        """Returns list of Gene objects from TraSH formatted file"""
        gene_list = []; last_gene = ""; reads = []; raw_data = []; tosses = ""; orf="";
        count = 1
        for line in open(self.path, "r"):
            if line.startswith("#"): continue
            if line.startswith("\n"): continue
            tmp = line.split()
            if not line.startswith(" "):
                if raw_data and orf != "non-coding":
                    gene_list.append( Gene(orf, name, tosses, raw_data, reads, -1, -1,"") )
                elif orf == "non-coding" and self.nc:
                    #print raw_data
                    gene_list.append(Gene("non_coding_%d" % count,"-",tosses,raw_data,reads,-1,-1,""))
                    count +=1
                orf = tmp[0]; name = tmp[1]; raw_data = []; reads = []; tosses = "";
                continue
            if tmp[1] != "TA": continue
            start = int(tmp[0])
            if tmp[2].startswith("()") or tmp[2].startswith("|"): pos = 0.00
            else:
                pos = float(tmp[2][1:-1])
            row = [map(int,x.split()) for x in line[:-2].strip().split("|")[1:]]
            
            read = 0
            for i in self.lanes:
                if not self.bs:
                    read += sum(row[i-1])
                else:
                    if row[i-1][0] > 0 and row[i-1][1] > 0:
                        read += sum(row[i-1])
                    
            if read >= self.min_read: tosses+="1"
            else: tosses+="0"
            raw_data.append( (int(tmp[0]), pos, "TA", read)   )
            reads.append(read)
        if orf == "non-coding" and self.nc:
            gene_list.append(Gene("non_coding_%d" % count,"-",tosses,raw_data,reads,-1,-1,""))
        
        return(gene_list)


    #@classmethod
    def fromWigFile(self, hash, orf2name):
        """Returns list of Gene objects from Wig formatted file"""
        gene2wig = {}
        genes = []
        nc_count = 1
        for line in open(self.path):
            if line.startswith("#"): continue
            if line.startswith("variable"): continue
            if line.startswith("location"): continue
            loc,read = line.split(); loc = int(loc); read = int(read);
            gene = hash.get(loc, "non_coding_%d" % nc_count)
            if gene not in gene2wig:
                if gene.startswith("non_coding"):
                    nc_count+=1
                genes.append(gene)
                gene2wig[gene] = []
            gene2wig[gene].append((loc, read))

        gene_list = []
        for gene in genes:
            name = orf2name.get(gene, "-")
            if not self.nc and gene.startswith("non"): continue
            raw_reads = []
            tosses = ""
            reads = []
            for loc,read in gene2wig[gene]:
                if read >= self.min_read: tosses +="1"
                else: tosses +="0"
                raw_reads.append((loc, -1.00, "TA", read))
                reads.append(read)
            gene_list.append(Gene(gene, name, tosses, raw_reads, reads, -1, -1, ""))

        return gene_list 
            
   
    #@classmethod
    def fromIGVFile(self, hash, orf2name):
        """Returns list of Gene objects from IGV formatted file"""
        gene2igv = {}
        genes = []
        nc_count = 1
        for line in open(self.path):
            if line.startswith("#"): continue
            if line.startswith("Chromosome"):
                for i,h  in enumerate(line.strip().split("\t")):
                    if h.lower() == "start":
                        s_pos = i
                    if h.lower() == "reads":
                        r_pos = i
                continue
            tmp = line.strip().split("\t"); loc = int(tmp[s_pos]); read = int(tmp[r_pos])
            gene = hash.get(loc, "non_coding_%d" % nc_count)
            if gene not in gene2igv:
                if gene.startswith("non_coding"):
                    nc_count+=1
                genes.append(gene)
                gene2igv[gene] = []
            gene2igv[gene].append((loc, read))

        gene_list = []
        for gene in genes:
            name = orf2name.get(gene, "-")
            if not self.nc and gene.startswith("non"): continue
            raw_reads = []
            tosses = ""
            reads = []
            for loc,read in gene2igv[gene]:
                if read >= self.min_read: tosses +="1"
                else: tosses +="0"
                raw_reads.append((loc, -1.00, "TA", read))
                reads.append(read)
            gene_list.append(Gene(gene, name, tosses, raw_reads, reads, -1, -1, ""))

        return gene_list


        
    def k(self):
        """Returns numpy array with the number of insertions, 'k', for each gene FLORF"""
        G = len(self.genes)
        K = numpy.zeros(G, dtype=int)
        for i in xrange(G):
            K[i] = self.genes[i].k
        return(K)


    def n(self):
        """Returns numpy array with total number of TA sites, 'n', for each gene"""
        G = len(self.genes)
        N = numpy.zeros(G, dtype=int)
        for i in range(G):
            N[i] = self.genes[i].n
        return(N)


    def r(self):
        """Returns numpy array with total number of TA sites, 'r', for each gene"""
        G = len(self.genes)
        R = numpy.zeros(G, dtype=int)
        for i in xrange(G):
            R[i] = self.genes[i].r
        return(R)


    def reads(self):
        """Returns numpy array of lists containing the read counts for each gene"""
        all_reads = []
        G = len(self.genes)
        for i in xrange(G):
            all_reads.extend(self.genes[i].reads)
        return(numpy.array(all_reads))


    def theta(self):
        """Returns numpy array of insertion frequencies, 'theta', for each gene"""
        G = len(self.genes)
        theta = numpy.zeros(G)
        for i in xrange(G):
            theta[i] = self.genes[i].theta()
        return(theta)


    def phi(self):
        """Returns numpy array of non-insertion frequency, 'phi', for each gene"""
        return(1.0 - self.theta())


    def global_theta(self):
        """Returns global insertion frequency, of the library"""
        return( float(sum(self.k()))/sum(self.n()) )


    def global_phi(self):
        """Returns global non-insertion frequency, of the library"""
        return(1.0 - self.global_theta())


    def global_reads(self):
        """Returns total reads among the library"""
        reads_total = 0
        for g in self.genes:
            reads_total += g.total_reads()
        return(reads_total)


    def mid_reads(self):
        """Returns list of reads in middle region (5%-80%) in all genes"""
        mid_list = []
        for g in self.genes:
            mid_list.extend(g.mid_reads())
        return(numpy.array(mid_list))


    def tosses(self):
        """Returns list of bernoulli trials, 'tosses', representing insertions in the gene"""
        reads = self.reads()
        t = []
        for r in reads:
            if r >= self.min_read: t.append(1)
            else: t.append(0)   
        return(t)


    def runs(self):
        """Returns list of all runs in a gene"""
        reads = self.reads()
        runs = []; current_r = 0;
        for read in reads:
            if read >= self.min_read:
                if current_r > 0: runs.append(current_r)
                current_r = 0
            else:
                current_r += 1
        if current_r > 0:
            runs.append(current_r)
        return(runs)


    def sites(self):
        """Returns list of TA sites in the gene"""
        sites = []
        for gene in self.genes:
            for read in gene.raw_data:
                sites.append(read[0])
        return(sites)


    def site2reads(self):
        """Returns dictionary of ta sites to reads in the gene"""
        ta_sites = {}
        for gene in self.genes:
            for read in gene.raw_data:
                ta_sites[read[0]] = read[-1]
        return(ta_sites)


    def getpath(self):
        """Return path to data file"""
        return(self.path)



def fdr_thresholds(Z_raw, ALPHA=0.05):
    Z = sorted(Z_raw)[::-1]
    W = [1.0 - z for z in Z]
    N = len(Z)

    ess_threshold = 1.0

    INDEX = range(3, N+1)
    count = 0
    for i in INDEX:
        count +=1
        wi = 1 - Z[i-1]
        ai_n = (ALPHA*i)/N
        mean_wi = sum(W[0:i-2])/float(len(W[0:i-2]))
        delta_w = wi - mean_wi
        if delta_w > ai_n:
            ess_threshold = Z[i-1]
            break

    noness_threshold = 0
    count = 0
    INDEX = range(0, N+1)
    INDEX.sort(reverse=True)
    for i in INDEX:
        wi = Z[N-i+1]
        ai_n = (ALPHA*i)/N
        mean_wi = sum(Z[N-i+1:])/float(len(Z[N-i+1:]))
        delta_w = Z[N-i+1] - mean_wi
        count +=1
        if ai_n > delta_w:
            break
        noness_threshold = Z[N-i]

    return(ess_threshold, noness_threshold)


def hash_gff_genes(path):
    hash = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        if tmp[2] != "gene": continue
        features = dict([tuple(f.split("=")) for f in tmp[8].split(";")])
        start, end = int(tmp[3]), int(tmp[4])
        for i in range(start,end+1):
            #if i not in hash:
            #    hash[i] = features.get("ID", "missing")
            hash[i] = features.get("ID", "missing")
    return hash


def hash_prot_genes(path):
    hash = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        start, end = int(tmp[1]), int(tmp[2])
        for i in range(start,end+1):
            #if i not in hash:
            #    hash[i] = tmp[8]
            hash[i] = tmp[8]
    return hash


def get_gff_names(path):
    orf2name = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.split("\t")
        if tmp[2] != "gene": continue
        features = dict([tuple(f.split("=")) for f in tmp[8].strip().split(";")])
        orf2name[features.get("ID", "missing")] = features.get("Name", "-")
    return(orf2name)


def get_prot_names(path):
    orf2name = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        orf2name[tmp[8]] = tmp[7]
    return(orf2name)


def maxrun(lst,item="0"):
    best = 0
    i,n = 0,len(lst)
    while i<n:
        if lst[i]==item:
            j = i+1
            while j<n and lst[j]==item: j += 1
            r = j-i
            if r>best: best = r
            i = j
        else: i += 1
    return best


def isNumber(num):
    try:
        x = float(num)
        return(True)
    except:
        return(False)


def getR1(n):
    """Small Correction term. Defaults to 0.000016 for now"""
    return(0.000016)


def getR2(n):
    """Small Correction term. Defaults to 0.00006 for now"""
    return(0.00006)


def getE1(n):
    """Small Correction term. Defaults to 0.01 for now"""
    return(0.01)

    
def getE2(n):
    """Small Correction term. Defaults to 0.01 for now"""
    return(0.01)


def getGamma():
    """Euler-Mascheroni constant ~ 0.577215664901 """
    return(0.5772156649015328606)

    
def ExpectedRuns(n,p):
    """ER_n =  log(1/p)(nq) + gamma/ln(1/p) -1/2 + r1(n) + E1(n) (Schilling, 1990)"""    
    q = 1-p
    gamma = getGamma()
    r1 = getR1(n)
    E1 = getE1(n)
    A = math.log(n*q,1.0/p)
    B = gamma/math.log(1.0/p)
    ER = A + B -0.5 + r1 + E1
    return ER 
    

def VarR(n,p):
    """ VarR_n =  (pi^2)/(6*ln(1/p)^2) + 1/12 + r2(n) + E2(n) (Schilling, 1990)"""
    r2 = getR2(n)
    E2 = getE2(n)    
    A = math.pow(math.pi,2.0)/(6* math.pow(math.log(1.0/p),2.0))
    V = A + 1/12.0 + r2 + E2
    return V

    
def Gumbel(x,u,B):
    """CDF of the Gumbel distribution
     e^(-e^( (u-x)/B)) """
    return (math.exp( -1 * math.exp((u-x)/B )))
    

def trash_analysis(trash_data, p, gene2name = {}, gene2other_ess ={}):
    q = 1 - p
    results = {}
    gene_list = trash_data.keys()

    #for gene ORF in the gene2info...
    for gene in gene_list:
        
        insert_count = trash_data[gene][0]
        n = trash_data[gene][1]
        maxrun = trash_data[gene][2]

        if n == 0:
            results[gene] = [gene, gene2name.get(gene,"-"), insert_count,n,maxrun,"-","-", "-",gene2other_ess.get(gene,"-")]
        else:
            B = 1/math.log(1/p)
            u = math.log(n*q,1/p)
            results[gene] =  [gene, gene2name.get(gene,"-"), insert_count,n,maxrun, "%1.3f" % ExpectedRuns(n, p), "%1.3f" % math.pow(VarR(n, p),0.5),  "%1.5f" %(1-Gumbel(maxrun,u,B)),gene2other_ess.get(gene,"no-data")]
        
    return(results)



