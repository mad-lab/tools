import sys
import math
import numpy
import scipy.stats
import transit_tools

class Gene:
    """This is a longer explanation, which may include math with latex syntax
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


    Attributes:
        orf (str): ORF ID. Must be unique.
        name (str): Human readable name of the ORF.
        reads (list): Read-counts data for the ORF.
        position (list): Position of TA sites for the ORF.
        start (int): Start coordinate for the ORF.
        end (int): End coordinate for the ORF.
        strand (str): Strand for the ORF.

    :Example:

        >>> import tnseq_tools
        >>> G = tnseq_tools.Gene("Rv0001", "dnaA", [[0,0,0,0,1,0,32]], start=1, end=1500, strand="+")
        >>> print G
        Rv0001 (dnaA): 2,7,4

        .. warning:: orf must be unique.
        .. seealso:: :class:`Genes`
    """

    def __init__(self, orf, name, desc, reads, position, start=0, end=0, strand=""):
        """Initializes the Gene object."""

        self.orf = orf
        self.name = name
        self.desc = desc
        self.start = start
        self.end = end
        self.strand = strand
        self.reads = numpy.array(reads)
        self.position = numpy.array(position, dtype=int)
        self.tosses = tossify(self.reads)
        self.runs = runs(self.tosses)
        self.k = int(numpy.sum(self.tosses))
        self.n = len(self.tosses)
        self.r = numpy.max(self.runs)
        self.s = self.getspan()
        self.t = self.getlength()

    def __getitem__(self, i):
        """Return read-counts at position i."""
        return self.reads[:, i]

    def __str__(self):
        """Return a string representation of the object."""
        return "%s\t(%s)\tk=%d\tn=%d\tr=%d\ttheta=%1.5f" % (self.orf, self.name, self.k, self.n, self.r, self.theta())


    def getspan(self):
        """Returns the span of the maxrun of the gene (i.e. number of nucleotides)."""
        if len(self.position) > 0:
            index = runindex(self.runs)
            maxii = numpy.argmax(self.runs)
            runstart = index[maxii]
            runend = runstart + max(self.runs) - 1
            return self.position[runend] - self.position[runstart] + 2
        else:
            return 0
        
    
    def getlength(self):
        """Returns the number of nucleotides spanned by the gene."""
        if len(self.position) > 0:
            return self.position[-1] - self.position[0] + 2
        return 0


    def theta(self):
        """Return the insertion density ("theta") for the gene."""
        if self.n:
            return float(self.k)/self.n
        else:
            return 0.0

    def phi(self):
        """ Return the non-insertion density ("phi") for the gene."""
        return 1.0 - self.theta()

    def total_reads(self):
        """ Return the total reads for the gene."""
        return numpy.sum(self.reads, 1)


    def calculate_span(self):
        """Caclulates the span based on the coordinates"""
        # TODO: Check if it works.
        runs = self.runs()
        if len(self.raw_data) > 0:
            runs = self.runs(include_pos=True) or [(0,0)]
            maxrun = max(runs)
            return self.get_TA_coord(maxrun[1] + maxrun[0]-1) + 2 - self.get_TA_coord(maxrun[1]) 
        else:
            return -1

    def calculate_length(self, raw_data):
        """Caclulates the length based on the coordinates"""
        # TODO: Check if it works.
        if len(self.raw_data) > 0:
            return self.raw_data[-1][0] + 2 - self.raw_data[0][0]
        else:
            return -1



class Genes:

    def __getitem__(self, i):
        """Defines __getitem__ method so that it works as dictionary and list."""
        if isinstance(i, int):
            return(self.genes[i])

        if isinstance(i, basestring):
            return self.genes[self.orf2index[i]]


    def __contains__(self, item):
        """Defines __contains__ to check if gene exists in the list."""
        return item in self.orf2index


    def __len__(self):
        """Defines __len__ returning number of genes."""
        return len(self.genes)


    def __str__(self):
        """Defines __str__ to print a generic str with the size of the list."""
        return "Genes Object (N=%d)" % len(self.genes)
    
    def __init__(self, wigList, protTable, norm="nonorm", minread=1, ignoreCodon = True, nterm=0.0, cterm=0.0, include_nc = False):
        """Initializes the gene list based on the list of wig files and a prot_table."""
        self.wigList = wigList
        self.protTable = protTable
        self.norm = norm
        self.minread = minread
        self.ignoreCodon = ignoreCodon
        self.nterm = nterm
        self.cterm = cterm
        self.include_nc = include_nc

        self.orf2index = {}
        self.genes = []
        

        orf2info = transit_tools.get_gene_info(self.protTable)
        (data, position) = transit_tools.get_data(self.wigList)
        hash = transit_tools.get_pos_hash(self.protTable)
        (data, factors) = transit_tools.normalize_data(data, norm, self.wigList, self.protTable)
        
        K,N = data.shape

        orf2posindex = {}
        visited_list = []
        for i in range(N):
            genes_with_coord = hash.get(position[i], [])
            for gene in genes_with_coord:
                if gene not in orf2posindex: visited_list.append(gene)
                if gene not in orf2posindex: orf2posindex[gene] = []

                name,desc,start,end,strand = orf2info.get(gene, ["", "", 0, 0, "+"])
                
                if strand == "+":
                    if self.ignoreCodon and position[i] > end - 3:
                        continue
                else:
                    if self.ignoreCodon and position[i] < start + 3:
                        continue

                if (position[i]-start)/float(end-start) < (self.nterm/100.0):
                    continue
                
                if (position[i]-start)/float(end-start) > ((100-self.cterm)/100.0):
                    continue

                orf2posindex[gene].append(i)

        count = 0
        for line in open(self.protTable):
            tmp = line.split("\t")
            gene = tmp[8]
            name,desc,start,end,strand = orf2info.get(gene, ["", "", 0, 0, "+"])
            posindex = orf2posindex.get(gene, [])
            if posindex:
                pos_start = orf2posindex[gene][0]
                pos_end = orf2posindex[gene][-1]
                self.genes.append(Gene(gene, name, desc, data[:, pos_start:pos_end+1], position[pos_start:pos_end+1], start, end, strand))
            else:
                self.genes.append(Gene(gene, name, desc, numpy.array([[]]), numpy.array([]), start, end, strand))
            self.orf2index[gene] = count
            count += 1


    def local_insertions(self):
        """Returns numpy array with the number of insertions, 'k', for each gene."""
        G = len(self.genes)
        K = numpy.zeros(G)
        for i in xrange(G):
            K[i] = self.genes[i].k
        return K


    def local_sites(self):
        """Returns numpy array with total number of TA sites, 'n', for each gene."""
        G = len(self.genes)
        N = numpy.zeros(G)
        for i in range(G):
            N[i] = self.genes[i].n
        return N


    def local_runs(self):
        """Returns numpy array with total number of TA sites, 'r', for each gene."""
        G = len(self.genes)
        R = numpy.zeros(G)
        for i in xrange(G):
            R[i] = self.genes[i].r
        return R 


    def local_reads(self):
        """Returns numpy array of lists containing the read counts for each gene."""
        all_reads = []
        G = len(self.genes)
        for i in xrange(G):
            all_reads.extend(self.genes[i].reads)
        return numpy.array(all_reads)


    def local_thetas(self):
        """Returns numpy array of insertion frequencies, 'theta', for each gene."""
        G = len(self.genes)
        theta = numpy.zeros(G)
        for i in xrange(G):
            theta[i] = self.genes[i].theta()
        return theta


    def local_phis(self):
        """Returns numpy array of non-insertion frequency, 'phi', for each gene."""
        return 1.0 - self.theta()


    ######

    def global_insertion(self):
        """Returns total number of insertions, i.e. sum of 'k' over all genes."""
        G = len(self.genes)
        total = 0
        for i in xrange(G):
            total += self.genes[i].k
        return total

    def global_sites(self):
        """Returns total number of sites, i.e. sum of 'n' over all genes."""
        G = len(self.genes)
        total = 0
        for i in xrange(G):
            total += self.genes[i].n
        return total

    def global_run(self):
        """Returns the run assuming all genes were concatenated together."""
        return maxrun(self.tosses())

    def global_reads(self):
        """Returns the reads among the library."""
        return self.data

    def global_theta(self):
        """Returns global insertion frequency, of the library."""
        return float(self.global_insertion())/self.global_sites()

    def global_phi(self):
        """Returns global non-insertion frequency, of the library."""
        return 1.0 - self.global_theta()

    def total_reads(self):
        """Returns total reads among the library."""
        reads_total = 0
        for g in self.genes:
            reads_total += g.total_reads()
        return reads_total

    def tosses(self):
        """Returns list of bernoulli trials, 'tosses', representing insertions in the gene."""
        all_tosses = []
        for g in self.genes:
            all_tosses.extend(g.tosses)
        return all_tosses




def tossify(data):
    """Reduces the data into Bernoulli trials (or 'tosses') based on whether counts were observed or not."""
    K,N = data.shape
    reduced = numpy.sum(data,0)
    return numpy.zeros(N) + (numpy.sum(data, 0) > 0)


def runs(data):
    """Return list of all the runs of zero, given a list of counts."""
    if len(data) == 0: return [0]
    runs = []
    current_r = 0
    for read in data:
        if read > 0: # If ending a run of zeros
            #if current_r > 0: # If we were in a run, add to list
            runs.append(current_r)
            current_r = 0
        else:
            current_r += 1
    # If we ended in a run, add it
    if current_r > 0:
        runs.append(current_r)
    return runs

def runindex(runs):
    """Returns a list of the indexes of the start of the runs; complements runs()."""
    index = 0
    index_list = []
    runindex = 0
    for r in runs:
        for i in range(r):
            if i == 0:
                runindex = index
            index+=1
        if r == 0:
            runindex = index
            index+=1
        index_list.append(runindex)
    return index_list


def get_data(wig_list):
    """ Returns a tuple of (data, position) containing a matrix of raw read counts, and list of coordinates. """
    K = len(wig_list)
    T = 0
    for line in open(wig_list[0]):
        if line.startswith("#"): continue
        if line.startswith("location"): continue
        if line.startswith("variable"): continue
        T+=1

    data = numpy.zeros((K,T))
    position = numpy.zeros(T)
    for j,path in enumerate(wig_list):
        reads = []
        i = 0
        for line in open(path):
            if line.startswith("#"): continue
            if line.startswith("location"): continue
            if line.startswith("variable"): continue
            tmp = line.split()
            pos = int(tmp[0])
            rd = float(tmp[1])
            data[j,i] = rd
            position[i] = pos
            i+=1
    return (data, position)


def get_pos_hash(path):
    """Returns a dictionary that maps coordinates to a list of genes that occur at that coordinate."""
    hash = {}
    for line in open(path):
        if line.startswith("#"): continue
        tmp = line.strip().split("\t")
        orf = tmp[8]
        start = int(tmp[1])
        end = int(tmp[2])
        for pos in range(start, end+1):
            if pos not in hash: hash[pos] = []
            hash[pos].append(orf)
    return hash



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


def maxrun(lst,item=0):
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





if __name__ == "__main__":

    G = Genes(sys.argv[1].split(","), sys.argv[2])
    print "#Insertion: %s" % G.global_insertion()
    print "#Sites: %s" % G.global_sites()
    print "#Run: %s" % G.global_run()
    print "#Theta: %1.4f" % G.global_theta()
    print "#Phi: %1.4f" % G.global_phi()
    
    orf = "Rv0016c"
    print G[orf]
    print G[orf].reads
    print G[orf].position
    print G[orf].t


    for gene in G:
        k = gene.k
        n = gene.n
        r = gene.r
        s = gene.s
        t = gene.t
        pos = gene.position
        #print "%s" % gene
        print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (gene.orf, gene.desc, k, n, r, s, t)



