import math
import sys


#GOPATH = "/pacific/home/michael.dejesus/gene_ontology.1_2.obo"
GOPATH = "/pacific/home/michael.dejesus/FUNCTIONS/gene_ontology_edit.obo"

RV2GOPATH = "/pacific/home/michael.dejesus/FUNCTIONS/rv2go.dat"


class GoTerm:
    def __init__(self, id, name="", desc="", parents=[], genes=[]):
        self.id = id
        self.name = name
        self.desc = desc
        self.parents = parents
        self.genes = genes
        self.height = -1

    def __str__(self):
        return( "%s (%s)" % (self.id, self.name))



class GoTerms:

    def __getitem__(self, i):
        
        if isinstance(i, int):
            return(self.terms[i])
        if isinstance(i, basestring):
            for j in range(len(self.terms)):
                if self.terms[j].id == i:
                    return(self.terms[j])

    def __contains__(self, item):
        for term in self.terms:
            if term.id == item: return(True)
        return(False)

    def __len__(self):
        return( len(self.terms) )

    def __str__(self):
        return("GoTerm Object (N=%d)" % len(self.terms))

    def add(self, term):
        self.terms.append(term)

    #@classmethod
    #def fromtrashfile(self, path, M=MINIMUM_READ, L=LANES):
    def __init__(self, P):
        self.terms = []
        self.children = {"GO:0000000":["GO:0008150", "GO:0005575", "GO:0003674"]}
        self.parents = {}
        self.root = "GO:0000000"
        path = P or []
        id = ""; name = ""; desc = ""; parents = [];
        #print min_read, lanes, term_list, reads, tosses, self.terms
        file = open(path, "r")
        line = file.readline()
        while line:
            if line.startswith("\n"):
                if id != "":
                    self.terms.append(GoTerm(id, name, desc, parents, []))
                    for p in parents:
                        if p not in self.children: self.children[p] = []
                        self.children[p].append(id)
                    self.parents[id] = parents
                id = ""; name = ""; desc = ""; parents = [];
            if line.startswith("id:"): id = line.split()[1].strip()
            #print id
            if line.startswith("name:"): name = line.split(" ", 1)[1].strip()
            if line.startswith("def:"): desc = line.split(" ", 1)[1].strip()
            if line.startswith("is_a:"): parents.append( line.split(" ")[1].strip() )
            if line.startswith("relationship: part_of "): parents.append( line.split(" ")[2].strip() )

            line = file.readline()



    def full_children(self, id):
    
        full_list = []
        fringe = [c for c in self.children.get(id, [])]
        while fringe:
            term = fringe.pop()
            full_list.append(term)
            for c in self.children.get(term,[]):
                fringe.append(c)

        return(full_list)
        
    def full_parents(self, id):

        full_list = []
        fringe = [c for c in self.parents.get(id,[])]
        while fringe:
            term = fringe.pop()
            full_list.append(term)
            for c in self.parents.get(term,[]):
                fringe.append(c)
        return(full_list)

    def nth_ancestor(self, id, n):
        n_parents = {1:self.parents.get(id,[])}
        k = 1
        while k < n:
            k_parents = []
            for p in n_parents[k]:
                k_parents.extend(self.parents.get(p,[]))
            k+=1
            n_parents[k] = k_parents
    
        return(list(set(n_parents[n])))

    def nth_descendant(self,id,n):
        n_children = {1:self.children.get(id,[])}
        k = 1
        while k < n:
            k_children = []
            for p in n_children[k]:
                k_children.extend(self.children.get(p,[]))
            k+=1
            n_children[k] = k_children

        return(list(set(n_children[n])))

    def get_height(self, id):
        h = 0
        children = [self.root]
        while not id in children and len(children) > 0:
            next = []
            for c in children:
                next.extend(self.children.get(c,[]))
            children = next
            h+=1
        return(h)
        
