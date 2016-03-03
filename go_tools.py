import sys
import math
import random

GOPATH = "/pacific/HomeFrozen/michael.dejesus/FUNCTIONS/gene_ontology_edit.obo"
RV2GOPATH = "/pacific/HomeFrozen/michael.dejesus/FUNCTIONS/rv2go.dat"


class GoTerm:
    """Class representing a gene-ontology term.

    Attributes:
        id (str): A string representing the GO id. eg. "GO:000000".
        name (str): A string representing the 'name' or small description.
        desc (str): A string representing the 'description' of the GO term.
        parents (list): A list containing the ids of its parents.
        height (integer): A integer representing it's height in the graph.
    
    """
    def __init__(self, id, name="", desc="", parents=[], children=[]):
        """Initializes a GO Term with the provided attribute values."""
        self.id = id
        self.name = name
        self.desc = desc
        self.parents = parents
        self.children = children
        self.height = float('inf')

    def __str__(self):
        """Prints the id and name of the GO term so it looks nice."""
        return( "%s - %s" % (self.id, self.name))



class GoTerms:
    """Reads in a gene-ontology .obo file and creates a list of GoTerm objects.

    Faciliates handling GO terms, by reading in .obo files and saving their
    name, description. etc. Contains methods for finding their ancestors
    and descendants.

    Attributes:
        terms: A list of GoTerm objects for each go term.
        children: Dictionary going from GO id to its children.
        parents: Dictionary going from GO id to its parents.
        root: String defining the ID of the top level GO term.

    """

    def __getitem__(self, i):
        """Defines the getitem function to work as both a dictionary
        and/or a list. If a string is passed, looks for the GO term
        matching that id. (Not very efficient). If an integer is passed
        returns the id in that position.
        """
        if isinstance(i, int):
            return self.terms[i]

        if isinstance(i, basestring):
            j = self.term2index[i]
            return self.terms[j]

    def __contains__(self, item):
        """Defines the contain function to check if an id exists."""
        return item in self.term2index

    def __len__(self):
        """Defines the length function, returning the number of terms."""
        return( len(self.terms) )

    def __str__(self):
        """Defines the string function to make it print nice."""
        return("GoTerm Object (N=%d)" % len(self.terms))

    def __init__(self, path):
        """Initializes the GoTerms object given the required .obo file."""
        self.terms = []
        self.root = "GO:0000000"
        self.term2index = {self.root:0}
        self.children = {self.root:["GO:0008150", "GO:0005575", "GO:0003674"]}
        self.parents = {}
        id = name =  desc = ""
        parents = []
        file = open(path, "r")
        line = file.readline()
        while line:
            if line.startswith("\n"):
                if id != "":
                    self.term2index[id] = len(self.terms)
                    self.terms.append(GoTerm(id, name, desc, parents))
                    for p in parents:
                        if p not in self.children: self.children[p] = []
                        self.children[p].append(id)
                    self.parents[id] = parents
                id = name = desc = ""
                parents = []
            if line.startswith("id:"):
                id = line.split()[1].strip()
                if not id.startswith("GO"): id = ""
            if line.startswith("name:"): name = line.split(" ", 1)[1].strip()
            if line.startswith("def:"): desc = line.split(" ", 1)[1].strip()
            if line.startswith("is_a:"): parents.append( line.split(" ")[1].strip() )
            if line.startswith("relationship: part_of "): parents.append( line.split(" ")[2].strip() )

            line = file.readline()
        file.close()
        self.update_height()
        self.update_children()


    def full_children(self, id):
        """Returns a list of ALL the children of a specific go term.

        Args:
            id (str): String with the id you want children from.

        Returns:
            list. Includes all descdendants.
        """
        full_list = []
        fringe = [c for c in self.children.get(id, [])]
        while fringe:
            term = fringe.pop()
            full_list.append(term)
            for c in self.children.get(term,[]):
                fringe.append(c)

        return(full_list)
        
    def full_parents(self, id):
        """Returns a list of aLL the ancestors of a specific go term.

        Args:
            id (str): String with the id you want parents from.

        Returns:
            list. Includes all ancestors.

        """
        full_list = []
        fringe = [c for c in self.parents.get(id,[])]
        while fringe:
            term = fringe.pop()
            full_list.append(term)
            for c in self.parents.get(term,[]):
                fringe.append(c)
        return(full_list)

    def nth_ancestor(self, id, n):
        """Returns list of the ancestors up to a specific  height.

        Args:
            id (str): String of the id you want ancestors from.
            n (int): Integer of the max height height you want ancestors from.

        Returns:
            list. 
        """
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
        """Returns list of the descendants up to a specific height.

        Args:
            id (str): String of the id you want descendants from.
            n (int): Integer of the max height height you want descendants from.
        
        Returns:
            list.
        """
        n_children = {1:self.children.get(id,[])}
        k = 1
        while k < n:
            k_children = []
            for p in n_children[k]:
                k_children.extend(self.children.get(p,[]))
            k+=1
            n_children[k] = k_children

        return(list(set(n_children[n])))

    def update_height(self):
        """Updates the height of all the terms."""
        h = 0
        children = [self.root]
        self[self.root].height = h
        while not id in children and len(children) > 0:
            next = []
            for c in children:
                self[c].height = min(h, self[c].height)
                next.extend(self.children.get(c,[]))
            children = next
            h+=1

    def update_children(self):
        """Updates the children of all the terms."""
        children = [self.root]
        self[self.root].children = []
        while not id in children and len(children) > 0:
            next = []
            for c in children:
                self[c].children = self.children.get(c, [])
                next.extend(self.children.get(c,[]))
            children = next            


    def get_height(self, id):
        """Returns the height of the term with the given id.

        Args:
            id (str): String representing the id of the go term.
        
        Returns:
            int. 
        """
        return self[id].height



if __name__ == "__main__":


    # Example on how to use the GoTerms Class
    # GOPATH is a path to the desired .obo file

    print "Using .obo file: %s" % GOPATH
    G = GoTerms(GOPATH)
    N = len(G)
    print ""
    print G
    print "Size: %d" % N
    print "First Term: %s" % G[0]
    print "Last Term: %s" % G[N-1]

    ii = random.randint(0,N+1)
    ID = G[ii].id
    print "Random Term: %s" % G[ii]
    print "Height:     %s" % G.get_height(ID)
    print "Children: %s" % ",".join(G.full_children(ID))




