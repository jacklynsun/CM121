import copy
import string
#from collections import defaultdict
from itertools import product
import random
#import graphviz
#from graphviz import Digraph
import numpy as np

#generate random genome of size 1000
def randomGenome():
    randGenome = "".join(random.choice("ACGT") for _ in range(1000))
    return randGenome

#print out a pdf picture of graph for visuals, only input is X for demonstration
class printGraph:
    #constructor
    def __init__(self):
        self.place = 0
        self.preNodes = []
        self.postNodes = []

    #make Kmers: using genome read, make k-mers of size X and return array of k-mers
    def makeKmers(self, read, X):
        ans = []
        for i in range(len(read) - X + 1):
            ans.append(read[i:i + X].upper())

        return ans

    #makes Nodes: using kmers as input, splice kmers into prefix and suffix nodes
    def makeNodes(self, kmers):
        nodes = {}
        for kmer in kmers: #iterate through kmers
            prefix = kmer[:-1] #array splice for prefix
            self.preNodes.append(prefix)
            suffix = kmer[1:] #array splice for suffix
            self.postNodes.append(suffix)
            if prefix in nodes:
                nodes[prefix].append(suffix)
            else:
                nodes[prefix] = []
                nodes[prefix].append(suffix)
        return nodes

    #using graphviz to create graph and display for visuals
    def makeGraph(self, nodes):
        edges = []
        dot = Digraph(format='pdf')
        for k_1mer in nodes:
            dot.node(k_1mer, k_1mer)
        for head in nodes:
            tails = nodes[head]
            for tail in tails:
                edge = head[0] + tail  # edge label
                edges.append(edge)
                dot.edge(head, tail, edge)

        dot.view()
        return edges

    def generateDeBruijn(self, X):
        #read = randomGenome()
        #print(read)
        read = "TAATGTACCATGGTAGATGTT" #example genome read, you can input the randomGenome as well, graph just gets very messy
        kmers = self.makeKmers(read, X)
        nodes = self.makeNodes(kmers)
        return self.makeGraph(nodes)

#temp = printGraph()
#temp.generateDeBruijn(3) #k-mer length 3 as example


#Generating a DB graph with X, Y, Z as input
class DeBruijnGraph:

    #create k-mers from read
    @staticmethod
    def chop(read, X):
        for i in range(len(read) - (X - 1)):
            yield (read[i:i + X], read[i:i + X - 1], read[i + 1:i + X])

    class Node:
        #constructor
        def __init__(self, km1mer):
            self.km1mer = km1mer
            self.nin = 0
            self.nout = 0

        def isSemiBalanced(self):
            return abs(self.nin - self.nout) == 1

        def isBalanced(self):
            return self.nin == self.nout

        def __hash__(self):
            return hash(self.km1mer)

        def __str__(self):
            return self.km1mer

    #constructor to build db graph
    def __init__(self, strIter, k, circularize=False):
        #self.walkCount = 0
        self.G = {}  #map nodes to its adjacent nodes
        self.nodes = {}  #maps k-1-mers to Node objects
        for st in strIter:
            if circularize:
                st += st[:k - 1]
            for kmer, km1L, km1R in self.chop(st, k):
                nodeL, nodeR = None, None
                if km1L in self.nodes:
                    nodeL = self.nodes[km1L]
                else:
                    nodeL = self.nodes[km1L] = self.Node(km1L)
                if km1R in self.nodes:
                    nodeR = self.nodes[km1R]
                else:
                    nodeR = self.nodes[km1R] = self.Node(km1R)
                nodeL.nout += 1
                nodeR.nin += 1
                self.G.setdefault(nodeL, []).append(nodeR)

        # Iterate over nodes; tally # balanced, semi-balanced, neither
        self.nsemi, self.nbal, self.nneither = 0, 0, 0
        # Keep track of head and tail nodes in the case of a graph with
        # Eularian walk (not cycle)
        self.head, self.tail = None, None
        for node in iter(self.nodes.values()):
            if node.isBalanced():
                self.nbal += 1
            elif node.isSemiBalanced():
                if node.nin == node.nout + 1:
                    self.tail = node
                if node.nin == node.nout - 1:
                    self.head = node
                self.nsemi += 1
            else:
                self.nneither += 1

    #def walkCounter(self):
    #    self.walkCount = self.walkCount + 1
    #    return (self.walkCount)

    #returns number of nodes
    def nnodes(self):
        return len(self.nodes)

    #returns number of edges
    def nedges(self):
        return len(self.G)

    #returns true if eulerian walk exists
    def hasEulerianWalk(self):
        return self.nneither == 0 and self.nsemi == 2

    #returns true if eulerian cycle exists
    def hasEulerianCycle(self):
        return self.nneither == 0 and self.nsemi == 0

    #returns true if eulerian cycle or walk exists
    def isEulerian(self):
        return self.hasEulerianWalk() or self.hasEulerianCycle()

    #traverse and return nodes travelled in eulerian walk/cycle and print out
    def eulerianWalkOrCycle(self):
        assert self.isEulerian()
        g = self.G
        if self.hasEulerianWalk():
            g = g.copy()
            g.setdefault(self.tail, []).append(self.head)
        # graph g has an Eulerian cycle
        tour = []
        src = next(iter(g.keys()))  #pick any starting node

        def __visit(n):
            while len(g[n]) > 0:

                listSize = len(g[n])
                x = random.randint(0,listSize - 1)
                dst = g[n].pop(x) #remember it's pop(x) for random node in list to tour
                __visit(dst)
            tour.append(n)

        __visit(src)
        tour = tour[::-1][:-1]  #reverse, take all but last node to get eulerian path

        if self.hasEulerianWalk(): #if we find a eulerian walk
            #adjust list of nodes such that starts at head and ends at tail
            sti = tour.index(self.head)
            tour = tour[sti:] + tour[:sti]
        elif not self.hasEulerianWalk(): #else if we do not find a eulerian walk, look again
            #self.walkCounter(self)
            self.eulerianWalkOrCycle()
        # Return node list
        return list(map(str, tour))

#create repeats of length Y randomly and randomly place them in different locations of the genome
def repeats(read, Y, Z):
    subRead = []
    subRead = read[:(Y)]  #take the last y items of read
    for i in range(Z): #Z number of repeats
      insertLoc = random.randint(0, len(read) - 1) #random insertion into genome
      read = read[0:insertLoc] + subRead + read[insertLoc:len(read)] #insert repeat into genome

    return read

#user input
X = 3
Y = 2
Z = 2
genome = "TAATGCCATGGGATGTT"

#running the function
inputString = repeats(genome,Y,Z)
print(inputString)
g = DeBruijnGraph(inputString, X)
walk = DeBruijnGraph([inputString], X).eulerianWalkOrCycle()
walk[0] + ''.join(map(lambda x: x[-1], walk[1:]))
print(walk) #outputs node traversal/ euler path traversal
