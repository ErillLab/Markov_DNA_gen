import random
import numpy as np

class MCM:

    def __init__(self, k=1):
        """
        Markov Chain Model constructor.

        Args:
            k (int): Order of the Markov model.
        """
        # Markov model order
        self.k = k

        # Transition matrix
        # key: Represents the state, a sequence of n nucleotides
        # value: Transition probabilities for the state in form of dict.
        #         Each entry represents the probability to make a transition
        #         to a nucleotide, being the nucleotide the key
        self.transition = {}

        # Transition matrix of auxiliar Markov model
        # Similar structure as self.transition
        self.aux_transition = [{} for _ in range(self.k)]

    def train(self, seq):
        """
        Trains the model using a DNA sequence, computing the transition matrix of the model.

        An auxiliar Markov model is trained, with an order in which there is no state without 
        transitions. In case of arriving to some state withput transitions, this auxiliar model
        is used to generate the next transition.

        Args:
            seq (string): DNA sequence to train the model.
        """
        self.size = len(seq)
        seq = seq.upper()
        self.seq = seq
        self.transition = {}
        for i in range(0, len(seq)-self.k):
            kmer = seq[i: i+self.k]
            if not kmer in self.transition:
                self.transition[kmer] = {}
                self.transition[kmer]["A"] = 0
                self.transition[kmer]["C"] = 0
                self.transition[kmer]["G"] = 0
                self.transition[kmer]["T"] = 0
                self.transition[kmer]["tot"] = 0
            self.transition[kmer][seq[i+self.k]] += 1
            self.transition[kmer]["tot"] += 1

            ## Computes all posible models
            for k in range(0, self.k):
                if k == 0:
                    aux_kmer = ""
                else:
                    aux_kmer = kmer[-k:]
                if not aux_kmer in self.aux_transition[k]:
                    self.aux_transition[k][aux_kmer] = {}
                    self.aux_transition[k][aux_kmer]["A"] = 0
                    self.aux_transition[k][aux_kmer]["C"] = 0
                    self.aux_transition[k][aux_kmer]["G"] = 0
                    self.aux_transition[k][aux_kmer]["T"] = 0
                    self.aux_transition[k][aux_kmer]["tot"] = 0
                self.aux_transition[k][aux_kmer][seq[i+self.k]] += 1
                self.aux_transition[k][aux_kmer]["tot"] += 1

        for key in self.transition.keys():
            self.transition[key]["A"] /= self.transition[key]["tot"]
            self.transition[key]["C"] /= self.transition[key]["tot"]
            self.transition[key]["G"] /= self.transition[key]["tot"]
            self.transition[key]["T"] /= self.transition[key]["tot"]

        ## Searches for the model with the highest order that has no states without transitions
        for k in range(self.k-1, -1, -1):
            if len(self.aux_transition[k]) == 4**(k):
                self.aux_transition = self.aux_transition[k]
                for key in self.aux_transition.keys():
                    self.aux_transition[key]["A"] /= self.aux_transition[key]["tot"]
                    self.aux_transition[key]["C"] /= self.aux_transition[key]["tot"]
                    self.aux_transition[key]["G"] /= self.aux_transition[key]["tot"]
                    self.aux_transition[key]["T"] /= self.aux_transition[key]["tot"]
                self.aux_k = k
                break

    def sample(self, state):
        """
        Samples a new nucleotide given an array of probabilities.

        Args:
            state: Dictionary of probabilities for nucleotides.
        """
        if state in self.transition.keys():
            probs = self.transition[state]
        else:
            probs = self.aux_transition[state[-self.aux_k:]]
        return np.random.choice(["A","C","G","T"], p=[probs["A"], probs["C"], probs["G"], probs["T"]])

    def generate(self, size=None, N=1):
        """
        Generates a new sequence of nucleotides using the trained model.

        Args:
            size (int): Length of the sequence to be generated.

            N (int): Number of sequences to be generated.
        """
        if size == None:
            size = self.size

        res = []
        for _ in range(N):
            if self.k == 0:
                probs = self.transition[""]
                seq = np.random.choice(["A","C","G","T"], p=[probs["A"], probs["C"], probs["G"], probs["T"]], size=n)
                seq = "".join(seq)
            else:
                seq = self.initial_state()
                state = seq
                for _ in range(size-len(state)):
                    nuc = self.sample(state)
                    seq += nuc
                    state = seq[-self.k:]
            res.append(seq)
        return res
    
    def generator(self, size=None, N=None):
        """
        Generator for sampling new sequences of nucleotides using the trained model.

        Args:
            size (int):
                Length of each generate sequence. If None, the length of the original
                sequence used to train the model will be used.

            N (int):
                Number of sequences to be generated. If None, the generator is infinite.
        """
        if size == None:
            size = self.size

        if N == None:
            while True:
                seq = self.generate(size=size)[0]
                yield seq
        else:
            for _ in range(N):
                seq = self.generate(size=size)[0]
                yield seq

    def initial_state(self):
        """
        Computes the initial k-mer of the sequence sampling from the list
        of appeared k-mers in the training sequence.
        """
        keys = list(self.transition.keys())
        return random.choice(keys)
