import random
import numpy as np

class MCM:

    NUCLEOTIDES = ["A", "C", "G", "T"]
    NONE = 0
    CIRCULAR = 1
    AUXILIAR = 2
    RANDOM = 3

    def __init__(self, n=-1, mode=0):
        """
        Markov Chain Model constructor.

        Args:
            n (int): Order of the Markov model.

            mode (int): Can be 0, 1(CIRCULAR), 2(AUXILIAR), 3(RANDOM). 
                This mode is use in the generation of sequences, in case 
                of arriving to some state without transitions.

                - CIRCULAR: Duplicates the training sequence.
                - AUXILIAR: An auxiliar Markov model is trained, with an 
                  order in which there are no any states without transitions. 
                  In case of arriving to some state withput transitions, this 
                  auxiliar model is used to generate the next transition.
                - RANDOM: The transition probabilities are generated randomly.
        """
        # Markov model order
        self.n = n

        # Transition matrix
        self.transition = {}

        # Auxiliar mode
        self.mode = mode
        if self.mode == self.AUXILIAR:
            # Transition matrix of auxiliar Markov model
            self.aux_transition = [{} for _ in range(n-1)]

    def train(self, seq):
        """
        Trains the model using a DNA sequence, computing the transition matrix of the model.

        Args:
            seq (string): DNA sequence to train the model.
        """
        seq = seq.upper()
        if self.mode == self.CIRCULAR: 
            # Duplicates the sequence in order to avoid states without transitions
            seq += seq
        self.seq = seq
        self.transition = {}
        for i in range(0, len(seq)-self.n):
            kmer = seq[i: i+self.n]
            if not kmer in self.transition:
                self.transition[kmer] = {}
                self.transition[kmer]["A"] = 0
                self.transition[kmer]["C"] = 0
                self.transition[kmer]["G"] = 0
                self.transition[kmer]["T"] = 0
                self.transition[kmer]["tot"] = 0
            self.transition[kmer][seq[i+self.n]] += 1
            self.transition[kmer]["tot"] += 1
            if self.mode == self.AUXILIAR: 
                ## Computes all posible models
                for k in range(0, self.n-1):
                    aux_kmer = kmer[-k-1:]
                    if not kmer in self.aux_transition[k]:
                        self.aux_transition[k][aux_kmer] = {}
                        self.aux_transition[k][aux_kmer]["A"] = 0
                        self.aux_transition[k][aux_kmer]["C"] = 0
                        self.aux_transition[k][aux_kmer]["G"] = 0
                        self.aux_transition[k][aux_kmer]["T"] = 0
                        self.aux_transition[k][aux_kmer]["tot"] = 0
                    self.aux_transition[k][aux_kmer][seq[i+self.n]] += 1
                    self.aux_transition[k][aux_kmer]["tot"] += 1
        for key in self.transition.keys():
            self.transition[key]["A"] /= self.transition[key]["tot"]
            self.transition[key]["C"] /= self.transition[key]["tot"]
            self.transition[key]["G"] /= self.transition[key]["tot"]
            self.transition[key]["T"] /= self.transition[key]["tot"]

        if self.mode == self.AUXILIAR:
            ## Computes the transition matrix for the auxiliar Markov model
            for k in range(self.n-1-1, -1, -1):
                ## Searches for the model with the highest order that has no states without transitions
                print(len(self.aux_transition[k]), 4**(k+1), k+1)
                if len(self.aux_transition[k]) == 4**(k+1):
                    # print("Auxiliar transition matrix of order", k+1)
                    self.aux_transition = self.aux_transition[k]
                    for key in self.aux_transition.keys():
                        self.aux_transition[key]["A"] /= self.aux_transition[key]["tot"]
                        self.aux_transition[key]["C"] /= self.aux_transition[key]["tot"]
                        self.aux_transition[key]["G"] /= self.aux_transition[key]["tot"]
                        self.aux_transition[key]["T"] /= self.aux_transition[key]["tot"]
                    self.aux_k = k+1
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
            probs = self.auxiliar_transition(state)
        res = np.random.multinomial(1, [probs["A"], probs["C"], probs["G"], probs["T"]])
        if res[0] == 1:
            return "A"
        if res[1] == 1:
            return "C"
        if res[2] == 1:
            return "G"
        return "T"

    def generate(self, n, N=1):
        """
        Generates a new sequence of nucleotides using the trained model.

        Args:
            n (int): Length of the sequence to be generated.

            N (int): Number of sequences to be generated.
        """
        res = []
        for _ in range(N):
            seq = self.initial_state()
            state = seq
            for _ in range(n):
                nuc = self.sample(state)
                seq += nuc
                state = seq[-self.n:]
            res.append(seq)
        return res
    
    def generator(self, n, N=1):
        """
        Generator for sampling new sequences of nucleotides using the trained model.

        Args:
            n (int): Length of the sequence to be generated.

            N (int): Number of sequences to be generated.
        """
        for _ in range(N):
            seq = self.initial_state()
            state = seq
            for _ in range(n):
                nuc = self.sample(state)
                seq += nuc
                state = seq[-self.n:]
            yield seq

    def initial_state(self):
        """
        Computes the initial k-mer of the sequence sampling from the list
        of appeared k-mers in the training sequence.
        """
        keys = list(self.transition.keys())
        return random.choice(keys)
    
    def auxiliar_transition(self, state):
        """
        Generates auxiliar transition probabilities with the specified
        method by the self.mode parameter.

        Args:
            state (string): The current state to in case of using the 
            auxiliar Markov model.
        """
        if self.mode == self.AUXILIAR:
            # Uses the transition probabilities of the auxiliar model trained.
            # The state used for the auxiliar model has less nucleotides, since
            # the trained model has a lower order -> state[-self.aux_k:]
            probs = self.aux_transition[state[-self.aux_k:]]
        elif self.mode == self.RANDOM:
            # Generates auxiliar transition probabilities randomly
            a = random.randint()
            c = random.randint()
            g = random.randint()
            t = random.randint()
            tot = a + c + g + t
            probs = {"A:": a/tot, "C": c/tot, "G": g/tot, "T": t/tot}
        else:
            print("ERROR: Undefined transitions. (%s)"%state)
            exit(-1)
        return probs
