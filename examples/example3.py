from Markov_DNA import MCM

if __name__ == "__main__":
    seq = "TGAGGACTTTaggatAGGATTTTGTCATCATCAAAAGACATTCTTGTAAATTATATGAAACGTCGGT"
    
    mcm = MCM(4, MCM.AUXILIAR)

    mcm.train(seq)

    for new_seq in mcm.generator(len(seq), 10):
        print(new_seq)