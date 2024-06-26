# Markov_DNA_gen

A Markov Model DNA sequence generator to generate pseudo-replicate sequences based on an input sequence.

## Installation

``` bash
pip3 install .
```

## User guide

Import module:

```python
    from Markov_DNA import MCM
```

Train model:

```python
    mcm = MCM(k, mode=mode)
    mcm.train(seq)
```

Where:
  - k is the model order
  - seq is the reference sequence to train the model
  - mode is the selected mode that determines how the transition probabilities are computed in case of arriving to null transition states.
    - NONE (0): Normal mode, no additional transition probabilities are computed, in case of encountering null transition states the program will fail.
    - CIRCULAR (1): Duplicates the input sequence so that there are no null transitions.
    - AUXILIAR (2): Trains an alternative model with a lower order in which are present all possible states
    - RANDOM (3): The auxiliar transition probabilities are computed randomly.

Generate sequences:

```python
    seqs = mcm.generate(l, N=N)
```

Where:
  - l is the length of the sequence to be generated.
  - N is the number of sequences to be generated.

## Example

The example.py file is an example that generates 10 sequences from a base sequence. It can be executed as follows:

```bash
    python3 example.py
```

## Authors

Erik Blázquez Fernández (erikblazfer@outlook.es)