# Markov_DNA_gen

A Markov Model DNA sequence generator to generate pseudo-replicate sequences based on an input sequence.

## Requirements

- NumPy

## Installation

It can be installed using pip:

``` bash
pip3 install .
```

Or using conda:

``` bash
conda install erikblazfer::markov_dna
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

Can be done in two different ways, using a generator in an iterative way:

```python
    for new_seq in mcm.generator(l, N=N):
        print(new_seq)
```

Or by calling a function that generates a list of sequences:

```python
    seqs = mcm.generate(l, N=N)
```

Where:
  - l is the length of the sequence to be generated.
  - N is the number of sequences to be generated.

The advantage of the first method is that you do not need to keep all the sequences in memory, while the second one allows you to obtain a list of sequences directly.

## Example

There are some examples located in the examples folder. To execute them, the following command can be used:

```bash
    python3 examples/exampleX.py
```

Where X can be [1..4] depending on the genereting method.

- Example 1: Shows how to use the function that returns the list of sequences.
- Example 2: Shows how to use the iterative form of the generator generator.
- Example 3: An example using the ```AUXILIAR``` mode for empty transitions.
- Example 4: An example using the ```RANDOM``` mode for empty transitions.

## Authors

Erik Blázquez Fernández (erikblazfer@outlook.es)