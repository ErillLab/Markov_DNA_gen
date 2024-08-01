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
    mcm = MCM(k)
    mcm.train(seq)
```

where:
  - k is the model order
  - seq is the reference sequence to train the model

Generate sequences:

Can be done in two different ways, using a generator in an iterative way:

```python
    for new_seq in mcm.generator(size=l, N=N):
        print(new_seq)
```

Or by calling a function that generates a list of sequences:

```python
    seqs = mcm.generate(size=l, N=N)
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

Where X can be [1..2] depending on the genereting method.

- Example 1: Shows how to use the function that returns the list of sequences.
- Example 2: Shows how to use the iterative form of the generator generator.

## Authors

Erik Blázquez Fernández (erikblazfer@outlook.es)