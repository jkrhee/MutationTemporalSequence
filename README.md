# MutationTemporalSequence


## Prerequsition

To run this program, 'numpy' package is required. If this package is not on your system, please install it using `pip` or `easy-install`.


### Input data

This program requires 2 input data, gene_names.txt and input data for analysis.

The columns of the input data are composed of gene name, chromsome #, coordinate, ref, alt, sample id, CCF distribution [0.01, 1].


## Running this tool

Please type:

```
  python mutTempSeq.py [input_data] [output]
```

Eg.
```
  python mutTempSeq.py example_data result.txt
```

