rtm
===

Regularized Topic Models

Implements the regularized topic models in Newman, Bonilla and Buntine [1]

Authors: 
-----
David Newman (newman@uci.edu)
Edwin V. Bonilla (edwinbonilla@gmail.com)


Requirements:
-----------
This is a self-contained package. You will need to compile the necessary
mex files: 
>> mex gibbsmex_lda.c
>> mex gibbsmex_semi.c


Main Contents
--------------
1. regularized_lda.m    : Learns an RTM model. 
2. make_reg_matrix.m    : Builds regularization matrices 
3. run_regularized_lda.m: An example of how to run the algorithms
 
Additionally, we provide in a separate file rtm_data.zip the "climate" 
dataset used in [1] as a test example:
a) Ndw.txt: File containing the corpus data
(b) vocab.txt: File containing the vocabulary

References
-----------
[1] David Newman and Edwin V. Bonilla and Wray Buntine.
Improving Topic Coherence with Regularized Topic Models.
Advances in Neural Information Processing Systems 24: NIPS'2011

