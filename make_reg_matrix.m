%%% Builds regularization matrices
%%% David Newman (newman@uci.edu)
%%% Edwin V. Bonilla (edwin.bonilla@nicta.com.au)
%%% Last udpdate: 15/03/2012

clear all; 

vocab_fname = 'vocab.txt'; 


%% QuadReg
regmatrix_fname  = 'w1w2pmi.txt';
make_Sdiag1(regmatrix_fname, vocab_fname , 'Sdiag1.mat');

%% ConvReg
regmatrix_fname  = 'w1w2cooc.txt';
make_Scon1(regmatrix_fname, vocab_fname, 'Scon1.mat');




