%%% Regularized LDA
%%% Newman, Bonilla, Buntine
%%% NIPS 2011
%%% David Newman (newman@uci.edu)
%%% Edwin V. Bonilla (edwin.bonilla@nicta.com.au)
%%% Last udpdate: 15/03/2012
% This runs an example of regularized LDA

clear all;

%% Data files
corpus_fname = 'Ndw.txt';
vocab_fname  = 'vocab.txt';


%% Parameter setting
config.T            = 10;
config.beta         = 0.01;
config.alpha        = -999;  %% later set alpha = 0.05 * N / (D*T)
config.gibbs_iter   = 500;
config.lag_iter     = 50;
config.reg_iter     = 10;
config.nu           = 0.5;


%% Run baselie LDA
reg = 0; regmatrix_fname = [];
rand('state', 7); 
[Nwt, Ndt, PHIwt] = regularized_lda(corpus_fname,vocab_fname, reg, regmatrix_fname, config);
pause;

%% run Quad-Reg LDA
reg = 1;  regmatrix_fname   = 'Sdiag1.mat';
rand('state', 7);
[Nwt, Ndt, PHIwt] = regularized_lda(corpus_fname, vocab_fname, reg, regmatrix_fname, config);
pause;


%% run Conv-Reg LDA
reg = 2;  regmatrix_fname   = 'Scon1.mat';
rand('state', 7); 
[Nwt, Ndt, PHIwt] =  regularized_lda(corpus_fname, vocab_fname, reg, regmatrix_fname, config);



