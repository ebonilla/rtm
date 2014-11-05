function [Nwt, Ndt, PHIwt] = regularized_lda(corpus_fname,vocab_fname, reg, regmatrix_fname, config)
%%% Regularized LDA
%%% Newman, Bonilla, Buntine
%%% NIPS 2011
%%% David Newman (newman@uci.edu)
%%% Edwin V. Bonilla (edwin.bonilla@nicta.com.au)
%%% Last Update: 15/03/2012
%
% INPUT:
%   - corpus_fname: Name of the text file containing the corpus data
%      It must have the following format:
%       * First three lines: Number of documents, size of vocabulary and       
%         total number of woeds
%       * Rest: three columns containing document-id, word-id, word-count
%   - vocab_fname: Name of the text file containing the vocabulary
%   - reg: Regularization method:
%       * reg = 0;   run standard LDA
%       * reg = 1;   run Quad-Reg LDA
%       * reg = 2;   run Conv-Reg LDA
%   - reg_matrix: Name of the matlab file containg the regularization 
%      matrix 'S'. See make_reg_matrix.m  details.
%   - config: A structure containing the parameters for the models:
%       .T          : Number of topics
%       .beta       : Dirichlet hyerparameter over words
%       .alpha      : Dirichlet hyerparameter over topics
%       .nu         : Regularization parameter for QuadReg
%       .gibbs_iter : Number of Gibbs iterations  
%       .lag_iter   : Interval on which the regularization is applied   
%       .reg_iter   : Regularization iterations  
%
% OUTPUT:
%   - Nwt: Matrix of counts on word w assigned to topic t
%   - Ndt: Matrix of counts on document d assignet to topic t
%   - PHIwt: Word topic distributions

beta         = config.beta;
alpha        = config.alpha;  %% later set alpha = 0.05 * N / (D*T)
gibbs_iter   = config.gibbs_iter;
lag_iter     = config.lag_iter;
reg_iter     = config.reg_iter;
T            = config.T;
nu           = config.nu;

%--------------------------------------------
% read corpus
%--------------------------------------------
[num] = textread(corpus_fname,'%d',3);
[did,wid,cnt] = textread(corpus_fname,'%d%d%d','headerlines',3);
D = num(1);         % number of docs
W = num(2);         % size of vocab
N = sum(cnt);       % total number of words
[word] = textread(vocab_fname,'%s');
assert(length(word)==W);

alpha = 0.05 * N / (D*T);

%--------------------------------------------
% reg
%--------------------------------------------
if (reg>0)
    load(regmatrix_fname, 'S');
  if (reg==1)
    nu = nu * (N/T);
  end
end

%--------------------------------------------
% allocate memory
%--------------------------------------------
w = zeros(N,1);
d = zeros(N,1);
z = zeros(N,1);
Nwt = zeros(W,T);
Ndt = zeros(D,T);

%--------------------------------------------
% fill w and d
%--------------------------------------------
count = 1;
for j = 1:length(cnt)
  for i = count:count+cnt(j)-1 
    w(i) = wid(j);
    d(i) = did(j);
  end
  count = count + cnt(j);
end
assert(max(w)==max(wid))
assert(max(d)==max(did))
assert(count-1==N)

%--------------------------------------------
% random initial assignment
%--------------------------------------------

z = floor(T*rand(N,1)) + 1;
for n = 1:N
  Nwt(w(n),z(n)) = Nwt(w(n),z(n)) + 1;
  Ndt(d(n),z(n)) = Ndt(d(n),z(n)) + 1;
end
Nt = sum(Nwt,1);
Ntchk = sum(Ndt,1);
assert(norm(Nt-Ntchk)==0);
assert(sum(Nt)==N);

PHIwt    = zeros(W,T);
GAMMA    = beta;
psi_val0 = zeros(W,1) + GAMMA;
for t = 1 : T
  PHIwt(:,t) = Nwt(:,t) + GAMMA;
  PHIwt(:,t) = PHIwt(:,t)/sum(PHIwt(:,t));
end

%--------------------------------------------
% Gibbs sampler
%--------------------------------------------
for iter = 1 : gibbs_iter

  Nt = sum(Nwt,1);
  if (reg==0), [z, Nwt, Ndt] = gibbsmex_lda (z, Nwt, Ndt, Nt, w, d);
  else        [z, Nwt, Ndt] = gibbsmex_semi(z, Nwt, Ndt,     w, d, alpha, PHIwt);
  end
  Nt = sum(Nwt,1);

  %%% Quad-Reg LDA
  if ( reg==1 && rem(iter,lag_iter)==0)
    for t = 1:T
      phi_t = Nwt(:,t) + 0.01;
      phi_t = phi_t / sum(phi_t);
      for kk = 1 : reg_iter
        num         = S * phi_t;
        den         = phi_t' * num;
        beta_t      = 2 * nu * phi_t.*num / den;
        phi_t = Nwt(:,t) + beta_t;
        phi_t = phi_t / sum(phi_t);
      end
      PHIwt(:,t) = phi_t;
      PHIwt(:,t) = PHIwt(:,t)/sum(PHIwt(:,t));
    end
  end
  
  %%% Conv-Reg LDA
  if ( reg==2 && rem(iter,lag_iter)==0)
    for t = 1:T
      Nw   = Nwt(:,t);
      psi_val = psi_val0/sum(psi_val0);
      for kk = 1 : reg_iter
        psi_val = (S*(Nw./(S*psi_val))).*psi_val + GAMMA;
        psi_val = psi_val/sum(psi_val);
      end
      PHIwt(:,t) = full(S*psi_val);
      PHIwt(:,t) = PHIwt(:,t)/sum(PHIwt(:,t));
    end
  end

  %%% print topics
  if ( rem(iter,lag_iter)==0)
    fprintf('\n\niter %d\n', iter);
    switch reg
      case 0, fprintf('>>> LDA\n')
      case 1, fprintf('>>> Quad-Reg LDA\n')
      case 2, fprintf('>>> Conv-Reg LDA\n')
    end
    if (reg==0)
      PHIwt = Nwt+beta;
    else
      PHIwt = (N/T)*PHIwt;
    end
    for t = 1:T
      fprintf('[t%d] ', t);
      [xsort,isort] = sort(-PHIwt(:,t));
      for ww = 1:8
        www = isort(ww);
        fprintf('%s(%d) ', word{www}, Nwt(www,t) );
      end
      fprintf('\n');
    end
  end
  
  %fprintf('iter %d\n',iter);
  
end % gibbs_iter

