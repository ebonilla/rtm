function make_Sdiag1(regmatrix_fname, vocab_fname, output_fname)
%%% make regularization matrix S for QuadReg method
%%% David Newman (newman@uci.edu)
%%% Edwin V. Bonilla (edwin.bonilla@nicta.com.au)
%%% Last udpdate: 15/03/2012
%
% INPUT:
%   - regmatrix_fname: Text file name of regularization matrix.
%       The format is: word_id1, word_id2, pmi_value
%   - vocab_fname: Name of the text file containing the vocabulary
%   - output_fname: Name of matlab file where to save the output matrix
%
% Saves output matrix 'S' in output_fname

word = textread(vocab_fname,'%s');
W = length(word);

[w1,w2,ct] = textread(regmatrix_fname,'%d%d%f','headerlines',3);
  ct(find(ct<0))=0; %%% zero PMIs that are negative
  S = sparse(w1,w2,ct,W,W);
  S = S + S';
  for w=1:W
    S(w,w)=1;
  end
save(output_fname, 'S');




