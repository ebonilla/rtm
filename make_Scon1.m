function make_Scon1(regmatrix_fname, vocab_fname, output_fname)
%%% make regularization matrix S for convReg method
%%% David Newman (newman@uci.edu)
%%% Edwin V. Bonilla (edwin.bonilla@nicta.com.au)
%%% Last udpdate: 15/03/2012
%
% INPUT:
%   - regmatrix_fname: Text file name of regularization matrix.
%       The format is: word_id1, word_id2, count_value
%   - vocab_fname: Name of the text file containing the vocabulary
%   - output_fname: Name of matlab file where to save the output matrix
%
% Saves output matrix 'S' in output_fname


word = textread(vocab_fname,'%s');
W = length(word);

[w1,w2,ct] = textread(regmatrix_fname,'%d%d%d','headerlines',3);

idx_diag = (1 : W)';
idx_nonzero = w1(w1==w2);
idx_diag(idx_nonzero) = [];

% adding extra diagonal elements
% w11 = [w1; idx_diag];
% w22 = [w2; idx_diag];
% ct = [ct; 0.5*ones(length(idx_diag),1)]; 
S = sparse(w1,w2,ct,W,W);
fprintf('Normalizing Co-occurence ...');
for i = 1 : length(w1)
  row = w1(i);
  col = w2(i);
  if (row ~= col)
    val = S(row,row) + S(col,col) - S(row,col);
    S(row,col) = S(row,col)/val;
  end 
end
diag_idx = sub2ind([W W], idx_nonzero, idx_nonzero);
S(diag_idx) = 0.5;
S = S + S';

for w=1:W
  S(w,w)=1;
end

save(output_fname, 'S');


