%%% make regularization matrix S for conv-dir-reg

clear all

word = textread('vocab.txt','%s');
W = length(word)

[w1,w2,ct] = textread('w1w2cooc.txt','%d%d%d','headerlines',3);

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

save Scon.mat S
