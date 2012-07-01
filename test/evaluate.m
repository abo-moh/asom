function [difference, sensitivity, specificity] = evaluate(i_eval_bw, i_reference_bw)
  n_eval_white = sum(i_eval_bw(:));
  n_reference_white = sum(i_reference_bw(:));
  
  i_tp_bw = i_eval_bw & i_reference_bw;
  i_fn_bw = ~i_eval_bw & i_reference_bw;
  
  i_tn_bw = ~i_eval_bw & ~i_reference_bw;
  i_fp_bw = i_eval_bw & ~i_reference_bw;
  
  %figure('name', 'Evaluation');

  %subplot(1, 2, 1);
  %imshow(i_fn_bw);
  %title('Falske negativer');

  %subplot(1, 2, 2);
  %imshow(i_fp_bw);
  %title('Falske positiver');
  
  n_tp = sum(i_tp_bw(:));
  n_fn = sum(i_fn_bw(:));
  n_tn = sum(i_tn_bw(:));
  n_fp = sum(i_fp_bw(:));

  difference = n_eval_white - n_reference_white;

  % sensitivity is found as the percentage of true positives amongst
  % true positives and false negatives
  sensitivity = n_tp / (n_tp + n_fn);

  % and we find specificity as the percentage of true negatives amongst
  % true negatives and false positives
  specificity = n_tn / (n_tn + n_fp);
end
