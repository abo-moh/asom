function [pse, psp, bse, bsp, tse, tsp, nancount] = standard(dirpath, ground)
  diropen = dir(dirpath);
  pse = size(length(diropen));
  psp = size(length(diropen));
  bse = size(length(diropen));
  bsp = size(length(diropen));
  tse = size(length(diropen));
  tsp = size(length(diropen));
  nancount = 0;
  for i = 1:length(diropen)
    filename = diropen(i).name;
    if ~strcmp(filename, '.') && ~strcmp(filename, '..')
      %disp(['evaluating ' filename]);
      [pathstr name ext]= fileparts([dirpath diropen(i).name]);

      % find the mask
      mask_filename = [ground 'mask' name(2:length(name))];
      mask = imread(mask_filename, 'tif');

      % load the mat file
      mat = load([dirpath filename]);
      
      if isfield(mat, 'pec_a') && isfield(mat', 'pec_b')
        if isnan(mat.pec_a) || isnan(mat.pec_b)
          nancount = nancount + 1;
          continue
        end
      end
      
      % display real image
      %figure();
      %imshow(mat.sa_bound, [min(mat.sa_bound(:)) max(mat.sa_bound(:))]);
      
      % create a figure
      %figure('name', 'Test');

      %subplot(2, 2, 1);
      %imshow(mat.sa_bound, [min(mat.sa_bound(:)) max(mat.sa_bound(:))]);
      %title('Eval');

      %subplot(2, 2, 2);
      %imshow(mask, [min(mask(:)) max(mask(:))]);
      %title('Reference');

      %subplot(2, 2, 3);
      %imshow(eval_pectoral);
      %title('Eval Pectoral');

      %subplot(2, 2, 4);
      %imshow(mask_pectoral);
      %title('Reference Pectoral');

      % run the pectoral evaluation
      eval_pectoral = mat.sa_bound == 2;
      mask_pectoral = mask == 2;

      [p_difference p_sensitivity p_specificity] = evaluate(eval_pectoral, mask_pectoral);
      if isnan(p_sensitivity) || isnan(p_specificity)
        disp(filename);
        figure(); imshow(mask, [min(mask(:)) max(mask(:))]);
        disp(min(mask(:)));
        disp(max(mask(:)));
        continue;
      end
      pse(i) = p_sensitivity;
      psp(i) = p_specificity;
      
      % run the breast evaluation
      eval_breast = mat.sa_bound == 1;
      mask_breast = mask == 1;

      [b_difference b_sensitivity b_specificity] = evaluate(eval_breast, mask_breast);
      if isnan(b_sensitivity) || isnan(b_specificity)
        disp(filename);
        figure(); imshow(mask, [min(mask(:)) max(mask(:))]);
        disp(min(mask(:)));
        disp(max(mask(:)));
        continue;
      end
      bse(i) = b_sensitivity;
      bsp(i) = b_specificity;
      
      % run the total evaluation
      eval_total = eval_pectoral | eval_breast;
      mask_total = mask_pectoral | mask_breast;
 
      [t_difference t_sensitivity t_specificity] = evaluate(eval_total, mask_total);
      if isnan(t_sensitivity) || isnan(t_specificity)
        disp(filename);
        figure(); imshow(mask, [min(mask(:)) max(mask(:))]);
        disp(min(mask(:)));
        disp(max(mask(:)));
        continue;
      end
      tse(i) = t_sensitivity;
      tsp(i) = t_specificity;
    end
  end
end
