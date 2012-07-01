function [difference, sensitivity, specificity] = evaluate(candidate, reference)
  % create a new figure
  figure('name', 'Evaluation');

  % read the candidate image
  i_candidate = imread(candidate, 'tif');

  % here we call a function which returns the slope and offset of the
  % linear function that bounds the pectoral muscle
  a = -2;
  b = 950;

  % fetch information on the candidate image, we use the pectoral
  % function to draw a polygon mask, specifically a polygon limited to
  % the top left hand corner, and a linear function.
  [w, h] = size(i_candidate);
  i_candidate_bw = pectoral(w, h, a, b);

  % plotting
  subplot(2, 2, 1);
  imshow(i_candidate_bw);
  title('Candidate');

  % our masks are oriented to the right hand side, so we flip em to match
  % the left hand side candidates
  i_reference = fliplr(imread(reference, 'tif'));

  % in the masks, white pixels equal 2, therefore, we'd like a matrix of
  % those cells exclusively
  i_reference_bw = i_reference == 2;

  % plotting
  subplot(2, 2, 2);
  imshow(i_reference_bw);
  title('Reference');

  % summerize the amount of white pixels in the candidate and reference
  % images
  n_candidate_white = sum(i_candidate_bw(:));
  n_reference_white = sum(i_reference_bw(:));

  % we know from Woods, Sallam and Bowyer that sensitivity is false
  % negatives
  i_sensitivity_bw = ~i_candidate_bw & i_reference_bw;

  % plotting
  subplot(2, 2, 3);
  imshow(i_sensitivity_bw);
  title('Sensitivity');

  % we know from Woods, Sallam and Bowyer that specificity is false
  % positives
  i_specificity_bw = i_candidate_bw & ~i_reference_bw;

  % plotting
  subplot(2, 2, 4);
  imshow(i_specificity_bw);
  title('Specificity');

  % the difference isn't really that interesting, but gives us an idea on
  % the scope
  difference = n_candidate_white - n_reference_white;

  % we compute sensitivity as a percentage of false negatives in the
  % reference muscle space
  sensitivity = 1.0 - (sum(i_sensitivity_bw(:)) / n_reference_white);

  % and we find specificity as the percetange of false positives in the
  % candidate muscle space
  specificity = 1.0 - (sum(i_specificity_bw(:)) / n_candidate_white);
end
