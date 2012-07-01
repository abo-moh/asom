% this function performs evaluation on the results given in the folders,
% karssemeijer, canny and blackbox.

% pectoral sensitivity, pectoral specificity, breast sensitivity, breast
% specificity, total sensitivity and total specificity
[pse psp bse bsp tse tsp nancount] = standard('karssemeijer/', 'ground/');
disp('karssemeijer');
disp(nancount);
figure('name', 'Karssemeijer PSE'); imhist(pse);
figure('name', 'Karssemeijer PSP'); imhist(psp);
figure('name', 'Karssemeijer BSE'); imhist(bse);
figure('name', 'Karssemeijer BSP'); imhist(bsp);
figure('name', 'Karssemeijer TSE'); imhist(tse);
figure('name', 'Karssemeijer TSP'); imhist(tsp);

disp('mk pse mean');
disp(mean(pse));
disp('mk pse var');
disp(var(pse));

disp('mk psp mean');
disp(mean(psp));
disp('mk psp var');
disp(var(psp));

disp('mk bse mean');
disp(mean(bse));
disp('mk bse var');
disp(var(bse));

disp('mk bsp mean');
disp(mean(bsp));
disp('mk psp var');
disp(var(bsp));

disp('mk tse mean');
disp(mean(tse));
disp('mk tse var');
disp(var(tse));

disp('mk tsp mean');
disp(mean(tsp));
disp('mk tsp var');
disp(var(tsp));

[pse psp bse bsp tse tsp nancount] = standard('canny/', 'ground/');
disp('canny');
disp(nancount);
figure('name', 'Canny PSE'); imhist(pse);
figure('name', 'Canny PSP'); imhist(psp);
figure('name', 'Canny BSE'); imhist(bse);
figure('name', 'Canny BSP'); imhist(bsp);
figure('name', 'Canny TSE'); imhist(tse);
figure('name', 'Canny TSP'); imhist(tsp);

disp('mc pse mean');
disp(mean(pse));
disp('mc pse var');
disp(var(pse));

disp('mc psp mean');
disp(mean(psp));
disp('mc psp var');
disp(var(psp));

disp('mc bse mean');
disp(mean(bse));
disp('mc bse var');
disp(var(bse));

disp('mc bsp mean');
disp(mean(bsp));
disp('mc bsp var');
disp(var(bsp));

disp('mc tse mean');
disp(mean(tse));
disp('mc tse var');
disp(var(tse));

disp('mc tsp mean');
disp(mean(tsp));
disp('mc tsp var');
disp(var(tsp));

[pse psp bse bsp tse tsp nancount] = standard('blackbox/', 'ground/');
disp('blackbox');
disp(nancount);
figure('name', 'Blackbox PSE'); imhist(pse);
figure('name', 'Blackbox PSP'); imhist(psp);
figure('name', 'Blackbox BSE'); imhist(bse);
figure('name', 'Blackbox BSP'); imhist(bsp);
figure('name', 'Blackbox TSE'); imhist(tse);
figure('name', 'Blackbox TSP'); imhist(tsp);

disp('mb pse mean');
disp(mean(pse));
disp('mb pse var');
disp(var(pse));

disp('mb psp mean');
disp(mean(psp));
disp('mb psp var');
disp(var(psp));

disp('mb bse mean');
disp(mean(bse));
disp('mb bse var');
disp(var(bse));

disp('mb bsp mean');
disp(mean(bsp));
disp('mb psp var');
disp(var(bsp));

disp('mb tse mean');
disp(mean(tse));
disp('mb tse var');
disp(var(tse));

disp('mb tsp mean');
disp(mean(tsp));
disp('mb tsp var');
disp(var(tsp));
