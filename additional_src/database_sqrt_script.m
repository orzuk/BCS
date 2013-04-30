% a  one-time script that converts the database to single

X  = load(fullfile(data_dir, 'data.mat')); % load data
X.predbinSqr = single(X.predbinSqr);

