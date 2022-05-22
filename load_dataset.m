function [X, Y]=load_dataset(dataset_num)
    if dataset_num==1
        load data/Yale_1024.mat
        Y=y;
    end
end