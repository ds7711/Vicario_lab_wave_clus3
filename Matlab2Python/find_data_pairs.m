function [fn_pairs, unpaired_files] = find_data_pairs(filenames)

% double loop 
num_names = length(filenames);
split_kw = '_';

fn_pairs = cell(10000, 2);
unpaired_files = cell(10000, 1);
unpaired_count = 0;
pair_count = 0;
paired_idxes = zeros(10000);
for iii = 1 : (num_names-1)
    if any(iii == paired_idxes(:))
        continue;
    end
    st_fn = filenames{iii};
    find_pair_flag = 0;
    st_fn_splited = strsplit(st_fn, split_kw);
    for jjj = (iii+1) : num_names
        
        if any(jjj == paired_idxes(:))
            continue;
        end
        
        ed_fn = filenames{jjj};
        ed_fn_splited = strsplit(ed_fn, split_kw);
        try 
            cmp_result = strcmp(st_fn_splited, ed_fn_splited);
        catch
            cmp_result = zeros(4, 1);
        end
        cmp_result(2) = ~ cmp_result(2);
        if all(cmp_result)
            % store the fn pairs and skip to the next one
            paired_idxes(pair_count * 2 + 1) = iii;
            paired_idxes(pair_count * 2 + 2) = jjj;
            pair_count = pair_count + 1;
            fn_pairs{pair_count, 1} = st_fn;
            fn_pairs{pair_count, 2} = ed_fn;
            find_pair_flag = 1;
            break;
        end
    end
    if ~ find_pair_flag
        unpaired_count = unpaired_count + 1;
        unpaired_files{unpaired_count} = st_fn;
    end
end
if pair_count ~= 0
    fn_pairs = fn_pairs(1:pair_count, 1:2);
else
    fn_pairs = cell(1, 2);
end
if unpaired_count ~= 0
    unpaired_files = unpaired_files(1:unpaired_count);
else
    unpaired_files = cell(1, 1);
end
end