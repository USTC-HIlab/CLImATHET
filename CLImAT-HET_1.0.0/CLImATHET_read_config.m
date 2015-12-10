function CLImAT_read_config(configfile)
%2014/06/12 by Zhenhua

global td_th
global gc_th
global map_th

% set parameters to default values
td_th = [10 300];
gc_th = [0.1 0.9];
map_th = [0 0.98];

fp_conf = fopen(configfile, 'r');
if fp_conf == -1
    error('Can not open the config file, please check it again!');
end
while 1
    confline = fgetl(fp_conf);
    if confline == -1
        break;
    elseif isempty(confline)
        continue;
    end
    if confline(1) == '/' || confline(1) == '#'
        continue;
    end
    result = regexp(confline,'minDepth\s*=\s*(\S+)','tokens','once');
    if ~isempty(result)
        td_th(1) = str2double(result{1});
    end
    result = regexp(confline,'maxDepth\s*=\s*(\S+)','tokens','once');
    if ~isempty(result)
        td_th(2) = str2double(result{1});
    end
    result = regexp(confline,'minGC\s*=\s*(\S+)','tokens','once');
    if ~isempty(result)
        gc_th(1) = str2double(result{1});
    end
    result = regexp(confline,'maxGC\s*=\s*(\S+)','tokens','once');
    if ~isempty(result)
        gc_th(2) = str2double(result{1});
    end
    result = regexp(confline,'minMapScore\s*=\s*(\S+)','tokens','once');
    if ~isempty(result)
        map_th(1) = str2double(result{1});
    end
    result = regexp(confline,'maxMapScore\s*=\s*(\S+)','tokens','once');
    if ~isempty(result)
        map_th(2) = str2double(result{1});
    end
end
fclose(fp_conf);

end