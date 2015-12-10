function corrected_rc = CLImATHET_GC_MAP_correction(data_rc,data_gc,data_map)
% 10/11/2014 by Zhenhua
% this function is used to correct read counts for GC-content and mappability

corrected_rc = data_rc;

rc_med = median(data_rc);
int_map = floor(data_map*100);
int_gc = floor(data_gc*100);
int_map_u = unique(int_map);
int_gc_u = unique(int_gc);

for i = 1:length(int_map_u)
    for j = 1:length(int_gc_u)
        tv = int_map == int_map_u(i) & int_gc == int_gc_u(j);
        if sum(tv) > 0
            loc_med = median(data_rc(tv));
            if loc_med > 0
                corrected_rc(tv) = data_rc(tv)*rc_med/loc_med;
            end
        end
    end
end

end