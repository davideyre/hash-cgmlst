tail -n 1 -q replicates_output/*/*_cgmlst.statlog  > tmp_assembly_stats.txt
cat replicates_output/*/*_cgmlst.statlog | head -1 > tmp_header.txt
cat tmp_header.txt tmp_assembly_stats.txt > assembly_stats.txt
rm tmp_assembly_stats.txt tmp_header.txt
