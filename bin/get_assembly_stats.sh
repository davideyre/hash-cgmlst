#run in root output folder, e.g. cd comparison_study_data/replicates_output
# get assembly stats
tail -n 1 -q */*_cgmlst.statlog  > tmp_assembly_stats.txt
cat */*_cgmlst.statlog | head -1 > tmp_header.txt
cat tmp_header.txt tmp_assembly_stats.txt > assembly_stats.txt
rm tmp_assembly_stats.txt tmp_header.txt

#get mlst
tail -n 1 -q */*_mlst.txt  > tmp_mlst.txt
cat */*_mlst.txt | head -1 > tmp_header.txt
cat tmp_header.txt tmp_mlst.txt > mlst_summary.txt
rm tmp_mlst.txt tmp_header.txt
