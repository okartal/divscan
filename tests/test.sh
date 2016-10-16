# test help
shannon -h
# test without options
shannon -i CpG/*.gz -c MT
# test without options with custom output
shannon -o MT.temp.bed -c MT -i CpG/*.gz
# test query without metadata
shannon -q "stage=='MI'" -c MT -i CpG/*.gz
# test query with metadata
shannon -q "stage=='MI'" -m metadata.csv -c MT -i CpG/*.gz
# test with query and groupby
shannon --query "stage=='MI'" --metadata metadata.csv --groupby patient age --contig MT --input CpG/*.gz
