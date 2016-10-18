# test help
shannon -h
# test without options
shannon -o mt.bed -c MT -i metadata.csv
# test query
shannon -q "stage=='MI'" -c MT -i metadata.csv -o mt.stage_MI.bed
# test with query and groupby
shannon -q "stage=='MI'" -i metadata.csv -g patient age -c MT -o mt.stage_MI.bed
