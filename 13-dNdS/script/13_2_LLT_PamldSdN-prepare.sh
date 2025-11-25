for i in $(seq 1 1577); do
    find CDS_${i}_core -type f \
        > comparisons_CDS_${i}_core_temp.txt

    awk -v i=$i '{print $1 "\ttemp/CDS_" i "_core/HEL_BA9262AA_AS.aln"}' \
        comparisons_CDS_${i}_core_temp.txt \
        > comparisons_CDS_${i}_core.txt
done
