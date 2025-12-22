
count_reads() {
    for fq in "$@"; do
        if [[ "$fq" == *.gz ]]; then
            echo $(( $(zcat "$fq" | wc -l) / 4 ))
        else
            echo $(( $(wc -l < "$fq") / 4 ))
        fi
    done
}

pair_fq_count_check () {
    fastq1="$1"
    fastq2="$2"
    output_file="$3"

    reads1=$(count_reads "$fastq1")
    reads2=$(count_reads "$fastq2")
    echo $fastq1 $reads1
    echo $fastq2 $reads2

    if [[ "$reads1" -eq "$reads2" ]]; then
        echo "$reads1" > "$output_file"
    fi
}
