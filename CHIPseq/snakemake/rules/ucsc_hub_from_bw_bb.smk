
rule ucsc_hub_from_bw_bb:
    input:
        bw=expand("{output_dir}/ucsc/{sample}_normalized.bw", sample=sample_names,output_dir=output_dir),
        peaks=expand("{output_dir}/ucsc/{case}.peaks.bb", case=cases,output_dir=output_dir)
    output:
        hub_dir=directory("/data/wade/linl7/{identifier}/mm10/"),
        hub_file="/data/wade/linl7/{identifier}/hub.txt",
        genomes_file="/data/wade/linl7/{identifier}/genomes.txt",
        trackdb_file="/data/wade/linl7/{identifier}/mm10/trackDb.txt"
    params:
        sorted_sample_names = sorted(sample_names),
        sorted_cases_names = sorted(cases)
    shell:
        """
        #mkdir -p {output.hub_dir}
        cp {input.bw} {input.peaks} {output.hub_dir}
    
        # Create the main hub file
        echo "hub {identifier}" > {output.hub_file}
        echo "shortLabel {identifier}" >> {output.hub_file}
        echo "longLabel ChIP-seq Data Hub" >> {output.hub_file}
        echo "genomesFile genomes.txt" >> {output.hub_file}
        echo "email long.lin@nih.gov" >> {output.hub_file}
        
        # Specify the genome version
        echo "genome mm10" > {output.genomes_file}
        echo "trackDb mm10/trackDb.txt" >> {output.genomes_file}
        
        # Create trackDb.txt with all samples
        echo -e "track cont\ncontainer multiWig\nshortLabel {identifier}\nlongLabel {identifier}\ntype bigWig\nvisibility full\nalwaysZero on\nautoScale on\nshowSubtrackColorOnUi on\naggregate none\ncolor 1,1,1\npriority 1\n " > {output.trackdb_file}
        
        index=1


        for sample in {params.sorted_sample_names};
        do 
            echo "track $sample " >> {output.trackdb_file}
            echo "bigDataUrl ${{sample}}_normalized.bw" >> {output.trackdb_file}
            echo "shortLabel $sample" >> {output.trackdb_file}
            echo "longLabel $sample" >> {output.trackdb_file}
            echo "type bigWig" >> {output.trackdb_file}
            echo "visibility full" >> {output.trackdb_file}
            echo "alwaysZero on" >> {output.trackdb_file}
            echo -e "autoScale on\nmaxHeightPixels 1:60:9999\ncolor 20,20,255\nparent cont" >> {output.trackdb_file}
            ((index++))
            echo "priority $index" >> {output.trackdb_file} 
            echo "" >> {output.trackdb_file}
        done

        for case in {params.sorted_cases_names}
        do
            echo "track ${{case}}" >> {output.trackdb_file}
            echo "bigDataUrl ${{case}}.peaks.bb" >> {output.trackdb_file}
            echo "shortLabel $case" >> {output.trackdb_file}
            echo "longLabel $case" >> {output.trackdb_file}
            echo "type bigBed" >> {output.trackdb_file}
            echo "visibility dense" >> {output.trackdb_file}
            echo -e "color 200,0,0" >> {output.trackdb_file}
            ((index++))
            echo "priority $index" >> {output.trackdb_file}
            echo "" >> {output.trackdb_file}
        done
        """

