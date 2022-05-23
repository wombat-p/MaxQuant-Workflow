process PREPARE_FILES { 
    publishDir "${params.outdir}/sdrfmerge"
    label 'process_medium'
    conda (params.enable_conda ? "conda-forge::python-3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "wombatp/maxquant-pipeline:dev"
    }


    input:
      path sdrf
      path parameters
      path exp_design
      path raws
      path map
     

    output:
      path "sdrf_local.tsv" , includeInputs:true         , emit: sdrf_local
      path "exp_design.txt" , includeInputs:true	    , emit: exp_design
      path "params.yml" , includeInputs:true		    , emit: params
      path "*.version.txt" , includeInputs:true          , emit: version
      path "{*.raw,*.RAW}" , includeInputs:true	    , emit: raws


    script:
    """
    if [[ "$sdrf" != "no_sdrf" ]] 
    then    
        if [[ "$sdrf" != "sdrf_local.tsv" ]]
        then
	    cp "${sdrf}" sdrf_local.tsv
        fi
    fi	
    if [[ "$exp_design" == "no_exp_design" ]] 
    then
        if [[ "$sdrf" == "no_sdrf" ]] 
        then
            printf "raw_file\texp_condition" >> exp_design_in.txt
	    for a in "$raws"
	    do
	        printf "\n\$a\tA"
	    done
        else 
	    $baseDir/scripts/sdrf2exp_design.py
        fi
    elif [[ "$exp_design" != "exp_design_in.txt" ]]
    then
	cp "${exp_design}" exp_design_in.txt
    fi
    if [[ "$sdrf" == "no_sdrf" ]] 
    then
	$baseDir/scripts/exp_design2sdrf.py
    fi
    if [[ "$raws" == "no_raws" ]]
    then
        # Download all files from column file uri		
	for a in \$(awk -F '\t' -v column_val='comment[file uri]' '{ if (NR==1) {val=-1; for(i=1;i<=NF;i++) { if (\$i == column_val) {val=i;}}} if(val != -1) { if (NR!=1) print \$val} } ' "$sdrf")
	do
            echo "Downloading \$a\n"
	    wget "\$a"
        done
    fi

    if [[ "$parameters" == "no_params" ]]
    then
	printf "params:\n  None:  \nrawfiles: None\nfastafile: None" >  params.yml
    elif [[ "$parameters" != "params.yml" ]] 
    then
        cp "${parameters}" params.yml
    fi
    if [[ "$map" != "params2sdrf.yml" ]]
    then
        cp "${map}" params2sdrf.yml
    fi
    echo "See workflow version" > prepare_files.version.txt
    """
}
