configfile: 'config/config.yml'

ncores = config['ncores']

rule final:
    input:
        expand("results/{cytokine}/summary_report.html", cytokine=cytokine)

subworkflow otherworkflow:
    workdir:
        "../snakemake-cytokine-analysis"
    snakefile:
        "../snakemake-cytokine-analysis/Snakefile"
    configfile:
        "../snakemake-cytokine-analysis/config/config.yml"

rule a:
    input:
        otherworkflow("test.txt")
    output: ...
    shell:  ...
