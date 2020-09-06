import os

samples = os.listdir("pdfs/")
samples = [x.replace('.pdf', '') for x in samples]

rule all:
    input:
        "comp_results/entities_img_dataframe.tsv",
        "comp_results/gn_dataframe.tsv"

rule pdftotext:
    input:
        "pdfs/{sample}.pdf"
    output:
        "txt/{sample}.txt"
    shell:
        "pdftotext -layout {input} {output}"

rule gnfinder:
    input:
        "txt/{sample}.txt"
    output:
        "gn_txt/{sample}_gn.txt"
    shell:
        "gnfinder find -c -l eng -s 4,12 {input} > {output}"

rule oscar:
    input:
        "pdfs/{sample}.pdf"
    output:
        "oscar_json/{sample}.json"
    shell:
        "oscarpdf2json {input} > {output}"

rule osra:
    input:
        "pdfs/{sample}.pdf"
    output:
        "osra_txt/{sample}.txt"
    shell:
        "osra {input} -w {output}"

rule postprocessing:
    input:
         gnfinder=expand("gn_txt/{sample}_gn.txt", sample=samples),
         oscar=expand("oscar_json/{sample}.json", sample=samples),
         osra=expand("osra_txt/{sample}.txt", sample=samples)
    params:
         consumer_key = config["csk"]
    output:
        entities="comp_results/entities_img_dataframe.tsv",
        gnames="comp_results/gn_dataframe.tsv"
    run:
        oscar_lst = input.oscar[0]
        gn_lst = input.gnfinder[0]
        osra_lst = input.osra[0]
        shell("postprocessing_script {oscar_lst} {gn_lst} {osra_lst} {params.consumer_key} {output.entities} {output.gnames}")

