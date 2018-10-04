# Snakemake hemophilia
localrules: all


rule all:
  input:
    expand(directory("output/{dataset}/{gatk}/report"),gatk=config["GATK"],dataset=config["datasets"].keys())

rule splithemophilia:
  input:
    R1="input/{dataset}/Undetermined_S0_L00{lane}_R1_001.fastq.gz",
    I="input/{dataset}/Undetermined_S0_L00{lane}_I1_001.fastq.gz",
    R2="input/{dataset}/Undetermined_S0_L00{lane}_R2_001.fastq.gz",
    lst="input/{dataset}/sample_index.lst"
  output:
    bam="output/{dataset}/mapping/sample_l{lane}.bam"
  log:
    "output/{dataset}/mapping/processing_stats_l{lane}.log"
  conda: "envs/rules.yml"
  shell:"""
    ( paste <( zcat {input.R1} | cut -c 1-120 ) \
      <( zcat {input.I} ) \
      <( zcat {input.R2} | cut -c 1-120 ) | \
    awk '{{ count+=1; if ((count == 1) || (count == 3)) {{ print $1 }} else {{ print $1$2$3 }}; if (count == 4) {{ count=0 }} }}' | \
    scripts/pipeline2.0/SplitFastQdoubleIndexBAM.py --bases_after_index=ATCTCGTATGCCGTCTTCTGCTTG --bases_after_2ndindex='' -l 10 -m 0 -s 131 --summary -i {input.lst} -q 10 -p --remove | scripts/pipeline2.0/MergeTrimReadsBAM.py --mergeoverlap -p \
    > {output.bam} ) 2> {log}
    """

def getLaneBAMs(wc):
  return(expand("output/{dataset}/mapping/sample_l{lane}.bam", dataset=wc.dataset,lane=config["datasets"][wc.dataset]["lanes"]))

rule mergebam:
  input: getLaneBAMs
  output: "output/{dataset}/mapping/sample.bam"
  conda: "envs/rules.yml"
  shell: "samtools merge -c {output} {input}"

rule reheadering:
  input:sam="input/sam_header_hg19_1000g.sam",
        lst="input/{dataset}/sample_index.lst"
  output:"input/{dataset}/new_header.sam"
  conda: "envs/rules.yml"
  shell:"""
        ( cat {input.sam}; tail -n +2 {input.lst}  | awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ print "@RG","ID:"$2,"PL:Illumina","LB:"$2,"SM:"$2}}'     ) > {output}
        """
def loadSamples(wc):
    file = open("input/%s/sample_index.lst" % wc.dataset, "r")
    output=[]
    for line in file:
      if line.startswith("#"):
          continue
      output.append(line.split("\t")[1].strip())
    return(output)

rule bysample:
  input:
    bam="output/{dataset}/mapping/sample.bam"
  output:
    "output/{dataset}/mapping/by_sample/{plate}.bam"
  params:
    plate="{plate}"
  conda: "envs/rules.yml"
  shell: """
    samtools view -u -F 513 -r {params.plate} {input.bam} | scripts/pipeline2.0/FilterBAM.py -q --qual_number 5 --qual_cutoff=15 -p > {output}
    """

rule aligning:
  input:
    bam="output/{dataset}/mapping/by_sample/{plate}.bam",
    fasta=config["references"]["bwa"],
    design="input/{dataset}/hemomips_design.txt",
    new_header="input/{dataset}/new_header.sam"
  output: "output/{dataset}/mapping/aligned/{plate}.bam"
  conda: "envs/rules.yml"
  shell:"""
    bwa mem -L 80 -M -C {input.fasta} <( samtools view -F 513 {input.bam} | awk 'BEGIN{{ OFS="\\n"; FS="\\t" }}{{ print "@"$1"\\t"$12"\\t"$13"\\t"$14,$10,"+",$11 }}' ) | samtools view -u - | samtools sort - | scripts/pipeline2.0/TrimMIParms.py -d {input.design} -p | samtools reheader {input.new_header} - | samtools sort -o {output} -
    """

rule indexing:
  input: "output/{dataset}/mapping/aligned/{plate}.bam"
  output: "output/{dataset}/mapping/aligned/{plate}.bai"
  conda: "envs/rules.yml"
  shell: "samtools index {input} {output}"

def sampleBamsAligned(wc):
    return (expand("output/{dataset}/mapping/aligned/{plate}.bam", dataset=wc.dataset, plate=loadSamples(wc)))
def sampleBamsAlignedIdx(wc):
    return (expand("output/{dataset}/mapping/aligned/{plate}.bai", dataset=wc.dataset, plate=loadSamples(wc)))

rule samplesexcheck:
  input:
    lst="input/{dataset}/sample_index.lst",
    bams=sampleBamsAligned,
    idx=sampleBamsAlignedIdx
  output:"output/{dataset}/samples_sex_check.txt"
  conda: "envs/rules.yml"
  shell: "( for i in {input.bams}; do echo $( basename $i ) $(samtools view $i Y | wc -l) $(samtools view -F u $i | wc -l); done )> {output}"

rule inversionmips:
  input:bam="output/{dataset}/mapping/sample.bam",
        sam="input/{dataset}/new_header.sam",
        inv=config["references"]["inv"]
  output: "output/{dataset}/mapping/inversion_mips/{plate}.bam"
  params:
    plate="{plate}"
  conda: "envs/rules.yml"
  shell: """
    (   grep "@RG" {input.sam}; \
        bwa mem -M -L 80 -C {input.inv} <(
            samtools view -r {params.plate} -F 1 {input.bam} | awk 'BEGIN{{ OFS="\\n" }}{{ if (length($10) >= 75) {{ print "@"$1" "$12"\\t"$13"\\t"$14,$10,"+",$11 }} }}' \
        ) | awk '{{ if (($0 ~ /^@/) || ($3 ~ /^inv/)) print }}'; \
        bwa mem -M -L 80 -p -C {input.inv} <( \
            samtools view -r {params.plate} -f 1 {input.bam} | awk 'BEGIN{{ OFS="\\n" }}{{ print "@"$1" "$12"\\t"$13"\\t"$14,$10,"+",$11 }}' \
        ) | awk '{{ if (($0 !~ /^@/) && ($3 ~ /^inv/)) print }}' \
     ) | samtools view -b -F 768 - | samtools sort -O bam -o {output} -
    """
def sampleBamsInversion(wc):
    return (expand("output/{dataset}/mapping/inversion_mips/{plate}.bam", dataset=wc.dataset, plate=loadSamples(wc)))

rule inversionsum:
  input:
    lst="input/{dataset}/sample_index.lst",
    bam=sampleBamsInversion
  output:"output/{dataset}/mapping/inversion_mips/inversion_summary_counts.txt"
  conda: "envs/rules.yml"
  shell: """
    ( for i in $(tail -n +2 {input.lst} | cut -f 2 ); do echo $i $( ( samtools view -F 513 {input.bam} | awk 'BEGIN{{ FS="\\t" }}{{ split($12,a,":"); if (($6 !~ /S/) && (a[1] == "NM") && (a[3] <= 10)) {{ print $3 }} }}'; samtools view -f 2 -F 512 {input.bam} | awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ split($12,a,":"); if (($6 !~ /S/) && (a[1] == "NM") && (a[3] <= 10)) {{ print $1,$3 }} }}' | sort | uniq -c | awk '{{ if ($1 == 2) print $3 }}' ) | sort | uniq -c | awk '{{ print $1":"$2 }}' ); done )> {output}
    """

############################################
# Realignment and Variant Calling
############################################

# GATK4
rule gatk4_HTcaller:
  input:
    bamin=sampleBamsAligned,
    baminIdx=sampleBamsAlignedIdx,
    targets="input/{dataset}/targets.intervals",
    fasta=config["references"]["fasta"]
  output:
    bamout="output/{dataset}/mapping/gatk4/realign_all_samples.bam",
    vcf="output/{dataset}/mapping/gatk4/bam.vcf.gz"
  conda:"envs/python3gatk4.yml"
  shell: "gatk HaplotypeCaller -R {input.fasta} -L {input.targets} $(ls -1 {input.bamin} | xargs -n 1 echo -I ) --output-mode EMIT_ALL_SITES -bamout {output.bamout} -O {output.vcf} --disable-optimizations"

rule gatk4_gvcfs:
  input:
    lst="input/{dataset}/sample_index.lst",
    fasta=config["references"]["fasta"],
    targets="input/{dataset}/targets.intervals",
    bam="output/{dataset}/mapping/aligned/{plate}.bam",
    idx="output/{dataset}/mapping/aligned/{plate}.bai",
  output:
    vcfgz="output/{dataset}/mapping/gatk4/gvcf/{plate}.g.vcf.gz"
  params:
    plate="{plate}"
  conda:"envs/python3gatk4.yml"
  shell:"""
    gatk --java-options "-Xmx8G" HaplotypeCaller --max-reads-per-alignment-start 0 --disable-optimizations \
    --sample-name {params.plate} -ERC GVCF -R {input.fasta} -I {input.bam} -L {input.targets} -O {output.vcfgz}
    """

def sampleBamsInversion(wc):
    return (expand("output/{dataset}/mapping/gatk4/gvcf/{plate}.g.vcf.gz", dataset=wc.dataset, plate=loadSamples(wc)))

rule gatk4_combine:
  input:fasta=config["references"]["fasta"],
        gvcf=sampleBamsInversion
  output:bamout="output/{dataset}/mapping/gatk4/realign_all_samples.all_sites.vcf.gz"
  conda:"python3gatk4.yml"
  shell:"gatk CombineGVCFs --break-bands-at-multiples-of 1 -R {input.fasta} $(ls -1 {input.gvcf} | xargs -n 1 echo -V ) -O {output.bamout}"


rule gatk4_genotype:
  input:
    fasta=config["references"]["fasta"],
    vcf="output/{dataset}/mapping/gatk4/realign_all_samples.all_sites.vcf.gz"
  output:"output/{dataset}/mapping/gatk4/realign_all_samples.vcf.gz"
  conda:"envs/python3gatk4.yml"
  shell:"gatk GenotypeGVCFs -R {input.fasta} -V {input.vcf} -O {output}"

# GATK3

rule gatk3_realign:
  input:
    bam=sampleBamsAligned
  output:"output/{dataset}/mapping/gatk3/realign_all_samples.bam"
  shell: "java -Xmx8G -jar tools/GATK322/GenomeAnalysisTK.jar -T IndelRealigner -R {fa} -DBQ 3 -filterNoBases -maxReads 1500000 -maxInMemory 1500000 -targetIntervals targets_split.intervals $(ls -1 {input.bam} | xargs -n 1 echo -I ) -o {output} -dt BY_SAMPLE -dcov 500"

rule gatk3_genotyping:
  input:bam="output/{dataset}/mapping/gatk3/realign_all_samples.bam",
        targets="input/{dataset}/targets.intervals",
        fasta=config["references"]["fasta"]
  output:"output/{dataset}/mapping/gatk3/realign_all_samples.all_sites.vcf.gz"
  conda: "envs/rules.yml"
  shell: "java -Xmx6G -jar tools/GATK3446/GenomeAnalysisTK.jar -T UnifiedGenotyper -R {input.fasta} -I {input.bam} -L {input.targets} -o >( bgzip -c > {output} ) -glm BOTH -rf BadCigar --max_alternate_alleles 15 --output_mode EMIT_ALL_SITES -dt NONE"

rule gatk3_subsetting:
    input:"output/{dataset}/mapping/gatk3/realign_all_samples.all_sites.vcf.gz"
    output:"output/{dataset}/mapping/gatk3/realign_all_samples.vcf.gz"
    conda: "envs/rules.yml"
    shell: """
      zcat {input} | awk 'BEGIN{{ FS="\\t" }}{{ if ($1 ~ /^#/) {{ print }} else {{ if ($5 != ".") print }} }}' | bgzip -c > {output}"""


rule gatk3_tabixing:
    input:"output/{dataset}/mapping/gatk3/realign_all_samples.vcf.gz"
    output:"output/{dataset}/mapping/gatk3/realign_all_samples.vcf.idx"
    conda: "envs/rules.yml"
    shell: "tabix -p vcf {input} > {output}"

rule gatk3_tabixing2:
    input:"output/{dataset}/mapping/gatk3/realign_all_samples.all_sites.vcf.gz"
    output:"output/{dataset}/mapping/gatk3/realign_all_samples.all_sites.vcf.idx"
    conda: "envs/rules.yml"
    shell: "tabix -p vcf {input} > {output}"


##############################################
# Statistics and Check
##############################################


rule MIPstats:
  input:"output/{dataset}/mapping/{gatk}/realign_all_samples.bam"
  output:"output/{dataset}/mapping/{gatk}/realign_all_samples.MIPstats.tsv"
  conda: "envs/rules.yml"
  shell:"scripts/pipeline2.0/MIPstats.py {input} -o {output}"

rule checkPileUp:
  input: bam="output/{dataset}/mapping/{gatk}/realign_all_samples.bam",
         sites="output/{dataset}/mapping/{gatk}/realign_all_samples.vcf.gz"
  output:"output/{dataset}/mapping/{gatk}/realign_all_samples.indel_check.txt"
  conda: "envs/rules.yml"
  shell: "scripts/processing/checkPileUpAtInDels.py -b {input.bam} -s {input.sites} -o {output}"

##############################################
# VEP
##############################################


rule VEP:
  input:vcf="output/{dataset}/mapping/{gatk}/realign_all_samples.vcf.gz",
        fasta2=config["references"]["fasta2"],
        veppath=config["tools"]["vep"],
        vepversion=config["parameters"]["vep-version"]
  output:"output/{dataset}/mapping/{gatk}/realign_all_samples.vep.tsv.gz"
  shell: """zcat {input.vcf} | scripts/processing/VCF2vepVCF.py | perl {input.veppath} --no_stats --fasta {input.fasta2} --quiet --buffer 2000 --cache --offline --species homo_sapiens --db_version={input.vepversion} --format vcf --symbol --hgvs --regulatory --gmaf --sift b --polyphen b --ccds --domains --numbers --canonical --shift_hgvs --output_file >( awk 'BEGIN{{ FS="\\t"; OFS="\\t"; }}{{ if ($1 ~ /^##/) {{ print }} else if ($1 ~ /^#/) {{ sub("^#","",$0); print "#Chrom","Start","End",$0 }} else {{ split($2,a,":"); if (a[2] ~ /-/) {{ split(a[2],b,"-"); print a[1],b[1],b[2],$0 }} else {{ print a[1],a[2],a[2],$0 }} }} }}' | grep -E "(^#|ENST00000360256|ENST00000218099)" | bgzip -c > {output} ) --force_overwrite
  """

##############################################
# Final summary report
##############################################


rule summaryreport_gatk4:
  input:vcf="output/{dataset}/mapping/{gatk}/realign_all_samples.all_sites.vcf.gz",
        vep="output/{dataset}/mapping/{gatk}/realign_all_samples.vep.tsv.gz",
        inv="output/{dataset}/mapping/inversion_mips/inversion_summary_counts.txt",
        sex="output/{dataset}/samples_sex_check.txt",
        target="input/{dataset}/target_coords.bed",
        mips="output/{dataset}/mapping/{gatk}/realign_all_samples.MIPstats.tsv",
        indel="output/{dataset}/mapping/{gatk}/realign_all_samples.indel_check.txt",
        hemomips="input/{dataset}/hemomips_design.txt",
        tg=config["references"]["annotation"]
  output:
        folder=directory("output/{dataset}/{gatk}/report")
  conda: "envs/rules.yml"
  shell: """scripts/processing/summary_report.py --vcf {input.vcf} --vep {input.vep} --inversions {input.inv} --sample_sex {input.sex} --target {input.target} --mipstats {input.mips} --indelCheck {input.indel} --design {input.hemomips} --TG {input.tg} && mv report/ {output.folder}"""
