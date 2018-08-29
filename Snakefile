# Snakemake hemophilia

rule all:
  input:
    expand("output/{dataset}/{gatk}/report/summary.html",gatk=config["GATK"],dataset=config["datasets"].keys())

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
  conda: "hemoMIPs.yml"
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
  conda: "hemoMIPs.yml"
  shell: "samtools merge -c {output} {input}"


def loadSamples(wc):
  file = file.read("input/%s/sample_index.lst" % wc.dataset)
  output=[]
  for line in file:
    output.append(line.split("\t")[1])
    return(output)

rule reheadering:
  input:sam="input/sam_header_hg19_1000g.sam",
        lst="input/{dataset}/sample_index.lst"
  output:"input/{dataset}/new_header.sam"
  conda: "hemoMIPs.yml"
  shell:"""
        ( cat {input.sam}; tail -n +2 {input.lst}  | awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ print "@RG","ID:"$2,"PL:Illumina","LB:"$2,"SM:"$2                }}'     ) > {output}
        """


rule bysample:
  input:bam="output/{dataset}/mapping/sample.bam",
        lst="input/{dataset}/sample_index.lst"
  output: expand("output/{{dataset}}/mapping/by_sample/{plate}.bam",plate=loadSamples)
  conda: "hemoMIPs.yml"
  shell: "for i in $( tail -n +2 {input.lst} | cut -f 2); do samtools view -u -F 513 -r ${{i}} {input.bam} | scripts/pipeline2.0/FilterBAM.py -q --qual_number 5 --qual_cutoff=15 -p > {output}"




rule aligning:
  input: 
    bam=expand("output/{{dataset}}/mapping/by_sample/{plate}.bam",plate=loadSamples),
    fasta=config["references"]["fasta"],
    design="input/hemomips_design.txt"
  output: expand("output/{{dataset}}/mapping/aligned/{plate}.bam",plate=loadSamples)
  conda: "hemoMIPs.yml"
  shell:"""
    bwa mem -L 80 -M -C {input.fasta} <( samtools view -F 513 {input.bam} | awk 'BEGIN{{ OFS="\\n"; FS="\\t" }}{{ print "@"$1"\\t"$12"\\t"$13"\\t"$14,$10,"+",$11 }}' ) | samtools view -u - | samtools sort - | scripts/pipeline2.0/TrimMIParms.py -d {input.design} -p | samtools reheader new_header.sam - | samtools sort -o {output} -
    """

rule indexing:
  input: expand("output/{{dataset}}/mapping/aligned/{plate}.bam",plate=loadSamples)
  output: expand("output/{{dataset}}/mapping/aligned/{plate}.bai",plate=loadSamples)
  conda: "hemoMIPs.yml"
  shell: "samtools index {input} {output}"

rule samplesexcheck:
  input: expand("output/{{dataset}}/mapping/aligned/{plate}.bam", plate=loadSamples)
  output:"output/{dataset}/samples_sex_check.txt"
  conda: "hemoMIPs.yml"
  shell: "( for i in {input}; do echo $( basename $i ) $(samtools view $i Y | wc -l) $(samtools view -F u $i | wc -l); done )> {output}"

rule inversionmips:
  input:bam="output/{dataset}/mapping/sample.bam",
        sam="input/{dataset}/new_header.sam",
        inv=expand("{reference}/hemomips_inv_ref.fa",reference=config["references"]["inv"])
  output:expand("output/{{dataset}}/mapping/inversion_mips/{plate}.bam", plate=loadSamples)
  conda: "hemoMIPs.yml"
  shell: """( grep "@RG" {input.sam}; bwa mem -M -L 80 -C {input.inv} <( samtools view -r ${{1}} -F 1 {input.bam} | awk 'BEGIN{{ OFS="\\n" }}{{ if (length($10) >= 75) {{ print "@"$1" "$12"\\t"$13"\\t"$14,$10,"+",$11 }} }}' ) | awk '{{ if (($0 ~ /^@/) || ($3 ~ /^inv/)) print }}';   bwa mem -M -L 80 -p -C {input.inv} <( samtools view -r ${{1}} -f 1 {input.bam} | awk 'BEGIN{{ OFS="\\n" }}{{ print "@"$1" "$12"\\t"$13"\\t"$14,$10,"+",$11 }}' ) | awk '{{ if (($0 !~ /^@/) && ($3 ~ /^inv/)) print }}' ) | samtools view -b -F 768 - | samtools sort -O bam -o {output} - """

rule inversionsum:
  input:lst="input/{dataset}/sample_index.lst", bam=expand("output/{{dataset}}/mapping/inversion_mips/{plate}.bam", plate=loadSamples)
  output:"output/{dataset}/inversion_mips/inversion_summary_counts.txt"
  conda: "hemoMIPs.yml"
  shell: """
    ( for i in $(tail -n +2 {input.lst} | cut -f 2 ); do echo $i $( ( samtools view -F 513 {input.bam} | awk 'BEGIN{{ FS="\\t" }}{{ split($12,a,":"); if (($6 !~ /S/) && (a[1] == "NM") && (a[3] <= 10)) {{ print $3 }} }}'; samtools view -f 2 -F 512 {input.bam} | awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ split($12,a,":"); if (($6 !~ /S/) && (a[1] == "NM") && (a[3] <= 10)) {{ print $1,$3 }} }}' | sort | uniq -c | awk '{{ if ($1 == 2) print $3 }}' ) | sort | uniq -c | awk '{{ print $1":"$2 }}' ); done )> {output}
    """



rule gatk3_align:
  input:bam=expand("output/{{dataset}}/mapping/aligned/{plate}.bam", plate=loadSamples)
  output:"output/{dataset}/mapping/gatk3/realign_all_samples.bam"
  shell: "java -Xmx8G -jar tools/GATK322/GenomeAnalysisTK.jar -T IndelRealigner -R {fa} -DBQ 3 -filterNoBases -maxReads 1500000 -maxInMemory 1500000 -targetIntervals targets_split.intervals $(ls -1 {input.bam} | xargs -n 1 echo -I ) -o {output} -dt BY_SAMPLE -dcov 500"


rule MIPstats:
  input:"output/{dataset}/mapping/gatk3/realign_all_samples.bam"
  output:"output/{dataset}/mapping/gatk3/realign_all_samples.MIPstats.tsv"
  shell:"scripts/pipeline2.0/MIPstats.py {input} -o {output}"


rule gatk3_genotyping:
  input:bam="output/{dataset}/mapping/gatk3/realign_all_samples.bam",
        targets="input/{dataset}/targets.intervals",
        fasta=config["references"]["fasta"]
  output:"output/{dataset}/mapping/gatk3/realign_all_samples.all_sites.vcf.gz"
  shell: "java -Xmx6G -jar scripts/GATK3446/GenomeAnalysisTK.jar -T UnifiedGenotyper -R {input.fasta} -I {input.bam} -L {input.targets} -o >( bgzip -c > {output} ) -glm BOTH -rf BadCigar --max_alternate_alleles 15 --output_mode EMIT_ALL_SITES -dt NONE"


rule subsetting:
  input:"output/{dataset}/mapping/gatk3/realign_all_samples.all_sites.vcf.gz"
  output:"output/{dataset}/mapping/gatk3/realign_all_samples.vcf.gz"
  shell: """
    zcat {input} | awk 'BEGIN{{ FS="\\t" }}{{ if ($1 ~ /^#/) {{ print }} else {{ if ($5 != ".") print }} }}' | bgzip -c > {output}"""


rule tabixing:
  input:"output/{dataset}/mapping/gatk3/realign_all_samples.vcf.gz"
  output:"output/{dataset}/mapping/gatk3/realign_all_samples.vcf.idx"
  shell: "tabix -p vcf {input} > {output}"

rule tabixing2:
  input:"output/{dataset}/mapping/gatk3/realign_all_samples.all_sites.vcf.gz"
  output:"output/{dataset}/mapping/gatk3/realign_all_samples.all_sites.vcf.idx"
  shell: "tabix -p vcf {input} > {output}"

rule VEP:
  input:"output/{dataset}/mapping/gatk3/realign_all_samples.vcf.gz"
  output:"output/{dataset}/mapping/gatk3/realign_all_samples.vep.tsv.gz"
  shell: "scripts/processing/VCF2VEP.sh {input} {output}"


rule checkPileUp:
  input: bam="output/{dataset}/mapping/gatk3/realign_all_samples.bam",
         sites="output/{dataset}/mapping/gatk3/realign_all_samples.vcf.gz"
  output:"output/{dataset}/mapping/gatk3/realign_all_samples.indel_check.txt"
  shell: "scripts/processing/checkPileUpAtInDels.py -b {input.bam} -s {input.sites} -o {output}"


rule summaryreport:
  input:vcf="output/{dataset}/mapping/gatk3/realign_all_samples.all_sites.vcf.gz",
        vep="output/{dataset}/mapping/gatk3/realign_all_samples.vep.tsv.gz",
        inv="output/{dataset}/inversion_mips/inversion_summary_counts.txt",
        sex="output/{dataset}/samples_sex_check.txt",
        target="input/target_coords.bed",
        mips="output/{dataset}/mapping/gatk3/realign_all_samples.MIPstats.tsv",
        indel="output/{dataset}/mapping/gatk3/realign_all_samples.indel_check.txt",
        hemomips="input/hemomips_design.txt",
        tg=config["references"]["annotation"]
  output:"output/{dataset}/gatk3/report/summary.html"
  shell: "scripts/processing/summary_report.py --vcf {input.vcf} --vep {input.vep} --inversions {input.inv} --sample_sex {input.sex} --target {input.target} --mipstats {input.mips} --indelCheck {input.indel} --design {input.hemomips} --TG {input.tg}"


rule gatk4_HTcaller:
  input:bamin=expand("output/{{dataset}}/mapping/aligned/{plate}.bam",plate=loadSamples),
        targets="input/{dataset}/targets.intervals",
        fasta=config["references"]["fasta"]
  output:bamout="output/{dataset}/mapping/gatk4/realign_all_samples.bam",
        vcf="output/{dataset}/mapping/gatk4/bam.vcf"
  conda:"python3gatk4.yml"
  shell: "gatk HaplotypeCaller -R {input.fasta} -L {input.targets} $(ls -1 {input.bamin} | xargs -n 1 echo -I ) --output-mode EMIT_ALL_SITES -bamout {output.bamout} -O {output.vcf} --disable-optimizations"

rule gatk4_gvcfs:
  input:lst="input/{dataset}/sample_index.lst",
        fasta=config["references"]["fasta"],
        targets="input/{dataset}/targets.intervals",
        bam=expand("output/{{dataset}}/mapping/aligned/{plate}.bam",plate=loadSamples)
  output:expand("output/{{dataset}}/mapping/gatk4/gvcf/{plate}.g.vcf",plate=loadSamples)
  conda:"python3gatk4.yml"
  shell:"""
    for i in $(less {input.lst} | cut -f 2 ); do (gatk HaplotypeCaller --max-reads-per-alignment-start 0 --disable-optimizations --sample-name ${{i}} -ERC GVCF -R {input.fasta} -I {input.bam} -L {input.targets} -O {output} ); done
    """

rule gatk4_combine:
  input:fasta=config["references"]["fasta"],
        gvcf=expand("output/{{dataset}}/mapping/gatk4/gvcf/{plate}.g.vcf",plate=loadSamples)
  output:bamout="output/{dataset}/mapping/gatk4/realign_all_samples.all_sites.vcf"
  conda:"python3gatk4.yml"
  shell:"gatk CombineGVCFs --break-bands-at-multiples-of 1 -R {input.fasta} $(ls -1 {input.gvcf} | xargs -n 1 echo -V ) -O {output}"
 

rule gatk4_genotype:
  input:fasta=config["references"]["fasta"],
        vcf="output/{dataset}/mapping/gatk4/realign_all_samples.all_sites.vcf"
  output:"output/{dataset}/mapping/gatk4/realign_all_samples.vcf"
  conda:"python3gatk4.yml"
  shell:"gatk GenotypeGVCFs -R {input.fasta} -V {input.vcf} -O {output}"

rule gatk4_gzip1:
  input:"output/{dataset}/mapping/gatk4/realign_all_samples.vcf"
  output:"output/{dataset}/mapping/gatk4/realign_all_samples.vcf.gz"
  shell: "bgzip {input} > {output}"

rule gatk4_gzip2:
  input:"output/{dataset}/mapping/gatk4/realign_all_samples.all_sites.vcf"
  output:"output/{dataset}/mapping/gatk4/realign_all_samples.all_sites.vcf.gz"
  shell: "bgzip {input} > {output}"

rule gatk4_MIPstats:
  input:"output/{dataset}/mapping/gatk4/realign_all_samples.bam"
  output:"output/{dataset}/mapping/gatk4/realign_all_samples.MIPstats.tsv"
  conda: "hemoMIPs.yml"
  shell: "/scripts/pipeline2.0/MIPstats.py {input} -o {output}"

rule VEP_gatk4:
  input:"output/{dataset}/mapping/gatk4/realign_all_samples.vcf.gz"
  output:"output/{dataset}/mapping/gatk4/realign_all_samples.vep.tsv.gz"
  conda: "hemoMIPs.yml"
  shell: "scripts/processing/VCF2VEP.sh {input} {output}"


rule checkPileUp_gatk4:
  input: bam="output/{dataset}/mapping/gatk4/realign_all_samples.bam",
         sites="output/{dataset}/mapping/gatk4/realign_all_samples.vcf.gz"
  output:"output/{dataset}/mapping/gatk4/realign_all_samples.indel_check.txt"
  conda: "hemoMIPs.yml"
  shell: "scripts/processing/checkPileUpAtInDels.py -b {input.bam} -s {input.sites} -o {output}"


rule summaryreport_gatk4:
  input:vcf="output/{dataset}/mapping/gatk4/realign_all_samples.all_sites.vcf.gz",
        vep="output/{dataset}/mapping/gatk4/realign_all_samples.vep.tsv.gz",
        inv="output/{dataset}/inversion_mips/inversion_summary_counts.txt",
        sex="output/{dataset}/samples_sex_check.txt",
        target="input/target_coords.bed",
        mips="output/{dataset}/mapping/gatk4/realign_all_samples.MIPstats.tsv",
        indel="output/{dataset}/mapping/gatk4/realign_all_samples.indel_check.txt",
        hemomips="input/hemomips_design.txt",
        tg=config["references"]["annotation"]
  output:"output/{dataset}/gatk4/report/summary.html"
  conda: "hemoMIPs.yml"
  shell: "scripts/processing/summary_report.py --vcf {input.vcf} --vep {input.vep} --inversions {input.inv} --sample_sex {input.sex} --target {input.target} --mipstats {input.mips} --indelCheck {input.indel} --design {input.hemomips} --TG {input.tg}"
