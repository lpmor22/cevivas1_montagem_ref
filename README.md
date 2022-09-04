# 1º Curso de Bioinformatica do CeVIVAS

## Roteiro de prática: tratamento das reads, montagem por referência e chamada de SNPs

### Para copiar e colar dados de/para o terminal Linux

* Copiar informação do terminal:
  * Selecionar com o mouse a informação a ser copiada
  * Copiar o dado com `Ctrl + Shift + C`
  * Colar no editor de texto de sua preferência com `Ctrl + V`

* Copiar para o terminal:
  * Copiar o dado do editor de texto de sua preferência com `Ctrl + C`
  * Clicar no terminal
  * Utilizar `Ctrl + Shift + V` para colar a informação no terminal

### Conectar à máquina virtual de cursos do Núcleo de Bioinformática e Biologia Computacional (NBBC) do Instituto Butantan

Para esta prática iremos acessar a máquina virtual via SSH (Secure Shell). É necessário trocar `usuario` pelo entregue pelo NCCB:

    ssh -p 2202 usuario@200.136.54.100

### Realizar análises de sequências utilizando o ViralFlow (via IGM_SARSCOV2)

Checar se o igm_sarscov2 está instalado:

    igm_sarscov2 -h

Checar se os ambientes condas estão instalados:

    conda env list

Rodar o pipeline ViralFlow (via IGM_SARSCOV2):

    igm_sarscov2 -w 1 -t 8 -p ARTIC_V3 -i cevivas/

## Como é dentro do pipeline a montagem?

Vamos testar com uma amostra SRR15365366 (https://www.ncbi.nlm.nih.gov/sra/SRR15365366). Esta amostra foi sequenciada utilizando Illumina NovaSeq 6000, paired-end, protocolo Illumina COVIDSeq com primers ARTIC versão 3.

- Ativar ambiente conda com as dependências necessárias para montagem:

    source activate igm-sars2_assembly

- Filtrar as leituras com qualidade PHRED >= 20 e manter leituras até o mínimo de 75 bp:

    fastp --cut_front --cut_tail --qualified_quality_phred 20 -l 75 -f 0 -t 0 -F 0 -T 0 \
      -i SRR15365366_1.fastq.gz -I SRR15365366_2.fastq.gz --adapter_fasta ARTIC_V3.primers.fasta \
      -o SRR15365366.trimado.R1.fastq.gz -O SRR15365366.trimado.R2.fastq.gz \
      -h SRR15365366.relatorio_fastp.html -j SRR15365366.relatorio_fastp.json

- Criar lista de index da sequência referência MN908947.3:

    bwa index MN908947.3.fasta

- XXX

    bwa mem MN908947.3.fasta SRR15365366.trimado.R1.fastq.gz SRR15365366.trimado.R2.fastq.gz -o SRR15365366.bam

- XXX

    samtools sort -o SRR15365366.sorted.bam SRR15365366.bam

- XXX

    samtools index SRR15365366.sorted.bam

- XXX

    samtools mpileup -d 50000 --reference MN908947.3.fasta -a -B SRR15365366.sorted.bam | \
      ivar variants -p SRR15365366 -q 30 -t 0.05

- XXX

    samtools mpileup -d 50000 --reference MN908947.3.fasta -a -B SRR15365366.sorted.bam | \
      ivar consensus -p SRR15365366 -q 30 -t 0 -m 10 -n N

- XXX

    sed -i -e 's/>.*/>'SRR15365366'/g' SRR15365366.fa

- XXX

    mafft --quiet --auto --keeplength --inputorder --6merpair --leavegappyregion \
      --addfragments SRR15365366.fa MN908947.3.fasta | seqkit grep -ip "MN908947.3"

- XXX

    conda deactivate

- XXX

    source activate igm-sars2_summary


- XXX

    echo "sample_id#num_total_reads#num_mapp_reads#avg_depth#depth_10x#depth_100x#depth_1000x#ref_cov#ncount#ncount_perc#pango_ver#pango_learn_ver#pango_lin#nextclade_ver#nextclade_lin" | tr '#' '\t' > sumario_montagem.txt

- XXX

    echo -n " SRR15365366""#" | tr '#' '\t' >> sumario_montagem.txt

- XXX

    samtools view -c SRR15365366.sorted.bam | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

- XXX

    samtools view -c -h -F 4 SRR15365366.sorted.bam | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

- XXX

    samtools depth SRR15365366.sorted.bam | awk '{sum+=$3} END {print sum/NR}' | | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

- XXX

    paste <(samtools depth SRR15365366.sorted.bam | awk '{if ($3 > '"10"') {print $0}}' | wc -l) <(fastalength MN908947.3.fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

- XXX

    paste <(samtools depth SRR15365366.sorted.bam | awk '{if ($3 > '"100"') {print $0}}' | wc -l) <(fastalength MN908947.3.fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

- XXX

    paste <(samtools depth SRR15365366.sorted.bam | awk '{if ($3 > '"1000"') {print $0}}' | wc -l) <(fastalength MN908947.3.fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

- XXX

    paste <(fastalength MN908947.3.fasta | awk '{print $1}') <(seqtk comp SRR15365366.fa | awk -F"\t" '{print $9}') | awk -F"\t" '{printf("%0.2f\n", ($1-$2)/$1*100)}') | tr '#' '\t' >> sumario_montagem.txt

- XXX

    seqtk comp SRR15365366.fa | awk -F"\t" '{print $9}') | tr '#' '\t' >> sumario_montagem.txt

- XXX

    paste <(seqtk comp SRR15365366.fa | awk -F"\t" '{print $9}') <(fastalength MN908947.3.fasta | awk '{print $1}')| awk -F"\t" '{printf("%0.2f\n", ($1/$2)*100)}') | tr '#' '\t' >> sumario_montagem.txt

- XXX

    pangolin SRR15365366.fa --outfile SRR15365366.pangolin.csv

- XXX

    cat SRR15365366.pangolin.csv | sed -n 2p | awk -F, '{print $10"\t"$9"\t"$2}' | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt


- XXX

    nextclade dataset get --name 'sars-cov-2' --output-dir nextclade_sc2

- XXX

    nextclade run --input-dataset nextclade_sc2 --output-tsv=SRR15365366.nextclade.tsv SRR15365366.fa

- XXX

    nextclade --version | awk '{print $2}' | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

- XXX

    cat SRR15365366.nextclade.tsv | sed -n 2p | awk -F"\t" '{print $2}' | awk '{printf $0}' >> sumario_montagem.txt

- XXX

    conda deactivate

- XXX

    fastcov.py -l SRR15365366.sorted.bam -o SRR15365366.coverage.pdf

- XXX

    conda deactivate
