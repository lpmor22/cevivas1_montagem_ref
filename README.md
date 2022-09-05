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

Vamos testar o pipeline com amostras do projeto `PRJNA686074` (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA686074).
Estas amostras foram sequenciadas utilizando Illumina NovaSeq 6000, paired-end, protocolo Illumina COVIDSeq com primers ARTIC versão 3.

    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/066/SRR15365366/SRR15365366_1.fastq.gz -o SRR15365366_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/066/SRR15365366/SRR15365366_2.fastq.gz -o SRR15365366_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/086/SRR15428386/SRR15428386_1.fastq.gz -o SRR15428386_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/086/SRR15428386/SRR15428386_2.fastq.gz -o SRR15428386_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/045/SRR15428345/SRR15428345_1.fastq.gz -o SRR15428345_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/045/SRR15428345/SRR15428345_2.fastq.gz -o SRR15428345_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/046/SRR15428346/SRR15428346_1.fastq.gz -o SRR15428346_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/046/SRR15428346/SRR15428346_2.fastq.gz -o SRR15428346_2.fastq.gz

Rodar o pipeline ViralFlow (via IGM_SARSCOV2):

    igm_sarscov2 -w 1 -t 8 -p ARTIC_V3 -i cevivas/

## Como é dentro do pipeline a montagem?

Vamos testar com a amostra `SRR15365366` do mesmo projeto `PRJNA686074` (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA686074).

Ativar ambiente conda com as dependências necessárias para montagem:

    source activate igm-sars2_assembly

Filtrar as leituras com qualidade PHRED >=20 e manter leituras até o mínimo de 75 bp:

    fastp --cut_front --cut_tail --qualified_quality_phred 20 -l 75 -f 0 -t 0 -F 0 -T 0 \
      -i SRR15365366_1.fastq.gz -I SRR15365366_2.fastq.gz --adapter_fasta ARTIC_V3.primers.fasta \
      -o SRR15365366.trimado.R1.fastq.gz -O SRR15365366.trimado.R2.fastq.gz \
      -h SRR15365366.relatorio_fastp.html -j SRR15365366.relatorio_fastp.json

Criar lista de index da sequência referência `MN908947.3`:

    bwa index MN908947.3.fasta

Mapear contra o genoma referência `MN908947.3`:

    bwa mem MN908947.3.fasta SRR15365366.trimado.R1.fastq.gz SRR15365366.trimado.R2.fastq.gz -o SRR15365366.bam

Organizar o mapeamento de acordo com as coordenadas do genoma referência `MN908947.3`:

    samtools sort -o SRR15365366.sorted.bam SRR15365366.bam

Criar lista de index do mapeamento:

    samtools index SRR15365366.sorted.bam

Criar arquivo contendo os SNPs em relação ao genoma referência `MN908947.3`:

    samtools mpileup -d 50000 --reference MN908947.3.fasta -a -B SRR15365366.sorted.bam | \
      ivar variants -p SRR15365366 -q 30 -t 0.05

Criar arquivo fasta com o genoma consenso - chamada de variante com profundidade 10x e inclusão de N em regiões <= 10x:

    samtools mpileup -d 50000 --reference MN908947.3.fasta -a -B SRR15365366.sorted.bam | \
      ivar consensus -p SRR15365366 -q 30 -t 0 -m 10 -n N

Mudar cabeçalho da sequência para o ID da amostra:

    sed -i -e 's/>.*/>'SRR15365366'/g' SRR15365366.fa

Alinhar contra o genoma referência `MN908947.3` para organizar os frames de leitura:

    mafft --quiet --auto --keeplength --inputorder --6merpair --leavegappyregion \
      --addfragments SRR15365366.fa MN908947.3.fasta | seqkit grep -ip "MN908947.3"

Desativar o ambiente conda:

    conda deactivate

Ativar ambiente conda com as dependências necessárias para obter métricas de montagem:

    source activate igm-sars2_summary

Criar arquivo de texto pra colocar as métricas de montagem:

    echo "id_amostra#num_leituras_total#num_leituras_mapeadas#medias_profundidade#profundidade_10x#profundidade_100x#profundidade_1000x#cobertura_referencia#contagem_n#contagem_n_porcent#pango_versao#pango_database#pango_linhagem#nextclade_versao#nextclade_clado" | tr '#' '\t' > sumario_montagem.txt

Adicionar ao arquivo do sumário o ID da amostra:

    echo -n "SRR15365366""#" | tr '#' '\t' >> sumario_montagem.txt

Obter número total de leituras geradas:

    samtools view -c SRR15365366.sorted.bam | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

Obter o número de leituras mapeadas contra o genoma referência `MN908947.3`:

    samtools view -c -h -F 4 SRR15365366.sorted.bam | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

Obter a média da profundidade:

    samtools depth SRR15365366.sorted.bam | awk '{sum+=$3} END {print sum/NR}' | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

Obter a porcentagem de profundidade 10x em relação à referência `MN908947.3`:

    paste <(samtools depth SRR15365366.sorted.bam | awk '{if ($3 > '"10"') {print $0}}' | wc -l) <(fastalength MN908947.3.fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

Obter a porcentagem de profundidade 100x em relação à referência `MN908947.3`:

    paste <(samtools depth SRR15365366.sorted.bam | awk '{if ($3 > '"100"') {print $0}}' | wc -l) <(fastalength MN908947.3.fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

Obter a porcentagem de profundidade 1000x em relação à referência `MN908947.3`:

    paste <(samtools depth SRR15365366.sorted.bam | awk '{if ($3 > '"1000"') {print $0}}' | wc -l) <(fastalength MN908947.3.fasta | awk '{print $1}') | awk -F"\t" '{printf("%0.2f\n", $1/$2*100)}' | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

Obter a cobertura em porcentagem em relação ao genoma referência `MN908947.3`:

    paste <(fastalength MN908947.3.fasta | awk '{print $1}') <(seqtk comp SRR15365366.fa | awk -F"\t" '{print $9}') | awk -F"\t" '{printf("%0.2f\n", ($1-$2)/$1*100)}') | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

Obter o número absoluto de Ns na sequência consenso:

    seqtk comp SRR15365366.fa | awk -F"\t" '{print $9}') | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

Obter a porcentagem de Ns na sequência consenso:

    paste <(seqtk comp SRR15365366.fa | awk -F"\t" '{print $9}') <(fastalength MN908947.3.fasta | awk '{print $1}')| awk -F"\t" '{printf("%0.2f\n", ($1/$2)*100)}') | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

Rodar classificação do pangolin:

    pangolin SRR15365366.fa --outfile SRR15365366.pangolin.csv

Obter as versões do pangolin, do banco de dados utilizado e linhagem do SARS-CoV-2:

    cat SRR15365366.pangolin.csv | sed -n 2p | awk -F, '{print $10"\t"$9"\t"$2}' | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

Criar dataset do nextclade para classificação:

    nextclade dataset get --name 'sars-cov-2' --output-dir nextclade_sc2

Rodar classificação do nextclade:

    nextclade run --input-dataset nextclade_sc2 --output-tsv=SRR15365366.nextclade.tsv SRR15365366.fa

Obter a versão do nextclade:

    nextclade --version | awk '{print $2}' | awk '{printf $0"#"}' | tr '#' '\t' >> sumario_montagem.txt

Obter o clado do SARS-CoV-2:

    cat SRR15365366.nextclade.tsv | sed -n 2p | awk -F"\t" '{print $2}' | awk '{printf $0}' >> sumario_montagem.txt

Obter o plot de cobertura e profundidade:

    fastcov.py -l SRR15365366.sorted.bam -o SRR15365366.coverage.pdf

Desativar o ambiente conda:

    conda deactivate
