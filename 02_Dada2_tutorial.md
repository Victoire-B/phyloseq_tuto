Dada2 tutorial
================

  - [Appel des séquences](#appel-des-séquences)
      - [Apprentissage des erreurs](#apprentissage-des-erreurs)
      - [Visualisation des modèles d’ereur du
        forward](#visualisation-des-modèles-dereur-du-forward)
      - [Exemple d’interférences](#exemple-dinterférences)
  - [Aligner les R1 et R2 en un
    contig](#aligner-les-r1-et-r2-en-un-contig)
  - [Construction d’une table
    d’observation](#construction-dune-table-dobservation)
  - [Inspection des longueurs de
    séquences](#inspection-des-longueurs-de-séquences)
  - [Chimères](#chimères)
      - [Enlever les chimères par méthode
        consensus](#enlever-les-chimères-par-méthode-consensus)
      - [Faire le ratio](#faire-le-ratio)
  - [Construction d’une table](#construction-dune-table)
      - [Assignation taxonomique n°2 Silva species
        assignement](#assignation-taxonomique-n2-silva-species-assignement)
      - [Evaluer la précision de l’assignation
        taxonomique](#evaluer-la-précision-de-lassignation-taxonomique)
  - [Conclusion](#conclusion)

# Appel des séquences

``` r
library("dada2")
```

    ## Loading required package: Rcpp

    ## Warning: multiple methods tables found for 'which'

On a appellé le package dada2

Les données ont été téléchargées On definit la variable Path pour
changer d’endroit dans l’arborescense du fichier une fois dézipé
MiSeq\_SOP :

``` r
path <- "~/MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

On voit les fichiers fastq

On crée des variables pour l’analyse du plot, les Read 1 vont dans la
variable fnFs et les Read2 dans la variable fnRs

  - la foction sample.names : extrait les noms des échantillons, tout en
    supposant que les noms de fichiers ont un format : NOM DE
    L’ÉCHANTILLON\_XXX.fastq

<!-- end list -->

``` r
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

On a récupéré seulement les read issu d’Illumina –\> on va évaluer la
qualité de chaque read On affiche les profils de qualité des Read1 (R1)
(= Forward) en avant du premier et du deuxième read issus de la variable
crée précédémment FnFs

``` r
plotQualityProfile(fnFs[1:2])
```

![](02_Dada2_tutorial_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

On affiche les profil de qualité des Read2 (R2) (= Reverse) en arrière
du premier et du deuxième read issus de la variable crée précédémment
FnRs

``` r
plotQualityProfile(fnRs[1:2])
```

![](02_Dada2_tutorial_files/figure-gfm/unnamed-chunk-5-1.png)<!-- --> \#
Filtration

On va filtrer les variables contenant les read 1 et 2. À partir de ces
filtrages on attribue des variable et des noms associées filtFs : pour
les read1 filtRs : pour les read2

``` r
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

Les paramètres de filtrage standards : -maxN=0 (Dada2 ne nécessite pas
de Ns), -truncQ=2 : on va tronquer pour les reads 1 a 240 et 160 pour
les reads 2 -rm.phix=TRUE -maxEE=2 -\> le nombre maximum d’“erreurs
attendues” autorisées dans une lecture, ce qui est un meilleur filtre
que la simple moyenne des scores de qualité

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

On obtient le nombre de bases (les paires de base) qui sont gardées. La
filtration est effectuée ici que pour les R1. Pour le reads F3D30 S188
R1 001 :7793 nt pour le read avant filtration et 7113 nt après
filtration. 680 nt ont été enlevé : les primers ont été supprimés et les
nucléotides pas bons.

## Apprentissage des erreurs

Il est possible d’avoir des erreurs, avec Dada2 on inspecte les
séquences. On utilise un modèle d’erreur paramétrique err pour les R1
et R2

Le But : Le modèle d’erreur de DADA2 permet identifier les positions
avec une forte probabilité d’erreur et par la suite changer avec la base
la plus probable. Cela veut dire celle quiest présente dans la séquence
la plus abondante

On crée les variables : -errF : recoit le modèle d’erreur paramétrique
par la fonction LearnErrors pour les R1 et R2 filtrés

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

## Visualisation des modèles d’ereur du forward

En abscisse on a la probabilité des mutations et en ordonnée le q score
Q30 : la probabilité que la base trouver sois la bonne

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](02_Dada2_tutorial_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Chaque mutation est possible, les taux d’erreur sont indiqués. - points
: les taux d’erreur observés pour chaque score de qualité. - ligne
noire: taux d’erreur estimé après convergence de l’algorithme
d’apprentissage. - ligne rouge: taux d’erreur attendu selon le
Q-score.

Les taux d’erreur estimés correspondent aux taux observés, et les taux
d’erreur diminuent avec l’augmentation de la qualité.

## Exemple d’interférences

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

Calcul pour les séquences composées de séquences uniques Dans notre
exemple avec les donnée de MiSeq\_SOP, l’échantillon 1 a 7113 read dont
1979 read unique les autres peuvent être retrouvées plusieus fois. Dans
les information supplémentaires on a les redondances, le clustering …

On fait pareil pour les reverses (Rs)

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

Ici, on inspecte la premiere étagère du premeir placard de dada2

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

On inspecte la deuxième étagère du premier placard de dada2 (deuxième
échantillon)

``` r
dadaFs[[2]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 113 sequence variants were inferred from 1639 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

Le filtrage des reads est fini.

# Aligner les R1 et R2 en un contig

Le but est de créer des contigs à partir des R1(forwards) et des R2
(reverse) Ici c’est posssible car c’est une amplification de la région
V4 de l’ARN 16S il y a donc un overlap avec read1 et read2

  - verbose : montrer les étapes avec du texte pendant qu’elles sont
    réalisées
  - objet merge: liste de data.frames de chaque échantillon. Chaque
    data.frame contient la séquence fusionnée, son abondance, et les
    indices des variantes de la séquence avant et arrière fusionnées.
  - mergePairs: suppression des read apparies qui se chevauchaient -\>
    réduction du bruit
  - head merge : regarder la première ligne

<!-- end list -->

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 6551 paired-reads (in 106 unique pairings) successfully merged out of 6907 (in 199 pairings) input.

    ## 5025 paired-reads (in 100 unique pairings) successfully merged out of 5188 (in 156 pairings) input.

    ## 4973 paired-reads (in 80 unique pairings) successfully merged out of 5268 (in 166 pairings) input.

    ## 2595 paired-reads (in 52 unique pairings) successfully merged out of 2756 (in 109 pairings) input.

    ## 2553 paired-reads (in 60 unique pairings) successfully merged out of 2785 (in 119 pairings) input.

    ## 3622 paired-reads (in 53 unique pairings) successfully merged out of 4103 (in 157 pairings) input.

    ## 6079 paired-reads (in 81 unique pairings) successfully merged out of 6515 (in 198 pairings) input.

    ## 3961 paired-reads (in 90 unique pairings) successfully merged out of 4384 (in 188 pairings) input.

    ## 14231 paired-reads (in 143 unique pairings) successfully merged out of 15358 (in 351 pairings) input.

    ## 10526 paired-reads (in 120 unique pairings) successfully merged out of 11166 (in 279 pairings) input.

    ## 11156 paired-reads (in 137 unique pairings) successfully merged out of 11799 (in 298 pairings) input.

    ## 4329 paired-reads (in 84 unique pairings) successfully merged out of 4788 (in 180 pairings) input.

    ## 17431 paired-reads (in 153 unique pairings) successfully merged out of 17812 (in 272 pairings) input.

    ## 5850 paired-reads (in 81 unique pairings) successfully merged out of 6095 (in 159 pairings) input.

    ## 3716 paired-reads (in 86 unique pairings) successfully merged out of 3894 (in 147 pairings) input.

    ## 6865 paired-reads (in 99 unique pairings) successfully merged out of 7193 (in 187 pairings) input.

    ## 4430 paired-reads (in 67 unique pairings) successfully merged out of 4605 (in 127 pairings) input.

    ## 4574 paired-reads (in 100 unique pairings) successfully merged out of 4736 (in 172 pairings) input.

    ## 6094 paired-reads (in 109 unique pairings) successfully merged out of 6314 (in 172 pairings) input.

    ## 4269 paired-reads (in 20 unique pairings) successfully merged out of 4281 (in 28 pairings) input.

``` r
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                       sequence
    ## 1 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACAGG
    ## 2 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAGG
    ## 3 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACCCCGGCCTGCCGTTGAAACTGGCGGCCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ## 4 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTTTTAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCCGGCCGTTGAAACTGGGGGCCTTGAGTGGGCGAGAAGAAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCCTTCCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCGAACAGG
    ## 5 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGACTCTCAAGTCAGCGGTCAAATCGCGGGGCTCAACCCCGTTCCGCCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAAGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTCATGCACGAAAGCGTGGGTATCGAACAGG
    ## 6 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1       579       1       1    148         0      0      1   TRUE
    ## 2       470       2       2    148         0      0      2   TRUE
    ## 3       449       3       4    148         0      0      1   TRUE
    ## 4       430       4       3    148         0      0      2   TRUE
    ## 5       345       5       6    148         0      0      1   TRUE
    ## 6       282       6       5    148         0      0      2   TRUE

# Construction d’une table d’observation

Table de variant de séquences d’amplicon (ASV), une version à plus haute
résolution de la table OTU produite par les méthodes traditionnelles.

On va partir de l’abondance, on crée la table à partir de notre objet
merger

  - objet seqtab : une table avec en ligne le nombre échantillon, en
    colonne les séquences elle-même à l’intérieur ou on observe la
    séquence dans l’échantillon
  - dim : récupérer ou définir la dimension d’un objet.

<!-- end list -->

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]  20 293

20 entrées = contigs et 293 colones = la séquence nucléotidique

Une matrice avec des lignes correspondant aux échantillons et des
colonnes correspondant aux variantes de séquences. Ce tableau contient
293 ASV, et les longueurs de nos séquences fusionnées se situent toutes
dans la plage prévue pour cet amplicon V4.

# Inspection des longueurs de séquences

``` r
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  88 196   6   2

# Chimères

Pendant l’amplification pendant PCR, on amplifie avec les primers Avec
le primer F(pour les forwards), on crée la séquence complémentaire, mais
pour x ou y, l’élongation s’arrête Comment se créer les chimères : On a
le 16S qui restera intact et un autre fragment simple brin plus court.
Lors du cycle suivant en PCR, on va avoir un des fragment 16S (rouge par
exemple) qui vont s’hydrider sur un autre fragment 16S comme un primer
et donc continuer élongation (verte par exemple) pour donner une
séquence hybride à la fin qui sera le début du rouge et a la fin du
vert

Ce phénomène est rare, mais lorsqu’il y a plein de séquence possible et
il faut les enlever du jeu de données Il n’est pas possible de les
détecter au niveau de la taille de la séquence. En revanche, on peut
regarder toutes les séquences rares dont le début correspond à une
séquence parent dans ce jeu donnée et la fin d’une autre séquence
parent

Appliquer à seqtab et transférer à une nouvelle variable seqtab.nochim

\-\> 1/5, les 293 uniques mais qui représente plus de séquence dans le
jeu de donnée

## Enlever les chimères par méthode consensus

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 61 bimeras out of 293 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]  20 232

61 chimères obtenu dans les 983 ASV

## Faire le ratio

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.964263

Il y a 0,96% de séquences chimériques dans notre jeu de données

# Construction d’une table

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##        input filtered denoisedF denoisedR merged nonchim
    ## F3D0    7793     7113      6996      6978   6551    6539
    ## F3D1    5869     5299      5227      5239   5025    5014
    ## F3D141  5958     5463      5339      5351   4973    4850
    ## F3D142  3183     2914      2799      2833   2595    2521
    ## F3D143  3178     2941      2822      2868   2553    2519
    ## F3D144  4827     4312      4146      4224   3622    3483

``` bash
cd $HOME
wget https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz
```

    ## --2020-12-03 15:14:07--  https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 137973851 (132M) [application/octet-stream]
    ## Saving to: ‘silva_nr99_v138_train_set.fa.gz.5’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 7.64M 17s
    ##     50K .......... .......... .......... .......... ..........  0% 10.2M 15s
    ##    100K .......... .......... .......... .......... ..........  0% 12.9M 13s
    ##    150K .......... .......... .......... .......... ..........  0% 46.9M 11s
    ##    200K .......... .......... .......... .......... ..........  0% 19.0M 10s
    ##    250K .......... .......... .......... .......... ..........  0% 42.4M 9s
    ##    300K .......... .......... .......... .......... ..........  0% 13.5M 9s
    ##    350K .......... .......... .......... .......... ..........  0%  114M 8s
    ##    400K .......... .......... .......... .......... ..........  0% 10.1M 9s
    ##    450K .......... .......... .......... .......... ..........  0% 55.0M 8s
    ##    500K .......... .......... .......... .......... ..........  0% 15.7M 8s
    ##    550K .......... .......... .......... .......... ..........  0% 75.2M 7s
    ##    600K .......... .......... .......... .......... ..........  0% 16.0M 7s
    ##    650K .......... .......... .......... .......... ..........  0%  105M 7s
    ##    700K .......... .......... .......... .......... ..........  0% 17.1M 7s
    ##    750K .......... .......... .......... .......... ..........  0% 76.6M 7s
    ##    800K .......... .......... .......... .......... ..........  0% 17.3M 7s
    ##    850K .......... .......... .......... .......... ..........  0% 97.1M 6s
    ##    900K .......... .......... .......... .......... ..........  0% 15.6M 7s
    ##    950K .......... .......... .......... .......... ..........  0% 94.2M 6s
    ##   1000K .......... .......... .......... .......... ..........  0% 13.6M 6s
    ##   1050K .......... .......... .......... .......... ..........  0% 88.2M 6s
    ##   1100K .......... .......... .......... .......... ..........  0% 14.9M 6s
    ##   1150K .......... .......... .......... .......... ..........  0% 15.3M 6s
    ##   1200K .......... .......... .......... .......... ..........  0% 72.9M 6s
    ##   1250K .......... .......... .......... .......... ..........  0% 75.6M 6s
    ##   1300K .......... .......... .......... .......... ..........  1% 10.6M 6s
    ##   1350K .......... .......... .......... .......... ..........  1% 96.7M 6s
    ##   1400K .......... .......... .......... .......... ..........  1% 17.9M 6s
    ##   1450K .......... .......... .......... .......... ..........  1% 29.6M 6s
    ##   1500K .......... .......... .......... .......... ..........  1% 22.6M 6s
    ##   1550K .......... .......... .......... .......... ..........  1% 37.9M 6s
    ##   1600K .......... .......... .......... .......... ..........  1% 19.8M 6s
    ##   1650K .......... .......... .......... .......... ..........  1% 92.3M 6s
    ##   1700K .......... .......... .......... .......... ..........  1% 54.4M 6s
    ##   1750K .......... .......... .......... .......... ..........  1% 15.9M 6s
    ##   1800K .......... .......... .......... .......... ..........  1% 14.6M 6s
    ##   1850K .......... .......... .......... .......... ..........  1% 84.3M 6s
    ##   1900K .......... .......... .......... .......... ..........  1% 20.2M 6s
    ##   1950K .......... .......... .......... .......... ..........  1% 30.1M 6s
    ##   2000K .......... .......... .......... .......... ..........  1% 36.4M 6s
    ##   2050K .......... .......... .......... .......... ..........  1% 20.2M 6s
    ##   2100K .......... .......... .......... .......... ..........  1% 14.1M 6s
    ##   2150K .......... .......... .......... .......... ..........  1% 68.4M 6s
    ##   2200K .......... .......... .......... .......... ..........  1% 93.8M 6s
    ##   2250K .......... .......... .......... .......... ..........  1% 14.0M 6s
    ##   2300K .......... .......... .......... .......... ..........  1% 44.1M 6s
    ##   2350K .......... .......... .......... .......... ..........  1% 24.2M 6s
    ##   2400K .......... .......... .......... .......... ..........  1% 61.9M 6s
    ##   2450K .......... .......... .......... .......... ..........  1% 14.6M 6s
    ##   2500K .......... .......... .......... .......... ..........  1% 92.0M 6s
    ##   2550K .......... .......... .......... .......... ..........  1% 13.3M 6s
    ##   2600K .......... .......... .......... .......... ..........  1% 62.5M 6s
    ##   2650K .......... .......... .......... .......... ..........  2% 20.3M 6s
    ##   2700K .......... .......... .......... .......... ..........  2% 62.3M 5s
    ##   2750K .......... .......... .......... .......... ..........  2% 14.1M 6s
    ##   2800K .......... .......... .......... .......... ..........  2% 69.8M 5s
    ##   2850K .......... .......... .......... .......... ..........  2%  111M 5s
    ##   2900K .......... .......... .......... .......... ..........  2% 18.7M 5s
    ##   2950K .......... .......... .......... .......... ..........  2% 82.8M 5s
    ##   3000K .......... .......... .......... .......... ..........  2% 19.2M 5s
    ##   3050K .......... .......... .......... .......... ..........  2% 78.4M 5s
    ##   3100K .......... .......... .......... .......... ..........  2% 22.3M 5s
    ##   3150K .......... .......... .......... .......... ..........  2% 49.1M 5s
    ##   3200K .......... .......... .......... .......... ..........  2% 84.5M 5s
    ##   3250K .......... .......... .......... .......... ..........  2% 13.5M 5s
    ##   3300K .......... .......... .......... .......... ..........  2% 86.1M 5s
    ##   3350K .......... .......... .......... .......... ..........  2% 13.8M 5s
    ##   3400K .......... .......... .......... .......... ..........  2% 94.0M 5s
    ##   3450K .......... .......... .......... .......... ..........  2% 14.3M 5s
    ##   3500K .......... .......... .......... .......... ..........  2% 68.4M 5s
    ##   3550K .......... .......... .......... .......... ..........  2%  118M 5s
    ##   3600K .......... .......... .......... .......... ..........  2% 14.2M 5s
    ##   3650K .......... .......... .......... .......... ..........  2% 97.2M 5s
    ##   3700K .......... .......... .......... .......... ..........  2% 21.0M 5s
    ##   3750K .......... .......... .......... .......... ..........  2% 79.4M 5s
    ##   3800K .......... .......... .......... .......... ..........  2% 22.6M 5s
    ##   3850K .......... .......... .......... .......... ..........  2% 26.9M 5s
    ##   3900K .......... .......... .......... .......... ..........  2% 99.0M 5s
    ##   3950K .......... .......... .......... .......... ..........  2% 12.7M 5s
    ##   4000K .......... .......... .......... .......... ..........  3% 84.5M 5s
    ##   4050K .......... .......... .......... .......... ..........  3%  117M 5s
    ##   4100K .......... .......... .......... .......... ..........  3% 13.8M 5s
    ##   4150K .......... .......... .......... .......... ..........  3%  115M 5s
    ##   4200K .......... .......... .......... .......... ..........  3% 14.4M 5s
    ##   4250K .......... .......... .......... .......... ..........  3% 69.0M 5s
    ##   4300K .......... .......... .......... .......... ..........  3%  109M 5s
    ##   4350K .......... .......... .......... .......... ..........  3% 15.7M 5s
    ##   4400K .......... .......... .......... .......... ..........  3% 82.7M 5s
    ##   4450K .......... .......... .......... .......... ..........  3% 14.2M 5s
    ##   4500K .......... .......... .......... .......... ..........  3% 62.1M 5s
    ##   4550K .......... .......... .......... .......... ..........  3%  102M 5s
    ##   4600K .......... .......... .......... .......... ..........  3% 16.7M 5s
    ##   4650K .......... .......... .......... .......... ..........  3% 54.5M 5s
    ##   4700K .......... .......... .......... .......... ..........  3%  102M 5s
    ##   4750K .......... .......... .......... .......... ..........  3% 12.2M 5s
    ##   4800K .......... .......... .......... .......... ..........  3% 87.3M 5s
    ##   4850K .......... .......... .......... .......... ..........  3% 12.9M 5s
    ##   4900K .......... .......... .......... .......... ..........  3%  106M 5s
    ##   4950K .......... .......... .......... .......... ..........  3%  114M 5s
    ##   5000K .......... .......... .......... .......... ..........  3% 14.4M 5s
    ##   5050K .......... .......... .......... .......... ..........  3%  107M 5s
    ##   5100K .......... .......... .......... .......... ..........  3% 14.4M 5s
    ##   5150K .......... .......... .......... .......... ..........  3% 89.4M 5s
    ##   5200K .......... .......... .......... .......... ..........  3% 98.8M 5s
    ##   5250K .......... .......... .......... .......... ..........  3% 14.0M 5s
    ##   5300K .......... .......... .......... .......... ..........  3% 98.5M 5s
    ##   5350K .......... .......... .......... .......... ..........  4% 18.3M 5s
    ##   5400K .......... .......... .......... .......... ..........  4% 73.5M 5s
    ##   5450K .......... .......... .......... .......... ..........  4%  109M 5s
    ##   5500K .......... .......... .......... .......... ..........  4% 17.2M 5s
    ##   5550K .......... .......... .......... .......... ..........  4% 73.3M 5s
    ##   5600K .......... .......... .......... .......... ..........  4% 17.4M 5s
    ##   5650K .......... .......... .......... .......... ..........  4% 97.9M 5s
    ##   5700K .......... .......... .......... .......... ..........  4% 67.6M 5s
    ##   5750K .......... .......... .......... .......... ..........  4% 13.4M 5s
    ##   5800K .......... .......... .......... .......... ..........  4%  100M 5s
    ##   5850K .......... .......... .......... .......... ..........  4% 17.0M 5s
    ##   5900K .......... .......... .......... .......... ..........  4% 50.5M 5s
    ##   5950K .......... .......... .......... .......... ..........  4% 69.5M 5s
    ##   6000K .......... .......... .......... .......... ..........  4% 20.2M 5s
    ##   6050K .......... .......... .......... .......... ..........  4% 70.7M 5s
    ##   6100K .......... .......... .......... .......... ..........  4% 98.5M 5s
    ##   6150K .......... .......... .......... .......... ..........  4% 24.3M 5s
    ##   6200K .......... .......... .......... .......... ..........  4% 35.4M 5s
    ##   6250K .......... .......... .......... .......... ..........  4% 25.7M 5s
    ##   6300K .......... .......... .......... .......... ..........  4% 29.4M 5s
    ##   6350K .......... .......... .......... .......... ..........  4%  121M 5s
    ##   6400K .......... .......... .......... .......... ..........  4% 31.9M 5s
    ##   6450K .......... .......... .......... .......... ..........  4% 33.0M 5s
    ##   6500K .......... .......... .......... .......... ..........  4% 31.8M 5s
    ##   6550K .......... .......... .......... .......... ..........  4% 22.5M 5s
    ##   6600K .......... .......... .......... .......... ..........  4% 66.2M 5s
    ##   6650K .......... .......... .......... .......... ..........  4% 76.9M 5s
    ##   6700K .......... .......... .......... .......... ..........  5% 19.2M 5s
    ##   6750K .......... .......... .......... .......... ..........  5% 72.4M 5s
    ##   6800K .......... .......... .......... .......... ..........  5% 18.4M 5s
    ##   6850K .......... .......... .......... .......... ..........  5% 58.8M 5s
    ##   6900K .......... .......... .......... .......... ..........  5%  118M 5s
    ##   6950K .......... .......... .......... .......... ..........  5% 19.3M 5s
    ##   7000K .......... .......... .......... .......... ..........  5% 86.7M 5s
    ##   7050K .......... .......... .......... .......... ..........  5%  143M 4s
    ##   7100K .......... .......... .......... .......... ..........  5% 17.6M 5s
    ##   7150K .......... .......... .......... .......... ..........  5% 95.8M 4s
    ##   7200K .......... .......... .......... .......... ..........  5% 15.9M 5s
    ##   7250K .......... .......... .......... .......... ..........  5%  100M 4s
    ##   7300K .......... .......... .......... .......... ..........  5%  109M 4s
    ##   7350K .......... .......... .......... .......... ..........  5% 16.0M 4s
    ##   7400K .......... .......... .......... .......... ..........  5%  116M 4s
    ##   7450K .......... .......... .......... .......... ..........  5% 86.0M 4s
    ##   7500K .......... .......... .......... .......... ..........  5% 18.1M 4s
    ##   7550K .......... .......... .......... .......... ..........  5%  101M 4s
    ##   7600K .......... .......... .......... .......... ..........  5% 99.7M 4s
    ##   7650K .......... .......... .......... .......... ..........  5% 13.1M 4s
    ##   7700K .......... .......... .......... .......... ..........  5%  103M 4s
    ##   7750K .......... .......... .......... .......... ..........  5%  151M 4s
    ##   7800K .......... .......... .......... .......... ..........  5% 12.9M 4s
    ##   7850K .......... .......... .......... .......... ..........  5%  128M 4s
    ##   7900K .......... .......... .......... .......... ..........  5% 13.3M 4s
    ##   7950K .......... .......... .......... .......... ..........  5%  123M 4s
    ##   8000K .......... .......... .......... .......... ..........  5%  120M 4s
    ##   8050K .......... .......... .......... .......... ..........  6% 16.9M 4s
    ##   8100K .......... .......... .......... .......... ..........  6% 72.9M 4s
    ##   8150K .......... .......... .......... .......... ..........  6%  123M 4s
    ##   8200K .......... .......... .......... .......... ..........  6% 16.5M 4s
    ##   8250K .......... .......... .......... .......... ..........  6%  103M 4s
    ##   8300K .......... .......... .......... .......... ..........  6%  120M 4s
    ##   8350K .......... .......... .......... .......... ..........  6% 13.8M 4s
    ##   8400K .......... .......... .......... .......... ..........  6%  107M 4s
    ##   8450K .......... .......... .......... .......... ..........  6% 60.8M 4s
    ##   8500K .......... .......... .......... .......... ..........  6% 28.4M 4s
    ##   8550K .......... .......... .......... .......... ..........  6% 86.0M 4s
    ##   8600K .......... .......... .......... .......... ..........  6% 26.2M 4s
    ##   8650K .......... .......... .......... .......... ..........  6% 32.1M 4s
    ##   8700K .......... .......... .......... .......... ..........  6% 55.5M 4s
    ##   8750K .......... .......... .......... .......... ..........  6% 50.4M 4s
    ##   8800K .......... .......... .......... .......... ..........  6% 30.1M 4s
    ##   8850K .......... .......... .......... .......... ..........  6% 44.7M 4s
    ##   8900K .......... .......... .......... .......... ..........  6% 34.6M 4s
    ##   8950K .......... .......... .......... .......... ..........  6% 77.8M 4s
    ##   9000K .......... .......... .......... .......... ..........  6% 28.3M 4s
    ##   9050K .......... .......... .......... .......... ..........  6% 43.8M 4s
    ##   9100K .......... .......... .......... .......... ..........  6% 27.5M 4s
    ##   9150K .......... .......... .......... .......... ..........  6% 32.3M 4s
    ##   9200K .......... .......... .......... .......... ..........  6% 85.8M 4s
    ##   9250K .......... .......... .......... .......... ..........  6% 32.6M 4s
    ##   9300K .......... .......... .......... .......... ..........  6% 29.6M 4s
    ##   9350K .......... .......... .......... .......... ..........  6% 33.8M 4s
    ##   9400K .......... .......... .......... .......... ..........  7%  104M 4s
    ##   9450K .......... .......... .......... .......... ..........  7% 28.1M 4s
    ##   9500K .......... .......... .......... .......... ..........  7% 26.6M 4s
    ##   9550K .......... .......... .......... .......... ..........  7% 30.8M 4s
    ##   9600K .......... .......... .......... .......... ..........  7% 27.8M 4s
    ##   9650K .......... .......... .......... .......... ..........  7% 29.4M 4s
    ##   9700K .......... .......... .......... .......... ..........  7%  112M 4s
    ##   9750K .......... .......... .......... .......... ..........  7% 20.5M 4s
    ##   9800K .......... .......... .......... .......... ..........  7% 46.9M 4s
    ##   9850K .......... .......... .......... .......... ..........  7% 20.0M 4s
    ##   9900K .......... .......... .......... .......... ..........  7% 88.7M 4s
    ##   9950K .......... .......... .......... .......... ..........  7% 88.7M 4s
    ##  10000K .......... .......... .......... .......... ..........  7% 17.5M 4s
    ##  10050K .......... .......... .......... .......... ..........  7% 73.8M 4s
    ##  10100K .......... .......... .......... .......... ..........  7% 36.3M 4s
    ##  10150K .......... .......... .......... .......... ..........  7% 39.2M 4s
    ##  10200K .......... .......... .......... .......... ..........  7% 56.0M 4s
    ##  10250K .......... .......... .......... .......... ..........  7% 28.0M 4s
    ##  10300K .......... .......... .......... .......... ..........  7% 33.6M 4s
    ##  10350K .......... .......... .......... .......... ..........  7%  107M 4s
    ##  10400K .......... .......... .......... .......... ..........  7% 26.5M 4s
    ##  10450K .......... .......... .......... .......... ..........  7% 48.4M 4s
    ##  10500K .......... .......... .......... .......... ..........  7% 29.4M 4s
    ##  10550K .......... .......... .......... .......... ..........  7% 31.9M 4s
    ##  10600K .......... .......... .......... .......... ..........  7% 90.1M 4s
    ##  10650K .......... .......... .......... .......... ..........  7% 36.1M 4s
    ##  10700K .......... .......... .......... .......... ..........  7% 23.0M 4s
    ##  10750K .......... .......... .......... .......... ..........  8% 31.4M 4s
    ##  10800K .......... .......... .......... .......... ..........  8% 29.3M 4s
    ##  10850K .......... .......... .......... .......... ..........  8% 96.9M 4s
    ##  10900K .......... .......... .......... .......... ..........  8% 46.4M 4s
    ##  10950K .......... .......... .......... .......... ..........  8% 22.3M 4s
    ##  11000K .......... .......... .......... .......... ..........  8% 49.7M 4s
    ##  11050K .......... .......... .......... .......... ..........  8% 19.5M 4s
    ##  11100K .......... .......... .......... .......... ..........  8% 83.5M 4s
    ##  11150K .......... .......... .......... .......... ..........  8% 69.7M 4s
    ##  11200K .......... .......... .......... .......... ..........  8% 22.0M 4s
    ##  11250K .......... .......... .......... .......... ..........  8% 45.4M 4s
    ##  11300K .......... .......... .......... .......... ..........  8% 21.4M 4s
    ##  11350K .......... .......... .......... .......... ..........  8% 94.3M 4s
    ##  11400K .......... .......... .......... .......... ..........  8% 46.1M 4s
    ##  11450K .......... .......... .......... .......... ..........  8% 20.2M 4s
    ##  11500K .......... .......... .......... .......... ..........  8% 62.4M 4s
    ##  11550K .......... .......... .......... .......... ..........  8% 92.7M 4s
    ##  11600K .......... .......... .......... .......... ..........  8% 18.6M 4s
    ##  11650K .......... .......... .......... .......... ..........  8% 85.8M 4s
    ##  11700K .......... .......... .......... .......... ..........  8% 15.1M 4s
    ##  11750K .......... .......... .......... .......... ..........  8% 88.1M 4s
    ##  11800K .......... .......... .......... .......... ..........  8% 14.1M 4s
    ##  11850K .......... .......... .......... .......... ..........  8% 95.6M 4s
    ##  11900K .......... .......... .......... .......... ..........  8% 92.6M 4s
    ##  11950K .......... .......... .......... .......... ..........  8% 13.4M 4s
    ##  12000K .......... .......... .......... .......... ..........  8% 73.1M 4s
    ##  12050K .......... .......... .......... .......... ..........  8%  117M 4s
    ##  12100K .......... .......... .......... .......... ..........  9% 12.7M 4s
    ##  12150K .......... .......... .......... .......... ..........  9% 98.9M 4s
    ##  12200K .......... .......... .......... .......... ..........  9% 86.2M 4s
    ##  12250K .......... .......... .......... .......... ..........  9% 28.9M 4s
    ##  12300K .......... .......... .......... .......... ..........  9% 51.8M 4s
    ##  12350K .......... .......... .......... .......... ..........  9%  101M 4s
    ##  12400K .......... .......... .......... .......... ..........  9% 23.2M 4s
    ##  12450K .......... .......... .......... .......... ..........  9% 39.8M 4s
    ##  12500K .......... .......... .......... .......... ..........  9% 76.3M 4s
    ##  12550K .......... .......... .......... .......... ..........  9% 14.8M 4s
    ##  12600K .......... .......... .......... .......... ..........  9% 90.5M 4s
    ##  12650K .......... .......... .......... .......... ..........  9% 14.1M 4s
    ##  12700K .......... .......... .......... .......... ..........  9%  103M 4s
    ##  12750K .......... .......... .......... .......... ..........  9% 93.6M 4s
    ##  12800K .......... .......... .......... .......... ..........  9% 14.2M 4s
    ##  12850K .......... .......... .......... .......... ..........  9%  106M 4s
    ##  12900K .......... .......... .......... .......... ..........  9% 14.6M 4s
    ##  12950K .......... .......... .......... .......... ..........  9% 76.8M 4s
    ##  13000K .......... .......... .......... .......... ..........  9% 89.4M 4s
    ##  13050K .......... .......... .......... .......... ..........  9% 16.4M 4s
    ##  13100K .......... .......... .......... .......... ..........  9% 82.2M 4s
    ##  13150K .......... .......... .......... .......... ..........  9% 86.2M 4s
    ##  13200K .......... .......... .......... .......... ..........  9% 26.2M 4s
    ##  13250K .......... .......... .......... .......... ..........  9% 77.5M 4s
    ##  13300K .......... .......... .......... .......... ..........  9% 92.8M 4s
    ##  13350K .......... .......... .......... .......... ..........  9% 19.7M 4s
    ##  13400K .......... .......... .......... .......... ..........  9% 86.7M 4s
    ##  13450K .......... .......... .......... .......... .......... 10% 85.4M 4s
    ##  13500K .......... .......... .......... .......... .......... 10% 20.4M 4s
    ##  13550K .......... .......... .......... .......... .......... 10% 61.9M 4s
    ##  13600K .......... .......... .......... .......... .......... 10% 20.2M 4s
    ##  13650K .......... .......... .......... .......... .......... 10% 86.1M 4s
    ##  13700K .......... .......... .......... .......... .......... 10% 64.2M 4s
    ##  13750K .......... .......... .......... .......... .......... 10% 17.1M 4s
    ##  13800K .......... .......... .......... .......... .......... 10% 78.2M 4s
    ##  13850K .......... .......... .......... .......... .......... 10% 20.7M 4s
    ##  13900K .......... .......... .......... .......... .......... 10%  101M 4s
    ##  13950K .......... .......... .......... .......... .......... 10% 47.8M 4s
    ##  14000K .......... .......... .......... .......... .......... 10% 21.5M 4s
    ##  14050K .......... .......... .......... .......... .......... 10% 56.6M 4s
    ##  14100K .......... .......... .......... .......... .......... 10% 30.0M 4s
    ##  14150K .......... .......... .......... .......... .......... 10% 27.0M 4s
    ##  14200K .......... .......... .......... .......... .......... 10% 82.9M 4s
    ##  14250K .......... .......... .......... .......... .......... 10% 49.8M 4s
    ##  14300K .......... .......... .......... .......... .......... 10% 21.5M 4s
    ##  14350K .......... .......... .......... .......... .......... 10% 21.3M 4s
    ##  14400K .......... .......... .......... .......... .......... 10% 84.4M 4s
    ##  14450K .......... .......... .......... .......... .......... 10% 82.1M 4s
    ##  14500K .......... .......... .......... .......... .......... 10% 60.7M 4s
    ##  14550K .......... .......... .......... .......... .......... 10% 17.6M 4s
    ##  14600K .......... .......... .......... .......... .......... 10% 76.4M 4s
    ##  14650K .......... .......... .......... .......... .......... 10%  101M 4s
    ##  14700K .......... .......... .......... .......... .......... 10% 14.8M 4s
    ##  14750K .......... .......... .......... .......... .......... 10%  104M 4s
    ##  14800K .......... .......... .......... .......... .......... 11% 15.3M 4s
    ##  14850K .......... .......... .......... .......... .......... 11% 66.2M 4s
    ##  14900K .......... .......... .......... .......... .......... 11% 85.6M 4s
    ##  14950K .......... .......... .......... .......... .......... 11%  114M 4s
    ##  15000K .......... .......... .......... .......... .......... 11% 23.0M 4s
    ##  15050K .......... .......... .......... .......... .......... 11% 72.8M 4s
    ##  15100K .......... .......... .......... .......... .......... 11% 27.8M 4s
    ##  15150K .......... .......... .......... .......... .......... 11% 80.4M 4s
    ##  15200K .......... .......... .......... .......... .......... 11% 36.4M 4s
    ##  15250K .......... .......... .......... .......... .......... 11% 26.1M 4s
    ##  15300K .......... .......... .......... .......... .......... 11% 32.2M 4s
    ##  15350K .......... .......... .......... .......... .......... 11% 85.2M 4s
    ##  15400K .......... .......... .......... .......... .......... 11% 23.3M 4s
    ##  15450K .......... .......... .......... .......... .......... 11% 42.2M 4s
    ##  15500K .......... .......... .......... .......... .......... 11% 20.6M 4s
    ##  15550K .......... .......... .......... .......... .......... 11% 98.4M 4s
    ##  15600K .......... .......... .......... .......... .......... 11% 38.0M 4s
    ##  15650K .......... .......... .......... .......... .......... 11% 28.1M 4s
    ##  15700K .......... .......... .......... .......... .......... 11% 31.9M 4s
    ##  15750K .......... .......... .......... .......... .......... 11% 38.1M 4s
    ##  15800K .......... .......... .......... .......... .......... 11% 77.3M 4s
    ##  15850K .......... .......... .......... .......... .......... 11% 24.5M 4s
    ##  15900K .......... .......... .......... .......... .......... 11% 32.7M 4s
    ##  15950K .......... .......... .......... .......... .......... 11% 24.9M 4s
    ##  16000K .......... .......... .......... .......... .......... 11% 41.5M 4s
    ##  16050K .......... .......... .......... .......... .......... 11% 67.8M 4s
    ##  16100K .......... .......... .......... .......... .......... 11% 28.1M 4s
    ##  16150K .......... .......... .......... .......... .......... 12% 31.5M 4s
    ##  16200K .......... .......... .......... .......... .......... 12% 29.0M 4s
    ##  16250K .......... .......... .......... .......... .......... 12% 25.3M 4s
    ##  16300K .......... .......... .......... .......... .......... 12% 74.4M 4s
    ##  16350K .......... .......... .......... .......... .......... 12% 43.3M 4s
    ##  16400K .......... .......... .......... .......... .......... 12% 24.4M 4s
    ##  16450K .......... .......... .......... .......... .......... 12% 78.7M 4s
    ##  16500K .......... .......... .......... .......... .......... 12% 74.4M 4s
    ##  16550K .......... .......... .......... .......... .......... 12% 24.6M 4s
    ##  16600K .......... .......... .......... .......... .......... 12% 37.4M 4s
    ##  16650K .......... .......... .......... .......... .......... 12% 23.3M 4s
    ##  16700K .......... .......... .......... .......... .......... 12% 87.7M 4s
    ##  16750K .......... .......... .......... .......... .......... 12% 65.9M 4s
    ##  16800K .......... .......... .......... .......... .......... 12% 19.0M 4s
    ##  16850K .......... .......... .......... .......... .......... 12% 73.1M 4s
    ##  16900K .......... .......... .......... .......... .......... 12% 12.9M 4s
    ##  16950K .......... .......... .......... .......... .......... 12% 53.1M 4s
    ##  17000K .......... .......... .......... .......... .......... 12% 72.0M 4s
    ##  17050K .......... .......... .......... .......... .......... 12% 94.3M 4s
    ##  17100K .......... .......... .......... .......... .......... 12% 15.2M 4s
    ##  17150K .......... .......... .......... .......... .......... 12% 67.9M 4s
    ##  17200K .......... .......... .......... .......... .......... 12% 83.9M 4s
    ##  17250K .......... .......... .......... .......... .......... 12% 95.3M 4s
    ##  17300K .......... .......... .......... .......... .......... 12% 26.0M 4s
    ##  17350K .......... .......... .......... .......... .......... 12% 36.6M 4s
    ##  17400K .......... .......... .......... .......... .......... 12% 81.9M 4s
    ##  17450K .......... .......... .......... .......... .......... 12% 85.3M 4s
    ##  17500K .......... .......... .......... .......... .......... 13% 21.4M 4s
    ##  17550K .......... .......... .......... .......... .......... 13% 46.3M 4s
    ##  17600K .......... .......... .......... .......... .......... 13% 86.2M 4s
    ##  17650K .......... .......... .......... .......... .......... 13% 79.6M 4s
    ##  17700K .......... .......... .......... .......... .......... 13% 34.8M 4s
    ##  17750K .......... .......... .......... .......... .......... 13% 33.7M 4s
    ##  17800K .......... .......... .......... .......... .......... 13% 83.9M 4s
    ##  17850K .......... .......... .......... .......... .......... 13% 39.1M 4s
    ##  17900K .......... .......... .......... .......... .......... 13%  101M 4s
    ##  17950K .......... .......... .......... .......... .......... 13% 26.1M 4s
    ##  18000K .......... .......... .......... .......... .......... 13% 95.6M 4s
    ##  18050K .......... .......... .......... .......... .......... 13% 85.0M 4s
    ##  18100K .......... .......... .......... .......... .......... 13% 18.2M 4s
    ##  18150K .......... .......... .......... .......... .......... 13% 90.5M 4s
    ##  18200K .......... .......... .......... .......... .......... 13% 81.3M 4s
    ##  18250K .......... .......... .......... .......... .......... 13% 96.2M 4s
    ##  18300K .......... .......... .......... .......... .......... 13% 17.8M 4s
    ##  18350K .......... .......... .......... .......... .......... 13% 63.5M 4s
    ##  18400K .......... .......... .......... .......... .......... 13% 86.3M 4s
    ##  18450K .......... .......... .......... .......... .......... 13% 95.5M 4s
    ##  18500K .......... .......... .......... .......... .......... 13% 19.6M 4s
    ##  18550K .......... .......... .......... .......... .......... 13% 80.9M 4s
    ##  18600K .......... .......... .......... .......... .......... 13% 76.0M 4s
    ##  18650K .......... .......... .......... .......... .......... 13% 90.4M 4s
    ##  18700K .......... .......... .......... .......... .......... 13% 29.9M 4s
    ##  18750K .......... .......... .......... .......... .......... 13% 54.8M 4s
    ##  18800K .......... .......... .......... .......... .......... 13% 75.9M 3s
    ##  18850K .......... .......... .......... .......... .......... 14% 95.2M 3s
    ##  18900K .......... .......... .......... .......... .......... 14% 28.3M 3s
    ##  18950K .......... .......... .......... .......... .......... 14% 62.0M 3s
    ##  19000K .......... .......... .......... .......... .......... 14% 72.2M 3s
    ##  19050K .......... .......... .......... .......... .......... 14% 86.6M 3s
    ##  19100K .......... .......... .......... .......... .......... 14% 39.1M 3s
    ##  19150K .......... .......... .......... .......... .......... 14% 35.3M 3s
    ##  19200K .......... .......... .......... .......... .......... 14% 87.2M 3s
    ##  19250K .......... .......... .......... .......... .......... 14% 33.2M 3s
    ##  19300K .......... .......... .......... .......... .......... 14% 39.0M 3s
    ##  19350K .......... .......... .......... .......... .......... 14% 93.9M 3s
    ##  19400K .......... .......... .......... .......... .......... 14% 77.2M 3s
    ##  19450K .......... .......... .......... .......... .......... 14% 42.3M 3s
    ##  19500K .......... .......... .......... .......... .......... 14% 24.6M 3s
    ##  19550K .......... .......... .......... .......... .......... 14% 97.9M 3s
    ##  19600K .......... .......... .......... .......... .......... 14% 91.6M 3s
    ##  19650K .......... .......... .......... .......... .......... 14% 57.8M 3s
    ##  19700K .......... .......... .......... .......... .......... 14% 25.7M 3s
    ##  19750K .......... .......... .......... .......... .......... 14% 57.7M 3s
    ##  19800K .......... .......... .......... .......... .......... 14% 26.6M 3s
    ##  19850K .......... .......... .......... .......... .......... 14% 84.0M 3s
    ##  19900K .......... .......... .......... .......... .......... 14% 25.6M 3s
    ##  19950K .......... .......... .......... .......... .......... 14%  107M 3s
    ##  20000K .......... .......... .......... .......... .......... 14% 69.4M 3s
    ##  20050K .......... .......... .......... .......... .......... 14% 89.4M 3s
    ##  20100K .......... .......... .......... .......... .......... 14% 20.3M 3s
    ##  20150K .......... .......... .......... .......... .......... 14%  103M 3s
    ##  20200K .......... .......... .......... .......... .......... 15% 96.2M 3s
    ##  20250K .......... .......... .......... .......... .......... 15%  126M 3s
    ##  20300K .......... .......... .......... .......... .......... 15% 24.4M 3s
    ##  20350K .......... .......... .......... .......... .......... 15% 70.6M 3s
    ##  20400K .......... .......... .......... .......... .......... 15% 81.3M 3s
    ##  20450K .......... .......... .......... .......... .......... 15% 24.2M 3s
    ##  20500K .......... .......... .......... .......... .......... 15% 60.4M 3s
    ##  20550K .......... .......... .......... .......... .......... 15% 90.5M 3s
    ##  20600K .......... .......... .......... .......... .......... 15% 84.1M 3s
    ##  20650K .......... .......... .......... .......... .......... 15% 30.1M 3s
    ##  20700K .......... .......... .......... .......... .......... 15% 43.8M 3s
    ##  20750K .......... .......... .......... .......... .......... 15% 95.1M 3s
    ##  20800K .......... .......... .......... .......... .......... 15% 51.7M 3s
    ##  20850K .......... .......... .......... .......... .......... 15% 38.9M 3s
    ##  20900K .......... .......... .......... .......... .......... 15% 43.3M 3s
    ##  20950K .......... .......... .......... .......... .......... 15% 65.8M 3s
    ##  21000K .......... .......... .......... .......... .......... 15% 83.9M 3s
    ##  21050K .......... .......... .......... .......... .......... 15% 47.3M 3s
    ##  21100K .......... .......... .......... .......... .......... 15% 20.6M 3s
    ##  21150K .......... .......... .......... .......... .......... 15% 83.8M 3s
    ##  21200K .......... .......... .......... .......... .......... 15% 93.9M 3s
    ##  21250K .......... .......... .......... .......... .......... 15%  113M 3s
    ##  21300K .......... .......... .......... .......... .......... 15% 23.1M 3s
    ##  21350K .......... .......... .......... .......... .......... 15% 81.3M 3s
    ##  21400K .......... .......... .......... .......... .......... 15% 80.8M 3s
    ##  21450K .......... .......... .......... .......... .......... 15% 33.9M 3s
    ##  21500K .......... .......... .......... .......... .......... 15% 83.1M 3s
    ##  21550K .......... .......... .......... .......... .......... 16% 18.9M 3s
    ##  21600K .......... .......... .......... .......... .......... 16% 93.1M 3s
    ##  21650K .......... .......... .......... .......... .......... 16%  109M 3s
    ##  21700K .......... .......... .......... .......... .......... 16% 96.1M 3s
    ##  21750K .......... .......... .......... .......... .......... 16% 20.4M 3s
    ##  21800K .......... .......... .......... .......... .......... 16% 99.6M 3s
    ##  21850K .......... .......... .......... .......... .......... 16% 89.2M 3s
    ##  21900K .......... .......... .......... .......... .......... 16% 80.7M 3s
    ##  21950K .......... .......... .......... .......... .......... 16% 38.2M 3s
    ##  22000K .......... .......... .......... .......... .......... 16% 70.1M 3s
    ##  22050K .......... .......... .......... .......... .......... 16% 37.7M 3s
    ##  22100K .......... .......... .......... .......... .......... 16% 31.8M 3s
    ##  22150K .......... .......... .......... .......... .......... 16% 53.6M 3s
    ##  22200K .......... .......... .......... .......... .......... 16% 93.8M 3s
    ##  22250K .......... .......... .......... .......... .......... 16% 72.8M 3s
    ##  22300K .......... .......... .......... .......... .......... 16% 21.0M 3s
    ##  22350K .......... .......... .......... .......... .......... 16% 75.7M 3s
    ##  22400K .......... .......... .......... .......... .......... 16% 73.1M 3s
    ##  22450K .......... .......... .......... .......... .......... 16%  102M 3s
    ##  22500K .......... .......... .......... .......... .......... 16% 26.0M 3s
    ##  22550K .......... .......... .......... .......... .......... 16% 43.1M 3s
    ##  22600K .......... .......... .......... .......... .......... 16% 79.2M 3s
    ##  22650K .......... .......... .......... .......... .......... 16%  110M 3s
    ##  22700K .......... .......... .......... .......... .......... 16% 24.0M 3s
    ##  22750K .......... .......... .......... .......... .......... 16% 29.9M 3s
    ##  22800K .......... .......... .......... .......... .......... 16% 96.9M 3s
    ##  22850K .......... .......... .......... .......... .......... 16%  105M 3s
    ##  22900K .......... .......... .......... .......... .......... 17% 42.2M 3s
    ##  22950K .......... .......... .......... .......... .......... 17% 34.1M 3s
    ##  23000K .......... .......... .......... .......... .......... 17% 60.2M 3s
    ##  23050K .......... .......... .......... .......... .......... 17% 74.7M 3s
    ##  23100K .......... .......... .......... .......... .......... 17% 25.4M 3s
    ##  23150K .......... .......... .......... .......... .......... 17% 77.3M 3s
    ##  23200K .......... .......... .......... .......... .......... 17% 71.8M 3s
    ##  23250K .......... .......... .......... .......... .......... 17% 86.3M 3s
    ##  23300K .......... .......... .......... .......... .......... 17% 27.6M 3s
    ##  23350K .......... .......... .......... .......... .......... 17% 78.8M 3s
    ##  23400K .......... .......... .......... .......... .......... 17% 59.6M 3s
    ##  23450K .......... .......... .......... .......... .......... 17% 94.7M 3s
    ##  23500K .......... .......... .......... .......... .......... 17% 60.7M 3s
    ##  23550K .......... .......... .......... .......... .......... 17% 21.9M 3s
    ##  23600K .......... .......... .......... .......... .......... 17% 78.6M 3s
    ##  23650K .......... .......... .......... .......... .......... 17%  115M 3s
    ##  23700K .......... .......... .......... .......... .......... 17% 45.9M 3s
    ##  23750K .......... .......... .......... .......... .......... 17% 32.0M 3s
    ##  23800K .......... .......... .......... .......... .......... 17% 81.6M 3s
    ##  23850K .......... .......... .......... .......... .......... 17%  102M 3s
    ##  23900K .......... .......... .......... .......... .......... 17% 39.9M 3s
    ##  23950K .......... .......... .......... .......... .......... 17% 33.9M 3s
    ##  24000K .......... .......... .......... .......... .......... 17% 74.9M 3s
    ##  24050K .......... .......... .......... .......... .......... 17% 99.2M 3s
    ##  24100K .......... .......... .......... .......... .......... 17% 43.9M 3s
    ##  24150K .......... .......... .......... .......... .......... 17% 17.3M 3s
    ##  24200K .......... .......... .......... .......... .......... 17% 97.0M 3s
    ##  24250K .......... .......... .......... .......... .......... 18% 89.3M 3s
    ##  24300K .......... .......... .......... .......... .......... 18% 94.5M 3s
    ##  24350K .......... .......... .......... .......... .......... 18% 20.8M 3s
    ##  24400K .......... .......... .......... .......... .......... 18% 60.4M 3s
    ##  24450K .......... .......... .......... .......... .......... 18% 97.6M 3s
    ##  24500K .......... .......... .......... .......... .......... 18% 89.8M 3s
    ##  24550K .......... .......... .......... .......... .......... 18% 17.1M 3s
    ##  24600K .......... .......... .......... .......... .......... 18% 78.1M 3s
    ##  24650K .......... .......... .......... .......... .......... 18% 93.6M 3s
    ##  24700K .......... .......... .......... .......... .......... 18% 17.8M 3s
    ##  24750K .......... .......... .......... .......... .......... 18% 91.9M 3s
    ##  24800K .......... .......... .......... .......... .......... 18% 71.4M 3s
    ##  24850K .......... .......... .......... .......... .......... 18% 85.3M 3s
    ##  24900K .......... .......... .......... .......... .......... 18% 26.3M 3s
    ##  24950K .......... .......... .......... .......... .......... 18% 49.3M 3s
    ##  25000K .......... .......... .......... .......... .......... 18% 66.7M 3s
    ##  25050K .......... .......... .......... .......... .......... 18%  106M 3s
    ##  25100K .......... .......... .......... .......... .......... 18% 32.1M 3s
    ##  25150K .......... .......... .......... .......... .......... 18% 67.7M 3s
    ##  25200K .......... .......... .......... .......... .......... 18% 76.2M 3s
    ##  25250K .......... .......... .......... .......... .......... 18% 97.3M 3s
    ##  25300K .......... .......... .......... .......... .......... 18% 31.6M 3s
    ##  25350K .......... .......... .......... .......... .......... 18% 34.8M 3s
    ##  25400K .......... .......... .......... .......... .......... 18% 64.1M 3s
    ##  25450K .......... .......... .......... .......... .......... 18% 72.9M 3s
    ##  25500K .......... .......... .......... .......... .......... 18% 52.5M 3s
    ##  25550K .......... .......... .......... .......... .......... 18% 45.4M 3s
    ##  25600K .......... .......... .......... .......... .......... 19% 40.3M 3s
    ##  25650K .......... .......... .......... .......... .......... 19% 71.1M 3s
    ##  25700K .......... .......... .......... .......... .......... 19% 49.6M 3s
    ##  25750K .......... .......... .......... .......... .......... 19% 87.3M 3s
    ##  25800K .......... .......... .......... .......... .......... 19% 39.6M 3s
    ##  25850K .......... .......... .......... .......... .......... 19% 51.7M 3s
    ##  25900K .......... .......... .......... .......... .......... 19% 73.0M 3s
    ##  25950K .......... .......... .......... .......... .......... 19%  100M 3s
    ##  26000K .......... .......... .......... .......... .......... 19% 33.3M 3s
    ##  26050K .......... .......... .......... .......... .......... 19% 48.8M 3s
    ##  26100K .......... .......... .......... .......... .......... 19% 47.0M 3s
    ##  26150K .......... .......... .......... .......... .......... 19% 73.4M 3s
    ##  26200K .......... .......... .......... .......... .......... 19% 38.9M 3s
    ##  26250K .......... .......... .......... .......... .......... 19% 57.1M 3s
    ##  26300K .......... .......... .......... .......... .......... 19% 61.6M 3s
    ##  26350K .......... .......... .......... .......... .......... 19% 36.6M 3s
    ##  26400K .......... .......... .......... .......... .......... 19% 37.1M 3s
    ##  26450K .......... .......... .......... .......... .......... 19% 98.0M 3s
    ##  26500K .......... .......... .......... .......... .......... 19% 64.1M 3s
    ##  26550K .......... .......... .......... .......... .......... 19% 55.6M 3s
    ##  26600K .......... .......... .......... .......... .......... 19% 30.3M 3s
    ##  26650K .......... .......... .......... .......... .......... 19%  106M 3s
    ##  26700K .......... .......... .......... .......... .......... 19% 49.9M 3s
    ##  26750K .......... .......... .......... .......... .......... 19%  101M 3s
    ##  26800K .......... .......... .......... .......... .......... 19% 32.5M 3s
    ##  26850K .......... .......... .......... .......... .......... 19%  106M 3s
    ##  26900K .......... .......... .......... .......... .......... 20% 26.3M 3s
    ##  26950K .......... .......... .......... .......... .......... 20% 89.0M 3s
    ##  27000K .......... .......... .......... .......... .......... 20% 24.1M 3s
    ##  27050K .......... .......... .......... .......... .......... 20% 41.2M 3s
    ##  27100K .......... .......... .......... .......... .......... 20% 43.6M 3s
    ##  27150K .......... .......... .......... .......... .......... 20%  127M 3s
    ##  27200K .......... .......... .......... .......... .......... 20% 35.6M 3s
    ##  27250K .......... .......... .......... .......... .......... 20% 80.6M 3s
    ##  27300K .......... .......... .......... .......... .......... 20% 60.0M 3s
    ##  27350K .......... .......... .......... .......... .......... 20% 90.3M 3s
    ##  27400K .......... .......... .......... .......... .......... 20% 16.8M 3s
    ##  27450K .......... .......... .......... .......... .......... 20% 84.9M 3s
    ##  27500K .......... .......... .......... .......... .......... 20% 99.4M 3s
    ##  27550K .......... .......... .......... .......... .......... 20% 77.8M 3s
    ##  27600K .......... .......... .......... .......... .......... 20% 43.1M 3s
    ##  27650K .......... .......... .......... .......... .......... 20% 50.4M 3s
    ##  27700K .......... .......... .......... .......... .......... 20% 57.6M 3s
    ##  27750K .......... .......... .......... .......... .......... 20% 47.9M 3s
    ##  27800K .......... .......... .......... .......... .......... 20% 42.1M 3s
    ##  27850K .......... .......... .......... .......... .......... 20% 84.9M 3s
    ##  27900K .......... .......... .......... .......... .......... 20% 39.7M 3s
    ##  27950K .......... .......... .......... .......... .......... 20% 56.4M 3s
    ##  28000K .......... .......... .......... .......... .......... 20% 42.4M 3s
    ##  28050K .......... .......... .......... .......... .......... 20%  104M 3s
    ##  28100K .......... .......... .......... .......... .......... 20% 50.2M 3s
    ##  28150K .......... .......... .......... .......... .......... 20% 40.0M 3s
    ##  28200K .......... .......... .......... .......... .......... 20% 53.5M 3s
    ##  28250K .......... .......... .......... .......... .......... 21% 52.5M 3s
    ##  28300K .......... .......... .......... .......... .......... 21% 36.4M 3s
    ##  28350K .......... .......... .......... .......... .......... 21% 88.2M 3s
    ##  28400K .......... .......... .......... .......... .......... 21% 60.7M 3s
    ##  28450K .......... .......... .......... .......... .......... 21% 48.5M 3s
    ##  28500K .......... .......... .......... .......... .......... 21% 36.9M 3s
    ##  28550K .......... .......... .......... .......... .......... 21%  115M 3s
    ##  28600K .......... .......... .......... .......... .......... 21% 42.0M 3s
    ##  28650K .......... .......... .......... .......... .......... 21% 94.5M 3s
    ##  28700K .......... .......... .......... .......... .......... 21% 28.0M 3s
    ##  28750K .......... .......... .......... .......... .......... 21% 38.5M 3s
    ##  28800K .......... .......... .......... .......... .......... 21% 96.1M 3s
    ##  28850K .......... .......... .......... .......... .......... 21% 71.3M 3s
    ##  28900K .......... .......... .......... .......... .......... 21% 72.8M 3s
    ##  28950K .......... .......... .......... .......... .......... 21% 17.9M 3s
    ##  29000K .......... .......... .......... .......... .......... 21% 70.5M 3s
    ##  29050K .......... .......... .......... .......... .......... 21%  122M 3s
    ##  29100K .......... .......... .......... .......... .......... 21%  101M 3s
    ##  29150K .......... .......... .......... .......... .......... 21% 18.3M 3s
    ##  29200K .......... .......... .......... .......... .......... 21%  105M 3s
    ##  29250K .......... .......... .......... .......... .......... 21%  106M 3s
    ##  29300K .......... .......... .......... .......... .......... 21%  111M 3s
    ##  29350K .......... .......... .......... .......... .......... 21% 18.2M 3s
    ##  29400K .......... .......... .......... .......... .......... 21% 96.2M 3s
    ##  29450K .......... .......... .......... .......... .......... 21% 64.3M 3s
    ##  29500K .......... .......... .......... .......... .......... 21%  117M 3s
    ##  29550K .......... .......... .......... .......... .......... 21% 24.1M 3s
    ##  29600K .......... .......... .......... .......... .......... 22% 80.7M 3s
    ##  29650K .......... .......... .......... .......... .......... 22% 81.0M 3s
    ##  29700K .......... .......... .......... .......... .......... 22% 16.3M 3s
    ##  29750K .......... .......... .......... .......... .......... 22% 91.6M 3s
    ##  29800K .......... .......... .......... .......... .......... 22%  102M 3s
    ##  29850K .......... .......... .......... .......... .......... 22%  115M 3s
    ##  29900K .......... .......... .......... .......... .......... 22% 17.7M 3s
    ##  29950K .......... .......... .......... .......... .......... 22% 92.3M 3s
    ##  30000K .......... .......... .......... .......... .......... 22%  108M 3s
    ##  30050K .......... .......... .......... .......... .......... 22% 91.8M 3s
    ##  30100K .......... .......... .......... .......... .......... 22% 22.6M 3s
    ##  30150K .......... .......... .......... .......... .......... 22%  136M 3s
    ##  30200K .......... .......... .......... .......... .......... 22% 33.3M 3s
    ##  30250K .......... .......... .......... .......... .......... 22% 89.8M 3s
    ##  30300K .......... .......... .......... .......... .......... 22% 38.9M 3s
    ##  30350K .......... .......... .......... .......... .......... 22% 23.8M 3s
    ##  30400K .......... .......... .......... .......... .......... 22% 97.2M 3s
    ##  30450K .......... .......... .......... .......... .......... 22%  133M 3s
    ##  30500K .......... .......... .......... .......... .......... 22% 46.1M 3s
    ##  30550K .......... .......... .......... .......... .......... 22% 23.2M 3s
    ##  30600K .......... .......... .......... .......... .......... 22% 74.1M 3s
    ##  30650K .......... .......... .......... .......... .......... 22% 73.6M 3s
    ##  30700K .......... .......... .......... .......... .......... 22% 82.7M 3s
    ##  30750K .......... .......... .......... .......... .......... 22% 21.7M 3s
    ##  30800K .......... .......... .......... .......... .......... 22% 71.5M 3s
    ##  30850K .......... .......... .......... .......... .......... 22% 82.6M 3s
    ##  30900K .......... .......... .......... .......... .......... 22% 23.4M 3s
    ##  30950K .......... .......... .......... .......... .......... 23% 74.4M 3s
    ##  31000K .......... .......... .......... .......... .......... 23% 66.5M 3s
    ##  31050K .......... .......... .......... .......... .......... 23% 77.1M 3s
    ##  31100K .......... .......... .......... .......... .......... 23% 29.2M 3s
    ##  31150K .......... .......... .......... .......... .......... 23% 68.5M 3s
    ##  31200K .......... .......... .......... .......... .......... 23% 63.2M 3s
    ##  31250K .......... .......... .......... .......... .......... 23% 95.8M 3s
    ##  31300K .......... .......... .......... .......... .......... 23% 27.1M 3s
    ##  31350K .......... .......... .......... .......... .......... 23% 89.5M 3s
    ##  31400K .......... .......... .......... .......... .......... 23% 69.6M 3s
    ##  31450K .......... .......... .......... .......... .......... 23% 95.5M 3s
    ##  31500K .......... .......... .......... .......... .......... 23% 21.7M 3s
    ##  31550K .......... .......... .......... .......... .......... 23% 65.1M 3s
    ##  31600K .......... .......... .......... .......... .......... 23% 79.1M 3s
    ##  31650K .......... .......... .......... .......... .......... 23% 92.2M 3s
    ##  31700K .......... .......... .......... .......... .......... 23% 23.6M 3s
    ##  31750K .......... .......... .......... .......... .......... 23% 86.4M 3s
    ##  31800K .......... .......... .......... .......... .......... 23% 77.8M 3s
    ##  31850K .......... .......... .......... .......... .......... 23% 85.7M 3s
    ##  31900K .......... .......... .......... .......... .......... 23% 25.9M 3s
    ##  31950K .......... .......... .......... .......... .......... 23% 87.6M 3s
    ##  32000K .......... .......... .......... .......... .......... 23% 72.3M 3s
    ##  32050K .......... .......... .......... .......... .......... 23% 88.3M 3s
    ##  32100K .......... .......... .......... .......... .......... 23% 21.5M 3s
    ##  32150K .......... .......... .......... .......... .......... 23% 94.4M 3s
    ##  32200K .......... .......... .......... .......... .......... 23% 49.1M 3s
    ##  32250K .......... .......... .......... .......... .......... 23% 34.8M 3s
    ##  32300K .......... .......... .......... .......... .......... 24% 57.5M 3s
    ##  32350K .......... .......... .......... .......... .......... 24% 84.8M 3s
    ##  32400K .......... .......... .......... .......... .......... 24% 70.1M 3s
    ##  32450K .......... .......... .......... .......... .......... 24% 32.2M 3s
    ##  32500K .......... .......... .......... .......... .......... 24% 69.3M 3s
    ##  32550K .......... .......... .......... .......... .......... 24% 95.4M 3s
    ##  32600K .......... .......... .......... .......... .......... 24% 74.0M 3s
    ##  32650K .......... .......... .......... .......... .......... 24% 21.1M 3s
    ##  32700K .......... .......... .......... .......... .......... 24% 57.2M 3s
    ##  32750K .......... .......... .......... .......... .......... 24% 94.6M 3s
    ##  32800K .......... .......... .......... .......... .......... 24% 87.6M 3s
    ##  32850K .......... .......... .......... .......... .......... 24% 29.3M 3s
    ##  32900K .......... .......... .......... .......... .......... 24% 64.8M 3s
    ##  32950K .......... .......... .......... .......... .......... 24% 80.2M 3s
    ##  33000K .......... .......... .......... .......... .......... 24% 30.5M 3s
    ##  33050K .......... .......... .......... .......... .......... 24% 65.1M 3s
    ##  33100K .......... .......... .......... .......... .......... 24% 66.6M 3s
    ##  33150K .......... .......... .......... .......... .......... 24%  100M 3s
    ##  33200K .......... .......... .......... .......... .......... 24% 25.3M 3s
    ##  33250K .......... .......... .......... .......... .......... 24% 71.2M 3s
    ##  33300K .......... .......... .......... .......... .......... 24% 70.5M 3s
    ##  33350K .......... .......... .......... .......... .......... 24% 79.8M 3s
    ##  33400K .......... .......... .......... .......... .......... 24% 22.2M 3s
    ##  33450K .......... .......... .......... .......... .......... 24% 74.7M 3s
    ##  33500K .......... .......... .......... .......... .......... 24% 62.6M 3s
    ##  33550K .......... .......... .......... .......... .......... 24% 97.0M 3s
    ##  33600K .......... .......... .......... .......... .......... 24% 32.5M 3s
    ##  33650K .......... .......... .......... .......... .......... 25% 66.4M 3s
    ##  33700K .......... .......... .......... .......... .......... 25% 87.9M 3s
    ##  33750K .......... .......... .......... .......... .......... 25% 82.8M 3s
    ##  33800K .......... .......... .......... .......... .......... 25% 30.1M 3s
    ##  33850K .......... .......... .......... .......... .......... 25% 76.5M 3s
    ##  33900K .......... .......... .......... .......... .......... 25% 72.2M 3s
    ##  33950K .......... .......... .......... .......... .......... 25% 32.2M 3s
    ##  34000K .......... .......... .......... .......... .......... 25% 68.8M 3s
    ##  34050K .......... .......... .......... .......... .......... 25% 9.63M 3s
    ##  34100K .......... .......... .......... .......... .......... 25% 41.0M 3s
    ##  34150K .......... .......... .......... .......... .......... 25% 6.53M 3s
    ##  34200K .......... .......... .......... .......... .......... 25% 67.5M 3s
    ##  34250K .......... .......... .......... .......... .......... 25% 18.8M 3s
    ##  34300K .......... .......... .......... .......... .......... 25% 3.12M 3s
    ##  34350K .......... .......... .......... .......... .......... 25% 65.8M 3s
    ##  34400K .......... .......... .......... .......... .......... 25% 19.8M 3s
    ##  34450K .......... .......... .......... .......... .......... 25% 89.2M 3s
    ##  34500K .......... .......... .......... .......... .......... 25% 93.9M 3s
    ##  34550K .......... .......... .......... .......... .......... 25% 19.5M 3s
    ##  34600K .......... .......... .......... .......... .......... 25% 82.1M 3s
    ##  34650K .......... .......... .......... .......... .......... 25% 82.6M 3s
    ##  34700K .......... .......... .......... .......... .......... 25% 88.3M 3s
    ##  34750K .......... .......... .......... .......... .......... 25% 79.7M 3s
    ##  34800K .......... .......... .......... .......... .......... 25% 40.0M 3s
    ##  34850K .......... .......... .......... .......... .......... 25% 59.0M 3s
    ##  34900K .......... .......... .......... .......... .......... 25% 80.8M 3s
    ##  34950K .......... .......... .......... .......... .......... 25% 80.6M 3s
    ##  35000K .......... .......... .......... .......... .......... 26% 36.8M 3s
    ##  35050K .......... .......... .......... .......... .......... 26% 55.7M 3s
    ##  35100K .......... .......... .......... .......... .......... 26% 50.6M 3s
    ##  35150K .......... .......... .......... .......... .......... 26% 77.8M 3s
    ##  35200K .......... .......... .......... .......... .......... 26% 59.0M 3s
    ##  35250K .......... .......... .......... .......... .......... 26% 69.5M 3s
    ##  35300K .......... .......... .......... .......... .......... 26% 40.6M 3s
    ##  35350K .......... .......... .......... .......... .......... 26% 32.2M 3s
    ##  35400K .......... .......... .......... .......... .......... 26% 71.8M 3s
    ##  35450K .......... .......... .......... .......... .......... 26% 95.8M 3s
    ##  35500K .......... .......... .......... .......... .......... 26% 67.5M 3s
    ##  35550K .......... .......... .......... .......... .......... 26% 34.0M 3s
    ##  35600K .......... .......... .......... .......... .......... 26% 71.0M 3s
    ##  35650K .......... .......... .......... .......... .......... 26% 64.9M 3s
    ##  35700K .......... .......... .......... .......... .......... 26% 48.9M 3s
    ##  35750K .......... .......... .......... .......... .......... 26% 32.1M 3s
    ##  35800K .......... .......... .......... .......... .......... 26% 81.4M 3s
    ##  35850K .......... .......... .......... .......... .......... 26% 90.5M 3s
    ##  35900K .......... .......... .......... .......... .......... 26% 45.7M 3s
    ##  35950K .......... .......... .......... .......... .......... 26% 35.1M 3s
    ##  36000K .......... .......... .......... .......... .......... 26% 71.5M 3s
    ##  36050K .......... .......... .......... .......... .......... 26%  101M 3s
    ##  36100K .......... .......... .......... .......... .......... 26% 46.9M 3s
    ##  36150K .......... .......... .......... .......... .......... 26% 35.9M 3s
    ##  36200K .......... .......... .......... .......... .......... 26% 81.9M 3s
    ##  36250K .......... .......... .......... .......... .......... 26% 40.5M 3s
    ##  36300K .......... .......... .......... .......... .......... 26% 94.2M 3s
    ##  36350K .......... .......... .......... .......... .......... 27% 31.4M 3s
    ##  36400K .......... .......... .......... .......... .......... 27% 38.6M 3s
    ##  36450K .......... .......... .......... .......... .......... 27% 54.7M 3s
    ##  36500K .......... .......... .......... .......... .......... 27% 24.1M 3s
    ##  36550K .......... .......... .......... .......... .......... 27% 78.6M 3s
    ##  36600K .......... .......... .......... .......... .......... 27% 85.4M 3s
    ##  36650K .......... .......... .......... .......... .......... 27% 35.6M 3s
    ##  36700K .......... .......... .......... .......... .......... 27% 64.0M 3s
    ##  36750K .......... .......... .......... .......... .......... 27% 30.1M 3s
    ##  36800K .......... .......... .......... .......... .......... 27% 67.0M 3s
    ##  36850K .......... .......... .......... .......... .......... 27% 58.7M 3s
    ##  36900K .......... .......... .......... .......... .......... 27% 69.4M 3s
    ##  36950K .......... .......... .......... .......... .......... 27% 96.7M 3s
    ##  37000K .......... .......... .......... .......... .......... 27%  106M 3s
    ##  37050K .......... .......... .......... .......... .......... 27% 18.9M 3s
    ##  37100K .......... .......... .......... .......... .......... 27% 83.7M 3s
    ##  37150K .......... .......... .......... .......... .......... 27% 36.4M 3s
    ##  37200K .......... .......... .......... .......... .......... 27%  104M 3s
    ##  37250K .......... .......... .......... .......... .......... 27%  110M 3s
    ##  37300K .......... .......... .......... .......... .......... 27% 77.6M 3s
    ##  37350K .......... .......... .......... .......... .......... 27% 26.8M 3s
    ##  37400K .......... .......... .......... .......... .......... 27%  101M 3s
    ##  37450K .......... .......... .......... .......... .......... 27% 95.4M 3s
    ##  37500K .......... .......... .......... .......... .......... 27% 96.1M 3s
    ##  37550K .......... .......... .......... .......... .......... 27%  112M 2s
    ##  37600K .......... .......... .......... .......... .......... 27% 14.3M 3s
    ##  37650K .......... .......... .......... .......... .......... 27% 68.9M 2s
    ##  37700K .......... .......... .......... .......... .......... 28% 94.4M 2s
    ##  37750K .......... .......... .......... .......... .......... 28% 99.5M 2s
    ##  37800K .......... .......... .......... .......... .......... 28%  103M 2s
    ##  37850K .......... .......... .......... .......... .......... 28% 28.1M 2s
    ##  37900K .......... .......... .......... .......... .......... 28% 77.1M 2s
    ##  37950K .......... .......... .......... .......... .......... 28% 84.2M 2s
    ##  38000K .......... .......... .......... .......... .......... 28% 82.2M 2s
    ##  38050K .......... .......... .......... .......... .......... 28%  106M 2s
    ##  38100K .......... .......... .......... .......... .......... 28% 47.3M 2s
    ##  38150K .......... .......... .......... .......... .......... 28% 31.0M 2s
    ##  38200K .......... .......... .......... .......... .......... 28% 69.6M 2s
    ##  38250K .......... .......... .......... .......... .......... 28% 75.5M 2s
    ##  38300K .......... .......... .......... .......... .......... 28%  112M 2s
    ##  38350K .......... .......... .......... .......... .......... 28% 90.0M 2s
    ##  38400K .......... .......... .......... .......... .......... 28% 35.7M 2s
    ##  38450K .......... .......... .......... .......... .......... 28% 82.4M 2s
    ##  38500K .......... .......... .......... .......... .......... 28% 29.7M 2s
    ##  38550K .......... .......... .......... .......... .......... 28% 89.5M 2s
    ##  38600K .......... .......... .......... .......... .......... 28%  111M 2s
    ##  38650K .......... .......... .......... .......... .......... 28%  115M 2s
    ##  38700K .......... .......... .......... .......... .......... 28% 67.8M 2s
    ##  38750K .......... .......... .......... .......... .......... 28% 29.2M 2s
    ##  38800K .......... .......... .......... .......... .......... 28% 81.2M 2s
    ##  38850K .......... .......... .......... .......... .......... 28%  104M 2s
    ##  38900K .......... .......... .......... .......... .......... 28%  116M 2s
    ##  38950K .......... .......... .......... .......... .......... 28%  103M 2s
    ##  39000K .......... .......... .......... .......... .......... 28% 21.1M 2s
    ##  39050K .......... .......... .......... .......... .......... 29% 65.9M 2s
    ##  39100K .......... .......... .......... .......... .......... 29% 74.9M 2s
    ##  39150K .......... .......... .......... .......... .......... 29%  109M 2s
    ##  39200K .......... .......... .......... .......... .......... 29% 89.9M 2s
    ##  39250K .......... .......... .......... .......... .......... 29% 38.1M 2s
    ##  39300K .......... .......... .......... .......... .......... 29% 48.0M 2s
    ##  39350K .......... .......... .......... .......... .......... 29%  114M 2s
    ##  39400K .......... .......... .......... .......... .......... 29% 90.6M 2s
    ##  39450K .......... .......... .......... .......... .......... 29% 87.0M 2s
    ##  39500K .......... .......... .......... .......... .......... 29% 20.1M 2s
    ##  39550K .......... .......... .......... .......... .......... 29%  106M 2s
    ##  39600K .......... .......... .......... .......... .......... 29%  102M 2s
    ##  39650K .......... .......... .......... .......... .......... 29%  122M 2s
    ##  39700K .......... .......... .......... .......... .......... 29%  102M 2s
    ##  39750K .......... .......... .......... .......... .......... 29% 25.0M 2s
    ##  39800K .......... .......... .......... .......... .......... 29% 88.1M 2s
    ##  39850K .......... .......... .......... .......... .......... 29% 90.4M 2s
    ##  39900K .......... .......... .......... .......... .......... 29% 81.3M 2s
    ##  39950K .......... .......... .......... .......... .......... 29%  113M 2s
    ##  40000K .......... .......... .......... .......... .......... 29% 43.5M 2s
    ##  40050K .......... .......... .......... .......... .......... 29% 35.3M 2s
    ##  40100K .......... .......... .......... .......... .......... 29% 92.7M 2s
    ##  40150K .......... .......... .......... .......... .......... 29% 43.4M 2s
    ##  40200K .......... .......... .......... .......... .......... 29% 98.1M 2s
    ##  40250K .......... .......... .......... .......... .......... 29% 99.5M 2s
    ##  40300K .......... .......... .......... .......... .......... 29%  105M 2s
    ##  40350K .......... .......... .......... .......... .......... 29% 42.6M 2s
    ##  40400K .......... .......... .......... .......... .......... 30% 27.1M 2s
    ##  40450K .......... .......... .......... .......... .......... 30% 98.0M 2s
    ##  40500K .......... .......... .......... .......... .......... 30%  104M 2s
    ##  40550K .......... .......... .......... .......... .......... 30%  104M 2s
    ##  40600K .......... .......... .......... .......... .......... 30%  120M 2s
    ##  40650K .......... .......... .......... .......... .......... 30% 26.7M 2s
    ##  40700K .......... .......... .......... .......... .......... 30% 87.9M 2s
    ##  40750K .......... .......... .......... .......... .......... 30%  112M 2s
    ##  40800K .......... .......... .......... .......... .......... 30% 87.2M 2s
    ##  40850K .......... .......... .......... .......... .......... 30% 84.4M 2s
    ##  40900K .......... .......... .......... .......... .......... 30% 38.5M 2s
    ##  40950K .......... .......... .......... .......... .......... 30% 75.2M 2s
    ##  41000K .......... .......... .......... .......... .......... 30% 72.4M 2s
    ##  41050K .......... .......... .......... .......... .......... 30% 96.7M 2s
    ##  41100K .......... .......... .......... .......... .......... 30% 66.1M 2s
    ##  41150K .......... .......... .......... .......... .......... 30% 16.6M 2s
    ##  41200K .......... .......... .......... .......... .......... 30% 92.0M 2s
    ##  41250K .......... .......... .......... .......... .......... 30% 60.7M 2s
    ##  41300K .......... .......... .......... .......... .......... 30% 80.3M 2s
    ##  41350K .......... .......... .......... .......... .......... 30%  100M 2s
    ##  41400K .......... .......... .......... .......... .......... 30% 32.1M 2s
    ##  41450K .......... .......... .......... .......... .......... 30% 62.0M 2s
    ##  41500K .......... .......... .......... .......... .......... 30% 56.9M 2s
    ##  41550K .......... .......... .......... .......... .......... 30% 69.6M 2s
    ##  41600K .......... .......... .......... .......... .......... 30%  101M 2s
    ##  41650K .......... .......... .......... .......... .......... 30% 90.9M 2s
    ##  41700K .......... .......... .......... .......... .......... 30% 79.0M 2s
    ##  41750K .......... .......... .......... .......... .......... 31%  110M 2s
    ##  41800K .......... .......... .......... .......... .......... 31% 37.9M 2s
    ##  41850K .......... .......... .......... .......... .......... 31% 89.7M 2s
    ##  41900K .......... .......... .......... .......... .......... 31% 41.8M 2s
    ##  41950K .......... .......... .......... .......... .......... 31% 65.6M 2s
    ##  42000K .......... .......... .......... .......... .......... 31% 94.5M 2s
    ##  42050K .......... .......... .......... .......... .......... 31% 42.5M 2s
    ##  42100K .......... .......... .......... .......... .......... 31% 79.9M 2s
    ##  42150K .......... .......... .......... .......... .......... 31% 16.0M 2s
    ##  42200K .......... .......... .......... .......... .......... 31% 63.1M 2s
    ##  42250K .......... .......... .......... .......... .......... 31% 92.5M 2s
    ##  42300K .......... .......... .......... .......... .......... 31% 95.1M 2s
    ##  42350K .......... .......... .......... .......... .......... 31%  105M 2s
    ##  42400K .......... .......... .......... .......... .......... 31% 21.0M 2s
    ##  42450K .......... .......... .......... .......... .......... 31% 65.9M 2s
    ##  42500K .......... .......... .......... .......... .......... 31% 80.5M 2s
    ##  42550K .......... .......... .......... .......... .......... 31% 67.2M 2s
    ##  42600K .......... .......... .......... .......... .......... 31%  103M 2s
    ##  42650K .......... .......... .......... .......... .......... 31%  120M 2s
    ##  42700K .......... .......... .......... .......... .......... 31% 61.9M 2s
    ##  42750K .......... .......... .......... .......... .......... 31% 30.0M 2s
    ##  42800K .......... .......... .......... .......... .......... 31% 87.4M 2s
    ##  42850K .......... .......... .......... .......... .......... 31% 55.7M 2s
    ##  42900K .......... .......... .......... .......... .......... 31% 87.7M 2s
    ##  42950K .......... .......... .......... .......... .......... 31%  106M 2s
    ##  43000K .......... .......... .......... .......... .......... 31% 59.1M 2s
    ##  43050K .......... .......... .......... .......... .......... 31% 27.0M 2s
    ##  43100K .......... .......... .......... .......... .......... 32% 90.7M 2s
    ##  43150K .......... .......... .......... .......... .......... 32% 88.6M 2s
    ##  43200K .......... .......... .......... .......... .......... 32% 82.6M 2s
    ##  43250K .......... .......... .......... .......... .......... 32%  103M 2s
    ##  43300K .......... .......... .......... .......... .......... 32% 24.8M 2s
    ##  43350K .......... .......... .......... .......... .......... 32% 92.1M 2s
    ##  43400K .......... .......... .......... .......... .......... 32% 41.6M 2s
    ##  43450K .......... .......... .......... .......... .......... 32% 88.5M 2s
    ##  43500K .......... .......... .......... .......... .......... 32% 86.1M 2s
    ##  43550K .......... .......... .......... .......... .......... 32% 81.4M 2s
    ##  43600K .......... .......... .......... .......... .......... 32%  103M 2s
    ##  43650K .......... .......... .......... .......... .......... 32% 30.0M 2s
    ##  43700K .......... .......... .......... .......... .......... 32% 67.7M 2s
    ##  43750K .......... .......... .......... .......... .......... 32% 90.4M 2s
    ##  43800K .......... .......... .......... .......... .......... 32% 91.2M 2s
    ##  43850K .......... .......... .......... .......... .......... 32% 35.5M 2s
    ##  43900K .......... .......... .......... .......... .......... 32% 67.0M 2s
    ##  43950K .......... .......... .......... .......... .......... 32% 98.0M 2s
    ##  44000K .......... .......... .......... .......... .......... 32% 70.2M 2s
    ##  44050K .......... .......... .......... .......... .......... 32% 93.1M 2s
    ##  44100K .......... .......... .......... .......... .......... 32% 85.2M 2s
    ##  44150K .......... .......... .......... .......... .......... 32% 68.9M 2s
    ##  44200K .......... .......... .......... .......... .......... 32% 76.7M 2s
    ##  44250K .......... .......... .......... .......... .......... 32% 67.3M 2s
    ##  44300K .......... .......... .......... .......... .......... 32% 60.2M 2s
    ##  44350K .......... .......... .......... .......... .......... 32% 92.7M 2s
    ##  44400K .......... .......... .......... .......... .......... 32% 50.5M 2s
    ##  44450K .......... .......... .......... .......... .......... 33% 63.8M 2s
    ##  44500K .......... .......... .......... .......... .......... 33% 86.8M 2s
    ##  44550K .......... .......... .......... .......... .......... 33% 86.5M 2s
    ##  44600K .......... .......... .......... .......... .......... 33% 67.5M 2s
    ##  44650K .......... .......... .......... .......... .......... 33% 51.6M 2s
    ##  44700K .......... .......... .......... .......... .......... 33% 48.4M 2s
    ##  44750K .......... .......... .......... .......... .......... 33% 33.9M 2s
    ##  44800K .......... .......... .......... .......... .......... 33% 75.8M 2s
    ##  44850K .......... .......... .......... .......... .......... 33%  104M 2s
    ##  44900K .......... .......... .......... .......... .......... 33% 86.0M 2s
    ##  44950K .......... .......... .......... .......... .......... 33% 83.7M 2s
    ##  45000K .......... .......... .......... .......... .......... 33% 32.2M 2s
    ##  45050K .......... .......... .......... .......... .......... 33% 66.1M 2s
    ##  45100K .......... .......... .......... .......... .......... 33%  101M 2s
    ##  45150K .......... .......... .......... .......... .......... 33% 51.5M 2s
    ##  45200K .......... .......... .......... .......... .......... 33% 70.8M 2s
    ##  45250K .......... .......... .......... .......... .......... 33% 86.1M 2s
    ##  45300K .......... .......... .......... .......... .......... 33% 61.1M 2s
    ##  45350K .......... .......... .......... .......... .......... 33% 42.8M 2s
    ##  45400K .......... .......... .......... .......... .......... 33% 78.9M 2s
    ##  45450K .......... .......... .......... .......... .......... 33% 64.5M 2s
    ##  45500K .......... .......... .......... .......... .......... 33% 76.4M 2s
    ##  45550K .......... .......... .......... .......... .......... 33% 88.8M 2s
    ##  45600K .......... .......... .......... .......... .......... 33% 71.0M 2s
    ##  45650K .......... .......... .......... .......... .......... 33% 51.6M 2s
    ##  45700K .......... .......... .......... .......... .......... 33% 78.7M 2s
    ##  45750K .......... .......... .......... .......... .......... 33% 70.2M 2s
    ##  45800K .......... .......... .......... .......... .......... 34% 48.3M 2s
    ##  45850K .......... .......... .......... .......... .......... 34% 83.1M 2s
    ##  45900K .......... .......... .......... .......... .......... 34% 70.6M 2s
    ##  45950K .......... .......... .......... .......... .......... 34% 76.6M 2s
    ##  46000K .......... .......... .......... .......... .......... 34% 18.3M 2s
    ##  46050K .......... .......... .......... .......... .......... 34% 51.6M 2s
    ##  46100K .......... .......... .......... .......... .......... 34% 97.5M 2s
    ##  46150K .......... .......... .......... .......... .......... 34% 86.9M 2s
    ##  46200K .......... .......... .......... .......... .......... 34%  107M 2s
    ##  46250K .......... .......... .......... .......... .......... 34%  106M 2s
    ##  46300K .......... .......... .......... .......... .......... 34% 83.7M 2s
    ##  46350K .......... .......... .......... .......... .......... 34%  112M 2s
    ##  46400K .......... .......... .......... .......... .......... 34% 67.9M 2s
    ##  46450K .......... .......... .......... .......... .......... 34% 96.6M 2s
    ##  46500K .......... .......... .......... .......... .......... 34% 86.7M 2s
    ##  46550K .......... .......... .......... .......... .......... 34% 40.4M 2s
    ##  46600K .......... .......... .......... .......... .......... 34% 47.6M 2s
    ##  46650K .......... .......... .......... .......... .......... 34% 65.7M 2s
    ##  46700K .......... .......... .......... .......... .......... 34% 61.5M 2s
    ##  46750K .......... .......... .......... .......... .......... 34% 63.2M 2s
    ##  46800K .......... .......... .......... .......... .......... 34%  114M 2s
    ##  46850K .......... .......... .......... .......... .......... 34% 64.0M 2s
    ##  46900K .......... .......... .......... .......... .......... 34% 51.7M 2s
    ##  46950K .......... .......... .......... .......... .......... 34% 70.9M 2s
    ##  47000K .......... .......... .......... .......... .......... 34% 59.8M 2s
    ##  47050K .......... .......... .......... .......... .......... 34%  134M 2s
    ##  47100K .......... .......... .......... .......... .......... 34% 63.8M 2s
    ##  47150K .......... .......... .......... .......... .......... 35%  109M 2s
    ##  47200K .......... .......... .......... .......... .......... 35% 35.6M 2s
    ##  47250K .......... .......... .......... .......... .......... 35% 54.6M 2s
    ##  47300K .......... .......... .......... .......... .......... 35% 54.9M 2s
    ##  47350K .......... .......... .......... .......... .......... 35%  136M 2s
    ##  47400K .......... .......... .......... .......... .......... 35%  109M 2s
    ##  47450K .......... .......... .......... .......... .......... 35% 40.6M 2s
    ##  47500K .......... .......... .......... .......... .......... 35%  105M 2s
    ##  47550K .......... .......... .......... .......... .......... 35% 73.6M 2s
    ##  47600K .......... .......... .......... .......... .......... 35% 64.3M 2s
    ##  47650K .......... .......... .......... .......... .......... 35%  138M 2s
    ##  47700K .......... .......... .......... .......... .......... 35% 27.3M 2s
    ##  47750K .......... .......... .......... .......... .......... 35% 76.1M 2s
    ##  47800K .......... .......... .......... .......... .......... 35% 67.2M 2s
    ##  47850K .......... .......... .......... .......... .......... 35%  132M 2s
    ##  47900K .......... .......... .......... .......... .......... 35% 89.4M 2s
    ##  47950K .......... .......... .......... .......... .......... 35%  130M 2s
    ##  48000K .......... .......... .......... .......... .......... 35% 45.9M 2s
    ##  48050K .......... .......... .......... .......... .......... 35% 73.7M 2s
    ##  48100K .......... .......... .......... .......... .......... 35% 73.2M 2s
    ##  48150K .......... .......... .......... .......... .......... 35%  113M 2s
    ##  48200K .......... .......... .......... .......... .......... 35% 45.5M 2s
    ##  48250K .......... .......... .......... .......... .......... 35% 72.8M 2s
    ##  48300K .......... .......... .......... .......... .......... 35% 63.7M 2s
    ##  48350K .......... .......... .......... .......... .......... 35%  115M 2s
    ##  48400K .......... .......... .......... .......... .......... 35% 70.0M 2s
    ##  48450K .......... .......... .......... .......... .......... 35% 55.2M 2s
    ##  48500K .......... .......... .......... .......... .......... 36% 65.2M 2s
    ##  48550K .......... .......... .......... .......... .......... 36% 63.8M 2s
    ##  48600K .......... .......... .......... .......... .......... 36%  110M 2s
    ##  48650K .......... .......... .......... .......... .......... 36% 38.5M 2s
    ##  48700K .......... .......... .......... .......... .......... 36%  116M 2s
    ##  48750K .......... .......... .......... .......... .......... 36% 33.5M 2s
    ##  48800K .......... .......... .......... .......... .......... 36% 72.8M 2s
    ##  48850K .......... .......... .......... .......... .......... 36%  117M 2s
    ##  48900K .......... .......... .......... .......... .......... 36% 35.9M 2s
    ##  48950K .......... .......... .......... .......... .......... 36%  130M 2s
    ##  49000K .......... .......... .......... .......... .......... 36% 83.2M 2s
    ##  49050K .......... .......... .......... .......... .......... 36% 76.7M 2s
    ##  49100K .......... .......... .......... .......... .......... 36%  114M 2s
    ##  49150K .......... .......... .......... .......... .......... 36% 21.0M 2s
    ##  49200K .......... .......... .......... .......... .......... 36% 39.0M 2s
    ##  49250K .......... .......... .......... .......... .......... 36% 49.9M 2s
    ##  49300K .......... .......... .......... .......... .......... 36% 38.2M 2s
    ##  49350K .......... .......... .......... .......... .......... 36% 48.3M 2s
    ##  49400K .......... .......... .......... .......... .......... 36% 53.4M 2s
    ##  49450K .......... .......... .......... .......... .......... 36% 51.3M 2s
    ##  49500K .......... .......... .......... .......... .......... 36% 55.8M 2s
    ##  49550K .......... .......... .......... .......... .......... 36% 6.92M 2s
    ##  49600K .......... .......... .......... .......... .......... 36% 80.1M 2s
    ##  49650K .......... .......... .......... .......... .......... 36% 94.4M 2s
    ##  49700K .......... .......... .......... .......... .......... 36% 70.1M 2s
    ##  49750K .......... .......... .......... .......... .......... 36%  101M 2s
    ##  49800K .......... .......... .......... .......... .......... 36% 81.8M 2s
    ##  49850K .......... .......... .......... .......... .......... 37% 51.3M 2s
    ##  49900K .......... .......... .......... .......... .......... 37% 53.8M 2s
    ##  49950K .......... .......... .......... .......... .......... 37% 45.5M 2s
    ##  50000K .......... .......... .......... .......... .......... 37% 96.8M 2s
    ##  50050K .......... .......... .......... .......... .......... 37%  108M 2s
    ##  50100K .......... .......... .......... .......... .......... 37% 62.8M 2s
    ##  50150K .......... .......... .......... .......... .......... 37% 88.9M 2s
    ##  50200K .......... .......... .......... .......... .......... 37%  101M 2s
    ##  50250K .......... .......... .......... .......... .......... 37% 95.9M 2s
    ##  50300K .......... .......... .......... .......... .......... 37% 48.2M 2s
    ##  50350K .......... .......... .......... .......... .......... 37%  106M 2s
    ##  50400K .......... .......... .......... .......... .......... 37% 99.6M 2s
    ##  50450K .......... .......... .......... .......... .......... 37% 18.3M 2s
    ##  50500K .......... .......... .......... .......... .......... 37% 50.9M 2s
    ##  50550K .......... .......... .......... .......... .......... 37% 35.5M 2s
    ##  50600K .......... .......... .......... .......... .......... 37%  104M 2s
    ##  50650K .......... .......... .......... .......... .......... 37% 83.2M 2s
    ##  50700K .......... .......... .......... .......... .......... 37%  108M 2s
    ##  50750K .......... .......... .......... .......... .......... 37% 52.0M 2s
    ##  50800K .......... .......... .......... .......... .......... 37%  103M 2s
    ##  50850K .......... .......... .......... .......... .......... 37% 33.8M 2s
    ##  50900K .......... .......... .......... .......... .......... 37%  104M 2s
    ##  50950K .......... .......... .......... .......... .......... 37%  117M 2s
    ##  51000K .......... .......... .......... .......... .......... 37% 63.9M 2s
    ##  51050K .......... .......... .......... .......... .......... 37% 55.8M 2s
    ##  51100K .......... .......... .......... .......... .......... 37% 78.8M 2s
    ##  51150K .......... .......... .......... .......... .......... 37% 55.6M 2s
    ##  51200K .......... .......... .......... .......... .......... 38%  110M 2s
    ##  51250K .......... .......... .......... .......... .......... 38% 59.9M 2s
    ##  51300K .......... .......... .......... .......... .......... 38% 52.9M 2s
    ##  51350K .......... .......... .......... .......... .......... 38% 71.2M 2s
    ##  51400K .......... .......... .......... .......... .......... 38% 73.3M 2s
    ##  51450K .......... .......... .......... .......... .......... 38%  106M 2s
    ##  51500K .......... .......... .......... .......... .......... 38% 69.0M 2s
    ##  51550K .......... .......... .......... .......... .......... 38% 55.4M 2s
    ##  51600K .......... .......... .......... .......... .......... 38% 55.3M 2s
    ##  51650K .......... .......... .......... .......... .......... 38% 69.2M 2s
    ##  51700K .......... .......... .......... .......... .......... 38%  106M 2s
    ##  51750K .......... .......... .......... .......... .......... 38% 60.0M 2s
    ##  51800K .......... .......... .......... .......... .......... 38% 55.6M 2s
    ##  51850K .......... .......... .......... .......... .......... 38% 49.0M 2s
    ##  51900K .......... .......... .......... .......... .......... 38% 67.8M 2s
    ##  51950K .......... .......... .......... .......... .......... 38% 71.5M 2s
    ##  52000K .......... .......... .......... .......... .......... 38%  112M 2s
    ##  52050K .......... .......... .......... .......... .......... 38% 59.6M 2s
    ##  52100K .......... .......... .......... .......... .......... 38% 54.5M 2s
    ##  52150K .......... .......... .......... .......... .......... 38%  119M 2s
    ##  52200K .......... .......... .......... .......... .......... 38% 29.7M 2s
    ##  52250K .......... .......... .......... .......... .......... 38% 65.7M 2s
    ##  52300K .......... .......... .......... .......... .......... 38% 71.9M 2s
    ##  52350K .......... .......... .......... .......... .......... 38%  109M 2s
    ##  52400K .......... .......... .......... .......... .......... 38% 32.9M 2s
    ##  52450K .......... .......... .......... .......... .......... 38% 96.1M 2s
    ##  52500K .......... .......... .......... .......... .......... 39% 63.8M 2s
    ##  52550K .......... .......... .......... .......... .......... 39%  126M 2s
    ##  52600K .......... .......... .......... .......... .......... 39%  106M 2s
    ##  52650K .......... .......... .......... .......... .......... 39% 30.0M 2s
    ##  52700K .......... .......... .......... .......... .......... 39% 54.6M 2s
    ##  52750K .......... .......... .......... .......... .......... 39% 86.6M 2s
    ##  52800K .......... .......... .......... .......... .......... 39%  117M 2s
    ##  52850K .......... .......... .......... .......... .......... 39%  135M 2s
    ##  52900K .......... .......... .......... .......... .......... 39% 53.2M 2s
    ##  52950K .......... .......... .......... .......... .......... 39% 44.5M 2s
    ##  53000K .......... .......... .......... .......... .......... 39% 57.7M 2s
    ##  53050K .......... .......... .......... .......... .......... 39%  102M 2s
    ##  53100K .......... .......... .......... .......... .......... 39%  116M 2s
    ##  53150K .......... .......... .......... .......... .......... 39% 82.5M 2s
    ##  53200K .......... .......... .......... .......... .......... 39% 50.0M 2s
    ##  53250K .......... .......... .......... .......... .......... 39% 82.9M 2s
    ##  53300K .......... .......... .......... .......... .......... 39% 60.0M 2s
    ##  53350K .......... .......... .......... .......... .......... 39%  110M 2s
    ##  53400K .......... .......... .......... .......... .......... 39% 62.0M 2s
    ##  53450K .......... .......... .......... .......... .......... 39% 45.2M 2s
    ##  53500K .......... .......... .......... .......... .......... 39% 84.2M 2s
    ##  53550K .......... .......... .......... .......... .......... 39% 80.2M 2s
    ##  53600K .......... .......... .......... .......... .......... 39% 87.2M 2s
    ##  53650K .......... .......... .......... .......... .......... 39% 76.9M 2s
    ##  53700K .......... .......... .......... .......... .......... 39% 58.5M 2s
    ##  53750K .......... .......... .......... .......... .......... 39% 62.8M 2s
    ##  53800K .......... .......... .......... .......... .......... 39%  114M 2s
    ##  53850K .......... .......... .......... .......... .......... 40% 79.5M 2s
    ##  53900K .......... .......... .......... .......... .......... 40% 67.1M 2s
    ##  53950K .......... .......... .......... .......... .......... 40% 89.5M 2s
    ##  54000K .......... .......... .......... .......... .......... 40% 64.7M 2s
    ##  54050K .......... .......... .......... .......... .......... 40% 56.8M 2s
    ##  54100K .......... .......... .......... .......... .......... 40% 82.3M 2s
    ##  54150K .......... .......... .......... .......... .......... 40% 72.0M 2s
    ##  54200K .......... .......... .......... .......... .......... 40% 74.9M 2s
    ##  54250K .......... .......... .......... .......... .......... 40% 83.5M 2s
    ##  54300K .......... .......... .......... .......... .......... 40% 65.4M 2s
    ##  54350K .......... .......... .......... .......... .......... 40% 64.2M 2s
    ##  54400K .......... .......... .......... .......... .......... 40% 63.9M 2s
    ##  54450K .......... .......... .......... .......... .......... 40%  108M 2s
    ##  54500K .......... .......... .......... .......... .......... 40% 66.1M 2s
    ##  54550K .......... .......... .......... .......... .......... 40% 60.3M 2s
    ##  54600K .......... .......... .......... .......... .......... 40% 53.4M 2s
    ##  54650K .......... .......... .......... .......... .......... 40% 54.5M 2s
    ##  54700K .......... .......... .......... .......... .......... 40% 74.5M 2s
    ##  54750K .......... .......... .......... .......... .......... 40%  145M 2s
    ##  54800K .......... .......... .......... .......... .......... 40% 40.3M 2s
    ##  54850K .......... .......... .......... .......... .......... 40%  119M 2s
    ##  54900K .......... .......... .......... .......... .......... 40% 76.8M 2s
    ##  54950K .......... .......... .......... .......... .......... 40% 71.3M 2s
    ##  55000K .......... .......... .......... .......... .......... 40%  105M 2s
    ##  55050K .......... .......... .......... .......... .......... 40% 19.0M 2s
    ##  55100K .......... .......... .......... .......... .......... 40% 64.9M 2s
    ##  55150K .......... .......... .......... .......... .......... 40% 78.6M 2s
    ##  55200K .......... .......... .......... .......... .......... 41%  135M 2s
    ##  55250K .......... .......... .......... .......... .......... 41%  118M 2s
    ##  55300K .......... .......... .......... .......... .......... 41% 26.2M 2s
    ##  55350K .......... .......... .......... .......... .......... 41% 41.3M 2s
    ##  55400K .......... .......... .......... .......... .......... 41%  105M 2s
    ##  55450K .......... .......... .......... .......... .......... 41% 93.1M 2s
    ##  55500K .......... .......... .......... .......... .......... 41%  126M 2s
    ##  55550K .......... .......... .......... .......... .......... 41%  106M 2s
    ##  55600K .......... .......... .......... .......... .......... 41% 73.0M 2s
    ##  55650K .......... .......... .......... .......... .......... 41%  143M 2s
    ##  55700K .......... .......... .......... .......... .......... 41% 25.7M 2s
    ##  55750K .......... .......... .......... .......... .......... 41%  113M 2s
    ##  55800K .......... .......... .......... .......... .......... 41% 75.3M 2s
    ##  55850K .......... .......... .......... .......... .......... 41% 79.2M 2s
    ##  55900K .......... .......... .......... .......... .......... 41%  129M 2s
    ##  55950K .......... .......... .......... .......... .......... 41%  140M 2s
    ##  56000K .......... .......... .......... .......... .......... 41% 29.0M 2s
    ##  56050K .......... .......... .......... .......... .......... 41% 53.2M 2s
    ##  56100K .......... .......... .......... .......... .......... 41% 76.0M 2s
    ##  56150K .......... .......... .......... .......... .......... 41% 79.2M 2s
    ##  56200K .......... .......... .......... .......... .......... 41% 92.5M 2s
    ##  56250K .......... .......... .......... .......... .......... 41% 76.8M 2s
    ##  56300K .......... .......... .......... .......... .......... 41% 37.7M 2s
    ##  56350K .......... .......... .......... .......... .......... 41% 80.8M 2s
    ##  56400K .......... .......... .......... .......... .......... 41%  107M 2s
    ##  56450K .......... .......... .......... .......... .......... 41% 70.9M 2s
    ##  56500K .......... .......... .......... .......... .......... 41%  116M 2s
    ##  56550K .......... .......... .......... .......... .......... 42% 25.8M 2s
    ##  56600K .......... .......... .......... .......... .......... 42% 88.5M 2s
    ##  56650K .......... .......... .......... .......... .......... 42% 49.2M 2s
    ##  56700K .......... .......... .......... .......... .......... 42% 94.7M 2s
    ##  56750K .......... .......... .......... .......... .......... 42%  137M 2s
    ##  56800K .......... .......... .......... .......... .......... 42% 63.0M 2s
    ##  56850K .......... .......... .......... .......... .......... 42% 37.8M 2s
    ##  56900K .......... .......... .......... .......... .......... 42% 99.9M 2s
    ##  56950K .......... .......... .......... .......... .......... 42% 56.9M 2s
    ##  57000K .......... .......... .......... .......... .......... 42% 65.9M 2s
    ##  57050K .......... .......... .......... .......... .......... 42%  106M 2s
    ##  57100K .......... .......... .......... .......... .......... 42%  126M 2s
    ##  57150K .......... .......... .......... .......... .......... 42% 63.2M 2s
    ##  57200K .......... .......... .......... .......... .......... 42% 45.5M 2s
    ##  57250K .......... .......... .......... .......... .......... 42% 72.7M 2s
    ##  57300K .......... .......... .......... .......... .......... 42% 59.7M 2s
    ##  57350K .......... .......... .......... .......... .......... 42% 87.7M 2s
    ##  57400K .......... .......... .......... .......... .......... 42%  106M 2s
    ##  57450K .......... .......... .......... .......... .......... 42% 56.5M 2s
    ##  57500K .......... .......... .......... .......... .......... 42% 77.7M 2s
    ##  57550K .......... .......... .......... .......... .......... 42%  125M 2s
    ##  57600K .......... .......... .......... .......... .......... 42% 34.7M 2s
    ##  57650K .......... .......... .......... .......... .......... 42%  131M 2s
    ##  57700K .......... .......... .......... .......... .......... 42% 62.8M 2s
    ##  57750K .......... .......... .......... .......... .......... 42% 66.7M 2s
    ##  57800K .......... .......... .......... .......... .......... 42%  113M 2s
    ##  57850K .......... .......... .......... .......... .......... 42% 29.4M 2s
    ##  57900K .......... .......... .......... .......... .......... 43% 61.5M 2s
    ##  57950K .......... .......... .......... .......... .......... 43% 63.7M 2s
    ##  58000K .......... .......... .......... .......... .......... 43% 61.6M 2s
    ##  58050K .......... .......... .......... .......... .......... 43% 79.8M 2s
    ##  58100K .......... .......... .......... .......... .......... 43% 63.8M 2s
    ##  58150K .......... .......... .......... .......... .......... 43%  145M 2s
    ##  58200K .......... .......... .......... .......... .......... 43% 48.9M 2s
    ##  58250K .......... .......... .......... .......... .......... 43% 70.5M 2s
    ##  58300K .......... .......... .......... .......... .......... 43% 84.7M 2s
    ##  58350K .......... .......... .......... .......... .......... 43% 98.1M 2s
    ##  58400K .......... .......... .......... .......... .......... 43%  119M 2s
    ##  58450K .......... .......... .......... .......... .......... 43% 47.2M 2s
    ##  58500K .......... .......... .......... .......... .......... 43% 53.6M 2s
    ##  58550K .......... .......... .......... .......... .......... 43% 35.4M 2s
    ##  58600K .......... .......... .......... .......... .......... 43% 56.6M 2s
    ##  58650K .......... .......... .......... .......... .......... 43% 86.4M 2s
    ##  58700K .......... .......... .......... .......... .......... 43% 86.9M 2s
    ##  58750K .......... .......... .......... .......... .......... 43% 91.7M 2s
    ##  58800K .......... .......... .......... .......... .......... 43% 76.2M 2s
    ##  58850K .......... .......... .......... .......... .......... 43% 69.5M 2s
    ##  58900K .......... .......... .......... .......... .......... 43% 67.3M 2s
    ##  58950K .......... .......... .......... .......... .......... 43% 97.9M 2s
    ##  59000K .......... .......... .......... .......... .......... 43% 33.9M 2s
    ##  59050K .......... .......... .......... .......... .......... 43% 78.3M 2s
    ##  59100K .......... .......... .......... .......... .......... 43% 91.0M 2s
    ##  59150K .......... .......... .......... .......... .......... 43% 71.0M 2s
    ##  59200K .......... .......... .......... .......... .......... 43% 74.2M 2s
    ##  59250K .......... .......... .......... .......... .......... 44% 76.7M 2s
    ##  59300K .......... .......... .......... .......... .......... 44% 51.0M 2s
    ##  59350K .......... .......... .......... .......... .......... 44% 79.5M 2s
    ##  59400K .......... .......... .......... .......... .......... 44% 69.5M 2s
    ##  59450K .......... .......... .......... .......... .......... 44% 84.4M 2s
    ##  59500K .......... .......... .......... .......... .......... 44% 68.3M 2s
    ##  59550K .......... .......... .......... .......... .......... 44% 19.9M 2s
    ##  59600K .......... .......... .......... .......... .......... 44% 82.3M 2s
    ##  59650K .......... .......... .......... .......... .......... 44% 78.4M 2s
    ##  59700K .......... .......... .......... .......... .......... 44% 86.4M 2s
    ##  59750K .......... .......... .......... .......... .......... 44% 78.5M 2s
    ##  59800K .......... .......... .......... .......... .......... 44% 80.0M 2s
    ##  59850K .......... .......... .......... .......... .......... 44%  102M 2s
    ##  59900K .......... .......... .......... .......... .......... 44% 22.2M 2s
    ##  59950K .......... .......... .......... .......... .......... 44%  104M 2s
    ##  60000K .......... .......... .......... .......... .......... 44% 77.8M 2s
    ##  60050K .......... .......... .......... .......... .......... 44% 84.0M 2s
    ##  60100K .......... .......... .......... .......... .......... 44% 86.6M 2s
    ##  60150K .......... .......... .......... .......... .......... 44% 77.4M 2s
    ##  60200K .......... .......... .......... .......... .......... 44% 31.3M 2s
    ##  60250K .......... .......... .......... .......... .......... 44% 79.5M 2s
    ##  60300K .......... .......... .......... .......... .......... 44% 69.5M 2s
    ##  60350K .......... .......... .......... .......... .......... 44% 92.4M 2s
    ##  60400K .......... .......... .......... .......... .......... 44% 87.2M 2s
    ##  60450K .......... .......... .......... .......... .......... 44%  102M 2s
    ##  60500K .......... .......... .......... .......... .......... 44% 74.0M 2s
    ##  60550K .......... .......... .......... .......... .......... 44% 60.7M 2s
    ##  60600K .......... .......... .......... .......... .......... 45% 67.5M 2s
    ##  60650K .......... .......... .......... .......... .......... 45% 81.7M 2s
    ##  60700K .......... .......... .......... .......... .......... 45% 74.0M 2s
    ##  60750K .......... .......... .......... .......... .......... 45% 97.5M 2s
    ##  60800K .......... .......... .......... .......... .......... 45% 65.4M 2s
    ##  60850K .......... .......... .......... .......... .......... 45% 85.2M 2s
    ##  60900K .......... .......... .......... .......... .......... 45% 71.4M 2s
    ##  60950K .......... .......... .......... .......... .......... 45% 76.6M 2s
    ##  61000K .......... .......... .......... .......... .......... 45% 95.9M 2s
    ##  61050K .......... .......... .......... .......... .......... 45% 89.5M 2s
    ##  61100K .......... .......... .......... .......... .......... 45% 69.0M 2s
    ##  61150K .......... .......... .......... .......... .......... 45% 71.6M 2s
    ##  61200K .......... .......... .......... .......... .......... 45% 80.8M 2s
    ##  61250K .......... .......... .......... .......... .......... 45% 99.8M 2s
    ##  61300K .......... .......... .......... .......... .......... 45% 73.2M 2s
    ##  61350K .......... .......... .......... .......... .......... 45% 85.7M 2s
    ##  61400K .......... .......... .......... .......... .......... 45% 58.7M 2s
    ##  61450K .......... .......... .......... .......... .......... 45% 70.7M 2s
    ##  61500K .......... .......... .......... .......... .......... 45% 78.9M 2s
    ##  61550K .......... .......... .......... .......... .......... 45% 80.1M 2s
    ##  61600K .......... .......... .......... .......... .......... 45% 70.7M 2s
    ##  61650K .......... .......... .......... .......... .......... 45%  116M 2s
    ##  61700K .......... .......... .......... .......... .......... 45% 72.5M 2s
    ##  61750K .......... .......... .......... .......... .......... 45% 84.1M 2s
    ##  61800K .......... .......... .......... .......... .......... 45% 87.8M 2s
    ##  61850K .......... .......... .......... .......... .......... 45% 79.1M 2s
    ##  61900K .......... .......... .......... .......... .......... 45% 94.2M 2s
    ##  61950K .......... .......... .......... .......... .......... 46% 88.3M 2s
    ##  62000K .......... .......... .......... .......... .......... 46% 88.7M 2s
    ##  62050K .......... .......... .......... .......... .......... 46% 78.6M 2s
    ##  62100K .......... .......... .......... .......... .......... 46% 77.8M 2s
    ##  62150K .......... .......... .......... .......... .......... 46% 66.9M 2s
    ##  62200K .......... .......... .......... .......... .......... 46% 93.9M 2s
    ##  62250K .......... .......... .......... .......... .......... 46% 89.9M 2s
    ##  62300K .......... .......... .......... .......... .......... 46% 86.7M 2s
    ##  62350K .......... .......... .......... .......... .......... 46% 56.9M 2s
    ##  62400K .......... .......... .......... .......... .......... 46% 78.2M 2s
    ##  62450K .......... .......... .......... .......... .......... 46% 68.2M 2s
    ##  62500K .......... .......... .......... .......... .......... 46%  103M 2s
    ##  62550K .......... .......... .......... .......... .......... 46%  107M 2s
    ##  62600K .......... .......... .......... .......... .......... 46% 16.8M 2s
    ##  62650K .......... .......... .......... .......... .......... 46% 80.2M 2s
    ##  62700K .......... .......... .......... .......... .......... 46% 45.6M 2s
    ##  62750K .......... .......... .......... .......... .......... 46% 69.8M 2s
    ##  62800K .......... .......... .......... .......... .......... 46% 64.2M 2s
    ##  62850K .......... .......... .......... .......... .......... 46% 90.5M 2s
    ##  62900K .......... .......... .......... .......... .......... 46% 74.5M 2s
    ##  62950K .......... .......... .......... .......... .......... 46% 64.9M 2s
    ##  63000K .......... .......... .......... .......... .......... 46% 70.9M 2s
    ##  63050K .......... .......... .......... .......... .......... 46% 89.8M 2s
    ##  63100K .......... .......... .......... .......... .......... 46% 53.9M 2s
    ##  63150K .......... .......... .......... .......... .......... 46% 85.8M 2s
    ##  63200K .......... .......... .......... .......... .......... 46% 75.3M 2s
    ##  63250K .......... .......... .......... .......... .......... 46% 97.4M 2s
    ##  63300K .......... .......... .......... .......... .......... 47% 86.4M 2s
    ##  63350K .......... .......... .......... .......... .......... 47% 86.9M 2s
    ##  63400K .......... .......... .......... .......... .......... 47% 87.0M 2s
    ##  63450K .......... .......... .......... .......... .......... 47% 79.2M 2s
    ##  63500K .......... .......... .......... .......... .......... 47% 61.6M 2s
    ##  63550K .......... .......... .......... .......... .......... 47% 74.5M 2s
    ##  63600K .......... .......... .......... .......... .......... 47% 86.6M 2s
    ##  63650K .......... .......... .......... .......... .......... 47% 81.3M 2s
    ##  63700K .......... .......... .......... .......... .......... 47% 86.9M 2s
    ##  63750K .......... .......... .......... .......... .......... 47%  106M 2s
    ##  63800K .......... .......... .......... .......... .......... 47% 44.8M 2s
    ##  63850K .......... .......... .......... .......... .......... 47% 84.1M 2s
    ##  63900K .......... .......... .......... .......... .......... 47%  100M 2s
    ##  63950K .......... .......... .......... .......... .......... 47% 94.6M 2s
    ##  64000K .......... .......... .......... .......... .......... 47% 81.0M 2s
    ##  64050K .......... .......... .......... .......... .......... 47% 99.7M 2s
    ##  64100K .......... .......... .......... .......... .......... 47% 93.6M 2s
    ##  64150K .......... .......... .......... .......... .......... 47% 79.1M 2s
    ##  64200K .......... .......... .......... .......... .......... 47% 64.2M 2s
    ##  64250K .......... .......... .......... .......... .......... 47% 48.2M 2s
    ##  64300K .......... .......... .......... .......... .......... 47% 91.5M 2s
    ##  64350K .......... .......... .......... .......... .......... 47%  116M 2s
    ##  64400K .......... .......... .......... .......... .......... 47% 73.0M 2s
    ##  64450K .......... .......... .......... .......... .......... 47% 32.8M 2s
    ##  64500K .......... .......... .......... .......... .......... 47% 93.2M 2s
    ##  64550K .......... .......... .......... .......... .......... 47%  105M 2s
    ##  64600K .......... .......... .......... .......... .......... 47% 71.2M 2s
    ##  64650K .......... .......... .......... .......... .......... 48%  104M 2s
    ##  64700K .......... .......... .......... .......... .......... 48% 67.8M 2s
    ##  64750K .......... .......... .......... .......... .......... 48% 97.4M 1s
    ##  64800K .......... .......... .......... .......... .......... 48% 76.7M 1s
    ##  64850K .......... .......... .......... .......... .......... 48%  107M 1s
    ##  64900K .......... .......... .......... .......... .......... 48% 65.1M 1s
    ##  64950K .......... .......... .......... .......... .......... 48% 85.2M 1s
    ##  65000K .......... .......... .......... .......... .......... 48% 91.4M 1s
    ##  65050K .......... .......... .......... .......... .......... 48%  103M 1s
    ##  65100K .......... .......... .......... .......... .......... 48% 79.0M 1s
    ##  65150K .......... .......... .......... .......... .......... 48% 69.6M 1s
    ##  65200K .......... .......... .......... .......... .......... 48% 66.8M 1s
    ##  65250K .......... .......... .......... .......... .......... 48% 85.4M 1s
    ##  65300K .......... .......... .......... .......... .......... 48% 85.6M 1s
    ##  65350K .......... .......... .......... .......... .......... 48% 87.8M 1s
    ##  65400K .......... .......... .......... .......... .......... 48% 79.8M 1s
    ##  65450K .......... .......... .......... .......... .......... 48% 98.5M 1s
    ##  65500K .......... .......... .......... .......... .......... 48% 44.8M 1s
    ##  65550K .......... .......... .......... .......... .......... 48% 91.6M 1s
    ##  65600K .......... .......... .......... .......... .......... 48% 88.1M 1s
    ##  65650K .......... .......... .......... .......... .......... 48% 86.3M 1s
    ##  65700K .......... .......... .......... .......... .......... 48% 79.1M 1s
    ##  65750K .......... .......... .......... .......... .......... 48%  124M 1s
    ##  65800K .......... .......... .......... .......... .......... 48% 89.3M 1s
    ##  65850K .......... .......... .......... .......... .......... 48% 83.0M 1s
    ##  65900K .......... .......... .......... .......... .......... 48% 69.8M 1s
    ##  65950K .......... .......... .......... .......... .......... 48% 75.4M 1s
    ##  66000K .......... .......... .......... .......... .......... 49% 82.8M 1s
    ##  66050K .......... .......... .......... .......... .......... 49%  101M 1s
    ##  66100K .......... .......... .......... .......... .......... 49% 83.4M 1s
    ##  66150K .......... .......... .......... .......... .......... 49% 59.7M 1s
    ##  66200K .......... .......... .......... .......... .......... 49% 83.4M 1s
    ##  66250K .......... .......... .......... .......... .......... 49% 77.0M 1s
    ##  66300K .......... .......... .......... .......... .......... 49%  110M 1s
    ##  66350K .......... .......... .......... .......... .......... 49% 64.8M 1s
    ##  66400K .......... .......... .......... .......... .......... 49% 89.0M 1s
    ##  66450K .......... .......... .......... .......... .......... 49%  108M 1s
    ##  66500K .......... .......... .......... .......... .......... 49% 82.8M 1s
    ##  66550K .......... .......... .......... .......... .......... 49% 53.2M 1s
    ##  66600K .......... .......... .......... .......... .......... 49% 85.7M 1s
    ##  66650K .......... .......... .......... .......... .......... 49%  108M 1s
    ##  66700K .......... .......... .......... .......... .......... 49% 99.5M 1s
    ##  66750K .......... .......... .......... .......... .......... 49%  125M 1s
    ##  66800K .......... .......... .......... .......... .......... 49% 23.1M 1s
    ##  66850K .......... .......... .......... .......... .......... 49% 93.8M 1s
    ##  66900K .......... .......... .......... .......... .......... 49% 87.0M 1s
    ##  66950K .......... .......... .......... .......... .......... 49%  105M 1s
    ##  67000K .......... .......... .......... .......... .......... 49% 88.7M 1s
    ##  67050K .......... .......... .......... .......... .......... 49%  110M 1s
    ##  67100K .......... .......... .......... .......... .......... 49% 32.9M 1s
    ##  67150K .......... .......... .......... .......... .......... 49%  115M 1s
    ##  67200K .......... .......... .......... .......... .......... 49% 90.6M 1s
    ##  67250K .......... .......... .......... .......... .......... 49% 41.3M 1s
    ##  67300K .......... .......... .......... .......... .......... 49% 71.8M 1s
    ##  67350K .......... .......... .......... .......... .......... 50% 94.3M 1s
    ##  67400K .......... .......... .......... .......... .......... 50%  102M 1s
    ##  67450K .......... .......... .......... .......... .......... 50%  136M 1s
    ##  67500K .......... .......... .......... .......... .......... 50%  109M 1s
    ##  67550K .......... .......... .......... .......... .......... 50% 33.6M 1s
    ##  67600K .......... .......... .......... .......... .......... 50% 69.5M 1s
    ##  67650K .......... .......... .......... .......... .......... 50%  114M 1s
    ##  67700K .......... .......... .......... .......... .......... 50% 79.5M 1s
    ##  67750K .......... .......... .......... .......... .......... 50% 94.5M 1s
    ##  67800K .......... .......... .......... .......... .......... 50% 87.3M 1s
    ##  67850K .......... .......... .......... .......... .......... 50%  108M 1s
    ##  67900K .......... .......... .......... .......... .......... 50% 44.1M 1s
    ##  67950K .......... .......... .......... .......... .......... 50% 95.2M 1s
    ##  68000K .......... .......... .......... .......... .......... 50% 94.3M 1s
    ##  68050K .......... .......... .......... .......... .......... 50% 36.3M 1s
    ##  68100K .......... .......... .......... .......... .......... 50% 63.1M 1s
    ##  68150K .......... .......... .......... .......... .......... 50%  116M 1s
    ##  68200K .......... .......... .......... .......... .......... 50% 87.9M 1s
    ##  68250K .......... .......... .......... .......... .......... 50%  102M 1s
    ##  68300K .......... .......... .......... .......... .......... 50%  111M 1s
    ##  68350K .......... .......... .......... .......... .......... 50%  125M 1s
    ##  68400K .......... .......... .......... .......... .......... 50% 83.2M 1s
    ##  68450K .......... .......... .......... .......... .......... 50%  107M 1s
    ##  68500K .......... .......... .......... .......... .......... 50% 44.0M 1s
    ##  68550K .......... .......... .......... .......... .......... 50% 50.3M 1s
    ##  68600K .......... .......... .......... .......... .......... 50% 97.3M 1s
    ##  68650K .......... .......... .......... .......... .......... 50%  108M 1s
    ##  68700K .......... .......... .......... .......... .......... 51% 90.9M 1s
    ##  68750K .......... .......... .......... .......... .......... 51% 87.3M 1s
    ##  68800K .......... .......... .......... .......... .......... 51% 92.6M 1s
    ##  68850K .......... .......... .......... .......... .......... 51%  118M 1s
    ##  68900K .......... .......... .......... .......... .......... 51% 35.0M 1s
    ##  68950K .......... .......... .......... .......... .......... 51% 97.7M 1s
    ##  69000K .......... .......... .......... .......... .......... 51% 80.9M 1s
    ##  69050K .......... .......... .......... .......... .......... 51% 86.3M 1s
    ##  69100K .......... .......... .......... .......... .......... 51% 73.8M 1s
    ##  69150K .......... .......... .......... .......... .......... 51%  132M 1s
    ##  69200K .......... .......... .......... .......... .......... 51% 84.5M 1s
    ##  69250K .......... .......... .......... .......... .......... 51%  110M 1s
    ##  69300K .......... .......... .......... .......... .......... 51% 95.8M 1s
    ##  69350K .......... .......... .......... .......... .......... 51% 99.0M 1s
    ##  69400K .......... .......... .......... .......... .......... 51% 87.8M 1s
    ##  69450K .......... .......... .......... .......... .......... 51% 83.1M 1s
    ##  69500K .......... .......... .......... .......... .......... 51% 20.0M 1s
    ##  69550K .......... .......... .......... .......... .......... 51%  113M 1s
    ##  69600K .......... .......... .......... .......... .......... 51%  105M 1s
    ##  69650K .......... .......... .......... .......... .......... 51%  111M 1s
    ##  69700K .......... .......... .......... .......... .......... 51% 84.7M 1s
    ##  69750K .......... .......... .......... .......... .......... 51%  138M 1s
    ##  69800K .......... .......... .......... .......... .......... 51%  130M 1s
    ##  69850K .......... .......... .......... .......... .......... 51% 17.1M 1s
    ##  69900K .......... .......... .......... .......... .......... 51% 97.7M 1s
    ##  69950K .......... .......... .......... .......... .......... 51% 51.6M 1s
    ##  70000K .......... .......... .......... .......... .......... 51% 87.5M 1s
    ##  70050K .......... .......... .......... .......... .......... 52% 97.9M 1s
    ##  70100K .......... .......... .......... .......... .......... 52%  124M 1s
    ##  70150K .......... .......... .......... .......... .......... 52%  109M 1s
    ##  70200K .......... .......... .......... .......... .......... 52% 53.4M 1s
    ##  70250K .......... .......... .......... .......... .......... 52%  121M 1s
    ##  70300K .......... .......... .......... .......... .......... 52% 26.7M 1s
    ##  70350K .......... .......... .......... .......... .......... 52%  137M 1s
    ##  70400K .......... .......... .......... .......... .......... 52% 42.3M 1s
    ##  70450K .......... .......... .......... .......... .......... 52% 70.2M 1s
    ##  70500K .......... .......... .......... .......... .......... 52% 80.9M 1s
    ##  70550K .......... .......... .......... .......... .......... 52% 92.8M 1s
    ##  70600K .......... .......... .......... .......... .......... 52% 22.0M 1s
    ##  70650K .......... .......... .......... .......... .......... 52% 71.2M 1s
    ##  70700K .......... .......... .......... .......... .......... 52% 68.3M 1s
    ##  70750K .......... .......... .......... .......... .......... 52% 59.5M 1s
    ##  70800K .......... .......... .......... .......... .......... 52% 93.3M 1s
    ##  70850K .......... .......... .......... .......... .......... 52% 82.3M 1s
    ##  70900K .......... .......... .......... .......... .......... 52% 79.3M 1s
    ##  70950K .......... .......... .......... .......... .......... 52% 79.7M 1s
    ##  71000K .......... .......... .......... .......... .......... 52% 91.2M 1s
    ##  71050K .......... .......... .......... .......... .......... 52% 42.7M 1s
    ##  71100K .......... .......... .......... .......... .......... 52% 58.5M 1s
    ##  71150K .......... .......... .......... .......... .......... 52% 71.5M 1s
    ##  71200K .......... .......... .......... .......... .......... 52% 74.0M 1s
    ##  71250K .......... .......... .......... .......... .......... 52% 75.6M 1s
    ##  71300K .......... .......... .......... .......... .......... 52% 90.5M 1s
    ##  71350K .......... .......... .......... .......... .......... 52% 67.4M 1s
    ##  71400K .......... .......... .......... .......... .......... 53% 44.9M 1s
    ##  71450K .......... .......... .......... .......... .......... 53% 52.3M 1s
    ##  71500K .......... .......... .......... .......... .......... 53% 82.3M 1s
    ##  71550K .......... .......... .......... .......... .......... 53% 71.9M 1s
    ##  71600K .......... .......... .......... .......... .......... 53% 71.9M 1s
    ##  71650K .......... .......... .......... .......... .......... 53% 79.3M 1s
    ##  71700K .......... .......... .......... .......... .......... 53% 88.2M 1s
    ##  71750K .......... .......... .......... .......... .......... 53% 88.5M 1s
    ##  71800K .......... .......... .......... .......... .......... 53% 77.4M 1s
    ##  71850K .......... .......... .......... .......... .......... 53% 65.8M 1s
    ##  71900K .......... .......... .......... .......... .......... 53% 68.5M 1s
    ##  71950K .......... .......... .......... .......... .......... 53% 77.8M 1s
    ##  72000K .......... .......... .......... .......... .......... 53% 87.7M 1s
    ##  72050K .......... .......... .......... .......... .......... 53% 77.7M 1s
    ##  72100K .......... .......... .......... .......... .......... 53% 69.5M 1s
    ##  72150K .......... .......... .......... .......... .......... 53% 80.2M 1s
    ##  72200K .......... .......... .......... .......... .......... 53% 76.4M 1s
    ##  72250K .......... .......... .......... .......... .......... 53% 75.5M 1s
    ##  72300K .......... .......... .......... .......... .......... 53% 64.6M 1s
    ##  72350K .......... .......... .......... .......... .......... 53% 87.4M 1s
    ##  72400K .......... .......... .......... .......... .......... 53% 91.4M 1s
    ##  72450K .......... .......... .......... .......... .......... 53% 94.1M 1s
    ##  72500K .......... .......... .......... .......... .......... 53% 24.2M 1s
    ##  72550K .......... .......... .......... .......... .......... 53% 86.7M 1s
    ##  72600K .......... .......... .......... .......... .......... 53% 86.8M 1s
    ##  72650K .......... .......... .......... .......... .......... 53% 79.1M 1s
    ##  72700K .......... .......... .......... .......... .......... 53% 91.0M 1s
    ##  72750K .......... .......... .......... .......... .......... 54%  110M 1s
    ##  72800K .......... .......... .......... .......... .......... 54% 25.4M 1s
    ##  72850K .......... .......... .......... .......... .......... 54% 88.7M 1s
    ##  72900K .......... .......... .......... .......... .......... 54% 61.1M 1s
    ##  72950K .......... .......... .......... .......... .......... 54% 64.9M 1s
    ##  73000K .......... .......... .......... .......... .......... 54% 92.7M 1s
    ##  73050K .......... .......... .......... .......... .......... 54% 94.0M 1s
    ##  73100K .......... .......... .......... .......... .......... 54% 61.9M 1s
    ##  73150K .......... .......... .......... .......... .......... 54% 95.1M 1s
    ##  73200K .......... .......... .......... .......... .......... 54% 45.6M 1s
    ##  73250K .......... .......... .......... .......... .......... 54% 72.5M 1s
    ##  73300K .......... .......... .......... .......... .......... 54% 76.3M 1s
    ##  73350K .......... .......... .......... .......... .......... 54% 82.3M 1s
    ##  73400K .......... .......... .......... .......... .......... 54% 84.5M 1s
    ##  73450K .......... .......... .......... .......... .......... 54%  105M 1s
    ##  73500K .......... .......... .......... .......... .......... 54% 81.0M 1s
    ##  73550K .......... .......... .......... .......... .......... 54% 71.5M 1s
    ##  73600K .......... .......... .......... .......... .......... 54% 88.1M 1s
    ##  73650K .......... .......... .......... .......... .......... 54% 69.4M 1s
    ##  73700K .......... .......... .......... .......... .......... 54% 93.5M 1s
    ##  73750K .......... .......... .......... .......... .......... 54% 97.0M 1s
    ##  73800K .......... .......... .......... .......... .......... 54% 85.5M 1s
    ##  73850K .......... .......... .......... .......... .......... 54% 79.8M 1s
    ##  73900K .......... .......... .......... .......... .......... 54% 78.3M 1s
    ##  73950K .......... .......... .......... .......... .......... 54% 34.5M 1s
    ##  74000K .......... .......... .......... .......... .......... 54% 69.9M 1s
    ##  74050K .......... .......... .......... .......... .......... 54% 93.2M 1s
    ##  74100K .......... .......... .......... .......... .......... 55% 66.0M 1s
    ##  74150K .......... .......... .......... .......... .......... 55% 85.9M 1s
    ##  74200K .......... .......... .......... .......... .......... 55% 84.5M 1s
    ##  74250K .......... .......... .......... .......... .......... 55%  106M 1s
    ##  74300K .......... .......... .......... .......... .......... 55% 23.8M 1s
    ##  74350K .......... .......... .......... .......... .......... 55% 81.3M 1s
    ##  74400K .......... .......... .......... .......... .......... 55% 87.4M 1s
    ##  74450K .......... .......... .......... .......... .......... 55% 91.0M 1s
    ##  74500K .......... .......... .......... .......... .......... 55% 90.8M 1s
    ##  74550K .......... .......... .......... .......... .......... 55%  114M 1s
    ##  74600K .......... .......... .......... .......... .......... 55% 49.7M 1s
    ##  74650K .......... .......... .......... .......... .......... 55% 91.3M 1s
    ##  74700K .......... .......... .......... .......... .......... 55% 64.3M 1s
    ##  74750K .......... .......... .......... .......... .......... 55% 29.3M 1s
    ##  74800K .......... .......... .......... .......... .......... 55% 37.8M 1s
    ##  74850K .......... .......... .......... .......... .......... 55% 94.6M 1s
    ##  74900K .......... .......... .......... .......... .......... 55% 45.2M 1s
    ##  74950K .......... .......... .......... .......... .......... 55% 94.1M 1s
    ##  75000K .......... .......... .......... .......... .......... 55% 93.0M 1s
    ##  75050K .......... .......... .......... .......... .......... 55% 72.4M 1s
    ##  75100K .......... .......... .......... .......... .......... 55% 70.6M 1s
    ##  75150K .......... .......... .......... .......... .......... 55% 74.9M 1s
    ##  75200K .......... .......... .......... .......... .......... 55% 68.5M 1s
    ##  75250K .......... .......... .......... .......... .......... 55% 74.7M 1s
    ##  75300K .......... .......... .......... .......... .......... 55% 77.2M 1s
    ##  75350K .......... .......... .......... .......... .......... 55% 79.1M 1s
    ##  75400K .......... .......... .......... .......... .......... 55% 89.5M 1s
    ##  75450K .......... .......... .......... .......... .......... 56% 72.8M 1s
    ##  75500K .......... .......... .......... .......... .......... 56% 69.3M 1s
    ##  75550K .......... .......... .......... .......... .......... 56% 91.8M 1s
    ##  75600K .......... .......... .......... .......... .......... 56% 74.6M 1s
    ##  75650K .......... .......... .......... .......... .......... 56% 92.1M 1s
    ##  75700K .......... .......... .......... .......... .......... 56% 91.1M 1s
    ##  75750K .......... .......... .......... .......... .......... 56%  100M 1s
    ##  75800K .......... .......... .......... .......... .......... 56% 96.2M 1s
    ##  75850K .......... .......... .......... .......... .......... 56% 74.7M 1s
    ##  75900K .......... .......... .......... .......... .......... 56% 83.9M 1s
    ##  75950K .......... .......... .......... .......... .......... 56% 92.7M 1s
    ##  76000K .......... .......... .......... .......... .......... 56% 89.9M 1s
    ##  76050K .......... .......... .......... .......... .......... 56%  102M 1s
    ##  76100K .......... .......... .......... .......... .......... 56% 99.7M 1s
    ##  76150K .......... .......... .......... .......... .......... 56% 98.1M 1s
    ##  76200K .......... .......... .......... .......... .......... 56% 81.7M 1s
    ##  76250K .......... .......... .......... .......... .......... 56%  101M 1s
    ##  76300K .......... .......... .......... .......... .......... 56%  107M 1s
    ##  76350K .......... .......... .......... .......... .......... 56% 92.0M 1s
    ##  76400K .......... .......... .......... .......... .......... 56% 89.6M 1s
    ##  76450K .......... .......... .......... .......... .......... 56% 89.6M 1s
    ##  76500K .......... .......... .......... .......... .......... 56%  101M 1s
    ##  76550K .......... .......... .......... .......... .......... 56%  111M 1s
    ##  76600K .......... .......... .......... .......... .......... 56% 17.1M 1s
    ##  76650K .......... .......... .......... .......... .......... 56% 88.2M 1s
    ##  76700K .......... .......... .......... .......... .......... 56% 48.5M 1s
    ##  76750K .......... .......... .......... .......... .......... 56% 67.9M 1s
    ##  76800K .......... .......... .......... .......... .......... 57% 75.6M 1s
    ##  76850K .......... .......... .......... .......... .......... 57% 92.4M 1s
    ##  76900K .......... .......... .......... .......... .......... 57%  100M 1s
    ##  76950K .......... .......... .......... .......... .......... 57%  102M 1s
    ##  77000K .......... .......... .......... .......... .......... 57%  101M 1s
    ##  77050K .......... .......... .......... .......... .......... 57% 76.8M 1s
    ##  77100K .......... .......... .......... .......... .......... 57% 88.2M 1s
    ##  77150K .......... .......... .......... .......... .......... 57% 42.5M 1s
    ##  77200K .......... .......... .......... .......... .......... 57% 83.8M 1s
    ##  77250K .......... .......... .......... .......... .......... 57% 72.5M 1s
    ##  77300K .......... .......... .......... .......... .......... 57% 99.8M 1s
    ##  77350K .......... .......... .......... .......... .......... 57%  106M 1s
    ##  77400K .......... .......... .......... .......... .......... 57%  102M 1s
    ##  77450K .......... .......... .......... .......... .......... 57%  115M 1s
    ##  77500K .......... .......... .......... .......... .......... 57% 23.9M 1s
    ##  77550K .......... .......... .......... .......... .......... 57% 83.9M 1s
    ##  77600K .......... .......... .......... .......... .......... 57%  106M 1s
    ##  77650K .......... .......... .......... .......... .......... 57%  115M 1s
    ##  77700K .......... .......... .......... .......... .......... 57% 37.7M 1s
    ##  77750K .......... .......... .......... .......... .......... 57%  125M 1s
    ##  77800K .......... .......... .......... .......... .......... 57% 97.3M 1s
    ##  77850K .......... .......... .......... .......... .......... 57%  105M 1s
    ##  77900K .......... .......... .......... .......... .......... 57% 83.3M 1s
    ##  77950K .......... .......... .......... .......... .......... 57%  112M 1s
    ##  78000K .......... .......... .......... .......... .......... 57% 28.5M 1s
    ##  78050K .......... .......... .......... .......... .......... 57% 85.0M 1s
    ##  78100K .......... .......... .......... .......... .......... 58% 90.7M 1s
    ##  78150K .......... .......... .......... .......... .......... 58% 82.3M 1s
    ##  78200K .......... .......... .......... .......... .......... 58% 94.8M 1s
    ##  78250K .......... .......... .......... .......... .......... 58%  117M 1s
    ##  78300K .......... .......... .......... .......... .......... 58%  115M 1s
    ##  78350K .......... .......... .......... .......... .......... 58% 60.2M 1s
    ##  78400K .......... .......... .......... .......... .......... 58% 39.6M 1s
    ##  78450K .......... .......... .......... .......... .......... 58%  108M 1s
    ##  78500K .......... .......... .......... .......... .......... 58%  106M 1s
    ##  78550K .......... .......... .......... .......... .......... 58%  129M 1s
    ##  78600K .......... .......... .......... .......... .......... 58%  105M 1s
    ##  78650K .......... .......... .......... .......... .......... 58%  118M 1s
    ##  78700K .......... .......... .......... .......... .......... 58% 68.7M 1s
    ##  78750K .......... .......... .......... .......... .......... 58% 48.1M 1s
    ##  78800K .......... .......... .......... .......... .......... 58% 47.9M 1s
    ##  78850K .......... .......... .......... .......... .......... 58% 61.0M 1s
    ##  78900K .......... .......... .......... .......... .......... 58%  102M 1s
    ##  78950K .......... .......... .......... .......... .......... 58%  142M 1s
    ##  79000K .......... .......... .......... .......... .......... 58% 85.7M 1s
    ##  79050K .......... .......... .......... .......... .......... 58%  134M 1s
    ##  79100K .......... .......... .......... .......... .......... 58%  105M 1s
    ##  79150K .......... .......... .......... .......... .......... 58% 57.1M 1s
    ##  79200K .......... .......... .......... .......... .......... 58% 63.4M 1s
    ##  79250K .......... .......... .......... .......... .......... 58% 50.0M 1s
    ##  79300K .......... .......... .......... .......... .......... 58% 79.0M 1s
    ##  79350K .......... .......... .......... .......... .......... 58% 87.8M 1s
    ##  79400K .......... .......... .......... .......... .......... 58%  108M 1s
    ##  79450K .......... .......... .......... .......... .......... 59%  106M 1s
    ##  79500K .......... .......... .......... .......... .......... 59%  104M 1s
    ##  79550K .......... .......... .......... .......... .......... 59%  104M 1s
    ##  79600K .......... .......... .......... .......... .......... 59% 85.9M 1s
    ##  79650K .......... .......... .......... .......... .......... 59% 90.8M 1s
    ##  79700K .......... .......... .......... .......... .......... 59% 90.8M 1s
    ##  79750K .......... .......... .......... .......... .......... 59% 99.1M 1s
    ##  79800K .......... .......... .......... .......... .......... 59% 93.8M 1s
    ##  79850K .......... .......... .......... .......... .......... 59% 35.7M 1s
    ##  79900K .......... .......... .......... .......... .......... 59% 82.0M 1s
    ##  79950K .......... .......... .......... .......... .......... 59%  135M 1s
    ##  80000K .......... .......... .......... .......... .......... 59% 95.4M 1s
    ##  80050K .......... .......... .......... .......... .......... 59% 69.6M 1s
    ##  80100K .......... .......... .......... .......... .......... 59% 71.8M 1s
    ##  80150K .......... .......... .......... .......... .......... 59% 92.9M 1s
    ##  80200K .......... .......... .......... .......... .......... 59%  112M 1s
    ##  80250K .......... .......... .......... .......... .......... 59% 32.9M 1s
    ##  80300K .......... .......... .......... .......... .......... 59%  125M 1s
    ##  80350K .......... .......... .......... .......... .......... 59% 97.2M 1s
    ##  80400K .......... .......... .......... .......... .......... 59% 98.6M 1s
    ##  80450K .......... .......... .......... .......... .......... 59%  147M 1s
    ##  80500K .......... .......... .......... .......... .......... 59% 88.4M 1s
    ##  80550K .......... .......... .......... .......... .......... 59%  111M 1s
    ##  80600K .......... .......... .......... .......... .......... 59% 71.9M 1s
    ##  80650K .......... .......... .......... .......... .......... 59% 63.8M 1s
    ##  80700K .......... .......... .......... .......... .......... 59% 38.8M 1s
    ##  80750K .......... .......... .......... .......... .......... 59%  139M 1s
    ##  80800K .......... .......... .......... .......... .......... 60% 94.5M 1s
    ##  80850K .......... .......... .......... .......... .......... 60%  112M 1s
    ##  80900K .......... .......... .......... .......... .......... 60%  138M 1s
    ##  80950K .......... .......... .......... .......... .......... 60% 21.1M 1s
    ##  81000K .......... .......... .......... .......... .......... 60% 97.8M 1s
    ##  81050K .......... .......... .......... .......... .......... 60%  135M 1s
    ##  81100K .......... .......... .......... .......... .......... 60%  113M 1s
    ##  81150K .......... .......... .......... .......... .......... 60%  123M 1s
    ##  81200K .......... .......... .......... .......... .......... 60%  126M 1s
    ##  81250K .......... .......... .......... .......... .......... 60% 27.3M 1s
    ##  81300K .......... .......... .......... .......... .......... 60% 98.6M 1s
    ##  81350K .......... .......... .......... .......... .......... 60% 45.1M 1s
    ##  81400K .......... .......... .......... .......... .......... 60% 96.5M 1s
    ##  81450K .......... .......... .......... .......... .......... 60% 95.8M 1s
    ##  81500K .......... .......... .......... .......... .......... 60%  108M 1s
    ##  81550K .......... .......... .......... .......... .......... 60%  107M 1s
    ##  81600K .......... .......... .......... .......... .......... 60%  101M 1s
    ##  81650K .......... .......... .......... .......... .......... 60%  117M 1s
    ##  81700K .......... .......... .......... .......... .......... 60% 77.1M 1s
    ##  81750K .......... .......... .......... .......... .......... 60% 66.2M 1s
    ##  81800K .......... .......... .......... .......... .......... 60%  121M 1s
    ##  81850K .......... .......... .......... .......... .......... 60% 57.5M 1s
    ##  81900K .......... .......... .......... .......... .......... 60% 96.5M 1s
    ##  81950K .......... .......... .......... .......... .......... 60% 93.1M 1s
    ##  82000K .......... .......... .......... .......... .......... 60% 60.2M 1s
    ##  82050K .......... .......... .......... .......... .......... 60%  131M 1s
    ##  82100K .......... .......... .......... .......... .......... 60% 98.3M 1s
    ##  82150K .......... .......... .......... .......... .......... 61% 69.4M 1s
    ##  82200K .......... .......... .......... .......... .......... 61% 82.4M 1s
    ##  82250K .......... .......... .......... .......... .......... 61% 76.3M 1s
    ##  82300K .......... .......... .......... .......... .......... 61%  118M 1s
    ##  82350K .......... .......... .......... .......... .......... 61% 63.3M 1s
    ##  82400K .......... .......... .......... .......... .......... 61% 95.6M 1s
    ##  82450K .......... .......... .......... .......... .......... 61% 69.1M 1s
    ##  82500K .......... .......... .......... .......... .......... 61%  101M 1s
    ##  82550K .......... .......... .......... .......... .......... 61% 44.0M 1s
    ##  82600K .......... .......... .......... .......... .......... 61% 49.0M 1s
    ##  82650K .......... .......... .......... .......... .......... 61% 58.4M 1s
    ##  82700K .......... .......... .......... .......... .......... 61% 50.0M 1s
    ##  82750K .......... .......... .......... .......... .......... 61% 60.4M 1s
    ##  82800K .......... .......... .......... .......... .......... 61% 72.0M 1s
    ##  82850K .......... .......... .......... .......... .......... 61% 86.6M 1s
    ##  82900K .......... .......... .......... .......... .......... 61%  112M 1s
    ##  82950K .......... .......... .......... .......... .......... 61%  113M 1s
    ##  83000K .......... .......... .......... .......... .......... 61% 32.0M 1s
    ##  83050K .......... .......... .......... .......... .......... 61% 90.4M 1s
    ##  83100K .......... .......... .......... .......... .......... 61% 62.5M 1s
    ##  83150K .......... .......... .......... .......... .......... 61% 79.8M 1s
    ##  83200K .......... .......... .......... .......... .......... 61% 89.0M 1s
    ##  83250K .......... .......... .......... .......... .......... 61%  120M 1s
    ##  83300K .......... .......... .......... .......... .......... 61%  105M 1s
    ##  83350K .......... .......... .......... .......... .......... 61% 44.6M 1s
    ##  83400K .......... .......... .......... .......... .......... 61% 83.9M 1s
    ##  83450K .......... .......... .......... .......... .......... 61% 26.3M 1s
    ##  83500K .......... .......... .......... .......... .......... 62%  106M 1s
    ##  83550K .......... .......... .......... .......... .......... 62% 53.1M 1s
    ##  83600K .......... .......... .......... .......... .......... 62% 67.4M 1s
    ##  83650K .......... .......... .......... .......... .......... 62% 64.9M 1s
    ##  83700K .......... .......... .......... .......... .......... 62% 78.0M 1s
    ##  83750K .......... .......... .......... .......... .......... 62% 56.4M 1s
    ##  83800K .......... .......... .......... .......... .......... 62% 53.7M 1s
    ##  83850K .......... .......... .......... .......... .......... 62% 58.9M 1s
    ##  83900K .......... .......... .......... .......... .......... 62% 66.0M 1s
    ##  83950K .......... .......... .......... .......... .......... 62%  133M 1s
    ##  84000K .......... .......... .......... .......... .......... 62%  102M 1s
    ##  84050K .......... .......... .......... .......... .......... 62% 69.4M 1s
    ##  84100K .......... .......... .......... .......... .......... 62% 79.6M 1s
    ##  84150K .......... .......... .......... .......... .......... 62% 42.5M 1s
    ##  84200K .......... .......... .......... .......... .......... 62% 51.7M 1s
    ##  84250K .......... .......... .......... .......... .......... 62%  123M 1s
    ##  84300K .......... .......... .......... .......... .......... 62% 47.0M 1s
    ##  84350K .......... .......... .......... .......... .......... 62% 70.8M 1s
    ##  84400K .......... .......... .......... .......... .......... 62% 93.3M 1s
    ##  84450K .......... .......... .......... .......... .......... 62%  116M 1s
    ##  84500K .......... .......... .......... .......... .......... 62% 43.2M 1s
    ##  84550K .......... .......... .......... .......... .......... 62% 72.1M 1s
    ##  84600K .......... .......... .......... .......... .......... 62% 67.1M 1s
    ##  84650K .......... .......... .......... .......... .......... 62% 85.2M 1s
    ##  84700K .......... .......... .......... .......... .......... 62% 88.0M 1s
    ##  84750K .......... .......... .......... .......... .......... 62% 66.4M 1s
    ##  84800K .......... .......... .......... .......... .......... 62% 68.5M 1s
    ##  84850K .......... .......... .......... .......... .......... 63% 62.0M 1s
    ##  84900K .......... .......... .......... .......... .......... 63% 44.0M 1s
    ##  84950K .......... .......... .......... .......... .......... 63% 76.2M 1s
    ##  85000K .......... .......... .......... .......... .......... 63% 40.3M 1s
    ##  85050K .......... .......... .......... .......... .......... 63%  128M 1s
    ##  85100K .......... .......... .......... .......... .......... 63% 72.0M 1s
    ##  85150K .......... .......... .......... .......... .......... 63% 55.4M 1s
    ##  85200K .......... .......... .......... .......... .......... 63% 60.9M 1s
    ##  85250K .......... .......... .......... .......... .......... 63% 43.2M 1s
    ##  85300K .......... .......... .......... .......... .......... 63%  125M 1s
    ##  85350K .......... .......... .......... .......... .......... 63%  117M 1s
    ##  85400K .......... .......... .......... .......... .......... 63%  116M 1s
    ##  85450K .......... .......... .......... .......... .......... 63% 72.0M 1s
    ##  85500K .......... .......... .......... .......... .......... 63%  136M 1s
    ##  85550K .......... .......... .......... .......... .......... 63%  104M 1s
    ##  85600K .......... .......... .......... .......... .......... 63% 54.9M 1s
    ##  85650K .......... .......... .......... .......... .......... 63%  124M 1s
    ##  85700K .......... .......... .......... .......... .......... 63% 75.8M 1s
    ##  85750K .......... .......... .......... .......... .......... 63% 91.1M 1s
    ##  85800K .......... .......... .......... .......... .......... 63%  102M 1s
    ##  85850K .......... .......... .......... .......... .......... 63% 72.7M 1s
    ##  85900K .......... .......... .......... .......... .......... 63% 90.2M 1s
    ##  85950K .......... .......... .......... .......... .......... 63% 64.5M 1s
    ##  86000K .......... .......... .......... .......... .......... 63% 34.9M 1s
    ##  86050K .......... .......... .......... .......... .......... 63%  128M 1s
    ##  86100K .......... .......... .......... .......... .......... 63% 93.0M 1s
    ##  86150K .......... .......... .......... .......... .......... 63%  104M 1s
    ##  86200K .......... .......... .......... .......... .......... 64%  121M 1s
    ##  86250K .......... .......... .......... .......... .......... 64% 98.4M 1s
    ##  86300K .......... .......... .......... .......... .......... 64%  128M 1s
    ##  86350K .......... .......... .......... .......... .......... 64% 40.6M 1s
    ##  86400K .......... .......... .......... .......... .......... 64% 48.7M 1s
    ##  86450K .......... .......... .......... .......... .......... 64% 90.3M 1s
    ##  86500K .......... .......... .......... .......... .......... 64%  110M 1s
    ##  86550K .......... .......... .......... .......... .......... 64%  108M 1s
    ##  86600K .......... .......... .......... .......... .......... 64%  121M 1s
    ##  86650K .......... .......... .......... .......... .......... 64% 81.7M 1s
    ##  86700K .......... .......... .......... .......... .......... 64% 78.2M 1s
    ##  86750K .......... .......... .......... .......... .......... 64% 89.8M 1s
    ##  86800K .......... .......... .......... .......... .......... 64% 49.0M 1s
    ##  86850K .......... .......... .......... .......... .......... 64% 66.9M 1s
    ##  86900K .......... .......... .......... .......... .......... 64% 83.9M 1s
    ##  86950K .......... .......... .......... .......... .......... 64%  134M 1s
    ##  87000K .......... .......... .......... .......... .......... 64% 64.5M 1s
    ##  87050K .......... .......... .......... .......... .......... 64% 86.0M 1s
    ##  87100K .......... .......... .......... .......... .......... 64% 77.4M 1s
    ##  87150K .......... .......... .......... .......... .......... 64% 76.8M 1s
    ##  87200K .......... .......... .......... .......... .......... 64%  125M 1s
    ##  87250K .......... .......... .......... .......... .......... 64% 61.7M 1s
    ##  87300K .......... .......... .......... .......... .......... 64%  134M 1s
    ##  87350K .......... .......... .......... .......... .......... 64% 95.8M 1s
    ##  87400K .......... .......... .......... .......... .......... 64%  132M 1s
    ##  87450K .......... .......... .......... .......... .......... 64%  118M 1s
    ##  87500K .......... .......... .......... .......... .......... 64%  129M 1s
    ##  87550K .......... .......... .......... .......... .......... 65%  137M 1s
    ##  87600K .......... .......... .......... .......... .......... 65% 83.5M 1s
    ##  87650K .......... .......... .......... .......... .......... 65%  152M 1s
    ##  87700K .......... .......... .......... .......... .......... 65%  133M 1s
    ##  87750K .......... .......... .......... .......... .......... 65%  139M 1s
    ##  87800K .......... .......... .......... .......... .......... 65%  156M 1s
    ##  87850K .......... .......... .......... .......... .......... 65% 35.7M 1s
    ##  87900K .......... .......... .......... .......... .......... 65%  121M 1s
    ##  87950K .......... .......... .......... .......... .......... 65%  107M 1s
    ##  88000K .......... .......... .......... .......... .......... 65%  135M 1s
    ##  88050K .......... .......... .......... .......... .......... 65%  126M 1s
    ##  88100K .......... .......... .......... .......... .......... 65% 41.6M 1s
    ##  88150K .......... .......... .......... .......... .......... 65% 25.8M 1s
    ##  88200K .......... .......... .......... .......... .......... 65% 49.3M 1s
    ##  88250K .......... .......... .......... .......... .......... 65% 86.5M 1s
    ##  88300K .......... .......... .......... .......... .......... 65%  131M 1s
    ##  88350K .......... .......... .......... .......... .......... 65%  105M 1s
    ##  88400K .......... .......... .......... .......... .......... 65%  125M 1s
    ##  88450K .......... .......... .......... .......... .......... 65%  149M 1s
    ##  88500K .......... .......... .......... .......... .......... 65% 20.9M 1s
    ##  88550K .......... .......... .......... .......... .......... 65% 88.8M 1s
    ##  88600K .......... .......... .......... .......... .......... 65%  113M 1s
    ##  88650K .......... .......... .......... .......... .......... 65% 86.5M 1s
    ##  88700K .......... .......... .......... .......... .......... 65%  157M 1s
    ##  88750K .......... .......... .......... .......... .......... 65%  154M 1s
    ##  88800K .......... .......... .......... .......... .......... 65%  114M 1s
    ##  88850K .......... .......... .......... .......... .......... 65% 37.9M 1s
    ##  88900K .......... .......... .......... .......... .......... 66% 99.7M 1s
    ##  88950K .......... .......... .......... .......... .......... 66% 50.6M 1s
    ##  89000K .......... .......... .......... .......... .......... 66% 52.8M 1s
    ##  89050K .......... .......... .......... .......... .......... 66% 86.7M 1s
    ##  89100K .......... .......... .......... .......... .......... 66%  118M 1s
    ##  89150K .......... .......... .......... .......... .......... 66%  101M 1s
    ##  89200K .......... .......... .......... .......... .......... 66%  100M 1s
    ##  89250K .......... .......... .......... .......... .......... 66%  124M 1s
    ##  89300K .......... .......... .......... .......... .......... 66% 99.2M 1s
    ##  89350K .......... .......... .......... .......... .......... 66% 93.2M 1s
    ##  89400K .......... .......... .......... .......... .......... 66% 54.7M 1s
    ##  89450K .......... .......... .......... .......... .......... 66% 29.7M 1s
    ##  89500K .......... .......... .......... .......... .......... 66%  132M 1s
    ##  89550K .......... .......... .......... .......... .......... 66%  143M 1s
    ##  89600K .......... .......... .......... .......... .......... 66%  169M 1s
    ##  89650K .......... .......... .......... .......... .......... 66% 18.5M 1s
    ##  89700K .......... .......... .......... .......... .......... 66%  141M 1s
    ##  89750K .......... .......... .......... .......... .......... 66%  167M 1s
    ##  89800K .......... .......... .......... .......... .......... 66%  103M 1s
    ##  89850K .......... .......... .......... .......... .......... 66% 63.8M 1s
    ##  89900K .......... .......... .......... .......... .......... 66%  115M 1s
    ##  89950K .......... .......... .......... .......... .......... 66% 43.7M 1s
    ##  90000K .......... .......... .......... .......... .......... 66% 44.6M 1s
    ##  90050K .......... .......... .......... .......... .......... 66% 46.0M 1s
    ##  90100K .......... .......... .......... .......... .......... 66% 64.9M 1s
    ##  90150K .......... .......... .......... .......... .......... 66%  145M 1s
    ##  90200K .......... .......... .......... .......... .......... 66%  134M 1s
    ##  90250K .......... .......... .......... .......... .......... 67%  114M 1s
    ##  90300K .......... .......... .......... .......... .......... 67% 94.3M 1s
    ##  90350K .......... .......... .......... .......... .......... 67% 80.2M 1s
    ##  90400K .......... .......... .......... .......... .......... 67%  135M 1s
    ##  90450K .......... .......... .......... .......... .......... 67% 28.5M 1s
    ##  90500K .......... .......... .......... .......... .......... 67% 29.9M 1s
    ##  90550K .......... .......... .......... .......... .......... 67% 36.2M 1s
    ##  90600K .......... .......... .......... .......... .......... 67% 49.7M 1s
    ##  90650K .......... .......... .......... .......... .......... 67% 50.8M 1s
    ##  90700K .......... .......... .......... .......... .......... 67% 64.2M 1s
    ##  90750K .......... .......... .......... .......... .......... 67% 23.4M 1s
    ##  90800K .......... .......... .......... .......... .......... 67% 97.7M 1s
    ##  90850K .......... .......... .......... .......... .......... 67% 78.6M 1s
    ##  90900K .......... .......... .......... .......... .......... 67% 81.6M 1s
    ##  90950K .......... .......... .......... .......... .......... 67% 77.7M 1s
    ##  91000K .......... .......... .......... .......... .......... 67% 67.7M 1s
    ##  91050K .......... .......... .......... .......... .......... 67% 59.8M 1s
    ##  91100K .......... .......... .......... .......... .......... 67% 43.3M 1s
    ##  91150K .......... .......... .......... .......... .......... 67% 69.8M 1s
    ##  91200K .......... .......... .......... .......... .......... 67% 36.6M 1s
    ##  91250K .......... .......... .......... .......... .......... 67% 58.6M 1s
    ##  91300K .......... .......... .......... .......... .......... 67% 50.5M 1s
    ##  91350K .......... .......... .......... .......... .......... 67% 59.2M 1s
    ##  91400K .......... .......... .......... .......... .......... 67% 72.3M 1s
    ##  91450K .......... .......... .......... .......... .......... 67% 47.2M 1s
    ##  91500K .......... .......... .......... .......... .......... 67% 85.3M 1s
    ##  91550K .......... .......... .......... .......... .......... 67% 25.8M 1s
    ##  91600K .......... .......... .......... .......... .......... 68% 47.2M 1s
    ##  91650K .......... .......... .......... .......... .......... 68% 87.2M 1s
    ##  91700K .......... .......... .......... .......... .......... 68% 70.7M 1s
    ##  91750K .......... .......... .......... .......... .......... 68% 73.0M 1s
    ##  91800K .......... .......... .......... .......... .......... 68% 92.1M 1s
    ##  91850K .......... .......... .......... .......... .......... 68%  122M 1s
    ##  91900K .......... .......... .......... .......... .......... 68% 30.4M 1s
    ##  91950K .......... .......... .......... .......... .......... 68%  109M 1s
    ##  92000K .......... .......... .......... .......... .......... 68% 59.8M 1s
    ##  92050K .......... .......... .......... .......... .......... 68% 82.0M 1s
    ##  92100K .......... .......... .......... .......... .......... 68% 62.9M 1s
    ##  92150K .......... .......... .......... .......... .......... 68%  107M 1s
    ##  92200K .......... .......... .......... .......... .......... 68% 76.4M 1s
    ##  92250K .......... .......... .......... .......... .......... 68% 66.1M 1s
    ##  92300K .......... .......... .......... .......... .......... 68% 59.9M 1s
    ##  92350K .......... .......... .......... .......... .......... 68% 36.3M 1s
    ##  92400K .......... .......... .......... .......... .......... 68% 95.4M 1s
    ##  92450K .......... .......... .......... .......... .......... 68% 99.8M 1s
    ##  92500K .......... .......... .......... .......... .......... 68% 59.0M 1s
    ##  92550K .......... .......... .......... .......... .......... 68%  120M 1s
    ##  92600K .......... .......... .......... .......... .......... 68% 48.9M 1s
    ##  92650K .......... .......... .......... .......... .......... 68% 67.9M 1s
    ##  92700K .......... .......... .......... .......... .......... 68% 45.4M 1s
    ##  92750K .......... .......... .......... .......... .......... 68% 98.7M 1s
    ##  92800K .......... .......... .......... .......... .......... 68%  118M 1s
    ##  92850K .......... .......... .......... .......... .......... 68% 37.4M 1s
    ##  92900K .......... .......... .......... .......... .......... 68% 78.6M 1s
    ##  92950K .......... .......... .......... .......... .......... 69%  127M 1s
    ##  93000K .......... .......... .......... .......... .......... 69% 60.6M 1s
    ##  93050K .......... .......... .......... .......... .......... 69% 84.3M 1s
    ##  93100K .......... .......... .......... .......... .......... 69% 91.1M 1s
    ##  93150K .......... .......... .......... .......... .......... 69%  130M 1s
    ##  93200K .......... .......... .......... .......... .......... 69% 86.5M 1s
    ##  93250K .......... .......... .......... .......... .......... 69% 98.3M 1s
    ##  93300K .......... .......... .......... .......... .......... 69% 71.3M 1s
    ##  93350K .......... .......... .......... .......... .......... 69% 66.2M 1s
    ##  93400K .......... .......... .......... .......... .......... 69% 74.7M 1s
    ##  93450K .......... .......... .......... .......... .......... 69%  135M 1s
    ##  93500K .......... .......... .......... .......... .......... 69%  114M 1s
    ##  93550K .......... .......... .......... .......... .......... 69% 43.6M 1s
    ##  93600K .......... .......... .......... .......... .......... 69% 67.2M 1s
    ##  93650K .......... .......... .......... .......... .......... 69% 90.7M 1s
    ##  93700K .......... .......... .......... .......... .......... 69% 78.8M 1s
    ##  93750K .......... .......... .......... .......... .......... 69% 99.5M 1s
    ##  93800K .......... .......... .......... .......... .......... 69% 38.4M 1s
    ##  93850K .......... .......... .......... .......... .......... 69% 54.8M 1s
    ##  93900K .......... .......... .......... .......... .......... 69% 81.1M 1s
    ##  93950K .......... .......... .......... .......... .......... 69% 92.2M 1s
    ##  94000K .......... .......... .......... .......... .......... 69% 59.2M 1s
    ##  94050K .......... .......... .......... .......... .......... 69% 51.9M 1s
    ##  94100K .......... .......... .......... .......... .......... 69% 90.4M 1s
    ##  94150K .......... .......... .......... .......... .......... 69%  103M 1s
    ##  94200K .......... .......... .......... .......... .......... 69% 94.2M 1s
    ##  94250K .......... .......... .......... .......... .......... 69%  100M 1s
    ##  94300K .......... .......... .......... .......... .......... 70% 49.5M 1s
    ##  94350K .......... .......... .......... .......... .......... 70%  107M 1s
    ##  94400K .......... .......... .......... .......... .......... 70% 76.3M 1s
    ##  94450K .......... .......... .......... .......... .......... 70%  119M 1s
    ##  94500K .......... .......... .......... .......... .......... 70%  106M 1s
    ##  94550K .......... .......... .......... .......... .......... 70% 94.0M 1s
    ##  94600K .......... .......... .......... .......... .......... 70% 58.0M 1s
    ##  94650K .......... .......... .......... .......... .......... 70% 67.0M 1s
    ##  94700K .......... .......... .......... .......... .......... 70% 76.6M 1s
    ##  94750K .......... .......... .......... .......... .......... 70%  122M 1s
    ##  94800K .......... .......... .......... .......... .......... 70%  105M 1s
    ##  94850K .......... .......... .......... .......... .......... 70% 94.7M 1s
    ##  94900K .......... .......... .......... .......... .......... 70%  121M 1s
    ##  94950K .......... .......... .......... .......... .......... 70% 62.4M 1s
    ##  95000K .......... .......... .......... .......... .......... 70% 73.2M 1s
    ##  95050K .......... .......... .......... .......... .......... 70%  139M 1s
    ##  95100K .......... .......... .......... .......... .......... 70%  117M 1s
    ##  95150K .......... .......... .......... .......... .......... 70%  118M 1s
    ##  95200K .......... .......... .......... .......... .......... 70% 95.7M 1s
    ##  95250K .......... .......... .......... .......... .......... 70% 80.9M 1s
    ##  95300K .......... .......... .......... .......... .......... 70% 73.9M 1s
    ##  95350K .......... .......... .......... .......... .......... 70% 47.6M 1s
    ##  95400K .......... .......... .......... .......... .......... 70%  115M 1s
    ##  95450K .......... .......... .......... .......... .......... 70%  101M 1s
    ##  95500K .......... .......... .......... .......... .......... 70%  108M 1s
    ##  95550K .......... .......... .......... .......... .......... 70%  133M 1s
    ##  95600K .......... .......... .......... .......... .......... 70% 73.1M 1s
    ##  95650K .......... .......... .......... .......... .......... 71%  126M 1s
    ##  95700K .......... .......... .......... .......... .......... 71% 74.0M 1s
    ##  95750K .......... .......... .......... .......... .......... 71% 74.6M 1s
    ##  95800K .......... .......... .......... .......... .......... 71% 53.1M 1s
    ##  95850K .......... .......... .......... .......... .......... 71% 78.2M 1s
    ##  95900K .......... .......... .......... .......... .......... 71%  123M 1s
    ##  95950K .......... .......... .......... .......... .......... 71%  144M 1s
    ##  96000K .......... .......... .......... .......... .......... 71%  102M 1s
    ##  96050K .......... .......... .......... .......... .......... 71% 65.9M 1s
    ##  96100K .......... .......... .......... .......... .......... 71% 58.5M 1s
    ##  96150K .......... .......... .......... .......... .......... 71% 95.2M 1s
    ##  96200K .......... .......... .......... .......... .......... 71% 79.4M 1s
    ##  96250K .......... .......... .......... .......... .......... 71% 89.2M 1s
    ##  96300K .......... .......... .......... .......... .......... 71% 58.0M 1s
    ##  96350K .......... .......... .......... .......... .......... 71% 83.3M 1s
    ##  96400K .......... .......... .......... .......... .......... 71%  126M 1s
    ##  96450K .......... .......... .......... .......... .......... 71% 59.5M 1s
    ##  96500K .......... .......... .......... .......... .......... 71% 77.9M 1s
    ##  96550K .......... .......... .......... .......... .......... 71%  118M 1s
    ##  96600K .......... .......... .......... .......... .......... 71% 93.5M 1s
    ##  96650K .......... .......... .......... .......... .......... 71% 82.8M 1s
    ##  96700K .......... .......... .......... .......... .......... 71% 98.2M 1s
    ##  96750K .......... .......... .......... .......... .......... 71% 57.2M 1s
    ##  96800K .......... .......... .......... .......... .......... 71%  126M 1s
    ##  96850K .......... .......... .......... .......... .......... 71% 63.2M 1s
    ##  96900K .......... .......... .......... .......... .......... 71% 89.4M 1s
    ##  96950K .......... .......... .......... .......... .......... 71%  101M 1s
    ##  97000K .......... .......... .......... .......... .......... 72% 69.4M 1s
    ##  97050K .......... .......... .......... .......... .......... 72%  177M 1s
    ##  97100K .......... .......... .......... .......... .......... 72% 77.8M 1s
    ##  97150K .......... .......... .......... .......... .......... 72% 77.7M 1s
    ##  97200K .......... .......... .......... .......... .......... 72%  111M 1s
    ##  97250K .......... .......... .......... .......... .......... 72%  139M 1s
    ##  97300K .......... .......... .......... .......... .......... 72% 92.9M 1s
    ##  97350K .......... .......... .......... .......... .......... 72%  124M 1s
    ##  97400K .......... .......... .......... .......... .......... 72% 61.6M 1s
    ##  97450K .......... .......... .......... .......... .......... 72% 92.7M 1s
    ##  97500K .......... .......... .......... .......... .......... 72% 96.1M 1s
    ##  97550K .......... .......... .......... .......... .......... 72%  111M 1s
    ##  97600K .......... .......... .......... .......... .......... 72% 54.2M 1s
    ##  97650K .......... .......... .......... .......... .......... 72%  140M 1s
    ##  97700K .......... .......... .......... .......... .......... 72%  101M 1s
    ##  97750K .......... .......... .......... .......... .......... 72%  137M 1s
    ##  97800K .......... .......... .......... .......... .......... 72%  139M 1s
    ##  97850K .......... .......... .......... .......... .......... 72% 73.2M 1s
    ##  97900K .......... .......... .......... .......... .......... 72% 94.0M 1s
    ##  97950K .......... .......... .......... .......... .......... 72% 62.6M 1s
    ##  98000K .......... .......... .......... .......... .......... 72% 93.8M 1s
    ##  98050K .......... .......... .......... .......... .......... 72% 72.9M 1s
    ##  98100K .......... .......... .......... .......... .......... 72% 77.3M 1s
    ##  98150K .......... .......... .......... .......... .......... 72%  137M 1s
    ##  98200K .......... .......... .......... .......... .......... 72% 74.9M 1s
    ##  98250K .......... .......... .......... .......... .......... 72% 70.8M 1s
    ##  98300K .......... .......... .......... .......... .......... 72%  126M 1s
    ##  98350K .......... .......... .......... .......... .......... 73% 72.4M 1s
    ##  98400K .......... .......... .......... .......... .......... 73%  147M 1s
    ##  98450K .......... .......... .......... .......... .......... 73%  180M 1s
    ##  98500K .......... .......... .......... .......... .......... 73% 95.6M 1s
    ##  98550K .......... .......... .......... .......... .......... 73% 96.7M 1s
    ##  98600K .......... .......... .......... .......... .......... 73% 99.5M 1s
    ##  98650K .......... .......... .......... .......... .......... 73% 68.7M 1s
    ##  98700K .......... .......... .......... .......... .......... 73%  148M 1s
    ##  98750K .......... .......... .......... .......... .......... 73% 69.8M 1s
    ##  98800K .......... .......... .......... .......... .......... 73% 97.6M 1s
    ##  98850K .......... .......... .......... .......... .......... 73%  149M 1s
    ##  98900K .......... .......... .......... .......... .......... 73%  125M 1s
    ##  98950K .......... .......... .......... .......... .......... 73% 69.7M 1s
    ##  99000K .......... .......... .......... .......... .......... 73% 78.7M 1s
    ##  99050K .......... .......... .......... .......... .......... 73%  148M 1s
    ##  99100K .......... .......... .......... .......... .......... 73% 82.8M 1s
    ##  99150K .......... .......... .......... .......... .......... 73%  126M 1s
    ##  99200K .......... .......... .......... .......... .......... 73%  161M 1s
    ##  99250K .......... .......... .......... .......... .......... 73%  111M 1s
    ##  99300K .......... .......... .......... .......... .......... 73% 47.2M 1s
    ##  99350K .......... .......... .......... .......... .......... 73%  139M 1s
    ##  99400K .......... .......... .......... .......... .......... 73% 49.4M 1s
    ##  99450K .......... .......... .......... .......... .......... 73%  160M 1s
    ##  99500K .......... .......... .......... .......... .......... 73% 86.3M 1s
    ##  99550K .......... .......... .......... .......... .......... 73%  119M 1s
    ##  99600K .......... .......... .......... .......... .......... 73%  142M 1s
    ##  99650K .......... .......... .......... .......... .......... 73% 88.4M 1s
    ##  99700K .......... .......... .......... .......... .......... 74% 75.1M 1s
    ##  99750K .......... .......... .......... .......... .......... 74%  146M 1s
    ##  99800K .......... .......... .......... .......... .......... 74%  150M 1s
    ##  99850K .......... .......... .......... .......... .......... 74% 78.3M 1s
    ##  99900K .......... .......... .......... .......... .......... 74%  164M 1s
    ##  99950K .......... .......... .......... .......... .......... 74% 75.1M 1s
    ## 100000K .......... .......... .......... .......... .......... 74% 52.9M 1s
    ## 100050K .......... .......... .......... .......... .......... 74% 59.4M 1s
    ## 100100K .......... .......... .......... .......... .......... 74%  109M 1s
    ## 100150K .......... .......... .......... .......... .......... 74% 88.5M 1s
    ## 100200K .......... .......... .......... .......... .......... 74%  131M 1s
    ## 100250K .......... .......... .......... .......... .......... 74%  121M 1s
    ## 100300K .......... .......... .......... .......... .......... 74%  122M 1s
    ## 100350K .......... .......... .......... .......... .......... 74%  119M 1s
    ## 100400K .......... .......... .......... .......... .......... 74% 89.3M 1s
    ## 100450K .......... .......... .......... .......... .......... 74%  123M 1s
    ## 100500K .......... .......... .......... .......... .......... 74% 39.5M 1s
    ## 100550K .......... .......... .......... .......... .......... 74%  118M 1s
    ## 100600K .......... .......... .......... .......... .......... 74%  125M 1s
    ## 100650K .......... .......... .......... .......... .......... 74% 75.2M 1s
    ## 100700K .......... .......... .......... .......... .......... 74% 88.4M 1s
    ## 100750K .......... .......... .......... .......... .......... 74%  131M 1s
    ## 100800K .......... .......... .......... .......... .......... 74% 72.8M 1s
    ## 100850K .......... .......... .......... .......... .......... 74% 93.4M 1s
    ## 100900K .......... .......... .......... .......... .......... 74% 96.3M 1s
    ## 100950K .......... .......... .......... .......... .......... 74%  114M 1s
    ## 101000K .......... .......... .......... .......... .......... 74% 99.1M 1s
    ## 101050K .......... .......... .......... .......... .......... 75%  147M 1s
    ## 101100K .......... .......... .......... .......... .......... 75%  137M 1s
    ## 101150K .......... .......... .......... .......... .......... 75%  147M 1s
    ## 101200K .......... .......... .......... .......... .......... 75%  154M 1s
    ## 101250K .......... .......... .......... .......... .......... 75%  121M 1s
    ## 101300K .......... .......... .......... .......... .......... 75%  131M 1s
    ## 101350K .......... .......... .......... .......... .......... 75%  150M 1s
    ## 101400K .......... .......... .......... .......... .......... 75% 67.1M 1s
    ## 101450K .......... .......... .......... .......... .......... 75%  122M 1s
    ## 101500K .......... .......... .......... .......... .......... 75% 93.2M 1s
    ## 101550K .......... .......... .......... .......... .......... 75%  157M 1s
    ## 101600K .......... .......... .......... .......... .......... 75%  140M 1s
    ## 101650K .......... .......... .......... .......... .......... 75%  162M 1s
    ## 101700K .......... .......... .......... .......... .......... 75%  109M 1s
    ## 101750K .......... .......... .......... .......... .......... 75%  171M 1s
    ## 101800K .......... .......... .......... .......... .......... 75%  131M 1s
    ## 101850K .......... .......... .......... .......... .......... 75%  188M 1s
    ## 101900K .......... .......... .......... .......... .......... 75% 45.0M 1s
    ## 101950K .......... .......... .......... .......... .......... 75% 27.7M 1s
    ## 102000K .......... .......... .......... .......... .......... 75% 49.1M 1s
    ## 102050K .......... .......... .......... .......... .......... 75% 64.7M 1s
    ## 102100K .......... .......... .......... .......... .......... 75% 28.7M 1s
    ## 102150K .......... .......... .......... .......... .......... 75%  100M 1s
    ## 102200K .......... .......... .......... .......... .......... 75% 59.1M 1s
    ## 102250K .......... .......... .......... .......... .......... 75% 79.5M 1s
    ## 102300K .......... .......... .......... .......... .......... 75% 36.6M 1s
    ## 102350K .......... .......... .......... .......... .......... 75% 43.6M 1s
    ## 102400K .......... .......... .......... .......... .......... 76% 63.9M 1s
    ## 102450K .......... .......... .......... .......... .......... 76% 60.1M 1s
    ## 102500K .......... .......... .......... .......... .......... 76% 33.2M 1s
    ## 102550K .......... .......... .......... .......... .......... 76% 62.6M 1s
    ## 102600K .......... .......... .......... .......... .......... 76% 52.2M 1s
    ## 102650K .......... .......... .......... .......... .......... 76% 39.9M 1s
    ## 102700K .......... .......... .......... .......... .......... 76% 53.4M 1s
    ## 102750K .......... .......... .......... .......... .......... 76% 49.2M 1s
    ## 102800K .......... .......... .......... .......... .......... 76%  140M 1s
    ## 102850K .......... .......... .......... .......... .......... 76% 91.1M 1s
    ## 102900K .......... .......... .......... .......... .......... 76% 63.0M 1s
    ## 102950K .......... .......... .......... .......... .......... 76%  143M 1s
    ## 103000K .......... .......... .......... .......... .......... 76%  165M 1s
    ## 103050K .......... .......... .......... .......... .......... 76%  105M 1s
    ## 103100K .......... .......... .......... .......... .......... 76%  143M 1s
    ## 103150K .......... .......... .......... .......... .......... 76%  101M 1s
    ## 103200K .......... .......... .......... .......... .......... 76% 74.1M 1s
    ## 103250K .......... .......... .......... .......... .......... 76%  158M 1s
    ## 103300K .......... .......... .......... .......... .......... 76% 91.0M 1s
    ## 103350K .......... .......... .......... .......... .......... 76% 53.8M 1s
    ## 103400K .......... .......... .......... .......... .......... 76% 98.4M 1s
    ## 103450K .......... .......... .......... .......... .......... 76% 97.8M 1s
    ## 103500K .......... .......... .......... .......... .......... 76%  119M 1s
    ## 103550K .......... .......... .......... .......... .......... 76% 74.3M 1s
    ## 103600K .......... .......... .......... .......... .......... 76% 82.8M 1s
    ## 103650K .......... .......... .......... .......... .......... 76% 85.1M 1s
    ## 103700K .......... .......... .......... .......... .......... 77% 65.1M 1s
    ## 103750K .......... .......... .......... .......... .......... 77% 64.6M 1s
    ## 103800K .......... .......... .......... .......... .......... 77%  101M 1s
    ## 103850K .......... .......... .......... .......... .......... 77% 68.7M 1s
    ## 103900K .......... .......... .......... .......... .......... 77%  147M 1s
    ## 103950K .......... .......... .......... .......... .......... 77%  159M 1s
    ## 104000K .......... .......... .......... .......... .......... 77%  113M 1s
    ## 104050K .......... .......... .......... .......... .......... 77% 95.7M 1s
    ## 104100K .......... .......... .......... .......... .......... 77% 76.8M 1s
    ## 104150K .......... .......... .......... .......... .......... 77% 52.3M 1s
    ## 104200K .......... .......... .......... .......... .......... 77% 55.8M 1s
    ## 104250K .......... .......... .......... .......... .......... 77%  138M 1s
    ## 104300K .......... .......... .......... .......... .......... 77%  111M 1s
    ## 104350K .......... .......... .......... .......... .......... 77%  118M 1s
    ## 104400K .......... .......... .......... .......... .......... 77% 79.7M 1s
    ## 104450K .......... .......... .......... .......... .......... 77%  131M 1s
    ## 104500K .......... .......... .......... .......... .......... 77% 64.7M 1s
    ## 104550K .......... .......... .......... .......... .......... 77% 71.2M 1s
    ## 104600K .......... .......... .......... .......... .......... 77% 85.7M 1s
    ## 104650K .......... .......... .......... .......... .......... 77%  108M 1s
    ## 104700K .......... .......... .......... .......... .......... 77%  133M 1s
    ## 104750K .......... .......... .......... .......... .......... 77% 65.4M 1s
    ## 104800K .......... .......... .......... .......... .......... 77%  127M 1s
    ## 104850K .......... .......... .......... .......... .......... 77% 78.1M 1s
    ## 104900K .......... .......... .......... .......... .......... 77% 31.3M 1s
    ## 104950K .......... .......... .......... .......... .......... 77%  137M 1s
    ## 105000K .......... .......... .......... .......... .......... 77%  161M 1s
    ## 105050K .......... .......... .......... .......... .......... 78%  149M 1s
    ## 105100K .......... .......... .......... .......... .......... 78% 76.0M 1s
    ## 105150K .......... .......... .......... .......... .......... 78%  148M 1s
    ## 105200K .......... .......... .......... .......... .......... 78%  131M 1s
    ## 105250K .......... .......... .......... .......... .......... 78%  124M 1s
    ## 105300K .......... .......... .......... .......... .......... 78%  156M 1s
    ## 105350K .......... .......... .......... .......... .......... 78%  124M 1s
    ## 105400K .......... .......... .......... .......... .......... 78%  129M 1s
    ## 105450K .......... .......... .......... .......... .......... 78%  141M 1s
    ## 105500K .......... .......... .......... .......... .......... 78% 89.8M 1s
    ## 105550K .......... .......... .......... .......... .......... 78%  157M 1s
    ## 105600K .......... .......... .......... .......... .......... 78%  130M 1s
    ## 105650K .......... .......... .......... .......... .......... 78% 42.7M 1s
    ## 105700K .......... .......... .......... .......... .......... 78%  112M 1s
    ## 105750K .......... .......... .......... .......... .......... 78%  143M 1s
    ## 105800K .......... .......... .......... .......... .......... 78%  139M 1s
    ## 105850K .......... .......... .......... .......... .......... 78%  155M 1s
    ## 105900K .......... .......... .......... .......... .......... 78%  127M 1s
    ## 105950K .......... .......... .......... .......... .......... 78%  130M 1s
    ## 106000K .......... .......... .......... .......... .......... 78%  156M 1s
    ## 106050K .......... .......... .......... .......... .......... 78%  148M 1s
    ## 106100K .......... .......... .......... .......... .......... 78%  159M 1s
    ## 106150K .......... .......... .......... .......... .......... 78%  154M 1s
    ## 106200K .......... .......... .......... .......... .......... 78%  117M 1s
    ## 106250K .......... .......... .......... .......... .......... 78%  112M 1s
    ## 106300K .......... .......... .......... .......... .......... 78% 69.5M 1s
    ## 106350K .......... .......... .......... .......... .......... 78%  151M 1s
    ## 106400K .......... .......... .......... .......... .......... 79%  103M 1s
    ## 106450K .......... .......... .......... .......... .......... 79% 80.3M 1s
    ## 106500K .......... .......... .......... .......... .......... 79% 87.8M 1s
    ## 106550K .......... .......... .......... .......... .......... 79% 84.1M 1s
    ## 106600K .......... .......... .......... .......... .......... 79% 56.0M 1s
    ## 106650K .......... .......... .......... .......... .......... 79%  161M 1s
    ## 106700K .......... .......... .......... .......... .......... 79% 80.1M 1s
    ## 106750K .......... .......... .......... .......... .......... 79% 98.4M 1s
    ## 106800K .......... .......... .......... .......... .......... 79%  153M 1s
    ## 106850K .......... .......... .......... .......... .......... 79%  150M 1s
    ## 106900K .......... .......... .......... .......... .......... 79%  153M 1s
    ## 106950K .......... .......... .......... .......... .......... 79%  116M 1s
    ## 107000K .......... .......... .......... .......... .......... 79%  150M 1s
    ## 107050K .......... .......... .......... .......... .......... 79%  162M 0s
    ## 107100K .......... .......... .......... .......... .......... 79%  157M 0s
    ## 107150K .......... .......... .......... .......... .......... 79%  133M 0s
    ## 107200K .......... .......... .......... .......... .......... 79%  131M 0s
    ## 107250K .......... .......... .......... .......... .......... 79%  103M 0s
    ## 107300K .......... .......... .......... .......... .......... 79%  119M 0s
    ## 107350K .......... .......... .......... .......... .......... 79% 28.2M 0s
    ## 107400K .......... .......... .......... .......... .......... 79%  131M 0s
    ## 107450K .......... .......... .......... .......... .......... 79%  123M 0s
    ## 107500K .......... .......... .......... .......... .......... 79%  125M 0s
    ## 107550K .......... .......... .......... .......... .......... 79%  145M 0s
    ## 107600K .......... .......... .......... .......... .......... 79%  166M 0s
    ## 107650K .......... .......... .......... .......... .......... 79% 27.5M 0s
    ## 107700K .......... .......... .......... .......... .......... 79% 40.6M 0s
    ## 107750K .......... .......... .......... .......... .......... 80% 41.2M 0s
    ## 107800K .......... .......... .......... .......... .......... 80% 84.0M 0s
    ## 107850K .......... .......... .......... .......... .......... 80% 44.0M 0s
    ## 107900K .......... .......... .......... .......... .......... 80% 69.7M 0s
    ## 107950K .......... .......... .......... .......... .......... 80%  126M 0s
    ## 108000K .......... .......... .......... .......... .......... 80%  140M 0s
    ## 108050K .......... .......... .......... .......... .......... 80%  145M 0s
    ## 108100K .......... .......... .......... .......... .......... 80% 14.0M 0s
    ## 108150K .......... .......... .......... .......... .......... 80% 44.8M 0s
    ## 108200K .......... .......... .......... .......... .......... 80% 59.7M 0s
    ## 108250K .......... .......... .......... .......... .......... 80% 63.5M 0s
    ## 108300K .......... .......... .......... .......... .......... 80% 58.8M 0s
    ## 108350K .......... .......... .......... .......... .......... 80% 74.4M 0s
    ## 108400K .......... .......... .......... .......... .......... 80%  115M 0s
    ## 108450K .......... .......... .......... .......... .......... 80%  165M 0s
    ## 108500K .......... .......... .......... .......... .......... 80%  123M 0s
    ## 108550K .......... .......... .......... .......... .......... 80% 59.4M 0s
    ## 108600K .......... .......... .......... .......... .......... 80% 83.3M 0s
    ## 108650K .......... .......... .......... .......... .......... 80%  127M 0s
    ## 108700K .......... .......... .......... .......... .......... 80%  141M 0s
    ## 108750K .......... .......... .......... .......... .......... 80% 99.1M 0s
    ## 108800K .......... .......... .......... .......... .......... 80% 65.9M 0s
    ## 108850K .......... .......... .......... .......... .......... 80%  163M 0s
    ## 108900K .......... .......... .......... .......... .......... 80% 25.6M 0s
    ## 108950K .......... .......... .......... .......... .......... 80%  116M 0s
    ## 109000K .......... .......... .......... .......... .......... 80% 40.6M 0s
    ## 109050K .......... .......... .......... .......... .......... 80%  155M 0s
    ## 109100K .......... .......... .......... .......... .......... 81% 58.5M 0s
    ## 109150K .......... .......... .......... .......... .......... 81%  167M 0s
    ## 109200K .......... .......... .......... .......... .......... 81%  112M 0s
    ## 109250K .......... .......... .......... .......... .......... 81%  138M 0s
    ## 109300K .......... .......... .......... .......... .......... 81% 69.6M 0s
    ## 109350K .......... .......... .......... .......... .......... 81%  173M 0s
    ## 109400K .......... .......... .......... .......... .......... 81% 10.5M 0s
    ## 109450K .......... .......... .......... .......... .......... 81%  151M 0s
    ## 109500K .......... .......... .......... .......... .......... 81% 79.1M 0s
    ## 109550K .......... .......... .......... .......... .......... 81% 76.0M 0s
    ## 109600K .......... .......... .......... .......... .......... 81% 43.6M 0s
    ## 109650K .......... .......... .......... .......... .......... 81%  124M 0s
    ## 109700K .......... .......... .......... .......... .......... 81% 93.8M 0s
    ## 109750K .......... .......... .......... .......... .......... 81%  107M 0s
    ## 109800K .......... .......... .......... .......... .......... 81%  131M 0s
    ## 109850K .......... .......... .......... .......... .......... 81%  171M 0s
    ## 109900K .......... .......... .......... .......... .......... 81% 93.5M 0s
    ## 109950K .......... .......... .......... .......... .......... 81%  140M 0s
    ## 110000K .......... .......... .......... .......... .......... 81% 74.7M 0s
    ## 110050K .......... .......... .......... .......... .......... 81% 84.6M 0s
    ## 110100K .......... .......... .......... .......... .......... 81% 36.8M 0s
    ## 110150K .......... .......... .......... .......... .......... 81% 48.4M 0s
    ## 110200K .......... .......... .......... .......... .......... 81%  115M 0s
    ## 110250K .......... .......... .......... .......... .......... 81%  144M 0s
    ## 110300K .......... .......... .......... .......... .......... 81%  105M 0s
    ## 110350K .......... .......... .......... .......... .......... 81%  137M 0s
    ## 110400K .......... .......... .......... .......... .......... 81% 68.1M 0s
    ## 110450K .......... .......... .......... .......... .......... 82%  103M 0s
    ## 110500K .......... .......... .......... .......... .......... 82% 66.6M 0s
    ## 110550K .......... .......... .......... .......... .......... 82%  134M 0s
    ## 110600K .......... .......... .......... .......... .......... 82% 53.9M 0s
    ## 110650K .......... .......... .......... .......... .......... 82% 85.4M 0s
    ## 110700K .......... .......... .......... .......... .......... 82% 89.5M 0s
    ## 110750K .......... .......... .......... .......... .......... 82% 87.6M 0s
    ## 110800K .......... .......... .......... .......... .......... 82% 55.3M 0s
    ## 110850K .......... .......... .......... .......... .......... 82% 87.6M 0s
    ## 110900K .......... .......... .......... .......... .......... 82%  136M 0s
    ## 110950K .......... .......... .......... .......... .......... 82% 83.6M 0s
    ## 111000K .......... .......... .......... .......... .......... 82%  135M 0s
    ## 111050K .......... .......... .......... .......... .......... 82% 70.4M 0s
    ## 111100K .......... .......... .......... .......... .......... 82% 76.6M 0s
    ## 111150K .......... .......... .......... .......... .......... 82%  122M 0s
    ## 111200K .......... .......... .......... .......... .......... 82% 60.0M 0s
    ## 111250K .......... .......... .......... .......... .......... 82%  162M 0s
    ## 111300K .......... .......... .......... .......... .......... 82%  132M 0s
    ## 111350K .......... .......... .......... .......... .......... 82% 82.3M 0s
    ## 111400K .......... .......... .......... .......... .......... 82% 95.6M 0s
    ## 111450K .......... .......... .......... .......... .......... 82% 40.1M 0s
    ## 111500K .......... .......... .......... .......... .......... 82%  128M 0s
    ## 111550K .......... .......... .......... .......... .......... 82%  174M 0s
    ## 111600K .......... .......... .......... .......... .......... 82%  134M 0s
    ## 111650K .......... .......... .......... .......... .......... 82% 17.2M 0s
    ## 111700K .......... .......... .......... .......... .......... 82%  114M 0s
    ## 111750K .......... .......... .......... .......... .......... 82%  133M 0s
    ## 111800K .......... .......... .......... .......... .......... 83%  102M 0s
    ## 111850K .......... .......... .......... .......... .......... 83%  136M 0s
    ## 111900K .......... .......... .......... .......... .......... 83% 87.2M 0s
    ## 111950K .......... .......... .......... .......... .......... 83%  146M 0s
    ## 112000K .......... .......... .......... .......... .......... 83%  138M 0s
    ## 112050K .......... .......... .......... .......... .......... 83% 55.9M 0s
    ## 112100K .......... .......... .......... .......... .......... 83% 51.4M 0s
    ## 112150K .......... .......... .......... .......... .......... 83% 33.6M 0s
    ## 112200K .......... .......... .......... .......... .......... 83% 24.8M 0s
    ## 112250K .......... .......... .......... .......... .......... 83% 47.1M 0s
    ## 112300K .......... .......... .......... .......... .......... 83% 47.2M 0s
    ## 112350K .......... .......... .......... .......... .......... 83% 89.4M 0s
    ## 112400K .......... .......... .......... .......... .......... 83% 50.2M 0s
    ## 112450K .......... .......... .......... .......... .......... 83% 45.7M 0s
    ## 112500K .......... .......... .......... .......... .......... 83% 58.6M 0s
    ## 112550K .......... .......... .......... .......... .......... 83% 99.8M 0s
    ## 112600K .......... .......... .......... .......... .......... 83% 41.3M 0s
    ## 112650K .......... .......... .......... .......... .......... 83% 90.8M 0s
    ## 112700K .......... .......... .......... .......... .......... 83% 61.8M 0s
    ## 112750K .......... .......... .......... .......... .......... 83% 46.6M 0s
    ## 112800K .......... .......... .......... .......... .......... 83% 52.1M 0s
    ## 112850K .......... .......... .......... .......... .......... 83%  113M 0s
    ## 112900K .......... .......... .......... .......... .......... 83% 36.4M 0s
    ## 112950K .......... .......... .......... .......... .......... 83% 41.1M 0s
    ## 113000K .......... .......... .......... .......... .......... 83% 60.1M 0s
    ## 113050K .......... .......... .......... .......... .......... 83% 58.4M 0s
    ## 113100K .......... .......... .......... .......... .......... 83% 60.8M 0s
    ## 113150K .......... .......... .......... .......... .......... 84% 37.7M 0s
    ## 113200K .......... .......... .......... .......... .......... 84% 26.8M 0s
    ## 113250K .......... .......... .......... .......... .......... 84% 76.3M 0s
    ## 113300K .......... .......... .......... .......... .......... 84% 85.7M 0s
    ## 113350K .......... .......... .......... .......... .......... 84% 48.1M 0s
    ## 113400K .......... .......... .......... .......... .......... 84% 71.7M 0s
    ## 113450K .......... .......... .......... .......... .......... 84% 39.8M 0s
    ## 113500K .......... .......... .......... .......... .......... 84% 49.4M 0s
    ## 113550K .......... .......... .......... .......... .......... 84% 43.4M 0s
    ## 113600K .......... .......... .......... .......... .......... 84%  103M 0s
    ## 113650K .......... .......... .......... .......... .......... 84% 57.4M 0s
    ## 113700K .......... .......... .......... .......... .......... 84% 76.8M 0s
    ## 113750K .......... .......... .......... .......... .......... 84% 53.1M 0s
    ## 113800K .......... .......... .......... .......... .......... 84% 40.1M 0s
    ## 113850K .......... .......... .......... .......... .......... 84% 43.3M 0s
    ## 113900K .......... .......... .......... .......... .......... 84% 58.8M 0s
    ## 113950K .......... .......... .......... .......... .......... 84% 54.4M 0s
    ## 114000K .......... .......... .......... .......... .......... 84% 76.5M 0s
    ## 114050K .......... .......... .......... .......... .......... 84% 55.0M 0s
    ## 114100K .......... .......... .......... .......... .......... 84% 81.6M 0s
    ## 114150K .......... .......... .......... .......... .......... 84% 70.7M 0s
    ## 114200K .......... .......... .......... .......... .......... 84% 92.2M 0s
    ## 114250K .......... .......... .......... .......... .......... 84%  107M 0s
    ## 114300K .......... .......... .......... .......... .......... 84% 83.6M 0s
    ## 114350K .......... .......... .......... .......... .......... 84% 43.4M 0s
    ## 114400K .......... .......... .......... .......... .......... 84% 44.9M 0s
    ## 114450K .......... .......... .......... .......... .......... 84%  106M 0s
    ## 114500K .......... .......... .......... .......... .......... 85% 68.9M 0s
    ## 114550K .......... .......... .......... .......... .......... 85% 54.9M 0s
    ## 114600K .......... .......... .......... .......... .......... 85% 33.1M 0s
    ## 114650K .......... .......... .......... .......... .......... 85% 90.7M 0s
    ## 114700K .......... .......... .......... .......... .......... 85% 80.8M 0s
    ## 114750K .......... .......... .......... .......... .......... 85% 60.6M 0s
    ## 114800K .......... .......... .......... .......... .......... 85% 55.0M 0s
    ## 114850K .......... .......... .......... .......... .......... 85% 91.2M 0s
    ## 114900K .......... .......... .......... .......... .......... 85% 90.2M 0s
    ## 114950K .......... .......... .......... .......... .......... 85% 81.8M 0s
    ## 115000K .......... .......... .......... .......... .......... 85% 78.5M 0s
    ## 115050K .......... .......... .......... .......... .......... 85% 82.3M 0s
    ## 115100K .......... .......... .......... .......... .......... 85% 69.2M 0s
    ## 115150K .......... .......... .......... .......... .......... 85% 48.1M 0s
    ## 115200K .......... .......... .......... .......... .......... 85% 72.2M 0s
    ## 115250K .......... .......... .......... .......... .......... 85% 46.4M 0s
    ## 115300K .......... .......... .......... .......... .......... 85% 50.0M 0s
    ## 115350K .......... .......... .......... .......... .......... 85% 45.1M 0s
    ## 115400K .......... .......... .......... .......... .......... 85% 76.3M 0s
    ## 115450K .......... .......... .......... .......... .......... 85% 52.9M 0s
    ## 115500K .......... .......... .......... .......... .......... 85% 69.3M 0s
    ## 115550K .......... .......... .......... .......... .......... 85% 57.9M 0s
    ## 115600K .......... .......... .......... .......... .......... 85% 94.6M 0s
    ## 115650K .......... .......... .......... .......... .......... 85% 70.4M 0s
    ## 115700K .......... .......... .......... .......... .......... 85% 76.1M 0s
    ## 115750K .......... .......... .......... .......... .......... 85%  102M 0s
    ## 115800K .......... .......... .......... .......... .......... 85% 77.4M 0s
    ## 115850K .......... .......... .......... .......... .......... 86% 75.9M 0s
    ## 115900K .......... .......... .......... .......... .......... 86% 42.7M 0s
    ## 115950K .......... .......... .......... .......... .......... 86% 59.9M 0s
    ## 116000K .......... .......... .......... .......... .......... 86%  113M 0s
    ## 116050K .......... .......... .......... .......... .......... 86%  100M 0s
    ## 116100K .......... .......... .......... .......... .......... 86% 80.0M 0s
    ## 116150K .......... .......... .......... .......... .......... 86%  113M 0s
    ## 116200K .......... .......... .......... .......... .......... 86% 46.7M 0s
    ## 116250K .......... .......... .......... .......... .......... 86% 72.0M 0s
    ## 116300K .......... .......... .......... .......... .......... 86% 69.3M 0s
    ## 116350K .......... .......... .......... .......... .......... 86%  109M 0s
    ## 116400K .......... .......... .......... .......... .......... 86% 86.0M 0s
    ## 116450K .......... .......... .......... .......... .......... 86%  109M 0s
    ## 116500K .......... .......... .......... .......... .......... 86% 91.1M 0s
    ## 116550K .......... .......... .......... .......... .......... 86% 61.8M 0s
    ## 116600K .......... .......... .......... .......... .......... 86% 86.1M 0s
    ## 116650K .......... .......... .......... .......... .......... 86% 41.6M 0s
    ## 116700K .......... .......... .......... .......... .......... 86% 98.1M 0s
    ## 116750K .......... .......... .......... .......... .......... 86%  121M 0s
    ## 116800K .......... .......... .......... .......... .......... 86% 93.9M 0s
    ## 116850K .......... .......... .......... .......... .......... 86% 87.1M 0s
    ## 116900K .......... .......... .......... .......... .......... 86% 44.1M 0s
    ## 116950K .......... .......... .......... .......... .......... 86%  105M 0s
    ## 117000K .......... .......... .......... .......... .......... 86% 41.7M 0s
    ## 117050K .......... .......... .......... .......... .......... 86%  129M 0s
    ## 117100K .......... .......... .......... .......... .......... 86% 76.3M 0s
    ## 117150K .......... .......... .......... .......... .......... 86%  137M 0s
    ## 117200K .......... .......... .......... .......... .......... 87% 71.0M 0s
    ## 117250K .......... .......... .......... .......... .......... 87% 54.4M 0s
    ## 117300K .......... .......... .......... .......... .......... 87% 55.1M 0s
    ## 117350K .......... .......... .......... .......... .......... 87%  117M 0s
    ## 117400K .......... .......... .......... .......... .......... 87% 42.0M 0s
    ## 117450K .......... .......... .......... .......... .......... 87% 88.8M 0s
    ## 117500K .......... .......... .......... .......... .......... 87% 68.3M 0s
    ## 117550K .......... .......... .......... .......... .......... 87% 75.6M 0s
    ## 117600K .......... .......... .......... .......... .......... 87% 97.0M 0s
    ## 117650K .......... .......... .......... .......... .......... 87%  107M 0s
    ## 117700K .......... .......... .......... .......... .......... 87% 66.7M 0s
    ## 117750K .......... .......... .......... .......... .......... 87% 78.2M 0s
    ## 117800K .......... .......... .......... .......... .......... 87% 58.3M 0s
    ## 117850K .......... .......... .......... .......... .......... 87% 65.6M 0s
    ## 117900K .......... .......... .......... .......... .......... 87% 53.7M 0s
    ## 117950K .......... .......... .......... .......... .......... 87%  128M 0s
    ## 118000K .......... .......... .......... .......... .......... 87% 68.1M 0s
    ## 118050K .......... .......... .......... .......... .......... 87%  123M 0s
    ## 118100K .......... .......... .......... .......... .......... 87% 82.6M 0s
    ## 118150K .......... .......... .......... .......... .......... 87%  125M 0s
    ## 118200K .......... .......... .......... .......... .......... 87% 55.8M 0s
    ## 118250K .......... .......... .......... .......... .......... 87%  106M 0s
    ## 118300K .......... .......... .......... .......... .......... 87% 97.8M 0s
    ## 118350K .......... .......... .......... .......... .......... 87%  124M 0s
    ## 118400K .......... .......... .......... .......... .......... 87% 92.3M 0s
    ## 118450K .......... .......... .......... .......... .......... 87% 97.0M 0s
    ## 118500K .......... .......... .......... .......... .......... 87%  111M 0s
    ## 118550K .......... .......... .......... .......... .......... 88%  131M 0s
    ## 118600K .......... .......... .......... .......... .......... 88% 61.8M 0s
    ## 118650K .......... .......... .......... .......... .......... 88% 78.1M 0s
    ## 118700K .......... .......... .......... .......... .......... 88%  108M 0s
    ## 118750K .......... .......... .......... .......... .......... 88% 61.1M 0s
    ## 118800K .......... .......... .......... .......... .......... 88% 71.0M 0s
    ## 118850K .......... .......... .......... .......... .......... 88% 72.1M 0s
    ## 118900K .......... .......... .......... .......... .......... 88%  102M 0s
    ## 118950K .......... .......... .......... .......... .......... 88% 41.1M 0s
    ## 119000K .......... .......... .......... .......... .......... 88%  107M 0s
    ## 119050K .......... .......... .......... .......... .......... 88% 54.6M 0s
    ## 119100K .......... .......... .......... .......... .......... 88% 65.3M 0s
    ## 119150K .......... .......... .......... .......... .......... 88%  143M 0s
    ## 119200K .......... .......... .......... .......... .......... 88% 88.5M 0s
    ## 119250K .......... .......... .......... .......... .......... 88% 85.6M 0s
    ## 119300K .......... .......... .......... .......... .......... 88% 56.7M 0s
    ## 119350K .......... .......... .......... .......... .......... 88% 88.4M 0s
    ## 119400K .......... .......... .......... .......... .......... 88% 77.8M 0s
    ## 119450K .......... .......... .......... .......... .......... 88%  159M 0s
    ## 119500K .......... .......... .......... .......... .......... 88% 39.5M 0s
    ## 119550K .......... .......... .......... .......... .......... 88%  123M 0s
    ## 119600K .......... .......... .......... .......... .......... 88% 43.8M 0s
    ## 119650K .......... .......... .......... .......... .......... 88%  115M 0s
    ## 119700K .......... .......... .......... .......... .......... 88%  108M 0s
    ## 119750K .......... .......... .......... .......... .......... 88%  120M 0s
    ## 119800K .......... .......... .......... .......... .......... 88% 70.8M 0s
    ## 119850K .......... .......... .......... .......... .......... 88%  103M 0s
    ## 119900K .......... .......... .......... .......... .......... 89% 78.7M 0s
    ## 119950K .......... .......... .......... .......... .......... 89% 99.3M 0s
    ## 120000K .......... .......... .......... .......... .......... 89%  105M 0s
    ## 120050K .......... .......... .......... .......... .......... 89% 88.9M 0s
    ## 120100K .......... .......... .......... .......... .......... 89% 81.0M 0s
    ## 120150K .......... .......... .......... .......... .......... 89% 85.9M 0s
    ## 120200K .......... .......... .......... .......... .......... 89%  114M 0s
    ## 120250K .......... .......... .......... .......... .......... 89% 48.8M 0s
    ## 120300K .......... .......... .......... .......... .......... 89% 70.5M 0s
    ## 120350K .......... .......... .......... .......... .......... 89% 87.2M 0s
    ## 120400K .......... .......... .......... .......... .......... 89% 53.5M 0s
    ## 120450K .......... .......... .......... .......... .......... 89%  134M 0s
    ## 120500K .......... .......... .......... .......... .......... 89%  109M 0s
    ## 120550K .......... .......... .......... .......... .......... 89%  118M 0s
    ## 120600K .......... .......... .......... .......... .......... 89%  106M 0s
    ## 120650K .......... .......... .......... .......... .......... 89% 69.4M 0s
    ## 120700K .......... .......... .......... .......... .......... 89%  110M 0s
    ## 120750K .......... .......... .......... .......... .......... 89%  135M 0s
    ## 120800K .......... .......... .......... .......... .......... 89%  132M 0s
    ## 120850K .......... .......... .......... .......... .......... 89%  156M 0s
    ## 120900K .......... .......... .......... .......... .......... 89%  129M 0s
    ## 120950K .......... .......... .......... .......... .......... 89% 70.1M 0s
    ## 121000K .......... .......... .......... .......... .......... 89% 93.2M 0s
    ## 121050K .......... .......... .......... .......... .......... 89% 69.6M 0s
    ## 121100K .......... .......... .......... .......... .......... 89%  107M 0s
    ## 121150K .......... .......... .......... .......... .......... 89% 89.1M 0s
    ## 121200K .......... .......... .......... .......... .......... 89%  134M 0s
    ## 121250K .......... .......... .......... .......... .......... 90%  126M 0s
    ## 121300K .......... .......... .......... .......... .......... 90% 79.2M 0s
    ## 121350K .......... .......... .......... .......... .......... 90%  130M 0s
    ## 121400K .......... .......... .......... .......... .......... 90%  108M 0s
    ## 121450K .......... .......... .......... .......... .......... 90% 62.9M 0s
    ## 121500K .......... .......... .......... .......... .......... 90% 84.2M 0s
    ## 121550K .......... .......... .......... .......... .......... 90% 38.5M 0s
    ## 121600K .......... .......... .......... .......... .......... 90%  130M 0s
    ## 121650K .......... .......... .......... .......... .......... 90% 30.5M 0s
    ## 121700K .......... .......... .......... .......... .......... 90% 25.8M 0s
    ## 121750K .......... .......... .......... .......... .......... 90% 77.8M 0s
    ## 121800K .......... .......... .......... .......... .......... 90% 47.4M 0s
    ## 121850K .......... .......... .......... .......... .......... 90% 53.2M 0s
    ## 121900K .......... .......... .......... .......... .......... 90% 54.1M 0s
    ## 121950K .......... .......... .......... .......... .......... 90% 69.9M 0s
    ## 122000K .......... .......... .......... .......... .......... 90% 59.0M 0s
    ## 122050K .......... .......... .......... .......... .......... 90% 85.1M 0s
    ## 122100K .......... .......... .......... .......... .......... 90% 96.9M 0s
    ## 122150K .......... .......... .......... .......... .......... 90% 83.8M 0s
    ## 122200K .......... .......... .......... .......... .......... 90% 91.3M 0s
    ## 122250K .......... .......... .......... .......... .......... 90% 63.4M 0s
    ## 122300K .......... .......... .......... .......... .......... 90% 66.9M 0s
    ## 122350K .......... .......... .......... .......... .......... 90% 96.7M 0s
    ## 122400K .......... .......... .......... .......... .......... 90% 96.8M 0s
    ## 122450K .......... .......... .......... .......... .......... 90% 57.5M 0s
    ## 122500K .......... .......... .......... .......... .......... 90% 81.1M 0s
    ## 122550K .......... .......... .......... .......... .......... 90%  107M 0s
    ## 122600K .......... .......... .......... .......... .......... 91% 61.7M 0s
    ## 122650K .......... .......... .......... .......... .......... 91% 95.6M 0s
    ## 122700K .......... .......... .......... .......... .......... 91% 33.7M 0s
    ## 122750K .......... .......... .......... .......... .......... 91% 60.8M 0s
    ## 122800K .......... .......... .......... .......... .......... 91% 46.7M 0s
    ## 122850K .......... .......... .......... .......... .......... 91% 59.1M 0s
    ## 122900K .......... .......... .......... .......... .......... 91% 71.1M 0s
    ## 122950K .......... .......... .......... .......... .......... 91% 65.5M 0s
    ## 123000K .......... .......... .......... .......... .......... 91% 56.4M 0s
    ## 123050K .......... .......... .......... .......... .......... 91% 98.0M 0s
    ## 123100K .......... .......... .......... .......... .......... 91% 80.9M 0s
    ## 123150K .......... .......... .......... .......... .......... 91%  109M 0s
    ## 123200K .......... .......... .......... .......... .......... 91% 63.2M 0s
    ## 123250K .......... .......... .......... .......... .......... 91%  105M 0s
    ## 123300K .......... .......... .......... .......... .......... 91% 64.6M 0s
    ## 123350K .......... .......... .......... .......... .......... 91%  111M 0s
    ## 123400K .......... .......... .......... .......... .......... 91% 80.2M 0s
    ## 123450K .......... .......... .......... .......... .......... 91% 88.9M 0s
    ## 123500K .......... .......... .......... .......... .......... 91%  101M 0s
    ## 123550K .......... .......... .......... .......... .......... 91% 99.3M 0s
    ## 123600K .......... .......... .......... .......... .......... 91%  104M 0s
    ## 123650K .......... .......... .......... .......... .......... 91%  125M 0s
    ## 123700K .......... .......... .......... .......... .......... 91%  108M 0s
    ## 123750K .......... .......... .......... .......... .......... 91% 92.3M 0s
    ## 123800K .......... .......... .......... .......... .......... 91% 64.8M 0s
    ## 123850K .......... .......... .......... .......... .......... 91% 42.7M 0s
    ## 123900K .......... .......... .......... .......... .......... 91% 96.6M 0s
    ## 123950K .......... .......... .......... .......... .......... 92%  118M 0s
    ## 124000K .......... .......... .......... .......... .......... 92%  121M 0s
    ## 124050K .......... .......... .......... .......... .......... 92%  106M 0s
    ## 124100K .......... .......... .......... .......... .......... 92%  111M 0s
    ## 124150K .......... .......... .......... .......... .......... 92% 93.0M 0s
    ## 124200K .......... .......... .......... .......... .......... 92% 75.0M 0s
    ## 124250K .......... .......... .......... .......... .......... 92%  122M 0s
    ## 124300K .......... .......... .......... .......... .......... 92% 97.0M 0s
    ## 124350K .......... .......... .......... .......... .......... 92% 64.6M 0s
    ## 124400K .......... .......... .......... .......... .......... 92% 83.4M 0s
    ## 124450K .......... .......... .......... .......... .......... 92% 56.0M 0s
    ## 124500K .......... .......... .......... .......... .......... 92% 58.1M 0s
    ## 124550K .......... .......... .......... .......... .......... 92% 86.6M 0s
    ## 124600K .......... .......... .......... .......... .......... 92%  108M 0s
    ## 124650K .......... .......... .......... .......... .......... 92%  114M 0s
    ## 124700K .......... .......... .......... .......... .......... 92% 52.5M 0s
    ## 124750K .......... .......... .......... .......... .......... 92%  112M 0s
    ## 124800K .......... .......... .......... .......... .......... 92% 86.9M 0s
    ## 124850K .......... .......... .......... .......... .......... 92% 81.9M 0s
    ## 124900K .......... .......... .......... .......... .......... 92% 81.8M 0s
    ## 124950K .......... .......... .......... .......... .......... 92%  112M 0s
    ## 125000K .......... .......... .......... .......... .......... 92% 33.3M 0s
    ## 125050K .......... .......... .......... .......... .......... 92%  114M 0s
    ## 125100K .......... .......... .......... .......... .......... 92% 74.2M 0s
    ## 125150K .......... .......... .......... .......... .......... 92% 67.1M 0s
    ## 125200K .......... .......... .......... .......... .......... 92% 74.8M 0s
    ## 125250K .......... .......... .......... .......... .......... 92%  133M 0s
    ## 125300K .......... .......... .......... .......... .......... 93%  123M 0s
    ## 125350K .......... .......... .......... .......... .......... 93%  138M 0s
    ## 125400K .......... .......... .......... .......... .......... 93%  106M 0s
    ## 125450K .......... .......... .......... .......... .......... 93% 52.0M 0s
    ## 125500K .......... .......... .......... .......... .......... 93% 84.1M 0s
    ## 125550K .......... .......... .......... .......... .......... 93% 85.0M 0s
    ## 125600K .......... .......... .......... .......... .......... 93% 41.4M 0s
    ## 125650K .......... .......... .......... .......... .......... 93%  103M 0s
    ## 125700K .......... .......... .......... .......... .......... 93%  107M 0s
    ## 125750K .......... .......... .......... .......... .......... 93% 74.1M 0s
    ## 125800K .......... .......... .......... .......... .......... 93% 99.8M 0s
    ## 125850K .......... .......... .......... .......... .......... 93% 76.7M 0s
    ## 125900K .......... .......... .......... .......... .......... 93% 96.0M 0s
    ## 125950K .......... .......... .......... .......... .......... 93% 78.2M 0s
    ## 126000K .......... .......... .......... .......... .......... 93% 79.1M 0s
    ## 126050K .......... .......... .......... .......... .......... 93% 50.9M 0s
    ## 126100K .......... .......... .......... .......... .......... 93% 77.8M 0s
    ## 126150K .......... .......... .......... .......... .......... 93%  129M 0s
    ## 126200K .......... .......... .......... .......... .......... 93%  111M 0s
    ## 126250K .......... .......... .......... .......... .......... 93%  116M 0s
    ## 126300K .......... .......... .......... .......... .......... 93%  117M 0s
    ## 126350K .......... .......... .......... .......... .......... 93%  111M 0s
    ## 126400K .......... .......... .......... .......... .......... 93%  123M 0s
    ## 126450K .......... .......... .......... .......... .......... 93%  130M 0s
    ## 126500K .......... .......... .......... .......... .......... 93%  126M 0s
    ## 126550K .......... .......... .......... .......... .......... 93%  130M 0s
    ## 126600K .......... .......... .......... .......... .......... 93% 39.5M 0s
    ## 126650K .......... .......... .......... .......... .......... 94% 79.4M 0s
    ## 126700K .......... .......... .......... .......... .......... 94%  134M 0s
    ## 126750K .......... .......... .......... .......... .......... 94%  126M 0s
    ## 126800K .......... .......... .......... .......... .......... 94%  126M 0s
    ## 126850K .......... .......... .......... .......... .......... 94% 99.9M 0s
    ## 126900K .......... .......... .......... .......... .......... 94%  103M 0s
    ## 126950K .......... .......... .......... .......... .......... 94%  136M 0s
    ## 127000K .......... .......... .......... .......... .......... 94%  139M 0s
    ## 127050K .......... .......... .......... .......... .......... 94%  125M 0s
    ## 127100K .......... .......... .......... .......... .......... 94% 78.2M 0s
    ## 127150K .......... .......... .......... .......... .......... 94%  108M 0s
    ## 127200K .......... .......... .......... .......... .......... 94% 67.8M 0s
    ## 127250K .......... .......... .......... .......... .......... 94%  104M 0s
    ## 127300K .......... .......... .......... .......... .......... 94% 44.1M 0s
    ## 127350K .......... .......... .......... .......... .......... 94% 66.3M 0s
    ## 127400K .......... .......... .......... .......... .......... 94%  133M 0s
    ## 127450K .......... .......... .......... .......... .......... 94%  138M 0s
    ## 127500K .......... .......... .......... .......... .......... 94%  114M 0s
    ## 127550K .......... .......... .......... .......... .......... 94% 66.9M 0s
    ## 127600K .......... .......... .......... .......... .......... 94% 84.9M 0s
    ## 127650K .......... .......... .......... .......... .......... 94% 54.3M 0s
    ## 127700K .......... .......... .......... .......... .......... 94% 76.1M 0s
    ## 127750K .......... .......... .......... .......... .......... 94% 72.8M 0s
    ## 127800K .......... .......... .......... .......... .......... 94% 70.4M 0s
    ## 127850K .......... .......... .......... .......... .......... 94%  152M 0s
    ## 127900K .......... .......... .......... .......... .......... 94% 74.6M 0s
    ## 127950K .......... .......... .......... .......... .......... 94%  113M 0s
    ## 128000K .......... .......... .......... .......... .......... 95% 68.7M 0s
    ## 128050K .......... .......... .......... .......... .......... 95%  153M 0s
    ## 128100K .......... .......... .......... .......... .......... 95% 88.4M 0s
    ## 128150K .......... .......... .......... .......... .......... 95% 61.2M 0s
    ## 128200K .......... .......... .......... .......... .......... 95% 78.0M 0s
    ## 128250K .......... .......... .......... .......... .......... 95% 30.7M 0s
    ## 128300K .......... .......... .......... .......... .......... 95%  122M 0s
    ## 128350K .......... .......... .......... .......... .......... 95%  146M 0s
    ## 128400K .......... .......... .......... .......... .......... 95%  115M 0s
    ## 128450K .......... .......... .......... .......... .......... 95%  139M 0s
    ## 128500K .......... .......... .......... .......... .......... 95%  124M 0s
    ## 128550K .......... .......... .......... .......... .......... 95%  175M 0s
    ## 128600K .......... .......... .......... .......... .......... 95%  138M 0s
    ## 128650K .......... .......... .......... .......... .......... 95%  141M 0s
    ## 128700K .......... .......... .......... .......... .......... 95% 39.6M 0s
    ## 128750K .......... .......... .......... .......... .......... 95%  129M 0s
    ## 128800K .......... .......... .......... .......... .......... 95%  128M 0s
    ## 128850K .......... .......... .......... .......... .......... 95%  105M 0s
    ## 128900K .......... .......... .......... .......... .......... 95% 89.5M 0s
    ## 128950K .......... .......... .......... .......... .......... 95%  141M 0s
    ## 129000K .......... .......... .......... .......... .......... 95%  127M 0s
    ## 129050K .......... .......... .......... .......... .......... 95%  153M 0s
    ## 129100K .......... .......... .......... .......... .......... 95%  126M 0s
    ## 129150K .......... .......... .......... .......... .......... 95% 94.1M 0s
    ## 129200K .......... .......... .......... .......... .......... 95%  117M 0s
    ## 129250K .......... .......... .......... .......... .......... 95%  129M 0s
    ## 129300K .......... .......... .......... .......... .......... 95%  131M 0s
    ## 129350K .......... .......... .......... .......... .......... 96%  144M 0s
    ## 129400K .......... .......... .......... .......... .......... 96% 49.8M 0s
    ## 129450K .......... .......... .......... .......... .......... 96%  150M 0s
    ## 129500K .......... .......... .......... .......... .......... 96%  140M 0s
    ## 129550K .......... .......... .......... .......... .......... 96%  141M 0s
    ## 129600K .......... .......... .......... .......... .......... 96%  142M 0s
    ## 129650K .......... .......... .......... .......... .......... 96%  151M 0s
    ## 129700K .......... .......... .......... .......... .......... 96%  146M 0s
    ## 129750K .......... .......... .......... .......... .......... 96%  169M 0s
    ## 129800K .......... .......... .......... .......... .......... 96%  125M 0s
    ## 129850K .......... .......... .......... .......... .......... 96%  132M 0s
    ## 129900K .......... .......... .......... .......... .......... 96%  144M 0s
    ## 129950K .......... .......... .......... .......... .......... 96%  183M 0s
    ## 130000K .......... .......... .......... .......... .......... 96% 42.4M 0s
    ## 130050K .......... .......... .......... .......... .......... 96% 90.2M 0s
    ## 130100K .......... .......... .......... .......... .......... 96% 36.5M 0s
    ## 130150K .......... .......... .......... .......... .......... 96% 93.8M 0s
    ## 130200K .......... .......... .......... .......... .......... 96% 74.3M 0s
    ## 130250K .......... .......... .......... .......... .......... 96% 64.7M 0s
    ## 130300K .......... .......... .......... .......... .......... 96% 51.9M 0s
    ## 130350K .......... .......... .......... .......... .......... 96% 59.9M 0s
    ## 130400K .......... .......... .......... .......... .......... 96% 39.3M 0s
    ## 130450K .......... .......... .......... .......... .......... 96% 48.0M 0s
    ## 130500K .......... .......... .......... .......... .......... 96% 46.9M 0s
    ## 130550K .......... .......... .......... .......... .......... 96% 52.9M 0s
    ## 130600K .......... .......... .......... .......... .......... 96% 79.6M 0s
    ## 130650K .......... .......... .......... .......... .......... 97%  126M 0s
    ## 130700K .......... .......... .......... .......... .......... 97% 38.5M 0s
    ## 130750K .......... .......... .......... .......... .......... 97%  107M 0s
    ## 130800K .......... .......... .......... .......... .......... 97% 36.6M 0s
    ## 130850K .......... .......... .......... .......... .......... 97% 23.6M 0s
    ## 130900K .......... .......... .......... .......... .......... 97% 78.3M 0s
    ## 130950K .......... .......... .......... .......... .......... 97% 31.9M 0s
    ## 131000K .......... .......... .......... .......... .......... 97% 45.9M 0s
    ## 131050K .......... .......... .......... .......... .......... 97% 56.8M 0s
    ## 131100K .......... .......... .......... .......... .......... 97% 43.4M 0s
    ## 131150K .......... .......... .......... .......... .......... 97%  104M 0s
    ## 131200K .......... .......... .......... .......... .......... 97% 72.0M 0s
    ## 131250K .......... .......... .......... .......... .......... 97% 37.9M 0s
    ## 131300K .......... .......... .......... .......... .......... 97% 36.6M 0s
    ## 131350K .......... .......... .......... .......... .......... 97% 66.1M 0s
    ## 131400K .......... .......... .......... .......... .......... 97% 83.3M 0s
    ## 131450K .......... .......... .......... .......... .......... 97%  124M 0s
    ## 131500K .......... .......... .......... .......... .......... 97% 59.7M 0s
    ## 131550K .......... .......... .......... .......... .......... 97% 69.8M 0s
    ## 131600K .......... .......... .......... .......... .......... 97% 63.3M 0s
    ## 131650K .......... .......... .......... .......... .......... 97% 24.5M 0s
    ## 131700K .......... .......... .......... .......... .......... 97% 69.1M 0s
    ## 131750K .......... .......... .......... .......... .......... 97%  103M 0s
    ## 131800K .......... .......... .......... .......... .......... 97% 83.4M 0s
    ## 131850K .......... .......... .......... .......... .......... 97%  130M 0s
    ## 131900K .......... .......... .......... .......... .......... 97% 62.5M 0s
    ## 131950K .......... .......... .......... .......... .......... 97% 55.8M 0s
    ## 132000K .......... .......... .......... .......... .......... 98% 25.9M 0s
    ## 132050K .......... .......... .......... .......... .......... 98%  110M 0s
    ## 132100K .......... .......... .......... .......... .......... 98% 75.1M 0s
    ## 132150K .......... .......... .......... .......... .......... 98%  101M 0s
    ## 132200K .......... .......... .......... .......... .......... 98%  102M 0s
    ## 132250K .......... .......... .......... .......... .......... 98% 98.3M 0s
    ## 132300K .......... .......... .......... .......... .......... 98%  120M 0s
    ## 132350K .......... .......... .......... .......... .......... 98% 63.6M 0s
    ## 132400K .......... .......... .......... .......... .......... 98%  121M 0s
    ## 132450K .......... .......... .......... .......... .......... 98%  118M 0s
    ## 132500K .......... .......... .......... .......... .......... 98% 57.6M 0s
    ## 132550K .......... .......... .......... .......... .......... 98% 96.0M 0s
    ## 132600K .......... .......... .......... .......... .......... 98%  105M 0s
    ## 132650K .......... .......... .......... .......... .......... 98%  119M 0s
    ## 132700K .......... .......... .......... .......... .......... 98%  128M 0s
    ## 132750K .......... .......... .......... .......... .......... 98%  120M 0s
    ## 132800K .......... .......... .......... .......... .......... 98%  114M 0s
    ## 132850K .......... .......... .......... .......... .......... 98%  108M 0s
    ## 132900K .......... .......... .......... .......... .......... 98% 55.6M 0s
    ## 132950K .......... .......... .......... .......... .......... 98% 77.9M 0s
    ## 133000K .......... .......... .......... .......... .......... 98% 82.2M 0s
    ## 133050K .......... .......... .......... .......... .......... 98% 67.9M 0s
    ## 133100K .......... .......... .......... .......... .......... 98% 45.7M 0s
    ## 133150K .......... .......... .......... .......... .......... 98% 47.9M 0s
    ## 133200K .......... .......... .......... .......... .......... 98% 69.5M 0s
    ## 133250K .......... .......... .......... .......... .......... 98% 74.4M 0s
    ## 133300K .......... .......... .......... .......... .......... 98%  115M 0s
    ## 133350K .......... .......... .......... .......... .......... 99%  103M 0s
    ## 133400K .......... .......... .......... .......... .......... 99%  114M 0s
    ## 133450K .......... .......... .......... .......... .......... 99%  101M 0s
    ## 133500K .......... .......... .......... .......... .......... 99%  114M 0s
    ## 133550K .......... .......... .......... .......... .......... 99%  116M 0s
    ## 133600K .......... .......... .......... .......... .......... 99%  107M 0s
    ## 133650K .......... .......... .......... .......... .......... 99%  107M 0s
    ## 133700K .......... .......... .......... .......... .......... 99%  124M 0s
    ## 133750K .......... .......... .......... .......... .......... 99%  126M 0s
    ## 133800K .......... .......... .......... .......... .......... 99%  115M 0s
    ## 133850K .......... .......... .......... .......... .......... 99%  114M 0s
    ## 133900K .......... .......... .......... .......... .......... 99%  117M 0s
    ## 133950K .......... .......... .......... .......... .......... 99%  121M 0s
    ## 134000K .......... .......... .......... .......... .......... 99%  116M 0s
    ## 134050K .......... .......... .......... .......... .......... 99%  118M 0s
    ## 134100K .......... .......... .......... .......... .......... 99%  109M 0s
    ## 134150K .......... .......... .......... .......... .......... 99%  136M 0s
    ## 134200K .......... .......... .......... .......... .......... 99%  111M 0s
    ## 134250K .......... .......... .......... .......... .......... 99%  106M 0s
    ## 134300K .......... .......... .......... .......... .......... 99%  121M 0s
    ## 134350K .......... .......... .......... .......... .......... 99%  135M 0s
    ## 134400K .......... .......... .......... .......... .......... 99%  129M 0s
    ## 134450K .......... .......... .......... .......... .......... 99%  142M 0s
    ## 134500K .......... .......... .......... .......... .......... 99%  114M 0s
    ## 134550K .......... .......... .......... .......... .......... 99%  117M 0s
    ## 134600K .......... .......... .......... .......... .......... 99%  117M 0s
    ## 134650K .......... .......... .......... .......... .......... 99%  131M 0s
    ## 134700K .......... .......... .......... ..........           100%  128M=2.3s
    ## 
    ## 2020-12-03 15:14:10 (57.2 MB/s) - ‘silva_nr99_v138_train_set.fa.gz.5’ saved [137973851/137973851]

\#Assignation taxonomique

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

## Assignation taxonomique n°2 Silva species assignement

``` bash
wget https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
```

    ## --2020-12-03 15:16:49--  https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 81840166 (78M) [application/octet-stream]
    ## Saving to: ‘silva_species_assignment_v138.fa.gz.2’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 5.46M 14s
    ##     50K .......... .......... .......... .......... ..........  0% 24.2M 9s
    ##    100K .......... .......... .......... .......... ..........  0% 8.43M 9s
    ##    150K .......... .......... .......... .......... ..........  0% 20.6M 8s
    ##    200K .......... .......... .......... .......... ..........  0% 31.2M 7s
    ##    250K .......... .......... .......... .......... ..........  0% 14.4M 6s
    ##    300K .......... .......... .......... .......... ..........  0% 30.4M 6s
    ##    350K .......... .......... .......... .......... ..........  0% 22.3M 6s
    ##    400K .......... .......... .......... .......... ..........  0% 36.8M 5s
    ##    450K .......... .......... .......... .......... ..........  0% 22.8M 5s
    ##    500K .......... .......... .......... .......... ..........  0% 24.4M 5s
    ##    550K .......... .......... .......... .......... ..........  0% 20.9M 5s
    ##    600K .......... .......... .......... .......... ..........  0% 45.9M 4s
    ##    650K .......... .......... .......... .......... ..........  0% 18.1M 4s
    ##    700K .......... .......... .......... .......... ..........  0% 73.6M 4s
    ##    750K .......... .......... .......... .......... ..........  1% 15.9M 4s
    ##    800K .......... .......... .......... .......... ..........  1% 7.78M 5s
    ##    850K .......... .......... .......... .......... ..........  1%  101M 4s
    ##    900K .......... .......... .......... .......... ..........  1% 2.78M 6s
    ##    950K .......... .......... .......... .......... ..........  1% 90.3M 5s
    ##   1000K .......... .......... .......... .......... ..........  1%  105M 5s
    ##   1050K .......... .......... .......... .......... ..........  1% 14.0M 5s
    ##   1100K .......... .......... .......... .......... ..........  1% 96.8M 5s
    ##   1150K .......... .......... .......... .......... ..........  1% 13.5M 5s
    ##   1200K .......... .......... .......... .......... ..........  1% 87.7M 5s
    ##   1250K .......... .......... .......... .......... ..........  1% 16.0M 5s
    ##   1300K .......... .......... .......... .......... ..........  1% 74.7M 5s
    ##   1350K .......... .......... .......... .......... ..........  1% 21.3M 5s
    ##   1400K .......... .......... .......... .......... ..........  1% 65.9M 5s
    ##   1450K .......... .......... .......... .......... ..........  1% 16.8M 5s
    ##   1500K .......... .......... .......... .......... ..........  1% 15.9M 5s
    ##   1550K .......... .......... .......... .......... ..........  2% 50.4M 4s
    ##   1600K .......... .......... .......... .......... ..........  2% 87.3M 4s
    ##   1650K .......... .......... .......... .......... ..........  2% 22.7M 4s
    ##   1700K .......... .......... .......... .......... ..........  2% 50.8M 4s
    ##   1750K .......... .......... .......... .......... ..........  2% 16.2M 4s
    ##   1800K .......... .......... .......... .......... ..........  2% 87.1M 4s
    ##   1850K .......... .......... .......... .......... ..........  2% 20.3M 4s
    ##   1900K .......... .......... .......... .......... ..........  2% 47.8M 4s
    ##   1950K .......... .......... .......... .......... ..........  2% 18.7M 4s
    ##   2000K .......... .......... .......... .......... ..........  2% 47.1M 4s
    ##   2050K .......... .......... .......... .......... ..........  2% 58.2M 4s
    ##   2100K .......... .......... .......... .......... ..........  2% 19.0M 4s
    ##   2150K .......... .......... .......... .......... ..........  2% 63.9M 4s
    ##   2200K .......... .......... .......... .......... ..........  2% 22.6M 4s
    ##   2250K .......... .......... .......... .......... ..........  2% 54.5M 4s
    ##   2300K .......... .......... .......... .......... ..........  2% 67.5M 4s
    ##   2350K .......... .......... .......... .......... ..........  3% 15.5M 4s
    ##   2400K .......... .......... .......... .......... ..........  3% 80.8M 4s
    ##   2450K .......... .......... .......... .......... ..........  3% 78.9M 4s
    ##   2500K .......... .......... .......... .......... ..........  3% 17.1M 4s
    ##   2550K .......... .......... .......... .......... ..........  3% 76.5M 4s
    ##   2600K .......... .......... .......... .......... ..........  3% 14.5M 4s
    ##   2650K .......... .......... .......... .......... ..........  3% 74.8M 4s
    ##   2700K .......... .......... .......... .......... ..........  3% 96.1M 4s
    ##   2750K .......... .......... .......... .......... ..........  3% 14.3M 4s
    ##   2800K .......... .......... .......... .......... ..........  3% 63.5M 3s
    ##   2850K .......... .......... .......... .......... ..........  3% 8.91M 4s
    ##   2900K .......... .......... .......... .......... ..........  3% 77.8M 4s
    ##   2950K .......... .......... .......... .......... ..........  3% 18.2M 4s
    ##   3000K .......... .......... .......... .......... ..........  3% 53.0M 4s
    ##   3050K .......... .......... .......... .......... ..........  3% 98.9M 3s
    ##   3100K .......... .......... .......... .......... ..........  3% 19.7M 3s
    ##   3150K .......... .......... .......... .......... ..........  4% 58.0M 3s
    ##   3200K .......... .......... .......... .......... ..........  4% 94.2M 3s
    ##   3250K .......... .......... .......... .......... ..........  4% 21.1M 3s
    ##   3300K .......... .......... .......... .......... ..........  4% 64.5M 3s
    ##   3350K .......... .......... .......... .......... ..........  4% 21.0M 3s
    ##   3400K .......... .......... .......... .......... ..........  4% 31.5M 3s
    ##   3450K .......... .......... .......... .......... ..........  4% 93.5M 3s
    ##   3500K .......... .......... .......... .......... ..........  4% 25.2M 3s
    ##   3550K .......... .......... .......... .......... ..........  4% 39.0M 3s
    ##   3600K .......... .......... .......... .......... ..........  4% 14.5M 3s
    ##   3650K .......... .......... .......... .......... ..........  4% 87.1M 3s
    ##   3700K .......... .......... .......... .......... ..........  4% 87.5M 3s
    ##   3750K .......... .......... .......... .......... ..........  4% 10.8M 3s
    ##   3800K .......... .......... .......... .......... ..........  4% 87.5M 3s
    ##   3850K .......... .......... .......... .......... ..........  4% 18.7M 3s
    ##   3900K .......... .......... .......... .......... ..........  4% 36.2M 3s
    ##   3950K .......... .......... .......... .......... ..........  5% 90.3M 3s
    ##   4000K .......... .......... .......... .......... ..........  5% 22.2M 3s
    ##   4050K .......... .......... .......... .......... ..........  5% 53.9M 3s
    ##   4100K .......... .......... .......... .......... ..........  5% 81.6M 3s
    ##   4150K .......... .......... .......... .......... ..........  5% 17.6M 3s
    ##   4200K .......... .......... .......... .......... ..........  5% 74.1M 3s
    ##   4250K .......... .......... .......... .......... ..........  5% 43.2M 3s
    ##   4300K .......... .......... .......... .......... ..........  5% 32.7M 3s
    ##   4350K .......... .......... .......... .......... ..........  5% 45.0M 3s
    ##   4400K .......... .......... .......... .......... ..........  5% 37.6M 3s
    ##   4450K .......... .......... .......... .......... ..........  5% 32.5M 3s
    ##   4500K .......... .......... .......... .......... ..........  5% 27.3M 3s
    ##   4550K .......... .......... .......... .......... ..........  5% 32.6M 3s
    ##   4600K .......... .......... .......... .......... ..........  5% 28.7M 3s
    ##   4650K .......... .......... .......... .......... ..........  5% 55.2M 3s
    ##   4700K .......... .......... .......... .......... ..........  5% 37.1M 3s
    ##   4750K .......... .......... .......... .......... ..........  6% 20.7M 3s
    ##   4800K .......... .......... .......... .......... ..........  6% 55.0M 3s
    ##   4850K .......... .......... .......... .......... ..........  6% 58.3M 3s
    ##   4900K .......... .......... .......... .......... ..........  6% 68.2M 3s
    ##   4950K .......... .......... .......... .......... ..........  6% 11.1M 3s
    ##   5000K .......... .......... .......... .......... ..........  6% 82.2M 3s
    ##   5050K .......... .......... .......... .......... ..........  6% 12.8M 3s
    ##   5100K .......... .......... .......... .......... ..........  6% 58.5M 3s
    ##   5150K .......... .......... .......... .......... ..........  6% 89.0M 3s
    ##   5200K .......... .......... .......... .......... ..........  6% 18.4M 3s
    ##   5250K .......... .......... .......... .......... ..........  6% 56.4M 3s
    ##   5300K .......... .......... .......... .......... ..........  6% 89.4M 3s
    ##   5350K .......... .......... .......... .......... ..........  6% 14.0M 3s
    ##   5400K .......... .......... .......... .......... ..........  6% 93.5M 3s
    ##   5450K .......... .......... .......... .......... ..........  6% 6.72M 3s
    ##   5500K .......... .......... .......... .......... ..........  6% 85.6M 3s
    ##   5550K .......... .......... .......... .......... ..........  7%  113M 3s
    ##   5600K .......... .......... .......... .......... ..........  7% 18.2M 3s
    ##   5650K .......... .......... .......... .......... ..........  7% 93.7M 3s
    ##   5700K .......... .......... .......... .......... ..........  7% 64.4M 3s
    ##   5750K .......... .......... .......... .......... ..........  7% 29.1M 3s
    ##   5800K .......... .......... .......... .......... ..........  7% 79.5M 3s
    ##   5850K .......... .......... .......... .......... ..........  7% 7.63M 3s
    ##   5900K .......... .......... .......... .......... ..........  7% 56.7M 3s
    ##   5950K .......... .......... .......... .......... ..........  7%  103M 3s
    ##   6000K .......... .......... .......... .......... ..........  7% 97.8M 3s
    ##   6050K .......... .......... .......... .......... ..........  7% 34.7M 3s
    ##   6100K .......... .......... .......... .......... ..........  7% 55.9M 3s
    ##   6150K .......... .......... .......... .......... ..........  7% 97.4M 3s
    ##   6200K .......... .......... .......... .......... ..........  7% 28.2M 3s
    ##   6250K .......... .......... .......... .......... ..........  7% 65.8M 3s
    ##   6300K .......... .......... .......... .......... ..........  7% 47.1M 3s
    ##   6350K .......... .......... .......... .......... ..........  8% 94.4M 3s
    ##   6400K .......... .......... .......... .......... ..........  8% 8.03M 3s
    ##   6450K .......... .......... .......... .......... ..........  8% 52.9M 3s
    ##   6500K .......... .......... .......... .......... ..........  8% 84.0M 3s
    ##   6550K .......... .......... .......... .......... ..........  8%  106M 3s
    ##   6600K .......... .......... .......... .......... ..........  8% 8.12M 3s
    ##   6650K .......... .......... .......... .......... ..........  8% 64.3M 3s
    ##   6700K .......... .......... .......... .......... ..........  8% 85.0M 3s
    ##   6750K .......... .......... .......... .......... ..........  8%  111M 3s
    ##   6800K .......... .......... .......... .......... ..........  8% 92.0M 3s
    ##   6850K .......... .......... .......... .......... ..........  8% 19.9M 3s
    ##   6900K .......... .......... .......... .......... ..........  8% 54.9M 3s
    ##   6950K .......... .......... .......... .......... ..........  8% 91.8M 3s
    ##   7000K .......... .......... .......... .......... ..........  8% 91.9M 3s
    ##   7050K .......... .......... .......... .......... ..........  8% 31.0M 3s
    ##   7100K .......... .......... .......... .......... ..........  8% 60.9M 3s
    ##   7150K .......... .......... .......... .......... ..........  9% 63.9M 3s
    ##   7200K .......... .......... .......... .......... ..........  9%  103M 3s
    ##   7250K .......... .......... .......... .......... ..........  9% 32.2M 3s
    ##   7300K .......... .......... .......... .......... ..........  9% 25.6M 3s
    ##   7350K .......... .......... .......... .......... ..........  9%  102M 3s
    ##   7400K .......... .......... .......... .......... ..........  9% 62.4M 3s
    ##   7450K .......... .......... .......... .......... ..........  9% 92.1M 3s
    ##   7500K .......... .......... .......... .......... ..........  9% 19.2M 3s
    ##   7550K .......... .......... .......... .......... ..........  9% 63.5M 3s
    ##   7600K .......... .......... .......... .......... ..........  9% 67.0M 3s
    ##   7650K .......... .......... .......... .......... ..........  9%  103M 3s
    ##   7700K .......... .......... .......... .......... ..........  9% 18.1M 3s
    ##   7750K .......... .......... .......... .......... ..........  9% 86.1M 3s
    ##   7800K .......... .......... .......... .......... ..........  9%  104M 3s
    ##   7850K .......... .......... .......... .......... ..........  9% 28.8M 3s
    ##   7900K .......... .......... .......... .......... ..........  9% 48.0M 3s
    ##   7950K .......... .......... .......... .......... .......... 10% 44.2M 3s
    ##   8000K .......... .......... .......... .......... .......... 10% 96.3M 3s
    ##   8050K .......... .......... .......... .......... .......... 10% 82.1M 3s
    ##   8100K .......... .......... .......... .......... .......... 10% 37.0M 3s
    ##   8150K .......... .......... .......... .......... .......... 10% 56.6M 3s
    ##   8200K .......... .......... .......... .......... .......... 10% 50.8M 2s
    ##   8250K .......... .......... .......... .......... .......... 10% 95.3M 2s
    ##   8300K .......... .......... .......... .......... .......... 10% 31.6M 2s
    ##   8350K .......... .......... .......... .......... .......... 10% 55.5M 2s
    ##   8400K .......... .......... .......... .......... .......... 10% 59.4M 2s
    ##   8450K .......... .......... .......... .......... .......... 10% 96.4M 2s
    ##   8500K .......... .......... .......... .......... .......... 10% 67.0M 2s
    ##   8550K .......... .......... .......... .......... .......... 10% 5.57M 2s
    ##   8600K .......... .......... .......... .......... .......... 10% 94.1M 2s
    ##   8650K .......... .......... .......... .......... .......... 10%  131M 2s
    ##   8700K .......... .......... .......... .......... .......... 10% 71.9M 2s
    ##   8750K .......... .......... .......... .......... .......... 11% 24.1M 2s
    ##   8800K .......... .......... .......... .......... .......... 11% 45.5M 2s
    ##   8850K .......... .......... .......... .......... .......... 11%  117M 2s
    ##   8900K .......... .......... .......... .......... .......... 11% 95.5M 2s
    ##   8950K .......... .......... .......... .......... .......... 11% 28.4M 2s
    ##   9000K .......... .......... .......... .......... .......... 11% 4.64M 2s
    ##   9050K .......... .......... .......... .......... .......... 11% 43.6M 2s
    ##   9100K .......... .......... .......... .......... .......... 11% 89.1M 2s
    ##   9150K .......... .......... .......... .......... .......... 11%  125M 2s
    ##   9200K .......... .......... .......... .......... .......... 11% 84.4M 2s
    ##   9250K .......... .......... .......... .......... .......... 11% 62.6M 2s
    ##   9300K .......... .......... .......... .......... .......... 11% 22.9M 2s
    ##   9350K .......... .......... .......... .......... .......... 11% 45.0M 2s
    ##   9400K .......... .......... .......... .......... .......... 11% 44.1M 2s
    ##   9450K .......... .......... .......... .......... .......... 11% 86.6M 2s
    ##   9500K .......... .......... .......... .......... .......... 11% 46.9M 2s
    ##   9550K .......... .......... .......... .......... .......... 12% 65.5M 2s
    ##   9600K .......... .......... .......... .......... .......... 12% 77.8M 2s
    ##   9650K .......... .......... .......... .......... .......... 12% 14.4M 2s
    ##   9700K .......... .......... .......... .......... .......... 12% 73.2M 2s
    ##   9750K .......... .......... .......... .......... .......... 12% 94.9M 2s
    ##   9800K .......... .......... .......... .......... .......... 12% 73.5M 2s
    ##   9850K .......... .......... .......... .......... .......... 12% 17.8M 2s
    ##   9900K .......... .......... .......... .......... .......... 12% 52.8M 2s
    ##   9950K .......... .......... .......... .......... .......... 12% 81.9M 2s
    ##  10000K .......... .......... .......... .......... .......... 12% 69.9M 2s
    ##  10050K .......... .......... .......... .......... .......... 12%  106M 2s
    ##  10100K .......... .......... .......... .......... .......... 12% 13.4M 2s
    ##  10150K .......... .......... .......... .......... .......... 12%  104M 2s
    ##  10200K .......... .......... .......... .......... .......... 12% 81.8M 2s
    ##  10250K .......... .......... .......... .......... .......... 12% 93.8M 2s
    ##  10300K .......... .......... .......... .......... .......... 12% 16.7M 2s
    ##  10350K .......... .......... .......... .......... .......... 13% 79.3M 2s
    ##  10400K .......... .......... .......... .......... .......... 13% 81.7M 2s
    ##  10450K .......... .......... .......... .......... .......... 13% 47.9M 2s
    ##  10500K .......... .......... .......... .......... .......... 13% 59.5M 2s
    ##  10550K .......... .......... .......... .......... .......... 13% 50.8M 2s
    ##  10600K .......... .......... .......... .......... .......... 13% 75.2M 2s
    ##  10650K .......... .......... .......... .......... .......... 13% 24.0M 2s
    ##  10700K .......... .......... .......... .......... .......... 13% 39.6M 2s
    ##  10750K .......... .......... .......... .......... .......... 13% 96.4M 2s
    ##  10800K .......... .......... .......... .......... .......... 13% 81.4M 2s
    ##  10850K .......... .......... .......... .......... .......... 13% 7.58M 2s
    ##  10900K .......... .......... .......... .......... .......... 13% 42.2M 2s
    ##  10950K .......... .......... .......... .......... .......... 13% 95.1M 2s
    ##  11000K .......... .......... .......... .......... .......... 13% 68.5M 2s
    ##  11050K .......... .......... .......... .......... .......... 13% 22.2M 2s
    ##  11100K .......... .......... .......... .......... .......... 13% 43.1M 2s
    ##  11150K .......... .......... .......... .......... .......... 14% 90.6M 2s
    ##  11200K .......... .......... .......... .......... .......... 14% 47.0M 2s
    ##  11250K .......... .......... .......... .......... .......... 14% 44.2M 2s
    ##  11300K .......... .......... .......... .......... .......... 14% 41.5M 2s
    ##  11350K .......... .......... .......... .......... .......... 14% 88.6M 2s
    ##  11400K .......... .......... .......... .......... .......... 14% 67.7M 2s
    ##  11450K .......... .......... .......... .......... .......... 14% 47.3M 2s
    ##  11500K .......... .......... .......... .......... .......... 14% 83.6M 2s
    ##  11550K .......... .......... .......... .......... .......... 14% 7.17M 2s
    ##  11600K .......... .......... .......... .......... .......... 14% 57.6M 2s
    ##  11650K .......... .......... .......... .......... .......... 14% 38.0M 2s
    ##  11700K .......... .......... .......... .......... .......... 14% 74.5M 2s
    ##  11750K .......... .......... .......... .......... .......... 14%  108M 2s
    ##  11800K .......... .......... .......... .......... .......... 14% 50.8M 2s
    ##  11850K .......... .......... .......... .......... .......... 14% 54.2M 2s
    ##  11900K .......... .......... .......... .......... .......... 14% 42.7M 2s
    ##  11950K .......... .......... .......... .......... .......... 15% 58.4M 2s
    ##  12000K .......... .......... .......... .......... .......... 15% 59.7M 2s
    ##  12050K .......... .......... .......... .......... .......... 15% 62.4M 2s
    ##  12100K .......... .......... .......... .......... .......... 15% 82.1M 2s
    ##  12150K .......... .......... .......... .......... .......... 15% 72.5M 2s
    ##  12200K .......... .......... .......... .......... .......... 15% 43.3M 2s
    ##  12250K .......... .......... .......... .......... .......... 15% 97.0M 2s
    ##  12300K .......... .......... .......... .......... .......... 15% 42.3M 2s
    ##  12350K .......... .......... .......... .......... .......... 15% 92.2M 2s
    ##  12400K .......... .......... .......... .......... .......... 15% 55.8M 2s
    ##  12450K .......... .......... .......... .......... .......... 15% 48.5M 2s
    ##  12500K .......... .......... .......... .......... .......... 15% 81.9M 2s
    ##  12550K .......... .......... .......... .......... .......... 15% 12.0M 2s
    ##  12600K .......... .......... .......... .......... .......... 15% 36.5M 2s
    ##  12650K .......... .......... .......... .......... .......... 15% 42.7M 2s
    ##  12700K .......... .......... .......... .......... .......... 15% 81.1M 2s
    ##  12750K .......... .......... .......... .......... .......... 16%  130M 2s
    ##  12800K .......... .......... .......... .......... .......... 16% 91.0M 2s
    ##  12850K .......... .......... .......... .......... .......... 16% 26.5M 2s
    ##  12900K .......... .......... .......... .......... .......... 16% 65.3M 2s
    ##  12950K .......... .......... .......... .......... .......... 16% 39.1M 2s
    ##  13000K .......... .......... .......... .......... .......... 16% 93.5M 2s
    ##  13050K .......... .......... .......... .......... .......... 16%  114M 2s
    ##  13100K .......... .......... .......... .......... .......... 16% 49.6M 2s
    ##  13150K .......... .......... .......... .......... .......... 16% 64.7M 2s
    ##  13200K .......... .......... .......... .......... .......... 16% 42.6M 2s
    ##  13250K .......... .......... .......... .......... .......... 16% 94.4M 2s
    ##  13300K .......... .......... .......... .......... .......... 16% 12.9M 2s
    ##  13350K .......... .......... .......... .......... .......... 16% 39.8M 2s
    ##  13400K .......... .......... .......... .......... .......... 16% 83.1M 2s
    ##  13450K .......... .......... .......... .......... .......... 16%  105M 2s
    ##  13500K .......... .......... .......... .......... .......... 16% 83.8M 2s
    ##  13550K .......... .......... .......... .......... .......... 17% 88.5M 2s
    ##  13600K .......... .......... .......... .......... .......... 17% 60.8M 2s
    ##  13650K .......... .......... .......... .......... .......... 17% 30.8M 2s
    ##  13700K .......... .......... .......... .......... .......... 17% 50.3M 2s
    ##  13750K .......... .......... .......... .......... .......... 17% 59.8M 2s
    ##  13800K .......... .......... .......... .......... .......... 17% 81.6M 2s
    ##  13850K .......... .......... .......... .......... .......... 17% 92.3M 2s
    ##  13900K .......... .......... .......... .......... .......... 17% 32.5M 2s
    ##  13950K .......... .......... .......... .......... .......... 17% 54.6M 2s
    ##  14000K .......... .......... .......... .......... .......... 17% 51.1M 2s
    ##  14050K .......... .......... .......... .......... .......... 17%  106M 2s
    ##  14100K .......... .......... .......... .......... .......... 17% 43.4M 2s
    ##  14150K .......... .......... .......... .......... .......... 17% 70.9M 2s
    ##  14200K .......... .......... .......... .......... .......... 17% 62.2M 2s
    ##  14250K .......... .......... .......... .......... .......... 17% 57.3M 2s
    ##  14300K .......... .......... .......... .......... .......... 17% 77.6M 2s
    ##  14350K .......... .......... .......... .......... .......... 18% 60.2M 2s
    ##  14400K .......... .......... .......... .......... .......... 18% 62.4M 2s
    ##  14450K .......... .......... .......... .......... .......... 18% 65.2M 2s
    ##  14500K .......... .......... .......... .......... .......... 18% 62.8M 2s
    ##  14550K .......... .......... .......... .......... .......... 18%  102M 2s
    ##  14600K .......... .......... .......... .......... .......... 18% 72.6M 2s
    ##  14650K .......... .......... .......... .......... .......... 18% 93.4M 2s
    ##  14700K .......... .......... .......... .......... .......... 18% 76.2M 2s
    ##  14750K .......... .......... .......... .......... .......... 18% 49.9M 2s
    ##  14800K .......... .......... .......... .......... .......... 18% 57.1M 2s
    ##  14850K .......... .......... .......... .......... .......... 18%  109M 2s
    ##  14900K .......... .......... .......... .......... .......... 18% 69.7M 2s
    ##  14950K .......... .......... .......... .......... .......... 18% 50.7M 2s
    ##  15000K .......... .......... .......... .......... .......... 18% 65.4M 2s
    ##  15050K .......... .......... .......... .......... .......... 18%  107M 2s
    ##  15100K .......... .......... .......... .......... .......... 18% 23.9M 2s
    ##  15150K .......... .......... .......... .......... .......... 19%  105M 2s
    ##  15200K .......... .......... .......... .......... .......... 19% 50.1M 2s
    ##  15250K .......... .......... .......... .......... .......... 19% 98.9M 2s
    ##  15300K .......... .......... .......... .......... .......... 19% 77.8M 2s
    ##  15350K .......... .......... .......... .......... .......... 19% 46.5M 2s
    ##  15400K .......... .......... .......... .......... .......... 19% 61.9M 2s
    ##  15450K .......... .......... .......... .......... .......... 19% 44.8M 2s
    ##  15500K .......... .......... .......... .......... .......... 19% 73.1M 2s
    ##  15550K .......... .......... .......... .......... .......... 19% 82.0M 2s
    ##  15600K .......... .......... .......... .......... .......... 19% 98.2M 2s
    ##  15650K .......... .......... .......... .......... .......... 19% 27.0M 2s
    ##  15700K .......... .......... .......... .......... .......... 19% 42.0M 2s
    ##  15750K .......... .......... .......... .......... .......... 19% 78.7M 2s
    ##  15800K .......... .......... .......... .......... .......... 19% 94.0M 2s
    ##  15850K .......... .......... .......... .......... .......... 19% 12.4M 2s
    ##  15900K .......... .......... .......... .......... .......... 19% 71.8M 2s
    ##  15950K .......... .......... .......... .......... .......... 20% 95.7M 2s
    ##  16000K .......... .......... .......... .......... .......... 20% 92.0M 2s
    ##  16050K .......... .......... .......... .......... .......... 20%  100M 2s
    ##  16100K .......... .......... .......... .......... .......... 20% 64.4M 2s
    ##  16150K .......... .......... .......... .......... .......... 20% 32.8M 2s
    ##  16200K .......... .......... .......... .......... .......... 20% 37.3M 2s
    ##  16250K .......... .......... .......... .......... .......... 20%  105M 2s
    ##  16300K .......... .......... .......... .......... .......... 20% 85.5M 2s
    ##  16350K .......... .......... .......... .......... .......... 20% 96.5M 2s
    ##  16400K .......... .......... .......... .......... .......... 20% 18.3M 2s
    ##  16450K .......... .......... .......... .......... .......... 20% 22.7M 2s
    ##  16500K .......... .......... .......... .......... .......... 20% 60.1M 2s
    ##  16550K .......... .......... .......... .......... .......... 20% 67.2M 2s
    ##  16600K .......... .......... .......... .......... .......... 20% 90.6M 2s
    ##  16650K .......... .......... .......... .......... .......... 20% 70.6M 2s
    ##  16700K .......... .......... .......... .......... .......... 20% 74.3M 2s
    ##  16750K .......... .......... .......... .......... .......... 21% 50.4M 2s
    ##  16800K .......... .......... .......... .......... .......... 21% 45.1M 2s
    ##  16850K .......... .......... .......... .......... .......... 21%  104M 2s
    ##  16900K .......... .......... .......... .......... .......... 21% 71.1M 2s
    ##  16950K .......... .......... .......... .......... .......... 21%  101M 2s
    ##  17000K .......... .......... .......... .......... .......... 21% 51.6M 2s
    ##  17050K .......... .......... .......... .......... .......... 21% 55.4M 2s
    ##  17100K .......... .......... .......... .......... .......... 21% 75.9M 2s
    ##  17150K .......... .......... .......... .......... .......... 21% 60.7M 2s
    ##  17200K .......... .......... .......... .......... .......... 21%  102M 2s
    ##  17250K .......... .......... .......... .......... .......... 21% 55.2M 2s
    ##  17300K .......... .......... .......... .......... .......... 21% 49.6M 2s
    ##  17350K .......... .......... .......... .......... .......... 21% 54.9M 2s
    ##  17400K .......... .......... .......... .......... .......... 21% 97.7M 2s
    ##  17450K .......... .......... .......... .......... .......... 21% 98.0M 2s
    ##  17500K .......... .......... .......... .......... .......... 21% 73.0M 2s
    ##  17550K .......... .......... .......... .......... .......... 22% 26.7M 2s
    ##  17600K .......... .......... .......... .......... .......... 22% 61.3M 2s
    ##  17650K .......... .......... .......... .......... .......... 22%  104M 2s
    ##  17700K .......... .......... .......... .......... .......... 22%  116M 2s
    ##  17750K .......... .......... .......... .......... .......... 22% 35.2M 2s
    ##  17800K .......... .......... .......... .......... .......... 22% 52.7M 2s
    ##  17850K .......... .......... .......... .......... .......... 22% 93.4M 2s
    ##  17900K .......... .......... .......... .......... .......... 22% 31.4M 2s
    ##  17950K .......... .......... .......... .......... .......... 22%  116M 2s
    ##  18000K .......... .......... .......... .......... .......... 22% 82.2M 2s
    ##  18050K .......... .......... .......... .......... .......... 22% 60.8M 2s
    ##  18100K .......... .......... .......... .......... .......... 22% 74.4M 2s
    ##  18150K .......... .......... .......... .......... .......... 22% 81.2M 2s
    ##  18200K .......... .......... .......... .......... .......... 22% 56.3M 2s
    ##  18250K .......... .......... .......... .......... .......... 22%  106M 2s
    ##  18300K .......... .......... .......... .......... .......... 22% 55.5M 2s
    ##  18350K .......... .......... .......... .......... .......... 23% 57.8M 2s
    ##  18400K .......... .......... .......... .......... .......... 23% 57.2M 2s
    ##  18450K .......... .......... .......... .......... .......... 23% 80.1M 2s
    ##  18500K .......... .......... .......... .......... .......... 23% 80.4M 2s
    ##  18550K .......... .......... .......... .......... .......... 23% 16.3M 2s
    ##  18600K .......... .......... .......... .......... .......... 23% 65.1M 2s
    ##  18650K .......... .......... .......... .......... .......... 23% 39.3M 2s
    ##  18700K .......... .......... .......... .......... .......... 23% 87.8M 2s
    ##  18750K .......... .......... .......... .......... .......... 23% 66.5M 2s
    ##  18800K .......... .......... .......... .......... .......... 23% 74.2M 2s
    ##  18850K .......... .......... .......... .......... .......... 23% 64.5M 2s
    ##  18900K .......... .......... .......... .......... .......... 23% 66.8M 2s
    ##  18950K .......... .......... .......... .......... .......... 23% 50.9M 2s
    ##  19000K .......... .......... .......... .......... .......... 23% 53.9M 2s
    ##  19050K .......... .......... .......... .......... .......... 23%  109M 2s
    ##  19100K .......... .......... .......... .......... .......... 23%  112M 2s
    ##  19150K .......... .......... .......... .......... .......... 24% 54.3M 2s
    ##  19200K .......... .......... .......... .......... .......... 24% 62.9M 2s
    ##  19250K .......... .......... .......... .......... .......... 24% 60.3M 2s
    ##  19300K .......... .......... .......... .......... .......... 24% 67.2M 2s
    ##  19350K .......... .......... .......... .......... .......... 24%  110M 2s
    ##  19400K .......... .......... .......... .......... .......... 24% 67.8M 2s
    ##  19450K .......... .......... .......... .......... .......... 24% 37.6M 2s
    ##  19500K .......... .......... .......... .......... .......... 24% 70.3M 2s
    ##  19550K .......... .......... .......... .......... .......... 24% 78.4M 2s
    ##  19600K .......... .......... .......... .......... .......... 24% 87.3M 2s
    ##  19650K .......... .......... .......... .......... .......... 24% 77.8M 2s
    ##  19700K .......... .......... .......... .......... .......... 24%  101M 2s
    ##  19750K .......... .......... .......... .......... .......... 24% 36.7M 2s
    ##  19800K .......... .......... .......... .......... .......... 24% 51.3M 2s
    ##  19850K .......... .......... .......... .......... .......... 24% 67.5M 2s
    ##  19900K .......... .......... .......... .......... .......... 24% 96.1M 2s
    ##  19950K .......... .......... .......... .......... .......... 25% 48.5M 2s
    ##  20000K .......... .......... .......... .......... .......... 25% 58.6M 2s
    ##  20050K .......... .......... .......... .......... .......... 25% 63.9M 2s
    ##  20100K .......... .......... .......... .......... .......... 25% 61.5M 2s
    ##  20150K .......... .......... .......... .......... .......... 25% 93.1M 2s
    ##  20200K .......... .......... .......... .......... .......... 25% 50.5M 2s
    ##  20250K .......... .......... .......... .......... .......... 25% 64.5M 2s
    ##  20300K .......... .......... .......... .......... .......... 25% 87.5M 2s
    ##  20350K .......... .......... .......... .......... .......... 25% 59.3M 2s
    ##  20400K .......... .......... .......... .......... .......... 25% 35.8M 2s
    ##  20450K .......... .......... .......... .......... .......... 25% 59.9M 2s
    ##  20500K .......... .......... .......... .......... .......... 25% 96.5M 2s
    ##  20550K .......... .......... .......... .......... .......... 25% 93.2M 2s
    ##  20600K .......... .......... .......... .......... .......... 25% 62.0M 2s
    ##  20650K .......... .......... .......... .......... .......... 25% 77.9M 2s
    ##  20700K .......... .......... .......... .......... .......... 25% 46.7M 2s
    ##  20750K .......... .......... .......... .......... .......... 26% 98.3M 2s
    ##  20800K .......... .......... .......... .......... .......... 26%  103M 2s
    ##  20850K .......... .......... .......... .......... .......... 26% 46.2M 2s
    ##  20900K .......... .......... .......... .......... .......... 26% 41.5M 2s
    ##  20950K .......... .......... .......... .......... .......... 26% 61.5M 2s
    ##  21000K .......... .......... .......... .......... .......... 26% 87.1M 2s
    ##  21050K .......... .......... .......... .......... .......... 26% 80.8M 2s
    ##  21100K .......... .......... .......... .......... .......... 26% 61.0M 2s
    ##  21150K .......... .......... .......... .......... .......... 26% 51.1M 2s
    ##  21200K .......... .......... .......... .......... .......... 26% 35.5M 2s
    ##  21250K .......... .......... .......... .......... .......... 26% 81.0M 2s
    ##  21300K .......... .......... .......... .......... .......... 26% 56.1M 2s
    ##  21350K .......... .......... .......... .......... .......... 26% 68.6M 2s
    ##  21400K .......... .......... .......... .......... .......... 26% 51.2M 2s
    ##  21450K .......... .......... .......... .......... .......... 26% 76.6M 2s
    ##  21500K .......... .......... .......... .......... .......... 26% 53.2M 2s
    ##  21550K .......... .......... .......... .......... .......... 27% 74.6M 2s
    ##  21600K .......... .......... .......... .......... .......... 27% 86.7M 1s
    ##  21650K .......... .......... .......... .......... .......... 27%  107M 1s
    ##  21700K .......... .......... .......... .......... .......... 27%  113M 1s
    ##  21750K .......... .......... .......... .......... .......... 27% 60.7M 1s
    ##  21800K .......... .......... .......... .......... .......... 27%  102M 1s
    ##  21850K .......... .......... .......... .......... .......... 27%  126M 1s
    ##  21900K .......... .......... .......... .......... .......... 27% 50.1M 1s
    ##  21950K .......... .......... .......... .......... .......... 27%  119M 1s
    ##  22000K .......... .......... .......... .......... .......... 27% 32.6M 1s
    ##  22050K .......... .......... .......... .......... .......... 27% 21.8M 1s
    ##  22100K .......... .......... .......... .......... .......... 27% 70.3M 1s
    ##  22150K .......... .......... .......... .......... .......... 27% 41.7M 1s
    ##  22200K .......... .......... .......... .......... .......... 27% 64.6M 1s
    ##  22250K .......... .......... .......... .......... .......... 27%  111M 1s
    ##  22300K .......... .......... .......... .......... .......... 27%  104M 1s
    ##  22350K .......... .......... .......... .......... .......... 28% 75.2M 1s
    ##  22400K .......... .......... .......... .......... .......... 28% 58.2M 1s
    ##  22450K .......... .......... .......... .......... .......... 28% 9.05M 1s
    ##  22500K .......... .......... .......... .......... .......... 28% 61.7M 1s
    ##  22550K .......... .......... .......... .......... .......... 28% 49.8M 1s
    ##  22600K .......... .......... .......... .......... .......... 28% 60.6M 1s
    ##  22650K .......... .......... .......... .......... .......... 28% 82.8M 1s
    ##  22700K .......... .......... .......... .......... .......... 28% 41.4M 1s
    ##  22750K .......... .......... .......... .......... .......... 28% 39.8M 1s
    ##  22800K .......... .......... .......... .......... .......... 28% 78.8M 1s
    ##  22850K .......... .......... .......... .......... .......... 28% 81.3M 1s
    ##  22900K .......... .......... .......... .......... .......... 28% 78.8M 1s
    ##  22950K .......... .......... .......... .......... .......... 28% 53.7M 1s
    ##  23000K .......... .......... .......... .......... .......... 28% 66.2M 1s
    ##  23050K .......... .......... .......... .......... .......... 28%  108M 1s
    ##  23100K .......... .......... .......... .......... .......... 28% 78.5M 1s
    ##  23150K .......... .......... .......... .......... .......... 29% 85.1M 1s
    ##  23200K .......... .......... .......... .......... .......... 29%  106M 1s
    ##  23250K .......... .......... .......... .......... .......... 29% 97.2M 1s
    ##  23300K .......... .......... .......... .......... .......... 29% 16.4M 1s
    ##  23350K .......... .......... .......... .......... .......... 29% 71.1M 1s
    ##  23400K .......... .......... .......... .......... .......... 29% 44.2M 1s
    ##  23450K .......... .......... .......... .......... .......... 29% 79.3M 1s
    ##  23500K .......... .......... .......... .......... .......... 29%  105M 1s
    ##  23550K .......... .......... .......... .......... .......... 29% 90.9M 1s
    ##  23600K .......... .......... .......... .......... .......... 29% 81.5M 1s
    ##  23650K .......... .......... .......... .......... .......... 29% 43.7M 1s
    ##  23700K .......... .......... .......... .......... .......... 29% 54.0M 1s
    ##  23750K .......... .......... .......... .......... .......... 29% 53.9M 1s
    ##  23800K .......... .......... .......... .......... .......... 29% 80.2M 1s
    ##  23850K .......... .......... .......... .......... .......... 29% 83.1M 1s
    ##  23900K .......... .......... .......... .......... .......... 29% 59.3M 1s
    ##  23950K .......... .......... .......... .......... .......... 30% 89.4M 1s
    ##  24000K .......... .......... .......... .......... .......... 30% 51.6M 1s
    ##  24050K .......... .......... .......... .......... .......... 30%  107M 1s
    ##  24100K .......... .......... .......... .......... .......... 30% 14.5M 1s
    ##  24150K .......... .......... .......... .......... .......... 30% 58.8M 1s
    ##  24200K .......... .......... .......... .......... .......... 30% 50.9M 1s
    ##  24250K .......... .......... .......... .......... .......... 30%  129M 1s
    ##  24300K .......... .......... .......... .......... .......... 30%  101M 1s
    ##  24350K .......... .......... .......... .......... .......... 30% 96.5M 1s
    ##  24400K .......... .......... .......... .......... .......... 30%  119M 1s
    ##  24450K .......... .......... .......... .......... .......... 30% 23.4M 1s
    ##  24500K .......... .......... .......... .......... .......... 30%  130M 1s
    ##  24550K .......... .......... .......... .......... .......... 30% 80.0M 1s
    ##  24600K .......... .......... .......... .......... .......... 30%  114M 1s
    ##  24650K .......... .......... .......... .......... .......... 30%  139M 1s
    ##  24700K .......... .......... .......... .......... .......... 30% 32.1M 1s
    ##  24750K .......... .......... .......... .......... .......... 31% 55.3M 1s
    ##  24800K .......... .......... .......... .......... .......... 31% 72.9M 1s
    ##  24850K .......... .......... .......... .......... .......... 31%  122M 1s
    ##  24900K .......... .......... .......... .......... .......... 31%  101M 1s
    ##  24950K .......... .......... .......... .......... .......... 31% 95.8M 1s
    ##  25000K .......... .......... .......... .......... .......... 31% 65.0M 1s
    ##  25050K .......... .......... .......... .......... .......... 31% 86.7M 1s
    ##  25100K .......... .......... .......... .......... .......... 31%  127M 1s
    ##  25150K .......... .......... .......... .......... .......... 31% 36.7M 1s
    ##  25200K .......... .......... .......... .......... .......... 31% 10.7M 1s
    ##  25250K .......... .......... .......... .......... .......... 31% 38.7M 1s
    ##  25300K .......... .......... .......... .......... .......... 31% 69.4M 1s
    ##  25350K .......... .......... .......... .......... .......... 31% 53.4M 1s
    ##  25400K .......... .......... .......... .......... .......... 31% 75.2M 1s
    ##  25450K .......... .......... .......... .......... .......... 31%  116M 1s
    ##  25500K .......... .......... .......... .......... .......... 31% 57.7M 1s
    ##  25550K .......... .......... .......... .......... .......... 32% 86.3M 1s
    ##  25600K .......... .......... .......... .......... .......... 32% 80.1M 1s
    ##  25650K .......... .......... .......... .......... .......... 32% 47.4M 1s
    ##  25700K .......... .......... .......... .......... .......... 32% 93.8M 1s
    ##  25750K .......... .......... .......... .......... .......... 32% 81.2M 1s
    ##  25800K .......... .......... .......... .......... .......... 32% 54.2M 1s
    ##  25850K .......... .......... .......... .......... .......... 32% 27.6M 1s
    ##  25900K .......... .......... .......... .......... .......... 32% 98.7M 1s
    ##  25950K .......... .......... .......... .......... .......... 32% 81.6M 1s
    ##  26000K .......... .......... .......... .......... .......... 32%  111M 1s
    ##  26050K .......... .......... .......... .......... .......... 32%  146M 1s
    ##  26100K .......... .......... .......... .......... .......... 32% 42.4M 1s
    ##  26150K .......... .......... .......... .......... .......... 32% 38.4M 1s
    ##  26200K .......... .......... .......... .......... .......... 32% 55.0M 1s
    ##  26250K .......... .......... .......... .......... .......... 32%  117M 1s
    ##  26300K .......... .......... .......... .......... .......... 32%  103M 1s
    ##  26350K .......... .......... .......... .......... .......... 33% 16.3M 1s
    ##  26400K .......... .......... .......... .......... .......... 33% 71.3M 1s
    ##  26450K .......... .......... .......... .......... .......... 33% 44.4M 1s
    ##  26500K .......... .......... .......... .......... .......... 33% 96.6M 1s
    ##  26550K .......... .......... .......... .......... .......... 33%  124M 1s
    ##  26600K .......... .......... .......... .......... .......... 33% 39.7M 1s
    ##  26650K .......... .......... .......... .......... .......... 33%  139M 1s
    ##  26700K .......... .......... .......... .......... .......... 33% 45.9M 1s
    ##  26750K .......... .......... .......... .......... .......... 33%  107M 1s
    ##  26800K .......... .......... .......... .......... .......... 33% 54.9M 1s
    ##  26850K .......... .......... .......... .......... .......... 33% 68.1M 1s
    ##  26900K .......... .......... .......... .......... .......... 33%  111M 1s
    ##  26950K .......... .......... .......... .......... .......... 33% 77.5M 1s
    ##  27000K .......... .......... .......... .......... .......... 33% 36.1M 1s
    ##  27050K .......... .......... .......... .......... .......... 33% 64.6M 1s
    ##  27100K .......... .......... .......... .......... .......... 33% 60.1M 1s
    ##  27150K .......... .......... .......... .......... .......... 34%  136M 1s
    ##  27200K .......... .......... .......... .......... .......... 34% 80.3M 1s
    ##  27250K .......... .......... .......... .......... .......... 34% 69.7M 1s
    ##  27300K .......... .......... .......... .......... .......... 34% 44.8M 1s
    ##  27350K .......... .......... .......... .......... .......... 34% 85.3M 1s
    ##  27400K .......... .......... .......... .......... .......... 34% 97.4M 1s
    ##  27450K .......... .......... .......... .......... .......... 34%  145M 1s
    ##  27500K .......... .......... .......... .......... .......... 34% 16.7M 1s
    ##  27550K .......... .......... .......... .......... .......... 34% 51.7M 1s
    ##  27600K .......... .......... .......... .......... .......... 34% 54.9M 1s
    ##  27650K .......... .......... .......... .......... .......... 34% 72.0M 1s
    ##  27700K .......... .......... .......... .......... .......... 34%  108M 1s
    ##  27750K .......... .......... .......... .......... .......... 34%  146M 1s
    ##  27800K .......... .......... .......... .......... .......... 34% 22.0M 1s
    ##  27850K .......... .......... .......... .......... .......... 34% 98.9M 1s
    ##  27900K .......... .......... .......... .......... .......... 34% 55.6M 1s
    ##  27950K .......... .......... .......... .......... .......... 35% 81.6M 1s
    ##  28000K .......... .......... .......... .......... .......... 35% 88.3M 1s
    ##  28050K .......... .......... .......... .......... .......... 35% 96.7M 1s
    ##  28100K .......... .......... .......... .......... .......... 35% 46.9M 1s
    ##  28150K .......... .......... .......... .......... .......... 35% 77.6M 1s
    ##  28200K .......... .......... .......... .......... .......... 35% 85.8M 1s
    ##  28250K .......... .......... .......... .......... .......... 35% 56.2M 1s
    ##  28300K .......... .......... .......... .......... .......... 35% 28.2M 1s
    ##  28350K .......... .......... .......... .......... .......... 35% 99.2M 1s
    ##  28400K .......... .......... .......... .......... .......... 35% 72.3M 1s
    ##  28450K .......... .......... .......... .......... .......... 35%  105M 1s
    ##  28500K .......... .......... .......... .......... .......... 35% 82.2M 1s
    ##  28550K .......... .......... .......... .......... .......... 35% 48.0M 1s
    ##  28600K .......... .......... .......... .......... .......... 35% 77.5M 1s
    ##  28650K .......... .......... .......... .......... .......... 35% 80.1M 1s
    ##  28700K .......... .......... .......... .......... .......... 35% 79.7M 1s
    ##  28750K .......... .......... .......... .......... .......... 36% 69.8M 1s
    ##  28800K .......... .......... .......... .......... .......... 36% 52.5M 1s
    ##  28850K .......... .......... .......... .......... .......... 36% 66.0M 1s
    ##  28900K .......... .......... .......... .......... .......... 36% 72.0M 1s
    ##  28950K .......... .......... .......... .......... .......... 36% 87.1M 1s
    ##  29000K .......... .......... .......... .......... .......... 36% 65.5M 1s
    ##  29050K .......... .......... .......... .......... .......... 36% 61.5M 1s
    ##  29100K .......... .......... .......... .......... .......... 36% 83.2M 1s
    ##  29150K .......... .......... .......... .......... .......... 36% 32.8M 1s
    ##  29200K .......... .......... .......... .......... .......... 36% 88.1M 1s
    ##  29250K .......... .......... .......... .......... .......... 36%  109M 1s
    ##  29300K .......... .......... .......... .......... .......... 36% 80.5M 1s
    ##  29350K .......... .......... .......... .......... .......... 36%  103M 1s
    ##  29400K .......... .......... .......... .......... .......... 36% 32.4M 1s
    ##  29450K .......... .......... .......... .......... .......... 36%  101M 1s
    ##  29500K .......... .......... .......... .......... .......... 36% 72.2M 1s
    ##  29550K .......... .......... .......... .......... .......... 37%  103M 1s
    ##  29600K .......... .......... .......... .......... .......... 37% 90.7M 1s
    ##  29650K .......... .......... .......... .......... .......... 37% 32.5M 1s
    ##  29700K .......... .......... .......... .......... .......... 37% 70.5M 1s
    ##  29750K .......... .......... .......... .......... .......... 37% 48.4M 1s
    ##  29800K .......... .......... .......... .......... .......... 37% 73.0M 1s
    ##  29850K .......... .......... .......... .......... .......... 37%  110M 1s
    ##  29900K .......... .......... .......... .......... .......... 37% 38.7M 1s
    ##  29950K .......... .......... .......... .......... .......... 37%  101M 1s
    ##  30000K .......... .......... .......... .......... .......... 37% 69.9M 1s
    ##  30050K .......... .......... .......... .......... .......... 37% 76.5M 1s
    ##  30100K .......... .......... .......... .......... .......... 37%  101M 1s
    ##  30150K .......... .......... .......... .......... .......... 37% 64.6M 1s
    ##  30200K .......... .......... .......... .......... .......... 37% 67.6M 1s
    ##  30250K .......... .......... .......... .......... .......... 37% 38.7M 1s
    ##  30300K .......... .......... .......... .......... .......... 37% 70.7M 1s
    ##  30350K .......... .......... .......... .......... .......... 38%  119M 1s
    ##  30400K .......... .......... .......... .......... .......... 38% 47.0M 1s
    ##  30450K .......... .......... .......... .......... .......... 38% 80.8M 1s
    ##  30500K .......... .......... .......... .......... .......... 38% 77.3M 1s
    ##  30550K .......... .......... .......... .......... .......... 38% 95.8M 1s
    ##  30600K .......... .......... .......... .......... .......... 38% 87.6M 1s
    ##  30650K .......... .......... .......... .......... .......... 38% 40.1M 1s
    ##  30700K .......... .......... .......... .......... .......... 38% 89.0M 1s
    ##  30750K .......... .......... .......... .......... .......... 38%  108M 1s
    ##  30800K .......... .......... .......... .......... .......... 38% 63.8M 1s
    ##  30850K .......... .......... .......... .......... .......... 38% 83.8M 1s
    ##  30900K .......... .......... .......... .......... .......... 38% 63.1M 1s
    ##  30950K .......... .......... .......... .......... .......... 38% 27.6M 1s
    ##  31000K .......... .......... .......... .......... .......... 38% 69.3M 1s
    ##  31050K .......... .......... .......... .......... .......... 38%  115M 1s
    ##  31100K .......... .......... .......... .......... .......... 38% 54.5M 1s
    ##  31150K .......... .......... .......... .......... .......... 39% 97.0M 1s
    ##  31200K .......... .......... .......... .......... .......... 39% 75.2M 1s
    ##  31250K .......... .......... .......... .......... .......... 39% 98.7M 1s
    ##  31300K .......... .......... .......... .......... .......... 39% 76.8M 1s
    ##  31350K .......... .......... .......... .......... .......... 39% 41.4M 1s
    ##  31400K .......... .......... .......... .......... .......... 39% 79.2M 1s
    ##  31450K .......... .......... .......... .......... .......... 39%  103M 1s
    ##  31500K .......... .......... .......... .......... .......... 39% 91.2M 1s
    ##  31550K .......... .......... .......... .......... .......... 39%  112M 1s
    ##  31600K .......... .......... .......... .......... .......... 39% 21.0M 1s
    ##  31650K .......... .......... .......... .......... .......... 39%  106M 1s
    ##  31700K .......... .......... .......... .......... .......... 39% 77.2M 1s
    ##  31750K .......... .......... .......... .......... .......... 39% 82.9M 1s
    ##  31800K .......... .......... .......... .......... .......... 39% 96.5M 1s
    ##  31850K .......... .......... .......... .......... .......... 39%  100M 1s
    ##  31900K .......... .......... .......... .......... .......... 39% 53.9M 1s
    ##  31950K .......... .......... .......... .......... .......... 40% 74.4M 1s
    ##  32000K .......... .......... .......... .......... .......... 40% 91.6M 1s
    ##  32050K .......... .......... .......... .......... .......... 40% 87.4M 1s
    ##  32100K .......... .......... .......... .......... .......... 40% 54.8M 1s
    ##  32150K .......... .......... .......... .......... .......... 40% 30.5M 1s
    ##  32200K .......... .......... .......... .......... .......... 40% 89.6M 1s
    ##  32250K .......... .......... .......... .......... .......... 40%  112M 1s
    ##  32300K .......... .......... .......... .......... .......... 40% 86.9M 1s
    ##  32350K .......... .......... .......... .......... .......... 40%  108M 1s
    ##  32400K .......... .......... .......... .......... .......... 40% 93.6M 1s
    ##  32450K .......... .......... .......... .......... .......... 40% 31.8M 1s
    ##  32500K .......... .......... .......... .......... .......... 40% 62.8M 1s
    ##  32550K .......... .......... .......... .......... .......... 40% 89.8M 1s
    ##  32600K .......... .......... .......... .......... .......... 40% 70.8M 1s
    ##  32650K .......... .......... .......... .......... .......... 40% 63.1M 1s
    ##  32700K .......... .......... .......... .......... .......... 40% 65.7M 1s
    ##  32750K .......... .......... .......... .......... .......... 41% 80.0M 1s
    ##  32800K .......... .......... .......... .......... .......... 41% 98.6M 1s
    ##  32850K .......... .......... .......... .......... .......... 41% 99.6M 1s
    ##  32900K .......... .......... .......... .......... .......... 41% 72.4M 1s
    ##  32950K .......... .......... .......... .......... .......... 41% 78.8M 1s
    ##  33000K .......... .......... .......... .......... .......... 41% 87.8M 1s
    ##  33050K .......... .......... .......... .......... .......... 41% 22.8M 1s
    ##  33100K .......... .......... .......... .......... .......... 41% 11.9M 1s
    ##  33150K .......... .......... .......... .......... .......... 41%  116M 1s
    ##  33200K .......... .......... .......... .......... .......... 41%  111M 1s
    ##  33250K .......... .......... .......... .......... .......... 41%  117M 1s
    ##  33300K .......... .......... .......... .......... .......... 41% 24.7M 1s
    ##  33350K .......... .......... .......... .......... .......... 41% 69.7M 1s
    ##  33400K .......... .......... .......... .......... .......... 41%  113M 1s
    ##  33450K .......... .......... .......... .......... .......... 41% 71.6M 1s
    ##  33500K .......... .......... .......... .......... .......... 41% 38.6M 1s
    ##  33550K .......... .......... .......... .......... .......... 42%  129M 1s
    ##  33600K .......... .......... .......... .......... .......... 42% 98.2M 1s
    ##  33650K .......... .......... .......... .......... .......... 42% 34.8M 1s
    ##  33700K .......... .......... .......... .......... .......... 42% 78.6M 1s
    ##  33750K .......... .......... .......... .......... .......... 42% 93.9M 1s
    ##  33800K .......... .......... .......... .......... .......... 42% 92.8M 1s
    ##  33850K .......... .......... .......... .......... .......... 42%  115M 1s
    ##  33900K .......... .......... .......... .......... .......... 42%  113M 1s
    ##  33950K .......... .......... .......... .......... .......... 42% 30.5M 1s
    ##  34000K .......... .......... .......... .......... .......... 42% 76.7M 1s
    ##  34050K .......... .......... .......... .......... .......... 42%  111M 1s
    ##  34100K .......... .......... .......... .......... .......... 42% 63.4M 1s
    ##  34150K .......... .......... .......... .......... .......... 42%  110M 1s
    ##  34200K .......... .......... .......... .......... .......... 42%  106M 1s
    ##  34250K .......... .......... .......... .......... .......... 42% 17.3M 1s
    ##  34300K .......... .......... .......... .......... .......... 42% 93.2M 1s
    ##  34350K .......... .......... .......... .......... .......... 43% 84.2M 1s
    ##  34400K .......... .......... .......... .......... .......... 43% 91.8M 1s
    ##  34450K .......... .......... .......... .......... .......... 43%  101M 1s
    ##  34500K .......... .......... .......... .......... .......... 43% 70.7M 1s
    ##  34550K .......... .......... .......... .......... .......... 43% 78.7M 1s
    ##  34600K .......... .......... .......... .......... .......... 43% 38.3M 1s
    ##  34650K .......... .......... .......... .......... .......... 43% 79.3M 1s
    ##  34700K .......... .......... .......... .......... .......... 43% 78.0M 1s
    ##  34750K .......... .......... .......... .......... .......... 43% 94.5M 1s
    ##  34800K .......... .......... .......... .......... .......... 43% 63.3M 1s
    ##  34850K .......... .......... .......... .......... .......... 43% 94.7M 1s
    ##  34900K .......... .......... .......... .......... .......... 43% 53.2M 1s
    ##  34950K .......... .......... .......... .......... .......... 43% 83.9M 1s
    ##  35000K .......... .......... .......... .......... .......... 43% 72.1M 1s
    ##  35050K .......... .......... .......... .......... .......... 43% 81.5M 1s
    ##  35100K .......... .......... .......... .......... .......... 43% 76.6M 1s
    ##  35150K .......... .......... .......... .......... .......... 44% 75.4M 1s
    ##  35200K .......... .......... .......... .......... .......... 44% 69.1M 1s
    ##  35250K .......... .......... .......... .......... .......... 44% 80.1M 1s
    ##  35300K .......... .......... .......... .......... .......... 44% 64.6M 1s
    ##  35350K .......... .......... .......... .......... .......... 44% 93.2M 1s
    ##  35400K .......... .......... .......... .......... .......... 44% 61.6M 1s
    ##  35450K .......... .......... .......... .......... .......... 44% 48.4M 1s
    ##  35500K .......... .......... .......... .......... .......... 44% 71.9M 1s
    ##  35550K .......... .......... .......... .......... .......... 44%  109M 1s
    ##  35600K .......... .......... .......... .......... .......... 44% 67.5M 1s
    ##  35650K .......... .......... .......... .......... .......... 44% 88.8M 1s
    ##  35700K .......... .......... .......... .......... .......... 44% 75.9M 1s
    ##  35750K .......... .......... .......... .......... .......... 44% 79.7M 1s
    ##  35800K .......... .......... .......... .......... .......... 44% 78.5M 1s
    ##  35850K .......... .......... .......... .......... .......... 44% 58.5M 1s
    ##  35900K .......... .......... .......... .......... .......... 44% 76.9M 1s
    ##  35950K .......... .......... .......... .......... .......... 45% 48.9M 1s
    ##  36000K .......... .......... .......... .......... .......... 45% 72.5M 1s
    ##  36050K .......... .......... .......... .......... .......... 45%  104M 1s
    ##  36100K .......... .......... .......... .......... .......... 45% 8.40M 1s
    ##  36150K .......... .......... .......... .......... .......... 45% 44.9M 1s
    ##  36200K .......... .......... .......... .......... .......... 45% 79.5M 1s
    ##  36250K .......... .......... .......... .......... .......... 45% 71.7M 1s
    ##  36300K .......... .......... .......... .......... .......... 45% 93.8M 1s
    ##  36350K .......... .......... .......... .......... .......... 45%  108M 1s
    ##  36400K .......... .......... .......... .......... .......... 45% 85.4M 1s
    ##  36450K .......... .......... .......... .......... .......... 45% 80.2M 1s
    ##  36500K .......... .......... .......... .......... .......... 45% 60.9M 1s
    ##  36550K .......... .......... .......... .......... .......... 45% 28.7M 1s
    ##  36600K .......... .......... .......... .......... .......... 45% 88.2M 1s
    ##  36650K .......... .......... .......... .......... .......... 45% 81.7M 1s
    ##  36700K .......... .......... .......... .......... .......... 45% 90.9M 1s
    ##  36750K .......... .......... .......... .......... .......... 46% 84.2M 1s
    ##  36800K .......... .......... .......... .......... .......... 46% 91.4M 1s
    ##  36850K .......... .......... .......... .......... .......... 46% 85.0M 1s
    ##  36900K .......... .......... .......... .......... .......... 46% 89.9M 1s
    ##  36950K .......... .......... .......... .......... .......... 46% 61.1M 1s
    ##  37000K .......... .......... .......... .......... .......... 46% 63.6M 1s
    ##  37050K .......... .......... .......... .......... .......... 46% 93.5M 1s
    ##  37100K .......... .......... .......... .......... .......... 46% 99.8M 1s
    ##  37150K .......... .......... .......... .......... .......... 46% 89.5M 1s
    ##  37200K .......... .......... .......... .......... .......... 46% 97.5M 1s
    ##  37250K .......... .......... .......... .......... .......... 46% 34.6M 1s
    ##  37300K .......... .......... .......... .......... .......... 46% 64.6M 1s
    ##  37350K .......... .......... .......... .......... .......... 46% 43.1M 1s
    ##  37400K .......... .......... .......... .......... .......... 46% 81.4M 1s
    ##  37450K .......... .......... .......... .......... .......... 46% 89.6M 1s
    ##  37500K .......... .......... .......... .......... .......... 46% 85.1M 1s
    ##  37550K .......... .......... .......... .......... .......... 47% 92.9M 1s
    ##  37600K .......... .......... .......... .......... .......... 47% 92.9M 1s
    ##  37650K .......... .......... .......... .......... .......... 47%  120M 1s
    ##  37700K .......... .......... .......... .......... .......... 47% 40.2M 1s
    ##  37750K .......... .......... .......... .......... .......... 47% 85.7M 1s
    ##  37800K .......... .......... .......... .......... .......... 47% 83.0M 1s
    ##  37850K .......... .......... .......... .......... .......... 47%  126M 1s
    ##  37900K .......... .......... .......... .......... .......... 47% 60.2M 1s
    ##  37950K .......... .......... .......... .......... .......... 47% 93.9M 1s
    ##  38000K .......... .......... .......... .......... .......... 47% 56.1M 1s
    ##  38050K .......... .......... .......... .......... .......... 47%  106M 1s
    ##  38100K .......... .......... .......... .......... .......... 47% 78.6M 1s
    ##  38150K .......... .......... .......... .......... .......... 47%  115M 1s
    ##  38200K .......... .......... .......... .......... .......... 47% 48.7M 1s
    ##  38250K .......... .......... .......... .......... .......... 47% 81.4M 1s
    ##  38300K .......... .......... .......... .......... .......... 47% 94.2M 1s
    ##  38350K .......... .......... .......... .......... .......... 48%  113M 1s
    ##  38400K .......... .......... .......... .......... .......... 48% 59.9M 1s
    ##  38450K .......... .......... .......... .......... .......... 48% 79.1M 1s
    ##  38500K .......... .......... .......... .......... .......... 48% 96.6M 1s
    ##  38550K .......... .......... .......... .......... .......... 48%  100M 1s
    ##  38600K .......... .......... .......... .......... .......... 48% 85.2M 1s
    ##  38650K .......... .......... .......... .......... .......... 48% 45.5M 1s
    ##  38700K .......... .......... .......... .......... .......... 48% 75.3M 1s
    ##  38750K .......... .......... .......... .......... .......... 48%  132M 1s
    ##  38800K .......... .......... .......... .......... .......... 48% 95.9M 1s
    ##  38850K .......... .......... .......... .......... .......... 48%  117M 1s
    ##  38900K .......... .......... .......... .......... .......... 48% 71.6M 1s
    ##  38950K .......... .......... .......... .......... .......... 48% 46.0M 1s
    ##  39000K .......... .......... .......... .......... .......... 48% 89.0M 1s
    ##  39050K .......... .......... .......... .......... .......... 48%  101M 1s
    ##  39100K .......... .......... .......... .......... .......... 48% 92.3M 1s
    ##  39150K .......... .......... .......... .......... .......... 49%  111M 1s
    ##  39200K .......... .......... .......... .......... .......... 49%  110M 1s
    ##  39250K .......... .......... .......... .......... .......... 49% 63.9M 1s
    ##  39300K .......... .......... .......... .......... .......... 49% 31.3M 1s
    ##  39350K .......... .......... .......... .......... .......... 49% 85.9M 1s
    ##  39400K .......... .......... .......... .......... .......... 49% 81.4M 1s
    ##  39450K .......... .......... .......... .......... .......... 49% 84.6M 1s
    ##  39500K .......... .......... .......... .......... .......... 49% 83.6M 1s
    ##  39550K .......... .......... .......... .......... .......... 49%  106M 1s
    ##  39600K .......... .......... .......... .......... .......... 49%  116M 1s
    ##  39650K .......... .......... .......... .......... .......... 49% 49.5M 1s
    ##  39700K .......... .......... .......... .......... .......... 49% 27.4M 1s
    ##  39750K .......... .......... .......... .......... .......... 49% 91.2M 1s
    ##  39800K .......... .......... .......... .......... .......... 49% 50.1M 1s
    ##  39850K .......... .......... .......... .......... .......... 49% 83.7M 1s
    ##  39900K .......... .......... .......... .......... .......... 49% 91.5M 1s
    ##  39950K .......... .......... .......... .......... .......... 50%  140M 1s
    ##  40000K .......... .......... .......... .......... .......... 50%  121M 1s
    ##  40050K .......... .......... .......... .......... .......... 50% 97.5M 1s
    ##  40100K .......... .......... .......... .......... .......... 50% 15.9M 1s
    ##  40150K .......... .......... .......... .......... .......... 50% 95.2M 1s
    ##  40200K .......... .......... .......... .......... .......... 50%  104M 1s
    ##  40250K .......... .......... .......... .......... .......... 50%  102M 1s
    ##  40300K .......... .......... .......... .......... .......... 50%  124M 1s
    ##  40350K .......... .......... .......... .......... .......... 50%  129M 1s
    ##  40400K .......... .......... .......... .......... .......... 50% 65.5M 1s
    ##  40450K .......... .......... .......... .......... .......... 50% 41.4M 1s
    ##  40500K .......... .......... .......... .......... .......... 50% 76.4M 1s
    ##  40550K .......... .......... .......... .......... .......... 50% 86.9M 1s
    ##  40600K .......... .......... .......... .......... .......... 50% 58.9M 1s
    ##  40650K .......... .......... .......... .......... .......... 50% 90.9M 1s
    ##  40700K .......... .......... .......... .......... .......... 50% 85.8M 1s
    ##  40750K .......... .......... .......... .......... .......... 51% 7.75M 1s
    ##  40800K .......... .......... .......... .......... .......... 51% 70.5M 1s
    ##  40850K .......... .......... .......... .......... .......... 51% 91.5M 1s
    ##  40900K .......... .......... .......... .......... .......... 51% 62.5M 1s
    ##  40950K .......... .......... .......... .......... .......... 51% 89.8M 1s
    ##  41000K .......... .......... .......... .......... .......... 51% 77.1M 1s
    ##  41050K .......... .......... .......... .......... .......... 51% 95.8M 1s
    ##  41100K .......... .......... .......... .......... .......... 51% 7.93M 1s
    ##  41150K .......... .......... .......... .......... .......... 51% 57.1M 1s
    ##  41200K .......... .......... .......... .......... .......... 51% 55.1M 1s
    ##  41250K .......... .......... .......... .......... .......... 51% 71.4M 1s
    ##  41300K .......... .......... .......... .......... .......... 51% 83.7M 1s
    ##  41350K .......... .......... .......... .......... .......... 51% 77.1M 1s
    ##  41400K .......... .......... .......... .......... .......... 51%  103M 1s
    ##  41450K .......... .......... .......... .......... .......... 51% 91.3M 1s
    ##  41500K .......... .......... .......... .......... .......... 51% 85.0M 1s
    ##  41550K .......... .......... .......... .......... .......... 52% 42.2M 1s
    ##  41600K .......... .......... .......... .......... .......... 52% 61.9M 1s
    ##  41650K .......... .......... .......... .......... .......... 52% 80.7M 1s
    ##  41700K .......... .......... .......... .......... .......... 52%  104M 1s
    ##  41750K .......... .......... .......... .......... .......... 52% 99.3M 1s
    ##  41800K .......... .......... .......... .......... .......... 52% 73.2M 1s
    ##  41850K .......... .......... .......... .......... .......... 52% 96.1M 1s
    ##  41900K .......... .......... .......... .......... .......... 52% 43.9M 1s
    ##  41950K .......... .......... .......... .......... .......... 52% 86.5M 1s
    ##  42000K .......... .......... .......... .......... .......... 52% 73.3M 1s
    ##  42050K .......... .......... .......... .......... .......... 52% 81.7M 1s
    ##  42100K .......... .......... .......... .......... .......... 52% 62.6M 1s
    ##  42150K .......... .......... .......... .......... .......... 52% 77.7M 1s
    ##  42200K .......... .......... .......... .......... .......... 52% 81.7M 1s
    ##  42250K .......... .......... .......... .......... .......... 52% 55.8M 1s
    ##  42300K .......... .......... .......... .......... .......... 52% 88.4M 1s
    ##  42350K .......... .......... .......... .......... .......... 53% 93.6M 1s
    ##  42400K .......... .......... .......... .......... .......... 53% 49.3M 1s
    ##  42450K .......... .......... .......... .......... .......... 53% 97.6M 1s
    ##  42500K .......... .......... .......... .......... .......... 53% 34.2M 1s
    ##  42550K .......... .......... .......... .......... .......... 53% 86.8M 1s
    ##  42600K .......... .......... .......... .......... .......... 53% 91.7M 1s
    ##  42650K .......... .......... .......... .......... .......... 53% 92.9M 1s
    ##  42700K .......... .......... .......... .......... .......... 53% 55.6M 1s
    ##  42750K .......... .......... .......... .......... .......... 53% 73.2M 1s
    ##  42800K .......... .......... .......... .......... .......... 53% 99.5M 1s
    ##  42850K .......... .......... .......... .......... .......... 53%  118M 1s
    ##  42900K .......... .......... .......... .......... .......... 53% 92.7M 1s
    ##  42950K .......... .......... .......... .......... .......... 53% 81.6M 1s
    ##  43000K .......... .......... .......... .......... .......... 53% 39.8M 1s
    ##  43050K .......... .......... .......... .......... .......... 53% 85.2M 1s
    ##  43100K .......... .......... .......... .......... .......... 53% 75.4M 1s
    ##  43150K .......... .......... .......... .......... .......... 54% 91.9M 1s
    ##  43200K .......... .......... .......... .......... .......... 54% 75.3M 1s
    ##  43250K .......... .......... .......... .......... .......... 54%  120M 1s
    ##  43300K .......... .......... .......... .......... .......... 54% 98.4M 1s
    ##  43350K .......... .......... .......... .......... .......... 54% 95.9M 1s
    ##  43400K .......... .......... .......... .......... .......... 54% 84.5M 1s
    ##  43450K .......... .......... .......... .......... .......... 54% 38.6M 1s
    ##  43500K .......... .......... .......... .......... .......... 54% 90.3M 1s
    ##  43550K .......... .......... .......... .......... .......... 54% 88.2M 1s
    ##  43600K .......... .......... .......... .......... .......... 54% 67.6M 1s
    ##  43650K .......... .......... .......... .......... .......... 54% 59.3M 1s
    ##  43700K .......... .......... .......... .......... .......... 54% 99.3M 1s
    ##  43750K .......... .......... .......... .......... .......... 54%  120M 1s
    ##  43800K .......... .......... .......... .......... .......... 54% 86.1M 1s
    ##  43850K .......... .......... .......... .......... .......... 54% 35.8M 1s
    ##  43900K .......... .......... .......... .......... .......... 54% 67.7M 1s
    ##  43950K .......... .......... .......... .......... .......... 55% 71.8M 1s
    ##  44000K .......... .......... .......... .......... .......... 55% 72.2M 1s
    ##  44050K .......... .......... .......... .......... .......... 55%  102M 1s
    ##  44100K .......... .......... .......... .......... .......... 55% 74.8M 1s
    ##  44150K .......... .......... .......... .......... .......... 55% 19.8M 1s
    ##  44200K .......... .......... .......... .......... .......... 55% 90.5M 1s
    ##  44250K .......... .......... .......... .......... .......... 55%  117M 1s
    ##  44300K .......... .......... .......... .......... .......... 55%  109M 1s
    ##  44350K .......... .......... .......... .......... .......... 55%  101M 1s
    ##  44400K .......... .......... .......... .......... .......... 55%  115M 1s
    ##  44450K .......... .......... .......... .......... .......... 55%  108M 1s
    ##  44500K .......... .......... .......... .......... .......... 55% 25.6M 1s
    ##  44550K .......... .......... .......... .......... .......... 55%  108M 1s
    ##  44600K .......... .......... .......... .......... .......... 55% 82.9M 1s
    ##  44650K .......... .......... .......... .......... .......... 55% 72.6M 1s
    ##  44700K .......... .......... .......... .......... .......... 55%  108M 1s
    ##  44750K .......... .......... .......... .......... .......... 56%  119M 1s
    ##  44800K .......... .......... .......... .......... .......... 56% 36.5M 1s
    ##  44850K .......... .......... .......... .......... .......... 56% 88.3M 1s
    ##  44900K .......... .......... .......... .......... .......... 56%  101M 1s
    ##  44950K .......... .......... .......... .......... .......... 56% 97.1M 1s
    ##  45000K .......... .......... .......... .......... .......... 56% 52.2M 1s
    ##  45050K .......... .......... .......... .......... .......... 56%  128M 1s
    ##  45100K .......... .......... .......... .......... .......... 56% 80.9M 1s
    ##  45150K .......... .......... .......... .......... .......... 56% 95.1M 1s
    ##  45200K .......... .......... .......... .......... .......... 56% 95.9M 1s
    ##  45250K .......... .......... .......... .......... .......... 56%  114M 1s
    ##  45300K .......... .......... .......... .......... .......... 56% 29.9M 1s
    ##  45350K .......... .......... .......... .......... .......... 56% 55.9M 1s
    ##  45400K .......... .......... .......... .......... .......... 56% 52.0M 1s
    ##  45450K .......... .......... .......... .......... .......... 56% 99.3M 1s
    ##  45500K .......... .......... .......... .......... .......... 56% 39.2M 1s
    ##  45550K .......... .......... .......... .......... .......... 57% 45.2M 1s
    ##  45600K .......... .......... .......... .......... .......... 57% 59.9M 1s
    ##  45650K .......... .......... .......... .......... .......... 57% 82.5M 1s
    ##  45700K .......... .......... .......... .......... .......... 57% 25.7M 1s
    ##  45750K .......... .......... .......... .......... .......... 57%  115M 1s
    ##  45800K .......... .......... .......... .......... .......... 57% 85.2M 1s
    ##  45850K .......... .......... .......... .......... .......... 57%  119M 1s
    ##  45900K .......... .......... .......... .......... .......... 57% 12.5M 1s
    ##  45950K .......... .......... .......... .......... .......... 57%  104M 1s
    ##  46000K .......... .......... .......... .......... .......... 57%  104M 1s
    ##  46050K .......... .......... .......... .......... .......... 57%  127M 1s
    ##  46100K .......... .......... .......... .......... .......... 57%  110M 1s
    ##  46150K .......... .......... .......... .......... .......... 57%  114M 1s
    ##  46200K .......... .......... .......... .......... .......... 57% 20.6M 1s
    ##  46250K .......... .......... .......... .......... .......... 57%  105M 1s
    ##  46300K .......... .......... .......... .......... .......... 57%  114M 1s
    ##  46350K .......... .......... .......... .......... .......... 58% 68.4M 1s
    ##  46400K .......... .......... .......... .......... .......... 58%  101M 1s
    ##  46450K .......... .......... .......... .......... .......... 58%  127M 1s
    ##  46500K .......... .......... .......... .......... .......... 58% 53.4M 1s
    ##  46550K .......... .......... .......... .......... .......... 58% 79.3M 1s
    ##  46600K .......... .......... .......... .......... .......... 58% 84.8M 1s
    ##  46650K .......... .......... .......... .......... .......... 58% 89.0M 1s
    ##  46700K .......... .......... .......... .......... .......... 58% 54.0M 1s
    ##  46750K .......... .......... .......... .......... .......... 58%  101M 1s
    ##  46800K .......... .......... .......... .......... .......... 58% 91.7M 1s
    ##  46850K .......... .......... .......... .......... .......... 58%  116M 1s
    ##  46900K .......... .......... .......... .......... .......... 58% 24.2M 1s
    ##  46950K .......... .......... .......... .......... .......... 58% 75.2M 1s
    ##  47000K .......... .......... .......... .......... .......... 58% 68.2M 1s
    ##  47050K .......... .......... .......... .......... .......... 58% 79.6M 1s
    ##  47100K .......... .......... .......... .......... .......... 58% 83.2M 1s
    ##  47150K .......... .......... .......... .......... .......... 59% 86.0M 1s
    ##  47200K .......... .......... .......... .......... .......... 59% 87.3M 1s
    ##  47250K .......... .......... .......... .......... .......... 59% 20.5M 1s
    ##  47300K .......... .......... .......... .......... .......... 59% 69.3M 1s
    ##  47350K .......... .......... .......... .......... .......... 59% 42.0M 1s
    ##  47400K .......... .......... .......... .......... .......... 59% 49.5M 1s
    ##  47450K .......... .......... .......... .......... .......... 59% 94.8M 1s
    ##  47500K .......... .......... .......... .......... .......... 59% 83.1M 1s
    ##  47550K .......... .......... .......... .......... .......... 59%  114M 1s
    ##  47600K .......... .......... .......... .......... .......... 59%  106M 1s
    ##  47650K .......... .......... .......... .......... .......... 59% 52.7M 1s
    ##  47700K .......... .......... .......... .......... .......... 59% 79.1M 1s
    ##  47750K .......... .......... .......... .......... .......... 59% 97.1M 1s
    ##  47800K .......... .......... .......... .......... .......... 59% 36.3M 1s
    ##  47850K .......... .......... .......... .......... .......... 59%  108M 1s
    ##  47900K .......... .......... .......... .......... .......... 59% 98.2M 1s
    ##  47950K .......... .......... .......... .......... .......... 60% 9.81M 1s
    ##  48000K .......... .......... .......... .......... .......... 60% 67.3M 1s
    ##  48050K .......... .......... .......... .......... .......... 60% 84.2M 1s
    ##  48100K .......... .......... .......... .......... .......... 60% 26.6M 1s
    ##  48150K .......... .......... .......... .......... .......... 60% 89.2M 1s
    ##  48200K .......... .......... .......... .......... .......... 60% 81.9M 1s
    ##  48250K .......... .......... .......... .......... .......... 60% 90.4M 1s
    ##  48300K .......... .......... .......... .......... .......... 60%  103M 1s
    ##  48350K .......... .......... .......... .......... .......... 60% 25.0M 1s
    ##  48400K .......... .......... .......... .......... .......... 60% 72.2M 1s
    ##  48450K .......... .......... .......... .......... .......... 60% 54.0M 1s
    ##  48500K .......... .......... .......... .......... .......... 60% 67.2M 1s
    ##  48550K .......... .......... .......... .......... .......... 60% 86.8M 1s
    ##  48600K .......... .......... .......... .......... .......... 60% 84.1M 1s
    ##  48650K .......... .......... .......... .......... .......... 60% 98.7M 1s
    ##  48700K .......... .......... .......... .......... .......... 60% 84.6M 1s
    ##  48750K .......... .......... .......... .......... .......... 61% 18.1M 1s
    ##  48800K .......... .......... .......... .......... .......... 61% 64.7M 1s
    ##  48850K .......... .......... .......... .......... .......... 61% 96.4M 1s
    ##  48900K .......... .......... .......... .......... .......... 61%  101M 1s
    ##  48950K .......... .......... .......... .......... .......... 61%  113M 1s
    ##  49000K .......... .......... .......... .......... .......... 61% 93.6M 1s
    ##  49050K .......... .......... .......... .......... .......... 61% 11.5M 1s
    ##  49100K .......... .......... .......... .......... .......... 61% 48.6M 1s
    ##  49150K .......... .......... .......... .......... .......... 61% 69.9M 1s
    ##  49200K .......... .......... .......... .......... .......... 61% 59.9M 1s
    ##  49250K .......... .......... .......... .......... .......... 61% 60.8M 1s
    ##  49300K .......... .......... .......... .......... .......... 61% 26.8M 1s
    ##  49350K .......... .......... .......... .......... .......... 61%  100M 1s
    ##  49400K .......... .......... .......... .......... .......... 61% 91.9M 1s
    ##  49450K .......... .......... .......... .......... .......... 61% 44.7M 1s
    ##  49500K .......... .......... .......... .......... .......... 61% 78.2M 1s
    ##  49550K .......... .......... .......... .......... .......... 62% 82.0M 1s
    ##  49600K .......... .......... .......... .......... .......... 62% 74.8M 1s
    ##  49650K .......... .......... .......... .......... .......... 62%  108M 1s
    ##  49700K .......... .......... .......... .......... .......... 62%  101M 1s
    ##  49750K .......... .......... .......... .......... .......... 62% 4.80M 1s
    ##  49800K .......... .......... .......... .......... .......... 62%  105M 1s
    ##  49850K .......... .......... .......... .......... .......... 62% 58.2M 1s
    ##  49900K .......... .......... .......... .......... .......... 62% 61.0M 1s
    ##  49950K .......... .......... .......... .......... .......... 62%  105M 1s
    ##  50000K .......... .......... .......... .......... .......... 62%  103M 1s
    ##  50050K .......... .......... .......... .......... .......... 62% 75.6M 1s
    ##  50100K .......... .......... .......... .......... .......... 62% 86.5M 1s
    ##  50150K .......... .......... .......... .......... .......... 62% 88.8M 1s
    ##  50200K .......... .......... .......... .......... .......... 62% 50.7M 1s
    ##  50250K .......... .......... .......... .......... .......... 62% 97.7M 1s
    ##  50300K .......... .......... .......... .......... .......... 62% 67.4M 1s
    ##  50350K .......... .......... .......... .......... .......... 63% 80.1M 1s
    ##  50400K .......... .......... .......... .......... .......... 63% 84.2M 1s
    ##  50450K .......... .......... .......... .......... .......... 63%  100M 1s
    ##  50500K .......... .......... .......... .......... .......... 63% 54.1M 1s
    ##  50550K .......... .......... .......... .......... .......... 63% 24.5M 1s
    ##  50600K .......... .......... .......... .......... .......... 63% 73.2M 1s
    ##  50650K .......... .......... .......... .......... .......... 63% 79.9M 1s
    ##  50700K .......... .......... .......... .......... .......... 63% 76.4M 1s
    ##  50750K .......... .......... .......... .......... .......... 63%  122M 1s
    ##  50800K .......... .......... .......... .......... .......... 63% 5.74M 1s
    ##  50850K .......... .......... .......... .......... .......... 63% 74.4M 1s
    ##  50900K .......... .......... .......... .......... .......... 63% 90.2M 1s
    ##  50950K .......... .......... .......... .......... .......... 63% 76.6M 1s
    ##  51000K .......... .......... .......... .......... .......... 63% 60.2M 1s
    ##  51050K .......... .......... .......... .......... .......... 63%  103M 1s
    ##  51100K .......... .......... .......... .......... .......... 63% 87.5M 1s
    ##  51150K .......... .......... .......... .......... .......... 64% 47.0M 1s
    ##  51200K .......... .......... .......... .......... .......... 64% 74.8M 1s
    ##  51250K .......... .......... .......... .......... .......... 64%  120M 1s
    ##  51300K .......... .......... .......... .......... .......... 64% 86.9M 1s
    ##  51350K .......... .......... .......... .......... .......... 64% 68.8M 1s
    ##  51400K .......... .......... .......... .......... .......... 64% 89.4M 1s
    ##  51450K .......... .......... .......... .......... .......... 64% 9.53M 1s
    ##  51500K .......... .......... .......... .......... .......... 64% 25.0M 1s
    ##  51550K .......... .......... .......... .......... .......... 64%  107M 1s
    ##  51600K .......... .......... .......... .......... .......... 64% 54.6M 1s
    ##  51650K .......... .......... .......... .......... .......... 64%  108M 1s
    ##  51700K .......... .......... .......... .......... .......... 64% 87.7M 1s
    ##  51750K .......... .......... .......... .......... .......... 64%  104M 1s
    ##  51800K .......... .......... .......... .......... .......... 64% 99.6M 1s
    ##  51850K .......... .......... .......... .......... .......... 64% 34.7M 1s
    ##  51900K .......... .......... .......... .......... .......... 65% 58.0M 1s
    ##  51950K .......... .......... .......... .......... .......... 65% 45.4M 1s
    ##  52000K .......... .......... .......... .......... .......... 65% 73.8M 1s
    ##  52050K .......... .......... .......... .......... .......... 65% 40.7M 1s
    ##  52100K .......... .......... .......... .......... .......... 65%  102M 1s
    ##  52150K .......... .......... .......... .......... .......... 65% 97.5M 1s
    ##  52200K .......... .......... .......... .......... .......... 65% 95.8M 1s
    ##  52250K .......... .......... .......... .......... .......... 65%  124M 1s
    ##  52300K .......... .......... .......... .......... .......... 65% 66.9M 1s
    ##  52350K .......... .......... .......... .......... .......... 65% 96.9M 1s
    ##  52400K .......... .......... .......... .......... .......... 65% 85.2M 1s
    ##  52450K .......... .......... .......... .......... .......... 65%  102M 1s
    ##  52500K .......... .......... .......... .......... .......... 65% 87.8M 1s
    ##  52550K .......... .......... .......... .......... .......... 65% 87.4M 1s
    ##  52600K .......... .......... .......... .......... .......... 65% 77.6M 1s
    ##  52650K .......... .......... .......... .......... .......... 65% 6.69M 1s
    ##  52700K .......... .......... .......... .......... .......... 66% 66.1M 1s
    ##  52750K .......... .......... .......... .......... .......... 66% 45.7M 1s
    ##  52800K .......... .......... .......... .......... .......... 66% 74.0M 1s
    ##  52850K .......... .......... .......... .......... .......... 66% 92.0M 1s
    ##  52900K .......... .......... .......... .......... .......... 66% 91.4M 1s
    ##  52950K .......... .......... .......... .......... .......... 66% 97.0M 1s
    ##  53000K .......... .......... .......... .......... .......... 66%  105M 1s
    ##  53050K .......... .......... .......... .......... .......... 66% 6.78M 1s
    ##  53100K .......... .......... .......... .......... .......... 66% 43.4M 1s
    ##  53150K .......... .......... .......... .......... .......... 66% 58.3M 1s
    ##  53200K .......... .......... .......... .......... .......... 66% 68.3M 1s
    ##  53250K .......... .......... .......... .......... .......... 66% 98.0M 1s
    ##  53300K .......... .......... .......... .......... .......... 66% 87.0M 1s
    ##  53350K .......... .......... .......... .......... .......... 66%  100M 1s
    ##  53400K .......... .......... .......... .......... .......... 66% 82.2M 1s
    ##  53450K .......... .......... .......... .......... .......... 66%  115M 1s
    ##  53500K .......... .......... .......... .......... .......... 67% 76.4M 1s
    ##  53550K .......... .......... .......... .......... .......... 67% 86.8M 1s
    ##  53600K .......... .......... .......... .......... .......... 67% 47.5M 1s
    ##  53650K .......... .......... .......... .......... .......... 67% 85.0M 1s
    ##  53700K .......... .......... .......... .......... .......... 67% 80.8M 1s
    ##  53750K .......... .......... .......... .......... .......... 67%  105M 1s
    ##  53800K .......... .......... .......... .......... .......... 67% 78.2M 1s
    ##  53850K .......... .......... .......... .......... .......... 67% 84.5M 1s
    ##  53900K .......... .......... .......... .......... .......... 67% 99.3M 1s
    ##  53950K .......... .......... .......... .......... .......... 67% 33.3M 1s
    ##  54000K .......... .......... .......... .......... .......... 67% 28.2M 1s
    ##  54050K .......... .......... .......... .......... .......... 67%  109M 1s
    ##  54100K .......... .......... .......... .......... .......... 67%  102M 1s
    ##  54150K .......... .......... .......... .......... .......... 67%  113M 1s
    ##  54200K .......... .......... .......... .......... .......... 67%  103M 1s
    ##  54250K .......... .......... .......... .......... .......... 67% 14.8M 1s
    ##  54300K .......... .......... .......... .......... .......... 68% 93.2M 1s
    ##  54350K .......... .......... .......... .......... .......... 68% 33.7M 1s
    ##  54400K .......... .......... .......... .......... .......... 68% 69.0M 1s
    ##  54450K .......... .......... .......... .......... .......... 68%  101M 1s
    ##  54500K .......... .......... .......... .......... .......... 68% 96.2M 1s
    ##  54550K .......... .......... .......... .......... .......... 68%  115M 1s
    ##  54600K .......... .......... .......... .......... .......... 68%  101M 1s
    ##  54650K .......... .......... .......... .......... .......... 68% 18.2M 1s
    ##  54700K .......... .......... .......... .......... .......... 68% 65.1M 1s
    ##  54750K .......... .......... .......... .......... .......... 68% 87.6M 1s
    ##  54800K .......... .......... .......... .......... .......... 68% 60.1M 1s
    ##  54850K .......... .......... .......... .......... .......... 68% 94.9M 1s
    ##  54900K .......... .......... .......... .......... .......... 68% 83.8M 1s
    ##  54950K .......... .......... .......... .......... .......... 68%  111M 1s
    ##  55000K .......... .......... .......... .......... .......... 68% 82.6M 1s
    ##  55050K .......... .......... .......... .......... .......... 68% 87.4M 1s
    ##  55100K .......... .......... .......... .......... .......... 69% 89.6M 1s
    ##  55150K .......... .......... .......... .......... .......... 69% 85.1M 1s
    ##  55200K .......... .......... .......... .......... .......... 69% 69.4M 1s
    ##  55250K .......... .......... .......... .......... .......... 69% 77.3M 1s
    ##  55300K .......... .......... .......... .......... .......... 69% 91.5M 1s
    ##  55350K .......... .......... .......... .......... .......... 69% 59.0M 1s
    ##  55400K .......... .......... .......... .......... .......... 69%  101M 0s
    ##  55450K .......... .......... .......... .......... .......... 69%  109M 0s
    ##  55500K .......... .......... .......... .......... .......... 69% 88.7M 0s
    ##  55550K .......... .......... .......... .......... .......... 69% 98.1M 0s
    ##  55600K .......... .......... .......... .......... .......... 69% 6.33M 0s
    ##  55650K .......... .......... .......... .......... .......... 69% 92.4M 0s
    ##  55700K .......... .......... .......... .......... .......... 69% 94.4M 0s
    ##  55750K .......... .......... .......... .......... .......... 69% 79.1M 0s
    ##  55800K .......... .......... .......... .......... .......... 69% 95.5M 0s
    ##  55850K .......... .......... .......... .......... .......... 69%  106M 0s
    ##  55900K .......... .......... .......... .......... .......... 70%  111M 0s
    ##  55950K .......... .......... .......... .......... .......... 70% 84.8M 0s
    ##  56000K .......... .......... .......... .......... .......... 70%  114M 0s
    ##  56050K .......... .......... .......... .......... .......... 70% 16.7M 0s
    ##  56100K .......... .......... .......... .......... .......... 70% 77.9M 0s
    ##  56150K .......... .......... .......... .......... .......... 70%  107M 0s
    ##  56200K .......... .......... .......... .......... .......... 70%  105M 0s
    ##  56250K .......... .......... .......... .......... .......... 70%  123M 0s
    ##  56300K .......... .......... .......... .......... .......... 70%  110M 0s
    ##  56350K .......... .......... .......... .......... .......... 70% 36.3M 0s
    ##  56400K .......... .......... .......... .......... .......... 70%  104M 0s
    ##  56450K .......... .......... .......... .......... .......... 70% 86.8M 0s
    ##  56500K .......... .......... .......... .......... .......... 70% 92.8M 0s
    ##  56550K .......... .......... .......... .......... .......... 70% 84.6M 0s
    ##  56600K .......... .......... .......... .......... .......... 70% 96.3M 0s
    ##  56650K .......... .......... .......... .......... .......... 70%  112M 0s
    ##  56700K .......... .......... .......... .......... .......... 71% 12.5M 0s
    ##  56750K .......... .......... .......... .......... .......... 71% 78.7M 0s
    ##  56800K .......... .......... .......... .......... .......... 71% 34.7M 0s
    ##  56850K .......... .......... .......... .......... .......... 71% 65.0M 0s
    ##  56900K .......... .......... .......... .......... .......... 71% 80.8M 0s
    ##  56950K .......... .......... .......... .......... .......... 71% 99.1M 0s
    ##  57000K .......... .......... .......... .......... .......... 71% 87.6M 0s
    ##  57050K .......... .......... .......... .......... .......... 71%  121M 0s
    ##  57100K .......... .......... .......... .......... .......... 71% 98.2M 0s
    ##  57150K .......... .......... .......... .......... .......... 71% 42.2M 0s
    ##  57200K .......... .......... .......... .......... .......... 71%  100M 0s
    ##  57250K .......... .......... .......... .......... .......... 71% 63.7M 0s
    ##  57300K .......... .......... .......... .......... .......... 71% 77.7M 0s
    ##  57350K .......... .......... .......... .......... .......... 71% 87.0M 0s
    ##  57400K .......... .......... .......... .......... .......... 71% 76.6M 0s
    ##  57450K .......... .......... .......... .......... .......... 71% 94.4M 0s
    ##  57500K .......... .......... .......... .......... .......... 72%  117M 0s
    ##  57550K .......... .......... .......... .......... .......... 72% 32.7M 0s
    ##  57600K .......... .......... .......... .......... .......... 72% 69.2M 0s
    ##  57650K .......... .......... .......... .......... .......... 72% 41.5M 0s
    ##  57700K .......... .......... .......... .......... .......... 72% 85.6M 0s
    ##  57750K .......... .......... .......... .......... .......... 72% 89.9M 0s
    ##  57800K .......... .......... .......... .......... .......... 72% 90.6M 0s
    ##  57850K .......... .......... .......... .......... .......... 72%  120M 0s
    ##  57900K .......... .......... .......... .......... .......... 72% 93.2M 0s
    ##  57950K .......... .......... .......... .......... .......... 72%  110M 0s
    ##  58000K .......... .......... .......... .......... .......... 72% 91.2M 0s
    ##  58050K .......... .......... .......... .......... .......... 72% 60.0M 0s
    ##  58100K .......... .......... .......... .......... .......... 72% 66.3M 0s
    ##  58150K .......... .......... .......... .......... .......... 72% 91.8M 0s
    ##  58200K .......... .......... .......... .......... .......... 72% 90.6M 0s
    ##  58250K .......... .......... .......... .......... .......... 72% 93.3M 0s
    ##  58300K .......... .......... .......... .......... .......... 73% 85.0M 0s
    ##  58350K .......... .......... .......... .......... .......... 73%  111M 0s
    ##  58400K .......... .......... .......... .......... .......... 73% 93.0M 0s
    ##  58450K .......... .......... .......... .......... .......... 73%  111M 0s
    ##  58500K .......... .......... .......... .......... .......... 73% 79.7M 0s
    ##  58550K .......... .......... .......... .......... .......... 73% 80.3M 0s
    ##  58600K .......... .......... .......... .......... .......... 73% 92.6M 0s
    ##  58650K .......... .......... .......... .......... .......... 73%  116M 0s
    ##  58700K .......... .......... .......... .......... .......... 73% 60.2M 0s
    ##  58750K .......... .......... .......... .......... .......... 73% 56.9M 0s
    ##  58800K .......... .......... .......... .......... .......... 73% 77.3M 0s
    ##  58850K .......... .......... .......... .......... .......... 73% 92.9M 0s
    ##  58900K .......... .......... .......... .......... .......... 73%  117M 0s
    ##  58950K .......... .......... .......... .......... .......... 73% 88.0M 0s
    ##  59000K .......... .......... .......... .......... .......... 73% 23.6M 0s
    ##  59050K .......... .......... .......... .......... .......... 73%  117M 0s
    ##  59100K .......... .......... .......... .......... .......... 74% 70.5M 0s
    ##  59150K .......... .......... .......... .......... .......... 74% 82.5M 0s
    ##  59200K .......... .......... .......... .......... .......... 74%  105M 0s
    ##  59250K .......... .......... .......... .......... .......... 74% 10.2M 0s
    ##  59300K .......... .......... .......... .......... .......... 74% 79.2M 0s
    ##  59350K .......... .......... .......... .......... .......... 74%  111M 0s
    ##  59400K .......... .......... .......... .......... .......... 74% 55.0M 0s
    ##  59450K .......... .......... .......... .......... .......... 74% 63.3M 0s
    ##  59500K .......... .......... .......... .......... .......... 74% 76.8M 0s
    ##  59550K .......... .......... .......... .......... .......... 74%  115M 0s
    ##  59600K .......... .......... .......... .......... .......... 74% 86.0M 0s
    ##  59650K .......... .......... .......... .......... .......... 74%  134M 0s
    ##  59700K .......... .......... .......... .......... .......... 74% 26.4M 0s
    ##  59750K .......... .......... .......... .......... .......... 74% 98.0M 0s
    ##  59800K .......... .......... .......... .......... .......... 74% 46.7M 0s
    ##  59850K .......... .......... .......... .......... .......... 74%  110M 0s
    ##  59900K .......... .......... .......... .......... .......... 75% 64.5M 0s
    ##  59950K .......... .......... .......... .......... .......... 75%  112M 0s
    ##  60000K .......... .......... .......... .......... .......... 75% 98.5M 0s
    ##  60050K .......... .......... .......... .......... .......... 75%  120M 0s
    ##  60100K .......... .......... .......... .......... .......... 75%  111M 0s
    ##  60150K .......... .......... .......... .......... .......... 75%  117M 0s
    ##  60200K .......... .......... .......... .......... .......... 75% 58.6M 0s
    ##  60250K .......... .......... .......... .......... .......... 75% 63.3M 0s
    ##  60300K .......... .......... .......... .......... .......... 75%  108M 0s
    ##  60350K .......... .......... .......... .......... .......... 75%  117M 0s
    ##  60400K .......... .......... .......... .......... .......... 75% 88.9M 0s
    ##  60450K .......... .......... .......... .......... .......... 75%  105M 0s
    ##  60500K .......... .......... .......... .......... .......... 75% 76.7M 0s
    ##  60550K .......... .......... .......... .......... .......... 75% 98.9M 0s
    ##  60600K .......... .......... .......... .......... .......... 75%  106M 0s
    ##  60650K .......... .......... .......... .......... .......... 75% 23.8M 0s
    ##  60700K .......... .......... .......... .......... .......... 76%  113M 0s
    ##  60750K .......... .......... .......... .......... .......... 76%  135M 0s
    ##  60800K .......... .......... .......... .......... .......... 76% 82.2M 0s
    ##  60850K .......... .......... .......... .......... .......... 76%  144M 0s
    ##  60900K .......... .......... .......... .......... .......... 76% 14.4M 0s
    ##  60950K .......... .......... .......... .......... .......... 76% 87.7M 0s
    ##  61000K .......... .......... .......... .......... .......... 76% 27.2M 0s
    ##  61050K .......... .......... .......... .......... .......... 76% 88.4M 0s
    ##  61100K .......... .......... .......... .......... .......... 76% 55.5M 0s
    ##  61150K .......... .......... .......... .......... .......... 76% 81.8M 0s
    ##  61200K .......... .......... .......... .......... .......... 76% 89.5M 0s
    ##  61250K .......... .......... .......... .......... .......... 76% 98.3M 0s
    ##  61300K .......... .......... .......... .......... .......... 76% 77.1M 0s
    ##  61350K .......... .......... .......... .......... .......... 76% 88.9M 0s
    ##  61400K .......... .......... .......... .......... .......... 76% 57.7M 0s
    ##  61450K .......... .......... .......... .......... .......... 76% 91.3M 0s
    ##  61500K .......... .......... .......... .......... .......... 77% 70.5M 0s
    ##  61550K .......... .......... .......... .......... .......... 77% 80.9M 0s
    ##  61600K .......... .......... .......... .......... .......... 77% 71.2M 0s
    ##  61650K .......... .......... .......... .......... .......... 77% 62.5M 0s
    ##  61700K .......... .......... .......... .......... .......... 77% 69.9M 0s
    ##  61750K .......... .......... .......... .......... .......... 77%  103M 0s
    ##  61800K .......... .......... .......... .......... .......... 77% 74.1M 0s
    ##  61850K .......... .......... .......... .......... .......... 77% 89.3M 0s
    ##  61900K .......... .......... .......... .......... .......... 77% 70.9M 0s
    ##  61950K .......... .......... .......... .......... .......... 77% 48.5M 0s
    ##  62000K .......... .......... .......... .......... .......... 77% 77.4M 0s
    ##  62050K .......... .......... .......... .......... .......... 77% 87.1M 0s
    ##  62100K .......... .......... .......... .......... .......... 77% 78.5M 0s
    ##  62150K .......... .......... .......... .......... .......... 77% 67.5M 0s
    ##  62200K .......... .......... .......... .......... .......... 77% 79.0M 0s
    ##  62250K .......... .......... .......... .......... .......... 77% 87.2M 0s
    ##  62300K .......... .......... .......... .......... .......... 78% 88.0M 0s
    ##  62350K .......... .......... .......... .......... .......... 78% 97.2M 0s
    ##  62400K .......... .......... .......... .......... .......... 78% 91.5M 0s
    ##  62450K .......... .......... .......... .......... .......... 78% 81.3M 0s
    ##  62500K .......... .......... .......... .......... .......... 78% 55.5M 0s
    ##  62550K .......... .......... .......... .......... .......... 78% 85.0M 0s
    ##  62600K .......... .......... .......... .......... .......... 78% 71.7M 0s
    ##  62650K .......... .......... .......... .......... .......... 78% 75.2M 0s
    ##  62700K .......... .......... .......... .......... .......... 78% 9.65M 0s
    ##  62750K .......... .......... .......... .......... .......... 78% 86.8M 0s
    ##  62800K .......... .......... .......... .......... .......... 78% 82.5M 0s
    ##  62850K .......... .......... .......... .......... .......... 78% 91.3M 0s
    ##  62900K .......... .......... .......... .......... .......... 78% 87.7M 0s
    ##  62950K .......... .......... .......... .......... .......... 78% 91.8M 0s
    ##  63000K .......... .......... .......... .......... .......... 78% 86.6M 0s
    ##  63050K .......... .......... .......... .......... .......... 78% 82.8M 0s
    ##  63100K .......... .......... .......... .......... .......... 79% 99.2M 0s
    ##  63150K .......... .......... .......... .......... .......... 79%  103M 0s
    ##  63200K .......... .......... .......... .......... .......... 79%  101M 0s
    ##  63250K .......... .......... .......... .......... .......... 79% 28.0M 0s
    ##  63300K .......... .......... .......... .......... .......... 79% 54.6M 0s
    ##  63350K .......... .......... .......... .......... .......... 79% 98.8M 0s
    ##  63400K .......... .......... .......... .......... .......... 79%  104M 0s
    ##  63450K .......... .......... .......... .......... .......... 79%  117M 0s
    ##  63500K .......... .......... .......... .......... .......... 79%  107M 0s
    ##  63550K .......... .......... .......... .......... .......... 79%  111M 0s
    ##  63600K .......... .......... .......... .......... .......... 79% 20.7M 0s
    ##  63650K .......... .......... .......... .......... .......... 79% 63.2M 0s
    ##  63700K .......... .......... .......... .......... .......... 79% 88.7M 0s
    ##  63750K .......... .......... .......... .......... .......... 79% 84.4M 0s
    ##  63800K .......... .......... .......... .......... .......... 79% 94.2M 0s
    ##  63850K .......... .......... .......... .......... .......... 79% 82.1M 0s
    ##  63900K .......... .......... .......... .......... .......... 80% 74.7M 0s
    ##  63950K .......... .......... .......... .......... .......... 80% 84.1M 0s
    ##  64000K .......... .......... .......... .......... .......... 80% 82.9M 0s
    ##  64050K .......... .......... .......... .......... .......... 80% 85.9M 0s
    ##  64100K .......... .......... .......... .......... .......... 80% 73.2M 0s
    ##  64150K .......... .......... .......... .......... .......... 80% 55.0M 0s
    ##  64200K .......... .......... .......... .......... .......... 80% 72.4M 0s
    ##  64250K .......... .......... .......... .......... .......... 80% 77.5M 0s
    ##  64300K .......... .......... .......... .......... .......... 80% 88.2M 0s
    ##  64350K .......... .......... .......... .......... .......... 80% 84.7M 0s
    ##  64400K .......... .......... .......... .......... .......... 80% 74.0M 0s
    ##  64450K .......... .......... .......... .......... .......... 80% 91.9M 0s
    ##  64500K .......... .......... .......... .......... .......... 80% 77.7M 0s
    ##  64550K .......... .......... .......... .......... .......... 80% 28.7M 0s
    ##  64600K .......... .......... .......... .......... .......... 80% 68.5M 0s
    ##  64650K .......... .......... .......... .......... .......... 80% 91.5M 0s
    ##  64700K .......... .......... .......... .......... .......... 81%  104M 0s
    ##  64750K .......... .......... .......... .......... .......... 81% 87.3M 0s
    ##  64800K .......... .......... .......... .......... .......... 81% 69.2M 0s
    ##  64850K .......... .......... .......... .......... .......... 81% 91.0M 0s
    ##  64900K .......... .......... .......... .......... .......... 81% 82.4M 0s
    ##  64950K .......... .......... .......... .......... .......... 81% 86.7M 0s
    ##  65000K .......... .......... .......... .......... .......... 81% 97.4M 0s
    ##  65050K .......... .......... .......... .......... .......... 81%  102M 0s
    ##  65100K .......... .......... .......... .......... .......... 81% 58.0M 0s
    ##  65150K .......... .......... .......... .......... .......... 81% 86.1M 0s
    ##  65200K .......... .......... .......... .......... .......... 81% 80.5M 0s
    ##  65250K .......... .......... .......... .......... .......... 81% 72.1M 0s
    ##  65300K .......... .......... .......... .......... .......... 81% 83.1M 0s
    ##  65350K .......... .......... .......... .......... .......... 81% 87.3M 0s
    ##  65400K .......... .......... .......... .......... .......... 81% 76.8M 0s
    ##  65450K .......... .......... .......... .......... .......... 81% 98.3M 0s
    ##  65500K .......... .......... .......... .......... .......... 82% 78.1M 0s
    ##  65550K .......... .......... .......... .......... .......... 82% 68.3M 0s
    ##  65600K .......... .......... .......... .......... .......... 82% 95.7M 0s
    ##  65650K .......... .......... .......... .......... .......... 82% 98.2M 0s
    ##  65700K .......... .......... .......... .......... .......... 82% 89.3M 0s
    ##  65750K .......... .......... .......... .......... .......... 82% 39.0M 0s
    ##  65800K .......... .......... .......... .......... .......... 82% 51.0M 0s
    ##  65850K .......... .......... .......... .......... .......... 82% 83.7M 0s
    ##  65900K .......... .......... .......... .......... .......... 82%  102M 0s
    ##  65950K .......... .......... .......... .......... .......... 82% 70.0M 0s
    ##  66000K .......... .......... .......... .......... .......... 82% 95.6M 0s
    ##  66050K .......... .......... .......... .......... .......... 82% 91.0M 0s
    ##  66100K .......... .......... .......... .......... .......... 82% 62.9M 0s
    ##  66150K .......... .......... .......... .......... .......... 82% 93.3M 0s
    ##  66200K .......... .......... .......... .......... .......... 82%  102M 0s
    ##  66250K .......... .......... .......... .......... .......... 82% 26.3M 0s
    ##  66300K .......... .......... .......... .......... .......... 83% 89.3M 0s
    ##  66350K .......... .......... .......... .......... .......... 83%  101M 0s
    ##  66400K .......... .......... .......... .......... .......... 83% 98.5M 0s
    ##  66450K .......... .......... .......... .......... .......... 83% 73.4M 0s
    ##  66500K .......... .......... .......... .......... .......... 83%  107M 0s
    ##  66550K .......... .......... .......... .......... .......... 83% 92.4M 0s
    ##  66600K .......... .......... .......... .......... .......... 83% 83.9M 0s
    ##  66650K .......... .......... .......... .......... .......... 83% 71.7M 0s
    ##  66700K .......... .......... .......... .......... .......... 83% 58.9M 0s
    ##  66750K .......... .......... .......... .......... .......... 83% 57.3M 0s
    ##  66800K .......... .......... .......... .......... .......... 83% 84.6M 0s
    ##  66850K .......... .......... .......... .......... .......... 83% 83.5M 0s
    ##  66900K .......... .......... .......... .......... .......... 83% 90.5M 0s
    ##  66950K .......... .......... .......... .......... .......... 83% 87.0M 0s
    ##  67000K .......... .......... .......... .......... .......... 83% 90.7M 0s
    ##  67050K .......... .......... .......... .......... .......... 83% 46.0M 0s
    ##  67100K .......... .......... .......... .......... .......... 84% 90.8M 0s
    ##  67150K .......... .......... .......... .......... .......... 84% 16.9M 0s
    ##  67200K .......... .......... .......... .......... .......... 84% 87.8M 0s
    ##  67250K .......... .......... .......... .......... .......... 84% 94.8M 0s
    ##  67300K .......... .......... .......... .......... .......... 84%  103M 0s
    ##  67350K .......... .......... .......... .......... .......... 84% 98.5M 0s
    ##  67400K .......... .......... .......... .......... .......... 84% 77.3M 0s
    ##  67450K .......... .......... .......... .......... .......... 84%  103M 0s
    ##  67500K .......... .......... .......... .......... .......... 84% 17.6M 0s
    ##  67550K .......... .......... .......... .......... .......... 84% 55.9M 0s
    ##  67600K .......... .......... .......... .......... .......... 84% 90.5M 0s
    ##  67650K .......... .......... .......... .......... .......... 84% 97.8M 0s
    ##  67700K .......... .......... .......... .......... .......... 84% 84.4M 0s
    ##  67750K .......... .......... .......... .......... .......... 84%  109M 0s
    ##  67800K .......... .......... .......... .......... .......... 84% 89.7M 0s
    ##  67850K .......... .......... .......... .......... .......... 84% 67.9M 0s
    ##  67900K .......... .......... .......... .......... .......... 85% 86.3M 0s
    ##  67950K .......... .......... .......... .......... .......... 85% 65.0M 0s
    ##  68000K .......... .......... .......... .......... .......... 85% 69.1M 0s
    ##  68050K .......... .......... .......... .......... .......... 85% 92.3M 0s
    ##  68100K .......... .......... .......... .......... .......... 85% 94.4M 0s
    ##  68150K .......... .......... .......... .......... .......... 85% 85.5M 0s
    ##  68200K .......... .......... .......... .......... .......... 85% 92.6M 0s
    ##  68250K .......... .......... .......... .......... .......... 85% 75.0M 0s
    ##  68300K .......... .......... .......... .......... .......... 85% 75.0M 0s
    ##  68350K .......... .......... .......... .......... .......... 85% 93.2M 0s
    ##  68400K .......... .......... .......... .......... .......... 85% 79.7M 0s
    ##  68450K .......... .......... .......... .......... .......... 85% 58.1M 0s
    ##  68500K .......... .......... .......... .......... .......... 85% 83.5M 0s
    ##  68550K .......... .......... .......... .......... .......... 85% 92.8M 0s
    ##  68600K .......... .......... .......... .......... .......... 85% 91.1M 0s
    ##  68650K .......... .......... .......... .......... .......... 85% 98.0M 0s
    ##  68700K .......... .......... .......... .......... .......... 86% 77.8M 0s
    ##  68750K .......... .......... .......... .......... .......... 86% 83.8M 0s
    ##  68800K .......... .......... .......... .......... .......... 86% 80.8M 0s
    ##  68850K .......... .......... .......... .......... .......... 86% 94.3M 0s
    ##  68900K .......... .......... .......... .......... .......... 86% 84.6M 0s
    ##  68950K .......... .......... .......... .......... .......... 86% 34.1M 0s
    ##  69000K .......... .......... .......... .......... .......... 86% 69.4M 0s
    ##  69050K .......... .......... .......... .......... .......... 86%  114M 0s
    ##  69100K .......... .......... .......... .......... .......... 86% 83.6M 0s
    ##  69150K .......... .......... .......... .......... .......... 86% 95.8M 0s
    ##  69200K .......... .......... .......... .......... .......... 86% 85.5M 0s
    ##  69250K .......... .......... .......... .......... .......... 86%  112M 0s
    ##  69300K .......... .......... .......... .......... .......... 86% 80.7M 0s
    ##  69350K .......... .......... .......... .......... .......... 86%  100M 0s
    ##  69400K .......... .......... .......... .......... .......... 86% 78.6M 0s
    ##  69450K .......... .......... .......... .......... .......... 86% 86.4M 0s
    ##  69500K .......... .......... .......... .......... .......... 87% 62.5M 0s
    ##  69550K .......... .......... .......... .......... .......... 87% 71.7M 0s
    ##  69600K .......... .......... .......... .......... .......... 87% 40.3M 0s
    ##  69650K .......... .......... .......... .......... .......... 87% 88.8M 0s
    ##  69700K .......... .......... .......... .......... .......... 87% 90.9M 0s
    ##  69750K .......... .......... .......... .......... .......... 87% 87.6M 0s
    ##  69800K .......... .......... .......... .......... .......... 87% 93.0M 0s
    ##  69850K .......... .......... .......... .......... .......... 87% 75.7M 0s
    ##  69900K .......... .......... .......... .......... .......... 87% 81.0M 0s
    ##  69950K .......... .......... .......... .......... .......... 87% 93.4M 0s
    ##  70000K .......... .......... .......... .......... .......... 87% 73.9M 0s
    ##  70050K .......... .......... .......... .......... .......... 87%  114M 0s
    ##  70100K .......... .......... .......... .......... .......... 87%  102M 0s
    ##  70150K .......... .......... .......... .......... .......... 87% 91.6M 0s
    ##  70200K .......... .......... .......... .......... .......... 87% 98.3M 0s
    ##  70250K .......... .......... .......... .......... .......... 87%  102M 0s
    ##  70300K .......... .......... .......... .......... .......... 88% 78.5M 0s
    ##  70350K .......... .......... .......... .......... .......... 88%  104M 0s
    ##  70400K .......... .......... .......... .......... .......... 88% 99.0M 0s
    ##  70450K .......... .......... .......... .......... .......... 88% 77.3M 0s
    ##  70500K .......... .......... .......... .......... .......... 88% 34.7M 0s
    ##  70550K .......... .......... .......... .......... .......... 88% 64.8M 0s
    ##  70600K .......... .......... .......... .......... .......... 88% 95.1M 0s
    ##  70650K .......... .......... .......... .......... .......... 88% 94.3M 0s
    ##  70700K .......... .......... .......... .......... .......... 88%  124M 0s
    ##  70750K .......... .......... .......... .......... .......... 88% 94.6M 0s
    ##  70800K .......... .......... .......... .......... .......... 88% 31.4M 0s
    ##  70850K .......... .......... .......... .......... .......... 88% 81.9M 0s
    ##  70900K .......... .......... .......... .......... .......... 88% 89.6M 0s
    ##  70950K .......... .......... .......... .......... .......... 88% 71.5M 0s
    ##  71000K .......... .......... .......... .......... .......... 88%  128M 0s
    ##  71050K .......... .......... .......... .......... .......... 88%  123M 0s
    ##  71100K .......... .......... .......... .......... .......... 89%  115M 0s
    ##  71150K .......... .......... .......... .......... .......... 89%  112M 0s
    ##  71200K .......... .......... .......... .......... .......... 89%  107M 0s
    ##  71250K .......... .......... .......... .......... .......... 89%  109M 0s
    ##  71300K .......... .......... .......... .......... .......... 89%  107M 0s
    ##  71350K .......... .......... .......... .......... .......... 89% 76.5M 0s
    ##  71400K .......... .......... .......... .......... .......... 89% 90.4M 0s
    ##  71450K .......... .......... .......... .......... .......... 89%  109M 0s
    ##  71500K .......... .......... .......... .......... .......... 89%  120M 0s
    ##  71550K .......... .......... .......... .......... .......... 89%  119M 0s
    ##  71600K .......... .......... .......... .......... .......... 89% 71.9M 0s
    ##  71650K .......... .......... .......... .......... .......... 89% 70.5M 0s
    ##  71700K .......... .......... .......... .......... .......... 89% 57.0M 0s
    ##  71750K .......... .......... .......... .......... .......... 89% 89.1M 0s
    ##  71800K .......... .......... .......... .......... .......... 89% 74.1M 0s
    ##  71850K .......... .......... .......... .......... .......... 89% 54.9M 0s
    ##  71900K .......... .......... .......... .......... .......... 90%  114M 0s
    ##  71950K .......... .......... .......... .......... .......... 90%  130M 0s
    ##  72000K .......... .......... .......... .......... .......... 90%  127M 0s
    ##  72050K .......... .......... .......... .......... .......... 90%  105M 0s
    ##  72100K .......... .......... .......... .......... .......... 90%  126M 0s
    ##  72150K .......... .......... .......... .......... .......... 90%  101M 0s
    ##  72200K .......... .......... .......... .......... .......... 90% 15.9M 0s
    ##  72250K .......... .......... .......... .......... .......... 90% 93.6M 0s
    ##  72300K .......... .......... .......... .......... .......... 90%  123M 0s
    ##  72350K .......... .......... .......... .......... .......... 90%  119M 0s
    ##  72400K .......... .......... .......... .......... .......... 90%  127M 0s
    ##  72450K .......... .......... .......... .......... .......... 90%  133M 0s
    ##  72500K .......... .......... .......... .......... .......... 90%  119M 0s
    ##  72550K .......... .......... .......... .......... .......... 90%  115M 0s
    ##  72600K .......... .......... .......... .......... .......... 90% 9.34M 0s
    ##  72650K .......... .......... .......... .......... .......... 90%  101M 0s
    ##  72700K .......... .......... .......... .......... .......... 91% 96.5M 0s
    ##  72750K .......... .......... .......... .......... .......... 91% 93.5M 0s
    ##  72800K .......... .......... .......... .......... .......... 91%  119M 0s
    ##  72850K .......... .......... .......... .......... .......... 91%  105M 0s
    ##  72900K .......... .......... .......... .......... .......... 91% 75.1M 0s
    ##  72950K .......... .......... .......... .......... .......... 91%  119M 0s
    ##  73000K .......... .......... .......... .......... .......... 91% 71.0M 0s
    ##  73050K .......... .......... .......... .......... .......... 91% 54.1M 0s
    ##  73100K .......... .......... .......... .......... .......... 91% 47.4M 0s
    ##  73150K .......... .......... .......... .......... .......... 91% 75.6M 0s
    ##  73200K .......... .......... .......... .......... .......... 91% 98.9M 0s
    ##  73250K .......... .......... .......... .......... .......... 91% 80.9M 0s
    ##  73300K .......... .......... .......... .......... .......... 91%  112M 0s
    ##  73350K .......... .......... .......... .......... .......... 91%  135M 0s
    ##  73400K .......... .......... .......... .......... .......... 91% 80.1M 0s
    ##  73450K .......... .......... .......... .......... .......... 91% 85.4M 0s
    ##  73500K .......... .......... .......... .......... .......... 92%  103M 0s
    ##  73550K .......... .......... .......... .......... .......... 92%  130M 0s
    ##  73600K .......... .......... .......... .......... .......... 92% 99.9M 0s
    ##  73650K .......... .......... .......... .......... .......... 92%  112M 0s
    ##  73700K .......... .......... .......... .......... .......... 92% 78.8M 0s
    ##  73750K .......... .......... .......... .......... .......... 92%  109M 0s
    ##  73800K .......... .......... .......... .......... .......... 92% 39.6M 0s
    ##  73850K .......... .......... .......... .......... .......... 92% 91.3M 0s
    ##  73900K .......... .......... .......... .......... .......... 92%  112M 0s
    ##  73950K .......... .......... .......... .......... .......... 92%  122M 0s
    ##  74000K .......... .......... .......... .......... .......... 92%  127M 0s
    ##  74050K .......... .......... .......... .......... .......... 92% 83.9M 0s
    ##  74100K .......... .......... .......... .......... .......... 92%  118M 0s
    ##  74150K .......... .......... .......... .......... .......... 92%  136M 0s
    ##  74200K .......... .......... .......... .......... .......... 92%  108M 0s
    ##  74250K .......... .......... .......... .......... .......... 92% 98.7M 0s
    ##  74300K .......... .......... .......... .......... .......... 93% 37.1M 0s
    ##  74350K .......... .......... .......... .......... .......... 93%  112M 0s
    ##  74400K .......... .......... .......... .......... .......... 93% 30.0M 0s
    ##  74450K .......... .......... .......... .......... .......... 93%  124M 0s
    ##  74500K .......... .......... .......... .......... .......... 93%  104M 0s
    ##  74550K .......... .......... .......... .......... .......... 93%  150M 0s
    ##  74600K .......... .......... .......... .......... .......... 93%  136M 0s
    ##  74650K .......... .......... .......... .......... .......... 93%  120M 0s
    ##  74700K .......... .......... .......... .......... .......... 93% 16.2M 0s
    ##  74750K .......... .......... .......... .......... .......... 93% 47.0M 0s
    ##  74800K .......... .......... .......... .......... .......... 93% 28.5M 0s
    ##  74850K .......... .......... .......... .......... .......... 93%  138M 0s
    ##  74900K .......... .......... .......... .......... .......... 93% 64.2M 0s
    ##  74950K .......... .......... .......... .......... .......... 93% 98.6M 0s
    ##  75000K .......... .......... .......... .......... .......... 93%  120M 0s
    ##  75050K .......... .......... .......... .......... .......... 93%  139M 0s
    ##  75100K .......... .......... .......... .......... .......... 94%  126M 0s
    ##  75150K .......... .......... .......... .......... .......... 94%  143M 0s
    ##  75200K .......... .......... .......... .......... .......... 94% 74.7M 0s
    ##  75250K .......... .......... .......... .......... .......... 94%  117M 0s
    ##  75300K .......... .......... .......... .......... .......... 94% 96.3M 0s
    ##  75350K .......... .......... .......... .......... .......... 94%  143M 0s
    ##  75400K .......... .......... .......... .......... .......... 94% 79.5M 0s
    ##  75450K .......... .......... .......... .......... .......... 94%  113M 0s
    ##  75500K .......... .......... .......... .......... .......... 94%  117M 0s
    ##  75550K .......... .......... .......... .......... .......... 94%  144M 0s
    ##  75600K .......... .......... .......... .......... .......... 94% 92.5M 0s
    ##  75650K .......... .......... .......... .......... .......... 94% 68.5M 0s
    ##  75700K .......... .......... .......... .......... .......... 94% 87.5M 0s
    ##  75750K .......... .......... .......... .......... .......... 94% 97.5M 0s
    ##  75800K .......... .......... .......... .......... .......... 94% 98.9M 0s
    ##  75850K .......... .......... .......... .......... .......... 94%  117M 0s
    ##  75900K .......... .......... .......... .......... .......... 95% 87.7M 0s
    ##  75950K .......... .......... .......... .......... .......... 95% 62.4M 0s
    ##  76000K .......... .......... .......... .......... .......... 95% 79.9M 0s
    ##  76050K .......... .......... .......... .......... .......... 95% 74.1M 0s
    ##  76100K .......... .......... .......... .......... .......... 95%  105M 0s
    ##  76150K .......... .......... .......... .......... .......... 95%  119M 0s
    ##  76200K .......... .......... .......... .......... .......... 95%  123M 0s
    ##  76250K .......... .......... .......... .......... .......... 95%  142M 0s
    ##  76300K .......... .......... .......... .......... .......... 95% 99.3M 0s
    ##  76350K .......... .......... .......... .......... .......... 95%  123M 0s
    ##  76400K .......... .......... .......... .......... .......... 95% 52.8M 0s
    ##  76450K .......... .......... .......... .......... .......... 95% 65.2M 0s
    ##  76500K .......... .......... .......... .......... .......... 95% 35.5M 0s
    ##  76550K .......... .......... .......... .......... .......... 95%  105M 0s
    ##  76600K .......... .......... .......... .......... .......... 95% 98.6M 0s
    ##  76650K .......... .......... .......... .......... .......... 95%  121M 0s
    ##  76700K .......... .......... .......... .......... .......... 96%  139M 0s
    ##  76750K .......... .......... .......... .......... .......... 96% 61.0M 0s
    ##  76800K .......... .......... .......... .......... .......... 96%  115M 0s
    ##  76850K .......... .......... .......... .......... .......... 96% 46.1M 0s
    ##  76900K .......... .......... .......... .......... .......... 96%  124M 0s
    ##  76950K .......... .......... .......... .......... .......... 96% 84.9M 0s
    ##  77000K .......... .......... .......... .......... .......... 96%  144M 0s
    ##  77050K .......... .......... .......... .......... .......... 96% 43.3M 0s
    ##  77100K .......... .......... .......... .......... .......... 96% 96.0M 0s
    ##  77150K .......... .......... .......... .......... .......... 96%  120M 0s
    ##  77200K .......... .......... .......... .......... .......... 96%  132M 0s
    ##  77250K .......... .......... .......... .......... .......... 96% 12.5M 0s
    ##  77300K .......... .......... .......... .......... .......... 96%  118M 0s
    ##  77350K .......... .......... .......... .......... .......... 96%  137M 0s
    ##  77400K .......... .......... .......... .......... .......... 96%  126M 0s
    ##  77450K .......... .......... .......... .......... .......... 96%  138M 0s
    ##  77500K .......... .......... .......... .......... .......... 97%  123M 0s
    ##  77550K .......... .......... .......... .......... .......... 97%  124M 0s
    ##  77600K .......... .......... .......... .......... .......... 97%  147M 0s
    ##  77650K .......... .......... .......... .......... .......... 97% 21.9M 0s
    ##  77700K .......... .......... .......... .......... .......... 97% 86.5M 0s
    ##  77750K .......... .......... .......... .......... .......... 97% 81.7M 0s
    ##  77800K .......... .......... .......... .......... .......... 97% 88.6M 0s
    ##  77850K .......... .......... .......... .......... .......... 97% 86.3M 0s
    ##  77900K .......... .......... .......... .......... .......... 97%  112M 0s
    ##  77950K .......... .......... .......... .......... .......... 97%  147M 0s
    ##  78000K .......... .......... .......... .......... .......... 97%  133M 0s
    ##  78050K .......... .......... .......... .......... .......... 97%  125M 0s
    ##  78100K .......... .......... .......... .......... .......... 97% 57.2M 0s
    ##  78150K .......... .......... .......... .......... .......... 97% 90.8M 0s
    ##  78200K .......... .......... .......... .......... .......... 97%  113M 0s
    ##  78250K .......... .......... .......... .......... .......... 97% 38.3M 0s
    ##  78300K .......... .......... .......... .......... .......... 98%  105M 0s
    ##  78350K .......... .......... .......... .......... .......... 98%  103M 0s
    ##  78400K .......... .......... .......... .......... .......... 98%  102M 0s
    ##  78450K .......... .......... .......... .......... .......... 98%  158M 0s
    ##  78500K .......... .......... .......... .......... .......... 98% 75.1M 0s
    ##  78550K .......... .......... .......... .......... .......... 98%  162M 0s
    ##  78600K .......... .......... .......... .......... .......... 98%  123M 0s
    ##  78650K .......... .......... .......... .......... .......... 98% 49.8M 0s
    ##  78700K .......... .......... .......... .......... .......... 98% 83.8M 0s
    ##  78750K .......... .......... .......... .......... .......... 98%  141M 0s
    ##  78800K .......... .......... .......... .......... .......... 98% 54.4M 0s
    ##  78850K .......... .......... .......... .......... .......... 98%  128M 0s
    ##  78900K .......... .......... .......... .......... .......... 98% 54.5M 0s
    ##  78950K .......... .......... .......... .......... .......... 98%  112M 0s
    ##  79000K .......... .......... .......... .......... .......... 98%  118M 0s
    ##  79050K .......... .......... .......... .......... .......... 98% 26.3M 0s
    ##  79100K .......... .......... .......... .......... .......... 99% 82.3M 0s
    ##  79150K .......... .......... .......... .......... .......... 99%  171M 0s
    ##  79200K .......... .......... .......... .......... .......... 99%  140M 0s
    ##  79250K .......... .......... .......... .......... .......... 99% 73.9M 0s
    ##  79300K .......... .......... .......... .......... .......... 99%  111M 0s
    ##  79350K .......... .......... .......... .......... .......... 99%  158M 0s
    ##  79400K .......... .......... .......... .......... .......... 99%  147M 0s
    ##  79450K .......... .......... .......... .......... .......... 99% 16.1M 0s
    ##  79500K .......... .......... .......... .......... .......... 99% 60.0M 0s
    ##  79550K .......... .......... .......... .......... .......... 99% 43.4M 0s
    ##  79600K .......... .......... .......... .......... .......... 99% 99.2M 0s
    ##  79650K .......... .......... .......... .......... .......... 99%  107M 0s
    ##  79700K .......... .......... .......... .......... .......... 99%  114M 0s
    ##  79750K .......... .......... .......... .......... .......... 99%  146M 0s
    ##  79800K .......... .......... .......... .......... .......... 99%  120M 0s
    ##  79850K .......... .......... .......... .......... .......... 99%  129M 0s
    ##  79900K .......... .......... ..                              100%  141M=1.5s
    ## 
    ## 2020-12-03 15:16:51 (53.0 MB/s) - ‘silva_species_assignment_v138.fa.gz.2’ saved [81840166/81840166]

``` r
taxa <- addSpecies(taxa, "~/silva_species_assignment_v138.fa.gz")
```

``` r
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum         Class         Order           Family          
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      Genus         Species
    ## [1,] NA            NA     
    ## [2,] NA            NA     
    ## [3,] NA            NA     
    ## [4,] NA            NA     
    ## [5,] "Bacteroides" NA     
    ## [6,] NA            NA

On a construit une table qui commence au niveau du reigne et qui va
jusqu’à l’espèce. Les résultats montrent que les séquences n’ont pas pu
être assignées jusqu’à l’espèce voire même au niveua du genre
hormispour la séquence 5.

## Evaluer la précision de l’assignation taxonomique

``` r
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

    ## DADA2 inferred 20 sample sequences present in the Mock community.

``` r
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

    ## Of those, 20 were exact matches to the expected reference sequences.

# Conclusion

``` r
save.image(file="02_Dada2_tutorial_FinalEnv")
```
