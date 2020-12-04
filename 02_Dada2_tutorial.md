Contrôle continu 1 : Analyse des données avec dada2
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

    ## --2020-12-04 08:43:57--  https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 137973851 (132M) [application/octet-stream]
    ## Saving to: ‘silva_nr99_v138_train_set.fa.gz.6’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 14.5M 9s
    ##     50K .......... .......... .......... .......... ..........  0% 13.9M 9s
    ##    100K .......... .......... .......... .......... ..........  0% 9.65M 11s
    ##    150K .......... .......... .......... .......... ..........  0% 28.6M 9s
    ##    200K .......... .......... .......... .......... ..........  0% 19.6M 9s
    ##    250K .......... .......... .......... .......... ..........  0% 54.6M 8s
    ##    300K .......... .......... .......... .......... ..........  0% 18.0M 8s
    ##    350K .......... .......... .......... .......... ..........  0% 51.1M 7s
    ##    400K .......... .......... .......... .......... ..........  0% 15.6M 7s
    ##    450K .......... .......... .......... .......... ..........  0% 12.5M 7s
    ##    500K .......... .......... .......... .......... ..........  0%  109M 7s
    ##    550K .......... .......... .......... .......... ..........  0%  137M 6s
    ##    600K .......... .......... .......... .......... ..........  0% 13.0M 7s
    ##    650K .......... .......... .......... .......... ..........  0%  135M 6s
    ##    700K .......... .......... .......... .......... ..........  0% 16.9M 6s
    ##    750K .......... .......... .......... .......... ..........  0% 10.4M 7s
    ##    800K .......... .......... .......... .......... ..........  0% 87.1M 6s
    ##    850K .......... .......... .......... .......... ..........  0% 13.0M 7s
    ##    900K .......... .......... .......... .......... ..........  0%  103M 6s
    ##    950K .......... .......... .......... .......... ..........  0% 13.6M 6s
    ##   1000K .......... .......... .......... .......... ..........  0% 80.1M 6s
    ##   1050K .......... .......... .......... .......... ..........  0% 16.0M 6s
    ##   1100K .......... .......... .......... .......... ..........  0% 96.5M 6s
    ##   1150K .......... .......... .......... .......... ..........  0% 13.3M 6s
    ##   1200K .......... .......... .......... .......... ..........  0%  109M 6s
    ##   1250K .......... .......... .......... .......... ..........  0% 14.5M 6s
    ##   1300K .......... .......... .......... .......... ..........  1%  127M 6s
    ##   1350K .......... .......... .......... .......... ..........  1% 11.6M 6s
    ##   1400K .......... .......... .......... .......... ..........  1%  123M 6s
    ##   1450K .......... .......... .......... .......... ..........  1% 14.2M 6s
    ##   1500K .......... .......... .......... .......... ..........  1% 79.4M 6s
    ##   1550K .......... .......... .......... .......... ..........  1% 17.8M 6s
    ##   1600K .......... .......... .......... .......... ..........  1% 44.3M 6s
    ##   1650K .......... .......... .......... .......... ..........  1% 14.4M 6s
    ##   1700K .......... .......... .......... .......... ..........  1% 81.2M 6s
    ##   1750K .......... .......... .......... .......... ..........  1% 21.6M 6s
    ##   1800K .......... .......... .......... .......... ..........  1% 30.6M 6s
    ##   1850K .......... .......... .......... .......... ..........  1%  100M 6s
    ##   1900K .......... .......... .......... .......... ..........  1% 19.1M 6s
    ##   1950K .......... .......... .......... .......... ..........  1% 44.9M 6s
    ##   2000K .......... .......... .......... .......... ..........  1% 41.1M 6s
    ##   2050K .......... .......... .......... .......... ..........  1% 35.0M 6s
    ##   2100K .......... .......... .......... .......... ..........  1% 29.6M 6s
    ##   2150K .......... .......... .......... .......... ..........  1% 20.6M 6s
    ##   2200K .......... .......... .......... .......... ..........  1% 75.1M 5s
    ##   2250K .......... .......... .......... .......... ..........  1% 73.8M 5s
    ##   2300K .......... .......... .......... .......... ..........  1% 20.9M 5s
    ##   2350K .......... .......... .......... .......... ..........  1% 48.9M 5s
    ##   2400K .......... .......... .......... .......... ..........  1% 16.7M 5s
    ##   2450K .......... .......... .......... .......... ..........  1%  121M 5s
    ##   2500K .......... .......... .......... .......... ..........  1% 74.6M 5s
    ##   2550K .......... .......... .......... .......... ..........  1% 25.3M 5s
    ##   2600K .......... .......... .......... .......... ..........  1% 16.8M 5s
    ##   2650K .......... .......... .......... .......... ..........  2% 97.9M 5s
    ##   2700K .......... .......... .......... .......... ..........  2% 17.2M 5s
    ##   2750K .......... .......... .......... .......... ..........  2% 91.3M 5s
    ##   2800K .......... .......... .......... .......... ..........  2%  105M 5s
    ##   2850K .......... .......... .......... .......... ..........  2% 16.7M 5s
    ##   2900K .......... .......... .......... .......... ..........  2% 70.9M 5s
    ##   2950K .......... .......... .......... .......... ..........  2%  122M 5s
    ##   3000K .......... .......... .......... .......... ..........  2% 14.3M 5s
    ##   3050K .......... .......... .......... .......... ..........  2%  123M 5s
    ##   3100K .......... .......... .......... .......... ..........  2% 16.1M 5s
    ##   3150K .......... .......... .......... .......... ..........  2% 87.9M 5s
    ##   3200K .......... .......... .......... .......... ..........  2%  111M 5s
    ##   3250K .......... .......... .......... .......... ..........  2% 15.8M 5s
    ##   3300K .......... .......... .......... .......... ..........  2% 84.5M 5s
    ##   3350K .......... .......... .......... .......... ..........  2%  117M 5s
    ##   3400K .......... .......... .......... .......... ..........  2% 17.4M 5s
    ##   3450K .......... .......... .......... .......... ..........  2% 87.2M 5s
    ##   3500K .......... .......... .......... .......... ..........  2% 63.3M 5s
    ##   3550K .......... .......... .......... .......... ..........  2% 27.1M 5s
    ##   3600K .......... .......... .......... .......... ..........  2% 28.6M 5s
    ##   3650K .......... .......... .......... .......... ..........  2%  122M 5s
    ##   3700K .......... .......... .......... .......... ..........  2% 46.7M 5s
    ##   3750K .......... .......... .......... .......... ..........  2% 18.1M 5s
    ##   3800K .......... .......... .......... .......... ..........  2% 83.2M 5s
    ##   3850K .......... .......... .......... .......... ..........  2%  131M 5s
    ##   3900K .......... .......... .......... .......... ..........  2% 14.2M 5s
    ##   3950K .......... .......... .......... .......... ..........  2%  105M 5s
    ##   4000K .......... .......... .......... .......... ..........  3% 13.1M 5s
    ##   4050K .......... .......... .......... .......... ..........  3%  130M 5s
    ##   4100K .......... .......... .......... .......... ..........  3%  100M 5s
    ##   4150K .......... .......... .......... .......... ..........  3%  132M 5s
    ##   4200K .......... .......... .......... .......... ..........  3% 20.2M 5s
    ##   4250K .......... .......... .......... .......... ..........  3% 65.0M 5s
    ##   4300K .......... .......... .......... .......... ..........  3% 16.8M 5s
    ##   4350K .......... .......... .......... .......... ..........  3%  144M 5s
    ##   4400K .......... .......... .......... .......... ..........  3% 49.8M 5s
    ##   4450K .......... .......... .......... .......... ..........  3% 16.9M 5s
    ##   4500K .......... .......... .......... .......... ..........  3% 95.9M 5s
    ##   4550K .......... .......... .......... .......... ..........  3% 27.2M 5s
    ##   4600K .......... .......... .......... .......... ..........  3% 98.1M 4s
    ##   4650K .......... .......... .......... .......... ..........  3% 29.6M 4s
    ##   4700K .......... .......... .......... .......... ..........  3% 23.9M 4s
    ##   4750K .......... .......... .......... .......... ..........  3% 7.23M 5s
    ##   4800K .......... .......... .......... .......... ..........  3% 62.3M 5s
    ##   4850K .......... .......... .......... .......... ..........  3%  102M 5s
    ##   4900K .......... .......... .......... .......... ..........  3%  119M 5s
    ##   4950K .......... .......... .......... .......... ..........  3% 8.26M 5s
    ##   5000K .......... .......... .......... .......... ..........  3% 54.1M 5s
    ##   5050K .......... .......... .......... .......... ..........  3%  110M 5s
    ##   5100K .......... .......... .......... .......... ..........  3%  110M 5s
    ##   5150K .......... .......... .......... .......... ..........  3% 18.5M 5s
    ##   5200K .......... .......... .......... .......... ..........  3% 6.06M 5s
    ##   5250K .......... .......... .......... .......... ..........  3%  126M 5s
    ##   5300K .......... .......... .......... .......... ..........  3%  126M 5s
    ##   5350K .......... .......... .......... .......... ..........  4%  112M 5s
    ##   5400K .......... .......... .......... .......... ..........  4% 19.7M 5s
    ##   5450K .......... .......... .......... .......... ..........  4% 24.4M 5s
    ##   5500K .......... .......... .......... .......... ..........  4% 68.2M 5s
    ##   5550K .......... .......... .......... .......... ..........  4% 85.5M 5s
    ##   5600K .......... .......... .......... .......... ..........  4% 88.7M 5s
    ##   5650K .......... .......... .......... .......... ..........  4% 20.0M 5s
    ##   5700K .......... .......... .......... .......... ..........  4% 65.3M 5s
    ##   5750K .......... .......... .......... .......... ..........  4% 83.7M 4s
    ##   5800K .......... .......... .......... .......... ..........  4% 80.0M 4s
    ##   5850K .......... .......... .......... .......... ..........  4% 36.1M 4s
    ##   5900K .......... .......... .......... .......... ..........  4% 67.5M 4s
    ##   5950K .......... .......... .......... .......... ..........  4% 49.2M 4s
    ##   6000K .......... .......... .......... .......... ..........  4% 17.9M 4s
    ##   6050K .......... .......... .......... .......... ..........  4% 77.4M 4s
    ##   6100K .......... .......... .......... .......... ..........  4% 70.9M 4s
    ##   6150K .......... .......... .......... .......... ..........  4% 93.2M 4s
    ##   6200K .......... .......... .......... .......... ..........  4% 12.1M 4s
    ##   6250K .......... .......... .......... .......... ..........  4% 69.2M 4s
    ##   6300K .......... .......... .......... .......... ..........  4% 90.7M 4s
    ##   6350K .......... .......... .......... .......... ..........  4%  102M 4s
    ##   6400K .......... .......... .......... .......... ..........  4% 22.2M 4s
    ##   6450K .......... .......... .......... .......... ..........  4% 65.5M 4s
    ##   6500K .......... .......... .......... .......... ..........  4% 53.3M 4s
    ##   6550K .......... .......... .......... .......... ..........  4%  104M 4s
    ##   6600K .......... .......... .......... .......... ..........  4% 89.0M 4s
    ##   6650K .......... .......... .......... .......... ..........  4% 19.2M 4s
    ##   6700K .......... .......... .......... .......... ..........  5% 31.6M 4s
    ##   6750K .......... .......... .......... .......... ..........  5% 44.3M 4s
    ##   6800K .......... .......... .......... .......... ..........  5% 53.0M 4s
    ##   6850K .......... .......... .......... .......... ..........  5% 94.9M 4s
    ##   6900K .......... .......... .......... .......... ..........  5% 53.7M 4s
    ##   6950K .......... .......... .......... .......... ..........  5% 69.8M 4s
    ##   7000K .......... .......... .......... .......... ..........  5% 90.9M 4s
    ##   7050K .......... .......... .......... .......... ..........  5% 75.0M 4s
    ##   7100K .......... .......... .......... .......... ..........  5% 59.9M 4s
    ##   7150K .......... .......... .......... .......... ..........  5% 31.5M 4s
    ##   7200K .......... .......... .......... .......... ..........  5% 83.3M 4s
    ##   7250K .......... .......... .......... .......... ..........  5% 82.2M 4s
    ##   7300K .......... .......... .......... .......... ..........  5% 39.7M 4s
    ##   7350K .......... .......... .......... .......... ..........  5% 39.0M 4s
    ##   7400K .......... .......... .......... .......... ..........  5% 75.7M 4s
    ##   7450K .......... .......... .......... .......... ..........  5% 93.4M 4s
    ##   7500K .......... .......... .......... .......... ..........  5% 46.0M 4s
    ##   7550K .......... .......... .......... .......... ..........  5% 26.5M 4s
    ##   7600K .......... .......... .......... .......... ..........  5% 11.9M 4s
    ##   7650K .......... .......... .......... .......... ..........  5% 56.8M 4s
    ##   7700K .......... .......... .......... .......... ..........  5% 86.1M 4s
    ##   7750K .......... .......... .......... .......... ..........  5% 98.0M 4s
    ##   7800K .......... .......... .......... .......... ..........  5% 19.8M 4s
    ##   7850K .......... .......... .......... .......... ..........  5% 88.5M 4s
    ##   7900K .......... .......... .......... .......... ..........  5% 97.8M 4s
    ##   7950K .......... .......... .......... .......... ..........  5% 84.8M 4s
    ##   8000K .......... .......... .......... .......... ..........  5% 18.2M 4s
    ##   8050K .......... .......... .......... .......... ..........  6% 81.4M 4s
    ##   8100K .......... .......... .......... .......... ..........  6% 81.5M 4s
    ##   8150K .......... .......... .......... .......... ..........  6%  103M 4s
    ##   8200K .......... .......... .......... .......... ..........  6% 16.4M 4s
    ##   8250K .......... .......... .......... .......... ..........  6% 72.6M 4s
    ##   8300K .......... .......... .......... .......... ..........  6%  100M 4s
    ##   8350K .......... .......... .......... .......... ..........  6% 15.5M 4s
    ##   8400K .......... .......... .......... .......... ..........  6% 58.2M 4s
    ##   8450K .......... .......... .......... .......... ..........  6% 70.3M 4s
    ##   8500K .......... .......... .......... .......... ..........  6% 94.0M 4s
    ##   8550K .......... .......... .......... .......... ..........  6% 35.6M 4s
    ##   8600K .......... .......... .......... .......... ..........  6% 29.0M 4s
    ##   8650K .......... .......... .......... .......... ..........  6% 79.4M 4s
    ##   8700K .......... .......... .......... .......... ..........  6% 80.0M 4s
    ##   8750K .......... .......... .......... .......... ..........  6%  116M 4s
    ##   8800K .......... .......... .......... .......... ..........  6% 31.6M 4s
    ##   8850K .......... .......... .......... .......... ..........  6% 75.1M 4s
    ##   8900K .......... .......... .......... .......... ..........  6% 31.5M 4s
    ##   8950K .......... .......... .......... .......... ..........  6% 90.4M 4s
    ##   9000K .......... .......... .......... .......... ..........  6% 30.6M 4s
    ##   9050K .......... .......... .......... .......... ..........  6%  109M 4s
    ##   9100K .......... .......... .......... .......... ..........  6% 56.5M 4s
    ##   9150K .......... .......... .......... .......... ..........  6% 53.7M 4s
    ##   9200K .......... .......... .......... .......... ..........  6% 34.0M 4s
    ##   9250K .......... .......... .......... .......... ..........  6% 51.9M 4s
    ##   9300K .......... .......... .......... .......... ..........  6% 72.8M 4s
    ##   9350K .......... .......... .......... .......... ..........  6% 81.1M 4s
    ##   9400K .......... .......... .......... .......... ..........  7% 26.0M 4s
    ##   9450K .......... .......... .......... .......... ..........  7% 67.8M 4s
    ##   9500K .......... .......... .......... .......... ..........  7% 88.0M 4s
    ##   9550K .......... .......... .......... .......... ..........  7%  101M 4s
    ##   9600K .......... .......... .......... .......... ..........  7% 18.2M 4s
    ##   9650K .......... .......... .......... .......... ..........  7%  119M 4s
    ##   9700K .......... .......... .......... .......... ..........  7% 77.3M 4s
    ##   9750K .......... .......... .......... .......... ..........  7% 95.4M 4s
    ##   9800K .......... .......... .......... .......... ..........  7% 23.9M 4s
    ##   9850K .......... .......... .......... .......... ..........  7% 86.8M 4s
    ##   9900K .......... .......... .......... .......... ..........  7% 84.8M 4s
    ##   9950K .......... .......... .......... .......... ..........  7% 28.1M 4s
    ##  10000K .......... .......... .......... .......... ..........  7% 73.6M 4s
    ##  10050K .......... .......... .......... .......... ..........  7% 86.8M 4s
    ##  10100K .......... .......... .......... .......... ..........  7% 42.9M 4s
    ##  10150K .......... .......... .......... .......... ..........  7% 42.6M 4s
    ##  10200K .......... .......... .......... .......... ..........  7% 46.4M 4s
    ##  10250K .......... .......... .......... .......... ..........  7%  102M 4s
    ##  10300K .......... .......... .......... .......... ..........  7% 68.3M 4s
    ##  10350K .......... .......... .......... .......... ..........  7% 35.3M 4s
    ##  10400K .......... .......... .......... .......... ..........  7% 39.0M 4s
    ##  10450K .......... .......... .......... .......... ..........  7% 71.5M 4s
    ##  10500K .......... .......... .......... .......... ..........  7% 95.4M 4s
    ##  10550K .......... .......... .......... .......... ..........  7% 35.3M 4s
    ##  10600K .......... .......... .......... .......... ..........  7% 46.4M 4s
    ##  10650K .......... .......... .......... .......... ..........  7% 68.5M 4s
    ##  10700K .......... .......... .......... .......... ..........  7% 34.9M 4s
    ##  10750K .......... .......... .......... .......... ..........  8% 72.7M 4s
    ##  10800K .......... .......... .......... .......... ..........  8% 41.0M 4s
    ##  10850K .......... .......... .......... .......... ..........  8% 63.6M 4s
    ##  10900K .......... .......... .......... .......... ..........  8% 31.7M 4s
    ##  10950K .......... .......... .......... .......... ..........  8% 87.2M 4s
    ##  11000K .......... .......... .......... .......... ..........  8% 32.2M 4s
    ##  11050K .......... .......... .......... .......... ..........  8% 94.1M 4s
    ##  11100K .......... .......... .......... .......... ..........  8% 38.4M 3s
    ##  11150K .......... .......... .......... .......... ..........  8% 55.4M 3s
    ##  11200K .......... .......... .......... .......... ..........  8% 72.3M 3s
    ##  11250K .......... .......... .......... .......... ..........  8% 78.0M 3s
    ##  11300K .......... .......... .......... .......... ..........  8% 26.4M 3s
    ##  11350K .......... .......... .......... .......... ..........  8% 43.5M 3s
    ##  11400K .......... .......... .......... .......... ..........  8% 69.2M 3s
    ##  11450K .......... .......... .......... .......... ..........  8% 97.2M 3s
    ##  11500K .......... .......... .......... .......... ..........  8% 36.9M 3s
    ##  11550K .......... .......... .......... .......... ..........  8% 52.4M 3s
    ##  11600K .......... .......... .......... .......... ..........  8% 37.4M 3s
    ##  11650K .......... .......... .......... .......... ..........  8% 71.6M 3s
    ##  11700K .......... .......... .......... .......... ..........  8%  100M 3s
    ##  11750K .......... .......... .......... .......... ..........  8% 40.0M 3s
    ##  11800K .......... .......... .......... .......... ..........  8% 29.6M 3s
    ##  11850K .......... .......... .......... .......... ..........  8% 83.5M 3s
    ##  11900K .......... .......... .......... .......... ..........  8% 81.5M 3s
    ##  11950K .......... .......... .......... .......... ..........  8%  107M 3s
    ##  12000K .......... .......... .......... .......... ..........  8% 24.5M 3s
    ##  12050K .......... .......... .......... .......... ..........  8% 77.6M 3s
    ##  12100K .......... .......... .......... .......... ..........  9% 80.4M 3s
    ##  12150K .......... .......... .......... .......... ..........  9% 96.8M 3s
    ##  12200K .......... .......... .......... .......... ..........  9% 21.1M 3s
    ##  12250K .......... .......... .......... .......... ..........  9% 79.7M 3s
    ##  12300K .......... .......... .......... .......... ..........  9% 97.2M 3s
    ##  12350K .......... .......... .......... .......... ..........  9% 27.2M 3s
    ##  12400K .......... .......... .......... .......... ..........  9% 64.7M 3s
    ##  12450K .......... .......... .......... .......... ..........  9% 51.6M 3s
    ##  12500K .......... .......... .......... .......... ..........  9% 89.5M 3s
    ##  12550K .......... .......... .......... .......... ..........  9% 90.7M 3s
    ##  12600K .......... .......... .......... .......... ..........  9%  108M 3s
    ##  12650K .......... .......... .......... .......... ..........  9% 28.5M 3s
    ##  12700K .......... .......... .......... .......... ..........  9% 84.0M 3s
    ##  12750K .......... .......... .......... .......... ..........  9% 57.2M 3s
    ##  12800K .......... .......... .......... .......... ..........  9% 28.7M 3s
    ##  12850K .......... .......... .......... .......... ..........  9% 86.7M 3s
    ##  12900K .......... .......... .......... .......... ..........  9% 89.7M 3s
    ##  12950K .......... .......... .......... .......... ..........  9% 68.9M 3s
    ##  13000K .......... .......... .......... .......... ..........  9% 40.1M 3s
    ##  13050K .......... .......... .......... .......... ..........  9% 30.7M 3s
    ##  13100K .......... .......... .......... .......... ..........  9% 72.7M 3s
    ##  13150K .......... .......... .......... .......... ..........  9% 99.8M 3s
    ##  13200K .......... .......... .......... .......... ..........  9% 76.6M 3s
    ##  13250K .......... .......... .......... .......... ..........  9% 21.9M 3s
    ##  13300K .......... .......... .......... .......... ..........  9% 63.5M 3s
    ##  13350K .......... .......... .......... .......... ..........  9% 53.4M 3s
    ##  13400K .......... .......... .......... .......... ..........  9%  100M 3s
    ##  13450K .......... .......... .......... .......... .......... 10% 30.9M 3s
    ##  13500K .......... .......... .......... .......... .......... 10% 76.2M 3s
    ##  13550K .......... .......... .......... .......... .......... 10%  102M 3s
    ##  13600K .......... .......... .......... .......... .......... 10%  100M 3s
    ##  13650K .......... .......... .......... .......... .......... 10% 21.1M 3s
    ##  13700K .......... .......... .......... .......... .......... 10% 31.1M 3s
    ##  13750K .......... .......... .......... .......... .......... 10%  120M 3s
    ##  13800K .......... .......... .......... .......... .......... 10% 68.6M 3s
    ##  13850K .......... .......... .......... .......... .......... 10% 49.6M 3s
    ##  13900K .......... .......... .......... .......... .......... 10% 23.4M 3s
    ##  13950K .......... .......... .......... .......... .......... 10%  111M 3s
    ##  14000K .......... .......... .......... .......... .......... 10% 88.6M 3s
    ##  14050K .......... .......... .......... .......... .......... 10%  108M 3s
    ##  14100K .......... .......... .......... .......... .......... 10% 20.2M 3s
    ##  14150K .......... .......... .......... .......... .......... 10%  107M 3s
    ##  14200K .......... .......... .......... .......... .......... 10% 96.8M 3s
    ##  14250K .......... .......... .......... .......... .......... 10% 98.7M 3s
    ##  14300K .......... .......... .......... .......... .......... 10% 21.7M 3s
    ##  14350K .......... .......... .......... .......... .......... 10% 90.6M 3s
    ##  14400K .......... .......... .......... .......... .......... 10% 99.1M 3s
    ##  14450K .......... .......... .......... .......... .......... 10% 22.2M 3s
    ##  14500K .......... .......... .......... .......... .......... 10% 66.9M 3s
    ##  14550K .......... .......... .......... .......... .......... 10% 84.3M 3s
    ##  14600K .......... .......... .......... .......... .......... 10% 97.2M 3s
    ##  14650K .......... .......... .......... .......... .......... 10%  102M 3s
    ##  14700K .......... .......... .......... .......... .......... 10% 29.6M 3s
    ##  14750K .......... .......... .......... .......... .......... 10% 62.7M 3s
    ##  14800K .......... .......... .......... .......... .......... 11%  102M 3s
    ##  14850K .......... .......... .......... .......... .......... 11%  126M 3s
    ##  14900K .......... .......... .......... .......... .......... 11% 22.5M 3s
    ##  14950K .......... .......... .......... .......... .......... 11%  113M 3s
    ##  15000K .......... .......... .......... .......... .......... 11% 63.6M 3s
    ##  15050K .......... .......... .......... .......... .......... 11% 79.9M 3s
    ##  15100K .......... .......... .......... .......... .......... 11% 17.7M 3s
    ##  15150K .......... .......... .......... .......... .......... 11% 66.9M 3s
    ##  15200K .......... .......... .......... .......... .......... 11%  105M 3s
    ##  15250K .......... .......... .......... .......... .......... 11%  111M 3s
    ##  15300K .......... .......... .......... .......... .......... 11% 19.8M 3s
    ##  15350K .......... .......... .......... .......... .......... 11% 80.7M 3s
    ##  15400K .......... .......... .......... .......... .......... 11% 71.6M 3s
    ##  15450K .......... .......... .......... .......... .......... 11%  117M 3s
    ##  15500K .......... .......... .......... .......... .......... 11% 92.9M 3s
    ##  15550K .......... .......... .......... .......... .......... 11% 25.4M 3s
    ##  15600K .......... .......... .......... .......... .......... 11% 75.2M 3s
    ##  15650K .......... .......... .......... .......... .......... 11% 79.7M 3s
    ##  15700K .......... .......... .......... .......... .......... 11%  101M 3s
    ##  15750K .......... .......... .......... .......... .......... 11% 35.6M 3s
    ##  15800K .......... .......... .......... .......... .......... 11% 24.4M 3s
    ##  15850K .......... .......... .......... .......... .......... 11%  104M 3s
    ##  15900K .......... .......... .......... .......... .......... 11%  120M 3s
    ##  15950K .......... .......... .......... .......... .......... 11%  122M 3s
    ##  16000K .......... .......... .......... .......... .......... 11% 21.2M 3s
    ##  16050K .......... .......... .......... .......... .......... 11% 88.0M 3s
    ##  16100K .......... .......... .......... .......... .......... 11% 96.8M 3s
    ##  16150K .......... .......... .......... .......... .......... 12% 20.7M 3s
    ##  16200K .......... .......... .......... .......... .......... 12% 74.6M 3s
    ##  16250K .......... .......... .......... .......... .......... 12% 86.8M 3s
    ##  16300K .......... .......... .......... .......... .......... 12% 78.6M 3s
    ##  16350K .......... .......... .......... .......... .......... 12%  123M 3s
    ##  16400K .......... .......... .......... .......... .......... 12%  113M 3s
    ##  16450K .......... .......... .......... .......... .......... 12% 28.0M 3s
    ##  16500K .......... .......... .......... .......... .......... 12% 62.2M 3s
    ##  16550K .......... .......... .......... .......... .......... 12% 54.7M 3s
    ##  16600K .......... .......... .......... .......... .......... 12% 83.7M 3s
    ##  16650K .......... .......... .......... .......... .......... 12%  117M 3s
    ##  16700K .......... .......... .......... .......... .......... 12% 35.2M 3s
    ##  16750K .......... .......... .......... .......... .......... 12% 41.0M 3s
    ##  16800K .......... .......... .......... .......... .......... 12% 86.8M 3s
    ##  16850K .......... .......... .......... .......... .......... 12% 63.6M 3s
    ##  16900K .......... .......... .......... .......... .......... 12% 92.4M 3s
    ##  16950K .......... .......... .......... .......... .......... 12%  110M 3s
    ##  17000K .......... .......... .......... .......... .......... 12% 39.4M 3s
    ##  17050K .......... .......... .......... .......... .......... 12% 32.1M 3s
    ##  17100K .......... .......... .......... .......... .......... 12% 81.3M 3s
    ##  17150K .......... .......... .......... .......... .......... 12%  119M 3s
    ##  17200K .......... .......... .......... .......... .......... 12% 76.4M 3s
    ##  17250K .......... .......... .......... .......... .......... 12% 63.9M 3s
    ##  17300K .......... .......... .......... .......... .......... 12%  107M 3s
    ##  17350K .......... .......... .......... .......... .......... 12% 32.7M 3s
    ##  17400K .......... .......... .......... .......... .......... 12% 93.6M 3s
    ##  17450K .......... .......... .......... .......... .......... 12% 61.5M 3s
    ##  17500K .......... .......... .......... .......... .......... 13% 80.7M 3s
    ##  17550K .......... .......... .......... .......... .......... 13% 36.6M 3s
    ##  17600K .......... .......... .......... .......... .......... 13% 82.6M 3s
    ##  17650K .......... .......... .......... .......... .......... 13% 82.1M 3s
    ##  17700K .......... .......... .......... .......... .......... 13%  107M 3s
    ##  17750K .......... .......... .......... .......... .......... 13% 91.3M 3s
    ##  17800K .......... .......... .......... .......... .......... 13% 91.6M 3s
    ##  17850K .......... .......... .......... .......... .......... 13% 34.9M 3s
    ##  17900K .......... .......... .......... .......... .......... 13% 86.0M 3s
    ##  17950K .......... .......... .......... .......... .......... 13% 29.5M 3s
    ##  18000K .......... .......... .......... .......... .......... 13% 99.9M 3s
    ##  18050K .......... .......... .......... .......... .......... 13%  134M 3s
    ##  18100K .......... .......... .......... .......... .......... 13%  103M 3s
    ##  18150K .......... .......... .......... .......... .......... 13% 70.0M 3s
    ##  18200K .......... .......... .......... .......... .......... 13% 66.8M 3s
    ##  18250K .......... .......... .......... .......... .......... 13% 38.1M 3s
    ##  18300K .......... .......... .......... .......... .......... 13% 71.6M 3s
    ##  18350K .......... .......... .......... .......... .......... 13%  106M 3s
    ##  18400K .......... .......... .......... .......... .......... 13% 56.3M 3s
    ##  18450K .......... .......... .......... .......... .......... 13% 72.1M 3s
    ##  18500K .......... .......... .......... .......... .......... 13% 93.0M 3s
    ##  18550K .......... .......... .......... .......... .......... 13% 34.7M 3s
    ##  18600K .......... .......... .......... .......... .......... 13%  116M 3s
    ##  18650K .......... .......... .......... .......... .......... 13% 29.3M 3s
    ##  18700K .......... .......... .......... .......... .......... 13%  119M 3s
    ##  18750K .......... .......... .......... .......... .......... 13%  127M 3s
    ##  18800K .......... .......... .......... .......... .......... 13% 60.6M 3s
    ##  18850K .......... .......... .......... .......... .......... 14% 91.1M 3s
    ##  18900K .......... .......... .......... .......... .......... 14% 35.4M 3s
    ##  18950K .......... .......... .......... .......... .......... 14% 88.2M 3s
    ##  19000K .......... .......... .......... .......... .......... 14%  107M 3s
    ##  19050K .......... .......... .......... .......... .......... 14% 52.4M 3s
    ##  19100K .......... .......... .......... .......... .......... 14% 88.7M 3s
    ##  19150K .......... .......... .......... .......... .......... 14% 33.4M 3s
    ##  19200K .......... .......... .......... .......... .......... 14% 83.3M 3s
    ##  19250K .......... .......... .......... .......... .......... 14%  121M 3s
    ##  19300K .......... .......... .......... .......... .......... 14% 36.1M 3s
    ##  19350K .......... .......... .......... .......... .......... 14%  144M 3s
    ##  19400K .......... .......... .......... .......... .......... 14% 56.4M 3s
    ##  19450K .......... .......... .......... .......... .......... 14% 72.1M 3s
    ##  19500K .......... .......... .......... .......... .......... 14% 97.7M 3s
    ##  19550K .......... .......... .......... .......... .......... 14%  103M 3s
    ##  19600K .......... .......... .......... .......... .......... 14% 47.1M 3s
    ##  19650K .......... .......... .......... .......... .......... 14%  137M 3s
    ##  19700K .......... .......... .......... .......... .......... 14% 37.8M 3s
    ##  19750K .......... .......... .......... .......... .......... 14% 89.5M 3s
    ##  19800K .......... .......... .......... .......... .......... 14% 67.2M 3s
    ##  19850K .......... .......... .......... .......... .......... 14% 87.2M 3s
    ##  19900K .......... .......... .......... .......... .......... 14%  104M 3s
    ##  19950K .......... .......... .......... .......... .......... 14% 40.0M 3s
    ##  20000K .......... .......... .......... .......... .......... 14% 63.4M 3s
    ##  20050K .......... .......... .......... .......... .......... 14%  114M 3s
    ##  20100K .......... .......... .......... .......... .......... 14%  115M 3s
    ##  20150K .......... .......... .......... .......... .......... 14% 87.5M 3s
    ##  20200K .......... .......... .......... .......... .......... 15% 14.0M 3s
    ##  20250K .......... .......... .......... .......... .......... 15%  119M 3s
    ##  20300K .......... .......... .......... .......... .......... 15%  107M 3s
    ##  20350K .......... .......... .......... .......... .......... 15%  130M 3s
    ##  20400K .......... .......... .......... .......... .......... 15%  126M 3s
    ##  20450K .......... .......... .......... .......... .......... 15% 20.5M 3s
    ##  20500K .......... .......... .......... .......... .......... 15% 92.1M 3s
    ##  20550K .......... .......... .......... .......... .......... 15% 79.2M 3s
    ##  20600K .......... .......... .......... .......... .......... 15% 85.2M 3s
    ##  20650K .......... .......... .......... .......... .......... 15%  141M 3s
    ##  20700K .......... .......... .......... .......... .......... 15% 25.4M 3s
    ##  20750K .......... .......... .......... .......... .......... 15% 98.4M 3s
    ##  20800K .......... .......... .......... .......... .......... 15% 55.6M 3s
    ##  20850K .......... .......... .......... .......... .......... 15%  102M 3s
    ##  20900K .......... .......... .......... .......... .......... 15%  122M 3s
    ##  20950K .......... .......... .......... .......... .......... 15% 29.0M 3s
    ##  21000K .......... .......... .......... .......... .......... 15% 69.8M 3s
    ##  21050K .......... .......... .......... .......... .......... 15%  110M 3s
    ##  21100K .......... .......... .......... .......... .......... 15%  100M 3s
    ##  21150K .......... .......... .......... .......... .......... 15%  114M 3s
    ##  21200K .......... .......... .......... .......... .......... 15% 28.9M 3s
    ##  21250K .......... .......... .......... .......... .......... 15% 81.9M 3s
    ##  21300K .......... .......... .......... .......... .......... 15% 87.1M 3s
    ##  21350K .......... .......... .......... .......... .......... 15%  114M 3s
    ##  21400K .......... .......... .......... .......... .......... 15% 68.3M 3s
    ##  21450K .......... .......... .......... .......... .......... 15% 42.0M 3s
    ##  21500K .......... .......... .......... .......... .......... 15%  107M 3s
    ##  21550K .......... .......... .......... .......... .......... 16% 64.2M 3s
    ##  21600K .......... .......... .......... .......... .......... 16% 73.3M 3s
    ##  21650K .......... .......... .......... .......... .......... 16% 99.6M 3s
    ##  21700K .......... .......... .......... .......... .......... 16% 45.9M 3s
    ##  21750K .......... .......... .......... .......... .......... 16% 52.8M 3s
    ##  21800K .......... .......... .......... .......... .......... 16% 30.9M 3s
    ##  21850K .......... .......... .......... .......... .......... 16%  104M 3s
    ##  21900K .......... .......... .......... .......... .......... 16% 82.7M 3s
    ##  21950K .......... .......... .......... .......... .......... 16%  138M 3s
    ##  22000K .......... .......... .......... .......... .......... 16% 19.0M 3s
    ##  22050K .......... .......... .......... .......... .......... 16%  135M 3s
    ##  22100K .......... .......... .......... .......... .......... 16% 75.9M 3s
    ##  22150K .......... .......... .......... .......... .......... 16% 91.9M 3s
    ##  22200K .......... .......... .......... .......... .......... 16% 87.9M 3s
    ##  22250K .......... .......... .......... .......... .......... 16%  127M 3s
    ##  22300K .......... .......... .......... .......... .......... 16% 45.4M 3s
    ##  22350K .......... .......... .......... .......... .......... 16% 30.1M 3s
    ##  22400K .......... .......... .......... .......... .......... 16% 79.1M 3s
    ##  22450K .......... .......... .......... .......... .......... 16%  110M 3s
    ##  22500K .......... .......... .......... .......... .......... 16% 95.5M 3s
    ##  22550K .......... .......... .......... .......... .......... 16% 67.2M 3s
    ##  22600K .......... .......... .......... .......... .......... 16% 33.9M 3s
    ##  22650K .......... .......... .......... .......... .......... 16%  118M 3s
    ##  22700K .......... .......... .......... .......... .......... 16% 75.2M 3s
    ##  22750K .......... .......... .......... .......... .......... 16% 87.1M 3s
    ##  22800K .......... .......... .......... .......... .......... 16%  106M 3s
    ##  22850K .......... .......... .......... .......... .......... 16% 36.0M 3s
    ##  22900K .......... .......... .......... .......... .......... 17% 89.2M 3s
    ##  22950K .......... .......... .......... .......... .......... 17% 80.5M 3s
    ##  23000K .......... .......... .......... .......... .......... 17% 72.1M 3s
    ##  23050K .......... .......... .......... .......... .......... 17% 63.8M 3s
    ##  23100K .......... .......... .......... .......... .......... 17% 17.6M 3s
    ##  23150K .......... .......... .......... .......... .......... 17%  120M 3s
    ##  23200K .......... .......... .......... .......... .......... 17%  139M 3s
    ##  23250K .......... .......... .......... .......... .......... 17% 56.9M 3s
    ##  23300K .......... .......... .......... .......... .......... 17%  101M 3s
    ##  23350K .......... .......... .......... .......... .......... 17% 32.9M 3s
    ##  23400K .......... .......... .......... .......... .......... 17% 98.0M 3s
    ##  23450K .......... .......... .......... .......... .......... 17% 54.2M 3s
    ##  23500K .......... .......... .......... .......... .......... 17% 58.8M 3s
    ##  23550K .......... .......... .......... .......... .......... 17% 82.2M 3s
    ##  23600K .......... .......... .......... .......... .......... 17% 91.7M 2s
    ##  23650K .......... .......... .......... .......... .......... 17% 63.3M 2s
    ##  23700K .......... .......... .......... .......... .......... 17%  133M 2s
    ##  23750K .......... .......... .......... .......... .......... 17% 75.7M 2s
    ##  23800K .......... .......... .......... .......... .......... 17% 62.7M 2s
    ##  23850K .......... .......... .......... .......... .......... 17% 13.7M 2s
    ##  23900K .......... .......... .......... .......... .......... 17% 91.3M 2s
    ##  23950K .......... .......... .......... .......... .......... 17%  153M 2s
    ##  24000K .......... .......... .......... .......... .......... 17%  133M 2s
    ##  24050K .......... .......... .......... .......... .......... 17%  129M 2s
    ##  24100K .......... .......... .......... .......... .......... 17% 20.9M 2s
    ##  24150K .......... .......... .......... .......... .......... 17%  106M 2s
    ##  24200K .......... .......... .......... .......... .......... 17% 37.5M 2s
    ##  24250K .......... .......... .......... .......... .......... 18%  111M 2s
    ##  24300K .......... .......... .......... .......... .......... 18%  134M 2s
    ##  24350K .......... .......... .......... .......... .......... 18% 66.0M 2s
    ##  24400K .......... .......... .......... .......... .......... 18%  134M 2s
    ##  24450K .......... .......... .......... .......... .......... 18% 26.6M 2s
    ##  24500K .......... .......... .......... .......... .......... 18%  107M 2s
    ##  24550K .......... .......... .......... .......... .......... 18%  148M 2s
    ##  24600K .......... .......... .......... .......... .......... 18%  123M 2s
    ##  24650K .......... .......... .......... .......... .......... 18%  114M 2s
    ##  24700K .......... .......... .......... .......... .......... 18% 24.8M 2s
    ##  24750K .......... .......... .......... .......... .......... 18%  128M 2s
    ##  24800K .......... .......... .......... .......... .......... 18% 77.1M 2s
    ##  24850K .......... .......... .......... .......... .......... 18% 99.5M 2s
    ##  24900K .......... .......... .......... .......... .......... 18%  123M 2s
    ##  24950K .......... .......... .......... .......... .......... 18% 39.7M 2s
    ##  25000K .......... .......... .......... .......... .......... 18% 54.8M 2s
    ##  25050K .......... .......... .......... .......... .......... 18% 28.3M 2s
    ##  25100K .......... .......... .......... .......... .......... 18%  122M 2s
    ##  25150K .......... .......... .......... .......... .......... 18%  109M 2s
    ##  25200K .......... .......... .......... .......... .......... 18%  149M 2s
    ##  25250K .......... .......... .......... .......... .......... 18% 54.6M 2s
    ##  25300K .......... .......... .......... .......... .......... 18% 42.4M 2s
    ##  25350K .......... .......... .......... .......... .......... 18% 47.0M 2s
    ##  25400K .......... .......... .......... .......... .......... 18% 62.5M 2s
    ##  25450K .......... .......... .......... .......... .......... 18% 98.7M 2s
    ##  25500K .......... .......... .......... .......... .......... 18%  106M 2s
    ##  25550K .......... .......... .......... .......... .......... 18% 35.1M 2s
    ##  25600K .......... .......... .......... .......... .......... 19% 67.7M 2s
    ##  25650K .......... .......... .......... .......... .......... 19%  163M 2s
    ##  25700K .......... .......... .......... .......... .......... 19%  142M 2s
    ##  25750K .......... .......... .......... .......... .......... 19%  101M 2s
    ##  25800K .......... .......... .......... .......... .......... 19%  120M 2s
    ##  25850K .......... .......... .......... .......... .......... 19% 27.8M 2s
    ##  25900K .......... .......... .......... .......... .......... 19%  129M 2s
    ##  25950K .......... .......... .......... .......... .......... 19%  120M 2s
    ##  26000K .......... .......... .......... .......... .......... 19% 92.1M 2s
    ##  26050K .......... .......... .......... .......... .......... 19% 27.2M 2s
    ##  26100K .......... .......... .......... .......... .......... 19% 76.6M 2s
    ##  26150K .......... .......... .......... .......... .......... 19% 38.2M 2s
    ##  26200K .......... .......... .......... .......... .......... 19% 35.8M 2s
    ##  26250K .......... .......... .......... .......... .......... 19%  153M 2s
    ##  26300K .......... .......... .......... .......... .......... 19%  136M 2s
    ##  26350K .......... .......... .......... .......... .......... 19% 96.7M 2s
    ##  26400K .......... .......... .......... .......... .......... 19% 56.3M 2s
    ##  26450K .......... .......... .......... .......... .......... 19%  111M 2s
    ##  26500K .......... .......... .......... .......... .......... 19% 55.9M 2s
    ##  26550K .......... .......... .......... .......... .......... 19%  107M 2s
    ##  26600K .......... .......... .......... .......... .......... 19% 48.6M 2s
    ##  26650K .......... .......... .......... .......... .......... 19% 56.3M 2s
    ##  26700K .......... .......... .......... .......... .......... 19% 46.2M 2s
    ##  26750K .......... .......... .......... .......... .......... 19%  160M 2s
    ##  26800K .......... .......... .......... .......... .......... 19%  134M 2s
    ##  26850K .......... .......... .......... .......... .......... 19% 60.0M 2s
    ##  26900K .......... .......... .......... .......... .......... 20% 54.3M 2s
    ##  26950K .......... .......... .......... .......... .......... 20% 26.3M 2s
    ##  27000K .......... .......... .......... .......... .......... 20%  115M 2s
    ##  27050K .......... .......... .......... .......... .......... 20%  184M 2s
    ##  27100K .......... .......... .......... .......... .......... 20% 45.1M 2s
    ##  27150K .......... .......... .......... .......... .......... 20%  171M 2s
    ##  27200K .......... .......... .......... .......... .......... 20% 46.5M 2s
    ##  27250K .......... .......... .......... .......... .......... 20% 64.5M 2s
    ##  27300K .......... .......... .......... .......... .......... 20%  117M 2s
    ##  27350K .......... .......... .......... .......... .......... 20% 44.1M 2s
    ##  27400K .......... .......... .......... .......... .......... 20% 89.8M 2s
    ##  27450K .......... .......... .......... .......... .......... 20% 30.6M 2s
    ##  27500K .......... .......... .......... .......... .......... 20% 83.7M 2s
    ##  27550K .......... .......... .......... .......... .......... 20%  189M 2s
    ##  27600K .......... .......... .......... .......... .......... 20% 62.0M 2s
    ##  27650K .......... .......... .......... .......... .......... 20%  116M 2s
    ##  27700K .......... .......... .......... .......... .......... 20% 38.2M 2s
    ##  27750K .......... .......... .......... .......... .......... 20%  125M 2s
    ##  27800K .......... .......... .......... .......... .......... 20% 72.4M 2s
    ##  27850K .......... .......... .......... .......... .......... 20% 34.6M 2s
    ##  27900K .......... .......... .......... .......... .......... 20%  106M 2s
    ##  27950K .......... .......... .......... .......... .......... 20%  134M 2s
    ##  28000K .......... .......... .......... .......... .......... 20% 50.7M 2s
    ##  28050K .......... .......... .......... .......... .......... 20%  140M 2s
    ##  28100K .......... .......... .......... .......... .......... 20% 33.6M 2s
    ##  28150K .......... .......... .......... .......... .......... 20% 57.0M 2s
    ##  28200K .......... .......... .......... .......... .......... 20%  140M 2s
    ##  28250K .......... .......... .......... .......... .......... 21% 60.9M 2s
    ##  28300K .......... .......... .......... .......... .......... 21% 68.7M 2s
    ##  28350K .......... .......... .......... .......... .......... 21%  173M 2s
    ##  28400K .......... .......... .......... .......... .......... 21% 33.8M 2s
    ##  28450K .......... .......... .......... .......... .......... 21%  107M 2s
    ##  28500K .......... .......... .......... .......... .......... 21% 56.7M 2s
    ##  28550K .......... .......... .......... .......... .......... 21% 50.0M 2s
    ##  28600K .......... .......... .......... .......... .......... 21%  142M 2s
    ##  28650K .......... .......... .......... .......... .......... 21%  167M 2s
    ##  28700K .......... .......... .......... .......... .......... 21% 41.2M 2s
    ##  28750K .......... .......... .......... .......... .......... 21% 35.2M 2s
    ##  28800K .......... .......... .......... .......... .......... 21% 23.6M 2s
    ##  28850K .......... .......... .......... .......... .......... 21% 44.7M 2s
    ##  28900K .......... .......... .......... .......... .......... 21%  134M 2s
    ##  28950K .......... .......... .......... .......... .......... 21%  193M 2s
    ##  29000K .......... .......... .......... .......... .......... 21%  167M 2s
    ##  29050K .......... .......... .......... .......... .......... 21% 22.0M 2s
    ##  29100K .......... .......... .......... .......... .......... 21% 58.6M 2s
    ##  29150K .......... .......... .......... .......... .......... 21% 61.3M 2s
    ##  29200K .......... .......... .......... .......... .......... 21%  157M 2s
    ##  29250K .......... .......... .......... .......... .......... 21%  169M 2s
    ##  29300K .......... .......... .......... .......... .......... 21% 59.6M 2s
    ##  29350K .......... .......... .......... .......... .......... 21%  115M 2s
    ##  29400K .......... .......... .......... .......... .......... 21% 24.4M 2s
    ##  29450K .......... .......... .......... .......... .......... 21%  149M 2s
    ##  29500K .......... .......... .......... .......... .......... 21% 83.8M 2s
    ##  29550K .......... .......... .......... .......... .......... 21% 90.0M 2s
    ##  29600K .......... .......... .......... .......... .......... 22%  138M 2s
    ##  29650K .......... .......... .......... .......... .......... 22% 83.1M 2s
    ##  29700K .......... .......... .......... .......... .......... 22% 39.5M 2s
    ##  29750K .......... .......... .......... .......... .......... 22% 51.9M 2s
    ##  29800K .......... .......... .......... .......... .......... 22%  162M 2s
    ##  29850K .......... .......... .......... .......... .......... 22% 36.2M 2s
    ##  29900K .......... .......... .......... .......... .......... 22%  141M 2s
    ##  29950K .......... .......... .......... .......... .......... 22% 56.3M 2s
    ##  30000K .......... .......... .......... .......... .......... 22%  134M 2s
    ##  30050K .......... .......... .......... .......... .......... 22%  124M 2s
    ##  30100K .......... .......... .......... .......... .......... 22% 28.9M 2s
    ##  30150K .......... .......... .......... .......... .......... 22% 68.5M 2s
    ##  30200K .......... .......... .......... .......... .......... 22% 48.6M 2s
    ##  30250K .......... .......... .......... .......... .......... 22%  172M 2s
    ##  30300K .......... .......... .......... .......... .......... 22%  163M 2s
    ##  30350K .......... .......... .......... .......... .......... 22% 24.9M 2s
    ##  30400K .......... .......... .......... .......... .......... 22%  130M 2s
    ##  30450K .......... .......... .......... .......... .......... 22% 45.6M 2s
    ##  30500K .......... .......... .......... .......... .......... 22%  101M 2s
    ##  30550K .......... .......... .......... .......... .......... 22%  157M 2s
    ##  30600K .......... .......... .......... .......... .......... 22% 85.6M 2s
    ##  30650K .......... .......... .......... .......... .......... 22% 49.3M 2s
    ##  30700K .......... .......... .......... .......... .......... 22% 40.9M 2s
    ##  30750K .......... .......... .......... .......... .......... 22% 76.6M 2s
    ##  30800K .......... .......... .......... .......... .......... 22% 84.0M 2s
    ##  30850K .......... .......... .......... .......... .......... 22%  107M 2s
    ##  30900K .......... .......... .......... .......... .......... 22%  141M 2s
    ##  30950K .......... .......... .......... .......... .......... 23% 44.9M 2s
    ##  31000K .......... .......... .......... .......... .......... 23%  141M 2s
    ##  31050K .......... .......... .......... .......... .......... 23% 47.0M 2s
    ##  31100K .......... .......... .......... .......... .......... 23% 70.2M 2s
    ##  31150K .......... .......... .......... .......... .......... 23%  149M 2s
    ##  31200K .......... .......... .......... .......... .......... 23% 41.3M 2s
    ##  31250K .......... .......... .......... .......... .......... 23%  138M 2s
    ##  31300K .......... .......... .......... .......... .......... 23% 50.9M 2s
    ##  31350K .......... .......... .......... .......... .......... 23% 65.7M 2s
    ##  31400K .......... .......... .......... .......... .......... 23%  119M 2s
    ##  31450K .......... .......... .......... .......... .......... 23%  149M 2s
    ##  31500K .......... .......... .......... .......... .......... 23% 56.3M 2s
    ##  31550K .......... .......... .......... .......... .......... 23% 51.3M 2s
    ##  31600K .......... .......... .......... .......... .......... 23% 53.7M 2s
    ##  31650K .......... .......... .......... .......... .......... 23% 62.9M 2s
    ##  31700K .......... .......... .......... .......... .......... 23%  109M 2s
    ##  31750K .......... .......... .......... .......... .......... 23% 61.2M 2s
    ##  31800K .......... .......... .......... .......... .......... 23% 55.0M 2s
    ##  31850K .......... .......... .......... .......... .......... 23% 25.1M 2s
    ##  31900K .......... .......... .......... .......... .......... 23% 45.7M 2s
    ##  31950K .......... .......... .......... .......... .......... 23%  123M 2s
    ##  32000K .......... .......... .......... .......... .......... 23%  141M 2s
    ##  32050K .......... .......... .......... .......... .......... 23%  137M 2s
    ##  32100K .......... .......... .......... .......... .......... 23% 10.8M 2s
    ##  32150K .......... .......... .......... .......... .......... 23% 92.5M 2s
    ##  32200K .......... .......... .......... .......... .......... 23% 65.0M 2s
    ##  32250K .......... .......... .......... .......... .......... 23%  111M 2s
    ##  32300K .......... .......... .......... .......... .......... 24%  143M 2s
    ##  32350K .......... .......... .......... .......... .......... 24% 32.5M 2s
    ##  32400K .......... .......... .......... .......... .......... 24%  161M 2s
    ##  32450K .......... .......... .......... .......... .......... 24%  143M 2s
    ##  32500K .......... .......... .......... .......... .......... 24% 32.5M 2s
    ##  32550K .......... .......... .......... .......... .......... 24%  162M 2s
    ##  32600K .......... .......... .......... .......... .......... 24% 42.1M 2s
    ##  32650K .......... .......... .......... .......... .......... 24% 91.4M 2s
    ##  32700K .......... .......... .......... .......... .......... 24%  163M 2s
    ##  32750K .......... .......... .......... .......... .......... 24% 46.3M 2s
    ##  32800K .......... .......... .......... .......... .......... 24%  116M 2s
    ##  32850K .......... .......... .......... .......... .......... 24% 38.8M 2s
    ##  32900K .......... .......... .......... .......... .......... 24%  112M 2s
    ##  32950K .......... .......... .......... .......... .......... 24%  127M 2s
    ##  33000K .......... .......... .......... .......... .......... 24% 26.0M 2s
    ##  33050K .......... .......... .......... .......... .......... 24%  177M 2s
    ##  33100K .......... .......... .......... .......... .......... 24%  103M 2s
    ##  33150K .......... .......... .......... .......... .......... 24%  136M 2s
    ##  33200K .......... .......... .......... .......... .......... 24% 90.5M 2s
    ##  33250K .......... .......... .......... .......... .......... 24% 29.5M 2s
    ##  33300K .......... .......... .......... .......... .......... 24%  130M 2s
    ##  33350K .......... .......... .......... .......... .......... 24%  168M 2s
    ##  33400K .......... .......... .......... .......... .......... 24% 42.5M 2s
    ##  33450K .......... .......... .......... .......... .......... 24%  183M 2s
    ##  33500K .......... .......... .......... .......... .......... 24% 31.6M 2s
    ##  33550K .......... .......... .......... .......... .......... 24%  129M 2s
    ##  33600K .......... .......... .......... .......... .......... 24%  122M 2s
    ##  33650K .......... .......... .......... .......... .......... 25%  110M 2s
    ##  33700K .......... .......... .......... .......... .......... 25% 82.7M 2s
    ##  33750K .......... .......... .......... .......... .......... 25% 26.3M 2s
    ##  33800K .......... .......... .......... .......... .......... 25% 66.0M 2s
    ##  33850K .......... .......... .......... .......... .......... 25% 74.3M 2s
    ##  33900K .......... .......... .......... .......... .......... 25% 97.3M 2s
    ##  33950K .......... .......... .......... .......... .......... 25% 85.4M 2s
    ##  34000K .......... .......... .......... .......... .......... 25% 15.8M 2s
    ##  34050K .......... .......... .......... .......... .......... 25% 83.0M 2s
    ##  34100K .......... .......... .......... .......... .......... 25% 17.1M 2s
    ##  34150K .......... .......... .......... .......... .......... 25% 67.9M 2s
    ##  34200K .......... .......... .......... .......... .......... 25% 70.4M 2s
    ##  34250K .......... .......... .......... .......... .......... 25% 79.5M 2s
    ##  34300K .......... .......... .......... .......... .......... 25% 93.9M 2s
    ##  34350K .......... .......... .......... .......... .......... 25%  117M 2s
    ##  34400K .......... .......... .......... .......... .......... 25% 37.3M 2s
    ##  34450K .......... .......... .......... .......... .......... 25% 87.4M 2s
    ##  34500K .......... .......... .......... .......... .......... 25%  114M 2s
    ##  34550K .......... .......... .......... .......... .......... 25% 85.2M 2s
    ##  34600K .......... .......... .......... .......... .......... 25% 95.0M 2s
    ##  34650K .......... .......... .......... .......... .......... 25%  119M 2s
    ##  34700K .......... .......... .......... .......... .......... 25% 44.2M 2s
    ##  34750K .......... .......... .......... .......... .......... 25% 57.9M 2s
    ##  34800K .......... .......... .......... .......... .......... 25% 51.5M 2s
    ##  34850K .......... .......... .......... .......... .......... 25% 82.3M 2s
    ##  34900K .......... .......... .......... .......... .......... 25% 81.1M 2s
    ##  34950K .......... .......... .......... .......... .......... 25% 75.0M 2s
    ##  35000K .......... .......... .......... .......... .......... 26% 32.9M 2s
    ##  35050K .......... .......... .......... .......... .......... 26% 76.3M 2s
    ##  35100K .......... .......... .......... .......... .......... 26% 75.3M 2s
    ##  35150K .......... .......... .......... .......... .......... 26% 95.6M 2s
    ##  35200K .......... .......... .......... .......... .......... 26% 87.2M 2s
    ##  35250K .......... .......... .......... .......... .......... 26% 62.1M 2s
    ##  35300K .......... .......... .......... .......... .......... 26% 77.9M 2s
    ##  35350K .......... .......... .......... .......... .......... 26% 81.7M 2s
    ##  35400K .......... .......... .......... .......... .......... 26% 60.0M 2s
    ##  35450K .......... .......... .......... .......... .......... 26% 95.3M 2s
    ##  35500K .......... .......... .......... .......... .......... 26% 79.5M 2s
    ##  35550K .......... .......... .......... .......... .......... 26% 85.8M 2s
    ##  35600K .......... .......... .......... .......... .......... 26% 63.4M 2s
    ##  35650K .......... .......... .......... .......... .......... 26% 62.9M 2s
    ##  35700K .......... .......... .......... .......... .......... 26% 57.8M 2s
    ##  35750K .......... .......... .......... .......... .......... 26%  104M 2s
    ##  35800K .......... .......... .......... .......... .......... 26% 48.4M 2s
    ##  35850K .......... .......... .......... .......... .......... 26% 32.0M 2s
    ##  35900K .......... .......... .......... .......... .......... 26% 79.6M 2s
    ##  35950K .......... .......... .......... .......... .......... 26% 61.2M 2s
    ##  36000K .......... .......... .......... .......... .......... 26% 80.7M 2s
    ##  36050K .......... .......... .......... .......... .......... 26%  123M 2s
    ##  36100K .......... .......... .......... .......... .......... 26% 45.5M 2s
    ##  36150K .......... .......... .......... .......... .......... 26% 73.7M 2s
    ##  36200K .......... .......... .......... .......... .......... 26% 75.2M 2s
    ##  36250K .......... .......... .......... .......... .......... 26% 51.1M 2s
    ##  36300K .......... .......... .......... .......... .......... 26% 69.3M 2s
    ##  36350K .......... .......... .......... .......... .......... 27%  115M 2s
    ##  36400K .......... .......... .......... .......... .......... 27% 55.3M 2s
    ##  36450K .......... .......... .......... .......... .......... 27% 90.4M 2s
    ##  36500K .......... .......... .......... .......... .......... 27% 62.5M 2s
    ##  36550K .......... .......... .......... .......... .......... 27% 58.6M 2s
    ##  36600K .......... .......... .......... .......... .......... 27% 66.2M 2s
    ##  36650K .......... .......... .......... .......... .......... 27%  130M 2s
    ##  36700K .......... .......... .......... .......... .......... 27% 71.0M 2s
    ##  36750K .......... .......... .......... .......... .......... 27% 74.9M 2s
    ##  36800K .......... .......... .......... .......... .......... 27% 30.7M 2s
    ##  36850K .......... .......... .......... .......... .......... 27% 89.7M 2s
    ##  36900K .......... .......... .......... .......... .......... 27% 86.9M 2s
    ##  36950K .......... .......... .......... .......... .......... 27%  124M 2s
    ##  37000K .......... .......... .......... .......... .......... 27% 93.6M 2s
    ##  37050K .......... .......... .......... .......... .......... 27% 18.5M 2s
    ##  37100K .......... .......... .......... .......... .......... 27% 90.8M 2s
    ##  37150K .......... .......... .......... .......... .......... 27% 86.3M 2s
    ##  37200K .......... .......... .......... .......... .......... 27% 90.7M 2s
    ##  37250K .......... .......... .......... .......... .......... 27% 99.3M 2s
    ##  37300K .......... .......... .......... .......... .......... 27% 44.3M 2s
    ##  37350K .......... .......... .......... .......... .......... 27% 40.4M 2s
    ##  37400K .......... .......... .......... .......... .......... 27%  106M 2s
    ##  37450K .......... .......... .......... .......... .......... 27% 60.1M 2s
    ##  37500K .......... .......... .......... .......... .......... 27%  103M 2s
    ##  37550K .......... .......... .......... .......... .......... 27% 70.3M 2s
    ##  37600K .......... .......... .......... .......... .......... 27% 84.1M 2s
    ##  37650K .......... .......... .......... .......... .......... 27% 52.8M 2s
    ##  37700K .......... .......... .......... .......... .......... 28% 47.6M 2s
    ##  37750K .......... .......... .......... .......... .......... 28% 54.1M 2s
    ##  37800K .......... .......... .......... .......... .......... 28% 68.4M 2s
    ##  37850K .......... .......... .......... .......... .......... 28%  112M 2s
    ##  37900K .......... .......... .......... .......... .......... 28%  119M 2s
    ##  37950K .......... .......... .......... .......... .......... 28% 50.3M 2s
    ##  38000K .......... .......... .......... .......... .......... 28% 53.3M 2s
    ##  38050K .......... .......... .......... .......... .......... 28% 61.6M 2s
    ##  38100K .......... .......... .......... .......... .......... 28% 86.2M 2s
    ##  38150K .......... .......... .......... .......... .......... 28% 84.0M 2s
    ##  38200K .......... .......... .......... .......... .......... 28%  112M 2s
    ##  38250K .......... .......... .......... .......... .......... 28% 50.3M 2s
    ##  38300K .......... .......... .......... .......... .......... 28% 91.4M 2s
    ##  38350K .......... .......... .......... .......... .......... 28% 60.1M 2s
    ##  38400K .......... .......... .......... .......... .......... 28% 88.1M 2s
    ##  38450K .......... .......... .......... .......... .......... 28% 96.6M 2s
    ##  38500K .......... .......... .......... .......... .......... 28% 32.5M 2s
    ##  38550K .......... .......... .......... .......... .......... 28%  122M 2s
    ##  38600K .......... .......... .......... .......... .......... 28% 57.7M 2s
    ##  38650K .......... .......... .......... .......... .......... 28% 84.5M 2s
    ##  38700K .......... .......... .......... .......... .......... 28%  104M 2s
    ##  38750K .......... .......... .......... .......... .......... 28% 70.7M 2s
    ##  38800K .......... .......... .......... .......... .......... 28% 48.1M 2s
    ##  38850K .......... .......... .......... .......... .......... 28% 42.0M 2s
    ##  38900K .......... .......... .......... .......... .......... 28% 80.3M 2s
    ##  38950K .......... .......... .......... .......... .......... 28% 92.1M 2s
    ##  39000K .......... .......... .......... .......... .......... 28% 96.8M 2s
    ##  39050K .......... .......... .......... .......... .......... 29%  110M 2s
    ##  39100K .......... .......... .......... .......... .......... 29% 76.6M 2s
    ##  39150K .......... .......... .......... .......... .......... 29% 75.4M 2s
    ##  39200K .......... .......... .......... .......... .......... 29% 92.6M 2s
    ##  39250K .......... .......... .......... .......... .......... 29% 25.1M 2s
    ##  39300K .......... .......... .......... .......... .......... 29%  104M 2s
    ##  39350K .......... .......... .......... .......... .......... 29% 99.0M 2s
    ##  39400K .......... .......... .......... .......... .......... 29%  106M 2s
    ##  39450K .......... .......... .......... .......... .......... 29% 29.5M 2s
    ##  39500K .......... .......... .......... .......... .......... 29% 43.1M 2s
    ##  39550K .......... .......... .......... .......... .......... 29% 78.3M 2s
    ##  39600K .......... .......... .......... .......... .......... 29% 69.0M 2s
    ##  39650K .......... .......... .......... .......... .......... 29%  127M 2s
    ##  39700K .......... .......... .......... .......... .......... 29% 85.0M 2s
    ##  39750K .......... .......... .......... .......... .......... 29%  135M 2s
    ##  39800K .......... .......... .......... .......... .......... 29% 40.4M 2s
    ##  39850K .......... .......... .......... .......... .......... 29% 31.3M 2s
    ##  39900K .......... .......... .......... .......... .......... 29% 94.5M 2s
    ##  39950K .......... .......... .......... .......... .......... 29% 88.7M 2s
    ##  40000K .......... .......... .......... .......... .......... 29% 90.6M 2s
    ##  40050K .......... .......... .......... .......... .......... 29%  122M 2s
    ##  40100K .......... .......... .......... .......... .......... 29% 99.9M 2s
    ##  40150K .......... .......... .......... .......... .......... 29% 58.1M 2s
    ##  40200K .......... .......... .......... .......... .......... 29%  110M 2s
    ##  40250K .......... .......... .......... .......... .......... 29% 82.0M 2s
    ##  40300K .......... .......... .......... .......... .......... 29% 44.1M 2s
    ##  40350K .......... .......... .......... .......... .......... 29%  118M 2s
    ##  40400K .......... .......... .......... .......... .......... 30%  108M 2s
    ##  40450K .......... .......... .......... .......... .......... 30%  101M 2s
    ##  40500K .......... .......... .......... .......... .......... 30% 43.0M 2s
    ##  40550K .......... .......... .......... .......... .......... 30%  125M 2s
    ##  40600K .......... .......... .......... .......... .......... 30% 75.0M 2s
    ##  40650K .......... .......... .......... .......... .......... 30% 25.3M 2s
    ##  40700K .......... .......... .......... .......... .......... 30% 86.9M 2s
    ##  40750K .......... .......... .......... .......... .......... 30% 93.3M 2s
    ##  40800K .......... .......... .......... .......... .......... 30%  119M 2s
    ##  40850K .......... .......... .......... .......... .......... 30%  104M 2s
    ##  40900K .......... .......... .......... .......... .......... 30%  101M 2s
    ##  40950K .......... .......... .......... .......... .......... 30%  139M 2s
    ##  41000K .......... .......... .......... .......... .......... 30% 41.0M 2s
    ##  41050K .......... .......... .......... .......... .......... 30% 63.8M 2s
    ##  41100K .......... .......... .......... .......... .......... 30% 33.3M 2s
    ##  41150K .......... .......... .......... .......... .......... 30%  133M 2s
    ##  41200K .......... .......... .......... .......... .......... 30%  113M 2s
    ##  41250K .......... .......... .......... .......... .......... 30%  124M 2s
    ##  41300K .......... .......... .......... .......... .......... 30% 54.0M 2s
    ##  41350K .......... .......... .......... .......... .......... 30% 56.2M 2s
    ##  41400K .......... .......... .......... .......... .......... 30% 79.9M 2s
    ##  41450K .......... .......... .......... .......... .......... 30%  116M 2s
    ##  41500K .......... .......... .......... .......... .......... 30% 59.1M 2s
    ##  41550K .......... .......... .......... .......... .......... 30% 99.1M 2s
    ##  41600K .......... .......... .......... .......... .......... 30% 96.8M 2s
    ##  41650K .......... .......... .......... .......... .......... 30%  110M 2s
    ##  41700K .......... .......... .......... .......... .......... 30% 76.1M 2s
    ##  41750K .......... .......... .......... .......... .......... 31% 69.9M 2s
    ##  41800K .......... .......... .......... .......... .......... 31% 88.3M 2s
    ##  41850K .......... .......... .......... .......... .......... 31% 89.7M 2s
    ##  41900K .......... .......... .......... .......... .......... 31% 65.4M 2s
    ##  41950K .......... .......... .......... .......... .......... 31% 92.4M 2s
    ##  42000K .......... .......... .......... .......... .......... 31% 79.8M 2s
    ##  42050K .......... .......... .......... .......... .......... 31% 80.1M 2s
    ##  42100K .......... .......... .......... .......... .......... 31% 96.3M 2s
    ##  42150K .......... .......... .......... .......... .......... 31% 78.1M 2s
    ##  42200K .......... .......... .......... .......... .......... 31% 50.3M 2s
    ##  42250K .......... .......... .......... .......... .......... 31%  114M 2s
    ##  42300K .......... .......... .......... .......... .......... 31% 95.7M 2s
    ##  42350K .......... .......... .......... .......... .......... 31%  110M 2s
    ##  42400K .......... .......... .......... .......... .......... 31% 98.0M 2s
    ##  42450K .......... .......... .......... .......... .......... 31% 47.1M 2s
    ##  42500K .......... .......... .......... .......... .......... 31% 36.7M 2s
    ##  42550K .......... .......... .......... .......... .......... 31%  111M 2s
    ##  42600K .......... .......... .......... .......... .......... 31%  101M 2s
    ##  42650K .......... .......... .......... .......... .......... 31% 86.0M 2s
    ##  42700K .......... .......... .......... .......... .......... 31% 91.4M 2s
    ##  42750K .......... .......... .......... .......... .......... 31% 23.6M 2s
    ##  42800K .......... .......... .......... .......... .......... 31% 92.6M 2s
    ##  42850K .......... .......... .......... .......... .......... 31%  125M 2s
    ##  42900K .......... .......... .......... .......... .......... 31% 77.7M 2s
    ##  42950K .......... .......... .......... .......... .......... 31% 40.5M 2s
    ##  43000K .......... .......... .......... .......... .......... 31% 76.8M 2s
    ##  43050K .......... .......... .......... .......... .......... 31%  114M 2s
    ##  43100K .......... .......... .......... .......... .......... 32% 38.4M 2s
    ##  43150K .......... .......... .......... .......... .......... 32%  114M 2s
    ##  43200K .......... .......... .......... .......... .......... 32%  117M 2s
    ##  43250K .......... .......... .......... .......... .......... 32%  118M 2s
    ##  43300K .......... .......... .......... .......... .......... 32% 96.6M 2s
    ##  43350K .......... .......... .......... .......... .......... 32% 71.1M 2s
    ##  43400K .......... .......... .......... .......... .......... 32% 38.2M 2s
    ##  43450K .......... .......... .......... .......... .......... 32% 61.3M 2s
    ##  43500K .......... .......... .......... .......... .......... 32%  116M 2s
    ##  43550K .......... .......... .......... .......... .......... 32% 91.1M 2s
    ##  43600K .......... .......... .......... .......... .......... 32% 57.5M 2s
    ##  43650K .......... .......... .......... .......... .......... 32%  116M 2s
    ##  43700K .......... .......... .......... .......... .......... 32%  118M 2s
    ##  43750K .......... .......... .......... .......... .......... 32% 64.0M 2s
    ##  43800K .......... .......... .......... .......... .......... 32%  102M 2s
    ##  43850K .......... .......... .......... .......... .......... 32% 80.0M 2s
    ##  43900K .......... .......... .......... .......... .......... 32% 41.5M 2s
    ##  43950K .......... .......... .......... .......... .......... 32% 99.4M 2s
    ##  44000K .......... .......... .......... .......... .......... 32%  120M 2s
    ##  44050K .......... .......... .......... .......... .......... 32%  101M 2s
    ##  44100K .......... .......... .......... .......... .......... 32% 72.4M 2s
    ##  44150K .......... .......... .......... .......... .......... 32% 81.3M 2s
    ##  44200K .......... .......... .......... .......... .......... 32% 99.1M 2s
    ##  44250K .......... .......... .......... .......... .......... 32%  112M 2s
    ##  44300K .......... .......... .......... .......... .......... 32% 67.4M 2s
    ##  44350K .......... .......... .......... .......... .......... 32% 89.8M 2s
    ##  44400K .......... .......... .......... .......... .......... 32%  111M 2s
    ##  44450K .......... .......... .......... .......... .......... 33%  130M 2s
    ##  44500K .......... .......... .......... .......... .......... 33% 34.7M 2s
    ##  44550K .......... .......... .......... .......... .......... 33% 34.5M 2s
    ##  44600K .......... .......... .......... .......... .......... 33% 92.1M 2s
    ##  44650K .......... .......... .......... .......... .......... 33%  116M 2s
    ##  44700K .......... .......... .......... .......... .......... 33%  123M 2s
    ##  44750K .......... .......... .......... .......... .......... 33%  139M 2s
    ##  44800K .......... .......... .......... .......... .......... 33% 22.1M 2s
    ##  44850K .......... .......... .......... .......... .......... 33% 39.5M 2s
    ##  44900K .......... .......... .......... .......... .......... 33% 57.1M 2s
    ##  44950K .......... .......... .......... .......... .......... 33% 92.8M 2s
    ##  45000K .......... .......... .......... .......... .......... 33% 96.0M 2s
    ##  45050K .......... .......... .......... .......... .......... 33%  139M 2s
    ##  45100K .......... .......... .......... .......... .......... 33% 57.7M 2s
    ##  45150K .......... .......... .......... .......... .......... 33%  116M 2s
    ##  45200K .......... .......... .......... .......... .......... 33% 41.1M 2s
    ##  45250K .......... .......... .......... .......... .......... 33% 51.2M 2s
    ##  45300K .......... .......... .......... .......... .......... 33%  102M 2s
    ##  45350K .......... .......... .......... .......... .......... 33% 53.0M 2s
    ##  45400K .......... .......... .......... .......... .......... 33%  115M 2s
    ##  45450K .......... .......... .......... .......... .......... 33%  135M 2s
    ##  45500K .......... .......... .......... .......... .......... 33%  115M 2s
    ##  45550K .......... .......... .......... .......... .......... 33%  129M 2s
    ##  45600K .......... .......... .......... .......... .......... 33% 54.4M 2s
    ##  45650K .......... .......... .......... .......... .......... 33% 97.0M 2s
    ##  45700K .......... .......... .......... .......... .......... 33% 67.0M 2s
    ##  45750K .......... .......... .......... .......... .......... 33%  124M 2s
    ##  45800K .......... .......... .......... .......... .......... 34% 91.9M 2s
    ##  45850K .......... .......... .......... .......... .......... 34%  130M 2s
    ##  45900K .......... .......... .......... .......... .......... 34% 40.5M 2s
    ##  45950K .......... .......... .......... .......... .......... 34% 98.2M 2s
    ##  46000K .......... .......... .......... .......... .......... 34%  121M 2s
    ##  46050K .......... .......... .......... .......... .......... 34%  102M 2s
    ##  46100K .......... .......... .......... .......... .......... 34%  117M 2s
    ##  46150K .......... .......... .......... .......... .......... 34% 80.3M 2s
    ##  46200K .......... .......... .......... .......... .......... 34% 26.6M 2s
    ##  46250K .......... .......... .......... .......... .......... 34% 23.9M 2s
    ##  46300K .......... .......... .......... .......... .......... 34%  113M 2s
    ##  46350K .......... .......... .......... .......... .......... 34%  155M 2s
    ##  46400K .......... .......... .......... .......... .......... 34%  119M 2s
    ##  46450K .......... .......... .......... .......... .......... 34%  145M 2s
    ##  46500K .......... .......... .......... .......... .......... 34% 20.2M 2s
    ##  46550K .......... .......... .......... .......... .......... 34% 98.8M 2s
    ##  46600K .......... .......... .......... .......... .......... 34% 34.0M 2s
    ##  46650K .......... .......... .......... .......... .......... 34% 60.3M 2s
    ##  46700K .......... .......... .......... .......... .......... 34% 95.3M 2s
    ##  46750K .......... .......... .......... .......... .......... 34% 70.9M 2s
    ##  46800K .......... .......... .......... .......... .......... 34%  114M 2s
    ##  46850K .......... .......... .......... .......... .......... 34%  145M 2s
    ##  46900K .......... .......... .......... .......... .......... 34% 69.9M 2s
    ##  46950K .......... .......... .......... .......... .......... 34% 77.6M 2s
    ##  47000K .......... .......... .......... .......... .......... 34% 84.8M 2s
    ##  47050K .......... .......... .......... .......... .......... 34% 77.8M 2s
    ##  47100K .......... .......... .......... .......... .......... 34%  122M 2s
    ##  47150K .......... .......... .......... .......... .......... 35%  118M 2s
    ##  47200K .......... .......... .......... .......... .......... 35% 44.0M 2s
    ##  47250K .......... .......... .......... .......... .......... 35%  112M 2s
    ##  47300K .......... .......... .......... .......... .......... 35% 92.5M 2s
    ##  47350K .......... .......... .......... .......... .......... 35% 55.5M 2s
    ##  47400K .......... .......... .......... .......... .......... 35% 86.4M 2s
    ##  47450K .......... .......... .......... .......... .......... 35% 70.2M 2s
    ##  47500K .......... .......... .......... .......... .......... 35%  105M 2s
    ##  47550K .......... .......... .......... .......... .......... 35% 30.7M 2s
    ##  47600K .......... .......... .......... .......... .......... 35%  107M 2s
    ##  47650K .......... .......... .......... .......... .......... 35%  161M 2s
    ##  47700K .......... .......... .......... .......... .......... 35%  133M 2s
    ##  47750K .......... .......... .......... .......... .......... 35%  155M 2s
    ##  47800K .......... .......... .......... .......... .......... 35% 25.2M 2s
    ##  47850K .......... .......... .......... .......... .......... 35%  113M 2s
    ##  47900K .......... .......... .......... .......... .......... 35%  122M 2s
    ##  47950K .......... .......... .......... .......... .......... 35%  125M 2s
    ##  48000K .......... .......... .......... .......... .......... 35% 99.6M 2s
    ##  48050K .......... .......... .......... .......... .......... 35%  126M 2s
    ##  48100K .......... .......... .......... .......... .......... 35% 31.0M 2s
    ##  48150K .......... .......... .......... .......... .......... 35% 94.0M 2s
    ##  48200K .......... .......... .......... .......... .......... 35% 94.3M 2s
    ##  48250K .......... .......... .......... .......... .......... 35%  122M 2s
    ##  48300K .......... .......... .......... .......... .......... 35% 43.8M 2s
    ##  48350K .......... .......... .......... .......... .......... 35% 44.7M 2s
    ##  48400K .......... .......... .......... .......... .......... 35%  145M 2s
    ##  48450K .......... .......... .......... .......... .......... 35% 93.4M 2s
    ##  48500K .......... .......... .......... .......... .......... 36%  104M 2s
    ##  48550K .......... .......... .......... .......... .......... 36%  167M 2s
    ##  48600K .......... .......... .......... .......... .......... 36% 62.4M 2s
    ##  48650K .......... .......... .......... .......... .......... 36%  120M 2s
    ##  48700K .......... .......... .......... .......... .......... 36% 36.4M 2s
    ##  48750K .......... .......... .......... .......... .......... 36% 77.6M 2s
    ##  48800K .......... .......... .......... .......... .......... 36% 54.3M 2s
    ##  48850K .......... .......... .......... .......... .......... 36%  112M 2s
    ##  48900K .......... .......... .......... .......... .......... 36%  117M 2s
    ##  48950K .......... .......... .......... .......... .......... 36%  133M 2s
    ##  49000K .......... .......... .......... .......... .......... 36% 47.4M 2s
    ##  49050K .......... .......... .......... .......... .......... 36% 20.5M 2s
    ##  49100K .......... .......... .......... .......... .......... 36% 93.7M 2s
    ##  49150K .......... .......... .......... .......... .......... 36%  160M 2s
    ##  49200K .......... .......... .......... .......... .......... 36%  147M 2s
    ##  49250K .......... .......... .......... .......... .......... 36%  162M 2s
    ##  49300K .......... .......... .......... .......... .......... 36%  155M 2s
    ##  49350K .......... .......... .......... .......... .......... 36% 16.1M 2s
    ##  49400K .......... .......... .......... .......... .......... 36% 59.1M 2s
    ##  49450K .......... .......... .......... .......... .......... 36% 68.8M 2s
    ##  49500K .......... .......... .......... .......... .......... 36%  111M 2s
    ##  49550K .......... .......... .......... .......... .......... 36%  116M 2s
    ##  49600K .......... .......... .......... .......... .......... 36%  159M 2s
    ##  49650K .......... .......... .......... .......... .......... 36%  149M 2s
    ##  49700K .......... .......... .......... .......... .......... 36% 25.6M 2s
    ##  49750K .......... .......... .......... .......... .......... 36% 77.7M 2s
    ##  49800K .......... .......... .......... .......... .......... 36%  108M 2s
    ##  49850K .......... .......... .......... .......... .......... 37% 46.7M 2s
    ##  49900K .......... .......... .......... .......... .......... 37%  136M 2s
    ##  49950K .......... .......... .......... .......... .......... 37%  165M 2s
    ##  50000K .......... .......... .......... .......... .......... 37% 85.5M 2s
    ##  50050K .......... .......... .......... .......... .......... 37% 53.8M 2s
    ##  50100K .......... .......... .......... .......... .......... 37%  139M 2s
    ##  50150K .......... .......... .......... .......... .......... 37% 52.2M 2s
    ##  50200K .......... .......... .......... .......... .......... 37%  122M 2s
    ##  50250K .......... .......... .......... .......... .......... 37% 90.2M 2s
    ##  50300K .......... .......... .......... .......... .......... 37%  121M 2s
    ##  50350K .......... .......... .......... .......... .......... 37% 88.7M 2s
    ##  50400K .......... .......... .......... .......... .......... 37%  122M 2s
    ##  50450K .......... .......... .......... .......... .......... 37% 48.0M 2s
    ##  50500K .......... .......... .......... .......... .......... 37% 87.7M 2s
    ##  50550K .......... .......... .......... .......... .......... 37%  163M 2s
    ##  50600K .......... .......... .......... .......... .......... 37% 87.7M 2s
    ##  50650K .......... .......... .......... .......... .......... 37% 50.3M 2s
    ##  50700K .......... .......... .......... .......... .......... 37%  138M 2s
    ##  50750K .......... .......... .......... .......... .......... 37% 56.3M 2s
    ##  50800K .......... .......... .......... .......... .......... 37%  149M 2s
    ##  50850K .......... .......... .......... .......... .......... 37%  136M 2s
    ##  50900K .......... .......... .......... .......... .......... 37%  134M 2s
    ##  50950K .......... .......... .......... .......... .......... 37% 48.6M 2s
    ##  51000K .......... .......... .......... .......... .......... 37%  107M 2s
    ##  51050K .......... .......... .......... .......... .......... 37%  157M 2s
    ##  51100K .......... .......... .......... .......... .......... 37% 45.2M 2s
    ##  51150K .......... .......... .......... .......... .......... 37%  130M 2s
    ##  51200K .......... .......... .......... .......... .......... 38%  149M 2s
    ##  51250K .......... .......... .......... .......... .......... 38% 58.2M 2s
    ##  51300K .......... .......... .......... .......... .......... 38%  129M 2s
    ##  51350K .......... .......... .......... .......... .......... 38% 53.1M 2s
    ##  51400K .......... .......... .......... .......... .......... 38%  101M 2s
    ##  51450K .......... .......... .......... .......... .......... 38% 66.1M 2s
    ##  51500K .......... .......... .......... .......... .......... 38%  139M 2s
    ##  51550K .......... .......... .......... .......... .......... 38% 90.3M 2s
    ##  51600K .......... .......... .......... .......... .......... 38% 98.3M 2s
    ##  51650K .......... .......... .......... .......... .......... 38% 68.0M 2s
    ##  51700K .......... .......... .......... .......... .......... 38% 96.4M 2s
    ##  51750K .......... .......... .......... .......... .......... 38%  123M 2s
    ##  51800K .......... .......... .......... .......... .......... 38% 61.2M 2s
    ##  51850K .......... .......... .......... .......... .......... 38% 33.4M 2s
    ##  51900K .......... .......... .......... .......... .......... 38% 45.0M 2s
    ##  51950K .......... .......... .......... .......... .......... 38% 66.8M 2s
    ##  52000K .......... .......... .......... .......... .......... 38% 73.7M 2s
    ##  52050K .......... .......... .......... .......... .......... 38%  114M 2s
    ##  52100K .......... .......... .......... .......... .......... 38%  131M 2s
    ##  52150K .......... .......... .......... .......... .......... 38%  117M 2s
    ##  52200K .......... .......... .......... .......... .......... 38% 64.6M 2s
    ##  52250K .......... .......... .......... .......... .......... 38%  105M 2s
    ##  52300K .......... .......... .......... .......... .......... 38% 64.9M 1s
    ##  52350K .......... .......... .......... .......... .......... 38% 69.4M 1s
    ##  52400K .......... .......... .......... .......... .......... 38% 73.1M 1s
    ##  52450K .......... .......... .......... .......... .......... 38%  115M 1s
    ##  52500K .......... .......... .......... .......... .......... 39%  124M 1s
    ##  52550K .......... .......... .......... .......... .......... 39% 51.3M 1s
    ##  52600K .......... .......... .......... .......... .......... 39%  115M 1s
    ##  52650K .......... .......... .......... .......... .......... 39% 31.4M 1s
    ##  52700K .......... .......... .......... .......... .......... 39% 90.9M 1s
    ##  52750K .......... .......... .......... .......... .......... 39% 38.9M 1s
    ##  52800K .......... .......... .......... .......... .......... 39% 50.9M 1s
    ##  52850K .......... .......... .......... .......... .......... 39% 79.7M 1s
    ##  52900K .......... .......... .......... .......... .......... 39% 62.8M 1s
    ##  52950K .......... .......... .......... .......... .......... 39% 92.0M 1s
    ##  53000K .......... .......... .......... .......... .......... 39% 48.8M 1s
    ##  53050K .......... .......... .......... .......... .......... 39% 80.8M 1s
    ##  53100K .......... .......... .......... .......... .......... 39% 45.6M 1s
    ##  53150K .......... .......... .......... .......... .......... 39% 44.9M 1s
    ##  53200K .......... .......... .......... .......... .......... 39%  133M 1s
    ##  53250K .......... .......... .......... .......... .......... 39% 68.7M 1s
    ##  53300K .......... .......... .......... .......... .......... 39% 45.2M 1s
    ##  53350K .......... .......... .......... .......... .......... 39% 62.2M 1s
    ##  53400K .......... .......... .......... .......... .......... 39% 54.0M 1s
    ##  53450K .......... .......... .......... .......... .......... 39% 52.9M 1s
    ##  53500K .......... .......... .......... .......... .......... 39% 94.7M 1s
    ##  53550K .......... .......... .......... .......... .......... 39% 90.3M 1s
    ##  53600K .......... .......... .......... .......... .......... 39% 57.1M 1s
    ##  53650K .......... .......... .......... .......... .......... 39%  101M 1s
    ##  53700K .......... .......... .......... .......... .......... 39% 59.5M 1s
    ##  53750K .......... .......... .......... .......... .......... 39% 76.1M 1s
    ##  53800K .......... .......... .......... .......... .......... 39% 48.8M 1s
    ##  53850K .......... .......... .......... .......... .......... 40%  138M 1s
    ##  53900K .......... .......... .......... .......... .......... 40% 46.9M 1s
    ##  53950K .......... .......... .......... .......... .......... 40% 76.7M 1s
    ##  54000K .......... .......... .......... .......... .......... 40% 71.8M 1s
    ##  54050K .......... .......... .......... .......... .......... 40% 76.1M 1s
    ##  54100K .......... .......... .......... .......... .......... 40%  106M 1s
    ##  54150K .......... .......... .......... .......... .......... 40% 39.6M 1s
    ##  54200K .......... .......... .......... .......... .......... 40%  106M 1s
    ##  54250K .......... .......... .......... .......... .......... 40% 59.0M 1s
    ##  54300K .......... .......... .......... .......... .......... 40% 64.3M 1s
    ##  54350K .......... .......... .......... .......... .......... 40% 67.3M 1s
    ##  54400K .......... .......... .......... .......... .......... 40% 16.3M 1s
    ##  54450K .......... .......... .......... .......... .......... 40%  140M 1s
    ##  54500K .......... .......... .......... .......... .......... 40%  124M 1s
    ##  54550K .......... .......... .......... .......... .......... 40%  141M 1s
    ##  54600K .......... .......... .......... .......... .......... 40%  150M 1s
    ##  54650K .......... .......... .......... .......... .......... 40%  138M 1s
    ##  54700K .......... .......... .......... .......... .......... 40% 34.6M 1s
    ##  54750K .......... .......... .......... .......... .......... 40% 17.7M 1s
    ##  54800K .......... .......... .......... .......... .......... 40% 76.1M 1s
    ##  54850K .......... .......... .......... .......... .......... 40%  130M 1s
    ##  54900K .......... .......... .......... .......... .......... 40%  121M 1s
    ##  54950K .......... .......... .......... .......... .......... 40%  130M 1s
    ##  55000K .......... .......... .......... .......... .......... 40%  132M 1s
    ##  55050K .......... .......... .......... .......... .......... 40% 36.3M 1s
    ##  55100K .......... .......... .......... .......... .......... 40% 38.9M 1s
    ##  55150K .......... .......... .......... .......... .......... 40% 46.1M 1s
    ##  55200K .......... .......... .......... .......... .......... 41% 47.8M 1s
    ##  55250K .......... .......... .......... .......... .......... 41% 88.7M 1s
    ##  55300K .......... .......... .......... .......... .......... 41% 95.4M 1s
    ##  55350K .......... .......... .......... .......... .......... 41%  101M 1s
    ##  55400K .......... .......... .......... .......... .......... 41% 98.0M 1s
    ##  55450K .......... .......... .......... .......... .......... 41% 90.6M 1s
    ##  55500K .......... .......... .......... .......... .......... 41% 43.5M 1s
    ##  55550K .......... .......... .......... .......... .......... 41%  144M 1s
    ##  55600K .......... .......... .......... .......... .......... 41% 18.0M 1s
    ##  55650K .......... .......... .......... .......... .......... 41%  134M 1s
    ##  55700K .......... .......... .......... .......... .......... 41%  111M 1s
    ##  55750K .......... .......... .......... .......... .......... 41%  103M 1s
    ##  55800K .......... .......... .......... .......... .......... 41%  128M 1s
    ##  55850K .......... .......... .......... .......... .......... 41%  110M 1s
    ##  55900K .......... .......... .......... .......... .......... 41% 28.1M 1s
    ##  55950K .......... .......... .......... .......... .......... 41%  116M 1s
    ##  56000K .......... .......... .......... .......... .......... 41%  101M 1s
    ##  56050K .......... .......... .......... .......... .......... 41%  127M 1s
    ##  56100K .......... .......... .......... .......... .......... 41%  124M 1s
    ##  56150K .......... .......... .......... .......... .......... 41% 50.9M 1s
    ##  56200K .......... .......... .......... .......... .......... 41% 35.7M 1s
    ##  56250K .......... .......... .......... .......... .......... 41% 69.6M 1s
    ##  56300K .......... .......... .......... .......... .......... 41% 45.8M 1s
    ##  56350K .......... .......... .......... .......... .......... 41% 86.2M 1s
    ##  56400K .......... .......... .......... .......... .......... 41% 42.1M 1s
    ##  56450K .......... .......... .......... .......... .......... 41%  117M 1s
    ##  56500K .......... .......... .......... .......... .......... 41%  164M 1s
    ##  56550K .......... .......... .......... .......... .......... 42% 77.0M 1s
    ##  56600K .......... .......... .......... .......... .......... 42%  125M 1s
    ##  56650K .......... .......... .......... .......... .......... 42% 84.0M 1s
    ##  56700K .......... .......... .......... .......... .......... 42% 34.5M 1s
    ##  56750K .......... .......... .......... .......... .......... 42%  147M 1s
    ##  56800K .......... .......... .......... .......... .......... 42%  100M 1s
    ##  56850K .......... .......... .......... .......... .......... 42% 37.6M 1s
    ##  56900K .......... .......... .......... .......... .......... 42%  110M 1s
    ##  56950K .......... .......... .......... .......... .......... 42%  128M 1s
    ##  57000K .......... .......... .......... .......... .......... 42% 83.3M 1s
    ##  57050K .......... .......... .......... .......... .......... 42% 23.1M 1s
    ##  57100K .......... .......... .......... .......... .......... 42% 55.7M 1s
    ##  57150K .......... .......... .......... .......... .......... 42% 75.9M 1s
    ##  57200K .......... .......... .......... .......... .......... 42% 42.3M 1s
    ##  57250K .......... .......... .......... .......... .......... 42% 91.2M 1s
    ##  57300K .......... .......... .......... .......... .......... 42%  119M 1s
    ##  57350K .......... .......... .......... .......... .......... 42%  161M 1s
    ##  57400K .......... .......... .......... .......... .......... 42% 91.1M 1s
    ##  57450K .......... .......... .......... .......... .......... 42%  145M 1s
    ##  57500K .......... .......... .......... .......... .......... 42% 31.1M 1s
    ##  57550K .......... .......... .......... .......... .......... 42% 90.9M 1s
    ##  57600K .......... .......... .......... .......... .......... 42% 72.4M 1s
    ##  57650K .......... .......... .......... .......... .......... 42% 77.4M 1s
    ##  57700K .......... .......... .......... .......... .......... 42%  116M 1s
    ##  57750K .......... .......... .......... .......... .......... 42%  104M 1s
    ##  57800K .......... .......... .......... .......... .......... 42%  100M 1s
    ##  57850K .......... .......... .......... .......... .......... 42% 68.8M 1s
    ##  57900K .......... .......... .......... .......... .......... 43%  123M 1s
    ##  57950K .......... .......... .......... .......... .......... 43%  110M 1s
    ##  58000K .......... .......... .......... .......... .......... 43% 70.5M 1s
    ##  58050K .......... .......... .......... .......... .......... 43% 95.1M 1s
    ##  58100K .......... .......... .......... .......... .......... 43% 84.7M 1s
    ##  58150K .......... .......... .......... .......... .......... 43%  126M 1s
    ##  58200K .......... .......... .......... .......... .......... 43% 62.9M 1s
    ##  58250K .......... .......... .......... .......... .......... 43% 83.7M 1s
    ##  58300K .......... .......... .......... .......... .......... 43%  103M 1s
    ##  58350K .......... .......... .......... .......... .......... 43%  149M 1s
    ##  58400K .......... .......... .......... .......... .......... 43%  119M 1s
    ##  58450K .......... .......... .......... .......... .......... 43% 43.2M 1s
    ##  58500K .......... .......... .......... .......... .......... 43%  108M 1s
    ##  58550K .......... .......... .......... .......... .......... 43%  127M 1s
    ##  58600K .......... .......... .......... .......... .......... 43%  115M 1s
    ##  58650K .......... .......... .......... .......... .......... 43%  132M 1s
    ##  58700K .......... .......... .......... .......... .......... 43% 93.9M 1s
    ##  58750K .......... .......... .......... .......... .......... 43% 77.6M 1s
    ##  58800K .......... .......... .......... .......... .......... 43% 59.0M 1s
    ##  58850K .......... .......... .......... .......... .......... 43%  140M 1s
    ##  58900K .......... .......... .......... .......... .......... 43% 94.5M 1s
    ##  58950K .......... .......... .......... .......... .......... 43% 54.0M 1s
    ##  59000K .......... .......... .......... .......... .......... 43% 88.9M 1s
    ##  59050K .......... .......... .......... .......... .......... 43%  124M 1s
    ##  59100K .......... .......... .......... .......... .......... 43%  111M 1s
    ##  59150K .......... .......... .......... .......... .......... 43%  100M 1s
    ##  59200K .......... .......... .......... .......... .......... 43% 85.6M 1s
    ##  59250K .......... .......... .......... .......... .......... 44%  136M 1s
    ##  59300K .......... .......... .......... .......... .......... 44% 51.8M 1s
    ##  59350K .......... .......... .......... .......... .......... 44%  109M 1s
    ##  59400K .......... .......... .......... .......... .......... 44% 83.6M 1s
    ##  59450K .......... .......... .......... .......... .......... 44% 89.5M 1s
    ##  59500K .......... .......... .......... .......... .......... 44%  106M 1s
    ##  59550K .......... .......... .......... .......... .......... 44% 71.8M 1s
    ##  59600K .......... .......... .......... .......... .......... 44%  114M 1s
    ##  59650K .......... .......... .......... .......... .......... 44% 34.6M 1s
    ##  59700K .......... .......... .......... .......... .......... 44% 51.5M 1s
    ##  59750K .......... .......... .......... .......... .......... 44%  152M 1s
    ##  59800K .......... .......... .......... .......... .......... 44% 83.2M 1s
    ##  59850K .......... .......... .......... .......... .......... 44%  117M 1s
    ##  59900K .......... .......... .......... .......... .......... 44%  115M 1s
    ##  59950K .......... .......... .......... .......... .......... 44% 66.1M 1s
    ##  60000K .......... .......... .......... .......... .......... 44% 63.6M 1s
    ##  60050K .......... .......... .......... .......... .......... 44% 68.6M 1s
    ##  60100K .......... .......... .......... .......... .......... 44%  100M 1s
    ##  60150K .......... .......... .......... .......... .......... 44% 96.8M 1s
    ##  60200K .......... .......... .......... .......... .......... 44% 71.7M 1s
    ##  60250K .......... .......... .......... .......... .......... 44%  145M 1s
    ##  60300K .......... .......... .......... .......... .......... 44% 60.5M 1s
    ##  60350K .......... .......... .......... .......... .......... 44% 44.7M 1s
    ##  60400K .......... .......... .......... .......... .......... 44% 92.5M 1s
    ##  60450K .......... .......... .......... .......... .......... 44%  126M 1s
    ##  60500K .......... .......... .......... .......... .......... 44%  139M 1s
    ##  60550K .......... .......... .......... .......... .......... 44% 45.5M 1s
    ##  60600K .......... .......... .......... .......... .......... 45%  130M 1s
    ##  60650K .......... .......... .......... .......... .......... 45%  129M 1s
    ##  60700K .......... .......... .......... .......... .......... 45% 99.9M 1s
    ##  60750K .......... .......... .......... .......... .......... 45% 72.3M 1s
    ##  60800K .......... .......... .......... .......... .......... 45%  105M 1s
    ##  60850K .......... .......... .......... .......... .......... 45% 58.1M 1s
    ##  60900K .......... .......... .......... .......... .......... 45% 49.9M 1s
    ##  60950K .......... .......... .......... .......... .......... 45%  157M 1s
    ##  61000K .......... .......... .......... .......... .......... 45% 75.5M 1s
    ##  61050K .......... .......... .......... .......... .......... 45%  110M 1s
    ##  61100K .......... .......... .......... .......... .......... 45% 73.6M 1s
    ##  61150K .......... .......... .......... .......... .......... 45%  103M 1s
    ##  61200K .......... .......... .......... .......... .......... 45%  105M 1s
    ##  61250K .......... .......... .......... .......... .......... 45% 63.1M 1s
    ##  61300K .......... .......... .......... .......... .......... 45%  128M 1s
    ##  61350K .......... .......... .......... .......... .......... 45%  107M 1s
    ##  61400K .......... .......... .......... .......... .......... 45%  127M 1s
    ##  61450K .......... .......... .......... .......... .......... 45%  130M 1s
    ##  61500K .......... .......... .......... .......... .......... 45% 46.8M 1s
    ##  61550K .......... .......... .......... .......... .......... 45% 67.0M 1s
    ##  61600K .......... .......... .......... .......... .......... 45% 45.5M 1s
    ##  61650K .......... .......... .......... .......... .......... 45% 81.8M 1s
    ##  61700K .......... .......... .......... .......... .......... 45% 95.6M 1s
    ##  61750K .......... .......... .......... .......... .......... 45%  129M 1s
    ##  61800K .......... .......... .......... .......... .......... 45%  126M 1s
    ##  61850K .......... .......... .......... .......... .......... 45% 92.1M 1s
    ##  61900K .......... .......... .......... .......... .......... 45%  123M 1s
    ##  61950K .......... .......... .......... .......... .......... 46% 47.4M 1s
    ##  62000K .......... .......... .......... .......... .......... 46% 51.5M 1s
    ##  62050K .......... .......... .......... .......... .......... 46% 93.8M 1s
    ##  62100K .......... .......... .......... .......... .......... 46% 84.3M 1s
    ##  62150K .......... .......... .......... .......... .......... 46%  144M 1s
    ##  62200K .......... .......... .......... .......... .......... 46%  111M 1s
    ##  62250K .......... .......... .......... .......... .......... 46% 70.8M 1s
    ##  62300K .......... .......... .......... .......... .......... 46% 60.2M 1s
    ##  62350K .......... .......... .......... .......... .......... 46% 89.4M 1s
    ##  62400K .......... .......... .......... .......... .......... 46% 60.8M 1s
    ##  62450K .......... .......... .......... .......... .......... 46%  165M 1s
    ##  62500K .......... .......... .......... .......... .......... 46%  103M 1s
    ##  62550K .......... .......... .......... .......... .......... 46%  110M 1s
    ##  62600K .......... .......... .......... .......... .......... 46%  134M 1s
    ##  62650K .......... .......... .......... .......... .......... 46% 46.2M 1s
    ##  62700K .......... .......... .......... .......... .......... 46% 29.0M 1s
    ##  62750K .......... .......... .......... .......... .......... 46%  116M 1s
    ##  62800K .......... .......... .......... .......... .......... 46%  148M 1s
    ##  62850K .......... .......... .......... .......... .......... 46%  137M 1s
    ##  62900K .......... .......... .......... .......... .......... 46%  133M 1s
    ##  62950K .......... .......... .......... .......... .......... 46% 20.5M 1s
    ##  63000K .......... .......... .......... .......... .......... 46% 51.7M 1s
    ##  63050K .......... .......... .......... .......... .......... 46% 54.3M 1s
    ##  63100K .......... .......... .......... .......... .......... 46%  127M 1s
    ##  63150K .......... .......... .......... .......... .......... 46% 71.4M 1s
    ##  63200K .......... .......... .......... .......... .......... 46% 81.5M 1s
    ##  63250K .......... .......... .......... .......... .......... 46%  158M 1s
    ##  63300K .......... .......... .......... .......... .......... 47% 44.9M 1s
    ##  63350K .......... .......... .......... .......... .......... 47%  141M 1s
    ##  63400K .......... .......... .......... .......... .......... 47%  104M 1s
    ##  63450K .......... .......... .......... .......... .......... 47% 52.3M 1s
    ##  63500K .......... .......... .......... .......... .......... 47%  111M 1s
    ##  63550K .......... .......... .......... .......... .......... 47%  128M 1s
    ##  63600K .......... .......... .......... .......... .......... 47%  146M 1s
    ##  63650K .......... .......... .......... .......... .......... 47% 32.0M 1s
    ##  63700K .......... .......... .......... .......... .......... 47%  122M 1s
    ##  63750K .......... .......... .......... .......... .......... 47% 52.1M 1s
    ##  63800K .......... .......... .......... .......... .......... 47% 57.8M 1s
    ##  63850K .......... .......... .......... .......... .......... 47%  103M 1s
    ##  63900K .......... .......... .......... .......... .......... 47%  139M 1s
    ##  63950K .......... .......... .......... .......... .......... 47%  123M 1s
    ##  64000K .......... .......... .......... .......... .......... 47%  110M 1s
    ##  64050K .......... .......... .......... .......... .......... 47%  123M 1s
    ##  64100K .......... .......... .......... .......... .......... 47% 32.7M 1s
    ##  64150K .......... .......... .......... .......... .......... 47%  112M 1s
    ##  64200K .......... .......... .......... .......... .......... 47%  153M 1s
    ##  64250K .......... .......... .......... .......... .......... 47%  121M 1s
    ##  64300K .......... .......... .......... .......... .......... 47%  103M 1s
    ##  64350K .......... .......... .......... .......... .......... 47%  104M 1s
    ##  64400K .......... .......... .......... .......... .......... 47% 16.2M 1s
    ##  64450K .......... .......... .......... .......... .......... 47%  149M 1s
    ##  64500K .......... .......... .......... .......... .......... 47%  115M 1s
    ##  64550K .......... .......... .......... .......... .......... 47%  124M 1s
    ##  64600K .......... .......... .......... .......... .......... 47%  133M 1s
    ##  64650K .......... .......... .......... .......... .......... 48%  121M 1s
    ##  64700K .......... .......... .......... .......... .......... 48%  146M 1s
    ##  64750K .......... .......... .......... .......... .......... 48%  150M 1s
    ##  64800K .......... .......... .......... .......... .......... 48% 23.5M 1s
    ##  64850K .......... .......... .......... .......... .......... 48% 61.6M 1s
    ##  64900K .......... .......... .......... .......... .......... 48%  110M 1s
    ##  64950K .......... .......... .......... .......... .......... 48% 51.7M 1s
    ##  65000K .......... .......... .......... .......... .......... 48%  116M 1s
    ##  65050K .......... .......... .......... .......... .......... 48%  163M 1s
    ##  65100K .......... .......... .......... .......... .......... 48%  137M 1s
    ##  65150K .......... .......... .......... .......... .......... 48% 76.6M 1s
    ##  65200K .......... .......... .......... .......... .......... 48% 95.8M 1s
    ##  65250K .......... .......... .......... .......... .......... 48%  111M 1s
    ##  65300K .......... .......... .......... .......... .......... 48%  112M 1s
    ##  65350K .......... .......... .......... .......... .......... 48% 92.1M 1s
    ##  65400K .......... .......... .......... .......... .......... 48%  106M 1s
    ##  65450K .......... .......... .......... .......... .......... 48%  120M 1s
    ##  65500K .......... .......... .......... .......... .......... 48% 30.5M 1s
    ##  65550K .......... .......... .......... .......... .......... 48%  124M 1s
    ##  65600K .......... .......... .......... .......... .......... 48% 43.9M 1s
    ##  65650K .......... .......... .......... .......... .......... 48%  128M 1s
    ##  65700K .......... .......... .......... .......... .......... 48% 21.1M 1s
    ##  65750K .......... .......... .......... .......... .......... 48%  113M 1s
    ##  65800K .......... .......... .......... .......... .......... 48%  104M 1s
    ##  65850K .......... .......... .......... .......... .......... 48%  137M 1s
    ##  65900K .......... .......... .......... .......... .......... 48%  141M 1s
    ##  65950K .......... .......... .......... .......... .......... 48%  136M 1s
    ##  66000K .......... .......... .......... .......... .......... 49%  129M 1s
    ##  66050K .......... .......... .......... .......... .......... 49% 91.4M 1s
    ##  66100K .......... .......... .......... .......... .......... 49% 44.9M 1s
    ##  66150K .......... .......... .......... .......... .......... 49% 28.3M 1s
    ##  66200K .......... .......... .......... .......... .......... 49% 63.9M 1s
    ##  66250K .......... .......... .......... .......... .......... 49% 60.5M 1s
    ##  66300K .......... .......... .......... .......... .......... 49%  105M 1s
    ##  66350K .......... .......... .......... .......... .......... 49%  117M 1s
    ##  66400K .......... .......... .......... .......... .......... 49%  142M 1s
    ##  66450K .......... .......... .......... .......... .......... 49%  119M 1s
    ##  66500K .......... .......... .......... .......... .......... 49% 83.3M 1s
    ##  66550K .......... .......... .......... .......... .......... 49% 42.4M 1s
    ##  66600K .......... .......... .......... .......... .......... 49% 88.4M 1s
    ##  66650K .......... .......... .......... .......... .......... 49%  134M 1s
    ##  66700K .......... .......... .......... .......... .......... 49%  109M 1s
    ##  66750K .......... .......... .......... .......... .......... 49% 79.6M 1s
    ##  66800K .......... .......... .......... .......... .......... 49% 63.8M 1s
    ##  66850K .......... .......... .......... .......... .......... 49%  111M 1s
    ##  66900K .......... .......... .......... .......... .......... 49%  124M 1s
    ##  66950K .......... .......... .......... .......... .......... 49% 30.9M 1s
    ##  67000K .......... .......... .......... .......... .......... 49% 90.2M 1s
    ##  67050K .......... .......... .......... .......... .......... 49%  133M 1s
    ##  67100K .......... .......... .......... .......... .......... 49%  107M 1s
    ##  67150K .......... .......... .......... .......... .......... 49%  138M 1s
    ##  67200K .......... .......... .......... .......... .......... 49%  141M 1s
    ##  67250K .......... .......... .......... .......... .......... 49% 34.1M 1s
    ##  67300K .......... .......... .......... .......... .......... 49% 21.6M 1s
    ##  67350K .......... .......... .......... .......... .......... 50%  120M 1s
    ##  67400K .......... .......... .......... .......... .......... 50%  131M 1s
    ##  67450K .......... .......... .......... .......... .......... 50%  163M 1s
    ##  67500K .......... .......... .......... .......... .......... 50%  123M 1s
    ##  67550K .......... .......... .......... .......... .......... 50%  162M 1s
    ##  67600K .......... .......... .......... .......... .......... 50%  124M 1s
    ##  67650K .......... .......... .......... .......... .......... 50% 29.3M 1s
    ##  67700K .......... .......... .......... .......... .......... 50% 35.6M 1s
    ##  67750K .......... .......... .......... .......... .......... 50% 60.3M 1s
    ##  67800K .......... .......... .......... .......... .......... 50%  113M 1s
    ##  67850K .......... .......... .......... .......... .......... 50%  111M 1s
    ##  67900K .......... .......... .......... .......... .......... 50%  110M 1s
    ##  67950K .......... .......... .......... .......... .......... 50% 99.8M 1s
    ##  68000K .......... .......... .......... .......... .......... 50%  106M 1s
    ##  68050K .......... .......... .......... .......... .......... 50%  174M 1s
    ##  68100K .......... .......... .......... .......... .......... 50% 33.2M 1s
    ##  68150K .......... .......... .......... .......... .......... 50%  138M 1s
    ##  68200K .......... .......... .......... .......... .......... 50%  118M 1s
    ##  68250K .......... .......... .......... .......... .......... 50%  108M 1s
    ##  68300K .......... .......... .......... .......... .......... 50% 35.5M 1s
    ##  68350K .......... .......... .......... .......... .......... 50%  165M 1s
    ##  68400K .......... .......... .......... .......... .......... 50% 38.5M 1s
    ##  68450K .......... .......... .......... .......... .......... 50% 81.4M 1s
    ##  68500K .......... .......... .......... .......... .......... 50%  136M 1s
    ##  68550K .......... .......... .......... .......... .......... 50%  155M 1s
    ##  68600K .......... .......... .......... .......... .......... 50% 12.8M 1s
    ##  68650K .......... .......... .......... .......... .......... 50%  138M 1s
    ##  68700K .......... .......... .......... .......... .......... 51% 40.2M 1s
    ##  68750K .......... .......... .......... .......... .......... 51%  124M 1s
    ##  68800K .......... .......... .......... .......... .......... 51%  111M 1s
    ##  68850K .......... .......... .......... .......... .......... 51%  100M 1s
    ##  68900K .......... .......... .......... .......... .......... 51%  133M 1s
    ##  68950K .......... .......... .......... .......... .......... 51%  150M 1s
    ##  69000K .......... .......... .......... .......... .......... 51%  118M 1s
    ##  69050K .......... .......... .......... .......... .......... 51% 41.8M 1s
    ##  69100K .......... .......... .......... .......... .......... 51% 63.6M 1s
    ##  69150K .......... .......... .......... .......... .......... 51% 82.4M 1s
    ##  69200K .......... .......... .......... .......... .......... 51% 44.2M 1s
    ##  69250K .......... .......... .......... .......... .......... 51%  158M 1s
    ##  69300K .......... .......... .......... .......... .......... 51% 93.0M 1s
    ##  69350K .......... .......... .......... .......... .......... 51% 92.9M 1s
    ##  69400K .......... .......... .......... .......... .......... 51%  106M 1s
    ##  69450K .......... .......... .......... .......... .......... 51%  142M 1s
    ##  69500K .......... .......... .......... .......... .......... 51% 18.0M 1s
    ##  69550K .......... .......... .......... .......... .......... 51%  154M 1s
    ##  69600K .......... .......... .......... .......... .......... 51% 89.9M 1s
    ##  69650K .......... .......... .......... .......... .......... 51% 86.8M 1s
    ##  69700K .......... .......... .......... .......... .......... 51%  123M 1s
    ##  69750K .......... .......... .......... .......... .......... 51%  149M 1s
    ##  69800K .......... .......... .......... .......... .......... 51%  146M 1s
    ##  69850K .......... .......... .......... .......... .......... 51% 52.1M 1s
    ##  69900K .......... .......... .......... .......... .......... 51% 57.8M 1s
    ##  69950K .......... .......... .......... .......... .......... 51% 19.4M 1s
    ##  70000K .......... .......... .......... .......... .......... 51%  143M 1s
    ##  70050K .......... .......... .......... .......... .......... 52% 53.0M 1s
    ##  70100K .......... .......... .......... .......... .......... 52%  132M 1s
    ##  70150K .......... .......... .......... .......... .......... 52%  144M 1s
    ##  70200K .......... .......... .......... .......... .......... 52%  120M 1s
    ##  70250K .......... .......... .......... .......... .......... 52%  131M 1s
    ##  70300K .......... .......... .......... .......... .......... 52% 49.0M 1s
    ##  70350K .......... .......... .......... .......... .......... 52% 69.2M 1s
    ##  70400K .......... .......... .......... .......... .......... 52% 81.5M 1s
    ##  70450K .......... .......... .......... .......... .......... 52% 69.0M 1s
    ##  70500K .......... .......... .......... .......... .......... 52% 66.5M 1s
    ##  70550K .......... .......... .......... .......... .......... 52% 55.9M 1s
    ##  70600K .......... .......... .......... .......... .......... 52% 69.7M 1s
    ##  70650K .......... .......... .......... .......... .......... 52% 31.5M 1s
    ##  70700K .......... .......... .......... .......... .......... 52% 78.1M 1s
    ##  70750K .......... .......... .......... .......... .......... 52%  110M 1s
    ##  70800K .......... .......... .......... .......... .......... 52% 30.9M 1s
    ##  70850K .......... .......... .......... .......... .......... 52% 60.9M 1s
    ##  70900K .......... .......... .......... .......... .......... 52% 74.5M 1s
    ##  70950K .......... .......... .......... .......... .......... 52% 97.2M 1s
    ##  71000K .......... .......... .......... .......... .......... 52% 91.6M 1s
    ##  71050K .......... .......... .......... .......... .......... 52%  112M 1s
    ##  71100K .......... .......... .......... .......... .......... 52% 78.8M 1s
    ##  71150K .......... .......... .......... .......... .......... 52% 92.7M 1s
    ##  71200K .......... .......... .......... .......... .......... 52% 80.2M 1s
    ##  71250K .......... .......... .......... .......... .......... 52% 86.0M 1s
    ##  71300K .......... .......... .......... .......... .......... 52% 69.0M 1s
    ##  71350K .......... .......... .......... .......... .......... 52% 77.1M 1s
    ##  71400K .......... .......... .......... .......... .......... 53% 71.0M 1s
    ##  71450K .......... .......... .......... .......... .......... 53% 75.0M 1s
    ##  71500K .......... .......... .......... .......... .......... 53% 62.9M 1s
    ##  71550K .......... .......... .......... .......... .......... 53% 68.3M 1s
    ##  71600K .......... .......... .......... .......... .......... 53% 76.9M 1s
    ##  71650K .......... .......... .......... .......... .......... 53%  103M 1s
    ##  71700K .......... .......... .......... .......... .......... 53% 73.2M 1s
    ##  71750K .......... .......... .......... .......... .......... 53% 60.9M 1s
    ##  71800K .......... .......... .......... .......... .......... 53% 71.4M 1s
    ##  71850K .......... .......... .......... .......... .......... 53% 79.8M 1s
    ##  71900K .......... .......... .......... .......... .......... 53% 65.6M 1s
    ##  71950K .......... .......... .......... .......... .......... 53% 93.7M 1s
    ##  72000K .......... .......... .......... .......... .......... 53% 46.8M 1s
    ##  72050K .......... .......... .......... .......... .......... 53% 87.6M 1s
    ##  72100K .......... .......... .......... .......... .......... 53% 71.6M 1s
    ##  72150K .......... .......... .......... .......... .......... 53%  113M 1s
    ##  72200K .......... .......... .......... .......... .......... 53% 93.0M 1s
    ##  72250K .......... .......... .......... .......... .......... 53% 89.8M 1s
    ##  72300K .......... .......... .......... .......... .......... 53% 82.1M 1s
    ##  72350K .......... .......... .......... .......... .......... 53% 95.3M 1s
    ##  72400K .......... .......... .......... .......... .......... 53% 85.6M 1s
    ##  72450K .......... .......... .......... .......... .......... 53%  103M 1s
    ##  72500K .......... .......... .......... .......... .......... 53% 61.7M 1s
    ##  72550K .......... .......... .......... .......... .......... 53% 77.7M 1s
    ##  72600K .......... .......... .......... .......... .......... 53% 47.7M 1s
    ##  72650K .......... .......... .......... .......... .......... 53% 98.6M 1s
    ##  72700K .......... .......... .......... .......... .......... 53% 84.9M 1s
    ##  72750K .......... .......... .......... .......... .......... 54%  113M 1s
    ##  72800K .......... .......... .......... .......... .......... 54% 79.7M 1s
    ##  72850K .......... .......... .......... .......... .......... 54% 80.9M 1s
    ##  72900K .......... .......... .......... .......... .......... 54%  102M 1s
    ##  72950K .......... .......... .......... .......... .......... 54%  107M 1s
    ##  73000K .......... .......... .......... .......... .......... 54% 93.9M 1s
    ##  73050K .......... .......... .......... .......... .......... 54%  118M 1s
    ##  73100K .......... .......... .......... .......... .......... 54% 75.9M 1s
    ##  73150K .......... .......... .......... .......... .......... 54% 83.6M 1s
    ##  73200K .......... .......... .......... .......... .......... 54% 70.3M 1s
    ##  73250K .......... .......... .......... .......... .......... 54%  116M 1s
    ##  73300K .......... .......... .......... .......... .......... 54% 82.4M 1s
    ##  73350K .......... .......... .......... .......... .......... 54%  119M 1s
    ##  73400K .......... .......... .......... .......... .......... 54% 80.3M 1s
    ##  73450K .......... .......... .......... .......... .......... 54% 88.5M 1s
    ##  73500K .......... .......... .......... .......... .......... 54% 93.1M 1s
    ##  73550K .......... .......... .......... .......... .......... 54%  115M 1s
    ##  73600K .......... .......... .......... .......... .......... 54% 91.3M 1s
    ##  73650K .......... .......... .......... .......... .......... 54% 43.5M 1s
    ##  73700K .......... .......... .......... .......... .......... 54% 71.6M 1s
    ##  73750K .......... .......... .......... .......... .......... 54% 80.1M 1s
    ##  73800K .......... .......... .......... .......... .......... 54% 88.5M 1s
    ##  73850K .......... .......... .......... .......... .......... 54%  106M 1s
    ##  73900K .......... .......... .......... .......... .......... 54%  101M 1s
    ##  73950K .......... .......... .......... .......... .......... 54% 78.2M 1s
    ##  74000K .......... .......... .......... .......... .......... 54% 81.3M 1s
    ##  74050K .......... .......... .......... .......... .......... 54%  105M 1s
    ##  74100K .......... .......... .......... .......... .......... 55%  101M 1s
    ##  74150K .......... .......... .......... .......... .......... 55% 65.5M 1s
    ##  74200K .......... .......... .......... .......... .......... 55% 71.1M 1s
    ##  74250K .......... .......... .......... .......... .......... 55% 39.7M 1s
    ##  74300K .......... .......... .......... .......... .......... 55% 79.4M 1s
    ##  74350K .......... .......... .......... .......... .......... 55% 99.8M 1s
    ##  74400K .......... .......... .......... .......... .......... 55% 98.6M 1s
    ##  74450K .......... .......... .......... .......... .......... 55% 73.6M 1s
    ##  74500K .......... .......... .......... .......... .......... 55%  115M 1s
    ##  74550K .......... .......... .......... .......... .......... 55% 91.6M 1s
    ##  74600K .......... .......... .......... .......... .......... 55% 77.0M 1s
    ##  74650K .......... .......... .......... .......... .......... 55% 91.2M 1s
    ##  74700K .......... .......... .......... .......... .......... 55%  101M 1s
    ##  74750K .......... .......... .......... .......... .......... 55% 84.8M 1s
    ##  74800K .......... .......... .......... .......... .......... 55% 94.4M 1s
    ##  74850K .......... .......... .......... .......... .......... 55%  114M 1s
    ##  74900K .......... .......... .......... .......... .......... 55% 81.4M 1s
    ##  74950K .......... .......... .......... .......... .......... 55%  132M 1s
    ##  75000K .......... .......... .......... .......... .......... 55% 67.3M 1s
    ##  75050K .......... .......... .......... .......... .......... 55% 89.1M 1s
    ##  75100K .......... .......... .......... .......... .......... 55%  105M 1s
    ##  75150K .......... .......... .......... .......... .......... 55%  106M 1s
    ##  75200K .......... .......... .......... .......... .......... 55% 93.6M 1s
    ##  75250K .......... .......... .......... .......... .......... 55%  109M 1s
    ##  75300K .......... .......... .......... .......... .......... 55% 91.2M 1s
    ##  75350K .......... .......... .......... .......... .......... 55% 91.2M 1s
    ##  75400K .......... .......... .......... .......... .......... 55% 53.5M 1s
    ##  75450K .......... .......... .......... .......... .......... 56% 75.2M 1s
    ##  75500K .......... .......... .......... .......... .......... 56% 85.5M 1s
    ##  75550K .......... .......... .......... .......... .......... 56% 60.4M 1s
    ##  75600K .......... .......... .......... .......... .......... 56% 94.2M 1s
    ##  75650K .......... .......... .......... .......... .......... 56% 97.3M 1s
    ##  75700K .......... .......... .......... .......... .......... 56% 55.0M 1s
    ##  75750K .......... .......... .......... .......... .......... 56% 89.1M 1s
    ##  75800K .......... .......... .......... .......... .......... 56%  106M 1s
    ##  75850K .......... .......... .......... .......... .......... 56% 78.8M 1s
    ##  75900K .......... .......... .......... .......... .......... 56% 90.0M 1s
    ##  75950K .......... .......... .......... .......... .......... 56% 84.0M 1s
    ##  76000K .......... .......... .......... .......... .......... 56%  130M 1s
    ##  76050K .......... .......... .......... .......... .......... 56%  124M 1s
    ##  76100K .......... .......... .......... .......... .......... 56% 78.3M 1s
    ##  76150K .......... .......... .......... .......... .......... 56%  101M 1s
    ##  76200K .......... .......... .......... .......... .......... 56% 69.1M 1s
    ##  76250K .......... .......... .......... .......... .......... 56% 86.3M 1s
    ##  76300K .......... .......... .......... .......... .......... 56% 86.2M 1s
    ##  76350K .......... .......... .......... .......... .......... 56%  110M 1s
    ##  76400K .......... .......... .......... .......... .......... 56% 90.7M 1s
    ##  76450K .......... .......... .......... .......... .......... 56% 89.7M 1s
    ##  76500K .......... .......... .......... .......... .......... 56% 82.8M 1s
    ##  76550K .......... .......... .......... .......... .......... 56% 77.0M 1s
    ##  76600K .......... .......... .......... .......... .......... 56%  110M 1s
    ##  76650K .......... .......... .......... .......... .......... 56%  106M 1s
    ##  76700K .......... .......... .......... .......... .......... 56%  106M 1s
    ##  76750K .......... .......... .......... .......... .......... 56%  126M 1s
    ##  76800K .......... .......... .......... .......... .......... 57% 21.5M 1s
    ##  76850K .......... .......... .......... .......... .......... 57% 85.9M 1s
    ##  76900K .......... .......... .......... .......... .......... 57% 93.4M 1s
    ##  76950K .......... .......... .......... .......... .......... 57%  101M 1s
    ##  77000K .......... .......... .......... .......... .......... 57%  106M 1s
    ##  77050K .......... .......... .......... .......... .......... 57%  121M 1s
    ##  77100K .......... .......... .......... .......... .......... 57%  108M 1s
    ##  77150K .......... .......... .......... .......... .......... 57%  140M 1s
    ##  77200K .......... .......... .......... .......... .......... 57% 39.4M 1s
    ##  77250K .......... .......... .......... .......... .......... 57%  121M 1s
    ##  77300K .......... .......... .......... .......... .......... 57% 28.3M 1s
    ##  77350K .......... .......... .......... .......... .......... 57% 95.1M 1s
    ##  77400K .......... .......... .......... .......... .......... 57% 85.8M 1s
    ##  77450K .......... .......... .......... .......... .......... 57% 57.6M 1s
    ##  77500K .......... .......... .......... .......... .......... 57% 98.5M 1s
    ##  77550K .......... .......... .......... .......... .......... 57% 96.1M 1s
    ##  77600K .......... .......... .......... .......... .......... 57%  110M 1s
    ##  77650K .......... .......... .......... .......... .......... 57% 93.8M 1s
    ##  77700K .......... .......... .......... .......... .......... 57% 89.1M 1s
    ##  77750K .......... .......... .......... .......... .......... 57%  141M 1s
    ##  77800K .......... .......... .......... .......... .......... 57% 83.6M 1s
    ##  77850K .......... .......... .......... .......... .......... 57%  134M 1s
    ##  77900K .......... .......... .......... .......... .......... 57% 87.7M 1s
    ##  77950K .......... .......... .......... .......... .......... 57% 91.3M 1s
    ##  78000K .......... .......... .......... .......... .......... 57% 69.0M 1s
    ##  78050K .......... .......... .......... .......... .......... 57%  146M 1s
    ##  78100K .......... .......... .......... .......... .......... 58% 88.5M 1s
    ##  78150K .......... .......... .......... .......... .......... 58%  121M 1s
    ##  78200K .......... .......... .......... .......... .......... 58%  106M 1s
    ##  78250K .......... .......... .......... .......... .......... 58%  120M 1s
    ##  78300K .......... .......... .......... .......... .......... 58% 50.7M 1s
    ##  78350K .......... .......... .......... .......... .......... 58%  127M 1s
    ##  78400K .......... .......... .......... .......... .......... 58% 90.9M 1s
    ##  78450K .......... .......... .......... .......... .......... 58% 64.9M 1s
    ##  78500K .......... .......... .......... .......... .......... 58% 86.2M 1s
    ##  78550K .......... .......... .......... .......... .......... 58% 93.9M 1s
    ##  78600K .......... .......... .......... .......... .......... 58%  103M 1s
    ##  78650K .......... .......... .......... .......... .......... 58%  123M 1s
    ##  78700K .......... .......... .......... .......... .......... 58%  113M 1s
    ##  78750K .......... .......... .......... .......... .......... 58%  103M 1s
    ##  78800K .......... .......... .......... .......... .......... 58% 50.5M 1s
    ##  78850K .......... .......... .......... .......... .......... 58% 88.9M 1s
    ##  78900K .......... .......... .......... .......... .......... 58% 71.9M 1s
    ##  78950K .......... .......... .......... .......... .......... 58%  129M 1s
    ##  79000K .......... .......... .......... .......... .......... 58%  110M 1s
    ##  79050K .......... .......... .......... .......... .......... 58%  112M 1s
    ##  79100K .......... .......... .......... .......... .......... 58% 87.9M 1s
    ##  79150K .......... .......... .......... .......... .......... 58% 99.2M 1s
    ##  79200K .......... .......... .......... .......... .......... 58%  113M 1s
    ##  79250K .......... .......... .......... .......... .......... 58%  105M 1s
    ##  79300K .......... .......... .......... .......... .......... 58%  106M 1s
    ##  79350K .......... .......... .......... .......... .......... 58%  110M 1s
    ##  79400K .......... .......... .......... .......... .......... 58% 93.8M 1s
    ##  79450K .......... .......... .......... .......... .......... 59% 99.5M 1s
    ##  79500K .......... .......... .......... .......... .......... 59%  103M 1s
    ##  79550K .......... .......... .......... .......... .......... 59%  137M 1s
    ##  79600K .......... .......... .......... .......... .......... 59%  140M 1s
    ##  79650K .......... .......... .......... .......... .......... 59% 50.4M 1s
    ##  79700K .......... .......... .......... .......... .......... 59%  126M 1s
    ##  79750K .......... .......... .......... .......... .......... 59% 35.9M 1s
    ##  79800K .......... .......... .......... .......... .......... 59% 86.7M 1s
    ##  79850K .......... .......... .......... .......... .......... 59% 40.0M 1s
    ##  79900K .......... .......... .......... .......... .......... 59% 24.3M 1s
    ##  79950K .......... .......... .......... .......... .......... 59%  138M 1s
    ##  80000K .......... .......... .......... .......... .......... 59%  121M 1s
    ##  80050K .......... .......... .......... .......... .......... 59%  136M 1s
    ##  80100K .......... .......... .......... .......... .......... 59% 96.8M 1s
    ##  80150K .......... .......... .......... .......... .......... 59%  141M 1s
    ##  80200K .......... .......... .......... .......... .......... 59%  111M 1s
    ##  80250K .......... .......... .......... .......... .......... 59%  159M 1s
    ##  80300K .......... .......... .......... .......... .......... 59% 64.7M 1s
    ##  80350K .......... .......... .......... .......... .......... 59%  102M 1s
    ##  80400K .......... .......... .......... .......... .......... 59% 38.2M 1s
    ##  80450K .......... .......... .......... .......... .......... 59% 76.1M 1s
    ##  80500K .......... .......... .......... .......... .......... 59% 87.6M 1s
    ##  80550K .......... .......... .......... .......... .......... 59%  117M 1s
    ##  80600K .......... .......... .......... .......... .......... 59% 36.8M 1s
    ##  80650K .......... .......... .......... .......... .......... 59%  100M 1s
    ##  80700K .......... .......... .......... .......... .......... 59% 89.7M 1s
    ##  80750K .......... .......... .......... .......... .......... 59%  138M 1s
    ##  80800K .......... .......... .......... .......... .......... 60% 58.2M 1s
    ##  80850K .......... .......... .......... .......... .......... 60% 33.8M 1s
    ##  80900K .......... .......... .......... .......... .......... 60% 76.2M 1s
    ##  80950K .......... .......... .......... .......... .......... 60% 59.5M 1s
    ##  81000K .......... .......... .......... .......... .......... 60% 81.4M 1s
    ##  81050K .......... .......... .......... .......... .......... 60%  127M 1s
    ##  81100K .......... .......... .......... .......... .......... 60%  105M 1s
    ##  81150K .......... .......... .......... .......... .......... 60% 91.3M 1s
    ##  81200K .......... .......... .......... .......... .......... 60% 39.3M 1s
    ##  81250K .......... .......... .......... .......... .......... 60%  103M 1s
    ##  81300K .......... .......... .......... .......... .......... 60% 80.2M 1s
    ##  81350K .......... .......... .......... .......... .......... 60% 48.3M 1s
    ##  81400K .......... .......... .......... .......... .......... 60% 91.6M 1s
    ##  81450K .......... .......... .......... .......... .......... 60%  130M 1s
    ##  81500K .......... .......... .......... .......... .......... 60% 48.0M 1s
    ##  81550K .......... .......... .......... .......... .......... 60% 94.4M 1s
    ##  81600K .......... .......... .......... .......... .......... 60% 49.8M 1s
    ##  81650K .......... .......... .......... .......... .......... 60%  111M 1s
    ##  81700K .......... .......... .......... .......... .......... 60% 97.3M 1s
    ##  81750K .......... .......... .......... .......... .......... 60% 44.5M 1s
    ##  81800K .......... .......... .......... .......... .......... 60%  102M 1s
    ##  81850K .......... .......... .......... .......... .......... 60% 70.5M 1s
    ##  81900K .......... .......... .......... .......... .......... 60% 92.4M 1s
    ##  81950K .......... .......... .......... .......... .......... 60%  134M 1s
    ##  82000K .......... .......... .......... .......... .......... 60% 86.1M 1s
    ##  82050K .......... .......... .......... .......... .......... 60% 86.3M 1s
    ##  82100K .......... .......... .......... .......... .......... 60% 40.8M 1s
    ##  82150K .......... .......... .......... .......... .......... 61% 93.5M 1s
    ##  82200K .......... .......... .......... .......... .......... 61%  117M 1s
    ##  82250K .......... .......... .......... .......... .......... 61% 57.4M 1s
    ##  82300K .......... .......... .......... .......... .......... 61% 79.5M 1s
    ##  82350K .......... .......... .......... .......... .......... 61%  132M 1s
    ##  82400K .......... .......... .......... .......... .......... 61% 94.2M 1s
    ##  82450K .......... .......... .......... .......... .......... 61%  103M 1s
    ##  82500K .......... .......... .......... .......... .......... 61% 44.8M 1s
    ##  82550K .......... .......... .......... .......... .......... 61%  155M 1s
    ##  82600K .......... .......... .......... .......... .......... 61% 79.6M 1s
    ##  82650K .......... .......... .......... .......... .......... 61% 78.4M 1s
    ##  82700K .......... .......... .......... .......... .......... 61% 84.1M 1s
    ##  82750K .......... .......... .......... .......... .......... 61% 80.3M 1s
    ##  82800K .......... .......... .......... .......... .......... 61% 54.8M 1s
    ##  82850K .......... .......... .......... .......... .......... 61%  141M 1s
    ##  82900K .......... .......... .......... .......... .......... 61% 88.2M 1s
    ##  82950K .......... .......... .......... .......... .......... 61%  123M 1s
    ##  83000K .......... .......... .......... .......... .......... 61% 43.9M 1s
    ##  83050K .......... .......... .......... .......... .......... 61%  100M 1s
    ##  83100K .......... .......... .......... .......... .......... 61%  137M 1s
    ##  83150K .......... .......... .......... .......... .......... 61% 59.5M 1s
    ##  83200K .......... .......... .......... .......... .......... 61%  111M 1s
    ##  83250K .......... .......... .......... .......... .......... 61% 75.6M 1s
    ##  83300K .......... .......... .......... .......... .......... 61%  124M 1s
    ##  83350K .......... .......... .......... .......... .......... 61% 91.8M 1s
    ##  83400K .......... .......... .......... .......... .......... 61% 56.4M 1s
    ##  83450K .......... .......... .......... .......... .......... 61%  107M 1s
    ##  83500K .......... .......... .......... .......... .......... 62% 95.1M 1s
    ##  83550K .......... .......... .......... .......... .......... 62% 68.4M 1s
    ##  83600K .......... .......... .......... .......... .......... 62%  110M 1s
    ##  83650K .......... .......... .......... .......... .......... 62% 64.5M 1s
    ##  83700K .......... .......... .......... .......... .......... 62%  119M 1s
    ##  83750K .......... .......... .......... .......... .......... 62% 66.7M 1s
    ##  83800K .......... .......... .......... .......... .......... 62% 89.7M 1s
    ##  83850K .......... .......... .......... .......... .......... 62%  125M 1s
    ##  83900K .......... .......... .......... .......... .......... 62% 79.6M 1s
    ##  83950K .......... .......... .......... .......... .......... 62% 73.8M 1s
    ##  84000K .......... .......... .......... .......... .......... 62% 66.6M 1s
    ##  84050K .......... .......... .......... .......... .......... 62%  101M 1s
    ##  84100K .......... .......... .......... .......... .......... 62%  104M 1s
    ##  84150K .......... .......... .......... .......... .......... 62% 56.1M 1s
    ##  84200K .......... .......... .......... .......... .......... 62% 87.1M 1s
    ##  84250K .......... .......... .......... .......... .......... 62%  157M 1s
    ##  84300K .......... .......... .......... .......... .......... 62%  116M 1s
    ##  84350K .......... .......... .......... .......... .......... 62% 89.0M 1s
    ##  84400K .......... .......... .......... .......... .......... 62% 56.4M 1s
    ##  84450K .......... .......... .......... .......... .......... 62%  106M 1s
    ##  84500K .......... .......... .......... .......... .......... 62%  139M 1s
    ##  84550K .......... .......... .......... .......... .......... 62% 73.0M 1s
    ##  84600K .......... .......... .......... .......... .......... 62%  117M 1s
    ##  84650K .......... .......... .......... .......... .......... 62% 84.1M 1s
    ##  84700K .......... .......... .......... .......... .......... 62%  112M 1s
    ##  84750K .......... .......... .......... .......... .......... 62%  129M 1s
    ##  84800K .......... .......... .......... .......... .......... 62% 53.4M 1s
    ##  84850K .......... .......... .......... .......... .......... 63%  111M 1s
    ##  84900K .......... .......... .......... .......... .......... 63% 83.7M 1s
    ##  84950K .......... .......... .......... .......... .......... 63% 92.0M 1s
    ##  85000K .......... .......... .......... .......... .......... 63%  141M 1s
    ##  85050K .......... .......... .......... .......... .......... 63% 88.6M 1s
    ##  85100K .......... .......... .......... .......... .......... 63% 87.8M 1s
    ##  85150K .......... .......... .......... .......... .......... 63% 66.6M 1s
    ##  85200K .......... .......... .......... .......... .......... 63%  102M 1s
    ##  85250K .......... .......... .......... .......... .......... 63%  102M 1s
    ##  85300K .......... .......... .......... .......... .......... 63% 71.1M 1s
    ##  85350K .......... .......... .......... .......... .......... 63%  118M 1s
    ##  85400K .......... .......... .......... .......... .......... 63% 86.6M 1s
    ##  85450K .......... .......... .......... .......... .......... 63% 82.4M 1s
    ##  85500K .......... .......... .......... .......... .......... 63%  128M 1s
    ##  85550K .......... .......... .......... .......... .......... 63% 78.2M 1s
    ##  85600K .......... .......... .......... .......... .......... 63%  103M 1s
    ##  85650K .......... .......... .......... .......... .......... 63%  115M 1s
    ##  85700K .......... .......... .......... .......... .......... 63% 92.8M 1s
    ##  85750K .......... .......... .......... .......... .......... 63% 97.9M 1s
    ##  85800K .......... .......... .......... .......... .......... 63% 71.1M 1s
    ##  85850K .......... .......... .......... .......... .......... 63% 99.1M 1s
    ##  85900K .......... .......... .......... .......... .......... 63%  108M 1s
    ##  85950K .......... .......... .......... .......... .......... 63%  120M 1s
    ##  86000K .......... .......... .......... .......... .......... 63% 90.4M 1s
    ##  86050K .......... .......... .......... .......... .......... 63% 88.3M 1s
    ##  86100K .......... .......... .......... .......... .......... 63% 79.5M 1s
    ##  86150K .......... .......... .......... .......... .......... 63% 98.4M 1s
    ##  86200K .......... .......... .......... .......... .......... 64%  121M 1s
    ##  86250K .......... .......... .......... .......... .......... 64%  128M 1s
    ##  86300K .......... .......... .......... .......... .......... 64% 78.8M 1s
    ##  86350K .......... .......... .......... .......... .......... 64% 88.2M 1s
    ##  86400K .......... .......... .......... .......... .......... 64%  157M 1s
    ##  86450K .......... .......... .......... .......... .......... 64% 96.9M 1s
    ##  86500K .......... .......... .......... .......... .......... 64%  100M 1s
    ##  86550K .......... .......... .......... .......... .......... 64%  121M 1s
    ##  86600K .......... .......... .......... .......... .......... 64% 87.4M 1s
    ##  86650K .......... .......... .......... .......... .......... 64%  100M 1s
    ##  86700K .......... .......... .......... .......... .......... 64% 78.3M 1s
    ##  86750K .......... .......... .......... .......... .......... 64%  106M 1s
    ##  86800K .......... .......... .......... .......... .......... 64% 88.3M 1s
    ##  86850K .......... .......... .......... .......... .......... 64% 92.1M 1s
    ##  86900K .......... .......... .......... .......... .......... 64%  122M 1s
    ##  86950K .......... .......... .......... .......... .......... 64%  126M 1s
    ##  87000K .......... .......... .......... .......... .......... 64% 90.2M 1s
    ##  87050K .......... .......... .......... .......... .......... 64%  135M 1s
    ##  87100K .......... .......... .......... .......... .......... 64% 86.6M 1s
    ##  87150K .......... .......... .......... .......... .......... 64% 93.7M 1s
    ##  87200K .......... .......... .......... .......... .......... 64% 98.5M 1s
    ##  87250K .......... .......... .......... .......... .......... 64% 91.4M 1s
    ##  87300K .......... .......... .......... .......... .......... 64%  108M 1s
    ##  87350K .......... .......... .......... .......... .......... 64%  147M 1s
    ##  87400K .......... .......... .......... .......... .......... 64%  100M 1s
    ##  87450K .......... .......... .......... .......... .......... 64% 87.4M 1s
    ##  87500K .......... .......... .......... .......... .......... 64%  125M 1s
    ##  87550K .......... .......... .......... .......... .......... 65% 40.0M 1s
    ##  87600K .......... .......... .......... .......... .......... 65%  139M 1s
    ##  87650K .......... .......... .......... .......... .......... 65% 22.6M 1s
    ##  87700K .......... .......... .......... .......... .......... 65%  152M 1s
    ##  87750K .......... .......... .......... .......... .......... 65%  122M 1s
    ##  87800K .......... .......... .......... .......... .......... 65%  146M 1s
    ##  87850K .......... .......... .......... .......... .......... 65%  159M 1s
    ##  87900K .......... .......... .......... .......... .......... 65%  150M 1s
    ##  87950K .......... .......... .......... .......... .......... 65%  165M 1s
    ##  88000K .......... .......... .......... .......... .......... 65%  140M 1s
    ##  88050K .......... .......... .......... .......... .......... 65% 36.5M 1s
    ##  88100K .......... .......... .......... .......... .......... 65%  127M 1s
    ##  88150K .......... .......... .......... .......... .......... 65%  130M 1s
    ##  88200K .......... .......... .......... .......... .......... 65%  130M 1s
    ##  88250K .......... .......... .......... .......... .......... 65%  133M 1s
    ##  88300K .......... .......... .......... .......... .......... 65%  160M 1s
    ##  88350K .......... .......... .......... .......... .......... 65%  163M 1s
    ##  88400K .......... .......... .......... .......... .......... 65%  138M 1s
    ##  88450K .......... .......... .......... .......... .......... 65% 36.0M 1s
    ##  88500K .......... .......... .......... .......... .......... 65% 45.9M 1s
    ##  88550K .......... .......... .......... .......... .......... 65% 55.9M 1s
    ##  88600K .......... .......... .......... .......... .......... 65% 47.5M 1s
    ##  88650K .......... .......... .......... .......... .......... 65% 83.7M 1s
    ##  88700K .......... .......... .......... .......... .......... 65%  147M 1s
    ##  88750K .......... .......... .......... .......... .......... 65%  143M 1s
    ##  88800K .......... .......... .......... .......... .......... 65%  150M 1s
    ##  88850K .......... .......... .......... .......... .......... 65%  151M 1s
    ##  88900K .......... .......... .......... .......... .......... 66% 15.1M 1s
    ##  88950K .......... .......... .......... .......... .......... 66%  144M 1s
    ##  89000K .......... .......... .......... .......... .......... 66%  159M 1s
    ##  89050K .......... .......... .......... .......... .......... 66%  143M 1s
    ##  89100K .......... .......... .......... .......... .......... 66%  152M 1s
    ##  89150K .......... .......... .......... .......... .......... 66%  127M 1s
    ##  89200K .......... .......... .......... .......... .......... 66% 20.6M 1s
    ##  89250K .......... .......... .......... .......... .......... 66%  149M 1s
    ##  89300K .......... .......... .......... .......... .......... 66%  156M 1s
    ##  89350K .......... .......... .......... .......... .......... 66% 84.7M 1s
    ##  89400K .......... .......... .......... .......... .......... 66%  155M 1s
    ##  89450K .......... .......... .......... .......... .......... 66% 88.2M 1s
    ##  89500K .......... .......... .......... .......... .......... 66%  128M 1s
    ##  89550K .......... .......... .......... .......... .......... 66%  155M 1s
    ##  89600K .......... .......... .......... .......... .......... 66% 18.3M 1s
    ##  89650K .......... .......... .......... .......... .......... 66%  142M 1s
    ##  89700K .......... .......... .......... .......... .......... 66%  123M 1s
    ##  89750K .......... .......... .......... .......... .......... 66%  110M 1s
    ##  89800K .......... .......... .......... .......... .......... 66%  130M 1s
    ##  89850K .......... .......... .......... .......... .......... 66%  111M 1s
    ##  89900K .......... .......... .......... .......... .......... 66%  119M 1s
    ##  89950K .......... .......... .......... .......... .......... 66%  130M 1s
    ##  90000K .......... .......... .......... .......... .......... 66% 25.0M 1s
    ##  90050K .......... .......... .......... .......... .......... 66% 70.2M 1s
    ##  90100K .......... .......... .......... .......... .......... 66% 53.2M 1s
    ##  90150K .......... .......... .......... .......... .......... 66% 52.1M 1s
    ##  90200K .......... .......... .......... .......... .......... 66%  141M 1s
    ##  90250K .......... .......... .......... .......... .......... 67%  114M 1s
    ##  90300K .......... .......... .......... .......... .......... 67%  144M 1s
    ##  90350K .......... .......... .......... .......... .......... 67%  158M 1s
    ##  90400K .......... .......... .......... .......... .......... 67%  145M 1s
    ##  90450K .......... .......... .......... .......... .......... 67%  154M 1s
    ##  90500K .......... .......... .......... .......... .......... 67%  136M 1s
    ##  90550K .......... .......... .......... .......... .......... 67% 65.9M 1s
    ##  90600K .......... .......... .......... .......... .......... 67% 64.2M 1s
    ##  90650K .......... .......... .......... .......... .......... 67% 61.5M 1s
    ##  90700K .......... .......... .......... .......... .......... 67% 44.8M 1s
    ##  90750K .......... .......... .......... .......... .......... 67% 86.0M 1s
    ##  90800K .......... .......... .......... .......... .......... 67%  121M 1s
    ##  90850K .......... .......... .......... .......... .......... 67%  140M 1s
    ##  90900K .......... .......... .......... .......... .......... 67% 62.9M 1s
    ##  90950K .......... .......... .......... .......... .......... 67%  101M 1s
    ##  91000K .......... .......... .......... .......... .......... 67%  122M 1s
    ##  91050K .......... .......... .......... .......... .......... 67% 46.0M 1s
    ##  91100K .......... .......... .......... .......... .......... 67%  124M 1s
    ##  91150K .......... .......... .......... .......... .......... 67% 59.7M 1s
    ##  91200K .......... .......... .......... .......... .......... 67%  115M 1s
    ##  91250K .......... .......... .......... .......... .......... 67% 40.0M 1s
    ##  91300K .......... .......... .......... .......... .......... 67% 73.7M 1s
    ##  91350K .......... .......... .......... .......... .......... 67%  129M 1s
    ##  91400K .......... .......... .......... .......... .......... 67%  148M 1s
    ##  91450K .......... .......... .......... .......... .......... 67%  149M 1s
    ##  91500K .......... .......... .......... .......... .......... 67%  114M 1s
    ##  91550K .......... .......... .......... .......... .......... 67% 90.9M 1s
    ##  91600K .......... .......... .......... .......... .......... 68% 32.4M 1s
    ##  91650K .......... .......... .......... .......... .......... 68% 88.3M 1s
    ##  91700K .......... .......... .......... .......... .......... 68%  116M 1s
    ##  91750K .......... .......... .......... .......... .......... 68%  133M 1s
    ##  91800K .......... .......... .......... .......... .......... 68%  139M 1s
    ##  91850K .......... .......... .......... .......... .......... 68%  126M 1s
    ##  91900K .......... .......... .......... .......... .......... 68%  117M 1s
    ##  91950K .......... .......... .......... .......... .......... 68%  110M 1s
    ##  92000K .......... .......... .......... .......... .......... 68%  124M 1s
    ##  92050K .......... .......... .......... .......... .......... 68% 45.4M 1s
    ##  92100K .......... .......... .......... .......... .......... 68%  113M 1s
    ##  92150K .......... .......... .......... .......... .......... 68%  121M 1s
    ##  92200K .......... .......... .......... .......... .......... 68%  147M 1s
    ##  92250K .......... .......... .......... .......... .......... 68%  162M 1s
    ##  92300K .......... .......... .......... .......... .......... 68%  111M 1s
    ##  92350K .......... .......... .......... .......... .......... 68% 57.4M 1s
    ##  92400K .......... .......... .......... .......... .......... 68%  111M 1s
    ##  92450K .......... .......... .......... .......... .......... 68% 58.7M 1s
    ##  92500K .......... .......... .......... .......... .......... 68% 86.3M 1s
    ##  92550K .......... .......... .......... .......... .......... 68%  125M 1s
    ##  92600K .......... .......... .......... .......... .......... 68% 64.5M 1s
    ##  92650K .......... .......... .......... .......... .......... 68%  140M 1s
    ##  92700K .......... .......... .......... .......... .......... 68%  135M 1s
    ##  92750K .......... .......... .......... .......... .......... 68%  107M 1s
    ##  92800K .......... .......... .......... .......... .......... 68%  135M 1s
    ##  92850K .......... .......... .......... .......... .......... 68%  144M 1s
    ##  92900K .......... .......... .......... .......... .......... 68% 21.3M 1s
    ##  92950K .......... .......... .......... .......... .......... 69%  139M 1s
    ##  93000K .......... .......... .......... .......... .......... 69%  155M 1s
    ##  93050K .......... .......... .......... .......... .......... 69%  161M 1s
    ##  93100K .......... .......... .......... .......... .......... 69%  105M 1s
    ##  93150K .......... .......... .......... .......... .......... 69%  172M 1s
    ##  93200K .......... .......... .......... .......... .......... 69%  116M 1s
    ##  93250K .......... .......... .......... .......... .......... 69%  149M 1s
    ##  93300K .......... .......... .......... .......... .......... 69% 53.0M 1s
    ##  93350K .......... .......... .......... .......... .......... 69% 52.3M 1s
    ##  93400K .......... .......... .......... .......... .......... 69%  126M 1s
    ##  93450K .......... .......... .......... .......... .......... 69% 72.9M 1s
    ##  93500K .......... .......... .......... .......... .......... 69%  153M 1s
    ##  93550K .......... .......... .......... .......... .......... 69%  149M 1s
    ##  93600K .......... .......... .......... .......... .......... 69%  120M 1s
    ##  93650K .......... .......... .......... .......... .......... 69%  174M 1s
    ##  93700K .......... .......... .......... .......... .......... 69% 64.4M 1s
    ##  93750K .......... .......... .......... .......... .......... 69% 52.5M 1s
    ##  93800K .......... .......... .......... .......... .......... 69%  102M 1s
    ##  93850K .......... .......... .......... .......... .......... 69%  113M 1s
    ##  93900K .......... .......... .......... .......... .......... 69% 26.6M 1s
    ##  93950K .......... .......... .......... .......... .......... 69%  108M 1s
    ##  94000K .......... .......... .......... .......... .......... 69%  133M 1s
    ##  94050K .......... .......... .......... .......... .......... 69%  167M 1s
    ##  94100K .......... .......... .......... .......... .......... 69%  153M 1s
    ##  94150K .......... .......... .......... .......... .......... 69%  151M 1s
    ##  94200K .......... .......... .......... .......... .......... 69%  152M 1s
    ##  94250K .......... .......... .......... .......... .......... 69%  165M 1s
    ##  94300K .......... .......... .......... .......... .......... 70% 21.4M 1s
    ##  94350K .......... .......... .......... .......... .......... 70% 78.6M 1s
    ##  94400K .......... .......... .......... .......... .......... 70% 56.5M 1s
    ##  94450K .......... .......... .......... .......... .......... 70% 85.3M 1s
    ##  94500K .......... .......... .......... .......... .......... 70% 87.6M 1s
    ##  94550K .......... .......... .......... .......... .......... 70%  107M 1s
    ##  94600K .......... .......... .......... .......... .......... 70%  113M 1s
    ##  94650K .......... .......... .......... .......... .......... 70%  150M 1s
    ##  94700K .......... .......... .......... .......... .......... 70% 86.7M 1s
    ##  94750K .......... .......... .......... .......... .......... 70%  113M 1s
    ##  94800K .......... .......... .......... .......... .......... 70%  128M 1s
    ##  94850K .......... .......... .......... .......... .......... 70%  129M 1s
    ##  94900K .......... .......... .......... .......... .......... 70%  126M 1s
    ##  94950K .......... .......... .......... .......... .......... 70%  152M 1s
    ##  95000K .......... .......... .......... .......... .......... 70%  122M 1s
    ##  95050K .......... .......... .......... .......... .......... 70%  139M 1s
    ##  95100K .......... .......... .......... .......... .......... 70% 56.9M 1s
    ##  95150K .......... .......... .......... .......... .......... 70%  116M 1s
    ##  95200K .......... .......... .......... .......... .......... 70%  125M 1s
    ##  95250K .......... .......... .......... .......... .......... 70% 64.1M 1s
    ##  95300K .......... .......... .......... .......... .......... 70% 79.3M 1s
    ##  95350K .......... .......... .......... .......... .......... 70% 84.4M 1s
    ##  95400K .......... .......... .......... .......... .......... 70% 49.5M 1s
    ##  95450K .......... .......... .......... .......... .......... 70%  170M 1s
    ##  95500K .......... .......... .......... .......... .......... 70% 76.5M 1s
    ##  95550K .......... .......... .......... .......... .......... 70% 84.3M 1s
    ##  95600K .......... .......... .......... .......... .......... 70%  133M 1s
    ##  95650K .......... .......... .......... .......... .......... 71%  148M 1s
    ##  95700K .......... .......... .......... .......... .......... 71% 52.1M 1s
    ##  95750K .......... .......... .......... .......... .......... 71%  164M 1s
    ##  95800K .......... .......... .......... .......... .......... 71%  107M 1s
    ##  95850K .......... .......... .......... .......... .......... 71%  140M 1s
    ##  95900K .......... .......... .......... .......... .......... 71% 95.1M 1s
    ##  95950K .......... .......... .......... .......... .......... 71%  127M 1s
    ##  96000K .......... .......... .......... .......... .......... 71% 99.5M 1s
    ##  96050K .......... .......... .......... .......... .......... 71% 40.0M 1s
    ##  96100K .......... .......... .......... .......... .......... 71%  121M 1s
    ##  96150K .......... .......... .......... .......... .......... 71%  161M 1s
    ##  96200K .......... .......... .......... .......... .......... 71% 85.7M 1s
    ##  96250K .......... .......... .......... .......... .......... 71% 90.1M 1s
    ##  96300K .......... .......... .......... .......... .......... 71% 57.6M 1s
    ##  96350K .......... .......... .......... .......... .......... 71%  161M 1s
    ##  96400K .......... .......... .......... .......... .......... 71% 85.2M 1s
    ##  96450K .......... .......... .......... .......... .......... 71%  148M 1s
    ##  96500K .......... .......... .......... .......... .......... 71%  123M 1s
    ##  96550K .......... .......... .......... .......... .......... 71%  145M 1s
    ##  96600K .......... .......... .......... .......... .......... 71% 94.2M 1s
    ##  96650K .......... .......... .......... .......... .......... 71% 79.0M 1s
    ##  96700K .......... .......... .......... .......... .......... 71% 88.5M 1s
    ##  96750K .......... .......... .......... .......... .......... 71%  144M 1s
    ##  96800K .......... .......... .......... .......... .......... 71% 53.6M 1s
    ##  96850K .......... .......... .......... .......... .......... 71%  133M 1s
    ##  96900K .......... .......... .......... .......... .......... 71% 65.2M 1s
    ##  96950K .......... .......... .......... .......... .......... 71% 68.5M 1s
    ##  97000K .......... .......... .......... .......... .......... 72%  121M 1s
    ##  97050K .......... .......... .......... .......... .......... 72%  127M 1s
    ##  97100K .......... .......... .......... .......... .......... 72%  157M 1s
    ##  97150K .......... .......... .......... .......... .......... 72%  154M 1s
    ##  97200K .......... .......... .......... .......... .......... 72%  125M 1s
    ##  97250K .......... .......... .......... .......... .......... 72%  151M 1s
    ##  97300K .......... .......... .......... .......... .......... 72%  104M 1s
    ##  97350K .......... .......... .......... .......... .......... 72%  117M 1s
    ##  97400K .......... .......... .......... .......... .......... 72%  107M 1s
    ##  97450K .......... .......... .......... .......... .......... 72% 50.3M 1s
    ##  97500K .......... .......... .......... .......... .......... 72%  131M 1s
    ##  97550K .......... .......... .......... .......... .......... 72%  167M 1s
    ##  97600K .......... .......... .......... .......... .......... 72% 79.4M 1s
    ##  97650K .......... .......... .......... .......... .......... 72% 59.4M 1s
    ##  97700K .......... .......... .......... .......... .......... 72% 76.5M 1s
    ##  97750K .......... .......... .......... .......... .......... 72% 84.8M 1s
    ##  97800K .......... .......... .......... .......... .......... 72%  155M 1s
    ##  97850K .......... .......... .......... .......... .......... 72% 48.4M 1s
    ##  97900K .......... .......... .......... .......... .......... 72% 58.5M 1s
    ##  97950K .......... .......... .......... .......... .......... 72%  181M 1s
    ##  98000K .......... .......... .......... .......... .......... 72%  151M 1s
    ##  98050K .......... .......... .......... .......... .......... 72%  117M 1s
    ##  98100K .......... .......... .......... .......... .......... 72%  154M 1s
    ##  98150K .......... .......... .......... .......... .......... 72%  152M 1s
    ##  98200K .......... .......... .......... .......... .......... 72% 49.4M 1s
    ##  98250K .......... .......... .......... .......... .......... 72% 50.8M 1s
    ##  98300K .......... .......... .......... .......... .......... 72%  144M 1s
    ##  98350K .......... .......... .......... .......... .......... 73% 97.0M 1s
    ##  98400K .......... .......... .......... .......... .......... 73%  110M 1s
    ##  98450K .......... .......... .......... .......... .......... 73%  182M 1s
    ##  98500K .......... .......... .......... .......... .......... 73% 51.1M 1s
    ##  98550K .......... .......... .......... .......... .......... 73%  152M 1s
    ##  98600K .......... .......... .......... .......... .......... 73%  129M 1s
    ##  98650K .......... .......... .......... .......... .......... 73% 93.6M 1s
    ##  98700K .......... .......... .......... .......... .......... 73%  117M 1s
    ##  98750K .......... .......... .......... .......... .......... 73% 64.4M 1s
    ##  98800K .......... .......... .......... .......... .......... 73%  100M 1s
    ##  98850K .......... .......... .......... .......... .......... 73% 67.4M 1s
    ##  98900K .......... .......... .......... .......... .......... 73%  158M 1s
    ##  98950K .......... .......... .......... .......... .......... 73%  129M 1s
    ##  99000K .......... .......... .......... .......... .......... 73%  163M 1s
    ##  99050K .......... .......... .......... .......... .......... 73%  119M 1s
    ##  99100K .......... .......... .......... .......... .......... 73%  128M 1s
    ##  99150K .......... .......... .......... .......... .......... 73% 60.9M 1s
    ##  99200K .......... .......... .......... .......... .......... 73% 24.5M 1s
    ##  99250K .......... .......... .......... .......... .......... 73% 96.2M 1s
    ##  99300K .......... .......... .......... .......... .......... 73%  158M 1s
    ##  99350K .......... .......... .......... .......... .......... 73%  152M 1s
    ##  99400K .......... .......... .......... .......... .......... 73%  184M 1s
    ##  99450K .......... .......... .......... .......... .......... 73% 99.6M 1s
    ##  99500K .......... .......... .......... .......... .......... 73%  119M 1s
    ##  99550K .......... .......... .......... .......... .......... 73%  174M 1s
    ##  99600K .......... .......... .......... .......... .......... 73%  142M 1s
    ##  99650K .......... .......... .......... .......... .......... 73% 39.5M 1s
    ##  99700K .......... .......... .......... .......... .......... 74%  150M 1s
    ##  99750K .......... .......... .......... .......... .......... 74% 48.3M 1s
    ##  99800K .......... .......... .......... .......... .......... 74% 89.4M 1s
    ##  99850K .......... .......... .......... .......... .......... 74%  102M 1s
    ##  99900K .......... .......... .......... .......... .......... 74%  108M 1s
    ##  99950K .......... .......... .......... .......... .......... 74% 72.1M 1s
    ## 100000K .......... .......... .......... .......... .......... 74%  147M 1s
    ## 100050K .......... .......... .......... .......... .......... 74%  143M 1s
    ## 100100K .......... .......... .......... .......... .......... 74%  133M 1s
    ## 100150K .......... .......... .......... .......... .......... 74%  130M 1s
    ## 100200K .......... .......... .......... .......... .......... 74% 61.8M 1s
    ## 100250K .......... .......... .......... .......... .......... 74% 79.6M 1s
    ## 100300K .......... .......... .......... .......... .......... 74% 90.9M 1s
    ## 100350K .......... .......... .......... .......... .......... 74% 52.8M 1s
    ## 100400K .......... .......... .......... .......... .......... 74% 79.7M 1s
    ## 100450K .......... .......... .......... .......... .......... 74%  132M 1s
    ## 100500K .......... .......... .......... .......... .......... 74%  124M 1s
    ## 100550K .......... .......... .......... .......... .......... 74%  145M 1s
    ## 100600K .......... .......... .......... .......... .......... 74% 83.1M 1s
    ## 100650K .......... .......... .......... .......... .......... 74%  121M 1s
    ## 100700K .......... .......... .......... .......... .......... 74%  106M 1s
    ## 100750K .......... .......... .......... .......... .......... 74% 70.3M 1s
    ## 100800K .......... .......... .......... .......... .......... 74%  129M 1s
    ## 100850K .......... .......... .......... .......... .......... 74% 57.6M 1s
    ## 100900K .......... .......... .......... .......... .......... 74%  132M 1s
    ## 100950K .......... .......... .......... .......... .......... 74% 55.3M 1s
    ## 101000K .......... .......... .......... .......... .......... 74% 91.2M 1s
    ## 101050K .......... .......... .......... .......... .......... 75%  137M 1s
    ## 101100K .......... .......... .......... .......... .......... 75% 70.3M 1s
    ## 101150K .......... .......... .......... .......... .......... 75%  139M 1s
    ## 101200K .......... .......... .......... .......... .......... 75%  107M 1s
    ## 101250K .......... .......... .......... .......... .......... 75%  147M 1s
    ## 101300K .......... .......... .......... .......... .......... 75%  131M 1s
    ## 101350K .......... .......... .......... .......... .......... 75% 99.4M 1s
    ## 101400K .......... .......... .......... .......... .......... 75% 90.0M 1s
    ## 101450K .......... .......... .......... .......... .......... 75% 94.6M 1s
    ## 101500K .......... .......... .......... .......... .......... 75%  123M 1s
    ## 101550K .......... .......... .......... .......... .......... 75%  115M 1s
    ## 101600K .......... .......... .......... .......... .......... 75% 88.1M 1s
    ## 101650K .......... .......... .......... .......... .......... 75% 95.1M 1s
    ## 101700K .......... .......... .......... .......... .......... 75%  127M 1s
    ## 101750K .......... .......... .......... .......... .......... 75%  114M 1s
    ## 101800K .......... .......... .......... .......... .......... 75% 82.8M 1s
    ## 101850K .......... .......... .......... .......... .......... 75% 84.0M 1s
    ## 101900K .......... .......... .......... .......... .......... 75% 72.1M 1s
    ## 101950K .......... .......... .......... .......... .......... 75%  109M 1s
    ## 102000K .......... .......... .......... .......... .......... 75%  174M 1s
    ## 102050K .......... .......... .......... .......... .......... 75%  120M 1s
    ## 102100K .......... .......... .......... .......... .......... 75% 63.4M 1s
    ## 102150K .......... .......... .......... .......... .......... 75%  123M 0s
    ## 102200K .......... .......... .......... .......... .......... 75%  136M 0s
    ## 102250K .......... .......... .......... .......... .......... 75% 88.9M 0s
    ## 102300K .......... .......... .......... .......... .......... 75% 21.6M 0s
    ## 102350K .......... .......... .......... .......... .......... 75% 30.1M 0s
    ## 102400K .......... .......... .......... .......... .......... 76%  126M 0s
    ## 102450K .......... .......... .......... .......... .......... 76%  120M 0s
    ## 102500K .......... .......... .......... .......... .......... 76%  126M 0s
    ## 102550K .......... .......... .......... .......... .......... 76% 99.5M 0s
    ## 102600K .......... .......... .......... .......... .......... 76%  135M 0s
    ## 102650K .......... .......... .......... .......... .......... 76%  162M 0s
    ## 102700K .......... .......... .......... .......... .......... 76% 77.2M 0s
    ## 102750K .......... .......... .......... .......... .......... 76% 31.8M 0s
    ## 102800K .......... .......... .......... .......... .......... 76% 96.7M 0s
    ## 102850K .......... .......... .......... .......... .......... 76% 60.2M 0s
    ## 102900K .......... .......... .......... .......... .......... 76%  112M 0s
    ## 102950K .......... .......... .......... .......... .......... 76% 76.7M 0s
    ## 103000K .......... .......... .......... .......... .......... 76%  125M 0s
    ## 103050K .......... .......... .......... .......... .......... 76%  177M 0s
    ## 103100K .......... .......... .......... .......... .......... 76%  116M 0s
    ## 103150K .......... .......... .......... .......... .......... 76%  154M 0s
    ## 103200K .......... .......... .......... .......... .......... 76%  148M 0s
    ## 103250K .......... .......... .......... .......... .......... 76% 53.1M 0s
    ## 103300K .......... .......... .......... .......... .......... 76% 75.3M 0s
    ## 103350K .......... .......... .......... .......... .......... 76% 58.8M 0s
    ## 103400K .......... .......... .......... .......... .......... 76% 71.9M 0s
    ## 103450K .......... .......... .......... .......... .......... 76%  101M 0s
    ## 103500K .......... .......... .......... .......... .......... 76% 65.2M 0s
    ## 103550K .......... .......... .......... .......... .......... 76% 67.9M 0s
    ## 103600K .......... .......... .......... .......... .......... 76% 40.4M 0s
    ## 103650K .......... .......... .......... .......... .......... 76%  191M 0s
    ## 103700K .......... .......... .......... .......... .......... 77%  110M 0s
    ## 103750K .......... .......... .......... .......... .......... 77%  131M 0s
    ## 103800K .......... .......... .......... .......... .......... 77%  140M 0s
    ## 103850K .......... .......... .......... .......... .......... 77%  192M 0s
    ## 103900K .......... .......... .......... .......... .......... 77%  163M 0s
    ## 103950K .......... .......... .......... .......... .......... 77% 63.9M 0s
    ## 104000K .......... .......... .......... .......... .......... 77%  127M 0s
    ## 104050K .......... .......... .......... .......... .......... 77% 23.3M 0s
    ## 104100K .......... .......... .......... .......... .......... 77% 94.7M 0s
    ## 104150K .......... .......... .......... .......... .......... 77%  161M 0s
    ## 104200K .......... .......... .......... .......... .......... 77% 25.9M 0s
    ## 104250K .......... .......... .......... .......... .......... 77%  127M 0s
    ## 104300K .......... .......... .......... .......... .......... 77%  115M 0s
    ## 104350K .......... .......... .......... .......... .......... 77%  123M 0s
    ## 104400K .......... .......... .......... .......... .......... 77%  118M 0s
    ## 104450K .......... .......... .......... .......... .......... 77%  164M 0s
    ## 104500K .......... .......... .......... .......... .......... 77% 37.1M 0s
    ## 104550K .......... .......... .......... .......... .......... 77%  148M 0s
    ## 104600K .......... .......... .......... .......... .......... 77%  130M 0s
    ## 104650K .......... .......... .......... .......... .......... 77% 32.9M 0s
    ## 104700K .......... .......... .......... .......... .......... 77%  119M 0s
    ## 104750K .......... .......... .......... .......... .......... 77%  116M 0s
    ## 104800K .......... .......... .......... .......... .......... 77%  128M 0s
    ## 104850K .......... .......... .......... .......... .......... 77%  171M 0s
    ## 104900K .......... .......... .......... .......... .......... 77%  151M 0s
    ## 104950K .......... .......... .......... .......... .......... 77%  160M 0s
    ## 105000K .......... .......... .......... .......... .......... 77%  147M 0s
    ## 105050K .......... .......... .......... .......... .......... 78% 63.0M 0s
    ## 105100K .......... .......... .......... .......... .......... 78% 41.9M 0s
    ## 105150K .......... .......... .......... .......... .......... 78% 98.5M 0s
    ## 105200K .......... .......... .......... .......... .......... 78% 76.8M 0s
    ## 105250K .......... .......... .......... .......... .......... 78%  160M 0s
    ## 105300K .......... .......... .......... .......... .......... 78%  159M 0s
    ## 105350K .......... .......... .......... .......... .......... 78% 5.66M 0s
    ## 105400K .......... .......... .......... .......... .......... 78% 91.0M 0s
    ## 105450K .......... .......... .......... .......... .......... 78% 45.1M 0s
    ## 105500K .......... .......... .......... .......... .......... 78% 72.3M 0s
    ## 105550K .......... .......... .......... .......... .......... 78%  136M 0s
    ## 105600K .......... .......... .......... .......... .......... 78% 41.1M 0s
    ## 105650K .......... .......... .......... .......... .......... 78%  124M 0s
    ## 105700K .......... .......... .......... .......... .......... 78%  146M 0s
    ## 105750K .......... .......... .......... .......... .......... 78%  164M 0s
    ## 105800K .......... .......... .......... .......... .......... 78%  104M 0s
    ## 105850K .......... .......... .......... .......... .......... 78%  157M 0s
    ## 105900K .......... .......... .......... .......... .......... 78%  154M 0s
    ## 105950K .......... .......... .......... .......... .......... 78%  145M 0s
    ## 106000K .......... .......... .......... .......... .......... 78%  142M 0s
    ## 106050K .......... .......... .......... .......... .......... 78% 22.7M 0s
    ## 106100K .......... .......... .......... .......... .......... 78% 56.5M 0s
    ## 106150K .......... .......... .......... .......... .......... 78% 49.7M 0s
    ## 106200K .......... .......... .......... .......... .......... 78% 63.7M 0s
    ## 106250K .......... .......... .......... .......... .......... 78%  143M 0s
    ## 106300K .......... .......... .......... .......... .......... 78% 88.3M 0s
    ## 106350K .......... .......... .......... .......... .......... 78%  153M 0s
    ## 106400K .......... .......... .......... .......... .......... 79%  134M 0s
    ## 106450K .......... .......... .......... .......... .......... 79%  135M 0s
    ## 106500K .......... .......... .......... .......... .......... 79% 54.9M 0s
    ## 106550K .......... .......... .......... .......... .......... 79%  157M 0s
    ## 106600K .......... .......... .......... .......... .......... 79% 67.8M 0s
    ## 106650K .......... .......... .......... .......... .......... 79%  130M 0s
    ## 106700K .......... .......... .......... .......... .......... 79%  109M 0s
    ## 106750K .......... .......... .......... .......... .......... 79% 80.8M 0s
    ## 106800K .......... .......... .......... .......... .......... 79% 20.6M 0s
    ## 106850K .......... .......... .......... .......... .......... 79% 78.7M 0s
    ## 106900K .......... .......... .......... .......... .......... 79%  115M 0s
    ## 106950K .......... .......... .......... .......... .......... 79%  139M 0s
    ## 107000K .......... .......... .......... .......... .......... 79%  145M 0s
    ## 107050K .......... .......... .......... .......... .......... 79%  112M 0s
    ## 107100K .......... .......... .......... .......... .......... 79%  129M 0s
    ## 107150K .......... .......... .......... .......... .......... 79%  151M 0s
    ## 107200K .......... .......... .......... .......... .......... 79%  154M 0s
    ## 107250K .......... .......... .......... .......... .......... 79% 23.9M 0s
    ## 107300K .......... .......... .......... .......... .......... 79% 35.3M 0s
    ## 107350K .......... .......... .......... .......... .......... 79%  192M 0s
    ## 107400K .......... .......... .......... .......... .......... 79% 27.8M 0s
    ## 107450K .......... .......... .......... .......... .......... 79% 44.6M 0s
    ## 107500K .......... .......... .......... .......... .......... 79% 40.0M 0s
    ## 107550K .......... .......... .......... .......... .......... 79%  119M 0s
    ## 107600K .......... .......... .......... .......... .......... 79% 88.1M 0s
    ## 107650K .......... .......... .......... .......... .......... 79%  122M 0s
    ## 107700K .......... .......... .......... .......... .......... 79% 62.1M 0s
    ## 107750K .......... .......... .......... .......... .......... 80%  114M 0s
    ## 107800K .......... .......... .......... .......... .......... 80% 72.1M 0s
    ## 107850K .......... .......... .......... .......... .......... 80%  126M 0s
    ## 107900K .......... .......... .......... .......... .......... 80%  114M 0s
    ## 107950K .......... .......... .......... .......... .......... 80%  130M 0s
    ## 108000K .......... .......... .......... .......... .......... 80% 88.3M 0s
    ## 108050K .......... .......... .......... .......... .......... 80% 89.5M 0s
    ## 108100K .......... .......... .......... .......... .......... 80%  110M 0s
    ## 108150K .......... .......... .......... .......... .......... 80%  102M 0s
    ## 108200K .......... .......... .......... .......... .......... 80% 81.8M 0s
    ## 108250K .......... .......... .......... .......... .......... 80%  104M 0s
    ## 108300K .......... .......... .......... .......... .......... 80% 98.7M 0s
    ## 108350K .......... .......... .......... .......... .......... 80%  113M 0s
    ## 108400K .......... .......... .......... .......... .......... 80% 98.8M 0s
    ## 108450K .......... .......... .......... .......... .......... 80%  108M 0s
    ## 108500K .......... .......... .......... .......... .......... 80%  110M 0s
    ## 108550K .......... .......... .......... .......... .......... 80%  124M 0s
    ## 108600K .......... .......... .......... .......... .......... 80% 89.0M 0s
    ## 108650K .......... .......... .......... .......... .......... 80%  102M 0s
    ## 108700K .......... .......... .......... .......... .......... 80% 41.7M 0s
    ## 108750K .......... .......... .......... .......... .......... 80%  110M 0s
    ## 108800K .......... .......... .......... .......... .......... 80%  106M 0s
    ## 108850K .......... .......... .......... .......... .......... 80%  102M 0s
    ## 108900K .......... .......... .......... .......... .......... 80%  113M 0s
    ## 108950K .......... .......... .......... .......... .......... 80%  102M 0s
    ## 109000K .......... .......... .......... .......... .......... 80% 88.6M 0s
    ## 109050K .......... .......... .......... .......... .......... 80%  111M 0s
    ## 109100K .......... .......... .......... .......... .......... 81% 79.3M 0s
    ## 109150K .......... .......... .......... .......... .......... 81% 92.4M 0s
    ## 109200K .......... .......... .......... .......... .......... 81%  103M 0s
    ## 109250K .......... .......... .......... .......... .......... 81%  100M 0s
    ## 109300K .......... .......... .......... .......... .......... 81%  101M 0s
    ## 109350K .......... .......... .......... .......... .......... 81%  113M 0s
    ## 109400K .......... .......... .......... .......... .......... 81% 89.2M 0s
    ## 109450K .......... .......... .......... .......... .......... 81%  109M 0s
    ## 109500K .......... .......... .......... .......... .......... 81%  106M 0s
    ## 109550K .......... .......... .......... .......... .......... 81% 88.3M 0s
    ## 109600K .......... .......... .......... .......... .......... 81% 94.8M 0s
    ## 109650K .......... .......... .......... .......... .......... 81% 96.7M 0s
    ## 109700K .......... .......... .......... .......... .......... 81% 58.1M 0s
    ## 109750K .......... .......... .......... .......... .......... 81%  114M 0s
    ## 109800K .......... .......... .......... .......... .......... 81%  116M 0s
    ## 109850K .......... .......... .......... .......... .......... 81%  113M 0s
    ## 109900K .......... .......... .......... .......... .......... 81%  108M 0s
    ## 109950K .......... .......... .......... .......... .......... 81%  147M 0s
    ## 110000K .......... .......... .......... .......... .......... 81%  109M 0s
    ## 110050K .......... .......... .......... .......... .......... 81% 82.8M 0s
    ## 110100K .......... .......... .......... .......... .......... 81%  124M 0s
    ## 110150K .......... .......... .......... .......... .......... 81%  140M 0s
    ## 110200K .......... .......... .......... .......... .......... 81% 93.3M 0s
    ## 110250K .......... .......... .......... .......... .......... 81%  130M 0s
    ## 110300K .......... .......... .......... .......... .......... 81% 95.7M 0s
    ## 110350K .......... .......... .......... .......... .......... 81%  108M 0s
    ## 110400K .......... .......... .......... .......... .......... 81% 97.8M 0s
    ## 110450K .......... .......... .......... .......... .......... 82%  125M 0s
    ## 110500K .......... .......... .......... .......... .......... 82% 90.1M 0s
    ## 110550K .......... .......... .......... .......... .......... 82%  117M 0s
    ## 110600K .......... .......... .......... .......... .......... 82% 52.4M 0s
    ## 110650K .......... .......... .......... .......... .......... 82%  127M 0s
    ## 110700K .......... .......... .......... .......... .......... 82%  108M 0s
    ## 110750K .......... .......... .......... .......... .......... 82% 99.5M 0s
    ## 110800K .......... .......... .......... .......... .......... 82%  106M 0s
    ## 110850K .......... .......... .......... .......... .......... 82%  146M 0s
    ## 110900K .......... .......... .......... .......... .......... 82% 83.5M 0s
    ## 110950K .......... .......... .......... .......... .......... 82%  121M 0s
    ## 111000K .......... .......... .......... .......... .......... 82% 89.4M 0s
    ## 111050K .......... .......... .......... .......... .......... 82%  137M 0s
    ## 111100K .......... .......... .......... .......... .......... 82%  129M 0s
    ## 111150K .......... .......... .......... .......... .......... 82% 86.7M 0s
    ## 111200K .......... .......... .......... .......... .......... 82% 70.3M 0s
    ## 111250K .......... .......... .......... .......... .......... 82% 39.1M 0s
    ## 111300K .......... .......... .......... .......... .......... 82%  105M 0s
    ## 111350K .......... .......... .......... .......... .......... 82%  138M 0s
    ## 111400K .......... .......... .......... .......... .......... 82% 66.9M 0s
    ## 111450K .......... .......... .......... .......... .......... 82% 86.1M 0s
    ## 111500K .......... .......... .......... .......... .......... 82%  100M 0s
    ## 111550K .......... .......... .......... .......... .......... 82%  129M 0s
    ## 111600K .......... .......... .......... .......... .......... 82%  130M 0s
    ## 111650K .......... .......... .......... .......... .......... 82%  101M 0s
    ## 111700K .......... .......... .......... .......... .......... 82%  143M 0s
    ## 111750K .......... .......... .......... .......... .......... 82% 45.0M 0s
    ## 111800K .......... .......... .......... .......... .......... 83% 77.9M 0s
    ## 111850K .......... .......... .......... .......... .......... 83%  103M 0s
    ## 111900K .......... .......... .......... .......... .......... 83% 96.1M 0s
    ## 111950K .......... .......... .......... .......... .......... 83% 67.1M 0s
    ## 112000K .......... .......... .......... .......... .......... 83%  132M 0s
    ## 112050K .......... .......... .......... .......... .......... 83%  119M 0s
    ## 112100K .......... .......... .......... .......... .......... 83% 78.9M 0s
    ## 112150K .......... .......... .......... .......... .......... 83% 93.6M 0s
    ## 112200K .......... .......... .......... .......... .......... 83% 28.3M 0s
    ## 112250K .......... .......... .......... .......... .......... 83% 87.4M 0s
    ## 112300K .......... .......... .......... .......... .......... 83% 86.7M 0s
    ## 112350K .......... .......... .......... .......... .......... 83% 90.0M 0s
    ## 112400K .......... .......... .......... .......... .......... 83% 52.9M 0s
    ## 112450K .......... .......... .......... .......... .......... 83%  102M 0s
    ## 112500K .......... .......... .......... .......... .......... 83%  112M 0s
    ## 112550K .......... .......... .......... .......... .......... 83% 50.3M 0s
    ## 112600K .......... .......... .......... .......... .......... 83% 76.4M 0s
    ## 112650K .......... .......... .......... .......... .......... 83%  118M 0s
    ## 112700K .......... .......... .......... .......... .......... 83%  122M 0s
    ## 112750K .......... .......... .......... .......... .......... 83% 72.4M 0s
    ## 112800K .......... .......... .......... .......... .......... 83% 92.5M 0s
    ## 112850K .......... .......... .......... .......... .......... 83% 93.3M 0s
    ## 112900K .......... .......... .......... .......... .......... 83%  122M 0s
    ## 112950K .......... .......... .......... .......... .......... 83% 88.1M 0s
    ## 113000K .......... .......... .......... .......... .......... 83%  125M 0s
    ## 113050K .......... .......... .......... .......... .......... 83%  109M 0s
    ## 113100K .......... .......... .......... .......... .......... 83% 95.8M 0s
    ## 113150K .......... .......... .......... .......... .......... 84% 43.7M 0s
    ## 113200K .......... .......... .......... .......... .......... 84% 87.0M 0s
    ## 113250K .......... .......... .......... .......... .......... 84%  120M 0s
    ## 113300K .......... .......... .......... .......... .......... 84%  103M 0s
    ## 113350K .......... .......... .......... .......... .......... 84% 69.9M 0s
    ## 113400K .......... .......... .......... .......... .......... 84% 88.6M 0s
    ## 113450K .......... .......... .......... .......... .......... 84%  127M 0s
    ## 113500K .......... .......... .......... .......... .......... 84%  102M 0s
    ## 113550K .......... .......... .......... .......... .......... 84%  142M 0s
    ## 113600K .......... .......... .......... .......... .......... 84% 58.7M 0s
    ## 113650K .......... .......... .......... .......... .......... 84% 78.7M 0s
    ## 113700K .......... .......... .......... .......... .......... 84% 94.4M 0s
    ## 113750K .......... .......... .......... .......... .......... 84%  131M 0s
    ## 113800K .......... .......... .......... .......... .......... 84% 47.7M 0s
    ## 113850K .......... .......... .......... .......... .......... 84% 98.8M 0s
    ## 113900K .......... .......... .......... .......... .......... 84% 73.3M 0s
    ## 113950K .......... .......... .......... .......... .......... 84%  108M 0s
    ## 114000K .......... .......... .......... .......... .......... 84% 59.0M 0s
    ## 114050K .......... .......... .......... .......... .......... 84% 84.3M 0s
    ## 114100K .......... .......... .......... .......... .......... 84% 89.7M 0s
    ## 114150K .......... .......... .......... .......... .......... 84%  107M 0s
    ## 114200K .......... .......... .......... .......... .......... 84% 78.8M 0s
    ## 114250K .......... .......... .......... .......... .......... 84%  110M 0s
    ## 114300K .......... .......... .......... .......... .......... 84%  121M 0s
    ## 114350K .......... .......... .......... .......... .......... 84% 87.4M 0s
    ## 114400K .......... .......... .......... .......... .......... 84%  140M 0s
    ## 114450K .......... .......... .......... .......... .......... 84%  119M 0s
    ## 114500K .......... .......... .......... .......... .......... 85% 96.9M 0s
    ## 114550K .......... .......... .......... .......... .......... 85%  103M 0s
    ## 114600K .......... .......... .......... .......... .......... 85% 80.6M 0s
    ## 114650K .......... .......... .......... .......... .......... 85% 61.0M 0s
    ## 114700K .......... .......... .......... .......... .......... 85% 92.8M 0s
    ## 114750K .......... .......... .......... .......... .......... 85%  108M 0s
    ## 114800K .......... .......... .......... .......... .......... 85% 82.8M 0s
    ## 114850K .......... .......... .......... .......... .......... 85%  120M 0s
    ## 114900K .......... .......... .......... .......... .......... 85% 82.5M 0s
    ## 114950K .......... .......... .......... .......... .......... 85%  119M 0s
    ## 115000K .......... .......... .......... .......... .......... 85%  138M 0s
    ## 115050K .......... .......... .......... .......... .......... 85% 87.5M 0s
    ## 115100K .......... .......... .......... .......... .......... 85%  102M 0s
    ## 115150K .......... .......... .......... .......... .......... 85% 90.9M 0s
    ## 115200K .......... .......... .......... .......... .......... 85%  118M 0s
    ## 115250K .......... .......... .......... .......... .......... 85% 56.0M 0s
    ## 115300K .......... .......... .......... .......... .......... 85%  130M 0s
    ## 115350K .......... .......... .......... .......... .......... 85% 58.8M 0s
    ## 115400K .......... .......... .......... .......... .......... 85%  111M 0s
    ## 115450K .......... .......... .......... .......... .......... 85%  161M 0s
    ## 115500K .......... .......... .......... .......... .......... 85% 84.3M 0s
    ## 115550K .......... .......... .......... .......... .......... 85%  131M 0s
    ## 115600K .......... .......... .......... .......... .......... 85%  122M 0s
    ## 115650K .......... .......... .......... .......... .......... 85%  134M 0s
    ## 115700K .......... .......... .......... .......... .......... 85% 58.7M 0s
    ## 115750K .......... .......... .......... .......... .......... 85%  113M 0s
    ## 115800K .......... .......... .......... .......... .......... 85%  120M 0s
    ## 115850K .......... .......... .......... .......... .......... 86%  106M 0s
    ## 115900K .......... .......... .......... .......... .......... 86%  138M 0s
    ## 115950K .......... .......... .......... .......... .......... 86% 67.7M 0s
    ## 116000K .......... .......... .......... .......... .......... 86% 70.3M 0s
    ## 116050K .......... .......... .......... .......... .......... 86% 94.5M 0s
    ## 116100K .......... .......... .......... .......... .......... 86% 92.9M 0s
    ## 116150K .......... .......... .......... .......... .......... 86%  130M 0s
    ## 116200K .......... .......... .......... .......... .......... 86% 93.0M 0s
    ## 116250K .......... .......... .......... .......... .......... 86%  127M 0s
    ## 116300K .......... .......... .......... .......... .......... 86%  120M 0s
    ## 116350K .......... .......... .......... .......... .......... 86%  124M 0s
    ## 116400K .......... .......... .......... .......... .......... 86%  130M 0s
    ## 116450K .......... .......... .......... .......... .......... 86%  104M 0s
    ## 116500K .......... .......... .......... .......... .......... 86%  100M 0s
    ## 116550K .......... .......... .......... .......... .......... 86% 70.5M 0s
    ## 116600K .......... .......... .......... .......... .......... 86%  106M 0s
    ## 116650K .......... .......... .......... .......... .......... 86%  110M 0s
    ## 116700K .......... .......... .......... .......... .......... 86%  108M 0s
    ## 116750K .......... .......... .......... .......... .......... 86%  134M 0s
    ## 116800K .......... .......... .......... .......... .......... 86%  107M 0s
    ## 116850K .......... .......... .......... .......... .......... 86%  128M 0s
    ## 116900K .......... .......... .......... .......... .......... 86%  148M 0s
    ## 116950K .......... .......... .......... .......... .......... 86%  104M 0s
    ## 117000K .......... .......... .......... .......... .......... 86% 22.2M 0s
    ## 117050K .......... .......... .......... .......... .......... 86%  125M 0s
    ## 117100K .......... .......... .......... .......... .......... 86%  105M 0s
    ## 117150K .......... .......... .......... .......... .......... 86%  161M 0s
    ## 117200K .......... .......... .......... .......... .......... 87%  119M 0s
    ## 117250K .......... .......... .......... .......... .......... 87%  104M 0s
    ## 117300K .......... .......... .......... .......... .......... 87%  123M 0s
    ## 117350K .......... .......... .......... .......... .......... 87%  133M 0s
    ## 117400K .......... .......... .......... .......... .......... 87%  147M 0s
    ## 117450K .......... .......... .......... .......... .......... 87%  120M 0s
    ## 117500K .......... .......... .......... .......... .......... 87% 62.4M 0s
    ## 117550K .......... .......... .......... .......... .......... 87% 39.0M 0s
    ## 117600K .......... .......... .......... .......... .......... 87% 96.0M 0s
    ## 117650K .......... .......... .......... .......... .......... 87%  111M 0s
    ## 117700K .......... .......... .......... .......... .......... 87%  124M 0s
    ## 117750K .......... .......... .......... .......... .......... 87%  157M 0s
    ## 117800K .......... .......... .......... .......... .......... 87%  113M 0s
    ## 117850K .......... .......... .......... .......... .......... 87%  126M 0s
    ## 117900K .......... .......... .......... .......... .......... 87%  107M 0s
    ## 117950K .......... .......... .......... .......... .......... 87%  105M 0s
    ## 118000K .......... .......... .......... .......... .......... 87%  152M 0s
    ## 118050K .......... .......... .......... .......... .......... 87%  163M 0s
    ## 118100K .......... .......... .......... .......... .......... 87% 19.0M 0s
    ## 118150K .......... .......... .......... .......... .......... 87% 77.4M 0s
    ## 118200K .......... .......... .......... .......... .......... 87% 89.6M 0s
    ## 118250K .......... .......... .......... .......... .......... 87%  152M 0s
    ## 118300K .......... .......... .......... .......... .......... 87%  141M 0s
    ## 118350K .......... .......... .......... .......... .......... 87%  146M 0s
    ## 118400K .......... .......... .......... .......... .......... 87% 98.4M 0s
    ## 118450K .......... .......... .......... .......... .......... 87%  140M 0s
    ## 118500K .......... .......... .......... .......... .......... 87% 95.6M 0s
    ## 118550K .......... .......... .......... .......... .......... 88% 63.5M 0s
    ## 118600K .......... .......... .......... .......... .......... 88% 67.5M 0s
    ## 118650K .......... .......... .......... .......... .......... 88% 64.6M 0s
    ## 118700K .......... .......... .......... .......... .......... 88% 61.4M 0s
    ## 118750K .......... .......... .......... .......... .......... 88% 87.6M 0s
    ## 118800K .......... .......... .......... .......... .......... 88% 82.4M 0s
    ## 118850K .......... .......... .......... .......... .......... 88% 56.7M 0s
    ## 118900K .......... .......... .......... .......... .......... 88% 40.8M 0s
    ## 118950K .......... .......... .......... .......... .......... 88% 93.0M 0s
    ## 119000K .......... .......... .......... .......... .......... 88% 84.0M 0s
    ## 119050K .......... .......... .......... .......... .......... 88% 69.1M 0s
    ## 119100K .......... .......... .......... .......... .......... 88% 63.0M 0s
    ## 119150K .......... .......... .......... .......... .......... 88% 84.6M 0s
    ## 119200K .......... .......... .......... .......... .......... 88% 65.9M 0s
    ## 119250K .......... .......... .......... .......... .......... 88% 88.2M 0s
    ## 119300K .......... .......... .......... .......... .......... 88% 90.6M 0s
    ## 119350K .......... .......... .......... .......... .......... 88% 63.2M 0s
    ## 119400K .......... .......... .......... .......... .......... 88% 87.6M 0s
    ## 119450K .......... .......... .......... .......... .......... 88% 71.8M 0s
    ## 119500K .......... .......... .......... .......... .......... 88% 81.8M 0s
    ## 119550K .......... .......... .......... .......... .......... 88% 93.0M 0s
    ## 119600K .......... .......... .......... .......... .......... 88% 92.1M 0s
    ## 119650K .......... .......... .......... .......... .......... 88% 42.4M 0s
    ## 119700K .......... .......... .......... .......... .......... 88% 79.5M 0s
    ## 119750K .......... .......... .......... .......... .......... 88%  107M 0s
    ## 119800K .......... .......... .......... .......... .......... 88% 95.1M 0s
    ## 119850K .......... .......... .......... .......... .......... 88% 73.1M 0s
    ## 119900K .......... .......... .......... .......... .......... 89% 84.0M 0s
    ## 119950K .......... .......... .......... .......... .......... 89% 78.9M 0s
    ## 120000K .......... .......... .......... .......... .......... 89% 92.3M 0s
    ## 120050K .......... .......... .......... .......... .......... 89%  102M 0s
    ## 120100K .......... .......... .......... .......... .......... 89% 56.9M 0s
    ## 120150K .......... .......... .......... .......... .......... 89% 64.2M 0s
    ## 120200K .......... .......... .......... .......... .......... 89% 94.8M 0s
    ## 120250K .......... .......... .......... .......... .......... 89%  106M 0s
    ## 120300K .......... .......... .......... .......... .......... 89%  106M 0s
    ## 120350K .......... .......... .......... .......... .......... 89% 36.5M 0s
    ## 120400K .......... .......... .......... .......... .......... 89% 89.3M 0s
    ## 120450K .......... .......... .......... .......... .......... 89% 58.5M 0s
    ## 120500K .......... .......... .......... .......... .......... 89% 83.0M 0s
    ## 120550K .......... .......... .......... .......... .......... 89% 89.0M 0s
    ## 120600K .......... .......... .......... .......... .......... 89% 58.4M 0s
    ## 120650K .......... .......... .......... .......... .......... 89% 71.1M 0s
    ## 120700K .......... .......... .......... .......... .......... 89% 72.0M 0s
    ## 120750K .......... .......... .......... .......... .......... 89% 73.0M 0s
    ## 120800K .......... .......... .......... .......... .......... 89% 50.5M 0s
    ## 120850K .......... .......... .......... .......... .......... 89% 92.0M 0s
    ## 120900K .......... .......... .......... .......... .......... 89%  105M 0s
    ## 120950K .......... .......... .......... .......... .......... 89% 81.6M 0s
    ## 121000K .......... .......... .......... .......... .......... 89% 88.3M 0s
    ## 121050K .......... .......... .......... .......... .......... 89% 74.5M 0s
    ## 121100K .......... .......... .......... .......... .......... 89% 85.8M 0s
    ## 121150K .......... .......... .......... .......... .......... 89%  132M 0s
    ## 121200K .......... .......... .......... .......... .......... 89% 72.4M 0s
    ## 121250K .......... .......... .......... .......... .......... 90% 94.6M 0s
    ## 121300K .......... .......... .......... .......... .......... 90% 78.8M 0s
    ## 121350K .......... .......... .......... .......... .......... 90% 86.3M 0s
    ## 121400K .......... .......... .......... .......... .......... 90%  111M 0s
    ## 121450K .......... .......... .......... .......... .......... 90% 90.4M 0s
    ## 121500K .......... .......... .......... .......... .......... 90% 95.8M 0s
    ## 121550K .......... .......... .......... .......... .......... 90%  124M 0s
    ## 121600K .......... .......... .......... .......... .......... 90% 99.9M 0s
    ## 121650K .......... .......... .......... .......... .......... 90% 84.0M 0s
    ## 121700K .......... .......... .......... .......... .......... 90% 88.8M 0s
    ## 121750K .......... .......... .......... .......... .......... 90%  102M 0s
    ## 121800K .......... .......... .......... .......... .......... 90% 58.3M 0s
    ## 121850K .......... .......... .......... .......... .......... 90% 89.2M 0s
    ## 121900K .......... .......... .......... .......... .......... 90% 89.6M 0s
    ## 121950K .......... .......... .......... .......... .......... 90% 82.6M 0s
    ## 122000K .......... .......... .......... .......... .......... 90% 96.0M 0s
    ## 122050K .......... .......... .......... .......... .......... 90% 99.0M 0s
    ## 122100K .......... .......... .......... .......... .......... 90%  106M 0s
    ## 122150K .......... .......... .......... .......... .......... 90%  106M 0s
    ## 122200K .......... .......... .......... .......... .......... 90% 97.2M 0s
    ## 122250K .......... .......... .......... .......... .......... 90% 96.5M 0s
    ## 122300K .......... .......... .......... .......... .......... 90%  106M 0s
    ## 122350K .......... .......... .......... .......... .......... 90% 95.4M 0s
    ## 122400K .......... .......... .......... .......... .......... 90% 92.4M 0s
    ## 122450K .......... .......... .......... .......... .......... 90% 94.0M 0s
    ## 122500K .......... .......... .......... .......... .......... 90% 55.9M 0s
    ## 122550K .......... .......... .......... .......... .......... 90% 91.9M 0s
    ## 122600K .......... .......... .......... .......... .......... 91%  105M 0s
    ## 122650K .......... .......... .......... .......... .......... 91%  102M 0s
    ## 122700K .......... .......... .......... .......... .......... 91% 66.6M 0s
    ## 122750K .......... .......... .......... .......... .......... 91%  127M 0s
    ## 122800K .......... .......... .......... .......... .......... 91%  108M 0s
    ## 122850K .......... .......... .......... .......... .......... 91%  134M 0s
    ## 122900K .......... .......... .......... .......... .......... 91% 92.8M 0s
    ## 122950K .......... .......... .......... .......... .......... 91% 86.5M 0s
    ## 123000K .......... .......... .......... .......... .......... 91% 75.6M 0s
    ## 123050K .......... .......... .......... .......... .......... 91% 94.9M 0s
    ## 123100K .......... .......... .......... .......... .......... 91% 93.0M 0s
    ## 123150K .......... .......... .......... .......... .......... 91%  100M 0s
    ## 123200K .......... .......... .......... .......... .......... 91% 75.4M 0s
    ## 123250K .......... .......... .......... .......... .......... 91%  140M 0s
    ## 123300K .......... .......... .......... .......... .......... 91% 98.4M 0s
    ## 123350K .......... .......... .......... .......... .......... 91% 90.5M 0s
    ## 123400K .......... .......... .......... .......... .......... 91%  101M 0s
    ## 123450K .......... .......... .......... .......... .......... 91%  128M 0s
    ## 123500K .......... .......... .......... .......... .......... 91% 79.4M 0s
    ## 123550K .......... .......... .......... .......... .......... 91%  138M 0s
    ## 123600K .......... .......... .......... .......... .......... 91% 58.0M 0s
    ## 123650K .......... .......... .......... .......... .......... 91%  100M 0s
    ## 123700K .......... .......... .......... .......... .......... 91%  103M 0s
    ## 123750K .......... .......... .......... .......... .......... 91%  125M 0s
    ## 123800K .......... .......... .......... .......... .......... 91% 89.6M 0s
    ## 123850K .......... .......... .......... .......... .......... 91%  156M 0s
    ## 123900K .......... .......... .......... .......... .......... 91%  112M 0s
    ## 123950K .......... .......... .......... .......... .......... 92%  112M 0s
    ## 124000K .......... .......... .......... .......... .......... 92% 93.4M 0s
    ## 124050K .......... .......... .......... .......... .......... 92%  125M 0s
    ## 124100K .......... .......... .......... .......... .......... 92%  130M 0s
    ## 124150K .......... .......... .......... .......... .......... 92%  126M 0s
    ## 124200K .......... .......... .......... .......... .......... 92%  116M 0s
    ## 124250K .......... .......... .......... .......... .......... 92% 89.5M 0s
    ## 124300K .......... .......... .......... .......... .......... 92% 87.3M 0s
    ## 124350K .......... .......... .......... .......... .......... 92%  142M 0s
    ## 124400K .......... .......... .......... .......... .......... 92% 77.5M 0s
    ## 124450K .......... .......... .......... .......... .......... 92%  125M 0s
    ## 124500K .......... .......... .......... .......... .......... 92%  117M 0s
    ## 124550K .......... .......... .......... .......... .......... 92% 93.3M 0s
    ## 124600K .......... .......... .......... .......... .......... 92% 78.4M 0s
    ## 124650K .......... .......... .......... .......... .......... 92%  130M 0s
    ## 124700K .......... .......... .......... .......... .......... 92%  109M 0s
    ## 124750K .......... .......... .......... .......... .......... 92%  125M 0s
    ## 124800K .......... .......... .......... .......... .......... 92%  107M 0s
    ## 124850K .......... .......... .......... .......... .......... 92%  111M 0s
    ## 124900K .......... .......... .......... .......... .......... 92% 97.9M 0s
    ## 124950K .......... .......... .......... .......... .......... 92%  105M 0s
    ## 125000K .......... .......... .......... .......... .......... 92%  109M 0s
    ## 125050K .......... .......... .......... .......... .......... 92% 95.3M 0s
    ## 125100K .......... .......... .......... .......... .......... 92%  103M 0s
    ## 125150K .......... .......... .......... .......... .......... 92%  108M 0s
    ## 125200K .......... .......... .......... .......... .......... 92%  121M 0s
    ## 125250K .......... .......... .......... .......... .......... 92%  109M 0s
    ## 125300K .......... .......... .......... .......... .......... 93%  125M 0s
    ## 125350K .......... .......... .......... .......... .......... 93%  135M 0s
    ## 125400K .......... .......... .......... .......... .......... 93% 88.5M 0s
    ## 125450K .......... .......... .......... .......... .......... 93%  106M 0s
    ## 125500K .......... .......... .......... .......... .......... 93%  128M 0s
    ## 125550K .......... .......... .......... .......... .......... 93% 93.5M 0s
    ## 125600K .......... .......... .......... .......... .......... 93%  109M 0s
    ## 125650K .......... .......... .......... .......... .......... 93% 99.2M 0s
    ## 125700K .......... .......... .......... .......... .......... 93%  127M 0s
    ## 125750K .......... .......... .......... .......... .......... 93% 81.4M 0s
    ## 125800K .......... .......... .......... .......... .......... 93%  136M 0s
    ## 125850K .......... .......... .......... .......... .......... 93%  129M 0s
    ## 125900K .......... .......... .......... .......... .......... 93% 96.5M 0s
    ## 125950K .......... .......... .......... .......... .......... 93%  139M 0s
    ## 126000K .......... .......... .......... .......... .......... 93%  104M 0s
    ## 126050K .......... .......... .......... .......... .......... 93%  130M 0s
    ## 126100K .......... .......... .......... .......... .......... 93%  106M 0s
    ## 126150K .......... .......... .......... .......... .......... 93%  128M 0s
    ## 126200K .......... .......... .......... .......... .......... 93%  129M 0s
    ## 126250K .......... .......... .......... .......... .......... 93%  109M 0s
    ## 126300K .......... .......... .......... .......... .......... 93%  139M 0s
    ## 126350K .......... .......... .......... .......... .......... 93%  157M 0s
    ## 126400K .......... .......... .......... .......... .......... 93% 89.4M 0s
    ## 126450K .......... .......... .......... .......... .......... 93%  153M 0s
    ## 126500K .......... .......... .......... .......... .......... 93%  126M 0s
    ## 126550K .......... .......... .......... .......... .......... 93% 90.4M 0s
    ## 126600K .......... .......... .......... .......... .......... 93%  146M 0s
    ## 126650K .......... .......... .......... .......... .......... 94%  127M 0s
    ## 126700K .......... .......... .......... .......... .......... 94%  105M 0s
    ## 126750K .......... .......... .......... .......... .......... 94%  112M 0s
    ## 126800K .......... .......... .......... .......... .......... 94%  104M 0s
    ## 126850K .......... .......... .......... .......... .......... 94%  151M 0s
    ## 126900K .......... .......... .......... .......... .......... 94% 99.4M 0s
    ## 126950K .......... .......... .......... .......... .......... 94%  176M 0s
    ## 127000K .......... .......... .......... .......... .......... 94%  115M 0s
    ## 127050K .......... .......... .......... .......... .......... 94% 91.8M 0s
    ## 127100K .......... .......... .......... .......... .......... 94% 95.6M 0s
    ## 127150K .......... .......... .......... .......... .......... 94%  158M 0s
    ## 127200K .......... .......... .......... .......... .......... 94%  113M 0s
    ## 127250K .......... .......... .......... .......... .......... 94%  136M 0s
    ## 127300K .......... .......... .......... .......... .......... 94%  125M 0s
    ## 127350K .......... .......... .......... .......... .......... 94%  133M 0s
    ## 127400K .......... .......... .......... .......... .......... 94%  109M 0s
    ## 127450K .......... .......... .......... .......... .......... 94%  172M 0s
    ## 127500K .......... .......... .......... .......... .......... 94%  129M 0s
    ## 127550K .......... .......... .......... .......... .......... 94%  135M 0s
    ## 127600K .......... .......... .......... .......... .......... 94%  108M 0s
    ## 127650K .......... .......... .......... .......... .......... 94% 94.7M 0s
    ## 127700K .......... .......... .......... .......... .......... 94%  116M 0s
    ## 127750K .......... .......... .......... .......... .......... 94%  126M 0s
    ## 127800K .......... .......... .......... .......... .......... 94%  123M 0s
    ## 127850K .......... .......... .......... .......... .......... 94%  130M 0s
    ## 127900K .......... .......... .......... .......... .......... 94% 96.4M 0s
    ## 127950K .......... .......... .......... .......... .......... 94%  112M 0s
    ## 128000K .......... .......... .......... .......... .......... 95%  114M 0s
    ## 128050K .......... .......... .......... .......... .......... 95%  145M 0s
    ## 128100K .......... .......... .......... .......... .......... 95%  117M 0s
    ## 128150K .......... .......... .......... .......... .......... 95%  175M 0s
    ## 128200K .......... .......... .......... .......... .......... 95%  105M 0s
    ## 128250K .......... .......... .......... .......... .......... 95%  126M 0s
    ## 128300K .......... .......... .......... .......... .......... 95%  123M 0s
    ## 128350K .......... .......... .......... .......... .......... 95%  141M 0s
    ## 128400K .......... .......... .......... .......... .......... 95%  126M 0s
    ## 128450K .......... .......... .......... .......... .......... 95%  173M 0s
    ## 128500K .......... .......... .......... .......... .......... 95%  108M 0s
    ## 128550K .......... .......... .......... .......... .......... 95%  132M 0s
    ## 128600K .......... .......... .......... .......... .......... 95%  118M 0s
    ## 128650K .......... .......... .......... .......... .......... 95%  143M 0s
    ## 128700K .......... .......... .......... .......... .......... 95%  139M 0s
    ## 128750K .......... .......... .......... .......... .......... 95%  137M 0s
    ## 128800K .......... .......... .......... .......... .......... 95%  134M 0s
    ## 128850K .......... .......... .......... .......... .......... 95%  152M 0s
    ## 128900K .......... .......... .......... .......... .......... 95% 79.2M 0s
    ## 128950K .......... .......... .......... .......... .......... 95%  111M 0s
    ## 129000K .......... .......... .......... .......... .......... 95%  112M 0s
    ## 129050K .......... .......... .......... .......... .......... 95%  152M 0s
    ## 129100K .......... .......... .......... .......... .......... 95% 97.9M 0s
    ## 129150K .......... .......... .......... .......... .......... 95%  147M 0s
    ## 129200K .......... .......... .......... .......... .......... 95%  124M 0s
    ## 129250K .......... .......... .......... .......... .......... 95%  110M 0s
    ## 129300K .......... .......... .......... .......... .......... 95%  109M 0s
    ## 129350K .......... .......... .......... .......... .......... 96%  124M 0s
    ## 129400K .......... .......... .......... .......... .......... 96%  101M 0s
    ## 129450K .......... .......... .......... .......... .......... 96%  147M 0s
    ## 129500K .......... .......... .......... .......... .......... 96%  117M 0s
    ## 129550K .......... .......... .......... .......... .......... 96%  131M 0s
    ## 129600K .......... .......... .......... .......... .......... 96%  121M 0s
    ## 129650K .......... .......... .......... .......... .......... 96%  152M 0s
    ## 129700K .......... .......... .......... .......... .......... 96%  142M 0s
    ## 129750K .......... .......... .......... .......... .......... 96%  122M 0s
    ## 129800K .......... .......... .......... .......... .......... 96%  130M 0s
    ## 129850K .......... .......... .......... .......... .......... 96%  144M 0s
    ## 129900K .......... .......... .......... .......... .......... 96%  107M 0s
    ## 129950K .......... .......... .......... .......... .......... 96%  159M 0s
    ## 130000K .......... .......... .......... .......... .......... 96%  128M 0s
    ## 130050K .......... .......... .......... .......... .......... 96%  151M 0s
    ## 130100K .......... .......... .......... .......... .......... 96%  122M 0s
    ## 130150K .......... .......... .......... .......... .......... 96%  118M 0s
    ## 130200K .......... .......... .......... .......... .......... 96%  106M 0s
    ## 130250K .......... .......... .......... .......... .......... 96% 8.60M 0s
    ## 130300K .......... .......... .......... .......... .......... 96% 91.1M 0s
    ## 130350K .......... .......... .......... .......... .......... 96%  112M 0s
    ## 130400K .......... .......... .......... .......... .......... 96%  117M 0s
    ## 130450K .......... .......... .......... .......... .......... 96%  137M 0s
    ## 130500K .......... .......... .......... .......... .......... 96%  147M 0s
    ## 130550K .......... .......... .......... .......... .......... 96%  199M 0s
    ## 130600K .......... .......... .......... .......... .......... 96%  138M 0s
    ## 130650K .......... .......... .......... .......... .......... 97%  139M 0s
    ## 130700K .......... .......... .......... .......... .......... 97%  149M 0s
    ## 130750K .......... .......... .......... .......... .......... 97%  183M 0s
    ## 130800K .......... .......... .......... .......... .......... 97% 42.1M 0s
    ## 130850K .......... .......... .......... .......... .......... 97%  124M 0s
    ## 130900K .......... .......... .......... .......... .......... 97% 40.9M 0s
    ## 130950K .......... .......... .......... .......... .......... 97% 89.9M 0s
    ## 131000K .......... .......... .......... .......... .......... 97% 69.3M 0s
    ## 131050K .......... .......... .......... .......... .......... 97%  176M 0s
    ## 131100K .......... .......... .......... .......... .......... 97%  126M 0s
    ## 131150K .......... .......... .......... .......... .......... 97% 35.2M 0s
    ## 131200K .......... .......... .......... .......... .......... 97%  126M 0s
    ## 131250K .......... .......... .......... .......... .......... 97%  152M 0s
    ## 131300K .......... .......... .......... .......... .......... 97%  136M 0s
    ## 131350K .......... .......... .......... .......... .......... 97%  111M 0s
    ## 131400K .......... .......... .......... .......... .......... 97% 89.2M 0s
    ## 131450K .......... .......... .......... .......... .......... 97%  149M 0s
    ## 131500K .......... .......... .......... .......... .......... 97%  123M 0s
    ## 131550K .......... .......... .......... .......... .......... 97%  119M 0s
    ## 131600K .......... .......... .......... .......... .......... 97%  107M 0s
    ## 131650K .......... .......... .......... .......... .......... 97%  153M 0s
    ## 131700K .......... .......... .......... .......... .......... 97%  131M 0s
    ## 131750K .......... .......... .......... .......... .......... 97%  132M 0s
    ## 131800K .......... .......... .......... .......... .......... 97%  130M 0s
    ## 131850K .......... .......... .......... .......... .......... 97% 51.2M 0s
    ## 131900K .......... .......... .......... .......... .......... 97% 63.5M 0s
    ## 131950K .......... .......... .......... .......... .......... 97% 29.4M 0s
    ## 132000K .......... .......... .......... .......... .......... 98% 95.9M 0s
    ## 132050K .......... .......... .......... .......... .......... 98%  148M 0s
    ## 132100K .......... .......... .......... .......... .......... 98% 89.8M 0s
    ## 132150K .......... .......... .......... .......... .......... 98%  124M 0s
    ## 132200K .......... .......... .......... .......... .......... 98% 81.5M 0s
    ## 132250K .......... .......... .......... .......... .......... 98%  183M 0s
    ## 132300K .......... .......... .......... .......... .......... 98%  119M 0s
    ## 132350K .......... .......... .......... .......... .......... 98%  106M 0s
    ## 132400K .......... .......... .......... .......... .......... 98% 80.5M 0s
    ## 132450K .......... .......... .......... .......... .......... 98%  101M 0s
    ## 132500K .......... .......... .......... .......... .......... 98% 37.2M 0s
    ## 132550K .......... .......... .......... .......... .......... 98%  100M 0s
    ## 132600K .......... .......... .......... .......... .......... 98% 73.5M 0s
    ## 132650K .......... .......... .......... .......... .......... 98%  111M 0s
    ## 132700K .......... .......... .......... .......... .......... 98%  107M 0s
    ## 132750K .......... .......... .......... .......... .......... 98%  118M 0s
    ## 132800K .......... .......... .......... .......... .......... 98% 77.2M 0s
    ## 132850K .......... .......... .......... .......... .......... 98%  120M 0s
    ## 132900K .......... .......... .......... .......... .......... 98%  121M 0s
    ## 132950K .......... .......... .......... .......... .......... 98% 86.2M 0s
    ## 133000K .......... .......... .......... .......... .......... 98%  109M 0s
    ## 133050K .......... .......... .......... .......... .......... 98% 59.8M 0s
    ## 133100K .......... .......... .......... .......... .......... 98%  115M 0s
    ## 133150K .......... .......... .......... .......... .......... 98% 72.0M 0s
    ## 133200K .......... .......... .......... .......... .......... 98% 65.9M 0s
    ## 133250K .......... .......... .......... .......... .......... 98% 95.9M 0s
    ## 133300K .......... .......... .......... .......... .......... 98%  113M 0s
    ## 133350K .......... .......... .......... .......... .......... 99%  137M 0s
    ## 133400K .......... .......... .......... .......... .......... 99% 51.8M 0s
    ## 133450K .......... .......... .......... .......... .......... 99%  146M 0s
    ## 133500K .......... .......... .......... .......... .......... 99%  112M 0s
    ## 133550K .......... .......... .......... .......... .......... 99%  132M 0s
    ## 133600K .......... .......... .......... .......... .......... 99%  137M 0s
    ## 133650K .......... .......... .......... .......... .......... 99% 90.1M 0s
    ## 133700K .......... .......... .......... .......... .......... 99% 88.1M 0s
    ## 133750K .......... .......... .......... .......... .......... 99% 99.4M 0s
    ## 133800K .......... .......... .......... .......... .......... 99%  146M 0s
    ## 133850K .......... .......... .......... .......... .......... 99% 57.3M 0s
    ## 133900K .......... .......... .......... .......... .......... 99% 73.2M 0s
    ## 133950K .......... .......... .......... .......... .......... 99% 99.3M 0s
    ## 134000K .......... .......... .......... .......... .......... 99%  108M 0s
    ## 134050K .......... .......... .......... .......... .......... 99%  170M 0s
    ## 134100K .......... .......... .......... .......... .......... 99% 91.2M 0s
    ## 134150K .......... .......... .......... .......... .......... 99% 68.7M 0s
    ## 134200K .......... .......... .......... .......... .......... 99%  102M 0s
    ## 134250K .......... .......... .......... .......... .......... 99%  118M 0s
    ## 134300K .......... .......... .......... .......... .......... 99%  117M 0s
    ## 134350K .......... .......... .......... .......... .......... 99%  125M 0s
    ## 134400K .......... .......... .......... .......... .......... 99% 50.1M 0s
    ## 134450K .......... .......... .......... .......... .......... 99%  127M 0s
    ## 134500K .......... .......... .......... .......... .......... 99%  120M 0s
    ## 134550K .......... .......... .......... .......... .......... 99%  108M 0s
    ## 134600K .......... .......... .......... .......... .......... 99%  102M 0s
    ## 134650K .......... .......... .......... .......... .......... 99%  193M 0s
    ## 134700K .......... .......... .......... ..........           100%  174M=1.9s
    ## 
    ## 2020-12-04 08:44:00 (68.1 MB/s) - ‘silva_nr99_v138_train_set.fa.gz.6’ saved [137973851/137973851]

\#Assignation taxonomique

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

## Assignation taxonomique n°2 Silva species assignement

``` bash
wget https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
```

    ## --2020-12-04 08:46:32--  https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 81840166 (78M) [application/octet-stream]
    ## Saving to: ‘silva_species_assignment_v138.fa.gz.3’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 7.80M 10s
    ##     50K .......... .......... .......... .......... ..........  0% 14.8M 8s
    ##    100K .......... .......... .......... .......... ..........  0% 5.24M 10s
    ##    150K .......... .......... .......... .......... ..........  0% 12.2M 9s
    ##    200K .......... .......... .......... .......... ..........  0% 68.2M 8s
    ##    250K .......... .......... .......... .......... ..........  0% 67.5M 6s
    ##    300K .......... .......... .......... .......... ..........  0% 87.9M 6s
    ##    350K .......... .......... .......... .......... ..........  0% 23.6M 5s
    ##    400K .......... .......... .......... .......... ..........  0% 69.5M 5s
    ##    450K .......... .......... .......... .......... ..........  0% 70.7M 4s
    ##    500K .......... .......... .......... .......... ..........  0% 82.0M 4s
    ##    550K .......... .......... .......... .......... ..........  0% 63.7M 4s
    ##    600K .......... .......... .......... .......... ..........  0% 64.2M 4s
    ##    650K .......... .......... .......... .......... ..........  0% 84.5M 4s
    ##    700K .......... .......... .......... .......... ..........  0% 80.8M 3s
    ##    750K .......... .......... .......... .......... ..........  1% 87.9M 3s
    ##    800K .......... .......... .......... .......... ..........  1% 67.5M 3s
    ##    850K .......... .......... .......... .......... ..........  1% 91.1M 3s
    ##    900K .......... .......... .......... .......... ..........  1% 74.7M 3s
    ##    950K .......... .......... .......... .......... ..........  1% 66.8M 3s
    ##   1000K .......... .......... .......... .......... ..........  1% 54.0M 3s
    ##   1050K .......... .......... .......... .......... ..........  1% 59.3M 3s
    ##   1100K .......... .......... .......... .......... ..........  1% 86.6M 3s
    ##   1150K .......... .......... .......... .......... ..........  1% 17.4M 3s
    ##   1200K .......... .......... .......... .......... ..........  1% 82.8M 3s
    ##   1250K .......... .......... .......... .......... ..........  1%  104M 2s
    ##   1300K .......... .......... .......... .......... ..........  1% 80.6M 2s
    ##   1350K .......... .......... .......... .......... ..........  1%  100M 2s
    ##   1400K .......... .......... .......... .......... ..........  1% 83.2M 2s
    ##   1450K .......... .......... .......... .......... ..........  1% 98.0M 2s
    ##   1500K .......... .......... .......... .......... ..........  1% 81.9M 2s
    ##   1550K .......... .......... .......... .......... ..........  2% 81.7M 2s
    ##   1600K .......... .......... .......... .......... ..........  2% 75.9M 2s
    ##   1650K .......... .......... .......... .......... ..........  2% 28.5M 2s
    ##   1700K .......... .......... .......... .......... ..........  2% 79.7M 2s
    ##   1750K .......... .......... .......... .......... ..........  2% 92.7M 2s
    ##   1800K .......... .......... .......... .......... ..........  2% 96.9M 2s
    ##   1850K .......... .......... .......... .......... ..........  2% 80.4M 2s
    ##   1900K .......... .......... .......... .......... ..........  2% 65.9M 2s
    ##   1950K .......... .......... .......... .......... ..........  2% 99.9M 2s
    ##   2000K .......... .......... .......... .......... ..........  2% 74.9M 2s
    ##   2050K .......... .......... .......... .......... ..........  2% 95.4M 2s
    ##   2100K .......... .......... .......... .......... ..........  2% 59.4M 2s
    ##   2150K .......... .......... .......... .......... ..........  2% 72.5M 2s
    ##   2200K .......... .......... .......... .......... ..........  2% 77.2M 2s
    ##   2250K .......... .......... .......... .......... ..........  2% 84.4M 2s
    ##   2300K .......... .......... .......... .......... ..........  2% 61.6M 2s
    ##   2350K .......... .......... .......... .......... ..........  3% 96.6M 2s
    ##   2400K .......... .......... .......... .......... ..........  3% 84.8M 2s
    ##   2450K .......... .......... .......... .......... ..........  3% 81.0M 2s
    ##   2500K .......... .......... .......... .......... ..........  3% 69.8M 2s
    ##   2550K .......... .......... .......... .......... ..........  3% 79.3M 2s
    ##   2600K .......... .......... .......... .......... ..........  3% 47.3M 2s
    ##   2650K .......... .......... .......... .......... ..........  3% 88.2M 2s
    ##   2700K .......... .......... .......... .......... ..........  3% 91.0M 2s
    ##   2750K .......... .......... .......... .......... ..........  3% 65.4M 2s
    ##   2800K .......... .......... .......... .......... ..........  3% 56.8M 2s
    ##   2850K .......... .......... .......... .......... ..........  3% 95.9M 2s
    ##   2900K .......... .......... .......... .......... ..........  3% 82.7M 2s
    ##   2950K .......... .......... .......... .......... ..........  3% 74.6M 2s
    ##   3000K .......... .......... .......... .......... ..........  3% 67.7M 2s
    ##   3050K .......... .......... .......... .......... ..........  3% 75.8M 2s
    ##   3100K .......... .......... .......... .......... ..........  3% 73.7M 2s
    ##   3150K .......... .......... .......... .......... ..........  4% 65.0M 2s
    ##   3200K .......... .......... .......... .......... ..........  4% 63.5M 2s
    ##   3250K .......... .......... .......... .......... ..........  4% 55.1M 2s
    ##   3300K .......... .......... .......... .......... ..........  4% 97.2M 2s
    ##   3350K .......... .......... .......... .......... ..........  4% 81.4M 2s
    ##   3400K .......... .......... .......... .......... ..........  4% 76.0M 2s
    ##   3450K .......... .......... .......... .......... ..........  4% 95.9M 2s
    ##   3500K .......... .......... .......... .......... ..........  4% 96.0M 2s
    ##   3550K .......... .......... .......... .......... ..........  4% 85.5M 2s
    ##   3600K .......... .......... .......... .......... ..........  4% 75.7M 1s
    ##   3650K .......... .......... .......... .......... ..........  4%  101M 1s
    ##   3700K .......... .......... .......... .......... ..........  4% 43.3M 1s
    ##   3750K .......... .......... .......... .......... ..........  4% 73.7M 1s
    ##   3800K .......... .......... .......... .......... ..........  4% 87.6M 1s
    ##   3850K .......... .......... .......... .......... ..........  4%  109M 1s
    ##   3900K .......... .......... .......... .......... ..........  4% 88.3M 1s
    ##   3950K .......... .......... .......... .......... ..........  5% 94.6M 1s
    ##   4000K .......... .......... .......... .......... ..........  5% 84.2M 1s
    ##   4050K .......... .......... .......... .......... ..........  5%  103M 1s
    ##   4100K .......... .......... .......... .......... ..........  5% 68.7M 1s
    ##   4150K .......... .......... .......... .......... ..........  5%  108M 1s
    ##   4200K .......... .......... .......... .......... ..........  5% 52.2M 1s
    ##   4250K .......... .......... .......... .......... ..........  5% 56.4M 1s
    ##   4300K .......... .......... .......... .......... ..........  5% 30.5M 1s
    ##   4350K .......... .......... .......... .......... ..........  5% 90.6M 1s
    ##   4400K .......... .......... .......... .......... ..........  5% 92.1M 1s
    ##   4450K .......... .......... .......... .......... ..........  5%  122M 1s
    ##   4500K .......... .......... .......... .......... ..........  5% 94.8M 1s
    ##   4550K .......... .......... .......... .......... ..........  5%  100M 1s
    ##   4600K .......... .......... .......... .......... ..........  5% 37.3M 1s
    ##   4650K .......... .......... .......... .......... ..........  5%  104M 1s
    ##   4700K .......... .......... .......... .......... ..........  5% 76.4M 1s
    ##   4750K .......... .......... .......... .......... ..........  6% 97.7M 1s
    ##   4800K .......... .......... .......... .......... ..........  6% 75.7M 1s
    ##   4850K .......... .......... .......... .......... ..........  6% 87.5M 1s
    ##   4900K .......... .......... .......... .......... ..........  6% 71.1M 1s
    ##   4950K .......... .......... .......... .......... ..........  6% 87.7M 1s
    ##   5000K .......... .......... .......... .......... ..........  6% 82.8M 1s
    ##   5050K .......... .......... .......... .......... ..........  6%  123M 1s
    ##   5100K .......... .......... .......... .......... ..........  6% 73.0M 1s
    ##   5150K .......... .......... .......... .......... ..........  6% 67.3M 1s
    ##   5200K .......... .......... .......... .......... ..........  6% 56.1M 1s
    ##   5250K .......... .......... .......... .......... ..........  6%  116M 1s
    ##   5300K .......... .......... .......... .......... ..........  6% 94.0M 1s
    ##   5350K .......... .......... .......... .......... ..........  6%  125M 1s
    ##   5400K .......... .......... .......... .......... ..........  6% 43.0M 1s
    ##   5450K .......... .......... .......... .......... ..........  6% 98.3M 1s
    ##   5500K .......... .......... .......... .......... ..........  6% 85.5M 1s
    ##   5550K .......... .......... .......... .......... ..........  7%  110M 1s
    ##   5600K .......... .......... .......... .......... ..........  7% 88.4M 1s
    ##   5650K .......... .......... .......... .......... ..........  7%  113M 1s
    ##   5700K .......... .......... .......... .......... ..........  7% 65.3M 1s
    ##   5750K .......... .......... .......... .......... ..........  7% 95.6M 1s
    ##   5800K .......... .......... .......... .......... ..........  7% 92.3M 1s
    ##   5850K .......... .......... .......... .......... ..........  7% 93.8M 1s
    ##   5900K .......... .......... .......... .......... ..........  7% 27.9M 1s
    ##   5950K .......... .......... .......... .......... ..........  7% 20.8M 1s
    ##   6000K .......... .......... .......... .......... ..........  7% 54.7M 1s
    ##   6050K .......... .......... .......... .......... ..........  7% 68.6M 1s
    ##   6100K .......... .......... .......... .......... ..........  7% 57.7M 1s
    ##   6150K .......... .......... .......... .......... ..........  7% 74.3M 1s
    ##   6200K .......... .......... .......... .......... ..........  7% 84.0M 1s
    ##   6250K .......... .......... .......... .......... ..........  7%  115M 1s
    ##   6300K .......... .......... .......... .......... ..........  7% 91.7M 1s
    ##   6350K .......... .......... .......... .......... ..........  8% 88.6M 1s
    ##   6400K .......... .......... .......... .......... ..........  8% 79.5M 1s
    ##   6450K .......... .......... .......... .......... ..........  8% 79.4M 1s
    ##   6500K .......... .......... .......... .......... ..........  8% 75.6M 1s
    ##   6550K .......... .......... .......... .......... ..........  8%  111M 1s
    ##   6600K .......... .......... .......... .......... ..........  8% 91.0M 1s
    ##   6650K .......... .......... .......... .......... ..........  8% 71.8M 1s
    ##   6700K .......... .......... .......... .......... ..........  8% 85.7M 1s
    ##   6750K .......... .......... .......... .......... ..........  8% 73.2M 1s
    ##   6800K .......... .......... .......... .......... ..........  8% 91.8M 1s
    ##   6850K .......... .......... .......... .......... ..........  8%  109M 1s
    ##   6900K .......... .......... .......... .......... ..........  8% 79.8M 1s
    ##   6950K .......... .......... .......... .......... ..........  8% 86.1M 1s
    ##   7000K .......... .......... .......... .......... ..........  8% 79.0M 1s
    ##   7050K .......... .......... .......... .......... ..........  8% 82.8M 1s
    ##   7100K .......... .......... .......... .......... ..........  8% 89.9M 1s
    ##   7150K .......... .......... .......... .......... ..........  9% 93.0M 1s
    ##   7200K .......... .......... .......... .......... ..........  9% 65.0M 1s
    ##   7250K .......... .......... .......... .......... ..........  9%  101M 1s
    ##   7300K .......... .......... .......... .......... ..........  9% 80.2M 1s
    ##   7350K .......... .......... .......... .......... ..........  9%  113M 1s
    ##   7400K .......... .......... .......... .......... ..........  9% 82.5M 1s
    ##   7450K .......... .......... .......... .......... ..........  9%  115M 1s
    ##   7500K .......... .......... .......... .......... ..........  9% 85.1M 1s
    ##   7550K .......... .......... .......... .......... ..........  9% 76.4M 1s
    ##   7600K .......... .......... .......... .......... ..........  9% 93.3M 1s
    ##   7650K .......... .......... .......... .......... ..........  9%  117M 1s
    ##   7700K .......... .......... .......... .......... ..........  9%  105M 1s
    ##   7750K .......... .......... .......... .......... ..........  9%  127M 1s
    ##   7800K .......... .......... .......... .......... ..........  9%  100M 1s
    ##   7850K .......... .......... .......... .......... ..........  9%  101M 1s
    ##   7900K .......... .......... .......... .......... ..........  9% 92.0M 1s
    ##   7950K .......... .......... .......... .......... .......... 10%  116M 1s
    ##   8000K .......... .......... .......... .......... .......... 10%  106M 1s
    ##   8050K .......... .......... .......... .......... .......... 10%  113M 1s
    ##   8100K .......... .......... .......... .......... .......... 10% 83.2M 1s
    ##   8150K .......... .......... .......... .......... .......... 10%  114M 1s
    ##   8200K .......... .......... .......... .......... .......... 10% 89.6M 1s
    ##   8250K .......... .......... .......... .......... .......... 10% 78.3M 1s
    ##   8300K .......... .......... .......... .......... .......... 10%  112M 1s
    ##   8350K .......... .......... .......... .......... .......... 10%  135M 1s
    ##   8400K .......... .......... .......... .......... .......... 10%  115M 1s
    ##   8450K .......... .......... .......... .......... .......... 10%  135M 1s
    ##   8500K .......... .......... .......... .......... .......... 10%  103M 1s
    ##   8550K .......... .......... .......... .......... .......... 10%  114M 1s
    ##   8600K .......... .......... .......... .......... .......... 10%  124M 1s
    ##   8650K .......... .......... .......... .......... .......... 10%  144M 1s
    ##   8700K .......... .......... .......... .......... .......... 10% 95.6M 1s
    ##   8750K .......... .......... .......... .......... .......... 11% 96.7M 1s
    ##   8800K .......... .......... .......... .......... .......... 11%  103M 1s
    ##   8850K .......... .......... .......... .......... .......... 11%  135M 1s
    ##   8900K .......... .......... .......... .......... .......... 11%  117M 1s
    ##   8950K .......... .......... .......... .......... .......... 11%  136M 1s
    ##   9000K .......... .......... .......... .......... .......... 11%  113M 1s
    ##   9050K .......... .......... .......... .......... .......... 11%  116M 1s
    ##   9100K .......... .......... .......... .......... .......... 11%  101M 1s
    ##   9150K .......... .......... .......... .......... .......... 11%  108M 1s
    ##   9200K .......... .......... .......... .......... .......... 11% 97.0M 1s
    ##   9250K .......... .......... .......... .......... .......... 11%  144M 1s
    ##   9300K .......... .......... .......... .......... .......... 11%  114M 1s
    ##   9350K .......... .......... .......... .......... .......... 11% 86.4M 1s
    ##   9400K .......... .......... .......... .......... .......... 11%  106M 1s
    ##   9450K .......... .......... .......... .......... .......... 11%  135M 1s
    ##   9500K .......... .......... .......... .......... .......... 11%  108M 1s
    ##   9550K .......... .......... .......... .......... .......... 12%  111M 1s
    ##   9600K .......... .......... .......... .......... .......... 12%  107M 1s
    ##   9650K .......... .......... .......... .......... .......... 12%  101M 1s
    ##   9700K .......... .......... .......... .......... .......... 12%  116M 1s
    ##   9750K .......... .......... .......... .......... .......... 12%  123M 1s
    ##   9800K .......... .......... .......... .......... .......... 12%  101M 1s
    ##   9850K .......... .......... .......... .......... .......... 12%  126M 1s
    ##   9900K .......... .......... .......... .......... .......... 12%  116M 1s
    ##   9950K .......... .......... .......... .......... .......... 12%  109M 1s
    ##  10000K .......... .......... .......... .......... .......... 12%  103M 1s
    ##  10050K .......... .......... .......... .......... .......... 12%  123M 1s
    ##  10100K .......... .......... .......... .......... .......... 12% 97.2M 1s
    ##  10150K .......... .......... .......... .......... .......... 12%  134M 1s
    ##  10200K .......... .......... .......... .......... .......... 12%  114M 1s
    ##  10250K .......... .......... .......... .......... .......... 12%  118M 1s
    ##  10300K .......... .......... .......... .......... .......... 12%  113M 1s
    ##  10350K .......... .......... .......... .......... .......... 13%  137M 1s
    ##  10400K .......... .......... .......... .......... .......... 13% 56.3M 1s
    ##  10450K .......... .......... .......... .......... .......... 13%  136M 1s
    ##  10500K .......... .......... .......... .......... .......... 13% 96.3M 1s
    ##  10550K .......... .......... .......... .......... .......... 13%  112M 1s
    ##  10600K .......... .......... .......... .......... .......... 13% 98.7M 1s
    ##  10650K .......... .......... .......... .......... .......... 13%  140M 1s
    ##  10700K .......... .......... .......... .......... .......... 13% 88.9M 1s
    ##  10750K .......... .......... .......... .......... .......... 13%  121M 1s
    ##  10800K .......... .......... .......... .......... .......... 13%  109M 1s
    ##  10850K .......... .......... .......... .......... .......... 13% 93.2M 1s
    ##  10900K .......... .......... .......... .......... .......... 13% 96.8M 1s
    ##  10950K .......... .......... .......... .......... .......... 13%  127M 1s
    ##  11000K .......... .......... .......... .......... .......... 13%  124M 1s
    ##  11050K .......... .......... .......... .......... .......... 13%  145M 1s
    ##  11100K .......... .......... .......... .......... .......... 13%  109M 1s
    ##  11150K .......... .......... .......... .......... .......... 14% 92.3M 1s
    ##  11200K .......... .......... .......... .......... .......... 14%  101M 1s
    ##  11250K .......... .......... .......... .......... .......... 14%  123M 1s
    ##  11300K .......... .......... .......... .......... .......... 14%  120M 1s
    ##  11350K .......... .......... .......... .......... .......... 14%  130M 1s
    ##  11400K .......... .......... .......... .......... .......... 14% 81.6M 1s
    ##  11450K .......... .......... .......... .......... .......... 14%  125M 1s
    ##  11500K .......... .......... .......... .......... .......... 14%  109M 1s
    ##  11550K .......... .......... .......... .......... .......... 14%  130M 1s
    ##  11600K .......... .......... .......... .......... .......... 14%  137M 1s
    ##  11650K .......... .......... .......... .......... .......... 14%  134M 1s
    ##  11700K .......... .......... .......... .......... .......... 14% 87.8M 1s
    ##  11750K .......... .......... .......... .......... .......... 14%  144M 1s
    ##  11800K .......... .......... .......... .......... .......... 14%  109M 1s
    ##  11850K .......... .......... .......... .......... .......... 14%  114M 1s
    ##  11900K .......... .......... .......... .......... .......... 14%  113M 1s
    ##  11950K .......... .......... .......... .......... .......... 15%  150M 1s
    ##  12000K .......... .......... .......... .......... .......... 15%  107M 1s
    ##  12050K .......... .......... .......... .......... .......... 15% 98.4M 1s
    ##  12100K .......... .......... .......... .......... .......... 15%  141M 1s
    ##  12150K .......... .......... .......... .......... .......... 15%  126M 1s
    ##  12200K .......... .......... .......... .......... .......... 15%  116M 1s
    ##  12250K .......... .......... .......... .......... .......... 15%  101M 1s
    ##  12300K .......... .......... .......... .......... .......... 15%  118M 1s
    ##  12350K .......... .......... .......... .......... .......... 15%  124M 1s
    ##  12400K .......... .......... .......... .......... .......... 15% 99.5M 1s
    ##  12450K .......... .......... .......... .......... .......... 15%  126M 1s
    ##  12500K .......... .......... .......... .......... .......... 15%  119M 1s
    ##  12550K .......... .......... .......... .......... .......... 15%  107M 1s
    ##  12600K .......... .......... .......... .......... .......... 15% 85.0M 1s
    ##  12650K .......... .......... .......... .......... .......... 15%  120M 1s
    ##  12700K .......... .......... .......... .......... .......... 15%  132M 1s
    ##  12750K .......... .......... .......... .......... .......... 16%  114M 1s
    ##  12800K .......... .......... .......... .......... .......... 16% 94.5M 1s
    ##  12850K .......... .......... .......... .......... .......... 16%  153M 1s
    ##  12900K .......... .......... .......... .......... .......... 16%  103M 1s
    ##  12950K .......... .......... .......... .......... .......... 16%  125M 1s
    ##  13000K .......... .......... .......... .......... .......... 16%  152M 1s
    ##  13050K .......... .......... .......... .......... .......... 16%  122M 1s
    ##  13100K .......... .......... .......... .......... .......... 16% 97.1M 1s
    ##  13150K .......... .......... .......... .......... .......... 16%  143M 1s
    ##  13200K .......... .......... .......... .......... .......... 16%  116M 1s
    ##  13250K .......... .......... .......... .......... .......... 16%  118M 1s
    ##  13300K .......... .......... .......... .......... .......... 16%  128M 1s
    ##  13350K .......... .......... .......... .......... .......... 16%  112M 1s
    ##  13400K .......... .......... .......... .......... .......... 16%  133M 1s
    ##  13450K .......... .......... .......... .......... .......... 16%  113M 1s
    ##  13500K .......... .......... .......... .......... .......... 16%  116M 1s
    ##  13550K .......... .......... .......... .......... .......... 17% 99.7M 1s
    ##  13600K .......... .......... .......... .......... .......... 17%  135M 1s
    ##  13650K .......... .......... .......... .......... .......... 17%  119M 1s
    ##  13700K .......... .......... .......... .......... .......... 17%  120M 1s
    ##  13750K .......... .......... .......... .......... .......... 17%  128M 1s
    ##  13800K .......... .......... .......... .......... .......... 17%  101M 1s
    ##  13850K .......... .......... .......... .......... .......... 17%  142M 1s
    ##  13900K .......... .......... .......... .......... .......... 17%  159M 1s
    ##  13950K .......... .......... .......... .......... .......... 17%  107M 1s
    ##  14000K .......... .......... .......... .......... .......... 17%  122M 1s
    ##  14050K .......... .......... .......... .......... .......... 17%  149M 1s
    ##  14100K .......... .......... .......... .......... .......... 17%  102M 1s
    ##  14150K .......... .......... .......... .......... .......... 17%  133M 1s
    ##  14200K .......... .......... .......... .......... .......... 17%  134M 1s
    ##  14250K .......... .......... .......... .......... .......... 17%  127M 1s
    ##  14300K .......... .......... .......... .......... .......... 17% 98.7M 1s
    ##  14350K .......... .......... .......... .......... .......... 18%  141M 1s
    ##  14400K .......... .......... .......... .......... .......... 18%  114M 1s
    ##  14450K .......... .......... .......... .......... .......... 18%  139M 1s
    ##  14500K .......... .......... .......... .......... .......... 18%  106M 1s
    ##  14550K .......... .......... .......... .......... .......... 18%  128M 1s
    ##  14600K .......... .......... .......... .......... .......... 18%  114M 1s
    ##  14650K .......... .......... .......... .......... .......... 18%  129M 1s
    ##  14700K .......... .......... .......... .......... .......... 18%  123M 1s
    ##  14750K .......... .......... .......... .......... .......... 18%  119M 1s
    ##  14800K .......... .......... .......... .......... .......... 18%  140M 1s
    ##  14850K .......... .......... .......... .......... .......... 18% 95.3M 1s
    ##  14900K .......... .......... .......... .......... .......... 18%  124M 1s
    ##  14950K .......... .......... .......... .......... .......... 18%  133M 1s
    ##  15000K .......... .......... .......... .......... .......... 18%  106M 1s
    ##  15050K .......... .......... .......... .......... .......... 18%  124M 1s
    ##  15100K .......... .......... .......... .......... .......... 18%  106M 1s
    ##  15150K .......... .......... .......... .......... .......... 19%  140M 1s
    ##  15200K .......... .......... .......... .......... .......... 19%  119M 1s
    ##  15250K .......... .......... .......... .......... .......... 19%  151M 1s
    ##  15300K .......... .......... .......... .......... .......... 19%  122M 1s
    ##  15350K .......... .......... .......... .......... .......... 19% 16.7M 1s
    ##  15400K .......... .......... .......... .......... .......... 19%  101M 1s
    ##  15450K .......... .......... .......... .......... .......... 19%  107M 1s
    ##  15500K .......... .......... .......... .......... .......... 19%  127M 1s
    ##  15550K .......... .......... .......... .......... .......... 19%  172M 1s
    ##  15600K .......... .......... .......... .......... .......... 19%  149M 1s
    ##  15650K .......... .......... .......... .......... .......... 19%  122M 1s
    ##  15700K .......... .......... .......... .......... .......... 19%  162M 1s
    ##  15750K .......... .......... .......... .......... .......... 19%  140M 1s
    ##  15800K .......... .......... .......... .......... .......... 19%  146M 1s
    ##  15850K .......... .......... .......... .......... .......... 19%  174M 1s
    ##  15900K .......... .......... .......... .......... .......... 19%  140M 1s
    ##  15950K .......... .......... .......... .......... .......... 20%  161M 1s
    ##  16000K .......... .......... .......... .......... .......... 20%  161M 1s
    ##  16050K .......... .......... .......... .......... .......... 20%  131M 1s
    ##  16100K .......... .......... .......... .......... .......... 20%  146M 1s
    ##  16150K .......... .......... .......... .......... .......... 20%  168M 1s
    ##  16200K .......... .......... .......... .......... .......... 20%  146M 1s
    ##  16250K .......... .......... .......... .......... .......... 20%  168M 1s
    ##  16300K .......... .......... .......... .......... .......... 20%  147M 1s
    ##  16350K .......... .......... .......... .......... .......... 20%  101M 1s
    ##  16400K .......... .......... .......... .......... .......... 20%  150M 1s
    ##  16450K .......... .......... .......... .......... .......... 20%  159M 1s
    ##  16500K .......... .......... .......... .......... .......... 20%  142M 1s
    ##  16550K .......... .......... .......... .......... .......... 20%  174M 1s
    ##  16600K .......... .......... .......... .......... .......... 20%  150M 1s
    ##  16650K .......... .......... .......... .......... .......... 20%  162M 1s
    ##  16700K .......... .......... .......... .......... .......... 20%  151M 1s
    ##  16750K .......... .......... .......... .......... .......... 21%  151M 1s
    ##  16800K .......... .......... .......... .......... .......... 21%  163M 1s
    ##  16850K .......... .......... .......... .......... .......... 21%  163M 1s
    ##  16900K .......... .......... .......... .......... .......... 21%  134M 1s
    ##  16950K .......... .......... .......... .......... .......... 21%  178M 1s
    ##  17000K .......... .......... .......... .......... .......... 21%  146M 1s
    ##  17050K .......... .......... .......... .......... .......... 21% 44.8M 1s
    ##  17100K .......... .......... .......... .......... .......... 21%  138M 1s
    ##  17150K .......... .......... .......... .......... .......... 21%  125M 1s
    ##  17200K .......... .......... .......... .......... .......... 21%  114M 1s
    ##  17250K .......... .......... .......... .......... .......... 21%  132M 1s
    ##  17300K .......... .......... .......... .......... .......... 21%  121M 1s
    ##  17350K .......... .......... .......... .......... .......... 21%  133M 1s
    ##  17400K .......... .......... .......... .......... .......... 21%  154M 1s
    ##  17450K .......... .......... .......... .......... .......... 21%  151M 1s
    ##  17500K .......... .......... .......... .......... .......... 21%  151M 1s
    ##  17550K .......... .......... .......... .......... .......... 22%  152M 1s
    ##  17600K .......... .......... .......... .......... .......... 22%  130M 1s
    ##  17650K .......... .......... .......... .......... .......... 22% 10.9M 1s
    ##  17700K .......... .......... .......... .......... .......... 22%  102M 1s
    ##  17750K .......... .......... .......... .......... .......... 22%  128M 1s
    ##  17800K .......... .......... .......... .......... .......... 22%  159M 1s
    ##  17850K .......... .......... .......... .......... .......... 22%  184M 1s
    ##  17900K .......... .......... .......... .......... .......... 22%  162M 1s
    ##  17950K .......... .......... .......... .......... .......... 22%  156M 1s
    ##  18000K .......... .......... .......... .......... .......... 22%  177M 1s
    ##  18050K .......... .......... .......... .......... .......... 22%  172M 1s
    ##  18100K .......... .......... .......... .......... .......... 22%  156M 1s
    ##  18150K .......... .......... .......... .......... .......... 22%  187M 1s
    ##  18200K .......... .......... .......... .......... .......... 22%  179M 1s
    ##  18250K .......... .......... .......... .......... .......... 22% 26.1M 1s
    ##  18300K .......... .......... .......... .......... .......... 22% 48.9M 1s
    ##  18350K .......... .......... .......... .......... .......... 23% 42.2M 1s
    ##  18400K .......... .......... .......... .......... .......... 23%  149M 1s
    ##  18450K .......... .......... .......... .......... .......... 23% 74.8M 1s
    ##  18500K .......... .......... .......... .......... .......... 23%  158M 1s
    ##  18550K .......... .......... .......... .......... .......... 23%  103M 1s
    ##  18600K .......... .......... .......... .......... .......... 23%  127M 1s
    ##  18650K .......... .......... .......... .......... .......... 23%  173M 1s
    ##  18700K .......... .......... .......... .......... .......... 23%  118M 1s
    ##  18750K .......... .......... .......... .......... .......... 23%  163M 1s
    ##  18800K .......... .......... .......... .......... .......... 23%  151M 1s
    ##  18850K .......... .......... .......... .......... .......... 23%  156M 1s
    ##  18900K .......... .......... .......... .......... .......... 23%  135M 1s
    ##  18950K .......... .......... .......... .......... .......... 23% 82.2M 1s
    ##  19000K .......... .......... .......... .......... .......... 23% 77.6M 1s
    ##  19050K .......... .......... .......... .......... .......... 23% 35.2M 1s
    ##  19100K .......... .......... .......... .......... .......... 23%  140M 1s
    ##  19150K .......... .......... .......... .......... .......... 24%  111M 1s
    ##  19200K .......... .......... .......... .......... .......... 24%  116M 1s
    ##  19250K .......... .......... .......... .......... .......... 24% 93.8M 1s
    ##  19300K .......... .......... .......... .......... .......... 24%  104M 1s
    ##  19350K .......... .......... .......... .......... .......... 24%  131M 1s
    ##  19400K .......... .......... .......... .......... .......... 24%  148M 1s
    ##  19450K .......... .......... .......... .......... .......... 24% 75.4M 1s
    ##  19500K .......... .......... .......... .......... .......... 24%  105M 1s
    ##  19550K .......... .......... .......... .......... .......... 24%  167M 1s
    ##  19600K .......... .......... .......... .......... .......... 24%  139M 1s
    ##  19650K .......... .......... .......... .......... .......... 24% 65.3M 1s
    ##  19700K .......... .......... .......... .......... .......... 24%  120M 1s
    ##  19750K .......... .......... .......... .......... .......... 24%  146M 1s
    ##  19800K .......... .......... .......... .......... .......... 24% 48.8M 1s
    ##  19850K .......... .......... .......... .......... .......... 24%  151M 1s
    ##  19900K .......... .......... .......... .......... .......... 24%  116M 1s
    ##  19950K .......... .......... .......... .......... .......... 25% 47.3M 1s
    ##  20000K .......... .......... .......... .......... .......... 25% 77.9M 1s
    ##  20050K .......... .......... .......... .......... .......... 25%  124M 1s
    ##  20100K .......... .......... .......... .......... .......... 25%  135M 1s
    ##  20150K .......... .......... .......... .......... .......... 25% 82.0M 1s
    ##  20200K .......... .......... .......... .......... .......... 25%  125M 1s
    ##  20250K .......... .......... .......... .......... .......... 25%  130M 1s
    ##  20300K .......... .......... .......... .......... .......... 25% 93.6M 1s
    ##  20350K .......... .......... .......... .......... .......... 25% 83.4M 1s
    ##  20400K .......... .......... .......... .......... .......... 25%  125M 1s
    ##  20450K .......... .......... .......... .......... .......... 25%  155M 1s
    ##  20500K .......... .......... .......... .......... .......... 25% 50.8M 1s
    ##  20550K .......... .......... .......... .......... .......... 25%  122M 1s
    ##  20600K .......... .......... .......... .......... .......... 25%  130M 1s
    ##  20650K .......... .......... .......... .......... .......... 25%  174M 1s
    ##  20700K .......... .......... .......... .......... .......... 25% 47.9M 1s
    ##  20750K .......... .......... .......... .......... .......... 26%  139M 1s
    ##  20800K .......... .......... .......... .......... .......... 26%  119M 1s
    ##  20850K .......... .......... .......... .......... .......... 26%  154M 1s
    ##  20900K .......... .......... .......... .......... .......... 26%  110M 1s
    ##  20950K .......... .......... .......... .......... .......... 26%  153M 1s
    ##  21000K .......... .......... .......... .......... .......... 26% 60.9M 1s
    ##  21050K .......... .......... .......... .......... .......... 26% 97.2M 1s
    ##  21100K .......... .......... .......... .......... .......... 26%  125M 1s
    ##  21150K .......... .......... .......... .......... .......... 26%  149M 1s
    ##  21200K .......... .......... .......... .......... .......... 26% 44.5M 1s
    ##  21250K .......... .......... .......... .......... .......... 26%  139M 1s
    ##  21300K .......... .......... .......... .......... .......... 26%  148M 1s
    ##  21350K .......... .......... .......... .......... .......... 26%  168M 1s
    ##  21400K .......... .......... .......... .......... .......... 26%  115M 1s
    ##  21450K .......... .......... .......... .......... .......... 26%  145M 1s
    ##  21500K .......... .......... .......... .......... .......... 26%  133M 1s
    ##  21550K .......... .......... .......... .......... .......... 27% 81.5M 1s
    ##  21600K .......... .......... .......... .......... .......... 27%  121M 1s
    ##  21650K .......... .......... .......... .......... .......... 27%  171M 1s
    ##  21700K .......... .......... .......... .......... .......... 27% 50.5M 1s
    ##  21750K .......... .......... .......... .......... .......... 27%  124M 1s
    ##  21800K .......... .......... .......... .......... .......... 27%  121M 1s
    ##  21850K .......... .......... .......... .......... .......... 27%  157M 1s
    ##  21900K .......... .......... .......... .......... .......... 27% 98.2M 1s
    ##  21950K .......... .......... .......... .......... .......... 27%  151M 1s
    ##  22000K .......... .......... .......... .......... .......... 27%  143M 1s
    ##  22050K .......... .......... .......... .......... .......... 27%  113M 1s
    ##  22100K .......... .......... .......... .......... .......... 27%  140M 1s
    ##  22150K .......... .......... .......... .......... .......... 27%  115M 1s
    ##  22200K .......... .......... .......... .......... .......... 27% 61.1M 1s
    ##  22250K .......... .......... .......... .......... .......... 27% 97.0M 1s
    ##  22300K .......... .......... .......... .......... .......... 27%  128M 1s
    ##  22350K .......... .......... .......... .......... .......... 28%  160M 1s
    ##  22400K .......... .......... .......... .......... .......... 28% 63.5M 1s
    ##  22450K .......... .......... .......... .......... .......... 28%  144M 1s
    ##  22500K .......... .......... .......... .......... .......... 28%  115M 1s
    ##  22550K .......... .......... .......... .......... .......... 28%  191M 1s
    ##  22600K .......... .......... .......... .......... .......... 28% 95.6M 1s
    ##  22650K .......... .......... .......... .......... .......... 28%  121M 1s
    ##  22700K .......... .......... .......... .......... .......... 28% 66.8M 1s
    ##  22750K .......... .......... .......... .......... .......... 28%  159M 1s
    ##  22800K .......... .......... .......... .......... .......... 28%  127M 1s
    ##  22850K .......... .......... .......... .......... .......... 28%  194M 1s
    ##  22900K .......... .......... .......... .......... .......... 28% 70.1M 1s
    ##  22950K .......... .......... .......... .......... .......... 28%  118M 1s
    ##  23000K .......... .......... .......... .......... .......... 28%  116M 1s
    ##  23050K .......... .......... .......... .......... .......... 28%  161M 1s
    ##  23100K .......... .......... .......... .......... .......... 28% 97.6M 1s
    ##  23150K .......... .......... .......... .......... .......... 29%  147M 1s
    ##  23200K .......... .......... .......... .......... .......... 29% 78.7M 1s
    ##  23250K .......... .......... .......... .......... .......... 29%  126M 1s
    ##  23300K .......... .......... .......... .......... .......... 29%  121M 1s
    ##  23350K .......... .......... .......... .......... .......... 29%  152M 1s
    ##  23400K .......... .......... .......... .......... .......... 29%  124M 1s
    ##  23450K .......... .......... .......... .......... .......... 29%  168M 1s
    ##  23500K .......... .......... .......... .......... .......... 29%  119M 1s
    ##  23550K .......... .......... .......... .......... .......... 29% 69.3M 1s
    ##  23600K .......... .......... .......... .......... .......... 29%  114M 1s
    ##  23650K .......... .......... .......... .......... .......... 29%  158M 1s
    ##  23700K .......... .......... .......... .......... .......... 29%  131M 1s
    ##  23750K .......... .......... .......... .......... .......... 29% 69.3M 1s
    ##  23800K .......... .......... .......... .......... .......... 29%  161M 1s
    ##  23850K .......... .......... .......... .......... .......... 29%  156M 1s
    ##  23900K .......... .......... .......... .......... .......... 29% 88.4M 1s
    ##  23950K .......... .......... .......... .......... .......... 30%  112M 1s
    ##  24000K .......... .......... .......... .......... .......... 30%  176M 1s
    ##  24050K .......... .......... .......... .......... .......... 30%  150M 1s
    ##  24100K .......... .......... .......... .......... .......... 30%  114M 1s
    ##  24150K .......... .......... .......... .......... .......... 30%  147M 1s
    ##  24200K .......... .......... .......... .......... .......... 30%  142M 1s
    ##  24250K .......... .......... .......... .......... .......... 30% 80.7M 1s
    ##  24300K .......... .......... .......... .......... .......... 30%  159M 1s
    ##  24350K .......... .......... .......... .......... .......... 30%  136M 1s
    ##  24400K .......... .......... .......... .......... .......... 30%  142M 1s
    ##  24450K .......... .......... .......... .......... .......... 30%  104M 1s
    ##  24500K .......... .......... .......... .......... .......... 30%  126M 1s
    ##  24550K .......... .......... .......... .......... .......... 30% 80.9M 1s
    ##  24600K .......... .......... .......... .......... .......... 30%  134M 1s
    ##  24650K .......... .......... .......... .......... .......... 30%  171M 1s
    ##  24700K .......... .......... .......... .......... .......... 30%  141M 1s
    ##  24750K .......... .......... .......... .......... .......... 31%  114M 1s
    ##  24800K .......... .......... .......... .......... .......... 31%  122M 1s
    ##  24850K .......... .......... .......... .......... .......... 31% 91.3M 1s
    ##  24900K .......... .......... .......... .......... .......... 31%  161M 1s
    ##  24950K .......... .......... .......... .......... .......... 31%  146M 1s
    ##  25000K .......... .......... .......... .......... .......... 31%  125M 1s
    ##  25050K .......... .......... .......... .......... .......... 31% 96.8M 1s
    ##  25100K .......... .......... .......... .......... .......... 31% 93.1M 1s
    ##  25150K .......... .......... .......... .......... .......... 31%  156M 1s
    ##  25200K .......... .......... .......... .......... .......... 31%  137M 1s
    ##  25250K .......... .......... .......... .......... .......... 31%  138M 1s
    ##  25300K .......... .......... .......... .......... .......... 31%  129M 1s
    ##  25350K .......... .......... .......... .......... .......... 31%  176M 1s
    ##  25400K .......... .......... .......... .......... .......... 31% 93.7M 1s
    ##  25450K .......... .......... .......... .......... .......... 31%  144M 1s
    ##  25500K .......... .......... .......... .......... .......... 31%  143M 1s
    ##  25550K .......... .......... .......... .......... .......... 32%  168M 1s
    ##  25600K .......... .......... .......... .......... .......... 32% 95.7M 1s
    ##  25650K .......... .......... .......... .......... .......... 32%  149M 1s
    ##  25700K .......... .......... .......... .......... .......... 32%  129M 1s
    ##  25750K .......... .......... .......... .......... .......... 32%  117M 1s
    ##  25800K .......... .......... .......... .......... .......... 32%  124M 1s
    ##  25850K .......... .......... .......... .......... .......... 32%  165M 1s
    ##  25900K .......... .......... .......... .......... .......... 32% 95.4M 1s
    ##  25950K .......... .......... .......... .......... .......... 32%  145M 1s
    ##  26000K .......... .......... .......... .......... .......... 32%  141M 1s
    ##  26050K .......... .......... .......... .......... .......... 32%  118M 1s
    ##  26100K .......... .......... .......... .......... .......... 32%  159M 1s
    ##  26150K .......... .......... .......... .......... .......... 32%  102M 1s
    ##  26200K .......... .......... .......... .......... .......... 32%  123M 1s
    ##  26250K .......... .......... .......... .......... .......... 32%  145M 1s
    ##  26300K .......... .......... .......... .......... .......... 32%  144M 1s
    ##  26350K .......... .......... .......... .......... .......... 33%  125M 1s
    ##  26400K .......... .......... .......... .......... .......... 33% 93.6M 1s
    ##  26450K .......... .......... .......... .......... .......... 33%  152M 1s
    ##  26500K .......... .......... .......... .......... .......... 33%  152M 1s
    ##  26550K .......... .......... .......... .......... .......... 33%  126M 1s
    ##  26600K .......... .......... .......... .......... .......... 33%  146M 1s
    ##  26650K .......... .......... .......... .......... .......... 33%  113M 1s
    ##  26700K .......... .......... .......... .......... .......... 33%  133M 1s
    ##  26750K .......... .......... .......... .......... .......... 33%  174M 1s
    ##  26800K .......... .......... .......... .......... .......... 33%  132M 1s
    ##  26850K .......... .......... .......... .......... .......... 33%  148M 1s
    ##  26900K .......... .......... .......... .......... .......... 33%  140M 1s
    ##  26950K .......... .......... .......... .......... .......... 33%  146M 1s
    ##  27000K .......... .......... .......... .......... .......... 33%  133M 1s
    ##  27050K .......... .......... .......... .......... .......... 33%  159M 1s
    ##  27100K .......... .......... .......... .......... .......... 33% 91.6M 1s
    ##  27150K .......... .......... .......... .......... .......... 34%  123M 1s
    ##  27200K .......... .......... .......... .......... .......... 34%  122M 1s
    ##  27250K .......... .......... .......... .......... .......... 34%  189M 1s
    ##  27300K .......... .......... .......... .......... .......... 34%  117M 1s
    ##  27350K .......... .......... .......... .......... .......... 34%  154M 1s
    ##  27400K .......... .......... .......... .......... .......... 34%  135M 1s
    ##  27450K .......... .......... .......... .......... .......... 34%  111M 1s
    ##  27500K .......... .......... .......... .......... .......... 34%  141M 1s
    ##  27550K .......... .......... .......... .......... .......... 34%  145M 1s
    ##  27600K .......... .......... .......... .......... .......... 34%  130M 1s
    ##  27650K .......... .......... .......... .......... .......... 34%  146M 1s
    ##  27700K .......... .......... .......... .......... .......... 34%  145M 1s
    ##  27750K .......... .......... .......... .......... .......... 34%  111M 1s
    ##  27800K .......... .......... .......... .......... .......... 34%  126M 1s
    ##  27850K .......... .......... .......... .......... .......... 34%  168M 1s
    ##  27900K .......... .......... .......... .......... .......... 34%  113M 1s
    ##  27950K .......... .......... .......... .......... .......... 35%  185M 1s
    ##  28000K .......... .......... .......... .......... .......... 35%  142M 1s
    ##  28050K .......... .......... .......... .......... .......... 35%  126M 1s
    ##  28100K .......... .......... .......... .......... .......... 35%  132M 1s
    ##  28150K .......... .......... .......... .......... .......... 35%  151M 1s
    ##  28200K .......... .......... .......... .......... .......... 35%  140M 1s
    ##  28250K .......... .......... .......... .......... .......... 35%  121M 1s
    ##  28300K .......... .......... .......... .......... .......... 35%  113M 1s
    ##  28350K .......... .......... .......... .......... .......... 35%  142M 1s
    ##  28400K .......... .......... .......... .......... .......... 35%  126M 1s
    ##  28450K .......... .......... .......... .......... .......... 35%  139M 1s
    ##  28500K .......... .......... .......... .......... .......... 35%  131M 1s
    ##  28550K .......... .......... .......... .......... .......... 35%  176M 1s
    ##  28600K .......... .......... .......... .......... .......... 35%  116M 1s
    ##  28650K .......... .......... .......... .......... .......... 35%  163M 1s
    ##  28700K .......... .......... .......... .......... .......... 35%  121M 1s
    ##  28750K .......... .......... .......... .......... .......... 36%  154M 1s
    ##  28800K .......... .......... .......... .......... .......... 36%  135M 1s
    ##  28850K .......... .......... .......... .......... .......... 36%  108M 1s
    ##  28900K .......... .......... .......... .......... .......... 36%  122M 1s
    ##  28950K .......... .......... .......... .......... .......... 36%  184M 1s
    ##  29000K .......... .......... .......... .......... .......... 36%  112M 1s
    ##  29050K .......... .......... .......... .......... .......... 36%  169M 1s
    ##  29100K .......... .......... .......... .......... .......... 36%  143M 1s
    ##  29150K .......... .......... .......... .......... .......... 36%  145M 1s
    ##  29200K .......... .......... .......... .......... .......... 36%  127M 1s
    ##  29250K .......... .......... .......... .......... .......... 36%  166M 1s
    ##  29300K .......... .......... .......... .......... .......... 36%  139M 1s
    ##  29350K .......... .......... .......... .......... .......... 36%  196M 1s
    ##  29400K .......... .......... .......... .......... .......... 36%  121M 1s
    ##  29450K .......... .......... .......... .......... .......... 36%  103M 1s
    ##  29500K .......... .......... .......... .......... .......... 36%  117M 1s
    ##  29550K .......... .......... .......... .......... .......... 37%  165M 1s
    ##  29600K .......... .......... .......... .......... .......... 37%  137M 1s
    ##  29650K .......... .......... .......... .......... .......... 37%  202M 1s
    ##  29700K .......... .......... .......... .......... .......... 37%  130M 1s
    ##  29750K .......... .......... .......... .......... .......... 37%  133M 1s
    ##  29800K .......... .......... .......... .......... .......... 37%  116M 1s
    ##  29850K .......... .......... .......... .......... .......... 37%  166M 1s
    ##  29900K .......... .......... .......... .......... .......... 37%  140M 1s
    ##  29950K .......... .......... .......... .......... .......... 37%  180M 1s
    ##  30000K .......... .......... .......... .......... .......... 37%  131M 1s
    ##  30050K .......... .......... .......... .......... .......... 37%  121M 1s
    ##  30100K .......... .......... .......... .......... .......... 37%  109M 1s
    ##  30150K .......... .......... .......... .......... .......... 37%  133M 1s
    ##  30200K .......... .......... .......... .......... .......... 37%  122M 1s
    ##  30250K .......... .......... .......... .......... .......... 37%  187M 1s
    ##  30300K .......... .......... .......... .......... .......... 37%  129M 1s
    ##  30350K .......... .......... .......... .......... .......... 38%  144M 1s
    ##  30400K .......... .......... .......... .......... .......... 38%  115M 1s
    ##  30450K .......... .......... .......... .......... .......... 38%  145M 1s
    ##  30500K .......... .......... .......... .......... .......... 38%  131M 1s
    ##  30550K .......... .......... .......... .......... .......... 38%  176M 1s
    ##  30600K .......... .......... .......... .......... .......... 38%  118M 1s
    ##  30650K .......... .......... .......... .......... .......... 38%  132M 1s
    ##  30700K .......... .......... .......... .......... .......... 38%  101M 1s
    ##  30750K .......... .......... .......... .......... .......... 38%  146M 1s
    ##  30800K .......... .......... .......... .......... .......... 38%  131M 1s
    ##  30850K .......... .......... .......... .......... .......... 38%  168M 1s
    ##  30900K .......... .......... .......... .......... .......... 38%  118M 1s
    ##  30950K .......... .......... .......... .......... .......... 38%  119M 1s
    ##  31000K .......... .......... .......... .......... .......... 38%  116M 1s
    ##  31050K .......... .......... .......... .......... .......... 38%  141M 1s
    ##  31100K .......... .......... .......... .......... .......... 38%  128M 1s
    ##  31150K .......... .......... .......... .......... .......... 39%  168M 1s
    ##  31200K .......... .......... .......... .......... .......... 39%  128M 1s
    ##  31250K .......... .......... .......... .......... .......... 39%  121M 1s
    ##  31300K .......... .......... .......... .......... .......... 39%  142M 1s
    ##  31350K .......... .......... .......... .......... .......... 39%  140M 1s
    ##  31400K .......... .......... .......... .......... .......... 39%  133M 1s
    ##  31450K .......... .......... .......... .......... .......... 39%  194M 1s
    ##  31500K .......... .......... .......... .......... .......... 39%  153M 1s
    ##  31550K .......... .......... .......... .......... .......... 39%  153M 1s
    ##  31600K .......... .......... .......... .......... .......... 39%  149M 1s
    ##  31650K .......... .......... .......... .......... .......... 39%  134M 1s
    ##  31700K .......... .......... .......... .......... .......... 39% 45.9M 1s
    ##  31750K .......... .......... .......... .......... .......... 39%  177M 1s
    ##  31800K .......... .......... .......... .......... .......... 39%  140M 1s
    ##  31850K .......... .......... .......... .......... .......... 39%  139M 0s
    ##  31900K .......... .......... .......... .......... .......... 39%  150M 0s
    ##  31950K .......... .......... .......... .......... .......... 40%  176M 0s
    ##  32000K .......... .......... .......... .......... .......... 40%  159M 0s
    ##  32050K .......... .......... .......... .......... .......... 40%  166M 0s
    ##  32100K .......... .......... .......... .......... .......... 40%  150M 0s
    ##  32150K .......... .......... .......... .......... .......... 40%  147M 0s
    ##  32200K .......... .......... .......... .......... .......... 40%  123M 0s
    ##  32250K .......... .......... .......... .......... .......... 40%  184M 0s
    ##  32300K .......... .......... .......... .......... .......... 40% 47.2M 0s
    ##  32350K .......... .......... .......... .......... .......... 40%  190M 0s
    ##  32400K .......... .......... .......... .......... .......... 40% 57.8M 0s
    ##  32450K .......... .......... .......... .......... .......... 40% 89.5M 0s
    ##  32500K .......... .......... .......... .......... .......... 40% 42.1M 0s
    ##  32550K .......... .......... .......... .......... .......... 40% 23.1M 0s
    ##  32600K .......... .......... .......... .......... .......... 40%  155M 0s
    ##  32650K .......... .......... .......... .......... .......... 40%  168M 0s
    ##  32700K .......... .......... .......... .......... .......... 40% 37.2M 0s
    ##  32750K .......... .......... .......... .......... .......... 41%  146M 0s
    ##  32800K .......... .......... .......... .......... .......... 41%  140M 0s
    ##  32850K .......... .......... .......... .......... .......... 41%  160M 0s
    ##  32900K .......... .......... .......... .......... .......... 41%  124M 0s
    ##  32950K .......... .......... .......... .......... .......... 41%  170M 0s
    ##  33000K .......... .......... .......... .......... .......... 41%  160M 0s
    ##  33050K .......... .......... .......... .......... .......... 41%  125M 0s
    ##  33100K .......... .......... .......... .......... .......... 41%  135M 0s
    ##  33150K .......... .......... .......... .......... .......... 41%  177M 0s
    ##  33200K .......... .......... .......... .......... .......... 41%  138M 0s
    ##  33250K .......... .......... .......... .......... .......... 41%  157M 0s
    ##  33300K .......... .......... .......... .......... .......... 41%  143M 0s
    ##  33350K .......... .......... .......... .......... .......... 41%  136M 0s
    ##  33400K .......... .......... .......... .......... .......... 41%  145M 0s
    ##  33450K .......... .......... .......... .......... .......... 41%  143M 0s
    ##  33500K .......... .......... .......... .......... .......... 41%  159M 0s
    ##  33550K .......... .......... .......... .......... .......... 42%  161M 0s
    ##  33600K .......... .......... .......... .......... .......... 42%  156M 0s
    ##  33650K .......... .......... .......... .......... .......... 42%  127M 0s
    ##  33700K .......... .......... .......... .......... .......... 42%  122M 0s
    ##  33750K .......... .......... .......... .......... .......... 42%  175M 0s
    ##  33800K .......... .......... .......... .......... .......... 42%  169M 0s
    ##  33850K .......... .......... .......... .......... .......... 42%  133M 0s
    ##  33900K .......... .......... .......... .......... .......... 42%  166M 0s
    ##  33950K .......... .......... .......... .......... .......... 42%  100M 0s
    ##  34000K .......... .......... .......... .......... .......... 42%  106M 0s
    ##  34050K .......... .......... .......... .......... .......... 42% 84.9M 0s
    ##  34100K .......... .......... .......... .......... .......... 42%  133M 0s
    ##  34150K .......... .......... .......... .......... .......... 42%  144M 0s
    ##  34200K .......... .......... .......... .......... .......... 42% 40.9M 0s
    ##  34250K .......... .......... .......... .......... .......... 42% 96.2M 0s
    ##  34300K .......... .......... .......... .......... .......... 42%  153M 0s
    ##  34350K .......... .......... .......... .......... .......... 43%  177M 0s
    ##  34400K .......... .......... .......... .......... .......... 43% 25.8M 0s
    ##  34450K .......... .......... .......... .......... .......... 43% 34.6M 0s
    ##  34500K .......... .......... .......... .......... .......... 43% 90.0M 0s
    ##  34550K .......... .......... .......... .......... .......... 43% 40.5M 0s
    ##  34600K .......... .......... .......... .......... .......... 43% 60.5M 0s
    ##  34650K .......... .......... .......... .......... .......... 43%  157M 0s
    ##  34700K .......... .......... .......... .......... .......... 43%  161M 0s
    ##  34750K .......... .......... .......... .......... .......... 43%  126M 0s
    ##  34800K .......... .......... .......... .......... .......... 43%  117M 0s
    ##  34850K .......... .......... .......... .......... .......... 43%  118M 0s
    ##  34900K .......... .......... .......... .......... .......... 43%  118M 0s
    ##  34950K .......... .......... .......... .......... .......... 43%  147M 0s
    ##  35000K .......... .......... .......... .......... .......... 43%  160M 0s
    ##  35050K .......... .......... .......... .......... .......... 43%  115M 0s
    ##  35100K .......... .......... .......... .......... .......... 43%  145M 0s
    ##  35150K .......... .......... .......... .......... .......... 44%  146M 0s
    ##  35200K .......... .......... .......... .......... .......... 44%  155M 0s
    ##  35250K .......... .......... .......... .......... .......... 44%  162M 0s
    ##  35300K .......... .......... .......... .......... .......... 44%  121M 0s
    ##  35350K .......... .......... .......... .......... .......... 44%  146M 0s
    ##  35400K .......... .......... .......... .......... .......... 44%  151M 0s
    ##  35450K .......... .......... .......... .......... .......... 44%  178M 0s
    ##  35500K .......... .......... .......... .......... .......... 44% 78.9M 0s
    ##  35550K .......... .......... .......... .......... .......... 44% 97.3M 0s
    ##  35600K .......... .......... .......... .......... .......... 44%  137M 0s
    ##  35650K .......... .......... .......... .......... .......... 44% 78.3M 0s
    ##  35700K .......... .......... .......... .......... .......... 44% 80.8M 0s
    ##  35750K .......... .......... .......... .......... .......... 44%  150M 0s
    ##  35800K .......... .......... .......... .......... .......... 44% 88.5M 0s
    ##  35850K .......... .......... .......... .......... .......... 44%  113M 0s
    ##  35900K .......... .......... .......... .......... .......... 44% 50.7M 0s
    ##  35950K .......... .......... .......... .......... .......... 45%  153M 0s
    ##  36000K .......... .......... .......... .......... .......... 45%  106M 0s
    ##  36050K .......... .......... .......... .......... .......... 45% 66.8M 0s
    ##  36100K .......... .......... .......... .......... .......... 45% 53.1M 0s
    ##  36150K .......... .......... .......... .......... .......... 45%  100M 0s
    ##  36200K .......... .......... .......... .......... .......... 45%  170M 0s
    ##  36250K .......... .......... .......... .......... .......... 45%  110M 0s
    ##  36300K .......... .......... .......... .......... .......... 45%  106M 0s
    ##  36350K .......... .......... .......... .......... .......... 45%  131M 0s
    ##  36400K .......... .......... .......... .......... .......... 45%  136M 0s
    ##  36450K .......... .......... .......... .......... .......... 45%  142M 0s
    ##  36500K .......... .......... .......... .......... .......... 45% 74.9M 0s
    ##  36550K .......... .......... .......... .......... .......... 45%  101M 0s
    ##  36600K .......... .......... .......... .......... .......... 45% 94.7M 0s
    ##  36650K .......... .......... .......... .......... .......... 45%  122M 0s
    ##  36700K .......... .......... .......... .......... .......... 45%  116M 0s
    ##  36750K .......... .......... .......... .......... .......... 46%  182M 0s
    ##  36800K .......... .......... .......... .......... .......... 46% 14.7M 0s
    ##  36850K .......... .......... .......... .......... .......... 46%  151M 0s
    ##  36900K .......... .......... .......... .......... .......... 46%  112M 0s
    ##  36950K .......... .......... .......... .......... .......... 46%  167M 0s
    ##  37000K .......... .......... .......... .......... .......... 46%  113M 0s
    ##  37050K .......... .......... .......... .......... .......... 46%  181M 0s
    ##  37100K .......... .......... .......... .......... .......... 46%  124M 0s
    ##  37150K .......... .......... .......... .......... .......... 46%  166M 0s
    ##  37200K .......... .......... .......... .......... .......... 46%  111M 0s
    ##  37250K .......... .......... .......... .......... .......... 46%  127M 0s
    ##  37300K .......... .......... .......... .......... .......... 46%  155M 0s
    ##  37350K .......... .......... .......... .......... .......... 46%  178M 0s
    ##  37400K .......... .......... .......... .......... .......... 46% 82.6M 0s
    ##  37450K .......... .......... .......... .......... .......... 46%  188M 0s
    ##  37500K .......... .......... .......... .......... .......... 46%  113M 0s
    ##  37550K .......... .......... .......... .......... .......... 47% 44.8M 0s
    ##  37600K .......... .......... .......... .......... .......... 47%  117M 0s
    ##  37650K .......... .......... .......... .......... .......... 47% 85.5M 0s
    ##  37700K .......... .......... .......... .......... .......... 47% 50.8M 0s
    ##  37750K .......... .......... .......... .......... .......... 47% 36.7M 0s
    ##  37800K .......... .......... .......... .......... .......... 47% 91.5M 0s
    ##  37850K .......... .......... .......... .......... .......... 47%  150M 0s
    ##  37900K .......... .......... .......... .......... .......... 47%  149M 0s
    ##  37950K .......... .......... .......... .......... .......... 47%  172M 0s
    ##  38000K .......... .......... .......... .......... .......... 47%  124M 0s
    ##  38050K .......... .......... .......... .......... .......... 47%  174M 0s
    ##  38100K .......... .......... .......... .......... .......... 47%  126M 0s
    ##  38150K .......... .......... .......... .......... .......... 47%  157M 0s
    ##  38200K .......... .......... .......... .......... .......... 47%  159M 0s
    ##  38250K .......... .......... .......... .......... .......... 47%  138M 0s
    ##  38300K .......... .......... .......... .......... .......... 47%  152M 0s
    ##  38350K .......... .......... .......... .......... .......... 48%  191M 0s
    ##  38400K .......... .......... .......... .......... .......... 48% 55.9M 0s
    ##  38450K .......... .......... .......... .......... .......... 48% 59.1M 0s
    ##  38500K .......... .......... .......... .......... .......... 48% 62.7M 0s
    ##  38550K .......... .......... .......... .......... .......... 48% 99.8M 0s
    ##  38600K .......... .......... .......... .......... .......... 48% 37.4M 0s
    ##  38650K .......... .......... .......... .......... .......... 48%  144M 0s
    ##  38700K .......... .......... .......... .......... .......... 48%  138M 0s
    ##  38750K .......... .......... .......... .......... .......... 48%  166M 0s
    ##  38800K .......... .......... .......... .......... .......... 48%  122M 0s
    ##  38850K .......... .......... .......... .......... .......... 48%  144M 0s
    ##  38900K .......... .......... .......... .......... .......... 48%  101M 0s
    ##  38950K .......... .......... .......... .......... .......... 48%  140M 0s
    ##  39000K .......... .......... .......... .......... .......... 48%  154M 0s
    ##  39050K .......... .......... .......... .......... .......... 48%  140M 0s
    ##  39100K .......... .......... .......... .......... .......... 48%  150M 0s
    ##  39150K .......... .......... .......... .......... .......... 49%  155M 0s
    ##  39200K .......... .......... .......... .......... .......... 49%  157M 0s
    ##  39250K .......... .......... .......... .......... .......... 49%  160M 0s
    ##  39300K .......... .......... .......... .......... .......... 49%  111M 0s
    ##  39350K .......... .......... .......... .......... .......... 49%  169M 0s
    ##  39400K .......... .......... .......... .......... .......... 49%  150M 0s
    ##  39450K .......... .......... .......... .......... .......... 49%  121M 0s
    ##  39500K .......... .......... .......... .......... .......... 49% 89.0M 0s
    ##  39550K .......... .......... .......... .......... .......... 49%  194M 0s
    ##  39600K .......... .......... .......... .......... .......... 49%  161M 0s
    ##  39650K .......... .......... .......... .......... .......... 49%  144M 0s
    ##  39700K .......... .......... .......... .......... .......... 49%  157M 0s
    ##  39750K .......... .......... .......... .......... .......... 49%  187M 0s
    ##  39800K .......... .......... .......... .......... .......... 49%  158M 0s
    ##  39850K .......... .......... .......... .......... .......... 49%  196M 0s
    ##  39900K .......... .......... .......... .......... .......... 49% 49.8M 0s
    ##  39950K .......... .......... .......... .......... .......... 50%  126M 0s
    ##  40000K .......... .......... .......... .......... .......... 50% 83.1M 0s
    ##  40050K .......... .......... .......... .......... .......... 50% 75.4M 0s
    ##  40100K .......... .......... .......... .......... .......... 50% 29.1M 0s
    ##  40150K .......... .......... .......... .......... .......... 50% 94.0M 0s
    ##  40200K .......... .......... .......... .......... .......... 50% 83.2M 0s
    ##  40250K .......... .......... .......... .......... .......... 50% 98.8M 0s
    ##  40300K .......... .......... .......... .......... .......... 50%  131M 0s
    ##  40350K .......... .......... .......... .......... .......... 50%  136M 0s
    ##  40400K .......... .......... .......... .......... .......... 50%  146M 0s
    ##  40450K .......... .......... .......... .......... .......... 50%  193M 0s
    ##  40500K .......... .......... .......... .......... .......... 50% 37.8M 0s
    ##  40550K .......... .......... .......... .......... .......... 50%  142M 0s
    ##  40600K .......... .......... .......... .......... .......... 50% 85.2M 0s
    ##  40650K .......... .......... .......... .......... .......... 50%  166M 0s
    ##  40700K .......... .......... .......... .......... .......... 50% 81.5M 0s
    ##  40750K .......... .......... .......... .......... .......... 51%  186M 0s
    ##  40800K .......... .......... .......... .......... .......... 51%  153M 0s
    ##  40850K .......... .......... .......... .......... .......... 51%  162M 0s
    ##  40900K .......... .......... .......... .......... .......... 51% 88.3M 0s
    ##  40950K .......... .......... .......... .......... .......... 51%  176M 0s
    ##  41000K .......... .......... .......... .......... .......... 51%  161M 0s
    ##  41050K .......... .......... .......... .......... .......... 51%  185M 0s
    ##  41100K .......... .......... .......... .......... .......... 51% 13.9M 0s
    ##  41150K .......... .......... .......... .......... .......... 51% 71.7M 0s
    ##  41200K .......... .......... .......... .......... .......... 51%  134M 0s
    ##  41250K .......... .......... .......... .......... .......... 51%  186M 0s
    ##  41300K .......... .......... .......... .......... .......... 51%  159M 0s
    ##  41350K .......... .......... .......... .......... .......... 51%  184M 0s
    ##  41400K .......... .......... .......... .......... .......... 51%  158M 0s
    ##  41450K .......... .......... .......... .......... .......... 51%  164M 0s
    ##  41500K .......... .......... .......... .......... .......... 51%  153M 0s
    ##  41550K .......... .......... .......... .......... .......... 52%  188M 0s
    ##  41600K .......... .......... .......... .......... .......... 52%  157M 0s
    ##  41650K .......... .......... .......... .......... .......... 52% 24.6M 0s
    ##  41700K .......... .......... .......... .......... .......... 52% 18.7M 0s
    ##  41750K .......... .......... .......... .......... .......... 52%  155M 0s
    ##  41800K .......... .......... .......... .......... .......... 52%  114M 0s
    ##  41850K .......... .......... .......... .......... .......... 52%  183M 0s
    ##  41900K .......... .......... .......... .......... .......... 52%  161M 0s
    ##  41950K .......... .......... .......... .......... .......... 52%  132M 0s
    ##  42000K .......... .......... .......... .......... .......... 52%  111M 0s
    ##  42050K .......... .......... .......... .......... .......... 52%  161M 0s
    ##  42100K .......... .......... .......... .......... .......... 52%  150M 0s
    ##  42150K .......... .......... .......... .......... .......... 52%  186M 0s
    ##  42200K .......... .......... .......... .......... .......... 52%  147M 0s
    ##  42250K .......... .......... .......... .......... .......... 52% 59.1M 0s
    ##  42300K .......... .......... .......... .......... .......... 52% 27.4M 0s
    ##  42350K .......... .......... .......... .......... .......... 53% 86.9M 0s
    ##  42400K .......... .......... .......... .......... .......... 53%  140M 0s
    ##  42450K .......... .......... .......... .......... .......... 53%  158M 0s
    ##  42500K .......... .......... .......... .......... .......... 53% 66.5M 0s
    ##  42550K .......... .......... .......... .......... .......... 53% 60.7M 0s
    ##  42600K .......... .......... .......... .......... .......... 53%  132M 0s
    ##  42650K .......... .......... .......... .......... .......... 53%  186M 0s
    ##  42700K .......... .......... .......... .......... .......... 53%  137M 0s
    ##  42750K .......... .......... .......... .......... .......... 53%  123M 0s
    ##  42800K .......... .......... .......... .......... .......... 53%  154M 0s
    ##  42850K .......... .......... .......... .......... .......... 53%  178M 0s
    ##  42900K .......... .......... .......... .......... .......... 53%  164M 0s
    ##  42950K .......... .......... .......... .......... .......... 53%  156M 0s
    ##  43000K .......... .......... .......... .......... .......... 53% 52.5M 0s
    ##  43050K .......... .......... .......... .......... .......... 53%  147M 0s
    ##  43100K .......... .......... .......... .......... .......... 53%  121M 0s
    ##  43150K .......... .......... .......... .......... .......... 54%  142M 0s
    ##  43200K .......... .......... .......... .......... .......... 54% 71.3M 0s
    ##  43250K .......... .......... .......... .......... .......... 54% 31.9M 0s
    ##  43300K .......... .......... .......... .......... .......... 54%  122M 0s
    ##  43350K .......... .......... .......... .......... .......... 54% 90.3M 0s
    ##  43400K .......... .......... .......... .......... .......... 54%  142M 0s
    ##  43450K .......... .......... .......... .......... .......... 54% 85.1M 0s
    ##  43500K .......... .......... .......... .......... .......... 54%  131M 0s
    ##  43550K .......... .......... .......... .......... .......... 54%  115M 0s
    ##  43600K .......... .......... .......... .......... .......... 54% 52.3M 0s
    ##  43650K .......... .......... .......... .......... .......... 54% 91.0M 0s
    ##  43700K .......... .......... .......... .......... .......... 54%  115M 0s
    ##  43750K .......... .......... .......... .......... .......... 54%  171M 0s
    ##  43800K .......... .......... .......... .......... .......... 54% 33.6M 0s
    ##  43850K .......... .......... .......... .......... .......... 54%  145M 0s
    ##  43900K .......... .......... .......... .......... .......... 54%  147M 0s
    ##  43950K .......... .......... .......... .......... .......... 55%  134M 0s
    ##  44000K .......... .......... .......... .......... .......... 55% 81.2M 0s
    ##  44050K .......... .......... .......... .......... .......... 55% 78.4M 0s
    ##  44100K .......... .......... .......... .......... .......... 55% 93.9M 0s
    ##  44150K .......... .......... .......... .......... .......... 55% 52.4M 0s
    ##  44200K .......... .......... .......... .......... .......... 55%  113M 0s
    ##  44250K .......... .......... .......... .......... .......... 55%  139M 0s
    ##  44300K .......... .......... .......... .......... .......... 55%  135M 0s
    ##  44350K .......... .......... .......... .......... .......... 55% 40.4M 0s
    ##  44400K .......... .......... .......... .......... .......... 55%  115M 0s
    ##  44450K .......... .......... .......... .......... .......... 55%  115M 0s
    ##  44500K .......... .......... .......... .......... .......... 55%  130M 0s
    ##  44550K .......... .......... .......... .......... .......... 55%  175M 0s
    ##  44600K .......... .......... .......... .......... .......... 55% 60.7M 0s
    ##  44650K .......... .......... .......... .......... .......... 55% 96.5M 0s
    ##  44700K .......... .......... .......... .......... .......... 55%  114M 0s
    ##  44750K .......... .......... .......... .......... .......... 56% 58.9M 0s
    ##  44800K .......... .......... .......... .......... .......... 56%  109M 0s
    ##  44850K .......... .......... .......... .......... .......... 56%  120M 0s
    ##  44900K .......... .......... .......... .......... .......... 56%  157M 0s
    ##  44950K .......... .......... .......... .......... .......... 56% 52.3M 0s
    ##  45000K .......... .......... .......... .......... .......... 56%  124M 0s
    ##  45050K .......... .......... .......... .......... .......... 56%  103M 0s
    ##  45100K .......... .......... .......... .......... .......... 56% 85.7M 0s
    ##  45150K .......... .......... .......... .......... .......... 56%  129M 0s
    ##  45200K .......... .......... .......... .......... .......... 56% 65.0M 0s
    ##  45250K .......... .......... .......... .......... .......... 56% 83.7M 0s
    ##  45300K .......... .......... .......... .......... .......... 56% 59.8M 0s
    ##  45350K .......... .......... .......... .......... .......... 56%  116M 0s
    ##  45400K .......... .......... .......... .......... .......... 56%  135M 0s
    ##  45450K .......... .......... .......... .......... .......... 56%  159M 0s
    ##  45500K .......... .......... .......... .......... .......... 56% 44.6M 0s
    ##  45550K .......... .......... .......... .......... .......... 57%  134M 0s
    ##  45600K .......... .......... .......... .......... .......... 57%  122M 0s
    ##  45650K .......... .......... .......... .......... .......... 57%  134M 0s
    ##  45700K .......... .......... .......... .......... .......... 57%  128M 0s
    ##  45750K .......... .......... .......... .......... .......... 57% 74.3M 0s
    ##  45800K .......... .......... .......... .......... .......... 57% 73.9M 0s
    ##  45850K .......... .......... .......... .......... .......... 57%  139M 0s
    ##  45900K .......... .......... .......... .......... .......... 57% 72.4M 0s
    ##  45950K .......... .......... .......... .......... .......... 57%  104M 0s
    ##  46000K .......... .......... .......... .......... .......... 57%  126M 0s
    ##  46050K .......... .......... .......... .......... .......... 57% 98.0M 0s
    ##  46100K .......... .......... .......... .......... .......... 57% 71.2M 0s
    ##  46150K .......... .......... .......... .......... .......... 57%  137M 0s
    ##  46200K .......... .......... .......... .......... .......... 57%  122M 0s
    ##  46250K .......... .......... .......... .......... .......... 57% 77.5M 0s
    ##  46300K .......... .......... .......... .......... .......... 57%  122M 0s
    ##  46350K .......... .......... .......... .......... .......... 58%  128M 0s
    ##  46400K .......... .......... .......... .......... .......... 58% 77.9M 0s
    ##  46450K .......... .......... .......... .......... .......... 58% 76.3M 0s
    ##  46500K .......... .......... .......... .......... .......... 58% 98.9M 0s
    ##  46550K .......... .......... .......... .......... .......... 58%  101M 0s
    ##  46600K .......... .......... .......... .......... .......... 58%  130M 0s
    ##  46650K .......... .......... .......... .......... .......... 58% 38.0M 0s
    ##  46700K .......... .......... .......... .......... .......... 58%  151M 0s
    ##  46750K .......... .......... .......... .......... .......... 58%  136M 0s
    ##  46800K .......... .......... .......... .......... .......... 58%  149M 0s
    ##  46850K .......... .......... .......... .......... .......... 58%  125M 0s
    ##  46900K .......... .......... .......... .......... .......... 58% 86.3M 0s
    ##  46950K .......... .......... .......... .......... .......... 58% 76.1M 0s
    ##  47000K .......... .......... .......... .......... .......... 58%  163M 0s
    ##  47050K .......... .......... .......... .......... .......... 58%  108M 0s
    ##  47100K .......... .......... .......... .......... .......... 58% 83.8M 0s
    ##  47150K .......... .......... .......... .......... .......... 59%  108M 0s
    ##  47200K .......... .......... .......... .......... .......... 59% 42.9M 0s
    ##  47250K .......... .......... .......... .......... .......... 59%  125M 0s
    ##  47300K .......... .......... .......... .......... .......... 59%  142M 0s
    ##  47350K .......... .......... .......... .......... .......... 59%  118M 0s
    ##  47400K .......... .......... .......... .......... .......... 59%  130M 0s
    ##  47450K .......... .......... .......... .......... .......... 59% 76.7M 0s
    ##  47500K .......... .......... .......... .......... .......... 59%  145M 0s
    ##  47550K .......... .......... .......... .......... .......... 59% 86.5M 0s
    ##  47600K .......... .......... .......... .......... .......... 59%  130M 0s
    ##  47650K .......... .......... .......... .......... .......... 59% 64.5M 0s
    ##  47700K .......... .......... .......... .......... .......... 59%  111M 0s
    ##  47750K .......... .......... .......... .......... .......... 59%  175M 0s
    ##  47800K .......... .......... .......... .......... .......... 59% 45.7M 0s
    ##  47850K .......... .......... .......... .......... .......... 59%  139M 0s
    ##  47900K .......... .......... .......... .......... .......... 59%  133M 0s
    ##  47950K .......... .......... .......... .......... .......... 60%  125M 0s
    ##  48000K .......... .......... .......... .......... .......... 60%  124M 0s
    ##  48050K .......... .......... .......... .......... .......... 60% 99.1M 0s
    ##  48100K .......... .......... .......... .......... .......... 60%  133M 0s
    ##  48150K .......... .......... .......... .......... .......... 60%  117M 0s
    ##  48200K .......... .......... .......... .......... .......... 60%  150M 0s
    ##  48250K .......... .......... .......... .......... .......... 60% 80.7M 0s
    ##  48300K .......... .......... .......... .......... .......... 60%  126M 0s
    ##  48350K .......... .......... .......... .......... .......... 60% 58.6M 0s
    ##  48400K .......... .......... .......... .......... .......... 60%  151M 0s
    ##  48450K .......... .......... .......... .......... .......... 60%  140M 0s
    ##  48500K .......... .......... .......... .......... .......... 60% 89.3M 0s
    ##  48550K .......... .......... .......... .......... .......... 60% 81.2M 0s
    ##  48600K .......... .......... .......... .......... .......... 60% 95.3M 0s
    ##  48650K .......... .......... .......... .......... .......... 60%  116M 0s
    ##  48700K .......... .......... .......... .......... .......... 60%  106M 0s
    ##  48750K .......... .......... .......... .......... .......... 61%  157M 0s
    ##  48800K .......... .......... .......... .......... .......... 61% 86.5M 0s
    ##  48850K .......... .......... .......... .......... .......... 61%  124M 0s
    ##  48900K .......... .......... .......... .......... .......... 61%  108M 0s
    ##  48950K .......... .......... .......... .......... .......... 61%  111M 0s
    ##  49000K .......... .......... .......... .......... .......... 61%  137M 0s
    ##  49050K .......... .......... .......... .......... .......... 61%  120M 0s
    ##  49100K .......... .......... .......... .......... .......... 61%  145M 0s
    ##  49150K .......... .......... .......... .......... .......... 61% 76.9M 0s
    ##  49200K .......... .......... .......... .......... .......... 61%  136M 0s
    ##  49250K .......... .......... .......... .......... .......... 61%  120M 0s
    ##  49300K .......... .......... .......... .......... .......... 61% 78.0M 0s
    ##  49350K .......... .......... .......... .......... .......... 61%  124M 0s
    ##  49400K .......... .......... .......... .......... .......... 61%  120M 0s
    ##  49450K .......... .......... .......... .......... .......... 61%  125M 0s
    ##  49500K .......... .......... .......... .......... .......... 61% 69.2M 0s
    ##  49550K .......... .......... .......... .......... .......... 62%  113M 0s
    ##  49600K .......... .......... .......... .......... .......... 62%  120M 0s
    ##  49650K .......... .......... .......... .......... .......... 62%  161M 0s
    ##  49700K .......... .......... .......... .......... .......... 62%  101M 0s
    ##  49750K .......... .......... .......... .......... .......... 62% 94.9M 0s
    ##  49800K .......... .......... .......... .......... .......... 62%  131M 0s
    ##  49850K .......... .......... .......... .......... .......... 62%  128M 0s
    ##  49900K .......... .......... .......... .......... .......... 62%  132M 0s
    ##  49950K .......... .......... .......... .......... .......... 62%  119M 0s
    ##  50000K .......... .......... .......... .......... .......... 62%  120M 0s
    ##  50050K .......... .......... .......... .......... .......... 62% 68.7M 0s
    ##  50100K .......... .......... .......... .......... .......... 62%  115M 0s
    ##  50150K .......... .......... .......... .......... .......... 62%  131M 0s
    ##  50200K .......... .......... .......... .......... .......... 62%  158M 0s
    ##  50250K .......... .......... .......... .......... .......... 62%  116M 0s
    ##  50300K .......... .......... .......... .......... .......... 62%  134M 0s
    ##  50350K .......... .......... .......... .......... .......... 63%  138M 0s
    ##  50400K .......... .......... .......... .......... .......... 63%  105M 0s
    ##  50450K .......... .......... .......... .......... .......... 63%  128M 0s
    ##  50500K .......... .......... .......... .......... .......... 63%  101M 0s
    ##  50550K .......... .......... .......... .......... .......... 63%  131M 0s
    ##  50600K .......... .......... .......... .......... .......... 63%  169M 0s
    ##  50650K .......... .......... .......... .......... .......... 63% 88.9M 0s
    ##  50700K .......... .......... .......... .......... .......... 63% 88.7M 0s
    ##  50750K .......... .......... .......... .......... .......... 63%  160M 0s
    ##  50800K .......... .......... .......... .......... .......... 63%  165M 0s
    ##  50850K .......... .......... .......... .......... .......... 63%  126M 0s
    ##  50900K .......... .......... .......... .......... .......... 63%  148M 0s
    ##  50950K .......... .......... .......... .......... .......... 63%  124M 0s
    ##  51000K .......... .......... .......... .......... .......... 63% 86.4M 0s
    ##  51050K .......... .......... .......... .......... .......... 63%  120M 0s
    ##  51100K .......... .......... .......... .......... .......... 63%  121M 0s
    ##  51150K .......... .......... .......... .......... .......... 64%  158M 0s
    ##  51200K .......... .......... .......... .......... .......... 64%  129M 0s
    ##  51250K .......... .......... .......... .......... .......... 64%  101M 0s
    ##  51300K .......... .......... .......... .......... .......... 64%  124M 0s
    ##  51350K .......... .......... .......... .......... .......... 64%  137M 0s
    ##  51400K .......... .......... .......... .......... .......... 64%  138M 0s
    ##  51450K .......... .......... .......... .......... .......... 64%  132M 0s
    ##  51500K .......... .......... .......... .......... .......... 64%  143M 0s
    ##  51550K .......... .......... .......... .......... .......... 64%  134M 0s
    ##  51600K .......... .......... .......... .......... .......... 64% 79.9M 0s
    ##  51650K .......... .......... .......... .......... .......... 64%  102M 0s
    ##  51700K .......... .......... .......... .......... .......... 64%  138M 0s
    ##  51750K .......... .......... .......... .......... .......... 64%  102M 0s
    ##  51800K .......... .......... .......... .......... .......... 64%  118M 0s
    ##  51850K .......... .......... .......... .......... .......... 64%  104M 0s
    ##  51900K .......... .......... .......... .......... .......... 65%  141M 0s
    ##  51950K .......... .......... .......... .......... .......... 65%  130M 0s
    ##  52000K .......... .......... .......... .......... .......... 65%  129M 0s
    ##  52050K .......... .......... .......... .......... .......... 65%  152M 0s
    ##  52100K .......... .......... .......... .......... .......... 65%  135M 0s
    ##  52150K .......... .......... .......... .......... .......... 65%  115M 0s
    ##  52200K .......... .......... .......... .......... .......... 65%  108M 0s
    ##  52250K .......... .......... .......... .......... .......... 65%  157M 0s
    ##  52300K .......... .......... .......... .......... .......... 65%  165M 0s
    ##  52350K .......... .......... .......... .......... .......... 65%  143M 0s
    ##  52400K .......... .......... .......... .......... .......... 65%  136M 0s
    ##  52450K .......... .......... .......... .......... .......... 65%  102M 0s
    ##  52500K .......... .......... .......... .......... .......... 65%  164M 0s
    ##  52550K .......... .......... .......... .......... .......... 65%  114M 0s
    ##  52600K .......... .......... .......... .......... .......... 65%  139M 0s
    ##  52650K .......... .......... .......... .......... .......... 65%  175M 0s
    ##  52700K .......... .......... .......... .......... .......... 66%  113M 0s
    ##  52750K .......... .......... .......... .......... .......... 66%  150M 0s
    ##  52800K .......... .......... .......... .......... .......... 66%  118M 0s
    ##  52850K .......... .......... .......... .......... .......... 66%  174M 0s
    ##  52900K .......... .......... .......... .......... .......... 66%  150M 0s
    ##  52950K .......... .......... .......... .......... .......... 66%  121M 0s
    ##  53000K .......... .......... .......... .......... .......... 66%  139M 0s
    ##  53050K .......... .......... .......... .......... .......... 66%  123M 0s
    ##  53100K .......... .......... .......... .......... .......... 66%  145M 0s
    ##  53150K .......... .......... .......... .......... .......... 66%  116M 0s
    ##  53200K .......... .......... .......... .......... .......... 66%  139M 0s
    ##  53250K .......... .......... .......... .......... .......... 66%  135M 0s
    ##  53300K .......... .......... .......... .......... .......... 66%  143M 0s
    ##  53350K .......... .......... .......... .......... .......... 66%  129M 0s
    ##  53400K .......... .......... .......... .......... .......... 66%  130M 0s
    ##  53450K .......... .......... .......... .......... .......... 66%  174M 0s
    ##  53500K .......... .......... .......... .......... .......... 67%  112M 0s
    ##  53550K .......... .......... .......... .......... .......... 67%  123M 0s
    ##  53600K .......... .......... .......... .......... .......... 67%  137M 0s
    ##  53650K .......... .......... .......... .......... .......... 67%  160M 0s
    ##  53700K .......... .......... .......... .......... .......... 67%  109M 0s
    ##  53750K .......... .......... .......... .......... .......... 67% 35.6M 0s
    ##  53800K .......... .......... .......... .......... .......... 67% 64.7M 0s
    ##  53850K .......... .......... .......... .......... .......... 67% 98.6M 0s
    ##  53900K .......... .......... .......... .......... .......... 67%  107M 0s
    ##  53950K .......... .......... .......... .......... .......... 67% 31.0M 0s
    ##  54000K .......... .......... .......... .......... .......... 67% 44.0M 0s
    ##  54050K .......... .......... .......... .......... .......... 67% 77.5M 0s
    ##  54100K .......... .......... .......... .......... .......... 67% 93.8M 0s
    ##  54150K .......... .......... .......... .......... .......... 67%  111M 0s
    ##  54200K .......... .......... .......... .......... .......... 67%  116M 0s
    ##  54250K .......... .......... .......... .......... .......... 67%  128M 0s
    ##  54300K .......... .......... .......... .......... .......... 68%  169M 0s
    ##  54350K .......... .......... .......... .......... .......... 68%  163M 0s
    ##  54400K .......... .......... .......... .......... .......... 68% 93.9M 0s
    ##  54450K .......... .......... .......... .......... .......... 68%  139M 0s
    ##  54500K .......... .......... .......... .......... .......... 68%  126M 0s
    ##  54550K .......... .......... .......... .......... .......... 68% 50.6M 0s
    ##  54600K .......... .......... .......... .......... .......... 68%  115M 0s
    ##  54650K .......... .......... .......... .......... .......... 68%  102M 0s
    ##  54700K .......... .......... .......... .......... .......... 68% 74.3M 0s
    ##  54750K .......... .......... .......... .......... .......... 68% 53.7M 0s
    ##  54800K .......... .......... .......... .......... .......... 68%  108M 0s
    ##  54850K .......... .......... .......... .......... .......... 68%  173M 0s
    ##  54900K .......... .......... .......... .......... .......... 68%  107M 0s
    ##  54950K .......... .......... .......... .......... .......... 68% 68.6M 0s
    ##  55000K .......... .......... .......... .......... .......... 68% 94.5M 0s
    ##  55050K .......... .......... .......... .......... .......... 68%  171M 0s
    ##  55100K .......... .......... .......... .......... .......... 69%  154M 0s
    ##  55150K .......... .......... .......... .......... .......... 69%  113M 0s
    ##  55200K .......... .......... .......... .......... .......... 69% 67.4M 0s
    ##  55250K .......... .......... .......... .......... .......... 69%  100M 0s
    ##  55300K .......... .......... .......... .......... .......... 69% 81.4M 0s
    ##  55350K .......... .......... .......... .......... .......... 69% 71.3M 0s
    ##  55400K .......... .......... .......... .......... .......... 69%  115M 0s
    ##  55450K .......... .......... .......... .......... .......... 69%  166M 0s
    ##  55500K .......... .......... .......... .......... .......... 69% 46.2M 0s
    ##  55550K .......... .......... .......... .......... .......... 69%  103M 0s
    ##  55600K .......... .......... .......... .......... .......... 69%  118M 0s
    ##  55650K .......... .......... .......... .......... .......... 69%  157M 0s
    ##  55700K .......... .......... .......... .......... .......... 69%  122M 0s
    ##  55750K .......... .......... .......... .......... .......... 69%  116M 0s
    ##  55800K .......... .......... .......... .......... .......... 69%  103M 0s
    ##  55850K .......... .......... .......... .......... .......... 69%  156M 0s
    ##  55900K .......... .......... .......... .......... .......... 70% 56.3M 0s
    ##  55950K .......... .......... .......... .......... .......... 70%  112M 0s
    ##  56000K .......... .......... .......... .......... .......... 70%  138M 0s
    ##  56050K .......... .......... .......... .......... .......... 70%  186M 0s
    ##  56100K .......... .......... .......... .......... .......... 70% 54.3M 0s
    ##  56150K .......... .......... .......... .......... .......... 70%  133M 0s
    ##  56200K .......... .......... .......... .......... .......... 70%  125M 0s
    ##  56250K .......... .......... .......... .......... .......... 70%  143M 0s
    ##  56300K .......... .......... .......... .......... .......... 70% 45.9M 0s
    ##  56350K .......... .......... .......... .......... .......... 70%  136M 0s
    ##  56400K .......... .......... .......... .......... .......... 70%  142M 0s
    ##  56450K .......... .......... .......... .......... .......... 70%  172M 0s
    ##  56500K .......... .......... .......... .......... .......... 70% 62.2M 0s
    ##  56550K .......... .......... .......... .......... .......... 70%  116M 0s
    ##  56600K .......... .......... .......... .......... .......... 70%  152M 0s
    ##  56650K .......... .......... .......... .......... .......... 70%  135M 0s
    ##  56700K .......... .......... .......... .......... .......... 71% 87.5M 0s
    ##  56750K .......... .......... .......... .......... .......... 71%  106M 0s
    ##  56800K .......... .......... .......... .......... .......... 71%  139M 0s
    ##  56850K .......... .......... .......... .......... .......... 71% 63.8M 0s
    ##  56900K .......... .......... .......... .......... .......... 71% 96.6M 0s
    ##  56950K .......... .......... .......... .......... .......... 71%  188M 0s
    ##  57000K .......... .......... .......... .......... .......... 71%  153M 0s
    ##  57050K .......... .......... .......... .......... .......... 71% 87.4M 0s
    ##  57100K .......... .......... .......... .......... .......... 71% 94.4M 0s
    ##  57150K .......... .......... .......... .......... .......... 71%  170M 0s
    ##  57200K .......... .......... .......... .......... .......... 71%  140M 0s
    ##  57250K .......... .......... .......... .......... .......... 71%  139M 0s
    ##  57300K .......... .......... .......... .......... .......... 71% 88.1M 0s
    ##  57350K .......... .......... .......... .......... .......... 71%  119M 0s
    ##  57400K .......... .......... .......... .......... .......... 71% 59.6M 0s
    ##  57450K .......... .......... .......... .......... .......... 71%  124M 0s
    ##  57500K .......... .......... .......... .......... .......... 72%  123M 0s
    ##  57550K .......... .......... .......... .......... .......... 72%  147M 0s
    ##  57600K .......... .......... .......... .......... .......... 72% 66.4M 0s
    ##  57650K .......... .......... .......... .......... .......... 72% 95.6M 0s
    ##  57700K .......... .......... .......... .......... .......... 72%  117M 0s
    ##  57750K .......... .......... .......... .......... .......... 72%  171M 0s
    ##  57800K .......... .......... .......... .......... .......... 72%  150M 0s
    ##  57850K .......... .......... .......... .......... .......... 72% 38.6M 0s
    ##  57900K .......... .......... .......... .......... .......... 72% 62.4M 0s
    ##  57950K .......... .......... .......... .......... .......... 72% 77.0M 0s
    ##  58000K .......... .......... .......... .......... .......... 72% 51.7M 0s
    ##  58050K .......... .......... .......... .......... .......... 72%  148M 0s
    ##  58100K .......... .......... .......... .......... .......... 72% 75.9M 0s
    ##  58150K .......... .......... .......... .......... .......... 72%  136M 0s
    ##  58200K .......... .......... .......... .......... .......... 72%  123M 0s
    ##  58250K .......... .......... .......... .......... .......... 72% 80.8M 0s
    ##  58300K .......... .......... .......... .......... .......... 73%  136M 0s
    ##  58350K .......... .......... .......... .......... .......... 73%  168M 0s
    ##  58400K .......... .......... .......... .......... .......... 73%  154M 0s
    ##  58450K .......... .......... .......... .......... .......... 73%  129M 0s
    ##  58500K .......... .......... .......... .......... .......... 73%  160M 0s
    ##  58550K .......... .......... .......... .......... .......... 73%  165M 0s
    ##  58600K .......... .......... .......... .......... .......... 73%  162M 0s
    ##  58650K .......... .......... .......... .......... .......... 73%  148M 0s
    ##  58700K .......... .......... .......... .......... .......... 73% 31.6M 0s
    ##  58750K .......... .......... .......... .......... .......... 73%  118M 0s
    ##  58800K .......... .......... .......... .......... .......... 73%  118M 0s
    ##  58850K .......... .......... .......... .......... .......... 73%  146M 0s
    ##  58900K .......... .......... .......... .......... .......... 73%  157M 0s
    ##  58950K .......... .......... .......... .......... .......... 73%  138M 0s
    ##  59000K .......... .......... .......... .......... .......... 73% 64.4M 0s
    ##  59050K .......... .......... .......... .......... .......... 73%  135M 0s
    ##  59100K .......... .......... .......... .......... .......... 74%  105M 0s
    ##  59150K .......... .......... .......... .......... .......... 74%  178M 0s
    ##  59200K .......... .......... .......... .......... .......... 74%  159M 0s
    ##  59250K .......... .......... .......... .......... .......... 74% 13.6M 0s
    ##  59300K .......... .......... .......... .......... .......... 74% 20.8M 0s
    ##  59350K .......... .......... .......... .......... .......... 74%  170M 0s
    ##  59400K .......... .......... .......... .......... .......... 74% 50.0M 0s
    ##  59450K .......... .......... .......... .......... .......... 74%  141M 0s
    ##  59500K .......... .......... .......... .......... .......... 74%  116M 0s
    ##  59550K .......... .......... .......... .......... .......... 74% 59.6M 0s
    ##  59600K .......... .......... .......... .......... .......... 74% 48.2M 0s
    ##  59650K .......... .......... .......... .......... .......... 74%  119M 0s
    ##  59700K .......... .......... .......... .......... .......... 74% 86.8M 0s
    ##  59750K .......... .......... .......... .......... .......... 74%  184M 0s
    ##  59800K .......... .......... .......... .......... .......... 74%  161M 0s
    ##  59850K .......... .......... .......... .......... .......... 74%  144M 0s
    ##  59900K .......... .......... .......... .......... .......... 75%  165M 0s
    ##  59950K .......... .......... .......... .......... .......... 75%  120M 0s
    ##  60000K .......... .......... .......... .......... .......... 75%  157M 0s
    ##  60050K .......... .......... .......... .......... .......... 75%  125M 0s
    ##  60100K .......... .......... .......... .......... .......... 75%  157M 0s
    ##  60150K .......... .......... .......... .......... .......... 75%  169M 0s
    ##  60200K .......... .......... .......... .......... .......... 75%  162M 0s
    ##  60250K .......... .......... .......... .......... .......... 75% 32.9M 0s
    ##  60300K .......... .......... .......... .......... .......... 75%  155M 0s
    ##  60350K .......... .......... .......... .......... .......... 75%  104M 0s
    ##  60400K .......... .......... .......... .......... .......... 75%  128M 0s
    ##  60450K .......... .......... .......... .......... .......... 75% 64.1M 0s
    ##  60500K .......... .......... .......... .......... .......... 75%  149M 0s
    ##  60550K .......... .......... .......... .......... .......... 75% 75.2M 0s
    ##  60600K .......... .......... .......... .......... .......... 75% 85.7M 0s
    ##  60650K .......... .......... .......... .......... .......... 75%  127M 0s
    ##  60700K .......... .......... .......... .......... .......... 76%  124M 0s
    ##  60750K .......... .......... .......... .......... .......... 76%  152M 0s
    ##  60800K .......... .......... .......... .......... .......... 76% 54.0M 0s
    ##  60850K .......... .......... .......... .......... .......... 76%  169M 0s
    ##  60900K .......... .......... .......... .......... .......... 76%  111M 0s
    ##  60950K .......... .......... .......... .......... .......... 76% 36.8M 0s
    ##  61000K .......... .......... .......... .......... .......... 76% 99.1M 0s
    ##  61050K .......... .......... .......... .......... .......... 76%  127M 0s
    ##  61100K .......... .......... .......... .......... .......... 76%  162M 0s
    ##  61150K .......... .......... .......... .......... .......... 76% 58.7M 0s
    ##  61200K .......... .......... .......... .......... .......... 76%  104M 0s
    ##  61250K .......... .......... .......... .......... .......... 76%  130M 0s
    ##  61300K .......... .......... .......... .......... .......... 76%  160M 0s
    ##  61350K .......... .......... .......... .......... .......... 76%  152M 0s
    ##  61400K .......... .......... .......... .......... .......... 76% 1.09M 0s
    ##  61450K .......... .......... .......... .......... .......... 76%  125M 0s
    ##  61500K .......... .......... .......... .......... .......... 77% 77.8M 0s
    ##  61550K .......... .......... .......... .......... .......... 77% 36.0M 0s
    ##  61600K .......... .......... .......... .......... .......... 77% 61.3M 0s
    ##  61650K .......... .......... .......... .......... .......... 77% 98.0M 0s
    ##  61700K .......... .......... .......... .......... .......... 77% 88.1M 0s
    ##  61750K .......... .......... .......... .......... .......... 77% 76.5M 0s
    ##  61800K .......... .......... .......... .......... .......... 77% 73.0M 0s
    ##  61850K .......... .......... .......... .......... .......... 77% 80.3M 0s
    ##  61900K .......... .......... .......... .......... .......... 77% 85.8M 0s
    ##  61950K .......... .......... .......... .......... .......... 77%  102M 0s
    ##  62000K .......... .......... .......... .......... .......... 77% 84.0M 0s
    ##  62050K .......... .......... .......... .......... .......... 77% 82.6M 0s
    ##  62100K .......... .......... .......... .......... .......... 77% 96.8M 0s
    ##  62150K .......... .......... .......... .......... .......... 77% 74.0M 0s
    ##  62200K .......... .......... .......... .......... .......... 77% 74.0M 0s
    ##  62250K .......... .......... .......... .......... .......... 77% 89.3M 0s
    ##  62300K .......... .......... .......... .......... .......... 78% 37.6M 0s
    ##  62350K .......... .......... .......... .......... .......... 78%  134M 0s
    ##  62400K .......... .......... .......... .......... .......... 78% 55.8M 0s
    ##  62450K .......... .......... .......... .......... .......... 78% 76.3M 0s
    ##  62500K .......... .......... .......... .......... .......... 78% 67.0M 0s
    ##  62550K .......... .......... .......... .......... .......... 78% 51.2M 0s
    ##  62600K .......... .......... .......... .......... .......... 78% 66.0M 0s
    ##  62650K .......... .......... .......... .......... .......... 78% 89.2M 0s
    ##  62700K .......... .......... .......... .......... .......... 78% 43.4M 0s
    ##  62750K .......... .......... .......... .......... .......... 78% 79.7M 0s
    ##  62800K .......... .......... .......... .......... .......... 78% 66.6M 0s
    ##  62850K .......... .......... .......... .......... .......... 78% 94.1M 0s
    ##  62900K .......... .......... .......... .......... .......... 78% 40.4M 0s
    ##  62950K .......... .......... .......... .......... .......... 78% 85.8M 0s
    ##  63000K .......... .......... .......... .......... .......... 78% 51.5M 0s
    ##  63050K .......... .......... .......... .......... .......... 78% 68.5M 0s
    ##  63100K .......... .......... .......... .......... .......... 79% 95.5M 0s
    ##  63150K .......... .......... .......... .......... .......... 79% 84.5M 0s
    ##  63200K .......... .......... .......... .......... .......... 79% 66.8M 0s
    ##  63250K .......... .......... .......... .......... .......... 79% 84.0M 0s
    ##  63300K .......... .......... .......... .......... .......... 79% 69.2M 0s
    ##  63350K .......... .......... .......... .......... .......... 79% 62.9M 0s
    ##  63400K .......... .......... .......... .......... .......... 79% 35.2M 0s
    ##  63450K .......... .......... .......... .......... .......... 79% 63.2M 0s
    ##  63500K .......... .......... .......... .......... .......... 79% 72.6M 0s
    ##  63550K .......... .......... .......... .......... .......... 79% 75.3M 0s
    ##  63600K .......... .......... .......... .......... .......... 79% 83.2M 0s
    ##  63650K .......... .......... .......... .......... .......... 79% 76.2M 0s
    ##  63700K .......... .......... .......... .......... .......... 79% 72.8M 0s
    ##  63750K .......... .......... .......... .......... .......... 79% 79.2M 0s
    ##  63800K .......... .......... .......... .......... .......... 79% 79.8M 0s
    ##  63850K .......... .......... .......... .......... .......... 79% 81.0M 0s
    ##  63900K .......... .......... .......... .......... .......... 80% 85.5M 0s
    ##  63950K .......... .......... .......... .......... .......... 80% 78.7M 0s
    ##  64000K .......... .......... .......... .......... .......... 80% 80.5M 0s
    ##  64050K .......... .......... .......... .......... .......... 80%  118M 0s
    ##  64100K .......... .......... .......... .......... .......... 80% 69.0M 0s
    ##  64150K .......... .......... .......... .......... .......... 80% 98.6M 0s
    ##  64200K .......... .......... .......... .......... .......... 80% 86.9M 0s
    ##  64250K .......... .......... .......... .......... .......... 80% 78.0M 0s
    ##  64300K .......... .......... .......... .......... .......... 80% 91.0M 0s
    ##  64350K .......... .......... .......... .......... .......... 80% 48.2M 0s
    ##  64400K .......... .......... .......... .......... .......... 80% 54.1M 0s
    ##  64450K .......... .......... .......... .......... .......... 80% 66.5M 0s
    ##  64500K .......... .......... .......... .......... .......... 80% 62.9M 0s
    ##  64550K .......... .......... .......... .......... .......... 80% 69.7M 0s
    ##  64600K .......... .......... .......... .......... .......... 80% 73.5M 0s
    ##  64650K .......... .......... .......... .......... .......... 80% 82.5M 0s
    ##  64700K .......... .......... .......... .......... .......... 81% 88.8M 0s
    ##  64750K .......... .......... .......... .......... .......... 81% 89.4M 0s
    ##  64800K .......... .......... .......... .......... .......... 81% 85.5M 0s
    ##  64850K .......... .......... .......... .......... .......... 81% 87.8M 0s
    ##  64900K .......... .......... .......... .......... .......... 81%  104M 0s
    ##  64950K .......... .......... .......... .......... .......... 81% 90.7M 0s
    ##  65000K .......... .......... .......... .......... .......... 81% 87.3M 0s
    ##  65050K .......... .......... .......... .......... .......... 81% 78.7M 0s
    ##  65100K .......... .......... .......... .......... .......... 81% 84.8M 0s
    ##  65150K .......... .......... .......... .......... .......... 81% 72.3M 0s
    ##  65200K .......... .......... .......... .......... .......... 81% 96.6M 0s
    ##  65250K .......... .......... .......... .......... .......... 81% 91.5M 0s
    ##  65300K .......... .......... .......... .......... .......... 81% 83.8M 0s
    ##  65350K .......... .......... .......... .......... .......... 81% 76.2M 0s
    ##  65400K .......... .......... .......... .......... .......... 81% 82.3M 0s
    ##  65450K .......... .......... .......... .......... .......... 81% 82.0M 0s
    ##  65500K .......... .......... .......... .......... .......... 82% 77.7M 0s
    ##  65550K .......... .......... .......... .......... .......... 82% 58.1M 0s
    ##  65600K .......... .......... .......... .......... .......... 82% 82.6M 0s
    ##  65650K .......... .......... .......... .......... .......... 82%  105M 0s
    ##  65700K .......... .......... .......... .......... .......... 82% 94.9M 0s
    ##  65750K .......... .......... .......... .......... .......... 82%  102M 0s
    ##  65800K .......... .......... .......... .......... .......... 82% 72.3M 0s
    ##  65850K .......... .......... .......... .......... .......... 82%  109M 0s
    ##  65900K .......... .......... .......... .......... .......... 82% 79.1M 0s
    ##  65950K .......... .......... .......... .......... .......... 82%  109M 0s
    ##  66000K .......... .......... .......... .......... .......... 82% 84.3M 0s
    ##  66050K .......... .......... .......... .......... .......... 82% 96.3M 0s
    ##  66100K .......... .......... .......... .......... .......... 82% 86.6M 0s
    ##  66150K .......... .......... .......... .......... .......... 82% 84.5M 0s
    ##  66200K .......... .......... .......... .......... .......... 82%  102M 0s
    ##  66250K .......... .......... .......... .......... .......... 82%  136M 0s
    ##  66300K .......... .......... .......... .......... .......... 83% 75.6M 0s
    ##  66350K .......... .......... .......... .......... .......... 83%  104M 0s
    ##  66400K .......... .......... .......... .......... .......... 83% 94.6M 0s
    ##  66450K .......... .......... .......... .......... .......... 83% 91.5M 0s
    ##  66500K .......... .......... .......... .......... .......... 83% 88.8M 0s
    ##  66550K .......... .......... .......... .......... .......... 83%  119M 0s
    ##  66600K .......... .......... .......... .......... .......... 83% 82.1M 0s
    ##  66650K .......... .......... .......... .......... .......... 83%  112M 0s
    ##  66700K .......... .......... .......... .......... .......... 83%  111M 0s
    ##  66750K .......... .......... .......... .......... .......... 83%  103M 0s
    ##  66800K .......... .......... .......... .......... .......... 83% 94.5M 0s
    ##  66850K .......... .......... .......... .......... .......... 83% 83.8M 0s
    ##  66900K .......... .......... .......... .......... .......... 83% 97.4M 0s
    ##  66950K .......... .......... .......... .......... .......... 83%  130M 0s
    ##  67000K .......... .......... .......... .......... .......... 83% 78.3M 0s
    ##  67050K .......... .......... .......... .......... .......... 83%  108M 0s
    ##  67100K .......... .......... .......... .......... .......... 84% 91.8M 0s
    ##  67150K .......... .......... .......... .......... .......... 84% 99.2M 0s
    ##  67200K .......... .......... .......... .......... .......... 84%  113M 0s
    ##  67250K .......... .......... .......... .......... .......... 84% 72.7M 0s
    ##  67300K .......... .......... .......... .......... .......... 84% 91.7M 0s
    ##  67350K .......... .......... .......... .......... .......... 84%  100M 0s
    ##  67400K .......... .......... .......... .......... .......... 84%  106M 0s
    ##  67450K .......... .......... .......... .......... .......... 84%  148M 0s
    ##  67500K .......... .......... .......... .......... .......... 84% 93.6M 0s
    ##  67550K .......... .......... .......... .......... .......... 84%  120M 0s
    ##  67600K .......... .......... .......... .......... .......... 84% 97.1M 0s
    ##  67650K .......... .......... .......... .......... .......... 84%  117M 0s
    ##  67700K .......... .......... .......... .......... .......... 84% 84.9M 0s
    ##  67750K .......... .......... .......... .......... .......... 84%  148M 0s
    ##  67800K .......... .......... .......... .......... .......... 84%  114M 0s
    ##  67850K .......... .......... .......... .......... .......... 84%  105M 0s
    ##  67900K .......... .......... .......... .......... .......... 85%  112M 0s
    ##  67950K .......... .......... .......... .......... .......... 85%  117M 0s
    ##  68000K .......... .......... .......... .......... .......... 85%  115M 0s
    ##  68050K .......... .......... .......... .......... .......... 85%  118M 0s
    ##  68100K .......... .......... .......... .......... .......... 85% 89.1M 0s
    ##  68150K .......... .......... .......... .......... .......... 85% 95.8M 0s
    ##  68200K .......... .......... .......... .......... .......... 85% 83.1M 0s
    ##  68250K .......... .......... .......... .......... .......... 85%  120M 0s
    ##  68300K .......... .......... .......... .......... .......... 85%  119M 0s
    ##  68350K .......... .......... .......... .......... .......... 85%  147M 0s
    ##  68400K .......... .......... .......... .......... .......... 85%  124M 0s
    ##  68450K .......... .......... .......... .......... .......... 85%  126M 0s
    ##  68500K .......... .......... .......... .......... .......... 85%  120M 0s
    ##  68550K .......... .......... .......... .......... .......... 85%  148M 0s
    ##  68600K .......... .......... .......... .......... .......... 85%  111M 0s
    ##  68650K .......... .......... .......... .......... .......... 85%  157M 0s
    ##  68700K .......... .......... .......... .......... .......... 86%  110M 0s
    ##  68750K .......... .......... .......... .......... .......... 86%  114M 0s
    ##  68800K .......... .......... .......... .......... .......... 86%  103M 0s
    ##  68850K .......... .......... .......... .......... .......... 86%  117M 0s
    ##  68900K .......... .......... .......... .......... .......... 86%  113M 0s
    ##  68950K .......... .......... .......... .......... .......... 86% 71.6M 0s
    ##  69000K .......... .......... .......... .......... .......... 86% 93.8M 0s
    ##  69050K .......... .......... .......... .......... .......... 86% 82.8M 0s
    ##  69100K .......... .......... .......... .......... .......... 86% 44.0M 0s
    ##  69150K .......... .......... .......... .......... .......... 86%  141M 0s
    ##  69200K .......... .......... .......... .......... .......... 86%  120M 0s
    ##  69250K .......... .......... .......... .......... .......... 86%  145M 0s
    ##  69300K .......... .......... .......... .......... .......... 86% 78.6M 0s
    ##  69350K .......... .......... .......... .......... .......... 86%  110M 0s
    ##  69400K .......... .......... .......... .......... .......... 86% 80.6M 0s
    ##  69450K .......... .......... .......... .......... .......... 86%  130M 0s
    ##  69500K .......... .......... .......... .......... .......... 87% 92.8M 0s
    ##  69550K .......... .......... .......... .......... .......... 87%  141M 0s
    ##  69600K .......... .......... .......... .......... .......... 87%  106M 0s
    ##  69650K .......... .......... .......... .......... .......... 87%  110M 0s
    ##  69700K .......... .......... .......... .......... .......... 87% 83.3M 0s
    ##  69750K .......... .......... .......... .......... .......... 87%  138M 0s
    ##  69800K .......... .......... .......... .......... .......... 87%  118M 0s
    ##  69850K .......... .......... .......... .......... .......... 87%  133M 0s
    ##  69900K .......... .......... .......... .......... .......... 87%  128M 0s
    ##  69950K .......... .......... .......... .......... .......... 87%  122M 0s
    ##  70000K .......... .......... .......... .......... .......... 87% 63.7M 0s
    ##  70050K .......... .......... .......... .......... .......... 87%  124M 0s
    ##  70100K .......... .......... .......... .......... .......... 87% 45.0M 0s
    ##  70150K .......... .......... .......... .......... .......... 87%  158M 0s
    ##  70200K .......... .......... .......... .......... .......... 87% 69.8M 0s
    ##  70250K .......... .......... .......... .......... .......... 87%  122M 0s
    ##  70300K .......... .......... .......... .......... .......... 88%  108M 0s
    ##  70350K .......... .......... .......... .......... .......... 88%  113M 0s
    ##  70400K .......... .......... .......... .......... .......... 88%  122M 0s
    ##  70450K .......... .......... .......... .......... .......... 88%  163M 0s
    ##  70500K .......... .......... .......... .......... .......... 88%  121M 0s
    ##  70550K .......... .......... .......... .......... .......... 88%  112M 0s
    ##  70600K .......... .......... .......... .......... .......... 88%  105M 0s
    ##  70650K .......... .......... .......... .......... .......... 88% 84.1M 0s
    ##  70700K .......... .......... .......... .......... .......... 88% 38.2M 0s
    ##  70750K .......... .......... .......... .......... .......... 88%  166M 0s
    ##  70800K .......... .......... .......... .......... .......... 88%  114M 0s
    ##  70850K .......... .......... .......... .......... .......... 88%  121M 0s
    ##  70900K .......... .......... .......... .......... .......... 88% 97.5M 0s
    ##  70950K .......... .......... .......... .......... .......... 88% 81.6M 0s
    ##  71000K .......... .......... .......... .......... .......... 88%  110M 0s
    ##  71050K .......... .......... .......... .......... .......... 88% 79.7M 0s
    ##  71100K .......... .......... .......... .......... .......... 89% 79.5M 0s
    ##  71150K .......... .......... .......... .......... .......... 89%  101M 0s
    ##  71200K .......... .......... .......... .......... .......... 89%  129M 0s
    ##  71250K .......... .......... .......... .......... .......... 89%  136M 0s
    ##  71300K .......... .......... .......... .......... .......... 89% 57.6M 0s
    ##  71350K .......... .......... .......... .......... .......... 89%  137M 0s
    ##  71400K .......... .......... .......... .......... .......... 89%  122M 0s
    ##  71450K .......... .......... .......... .......... .......... 89% 84.9M 0s
    ##  71500K .......... .......... .......... .......... .......... 89% 61.5M 0s
    ##  71550K .......... .......... .......... .......... .......... 89%  125M 0s
    ##  71600K .......... .......... .......... .......... .......... 89%  138M 0s
    ##  71650K .......... .......... .......... .......... .......... 89% 82.4M 0s
    ##  71700K .......... .......... .......... .......... .......... 89% 89.7M 0s
    ##  71750K .......... .......... .......... .......... .......... 89%  122M 0s
    ##  71800K .......... .......... .......... .......... .......... 89%  134M 0s
    ##  71850K .......... .......... .......... .......... .......... 89% 78.8M 0s
    ##  71900K .......... .......... .......... .......... .......... 90%  106M 0s
    ##  71950K .......... .......... .......... .......... .......... 90%  155M 0s
    ##  72000K .......... .......... .......... .......... .......... 90%  146M 0s
    ##  72050K .......... .......... .......... .......... .......... 90% 49.9M 0s
    ##  72100K .......... .......... .......... .......... .......... 90%  136M 0s
    ##  72150K .......... .......... .......... .......... .......... 90%  161M 0s
    ##  72200K .......... .......... .......... .......... .......... 90%  133M 0s
    ##  72250K .......... .......... .......... .......... .......... 90%  125M 0s
    ##  72300K .......... .......... .......... .......... .......... 90%  104M 0s
    ##  72350K .......... .......... .......... .......... .......... 90%  100M 0s
    ##  72400K .......... .......... .......... .......... .......... 90% 31.0M 0s
    ##  72450K .......... .......... .......... .......... .......... 90% 80.5M 0s
    ##  72500K .......... .......... .......... .......... .......... 90% 50.0M 0s
    ##  72550K .......... .......... .......... .......... .......... 90% 42.2M 0s
    ##  72600K .......... .......... .......... .......... .......... 90% 44.4M 0s
    ##  72650K .......... .......... .......... .......... .......... 90% 35.1M 0s
    ##  72700K .......... .......... .......... .......... .......... 91% 88.6M 0s
    ##  72750K .......... .......... .......... .......... .......... 91% 76.6M 0s
    ##  72800K .......... .......... .......... .......... .......... 91% 34.5M 0s
    ##  72850K .......... .......... .......... .......... .......... 91% 43.5M 0s
    ##  72900K .......... .......... .......... .......... .......... 91% 58.5M 0s
    ##  72950K .......... .......... .......... .......... .......... 91% 46.6M 0s
    ##  73000K .......... .......... .......... .......... .......... 91% 44.1M 0s
    ##  73050K .......... .......... .......... .......... .......... 91% 50.1M 0s
    ##  73100K .......... .......... .......... .......... .......... 91% 66.6M 0s
    ##  73150K .......... .......... .......... .......... .......... 91% 45.6M 0s
    ##  73200K .......... .......... .......... .......... .......... 91% 49.6M 0s
    ##  73250K .......... .......... .......... .......... .......... 91% 50.8M 0s
    ##  73300K .......... .......... .......... .......... .......... 91% 46.8M 0s
    ##  73350K .......... .......... .......... .......... .......... 91% 51.4M 0s
    ##  73400K .......... .......... .......... .......... .......... 91% 40.5M 0s
    ##  73450K .......... .......... .......... .......... .......... 91% 45.8M 0s
    ##  73500K .......... .......... .......... .......... .......... 92% 43.2M 0s
    ##  73550K .......... .......... .......... .......... .......... 92% 51.8M 0s
    ##  73600K .......... .......... .......... .......... .......... 92% 47.7M 0s
    ##  73650K .......... .......... .......... .......... .......... 92% 41.7M 0s
    ##  73700K .......... .......... .......... .......... .......... 92% 43.8M 0s
    ##  73750K .......... .......... .......... .......... .......... 92% 71.5M 0s
    ##  73800K .......... .......... .......... .......... .......... 92% 43.5M 0s
    ##  73850K .......... .......... .......... .......... .......... 92% 59.3M 0s
    ##  73900K .......... .......... .......... .......... .......... 92% 45.4M 0s
    ##  73950K .......... .......... .......... .......... .......... 92%  117M 0s
    ##  74000K .......... .......... .......... .......... .......... 92% 69.0M 0s
    ##  74050K .......... .......... .......... .......... .......... 92% 78.0M 0s
    ##  74100K .......... .......... .......... .......... .......... 92% 67.1M 0s
    ##  74150K .......... .......... .......... .......... .......... 92% 86.6M 0s
    ##  74200K .......... .......... .......... .......... .......... 92% 54.9M 0s
    ##  74250K .......... .......... .......... .......... .......... 92%  112M 0s
    ##  74300K .......... .......... .......... .......... .......... 93% 74.6M 0s
    ##  74350K .......... .......... .......... .......... .......... 93% 54.1M 0s
    ##  74400K .......... .......... .......... .......... .......... 93% 61.9M 0s
    ##  74450K .......... .......... .......... .......... .......... 93% 93.7M 0s
    ##  74500K .......... .......... .......... .......... .......... 93% 70.4M 0s
    ##  74550K .......... .......... .......... .......... .......... 93%  104M 0s
    ##  74600K .......... .......... .......... .......... .......... 93% 71.0M 0s
    ##  74650K .......... .......... .......... .......... .......... 93%  117M 0s
    ##  74700K .......... .......... .......... .......... .......... 93%  101M 0s
    ##  74750K .......... .......... .......... .......... .......... 93%  106M 0s
    ##  74800K .......... .......... .......... .......... .......... 93%  104M 0s
    ##  74850K .......... .......... .......... .......... .......... 93%  118M 0s
    ##  74900K .......... .......... .......... .......... .......... 93%  104M 0s
    ##  74950K .......... .......... .......... .......... .......... 93%  122M 0s
    ##  75000K .......... .......... .......... .......... .......... 93% 92.3M 0s
    ##  75050K .......... .......... .......... .......... .......... 93%  109M 0s
    ##  75100K .......... .......... .......... .......... .......... 94%  104M 0s
    ##  75150K .......... .......... .......... .......... .......... 94%  128M 0s
    ##  75200K .......... .......... .......... .......... .......... 94%  105M 0s
    ##  75250K .......... .......... .......... .......... .......... 94%  124M 0s
    ##  75300K .......... .......... .......... .......... .......... 94% 96.0M 0s
    ##  75350K .......... .......... .......... .......... .......... 94%  109M 0s
    ##  75400K .......... .......... .......... .......... .......... 94%  102M 0s
    ##  75450K .......... .......... .......... .......... .......... 94%  127M 0s
    ##  75500K .......... .......... .......... .......... .......... 94%  106M 0s
    ##  75550K .......... .......... .......... .......... .......... 94%  114M 0s
    ##  75600K .......... .......... .......... .......... .......... 94%  105M 0s
    ##  75650K .......... .......... .......... .......... .......... 94%  114M 0s
    ##  75700K .......... .......... .......... .......... .......... 94%  112M 0s
    ##  75750K .......... .......... .......... .......... .......... 94%  131M 0s
    ##  75800K .......... .......... .......... .......... .......... 94%  108M 0s
    ##  75850K .......... .......... .......... .......... .......... 94%  117M 0s
    ##  75900K .......... .......... .......... .......... .......... 95%  111M 0s
    ##  75950K .......... .......... .......... .......... .......... 95%  112M 0s
    ##  76000K .......... .......... .......... .......... .......... 95%  109M 0s
    ##  76050K .......... .......... .......... .......... .......... 95%  141M 0s
    ##  76100K .......... .......... .......... .......... .......... 95%  114M 0s
    ##  76150K .......... .......... .......... .......... .......... 95%  128M 0s
    ##  76200K .......... .......... .......... .......... .......... 95%  109M 0s
    ##  76250K .......... .......... .......... .......... .......... 95%  111M 0s
    ##  76300K .......... .......... .......... .......... .......... 95%  114M 0s
    ##  76350K .......... .......... .......... .......... .......... 95%  130M 0s
    ##  76400K .......... .......... .......... .......... .......... 95%  109M 0s
    ##  76450K .......... .......... .......... .......... .......... 95%  134M 0s
    ##  76500K .......... .......... .......... .......... .......... 95%  107M 0s
    ##  76550K .......... .......... .......... .......... .......... 95% 93.3M 0s
    ##  76600K .......... .......... .......... .......... .......... 95%  109M 0s
    ##  76650K .......... .......... .......... .......... .......... 95%  139M 0s
    ##  76700K .......... .......... .......... .......... .......... 96%  102M 0s
    ##  76750K .......... .......... .......... .......... .......... 96%  138M 0s
    ##  76800K .......... .......... .......... .......... .......... 96% 1.04M 0s
    ##  76850K .......... .......... .......... .......... .......... 96% 76.5M 0s
    ##  76900K .......... .......... .......... .......... .......... 96% 46.9M 0s
    ##  76950K .......... .......... .......... .......... .......... 96%  107M 0s
    ##  77000K .......... .......... .......... .......... .......... 96% 57.9M 0s
    ##  77050K .......... .......... .......... .......... .......... 96%  106M 0s
    ##  77100K .......... .......... .......... .......... .......... 96% 18.3M 0s
    ##  77150K .......... .......... .......... .......... .......... 96% 90.1M 0s
    ##  77200K .......... .......... .......... .......... .......... 96% 46.8M 0s
    ##  77250K .......... .......... .......... .......... .......... 96%  115M 0s
    ##  77300K .......... .......... .......... .......... .......... 96% 36.5M 0s
    ##  77350K .......... .......... .......... .......... .......... 96% 55.8M 0s
    ##  77400K .......... .......... .......... .......... .......... 96% 95.0M 0s
    ##  77450K .......... .......... .......... .......... .......... 96%  104M 0s
    ##  77500K .......... .......... .......... .......... .......... 97% 87.9M 0s
    ##  77550K .......... .......... .......... .......... .......... 97%  117M 0s
    ##  77600K .......... .......... .......... .......... .......... 97% 96.3M 0s
    ##  77650K .......... .......... .......... .......... .......... 97%  111M 0s
    ##  77700K .......... .......... .......... .......... .......... 97% 95.4M 0s
    ##  77750K .......... .......... .......... .......... .......... 97% 16.9M 0s
    ##  77800K .......... .......... .......... .......... .......... 97% 70.2M 0s
    ##  77850K .......... .......... .......... .......... .......... 97% 41.0M 0s
    ##  77900K .......... .......... .......... .......... .......... 97% 58.1M 0s
    ##  77950K .......... .......... .......... .......... .......... 97% 46.3M 0s
    ##  78000K .......... .......... .......... .......... .......... 97% 49.2M 0s
    ##  78050K .......... .......... .......... .......... .......... 97% 42.5M 0s
    ##  78100K .......... .......... .......... .......... .......... 97% 47.6M 0s
    ##  78150K .......... .......... .......... .......... .......... 97% 62.5M 0s
    ##  78200K .......... .......... .......... .......... .......... 97% 89.4M 0s
    ##  78250K .......... .......... .......... .......... .......... 97% 18.5M 0s
    ##  78300K .......... .......... .......... .......... .......... 98% 93.4M 0s
    ##  78350K .......... .......... .......... .......... .......... 98% 90.4M 0s
    ##  78400K .......... .......... .......... .......... .......... 98%  102M 0s
    ##  78450K .......... .......... .......... .......... .......... 98% 27.3M 0s
    ##  78500K .......... .......... .......... .......... .......... 98% 72.0M 0s
    ##  78550K .......... .......... .......... .......... .......... 98% 83.1M 0s
    ##  78600K .......... .......... .......... .......... .......... 98% 47.4M 0s
    ##  78650K .......... .......... .......... .......... .......... 98% 72.9M 0s
    ##  78700K .......... .......... .......... .......... .......... 98% 44.9M 0s
    ##  78750K .......... .......... .......... .......... .......... 98% 72.2M 0s
    ##  78800K .......... .......... .......... .......... .......... 98% 89.0M 0s
    ##  78850K .......... .......... .......... .......... .......... 98% 59.0M 0s
    ##  78900K .......... .......... .......... .......... .......... 98% 48.7M 0s
    ##  78950K .......... .......... .......... .......... .......... 98%  101M 0s
    ##  79000K .......... .......... .......... .......... .......... 98% 97.2M 0s
    ##  79050K .......... .......... .......... .......... .......... 98%  127M 0s
    ##  79100K .......... .......... .......... .......... .......... 99%  104M 0s
    ##  79150K .......... .......... .......... .......... .......... 99%  125M 0s
    ##  79200K .......... .......... .......... .......... .......... 99% 89.7M 0s
    ##  79250K .......... .......... .......... .......... .......... 99%  113M 0s
    ##  79300K .......... .......... .......... .......... .......... 99%  102M 0s
    ##  79350K .......... .......... .......... .......... .......... 99% 32.1M 0s
    ##  79400K .......... .......... .......... .......... .......... 99% 97.4M 0s
    ##  79450K .......... .......... .......... .......... .......... 99%  124M 0s
    ##  79500K .......... .......... .......... .......... .......... 99% 92.0M 0s
    ##  79550K .......... .......... .......... .......... .......... 99%  109M 0s
    ##  79600K .......... .......... .......... .......... .......... 99%  101M 0s
    ##  79650K .......... .......... .......... .......... .......... 99%  134M 0s
    ##  79700K .......... .......... .......... .......... .......... 99%  107M 0s
    ##  79750K .......... .......... .......... .......... .......... 99%  126M 0s
    ##  79800K .......... .......... .......... .......... .......... 99%  105M 0s
    ##  79850K .......... .......... .......... .......... .......... 99%  106M 0s
    ##  79900K .......... .......... ..                              100%  112M=0.9s
    ## 
    ## 2020-12-04 08:46:33 (82.2 MB/s) - ‘silva_species_assignment_v138.fa.gz.3’ saved [81840166/81840166]

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
