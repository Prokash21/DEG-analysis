# DEG-analysis
RNA seq analysis
R version 4.2.3 (2023-03-15 ucrt) -- "Shortstop Beagle"
Copyright (C) 2023 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

[Workspace loaded from ~/.RData]

> if (!requireNamespace("BiocManager", quietly = TRUE))
+   install.packages("BiocManager")
> url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
> utils::download.file(url, destfile = "GSE63310_RAW.tar", mode = "wb")
trying URL 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file'
Content type 'application/x-tar' length 1996800 bytes (1.9 MB)
downloaded 1.9 MB

> utils::untar("GSE63310_RAW.tar", exdir = ".")
> files_gz <- Sys.glob("GSM*txt.gz")
> for(f in files_gz)
+   R.utils::gunzip(f, overwrite = TRUE)
> library("limma")   # Linear models for differential expression
> library("Glimma")  # Interactive plots for exploration
> library("edgeR")   # Process count data from NGS experiments
> library("Mus.musculus") # Gene annotations for the Mus musculus genome
Loading required package: AnnotationDbi
Loading required package: stats4
Loading required package: BiocGenerics

Attaching package: ‘BiocGenerics’

The following object is masked from ‘package:limma’:

    plotMA

The following objects are masked from ‘package:stats’:

    IQR, mad, sd, var, xtabs

The following objects are masked from ‘package:base’:

    anyDuplicated, aperm, append, as.data.frame, basename, cbind, colnames, dirname, do.call,
    duplicated, eval, evalq, Filter, Find, get, grep, grepl, intersect, is.unsorted, lapply,
    Map, mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply, union, unique, unsplit,
    which.max, which.min

Loading required package: Biobase
Welcome to Bioconductor

    Vignettes contain introductory material; view with 'browseVignettes()'. To cite
    Bioconductor, see 'citation("Biobase")', and for packages 'citation("pkgname")'.

Loading required package: IRanges
Loading required package: S4Vectors

Attaching package: ‘S4Vectors’

The following objects are masked from ‘package:base’:

    expand.grid, I, unname


Attaching package: ‘IRanges’

The following object is masked from ‘package:grDevices’:

    windows

Loading required package: OrganismDbi
Loading required package: GenomicFeatures
Loading required package: GenomeInfoDb
Loading required package: GenomicRanges
Loading required package: GO.db

Loading required package: org.Mm.eg.db

Loading required package: TxDb.Mmusculus.UCSC.mm10.knownGene
> files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt",
+            "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt",
+            "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
+            "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt",
+            "GSM1545545_JMS9-P8c.txt")
> read.delim(files[1], nrow = 5)
   EntrezID GeneLength Count
1    497097       3634     1
2 100503874       3259     0
3 100038431       1634     0
4     19888       9747     0
5     20671       3130     1
> x <- readDGE(files, columns = c(1, 3))
> class(x)
[1] "DGEList"
attr(,"package")
[1] "edgeR"
> dim(x)
[1] 27179     9
> names(x)
[1] "samples" "counts" 
> x$samples
                                        files group lib.size norm.factors
GSM1545535_10_6_5_11 GSM1545535_10_6_5_11.txt     1 32863052            1
GSM1545536_9_6_5_11   GSM1545536_9_6_5_11.txt     1 35335491            1
GSM1545538_purep53     GSM1545538_purep53.txt     1 57160817            1
GSM1545539_JMS8-2       GSM1545539_JMS8-2.txt     1 51368625            1
GSM1545540_JMS8-3       GSM1545540_JMS8-3.txt     1 75795034            1
GSM1545541_JMS8-4       GSM1545541_JMS8-4.txt     1 60517657            1
GSM1545542_JMS8-5       GSM1545542_JMS8-5.txt     1 55086324            1
GSM1545544_JMS9-P7c   GSM1545544_JMS9-P7c.txt     1 21311068            1
GSM1545545_JMS9-P8c   GSM1545545_JMS9-P8c.txt     1 19958838            1
> samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
> samplenames
[1] "10_6_5_11" "9_6_5_11"  "purep53"   "JMS8-2"    "JMS8-3"    "JMS8-4"    "JMS8-5"    "JMS9-P7c" 
[9] "JMS9-P8c" 
> colnames(x) <- samplenames
> group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP",
+                      "Basal", "ML", "LP"))
> x$samples$group <- group
> lane <- as.factor(rep(c("L004", "L006", "L008"), c(3, 4, 2)))
> x$samples$lane <- lane
> x$samples
                             files group lib.size norm.factors lane
10_6_5_11 GSM1545535_10_6_5_11.txt    LP 32863052            1 L004
9_6_5_11   GSM1545536_9_6_5_11.txt    ML 35335491            1 L004
purep53     GSM1545538_purep53.txt Basal 57160817            1 L004
JMS8-2       GSM1545539_JMS8-2.txt Basal 51368625            1 L006
JMS8-3       GSM1545540_JMS8-3.txt    ML 75795034            1 L006
JMS8-4       GSM1545541_JMS8-4.txt    LP 60517657            1 L006
JMS8-5       GSM1545542_JMS8-5.txt Basal 55086324            1 L006
JMS9-P7c   GSM1545544_JMS9-P7c.txt    ML 21311068            1 L008
JMS9-P8c   GSM1545545_JMS9-P8c.txt    LP 19958838            1 L008
> View(x)
> View(x)
> head(x$counts)
           Samples
Tags        10_6_5_11 9_6_5_11 purep53 JMS8-2 JMS8-3 JMS8-4 JMS8-5 JMS9-P7c JMS9-P8c
  497097            1        2     342    526      3      3    535        2        0
  100503874         0        0       5      6      0      0      5        0        0
  100038431         0        0       0      0      0      0      1        0        0
  19888             0        1       0      0     17      2      0        1        0
  20671             1        1      76     40     33     14     98       18        8
  27395           431      771    1368   1268   1564    769    818      468      342
> dim(x$counts)
[1] 27179     9
> geneid <- rownames(x)
> genes <- select(Mus.musculus, keys = geneid,
+                 columns = c("SYMBOL", "TXCHROM"),
+                 keytype = "ENTREZID")
'select()' returned 1:many mapping between keys and columns
> head(genes)
   ENTREZID  SYMBOL TXCHROM
1    497097    Xkr4    chr1
2 100503874 Gm19938    <NA>
3 100038431 Gm10568    <NA>
4     19888     Rp1    chr1
5     20671   Sox17    chr1
6     27395  Mrpl15    chr1
> genes <- genes[!duplicated(genes$ENTREZID), ]
> x$genes <- genes
> x
An object of class "DGEList"
$samples
                             files group lib.size norm.factors lane
10_6_5_11 GSM1545535_10_6_5_11.txt    LP 32863052            1 L004
9_6_5_11   GSM1545536_9_6_5_11.txt    ML 35335491            1 L004
purep53     GSM1545538_purep53.txt Basal 57160817            1 L004
JMS8-2       GSM1545539_JMS8-2.txt Basal 51368625            1 L006
JMS8-3       GSM1545540_JMS8-3.txt    ML 75795034            1 L006
JMS8-4       GSM1545541_JMS8-4.txt    LP 60517657            1 L006
JMS8-5       GSM1545542_JMS8-5.txt Basal 55086324            1 L006
JMS9-P7c   GSM1545544_JMS9-P7c.txt    ML 21311068            1 L008
JMS9-P8c   GSM1545545_JMS9-P8c.txt    LP 19958838            1 L008

$counts
           Samples
Tags        10_6_5_11 9_6_5_11 purep53 JMS8-2 JMS8-3 JMS8-4 JMS8-5 JMS9-P7c JMS9-P8c
  497097            1        2     342    526      3      3    535        2        0
  100503874         0        0       5      6      0      0      5        0        0
  100038431         0        0       0      0      0      0      1        0        0
  19888             0        1       0      0     17      2      0        1        0
  20671             1        1      76     40     33     14     98       18        8
27174 more rows ...

$genes
   ENTREZID  SYMBOL TXCHROM
1    497097    Xkr4    chr1
2 100503874 Gm19938    <NA>
3 100038431 Gm10568    <NA>
4     19888     Rp1    chr1
5     20671   Sox17    chr1
27174 more rows ...

> cpm <- cpm(x)
> lcpm <- cpm(x, log = TRUE)
> table(rowSums(x$counts == 0) == 9)

FALSE  TRUE 
22026  5153 
> plotDensities(lcpm, legend = FALSE, main = "Before filtering")
> abline(v = 0, lty = 3)
> keep.exprs <- rowSums(cpm > 1) >= 3
> x <- x[keep.exprs, , keep.lib.sizes = FALSE]
> dim(x)
[1] 14165     9
> lcpm <- cpm(x, log=TRUE)
> plotDensities(lcpm, legend = FALSE, main = "After filtering")
> abline(v = 0, lty = 3)
> x <- calcNormFactors(x, method = "TMM")
> x$samples$norm.factors
[1] 0.8957309 1.0349196 1.0439552 1.0405040 1.0323599 0.9223424 0.9836603 1.0827381 0.9792607
> x2 <- x
> x2$samples$norm.factors <- 1
> x2$counts[,1] <- ceiling(x2$counts[, 1] * 0.05)
> x2$counts[,2] <- x2$counts[, 2] * 5
> lcpm2 <- cpm(x2, log = TRUE)
> boxplot(lcpm2, las = 2, main = "Before normalization")
> x2 <- calcNormFactors(x2)
> x2$samples$norm.factors
[1] 0.05472223 6.13059440 1.22927355 1.17051887 1.21487709 1.05622968 1.14587663 1.26129350 1.11702264
> lcpm2 <- cpm(x2, log = TRUE)
> boxplot(lcpm2, las=2, main = "After normalization")
> library("RColorBrewer")
> group
[1] LP    ML    Basal Basal ML    LP    Basal ML    LP   
Levels: Basal LP ML
> col.group <- group
> levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
> col.group <- as.character(col.group)
> lane
[1] L004 L004 L004 L006 L006 L006 L006 L008 L008
Levels: L004 L006 L008
> col.lane <- lane
> levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
> col.lane <- as.character(col.lane)
> plotMDS(lcpm, labels = group, col = col.group,
+         main = "group")
> plotMDS(lcpm, labels = lane, col = col.lane, dim = c(3, 4),
+         main = "lane")
> glMDSPlot(lcpm, labels = paste(group, lane, sep = "_"),
+           groups = x$samples[, c(2, 5)], launch = TRUE)
> design <- model.matrix(~0 + group + lane)
> colnames(design) <- gsub("group", "", colnames(design))
> design
  Basal LP ML laneL006 laneL008
1     0  1  0        0        0
2     0  0  1        0        0
3     1  0  0        0        0
4     1  0  0        1        0
5     0  0  1        1        0
6     0  1  0        1        0
7     1  0  0        1        0
8     0  0  1        0        1
9     0  1  0        0        1
attr(,"assign")
[1] 1 1 1 2 2
attr(,"contrasts")
attr(,"contrasts")$group
[1] "contr.treatment"

attr(,"contrasts")$lane
[1] "contr.treatment"

> contr.matrix <- makeContrasts(BasalvsLP = Basal - LP,
+                               BasalvsML = Basal - ML,
+                               LPvsML = LP - ML,
+                               levels = colnames(design))
> contr.matrix
          Contrasts
Levels     BasalvsLP BasalvsML LPvsML
  Basal            1         1      0
  LP              -1         0      1
  ML               0        -1     -1
  laneL006         0         0      0
  laneL008         0         0      0
> v <- voom(x, design, plot = TRUE)
> v
An object of class "EList"
$genes
   ENTREZID SYMBOL TXCHROM
1    497097   Xkr4    chr1
6     27395 Mrpl15    chr1
7     18777 Lypla1    chr1
9     21399  Tcea1    chr1
10    58175  Rgs20    chr1
14160 more rows ...

$targets
                             files group lib.size norm.factors lane
10_6_5_11 GSM1545535_10_6_5_11.txt    LP 29409426    0.8957309 L004
9_6_5_11   GSM1545536_9_6_5_11.txt    ML 36528591    1.0349196 L004
purep53     GSM1545538_purep53.txt Basal 59598629    1.0439552 L004
JMS8-2       GSM1545539_JMS8-2.txt Basal 53382070    1.0405040 L006
JMS8-3       GSM1545540_JMS8-3.txt    ML 78175314    1.0323599 L006
JMS8-4       GSM1545541_JMS8-4.txt    LP 55762781    0.9223424 L006
JMS8-5       GSM1545542_JMS8-5.txt Basal 54115150    0.9836603 L006
JMS9-P7c   GSM1545544_JMS9-P7c.txt    ML 23043111    1.0827381 L008
JMS9-P8c   GSM1545545_JMS9-P8c.txt    LP 19525423    0.9792607 L008

$E
        Samples
Tags     10_6_5_11  9_6_5_11   purep53    JMS8-2    JMS8-3    JMS8-4    JMS8-5  JMS9-P7c  JMS9-P8c
  497097 -4.293244 -3.869026  2.522753  3.302006 -4.481286 -3.993876  3.306782 -3.204336 -5.287282
  27395   3.875010  4.400568  4.521172  4.570624  4.322845  3.786547  3.918878  4.345642  4.132678
  18777   4.707695  5.559334  5.400569  5.171235  5.627798  5.081794  5.080061  5.757404  5.150470
  21399   4.784462  4.741999  5.374548  5.130925  4.848030  4.944024  5.158292  5.036933  4.987679
  58175   3.943567  3.294875 -1.767924 -1.880302  2.993289  3.357379 -2.114104  3.142621  3.523290
14160 more rows ...

$weights
          [,1]      [,2]      [,3]     [,4]      [,5]      [,6]      [,7]      [,8]      [,9]
[1,]  1.183974  1.183974 20.526779 20.97747  1.773562  1.217142 21.125740  1.183974  1.183974
[2,] 20.879554 26.561871 31.596323 29.66102 32.558344 26.745293 29.792090 21.900102 17.150677
[3,] 28.003202 33.695540 34.845507 34.45673 35.148529 33.550527 34.517259 31.440457 25.228325
[4,] 27.670233 29.595778 34.901302 34.43298 34.841349 33.159425 34.493456 26.136796 24.502247
[5,] 19.737381 18.658333  3.184207  2.62986 24.191635 24.014937  2.648747 13.149278 14.351930
14160 more rows ...

$design
  Basal LP ML laneL006 laneL008
1     0  1  0        0        0
2     0  0  1        0        0
3     1  0  0        0        0
4     1  0  0        1        0
5     0  0  1        1        0
6     0  1  0        1        0
7     1  0  0        1        0
8     0  0  1        0        1
9     0  1  0        0        1
attr(,"assign")
[1] 1 1 1 2 2
attr(,"contrasts")
attr(,"contrasts")$group
[1] "contr.treatment"

attr(,"contrasts")$lane
[1] "contr.treatment"


> vfit <- lmFit(v, design)
> vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
> efit <- eBayes(vfit)
> plotSA(efit, main = "Final model: Mean−variance trend")
> summary(decideTests(efit))
       BasalvsLP BasalvsML LPvsML
Down        4127      4338   2895
NotSig      5740      5655   8825
Up          4298      4172   2445
> tfit <- treat(vfit, lfc = 1)
> dt <- decideTests(tfit)
> summary(dt)
       BasalvsLP BasalvsML LPvsML
Down        1417      1512    203
NotSig     11030     10895  13780
Up          1718      1758    182
> head(dt)
TestResults matrix
        Contrasts
         BasalvsLP BasalvsML LPvsML
  497097         1         1      0
  27395          0         0      0
  18777          0         0      0
  21399          0         0      0
  58175         -1        -1      0
  108664         0         0      0
> de.common <- which(dt[, 1] != 0 & dt[, 2] != 0)
> length(de.common)
[1] 2409
> head(tfit$genes$SYMBOL[de.common], n = 20)
 [1] "Xkr4"     "Rgs20"    "Cpa6"     "Sulf1"    "Eya1"     "Msc"      "Sbspon"   "Pi15"     "Crispld1"
[10] "Kcnq5"    "Ptpn18"   "Arhgef4"  "Cracdl"   "Aff3"     "Npas2"    "Tbc1d8"   "Creg2"    "Il1r1"   
[19] "Il18r1"   "Il18rap" 
> vennDiagram(dt[, 1:2], circle.col = c("turquoise", "salmon"))
> write.fit(tfit, dt, file = "results.txt")
> basal.vs.lp <- topTreat(tfit, coef = 1, n = Inf)
> basal.vs.ml <- topTreat(tfit, coef = 2, n = Inf)
> head(basal.vs.lp)
       ENTREZID SYMBOL TXCHROM     logFC  AveExpr         t      P.Value    adj.P.Val
12759     12759    Clu   chr14 -5.442877 8.857907 -33.44429 3.990899e-10 2.703871e-06
53624     53624  Cldn7   chr11 -5.514605 6.296762 -32.94533 4.503694e-10 2.703871e-06
242505   242505  Rasef    chr4 -5.921741 5.119585 -31.77625 6.063249e-10 2.703871e-06
67451     67451   Pkp2   chr16 -5.724823 4.420495 -30.65370 8.010456e-10 2.703871e-06
228543   228543   Rhov    chr2 -6.253427 5.486640 -29.46244 1.112729e-09 2.703871e-06
70350     70350  Basp1   chr15 -6.073297 5.248349 -28.64890 1.380545e-09 2.703871e-06
> head(basal.vs.ml)
       ENTREZID  SYMBOL TXCHROM     logFC  AveExpr         t      P.Value    adj.P.Val
242505   242505   Rasef    chr4 -6.510470 5.119585 -35.49093 2.573575e-10 1.915485e-06
53624     53624   Cldn7   chr11 -5.469160 6.296762 -32.52520 4.978446e-10 1.915485e-06
12521     12521    Cd82    chr2 -4.667737 7.070963 -31.82187 5.796191e-10 1.915485e-06
71740     71740 Nectin4    chr1 -5.556046 5.166292 -31.29987 6.760578e-10 1.915485e-06
20661     20661   Sort1    chr3 -4.908119 6.705784 -31.23083 6.761331e-10 1.915485e-06
15375     15375   Foxa1   chr12 -5.753884 5.625064 -28.34612 1.487280e-09 2.281914e-06
> plotMD(tfit, column = 1, status = dt[, 1], main = colnames(tfit)[1],
+        xlim = c(-8, 13))
> glMDPlot(tfit, coef = 1, status = dt, main = colnames(tfit)[1],
+          side.main = "ENTREZID", counts = x$counts, groups = group,
+          launch = TRUE)
> library("gplots")

Attaching package: ‘gplots’

The following object is masked from ‘package:IRanges’:

    space

The following object is masked from ‘package:S4Vectors’:

    space

The following object is masked from ‘package:stats’:

    lowess

> basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
> i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
> mycol <- colorpanel(1000, "blue", "white", "red")
> heatmap.2(v$E[i, ], scale = "row",
+           labRow = v$genes$SYMBOL[i], labCol = group,
+           col = mycol, trace = "none", density.info = "none")
> load(system.file("extdata", "mouse_c2_v5p1.rda", package = "RNAseq123"))
> class(Mm.c2)
[1] "list"
> Mm.c2$KEGG_GLYCOLYSIS_GLUCONEOGENESIS
  [1] "100042025" "100042427" "103988"    "106557"    "109785"    "110695"    "11522"     "11529"    
  [9] "11532"     "11668"     "11669"     "11670"     "11671"     "11674"     "11676"     "12039"    
 [17] "12040"     "12183"     "13382"     "13806"     "13807"     "13808"     "14120"     "14121"    
 [25] "14377"     "14378"     "14433"     "14447"     "14751"     "14782"     "15275"     "15277"    
 [33] "15452"     "16497"     "16498"     "16499"     "16828"     "16832"     "16833"     "18534"    
 [41] "18597"     "18598"     "18641"     "18642"     "18648"     "18655"     "18663"     "18746"    
 [49] "18770"     "19378"     "20322"     "212032"    "215974"    "216019"    "21881"     "21991"    
 [57] "226041"    "230163"    "232223"    "232491"    "235339"    "26358"     "26462"     "26876"    
 [65] "26926"     "277333"    "319625"    "353204"    "380660"    "433182"    "50493"     "56012"    
 [73] "56421"     "56752"     "56847"     "58810"     "60525"     "621603"    "622339"    "638833"   
 [81] "640374"    "666258"    "666488"    "66681"     "667451"    "668435"    "669429"    "67689"    
 [89] "68263"     "68738"     "69117"     "70974"     "72157"     "72168"     "72535"     "73458"    
 [97] "74419"     "74551"     "78894"     "79459"     "83553"    
> idx <- ids2indices(Mm.c2, id = rownames(v))
> cam.BasalvsLP <- camera(v, idx, design, contrast = contr.matrix[, 1])
> head(cam.BasalvsLP, 5)
                                            NGenes Direction       PValue          FDR
LIM_MAMMARY_STEM_CELL_UP                       739        Up 1.134757e-18 5.360590e-15
LIM_MAMMARY_STEM_CELL_DN                       630      Down 1.569957e-15 3.708238e-12
ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER    163        Up 1.437987e-13 2.264351e-10
SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP         183        Up 2.181862e-13 2.576779e-10
LIM_MAMMARY_LUMINAL_PROGENITOR_UP               87      Down 6.734613e-13 6.362863e-10
> cam.BasalvsML <- camera(v, idx, design, contrast = contr.matrix[, 2])
> head(cam.BasalvsML, 5)
                                            NGenes Direction       PValue          FDR
LIM_MAMMARY_STEM_CELL_UP                       739        Up 5.090937e-23 2.404959e-19
LIM_MAMMARY_STEM_CELL_DN                       630      Down 5.132446e-19 1.212284e-15
LIM_MAMMARY_LUMINAL_MATURE_DN                  166        Up 8.875174e-16 1.397544e-12
LIM_MAMMARY_LUMINAL_MATURE_UP                  180      Down 6.287301e-13 7.425303e-10
ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER    163        Up 1.684323e-12 1.591348e-09
> cam.LPvsML <- camera(v, idx, design, contrast = contr.matrix[, 3])
> head(cam.LPvsML, 5)
                                        NGenes Direction       PValue          FDR
LIM_MAMMARY_LUMINAL_MATURE_UP              180      Down 8.497295e-14 3.401020e-10
LIM_MAMMARY_LUMINAL_MATURE_DN              166        Up 1.439890e-13 3.401020e-10
LIM_MAMMARY_LUMINAL_PROGENITOR_UP           87        Up 3.840915e-11 6.048160e-08
REACTOME_RESPIRATORY_ELECTRON_TRANSPORT     91      Down 2.655349e-08 3.135967e-05
NABA_CORE_MATRISOME                        222        Up 4.430361e-08 4.185805e-05
> barcodeplot(efit$t[, 3], index = idx$LIM_MAMMARY_LUMINAL_MATURE_UP,
+             index2 = idx$LIM_MAMMARY_LUMINAL_MATURE_DN,
+             main = "LPvsML")
> sessionInfo()
R version 4.2.3 (2023-03-15 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22621)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gplots_3.1.3                              RColorBrewer_1.1-3                       
 [3] Mus.musculus_1.3.1                        TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
 [5] org.Mm.eg.db_3.16.0                       GO.db_3.16.0                             
 [7] OrganismDbi_1.40.0                        GenomicFeatures_1.50.4                   
 [9] GenomicRanges_1.50.2                      GenomeInfoDb_1.34.9                      
[11] AnnotationDbi_1.60.2                      IRanges_2.32.0                           
[13] S4Vectors_0.36.2                          Biobase_2.58.0                           
[15] BiocGenerics_0.44.0                       edgeR_3.40.2                             
[17] Glimma_2.8.0                              limma_3.54.2                             

loaded via a namespace (and not attached):
 [1] bitops_1.0-7                matrixStats_0.63.0          bit64_4.0.5                
 [4] filelock_1.0.2              progress_1.2.2              httr_1.4.5                 
 [7] tools_4.2.3                 utf8_1.2.3                  R6_2.5.1                   
[10] KernSmooth_2.23-20          DBI_1.1.3                   colorspace_2.1-0           
[13] tidyselect_1.2.0            prettyunits_1.1.1           DESeq2_1.38.3              
[16] bit_4.0.5                   curl_5.0.0                  compiler_4.2.3             
[19] graph_1.76.0                cli_3.6.1                   xml2_1.3.4                 
[22] DelayedArray_0.24.0         rtracklayer_1.58.0          caTools_1.18.2             
[25] scales_1.2.1                RBGL_1.74.0                 rappdirs_0.3.3             
[28] stringr_1.5.0               digest_0.6.31               Rsamtools_2.14.0           
[31] R.utils_2.12.2              XVector_0.38.0              pkgconfig_2.0.3            
[34] htmltools_0.5.5             MatrixGenerics_1.10.0       dbplyr_2.3.2               
[37] fastmap_1.1.1               htmlwidgets_1.6.2           rlang_1.1.1                
[40] rstudioapi_0.14             RSQLite_2.3.1               BiocIO_1.8.0               
[43] generics_0.1.3              jsonlite_1.8.4              gtools_3.9.4               
[46] BiocParallel_1.32.6         dplyr_1.1.2                 R.oo_1.25.0                
[49] RCurl_1.98-1.12             magrittr_2.0.3              GenomeInfoDbData_1.2.9     
[52] Matrix_1.5-3                Rcpp_1.0.10                 munsell_0.5.0              
[55] fansi_1.0.4                 lifecycle_1.0.3             R.methodsS3_1.8.2          
[58] stringi_1.7.12              yaml_2.3.7                  SummarizedExperiment_1.28.0
[61] zlibbioc_1.44.0             BiocFileCache_2.6.1         grid_4.2.3                 
[64] blob_1.2.4                  parallel_4.2.3              crayon_1.5.2               
[67] lattice_0.20-45             Biostrings_2.66.0           annotate_1.76.0            
[70] hms_1.1.3                   KEGGREST_1.38.0             locfit_1.5-9.7             
[73] pillar_1.9.0                rjson_0.2.21                geneplotter_1.76.0         
[76] codetools_0.2-19            biomaRt_2.54.1              XML_3.99-0.14              
[79] glue_1.6.2                  BiocManager_1.30.20         png_0.1-8                  
[82] vctrs_0.6.2                 gtable_0.3.3                cachem_1.0.8               
[85] ggplot2_3.4.2               xtable_1.8-4                restfulr_0.0.15            
[88] tibble_3.2.1                GenomicAlignments_1.34.1    memoise_2.0.1              
> 
> sessionInfo()
R version 4.2.3 (2023-03-15 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22621)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gplots_3.1.3                              RColorBrewer_1.1-3                       
 [3] Mus.musculus_1.3.1                        TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
 [5] org.Mm.eg.db_3.16.0                       GO.db_3.16.0                             
 [7] OrganismDbi_1.40.0                        GenomicFeatures_1.50.4                   
 [9] GenomicRanges_1.50.2                      GenomeInfoDb_1.34.9                      
[11] AnnotationDbi_1.60.2                      IRanges_2.32.0                           
[13] S4Vectors_0.36.2                          Biobase_2.58.0                           
[15] BiocGenerics_0.44.0                       edgeR_3.40.2                             
[17] Glimma_2.8.0                              limma_3.54.2                             

loaded via a namespace (and not attached):
 [1] bitops_1.0-7                matrixStats_0.63.0          bit64_4.0.5                
 [4] filelock_1.0.2              progress_1.2.2              httr_1.4.5                 
 [7] tools_4.2.3                 utf8_1.2.3                  R6_2.5.1                   
[10] KernSmooth_2.23-20          DBI_1.1.3                   colorspace_2.1-0           
[13] tidyselect_1.2.0            prettyunits_1.1.1           DESeq2_1.38.3              
[16] bit_4.0.5                   curl_5.0.0                  compiler_4.2.3             
[19] graph_1.76.0                cli_3.6.1                   xml2_1.3.4                 
[22] DelayedArray_0.24.0         rtracklayer_1.58.0          caTools_1.18.2             
[25] scales_1.2.1                RBGL_1.74.0                 rappdirs_0.3.3             
[28] stringr_1.5.0               digest_0.6.31               Rsamtools_2.14.0           
[31] R.utils_2.12.2              XVector_0.38.0              pkgconfig_2.0.3            
[34] htmltools_0.5.5             MatrixGenerics_1.10.0       dbplyr_2.3.2               
[37] fastmap_1.1.1               htmlwidgets_1.6.2           rlang_1.1.1                
[40] rstudioapi_0.14             RSQLite_2.3.1               BiocIO_1.8.0               
[43] generics_0.1.3              jsonlite_1.8.4              gtools_3.9.4               
[46] BiocParallel_1.32.6         dplyr_1.1.2                 R.oo_1.25.0                
[49] RCurl_1.98-1.12             magrittr_2.0.3              GenomeInfoDbData_1.2.9     
[52] Matrix_1.5-3                Rcpp_1.0.10                 munsell_0.5.0              
[55] fansi_1.0.4                 lifecycle_1.0.3             R.methodsS3_1.8.2          
[58] stringi_1.7.12              yaml_2.3.7                  SummarizedExperiment_1.28.0
[61] zlibbioc_1.44.0             BiocFileCache_2.6.1         grid_4.2.3                 
[64] blob_1.2.4                  parallel_4.2.3              crayon_1.5.2               
[67] lattice_0.20-45             Biostrings_2.66.0           annotate_1.76.0            
[70] hms_1.1.3                   KEGGREST_1.38.0             locfit_1.5-9.7             
[73] pillar_1.9.0                rjson_0.2.21                geneplotter_1.76.0         
[76] codetools_0.2-19            biomaRt_2.54.1              XML_3.99-0.14              
[79] glue_1.6.2                  BiocManager_1.30.20         png_0.1-8                  
[82] vctrs_0.6.2                 gtable_0.3.3                cachem_1.0.8               
[85] ggplot2_3.4.2               xtable_1.8-4                restfulr_0.0.15            
[88] tibble_3.2.1                GenomicAlignments_1.34.1    memoise_2.0.1              
> 
> 
> # > sessionInfo()
> # > sessionInfo()
> # R version 3.6.3 (2020-02-29)
> if (!requireNamespace("BiocManager", quietly = TRUE))
