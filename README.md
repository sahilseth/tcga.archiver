tcga.archiver
=================

Creating metadata for TCGA lowpass data

### Installing this package
```
require(devtools)
install_github(repo="tcga.archiver", username="sahilseth")

```

### Reading in example params
```
require(tcga.archiver)
dcc_params = file.path(path.package("tcga.archiver"), "files/DCC.params")
opt_data <- read.delim(dcc_params, sep = "\t", header = TRUE, as.is = TRUE)
```

```
file_type = "tsv" ## or vcf for breakdancer
disease = 'HNSC'
platform = "IlluminaHiSeq_DNASeqC"
center = "hms.harvard.edu"


## these are the TSV files we upload:
tsvs = file.path(path.package("tcga.archiver"), 
  "files/TCGA-CN-5356-01A-01D-1431___120412_SN208_0283_AD0V1GACXX---TCGA-CN-5356-10A-01D-1431___120412_SN208_0284_BD0T51ACXX.tsv")

## and example sample
vec <- c(files = tsvs, 
  NORMALSAMPLEID = "TCGA-CN-5356-10A-01D-1431-02", TUMORSAMPLEID = "TCGA-CN-5356-01A-01D-1431-02", 
  NORMALANALYSISUUID = getCGHubStatus("TCGA-CN-5356-01A-01D-1431-02", 
    by = "legacy_sample_id", get = "both")$uuid,
  TUMORANALYSISUUID = getCGHubStatus("TCGA-CN-5356-10A-01D-1431-02", 
    by = "legacy_sample_id", get = "both")$uuid)

out <- get_sdrf_row(vec = vec, dat_level = 3 ,dat_batch = 0 ,dat_rev = 0 ,disease = disease, opt_data = opt_data,
  file_type = file_type, platform = platform, center = center)

```