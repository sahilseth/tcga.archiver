## compileDCC.R
## Author: Sahil Seth
## IACS, MD Anderson
## Contact: sseth@mdanderson.org
## Created: 12/07/2011
## Updated: 09/06/2012
## Title: This contains the main functions used to compile a archive for DCC submission.

## an os compatible command
get_md5 <- function(file){
  switch(Sys.info()[['sysname']],
         #Windows= {print("I'm a Windows PC.")},
         Linux  = {
           md5sum = strsplit(system(paste("md5sum ",file,sep=""),intern=TRUE)," ")[[1]][1];
           #print("I'm a penguin.")
           },
         Darwin = {
           md5sum = strsplit(system(sprintf("md5 %s", file), intern=TRUE), " = ")[[1]][2]
           #print("I'm a Mac.")
         }
  )
  return(md5sum)
}

get_sdrf_row <- function(vec, dat_level, dat_batch, dat_rev, disease, opt_data, file_type, center = 'hms.harvard.edu', platform){
  ## file_type <- tail(strsplit(file, "\\.")[[1]],1) #get file extension
  rowData <- opt_data[,file_type]
  names(rowData) <- rownames(opt_data)
  ## for each FILE:
  file <- unlist(vec["files"])
  cat("Working on ",file,"\n")
  rowData["curFile"] <- basename(unlist(vec["files"])) #the new name of the file, still to be copied over
  rowData["curArchName"] <- sprintf("%s_%s.%s.%s.%s.%s.0",center,disease,platform,dat_level,dat_batch,dat_rev)
  ## rowData["curFileMD5"] <- strsplit(grep(orgFile,md5Sums,value=TRUE)," ")[[1]][1] #to get the MD5
  rowData["curFileMD5"] <- get_md5(file)
  ## tumors
  tumData <- rowData
  ## tumData["barcode"] <- strsplit(basename(file),"_")[[1]][1]
  tumData["barcode"] <- vec["TUMORSAMPLEID"]
  ## tumData["extName"] <- getTCGAUUID(tumData["barcode"])
  tumData["extName"] <- vec["TUMORSAMPLEUUID"]
  tumData["labNm"] <- paste(tumData["barcode"],rowData["labNm"],sep="_")
  tumData["libPrepNm"] <- paste(tumData["barcode"],rowData["libPrepNm"],sep="_")
  tumData["clustGenNm"] <- paste(tumData["barcode"],rowData["clustGenNm"],sep="_")
  tumData["dnaSeqNm"] <- paste(tumData["barcode"],rowData["dnaSeqNm"],sep="_")
  tumData["seqAlnNm"] <- paste(tumData["barcode"],rowData["seqAlnNm"],sep="_")
  tumData["dnaSeqNm"] <- paste(tumData["barcode"],rowData["dnaSeqNm"],sep="_")
  #tmp <- getCGHubStatus(vec["TUMORANALYSISUUID"],by="analysis_id",get="all")
  tumData["seqAlnFile"] <- basename(as.character(vec["TUMORSAMPLEBAM"]))#tmp$filename
  tumData["seqAlnUUID"] <- vec["NORMALANALYSISUUID"]#tmp$uuid
  tumData["curDatTransName"] <- paste(tumData["barcode"],rowData["curDatTransName"],sep="_")
  ## cat(tmp$uuid,"\n")
  ## normals
  normData <- rowData
  ## normData["barcode"] <- strsplit(basename(file),"_")[[1]][2]
  normData["barcode"] <- vec["NORMALSAMPLEID"]
  ## normData["extName"] <- getTCGAUUID(normData["barcode"])
  normData["extName"] <- vec["NORMALSAMPLEUUID"]
  normData["labNm"] <- paste(normData["barcode"],rowData["labNm"],sep="_")
  normData["libPrepNm"] <- paste(normData["barcode"],rowData["libPrepNm"],sep="_")
  normData["clustGenNm"] <- paste(normData["barcode"],rowData["clustGenNm"],sep="_")
  normData["dnaSeqNm"] <- paste(normData["barcode"],rowData["dnaSeqNm"],sep="_")
  normData["seqAlnNm"] <- paste(normData["barcode"],rowData["seqAlnNm"],sep="_")
  normData["dnaSeqNm"] <- paste(normData["barcode"],rowData["dnaSeqNm"],sep="_")
  #tmp <- getCGHubStatus(vec["NORMALANALYSISUUID"],by="analysis_id",get="all")
  normData["seqAlnFile"] <- basename(as.character(vec["NORMALSAMPLEBAM"]))#tmp$filename
  normData["seqAlnUUID"] <- vec["NORMALANALYSISUUID"]#tmp$uuid
  normData["curDatTransName"] <- paste(normData["barcode"],rowData["curDatTransName"],sep="_")
  ## cat(tmp$uuid,"\n")
  return(c(tumData,normData))
}


createArchive <- function(mat=mat,disease,center="hms.harvard.edu",platform="IlluminaHiSeq_DNASeqC",file_type,type,
                          batchRev,dat_level=dat_level,dat_batch=dat_batch,dat_rev=dat_rev,wrkPath=wrkPath, opt_data){
  ## get old sdrf
  sdrf <- vector()
  rm <- mat[,"todo"]=="rm";matrm <- mat[rm,];mat <- mat[!rm,]
  files <- try(as.character(mat$files))
  if(batchRev$mage.rev!=0){## is this is not the first revision, update archive collumn, with latest revision (this one)
    inDir <- file.path(gccPath,"DCC/submissions",sprintf("%s_%s.%s.mage-tab.%s.%s.0",
                                                         center,disease,platform,batchRev$mage.batch,
                                                         batchRev$mage.rev-1))
    sdrf <- read.delim(file.path(inDir,sprintf("%s_%s.%s.sdrf.txt",center,disease,platform)),
                       sep = "\t", header = TRUE, as.is = TRUE)
    ## GET THE NEW ARCHIVE NAME, WHERE IT BELONGS TO A OLDER REVISION (OF THIS BATCH)
    if(TRUE){
      sdrf[,"Comment..TCGA.Archive.Name."] <- gsub(paste(dat_level,dat_batch,dat_rev-2,sep="."),
                                                   paste(dat_level,dat_batch,dat_rev-1,sep="."),
                                                   sdrf[,"Comment..TCGA.Archive.Name."])
    }
    ## REMOVE FILES FROM SDRF, THESE ARE THE FILES TO BE REMOVED
    sdrf <- sdrf[!sdrf[,"Derived.Data.File"] %in% matrm$newFiles,]
    ## remove VCFs from this SDRF:, stop doing this now !!
    ## sdrf <- sdrf[!sdrf[,"Comment..TCGA.File.Type."] =="vcf",]
  }
  ### GET NEW SDRF ROWS data files, for each of the rows in mat get SDRF table
  files <- mat$files;files2 <- mat$newFiles
  out <- mclapply(1:length(files), function(i){
    dcc.getSDRFTable(vec=mat[i,],dat_level=dat_level,dat_batch=dat_batch,dat_rev=dat_rev,disease=disease,opt_data=opt_data,
                     file_type=file_type)
  },mc.cores=1)
  tab <- matrix(unlist(out),byrow=TRUE,ncol=dim(opt_data)[1])
  ######### Working on Metadata ###############
  ## read old records for the disease and do rbind, then make respective folder and write SDRF
  inDir <- file.path(gccPath,"DCC/data",sprintf("%s_%s.%s.mage-tab",
                                                center,disease,platform)) #to get idf and description
  outDir <- file.path(wrkPath,sprintf("%s_%s.%s.mage-tab.%s.%s.0",
                                      center,disease,platform,batchRev$mage.batch, batchRev$mage.rev))
  if(file.exists(outDir)){file.remove(list.files(outDir,full.names=TRUE))}else{dir.create(outDir)}
  file.copy(list.files(inDir,full.name=TRUE), outDir) #copy description and idf from data folder
  ##############################################################################################################
  if(length(sdrf)<1){
    sdrf <- rbind(opt_data[,1],tab) # new stuff; tab has the new SDRF rows
  }else{
    sdrf <- rbind(opt_data[,1],as.matrix(sdrf),tab) # after things have cooled down
  }
  write.table(sdrf,file.path(outDir,sprintf("%s_%s.%s.sdrf.txt",center,disease,platform)),
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  ##Get new manifest
  setwd(outDir)
  system(paste("md5sum ", paste(list.files("."),collapse=" "),">> MANIFEST.txt"))
  ## remove rows from manifest
  setwd(dirname(outDir)); system(paste("tar -czvf ", outDir,".tar.gz ", basename(outDir),sep=""))
  system(paste("md5sum ", basename(outDir),".tar.gz > ",outDir,".tar.gz.md5",sep=""))
  ############ working on data ##################
  ##copy files to dataDir
  outDir <- file.path(wrkPath,sprintf("%s_%s.%s.%s.%s.%s.0",
                                      center,disease,platform,dat_level,dat_batch,dat_rev))
  ## incase of new archive, in and out will be the same
  inDir <- file.path(gsub("working","submissions",wrkPath),sprintf("%s_%s.%s.%s.%s.%s.0",
                                                                   center,disease,platform,dat_level,dat_batch,ifelse(dat_rev-1<0,0,dat_rev-1)))
  if(file.exists(outDir)){file.remove(list.files(outDir,full.names=TRUE))}else{dir.create(outDir)}
  system(sprintf("cp %s/MANIFEST.txt %s/MANIFEST.txt",inDir,outDir))
  file.copy(files,file.path(outDir,files2),overwrite=TRUE) #copy the files
  ## this may have duplicates
  setwd(outDir);system(paste("md5sum ", paste(list.files(".",pattern=".vcf|tsv"),collapse=" "),">> MANIFEST.txt"))
  ## remove the files from MANIFEST (those now unapproved)
  if(dim(matrm)[1]>0)
    system(sprintf("cat MANIFEST.txt | grep -v '%s' > MANIFEST.txt",paste(matrm$files,"|",sep="")))
  ## remove duplicates
  manifest <- read.table(sprintf("%s/MANIFEST.txt",outDir),sep=" ")[,c(1,3)]
  ## making tar gz file
  setwd(dirname(outDir))
  system(paste("tar -czvf ", outDir,".tar.gz ",basename(outDir),"/*",sep=""))
  system(paste("md5sum ", basename(outDir),".tar.gz > ",outDir,".tar.gz.md5",sep=""))
  ## check with QC Live
}


getCGHubStatus <- function(id, by="analysis_id", get="state"){ #get state, analysis_id, both
  require(XML);require(RCurl)
  url <- paste("https://cghub.ucsc.edu//cghub/metadata/analysisFull?",by,"=*",id,"*",sep="")
  ## xml <- system(sprintf("curl -s %s",url),intern=TRUE)
  xml <- system(sprintf("curl -s %s",url),intern=TRUE)
  state <- xmlValue(xmlRoot(xmlTreeParse(xml))[["Result"]][["state"]])
  uuid <- xmlValue(xmlRoot(xmlTreeParse(xml))[["Result"]][["analysis_id"]])
  filename <- try(xmlValue(xmlRoot(xmlTreeParse(xml))[["Result"]][["files"]][["file"]][["filename"]]))
  aliquotID <- try(xmlValue(xmlRoot(xmlTreeParse(xml))[["Result"]][["aliquot_id"]]))
  reason <- try(xmlValue(xmlRoot(xmlTreeParse(xml))[["Result"]][["reason"]]))
  if(get=="uuid"){
    ret <- uuid
  }else if(get=="state"){
    ret <- state
  }else if(get=="fileName"){
    ret <- filename
  }else if(get=="all" | get=="both")
    ret <- list(state=state,uuid=uuid,filename=filename,reason=reason,aliquotID=aliquotID)
  return(ret)
}
