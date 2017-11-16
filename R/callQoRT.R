#' @ Use QoRTs to get the count of subfeatures in each gene
#'
#' @ cmd1
#' @ inputfile:bam file, gene annotation file,
#' @ outfile
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
#' dir.name="/media/H_driver/Aimin_project/GOSJ_STAR_Bam/"
#'
#' file.name=dir(dir.name,recursive = TRUE,pattern="sorted.bam")
#' file.name.whole<-paste0(dir.name,file.name)
#' file.name.selected<-file.name.whole
#' file.name.selected.2<-as.list(file.name.selected)
#' names(file.name.selected.2)=sapply(strsplit(file.name.selected,split="\\/"),"[[",6)
#''
#'"
#' file.name.selected.3<-file.name.selected.2[-c(1,3)]
#'
#' cmd1="java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --noGzipOutput --keepMultiMapped"
#'
#' re.out<-lapply(file.name.selected.2,callQoRT,gtf_file=gtf1,runing_cmd=cmd1)
#'
#'
#' cmd2="java -Xmx5000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --noGzipOutput"
#'
#' gtf1="/media/H_driver/Aimin_project/GTF_Files/Homo_sapiens.GRCh38.84.processed.sorted.2.gtf"
#'
#' gtf1="/media/H_driver/Aimin_project/GTF_Files/Homo_sapiens.GRCh38.84.processed.sorted.2.chr4.only.gtf"
#'
#'
#' re.out<-lapply(file.name.selected.3,callQoRT,gtf_file=gtf1,runing_cmd=cmd2)
#'
#'
#'cmd3="java -Xmx5000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --noGzipOutput --testRun"
#'re.out<-lapply(file.name.selected.2,callQoRT,output.dir="test_chr4",gtf_file=gtf1,runing_cmd=cmd3)
#'
#'cmd33="java -Xmx5000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --noGzipOutput"
#'re.out<-lapply(file.name.selected.2,callQoRT,output.dir="chr4_only",gtf_file=gtf1,runing_cmd=cmd33)
#'
#'
#'data.chrom<-read.table("/media/H_driver/Aimin_project/SRR1660308_STAR_out.sorted.bam/QC.chromCount.txt",header=TRUE)
#'
#' chro=c(as.character(unique(data.chrom$CHROM)))
#'
#' class(chro)
#' chro.no.4<-chro[-which(chro=="chr4")]
#'
#' chro.list<-paste(chro.no.4, collapse = ',')
#'
#'
#' cmd4="java -Xmx5000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --noGzipOutput --dropChrom chr1,chr10,chr11,chr11_KI270721v1_random,chr12,chr13,chr14,chr14_GL000009v2_random,chr14_GL000194v1_random,chr14_GL000225v1_random,chr14_KI270722v1_random,chr14_KI270723v1_random,chr14_KI270726v1_random,chr15,chr15_KI270727v1_random,chr16,chr16_KI270728v1_random,chr17,chr17_GL000205v2_random,chr17_KI270729v1_random,chr18,chr19,chr1_KI270706v1_random,chr1_KI270707v1_random,chr1_KI270709v1_random,chr1_KI270710v1_random,chr1_KI270711v1_random,chr1_KI270712v1_random,chr1_KI270713v1_random,chr1_KI270714v1_random,chr2,chr20,chr21,chr22,chr22_KI270731v1_random,chr22_KI270732v1_random,chr22_KI270733v1_random,chr22_KI270734v1_random,chr22_KI270735v1_random,chr22_KI270736v1_random,chr22_KI270738v1_random,chr2_KI270716v1_random,chr3,chr3_GL000221v1_random,chr4_GL000008v2_random,chr5,chr6,chr7,chr8,chr9,chr9_KI270718v1_random,chr9_KI270719v1_random,chr9_KI270720v1_random,chrM,chrUn_GL000195v1,chrUn_GL000213v1,chrUn_GL000214v1,chrUn_GL000216v2,chrUn_GL000218v1,chrUn_GL000219v1,chrUn_GL000220v1,chrUn_GL000224v1,chrUn_KI270311v1,chrUn_KI270315v1,chrUn_KI270330v1,chrUn_KI270337v1,chrUn_KI270362v1,chrUn_KI270435v1,chrUn_KI270438v1,chrUn_KI270442v1,chrUn_KI270467v1,chrUn_KI270511v1,chrUn_KI270519v1,chrUn_KI270522v1,chrUn_KI270590v1,chrUn_KI270741v1,chrUn_KI270742v1,chrUn_KI270743v1,chrUn_KI270744v1,chrUn_KI270745v1,chrUn_KI270746v1,chrUn_KI270747v1,chrUn_KI270748v1,chrUn_KI270750v1,chrUn_KI270751v1,chrUn_KI270754v1,chrX,chrY"
#' re.out<-lapply(file.name.selected.2,callQoRT,output.dir="chr4_drop_other",gtf_file=gtf1,runing_cmd=cmd4)
#'
#' gtf.chr22="/media/H_driver/Aimin_project/GTF_Files/Homo_sapiens.GRCh38.84.processed.sorted.2.chr22.only.gtf"
#'
#' cmd5="java -Xmx5000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --noGzipOutput --dropChrom chr1,chr10,chr11,chr11_KI270721v1_random,chr12,chr13,chr14,chr14_GL000009v2_random,chr14_GL000194v1_random,chr14_GL000225v1_random,chr14_KI270722v1_random,chr14_KI270723v1_random,chr14_KI270726v1_random,chr15,chr15_KI270727v1_random,chr16,chr16_KI270728v1_random,chr17,chr17_GL000205v2_random,chr17_KI270729v1_random,chr18,chr19,chr1_KI270706v1_random,chr1_KI270707v1_random,chr1_KI270709v1_random,chr1_KI270710v1_random,chr1_KI270711v1_random,chr1_KI270712v1_random,chr1_KI270713v1_random,chr1_KI270714v1_random,chr2,chr20,chr21,chr22_KI270731v1_random,chr22_KI270732v1_random,chr22_KI270733v1_random,chr22_KI270734v1_random,chr22_KI270735v1_random,chr22_KI270736v1_random,chr22_KI270738v1_random,chr2_KI270716v1_random,chr3,chr3_GL000221v1_random,chr4,chr4_GL000008v2_random,chr5,chr6,chr7,chr8,chr9,chr9_KI270718v1_random,chr9_KI270719v1_random,chr9_KI270720v1_random,chrM,chrUn_GL000195v1,chrUn_GL000213v1,chrUn_GL000214v1,chrUn_GL000216v2,chrUn_GL000218v1,chrUn_GL000219v1,chrUn_GL000220v1,chrUn_GL000224v1,chrUn_KI270311v1,chrUn_KI270315v1,chrUn_KI270330v1,chrUn_KI270337v1,chrUn_KI270362v1,chrUn_KI270435v1,chrUn_KI270438v1,chrUn_KI270442v1,chrUn_KI270467v1,chrUn_KI270511v1,chrUn_KI270519v1,chrUn_KI270522v1,chrUn_KI270590v1,chrUn_KI270741v1,chrUn_KI270742v1,chrUn_KI270743v1,chrUn_KI270744v1,chrUn_KI270745v1,chrUn_KI270746v1,chrUn_KI270747v1,chrUn_KI270748v1,chrUn_KI270750v1,chrUn_KI270751v1,chrUn_KI270754v1,chrX,chrY"
#'
#' re.out<-lapply(file.name.selected.2,callQoRT,output.dir="chr22_drop_other",gtf_file=gtf.chr22,runing_cmd=cmd5)
#'
#'
#'
#'
#'

callQoRT<-function(input_file,runing_cmd,output.dir=NULL,gtf_file){

  inputfile=paste(input_file,gtf_file,sep=" ")

  if(is.null(output.dir)){
  outfile=paste("",sapply(strsplit(input_file,split="\\/"),"[[",2),sapply(strsplit(input_file,split="\\/"),"[[",3),
                sapply(strsplit(input_file,split="\\/"),"[[",4),sapply(strsplit(input_file,split="\\/"),"[[",6),sep="/")
}else{
  outfile=paste("",sapply(strsplit(input_file,split="\\/"),"[[",2),sapply(strsplit(input_file,split="\\/"),"[[",3),
                sapply(strsplit(input_file,split="\\/"),"[[",4),paste0(sapply(strsplit(input_file,split="\\/"),"[[",6),"_",output.dir),sep="/")
  }


  cmd2=paste(runing_cmd,inputfile,outfile,sep=" ")

  print(cmd2)


  system(cmd2, intern = TRUE, ignore.stderr = TRUE)

  return(cmd2)

}

# input.file.dir <- "/media/dropbox/Aimin_project/DI/IRFinderResults/"
#
# input.file.dir <- "/media/pegasus/aiminy_project/DI"

# runing_cmd= "java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --noGzipOutput --keepMultiMapped"

# output.dir= "/media/pegasus/aiminy_project/DI_Gene_Based"

# gtf_file = "/media/H_driver/Aimin_project/GTF_Files/Homo_sapiens.GRCh38.84.processed.sorted.2.gtf"
#
# GOSJ:::callQoRT2(input.file.dir,runing_cmd,output.dir,gtf_file)

callQoRT2<-function(input.file.dir,runing_cmd,output.dir,gtf_file){

  input.file <- list.files(input.file.dir,pattern = "Unsorted.bam", all.files = T,full.names = T, recursive = T,ignore.case = FALSE, include.dirs = T, no.. = FALSE)

  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

  cmdL <- lapply(input.file,function(u,gtf_file,runing_cmd,output.dir) {


    outfile <- file.path(output.dir,basename(dirname(u)))

    cmd=paste(runing_cmd,u,gtf_file,outfile,sep=" ")

    cat(cmd,"\n\n")

    system(cmd, intern = TRUE, ignore.stderr = TRUE)


  },gtf_file,runing_cmd,output.dir)


}







