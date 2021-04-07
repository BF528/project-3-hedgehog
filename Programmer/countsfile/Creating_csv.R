library(purrr)
library(data.table)

### Creating a directory file which points to the directory of text file location ### 
directory<-"/projectnb2/bf528/users/hedgehog/project_3/project-3-hedgehog/Programmer/countsfile"

###Change working directory###
setwd("/projectnb2/bf528/users/hedgehog/project_3/project-3-hedgehog/Programmer/countsfile/counts_only")

### Merging all text files###

cc <-c("SRR1177966.txt", "SRR1177969.txt","SRR1177970.txt","SRR1177993.txt", "SRR1177994.txt","SRR1177995.txt","SRR1177998.txt", "SRR1178001.txt","SRR1178003.txt")
class(cc)

allfiles<-file.path(directory, cc = c("SRR1177966.txt", "SRR1177969.txt","SRR1177970.txt","SRR1177993.txt", "SRR1177994.txt","SRR1177995.txt","SRR1177998.txt", "SRR1178001.txt","SRR1178003.txt"))
class(allfiles)
as.vector(allfiles)
class(allfiles)
cc_all<-sapply(allfiles,FUN = read.csv)

cc_all <- map_df(allfiles, fread) 
names(cc_all)[1] <- "Geneid"




### Changing the names header of the counts file to corresponding file name ### 
names(cc_all)[7:15] <- c("SRR1177966","SRR1177969","SRR1177970","SRR1177993", "SRR1177994","SRR1177995","SRR1177998","SRR1178001", "SRR1178003")
write.csv(cc_all, "allcounts.csv",row.names = FALSE)
str(cc_all)

