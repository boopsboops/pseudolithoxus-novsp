#clear memory
rm(list = ls())
#read in morpho table
tab1 <- read.table("/home/rupert/LaTeX/nhamunda-checklist/pseudolithoxus_morphometrics_table.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names="Measurement")#?read.table , 
#remove landmark number, keep abbreviations and convert to Matrix
tab <- subset(tab1, select=-Landmarks)

##check data is read in correctly
#class(tab)
#dim(tab)
#head(tab)
#tab$CTGA_14486
#tab["Standard_length",]
#rownames(tab)
#colnames(tab)


#subset body measurements, remove percents, and convert to matrix
body <- na.omit(tab[tab$Percents=="body",])
body2 <- as.matrix(subset(body, select=-Percents))
#selects the correct row with standard lengths
sl <- as.matrix(subset(tab["Standard_length",], select=-Percents))

#perform a 'sweep' over the columns to get the % of SL for each
slp <- sweep(body2, 2, sl, FUN="/") *100#?sweep

#plot
#boxplot(slp, use.cols = FALSE, cex.axis=0.6, las=2)#?boxplot?plot

#for head lengths
head <- na.omit(tab[tab$Percents=="head",])
head2 <- as.matrix(subset(head, select=-Percents))
hl <- as.matrix(subset(tab["Head_length",], select=-Percents))
hlp <- sweep(head2, 2, hl, FUN="/") *100
#boxplot(hlp, use.cols = FALSE, cex.axis=0.6, las=2)#?boxplot?plot

#combined boxplot
com <- rbind(slp, hlp)
#boxplot(com, use.cols = FALSE, cex.axis=0.6, las=2)#?boxplot?plot

#boxplot ordered by range
ran <- as.vector(apply(com, 1, sd))
ex <- as.data.frame(cbind(com,ran))
ex <- ex[order(ex$"ran", decreasing=TRUE),]
ex <- as.matrix(subset(ex, select=-ran))
par(mar=c(10,5,2,2))
boxplot(ex, use.cols = FALSE, cex.axis=0.6, las=2, ylab="percent of SL/HL sorted by standard devation")


#libraries to melt and use new boxplot
library(car)
library(reshape2)
#form data into 3 cols
mex <- melt(ex)#?melt
par(mar=c(10,5,2,2))
#boxplot with outliers annotated
Boxplot(mex$value~mex$Var1, data=mex, labels=as.character(mex$Var2), cex.axis=0.6, las=2)#?Boxplot #?par

#check the individual values
ex["Abdominal_length",]

##get some summaries
#drop percents
tab <- subset(tab, select=-Percents)
#apply summary over the rows
apply(tab, 1, summary)
apply(tab, 1, Mode)
apply(tab, 1, median)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


## create the final table for paper and reorganise the data into correct order

#get SL for all 
sl <- as.matrix(subset(tab["Standard_length",], select=-Percents))
#bind columns for all measurements
com <- rbind(sl, slp, hlp)

#run the summary functions over the rows
me <- apply(com, 1, mean)
sv <- apply(com, 1, sd)
mn <- apply(com, 1, min)
mx <- apply(com, 1, max)

#subset for the holotype only
slh <- sl[,"PN11"]
names(slh) <- rownames(sl)
holo <- c(slh, slp[,"PN11"], hlp[,"PN11"])

#make a vector for the percentage in column
pof <- c(NA, rep("SL", dim(slp)[1]), rep("HL", dim(hlp)[1]))

#create data frame of all 
fintab <- cbind(cbind(as.vector(na.omit(tab1$"Landmarks")), rownames(com)), pof, round(data.frame(cbind(as.numeric(holo), as.numeric(me), as.numeric(sv), as.numeric(mn), as.numeric(mx))), digits=1))
#label the column names
colnames(fintab) <- c("Landmarks", "Measurement", "Percentage in", "Holotype", "Mean", "SD", "Min.", "Max.")
#tidy up the entries
fintab$Landmarks <- as.factor(gsub("--", "‒", fintab$Landmarks))
fintab$Measurement <- gsub("--", "-", fintab$Measurement)
fintab$Measurement <- as.factor(gsub("_", " ", fintab$Measurement))
#check
print(fintab)

#save
write.table(fintab, file = "/media/1TB/auto_backed_up/LaTeX/nhamunda-checklist/final_morpho_table.csv", quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
