# This is example r script for punctaanalysis tool kit.# If you use batch program, you can use this script.# To use it, start R and change the working directory to where your pictures and data are stored, # and just copy everything here and paste it on R terminal and push return key.# Once you have run this script, you will get some pdf density graphs,# and ready to test data sets. ex. wtfixwidth, daffixwidth...# Then you can test with something like this ## t.test(wtfixwidth, daffixwidth)# # I don't know what kind of test is best. Making a choice is your job.#
# This script doesn't do log transformation.
# So, parameters involving length, such as width, distance, are shown in micrometer.
## Also, here is an example of density graph# plot(density(wtwidth), col=2, xlim=c(-2, 2), ylim=c(0,1))# par(new =T)# plot(density(dafwidth), xlim=c(-2, 2), ylim=c(0,1))# 
# percenthist
# You can make histogram with percent y-axis by percenthist function.
# ex. percenthist(wtfixwidth)
# right y-axis indicating percent.
#
# Also, you can specify x, y axis, range etc.
# percenthist(wtfixwidth, xlim=c(0,3),ylim=c(0,40), br=seq(0,5,by=0.3),border=1,col=2)
# xlim; Restrict x-axis range depicted in the graph. 
# ylim; Restrict y-axis range depicted in the graph by percent. But, bottom is fixed 0.
# br; Specify breakpoints between histogram cells. See help of hist function for detail.
# seq; Make a vector. above case from 0 to 5 by 0.3 increment. 
# When br assign by this way, it must cover whole data set. 
# In other words, if there is a data with large value, such as 10, you must inclued it;
# like seq(0,10,by=0.3) even xlim is 0 to 5.
# 
# I hope this example file helps you to learn R scripting.  20090326 Taizo######## read datacount <- read.table("count.txt", header = TRUE);width <- read.table("width.txt", header = TRUE);gap <- read.table("gap.txt", header = TRUE);distance <- read.table("distance.txt", header = TRUE);intensity <- read.table("intensity.txt", header = TRUE);lineardensity <- read.table("lineardensity.txt", header =TRUE);fixwidth <- read.table("fixwidth.txt", header =TRUE);fixvolume <- read.table("fixvolume.txt", header =TRUE);fixgap <- read.table("fixgap.txt", header =TRUE);#fixwidthlineardensity <- read.table("fixwidthlineardensity", header =TRUE)#fixgaplineardensity <- read.table("fixgaplineardensity", header =TRUE)###########  spadework #list up data types what you want to analyseanalysisdata<- c("width", "fixwidth","gap","fixvolume", "intensity", "lineardensity", "distance", "fixgap")datanum <- length(analysisdata)#Detect genotypes from files automaticallygenotype.list <- levels(count$genotype);numberofgenotype <- length(genotype.list);# dataframe operation. New vector data are made for each genotype# follwing makes genotypewidh, genotypegaps etc. (ex: wtwidth, e5ewidth, fc16witdh,....wtvolume, e5volme.)# "eval(parse(text = ....)))" run the text as command.for (j in 1:datanum){for (i in 1:numberofgenotype){# This is doing like this. wtwidth <- width[width$genotype == "wt",]eval(parse(text = paste(genotype.list[i], analysisdata[j], "<-", analysisdata[j], "[", analysisdata[j], "$genotype == genotype.list[", i, "], ]", sep ="")))# This is doing like this. wtwidth <- wtwidth$widtheval(parse(text = paste(genotype.list[i], analysisdata[j], "<-", genotype.list[i], analysisdata[j], "$",  analysisdata[j], sep = "")))}}######## density graphslibrary(KernSmooth)#list up data types what you want to analyse with density graph hereanalysisdata<- c("width", "fixwidth","gap","fixvolume", "intensity", "distance", "fixgap")datanum <- length(analysisdata)# genotypes that you want to analyse. everything or can be designate##################### You may list up genotypes what you want to analyse here.#genotype.list <- c("wt", "e719","e5", "fc16", "e1598e5", "e1598fc16");numberofgenotype <- length(genotype.list);for (j in 1:datanum){xmin <- NULLxmax <- NULLymin <- NULLymax <- NULLfor (i in 1:numberofgenotype){# This is doing something like this "wt <- bkde(wtwidth)"eval(parse(text = paste(genotype.list[i], "<- bkde(", genotype.list[i], analysisdata[j], ")", sep="")))# appending min and max values to use as axis depiction latereval(parse(text = paste("xmin <- append(xmin, min(", genotype.list[i], "$x))", sep="")))eval(parse(text = paste("xmax <- append(xmax, max(", genotype.list[i], "$x))", sep="")))eval(parse(text = paste("ymin <- append(ymin, min(", genotype.list[i], "$y))", sep="")))eval(parse(text = paste("ymax <- append(ymax, max(", genotype.list[i], "$y))", sep="")))}filename <- paste("density",analysisdata[j], ".pdf", sep="")pdf(file=filename)# prepare axis of graphfigtitle <- paste("Density graph of ", analysisdata[j], sep="")plot(0, main=figtitle, type = 'n', xlab = analysisdata[j], ylab = "kerneldensity", xlim = c(min(xmin), max(xmax)), ylim = c(min(ymin), max(ymax)))# drow lines for (i in 1:numberofgenotype){eval(parse(text = paste("lines(", genotype.list[i], ",type = 'l', col = ", i, ")", sep="")))}# write legendpar(usr = c(0, 1, 0, 1))legend(0.7, 1, genotype.list, col = 1:numberofgenotype, lty = 1)dev.off()}detach("package:KernSmooth")

# histogram with percent y-axis.
# only a few arguments are actually used.
# max of ylim is percent
# percenthist(wt2fixwidth, xlim=c(0,3),ylim=c(0,40), br=seq(0,5,by=0.3),border=1,col=2)
percenthist<-function(x, breaks = "Sturges",
	percent=TRUE,
     freq = NULL, probability = !freq,
     include.lowest = TRUE, right = TRUE,
     density = NULL, angle = 45, col = NULL, border = NULL,
     main = paste("Histogram of" , xname),
     xlim = NULL, ylim = NULL,
     xlab = xname, ylab,
     axes = TRUE, plot = TRUE, labels = FALSE,
      ...)
{
	xname<-deparse(substitute(x))
	ylim<-ylim
	xlim<-xlim
	histobj<-hist(x, breaks = breaks,
     include.lowest = include.lowest , right = right,
     plot = FALSE,
      ...)
      locbreaks<-histobj$breaks
      counts<-histobj$counts
      samplenum<-length(x)
      if(!is.null(ylim))
      {
      	ymax<-samplenum*ylim[2]/100
      }
      else
      {
      	ymax<-max(counts)
      }
       if(is.null(xlim))
      {
		locxlim<-range(locbreaks)
	}
	else
	{
		locxlim<-xlim
	}
      hist(x, breaks = locbreaks, freq=TRUE,xlab = paste(xname, "N =", samplenum), main = paste("Histogram of" , xname),xlim= locxlim,ylim=c(0,ymax),include.lowest = include.lowest , right = right,angle = angle, col = col, border = border)
	ticklab<-seq(0, ymax/samplenum*100,by=5)
	tickpos<-ticklab/100* samplenum
	axis(4,at= tickpos, labels= ticklab)
	
}