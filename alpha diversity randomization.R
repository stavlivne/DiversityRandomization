###randomization of alpha diversity
#import the data from the file
data <-read.csv(file.choose(), header=T,row.names=1)

#load R working library
library(vegan) 
library(plotrix)


#subset your data into seperate dataframes (a dataframe for each treatment combination)
#make sure to adjast the names according to your data
HONS <- t(data[,grep("HONS",names(data))])
HOS <- t(data[,grep("HOS",names(data))])
HNS <- t(data[,grep("HNS",names(data))])
HS <- t(data[,grep("HS",names(data))])

#calculating Simpson diversity, you can do this for any type of diversity calculation e.g. specnum/shannon/ Pielo's etc.

meanSim_HONS<-mean(diversity(HONS, index = "simpson", MARGIN = 1, base = exp(1)))
meanSim_HOS <-mean(diversity(HOS, index = "simpson", MARGIN = 1, base = exp(1)))
meanSim_HNS <-mean(diversity(HNS, index = "simpson", MARGIN = 1, base = exp(1)))
meanSim_HNS_HS  <-mean(diversity(HS, index = "simpson", MARGIN = 1, base = exp(1)))

names=c("Homogeneous natural","Homogeneous autoclaved","Heterogeneous natural","Heterogeneous autoclaved")
#combine all into one dataframe
mean_data <-cbind(meanSim_HONS,meanSim_HOS,meanSim_HNS,meanSim_HS)
colnames(mean_data)<-names

#stacked bar plot
plotTop<- max(mean_data)
plotBeta <-barplot(mean_data,names.arg = names,
                   ylab="mean alpha diversity",ylim=c(0,plotTop*2),legend.text = FALSE, col=grey.colors(2))

real_HONS_HOS <- abs(mean_data[1]-mean_data[2])
real_HONS_HNS <- abs(mean_data[1]-mean_data[3])
real_HONS_HS <-abs(mean_data[1]-mean_data[4])
real_HOS_HNS <-abs(mean_data[2]-mean_data[3])
real_HOS_HS <-abs(mean_data[2]-mean_data[4])
real_HNS_HS <-abs(mean_data[3]-mean_data[4])
real_results <- as.data.frame(cbind(real_HONS_HOS,real_HONS_HNS,real_HONS_HS,real_HOS_HNS,real_HOS_HS,real_HNS_HS))

####randomizing 
rand_results<-matrix(nrow=1,ncol=6)
colnames(rand_results)<- c("rand_HONS_HOS","rand_HONS_HNS","rand_HONS_HS","rand_HOS_HNS","rand_HOS_HS","rand_HNS_HS")

#replicate(10000, {
for (i in 1:10000){
  rand <- data[,sample(ncol(data))] #randomizing samples among treatments
  names(rand)<- names(data)
  rand_HONS <- t(rand[,grep("HONS",names(rand))])
  rand_HOS <- t(rand[,grep("HOS",names(rand))])
  rand_HNS <- t(rand[,grep("HNS",names(rand))])
  rand_HS <- t(rand[,grep("HS",names(rand))])
  
  #This is the part where you actually compute mean diversity
  #calculating the different components of beta diversity using betapart
  rand_meanSim_HONS <-mean(diversity(rand_HONS, index = "simpson", MARGIN = 1, base = exp(1)))
  rand_meanSim_HOS  <-mean(diversity(rand_HOS, index = "simpson", MARGIN = 1, base = exp(1)))
  rand_meanSim_HNS  <-mean(diversity(rand_HNS, index = "simpson", MARGIN = 1, base = exp(1)))
  rand_meanSim_HS   <-mean(diversity(rand_HS, index = "simpson", MARGIN = 1, base = exp(1)))
  
  names=c("Homogeneous natural","Homogeneous autoclaved","Heterogeneous natural","Heterogeneous autoclaved")
  #combine all into one dataframe
  rand_meanSim_data <-cbind(rand_meanSim_HONS,rand_meanSim_HOS,rand_meanSim_HNS,rand_meanSim_HS)
  colnames(rand_meanSim_data)<-names
  
  # calculating diff between rand treatments
  rand_HONS_HOS <- abs(rand_meanSim_data[1]-rand_meanSim_data[2])
  rand_HONS_HNS <- abs(rand_meanSim_data[1]-rand_meanSim_data[3])
  rand_HONS_HS <-abs(rand_meanSim_data[1]-rand_meanSim_data[4])
  rand_HOS_HNS <-abs(rand_meanSim_data[2]-rand_meanSim_data[3])
  rand_HOS_HS <-abs(rand_meanSim_data[2]-rand_meanSim_data[4])
  rand_HNS_HS <-abs(rand_meanSim_data[3]-rand_meanSim_data[4])
  rand_comp <- as.data.frame(cbind(rand_HONS_HOS,rand_HONS_HNS,rand_HONS_HS,rand_HOS_HNS,rand_HOS_HS,rand_HNS_HS))
  rand_results <- rbind(rand_results,rand_comp)
}
rand_results <- rand_results[2:10001,]
###calculating p value
#replicate the real data 10000 times
real_results <- real_results[rep(seq_len(nrow(real_results)), each=10000),]
#calculate the p.value as the nuber of times the real diff in ratio was equal to
#or larger than the randomize diff
p.values <-as.data.frame(colSums(rand_results>=real_results))/10000 
write.csv(p.values,file="Pvalues.csv")
write.csv(mean_data,file="diversity.csv")

