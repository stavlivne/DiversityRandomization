##Calculate real beta diversity
#Load data (abundance data)
#You can pre devide your data into seperate dataframes
#one for each treatment combination, and import each one seperatly
#or you can import all the data into one dataframe and seperate them in R

#import the data from the file
data <-read.csv(file.choose(), header=T,row.names=1)

#select R working library
library(vegan) 
library(plotrix)

# If you import all your data at ones you can subset into different combinations here
#make sure to adjast the values according to your data
#if you pre devided your data skip this stage
HONS <- t(rowSums(data[,grep("HONS",names(data))]))
HOS <- t(rowSums(data[,grep("HOS",names(data))]))
HNS <- t(rowSums(data[,grep("HNS",names(data))]))
HS <- t(rowSums(data[,grep("HS",names(data))]))

#This is the part where you actually compute beta diversity
#calculating the different components of beta diversity using betapart

H_HONS <-diversity(HONS, index = "shannon", MARGIN = 1, base = exp(1))
H_HOS <-diversity(HOS, index = "shannon", MARGIN = 1, base = exp(1))
H_HNS <-diversity(HNS, index = "shannon", MARGIN = 1, base = exp(1))
H_HS <-diversity(HS, index = "shannon", MARGIN = 1, base = exp(1))

S_HONS <-specnumber(HONS)
S_HOS <-specnumber(HOS)
S_HNS <-specnumber(HNS)
S_HS <-specnumber(HS)

Beta_HONS <-H_HONS/log(S_HONS)
Beta_HOS  <-H_HOS/log(S_HOS)
Beta_HNS  <-H_HNS/log(S_HNS)
Beta_HS   <-H_HS/log(S_HS)

names=c("Homo live","Homo sterile","Hetero live","Hetero sterile")
#combine all into one dataframe
mean_data <-cbind(Beta_HONS,Beta_HOS,Beta_HNS,Beta_HS)
colnames(mean_data)<-names

#stacked bar plot
plotTop<- max(mean_data)
plotBeta <-barplot(mean_data,names.arg = names,
                   ylab="Beta diversity",ylim=c(0,plotTop*2),legend.text = FALSE, col=grey.colors(2))

# calculating real diff between treatments
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

#replicate(1000, {
for (i in 1:1000){
  rand <- data[,sample(ncol(data))]
  names(rand)<- names(data)
  rand_HONS <- t(rowSums(rand[,grep("HONS",names(rand))]))
  rand_HOS <- t(rowSums(rand[,grep("HOS",names(rand))]))
  rand_HNS <- t(rowSums(rand[,grep("HNS",names(rand))]))
  rand_HS <- t(rowSums(rand[,grep("HS",names(rand))]))
  
  #This is the part where you actually compute mean diversity
  #calculating the different components of beta diversity using betapart
  rand_H_HONS <-diversity(rand_HONS, index = "shannon", MARGIN = 1, base = exp(1))
  rand_H_HOS <-diversity(rand_HOS, index = "shannon", MARGIN = 1, base = exp(1))
  rand_H_HNS <-diversity(rand_HNS, index = "shannon", MARGIN = 1, base = exp(1))
  rand_H_HS <-diversity(rand_HS, index = "shannon", MARGIN = 1, base = exp(1))
  
  rand_S_HONS <-specnumber(rand_HONS)
  rand_S_HOS <-specnumber(rand_HOS)
  rand_S_HNS <-specnumber(rand_HNS)
  rand_S_HS <-specnumber(rand_HS)
  
  rand_Beta_HONS <-rand_H_HONS/log(rand_S_HONS)
  rand_Beta_HOS  <-rand_H_HOS/log(rand_S_HOS)
  rand_Beta_HNS  <-rand_H_HNS/log(rand_S_HNS)
  rand_Beta_HS   <-rand_H_HS/log(rand_S_HS)
  
  names=c("Homo live","Homo sterile","Hetero live","Hetero sterile")
  #combine all into one dataframe
  rand_Beta_data <-cbind(rand_Beta_HONS,rand_Beta_HOS,rand_Beta_HNS,rand_Beta_HS)
  colnames(rand_Beta_data)<-names
  
  # calculating diff between rand treatments
  rand_HONS_HOS <- abs(rand_Beta_data[1]-rand_Beta_data[2])
  rand_HONS_HNS <- abs(rand_Beta_data[1]-rand_Beta_data[3])
  rand_HONS_HS <-abs(rand_Beta_data[1]-rand_Beta_data[4])
  rand_HOS_HNS <-abs(rand_Beta_data[2]-rand_Beta_data[3])
  rand_HOS_HS <-abs(rand_Beta_data[2]-rand_Beta_data[4])
  rand_HNS_HS <-abs(rand_Beta_data[3]-rand_Beta_data[4])
  rand_comp <- as.data.frame(cbind(rand_HONS_HOS,rand_HONS_HNS,rand_HONS_HS,rand_HOS_HNS,rand_HOS_HS,rand_HNS_HS))
  rand_results <- rbind(rand_results,rand_comp)
}
rand_results <- rand_results[2:1001,]
###calculating p value
#replicate the real data 1000 times
real_results <- real_results[rep(seq_len(nrow(real_results)), each=1000),]
#calculate the p.value as the nuber of times the real diff in ratio was equal to
#or larger than the randomize diff
p.values <-as.data.frame(colSums(rand_results>=real_results))/1000 
write.csv(p.values,file="Pvalues_betadiversity_Gil.csv")
write.csv(mean_data,file="betadiversity_data_gil.csv")
