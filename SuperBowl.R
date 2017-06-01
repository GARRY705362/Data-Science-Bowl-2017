#libraries used for the analysis
library(oro.dicom)
library(dplyr)
library(EBImage)
library(h2o)

#visualization of the images of a particular patient
Size=32  #Set the image size for reduction / expansion
setwd="D:/Predictive analysis/super bowl/SampleImages"
fileList=dir("D:/Predictive analysis/super bowl/SampleImages/00cba091fa4ad62cc3200a657aeb957e/", recursive=TRUE)
l=length(fileList) 
z=array(dim=c(length(fileList),Size,Size)) #defining an array
face=array(dim=c(length(fileList),Size*Size))
Hdata=array(dim=c(l,42))
par(mai=c(.05,.05,.05,.05))
par(mfrow=c(10,14)) 
#collapse the data and 
for (i in 1:l)
{ #l is length 
  a=paste("D:/Predictive analysis/super bowl/SampleImages/00cba091fa4ad62cc3200a657aeb957e/", fileList[i],sep="") 
  mydata=readDICOMFile(a) 
  y=resize(t(mydata$img), w=Size,h=Size) # resizing to 32x32
  image(y,col=gray(0:255/256), axes=FALSE, xlab="", ylab="") 
  z[i,,]=imageData(y) 
  face[i,]=imageData(y) 
  #removing from the memory
  rm(y)
  rm(h)
}

#making the train and test data. This section of code is for creating a flat file from sample images, but can be used
#to create flat files from stage 1(train) and stage 2(test) data.
LOP=list.files("D:/Predictive analysis/super bowl/SampleImages/")
for(Alpha in 1:length(LOP))
{  
 LOI=list.files(paste("D:/Predictive analysis/super bowl/SampleImages/",LOP[Alpha],sep = "")) 
 Beta=array(dim=c(length(LOI),32,32))
 Hdata=array(dim=c(length(LOI),3))
  for (i in 1:length(LOI))
  {
    DCM=paste("D:/Predictive analysis/super bowl/SampleImages/",LOP[Alpha],"/",LOI[i],sep="") #path
    mydata=readDICOMFile(DCM)
    newsize=resize(t(mydata$img), w=32,h=32) 
    Beta[i,,]=imageData(newsize)
    Gamma <- match(c("PatientID","SliceLocation"),mydata$hdr$name )
    temp=as.vector(c(mydata$hdr$value[Gamma],LOI[i]))
    Hdata[i,]=temp
    rm(newsize)
    rm(temp)
  }
  if(Alpha==1)
  {
    a=as.data.frame(z)
    b=as.data.frame(Hdata)
    colnames(b) <- c("PatientID","SliceLocation","image_id")
    Delta=cbind(b,a)
  }
  else
  {
    a=as.data.frame(z)
    b=as.data.frame(Hdata)
    colnames(b) <- c("PatientID","SliceLocation","image_id")
    Gamma <- cbind(b,a)
    Delta=rbind(Delta,Gamma)
  }
}

# writing the test and the train data into csv files
write.csv(Delta,"D:/Predictive analysis/super bowl/TrainData.csv") 
write.csv(Delta,"D:/Predictive analysis/super bowl/TestData.csv") 

#read TrainData set
TrainData <- read.csv("D:/Predictive analysis/super bowl/TrainData.csv",Hdata = TRUE)

# read Lablels for patients
LB <- read.csv("D:/Predictive analysis/super bowl/stage1_LB.csv",Hdata = TRUE)
colnames(LB)[1] <- "PatientID"
TrainData$Cancer.1 <- NULL
TrainData$Slicelocation<- NULL

#merge data with Labels 
TrainD1 <- TrainData
TrainD1 <- merge(TrainData,LB,by="PatientID")
TrainD1$SliceLocation <- NULL
TrainD1 <- aggregate(TrainD1[,c(-1,-2)],list(TrainD1$PatientID,TrainD1$Cancer),mean)

#data modifications
TrainD1$X<- NULL
TrainD1$X.1<- NULL
TrainD1$image_id <- NULL

#writing the new train data into a csv
write.csv(TrainD1,"D:/Predictive analysis/super bowl/TrainD1.csv")
DeltaNA <- subset(Delta,is.na(Delta$SliceLocation))

#data modifications
DeltaNA$SliceLocation <- NULL
DeltaNA$X <- NULL
DeltaNA$image_id<- NULL
DeltaNA <- aggregate(DeltaNA[,-1],list(DeltaNA$PatientID),mean)
colnames(DeltaNA)[1] <- c("PatientID")
TestD1<- subset(Delta,!is.na(Delta$SliceLocation))
TestD1$SliceLocation <- round(TestD1$SliceLocation/100)
TestD1$X <- NULL
TestD1$image_id <-  NULL
TestD1 <- aggregate(TestD1[,c(-1,-2)],list(TestD1$PatientID,TestD1$SliceLocation),mean)
colnames(TestD1)[1:2] <- c("PatientID","Slicelocation")
TestD1$Slicelocation <- as.factor(TestD1$Slicelocation)

#writing the new test data into a csv
write.csv(Delta,"D:/Predictive analysis/super bowl/TestD1.csv")

TrainData1 <- read.csv("D:/Predictive analysis/super bowl/TrainD1.csv",Hdata = TRUE)  


TestData1 <- read.csv("D:/Predictive analysis/super bowl/TestD1.csv")

#modelling
h2o.init(nthreads = -1)


Trainh2o <- as.h2o(TrainData1)
Testh2o <- as.h2o(TestData1)




#Trainh2o <- as.h2o(TrainData1)
#Testh2o <- as.h2o(TestData1)

#DL <- h2o.deeplearning(x=c(1,3:1026),y=2,TrainDataing_frame = Trainh2o,hidden = c(20,20),epochs = 100)
#best performance
RF <- h2o.randomForest(x=c(1,3:1026),y=2,TrainDataing_frame = Trainh2o,ntrees = 600,max_depth = 17,sample_rate = 0.95)
#tuning parameters
RF <- h2o.randomForest(x=c(1,3:1026),y=2,TrainDataing_frame = Trainh2o,ntrees = 750,max_depth = 20,sample_rate = 0.95)
RF <- h2o.randomForest(x=c(1,3:1026),y=2,TrainDataing_frame = Trainh2o,ntrees = 50,max_depth = 35,sample_rate = 0.95)

DL <- h2o.deeplearning(x=c(3:1025),y=2,TrainDataing_frame = Trainh2o,hidden = c(10),epochs = 10,standardize = TRUE)
#tuning parameters
#DL <- h2o.deeplearning(x=c(1,3:1026),y=2,TrainDataing_frame = Trainh2o,hidden = c(20,20),epochs = 100)
#DL <- h2o.deeplearning(x=c(3:1025),y=2,TrainDataing_frame = Trainh2o,ntrees = 100)



h2o.performance(model = RF)
h2o.performance(model = DL)

result <- h2o.predict(RF)
result <- as.data.frame(result)

result <- round(as.data.frame(result),digits = 1)
#Output <- cbind(as.data.frame(Testh2o)[,3],result)
Output <- cbind(test$PatientID,result)

#writing the output file
write.csv(Output,"D:/Predictive analysis/super bowl/Outputput.csv")

#shutting down h2o environment
h2o.shutdown()

