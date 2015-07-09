
#'@title Anderson-Darling Normality Test (norm.p)
#'@description Test for normality
#'@details Anderson-Darling Normality Test is used to determine whether a set of observations follows ‘Normal Distribution’.  The assumption of ‘Normality’ is widely used in Statistics in the areas of Inferential Statistics, Parametric Methods, Modeling (i.e.) from Simple Linear Regression to advanced Predictive Modeling Techniques.
#'@param Numeric column position x in the dataframe
#'@return P-value
#'@export


norm.p <- function(my.data) {
  
  my.matrix <- data.frame(matrix(data = 0,nrow = nrow(my.data),ncol = 7))
  
  my.matrix[,1] <- data.frame(sort(my.data[,1],decreasing = F))
  my.matrix[,2] <- pnorm(scale(my.matrix[,1],center = T,scale = T),0,1,T)
  my.matrix[,3] <- 1- (pnorm(scale(my.matrix[,1],center = T,scale = T),0,1,T))
  my.matrix[,4] <- sort(my.matrix[,3],decreasing = F)
  
  for (t in 1:nrow(my.matrix)) {
    my.matrix[t,5] <- (2*t-1) * (log(my.matrix[t,4])+log(my.matrix[t,2]))
    my.matrix[t,6] <- qnorm((t-3/8) / (nrow(my.matrix)+1/4),mean(my.matrix[,1]),sd(my.matrix[,1]))
    my.matrix[t,7] <- (t-3/8) / (nrow(my.matrix)+1/4)
  }
  
  si <- sum(my.matrix[,5])
  
  ad.stat <- -si / nrow(my.matrix)-nrow(my.matrix)
  
  ad.stat.star <- ad.stat * (1+0.75/nrow(my.matrix)+2.25/nrow(my.matrix)^2)
  
  if (ad.stat.star >= 0.6 && ad.stat.star < 13) {p1 <- exp(1.2937-5.709*ad.stat.star+0.0186*ad.stat.star^ 2) } else { p1<-0}
  
  if (ad.stat.star >= 0.34 && ad.stat.star <0.6) {p2 <- exp(0.9177-4.279*ad.stat.star-1.38*ad.stat.star^2)} else { p2<-0}
  
  if (ad.stat.star >= 0.2 && ad.stat.star <0.34) {p3 <- 1-exp(-8.318+42.796*ad.stat.star-59.938*ad.stat.star^2)} else {p3<-0}
  
  if (ad.stat.star < 0.2) {p4 <- 1-exp(-13.436+101.14*ad.stat.star-223.73*ad.stat.star^2)} else {p4<-0}
  
  if(max(p1,p2,p3,p4)<0.0005) {
    p.value <- 0.0005 } else { p.value <- max(p1,p2,p3,p4)}
  return(p.value)
  
}





#'@title Graphical Summary(gs)
#'@description Gives four types of graphical summary Histogram,Box plot,Time series plot and Auto correlation function plot

#'@details Histogram displays the values in Frequency Domain.  Histogram also displays information on Anderson-Darling Normality P-value, Mean (mu) & Standard Deviation (sigma).Time Series Plot, displays the values over Time Domain.  Time Plot displays values of Q1, Q2/Median & Q3.  These values are also displayed in the plot in the form of Horizontal Lines.Box-Whisker Plot, is a Five-Number summary.  Outliers and Skewness of the distribution can be assed using Box-Plot.Auto Correlation Function, shows the significant of each lags.
#'@param Numeric (column position x in the dataframe)
#'@return plots
#'@export

gs <- function(col.pos.x) {
  
  options(digits=5,scipen=10)
  
  mydata <- data.frame(col.pos.x)
  col.pos <- 1
  col.name <- colnames(mydata) [col.pos]
  
  par(mfcol=c(2,2))
  
  
  
  hist(mydata[,col.pos],col='tomato',main=paste('Histogram'),probability=T,xlab=paste('Nor.p.val',round(norm.p(mydata[col.pos]),4),'mu=',round(mean(mydata[,col.pos]),4),'sigma=',round(sd(mydata[,col.pos]),4)))
  
  rug(mydata[,col.pos])
  
  
  
  boxplot(mydata[,col.pos],col='purple',horizontal=T,main='Box-Whisker Plot',xlab=paste(' *Min=',min(mydata[,col.pos]),'Max=',max(mydata[,col.pos])))
  
  abline(v = mean(mydata[,col.pos]),col='red',lwd=2)
  
  xlab.v <- paste('Q1=',quantile(mydata[,col.pos],.25),'Q2/Median=',median(mydata[,col.pos]),'Q3=',quantile(mydata[,col.pos],.75))
  
  plot(mydata[,col.pos],type='b',main='Time Plot',ylab=col.name,xlab=xlab.v,col='purple')
  
  abline(h = median(mydata[,col.pos]),col='red')
  abline(h = quantile(mydata[,col.pos],.25),col='red')
  abline(h = quantile(mydata[,col.pos],.75),col='red')
  
  acf(mydata[,col.pos],main='Auto Correlation Function',col='red',type = 'correlation')
  
  title(paste('Graphical Summary of ',col.name),line = -1,outer = TRUE)
  
}






#'@title Deviations (my.deviate)
#'@description This functions addresses Deviation, MAD & MSD.
#'@param Four numeric parameters Parameter 1: It is a column of a Data Frame,Parameter 2: If the second parameter = 1, it then displays DEVIATIONS for each observation,Parameter 3: If the third parameter = 1, it then displays column required to compute MAD,Parameter 4: If the fourth parameter = 1, it then displays column required to compute MSD
#'@return Displays the computed columns
#'@export 
my.deviate <- function(col.pos.x,Deviate=0,MeanAbsDeviation=1,MeanSqDev=0) {
  
  col.pos <- as.numeric(1) # changes finalised
  mydata <- data.frame(col.pos.x) # changes finalised
  
  
  
  my.matrix <- data.frame(matrix(data = 0,nrow = nrow(mydata),ncol = 3))
  
  my.mean <- mean(mydata[,col.pos])
  
  for (t in 1:nrow(mydata)) {
    
    if(Deviate==1) {
      my.matrix[t,1] <- mydata[t,col.pos] - my.mean } # End of Deviate
    
    if(MeanAbsDeviation==1) {
      my.matrix[t,2] <- abs(mydata[t,col.pos] - my.mean) } # End of MAD
    
    if(MeanSqDev==1) {
      my.matrix[t,3] <- (mydata[t,col.pos] - my.mean) ^ 2 } # End of MSD
    
    
  }
  
  
  
  my.args <- array(data = 1,dim = c(6,1))
  my.args[3] <- MeanSqDev
  my.args[2] <- MeanAbsDeviation
  my.args[1] <- Deviate
  
  my.args[6] <- "cMSD"
  my.args[5] <- "cMAD"
  my.args[4] <- "cDEV"
  
  for (t in 1:3) {
    
    if (my.args[t]==1) {
      mydata <- data.frame(mydata,my.matrix[,t])
      colnames(mydata) [ncol(mydata)] <- paste(names(mydata[1]),"_",my.args[t+3],sep = "")
    }
    
  }
  
  return (mydata)
}


mad.d <- function(col.poss) {
  
  temp <- my.deviate(col.pos = col.poss,Deviate = 0,MeanAbsDeviation = 1,MeanSqDev = 0)
  
  return(sum(temp[,ncol(temp)]) / nrow(temp))
  
}


msd.d <- function(col.poss) {
  
  temp <- my.deviate(col.pos = col.poss,Deviate = 0,MeanAbsDeviation = 0,MeanSqDev = 1)
  
  return(sum(temp[,ncol(temp)]) / nrow(temp))
  
}




#'@title Creating Indicator Variables (create.iv)
#'@description Creates three types of indicator variables 
#'@details Indicator Variables is a method that is widely used in modeling particularly to handle data that are NOMINAL in nature.  It converts NOMINAL data into a meaningful numerical variables knows as ‘Indicator Variables’ (IV).  The other names for Indicator Variables are also known as Dummy Variables (DV).At times, Indicator Variables (IV) were rarely applied even on NOMINAL data, to convert them into meaningful numerical variables.Indicator Variables (IVs) is widely applied in Regression Techniques, Machine Learning Techniques & Multivariate Techniques such as Cluster Analysis, etc.
#'@param Numeric column position x in the dataframe and the type of indicator variable For REFERENCE CODING use ‘0’ ,ForEFFECT CODING use ‘-1’,For NO Y-INTERCEPT use ‘1’ and logical statement T or F to retain or drop the original column
#'@return Indicator Variables Columns
#'@export

create.iv <- function(col.pos.x,iv.type,retain.col=T) {
  
  mydata <- col.pos.x # recent insertion *parameter col.pos.x was col.pos
  col.pos <- 1 # recent insertion
  
  iv.colnam <- colnames(mydata) [col.pos]
  iv.uv <- data.frame(unique(mydata[,col.pos]))
  colnames(iv.uv) [1] <- colnames(mydata) [col.pos]
  
  
  if(iv.type==0) { # Base ZERO coding style
    
    temp.mat <- matrix(data = 0,nrow = nrow(mydata),ncol=(nrow(iv.uv)-1))
    temp.mat <- data.frame(temp.mat)
    
    for (i in 1:(nrow(iv.uv)-1)) {
      colnames(temp.mat) [i] <- paste(iv.colnam,i,sep = "_")
    }
    
    for (j in 1:nrow(mydata)) {
      
      for (k in 1:(nrow(iv.uv)-1)) {
        
        if(mydata[j,col.pos]==iv.uv[k+1,1]) {
          temp.mat[j,k] <- 1
        }
        
      }
      
    }
  } else {
    
    if(iv.type==1) { # n * n MATRIX coding style
      
      
      temp.mat <- matrix(data = 0,nrow = nrow(mydata),ncol=(nrow(iv.uv)-0))
      temp.mat <- data.frame(temp.mat)
      
      for (i in 1:(nrow(iv.uv)-0)) {
        colnames(temp.mat) [i] <- paste(iv.colnam,iv.uv[i,1],sep = "_")
      }
      
      for (j in 1:nrow(mydata)) {
        
        for (k in 1:(nrow(iv.uv)-0)) {
          
          if(mydata[j,col.pos]==iv.uv[k+0,1]) {
            temp.mat[j,k] <- 1
          }
          
        }
        
      } } else {
        
        
        if(iv.type==-1) { # Average coding style
          
          temp.mat <- matrix(data = 0 ,nrow = nrow(mydata),ncol=(nrow(iv.uv)-1))
          temp.mat <- data.frame(temp.mat)
          
          for (i in 1:(nrow(iv.uv)-1)) {
            colnames(temp.mat) [i] <- paste(iv.colnam,i,sep = "_")
          }
          
          for (j in 1:nrow(mydata)) {
            
            for (k in 1:(nrow(iv.uv)-1)) {
              
              if(mydata[j,col.pos]==iv.uv[1,1])   {
                temp.mat[j,k] <- -1
              }   else {                             
                
                
                if(mydata[j,col.pos]==iv.uv[k+1,1]) {
                  temp.mat[j,k] <- 1               
                } 
              }
              
              
              
              
            }
            
          }
        } }
  }
  
  
  mydata <- cbind(mydata,temp.mat) 
  
  
  if(retain.col ==T)
    return(mydata)
  else
    return(mydata[,-1])
  
}



# xTOz
#'@title XtoZ (stdZ)
#'@description Standardize 
#'@details Transforms originalvalue into its respective "Z" value
#'@param Numeric column position x in the dataframe 
#'@return Returns standardized column
#'@export

stdZ <- function(pos.x,replace=0) {
  
  mydata <- pos.x # recent insertion
  pos <- 1 # recent insertion
  
  f.o <- scale(mydata[,pos],center = T,scale = T)
  
  if (replace==0) {
    mydata <- data.frame(mydata,"TEMP"=f.o)
    colnames(mydata) [ncol(mydata)] <- paste(colnames(mydata) [pos],".Z",sep = "")
  } else {
    
    mydata[,pos] <- f.o
    colnames(mydata) [pos] <- paste(colnames(mydata) [pos],".Z",sep = "")
  }
  
  return(mydata)
  
}



#'@title  Skewness (ske)
#'@description Measure of Descriptive Statistics
#'@details Skewness addresses the characteristics of symmetry.
#'@param Numeric column position x in the dataframe
#'@return Gives the skewness value
#'@export

ske <- function(column.pos) {
  
  val <- 0
  
  temp <- my.deviate(col.pos=column.pos,Deviate = 1,MeanAbsDeviation = 0,MeanSqDev = 0)
  tot.col <- ncol(temp)
  temp.sd <- sd(temp[,tot.col])
  tot.row <- nrow(temp)
  
  for(t in 1:tot.row) {
    val <- val + ((temp[t,tot.col] / temp.sd) ^ 3)
  } 
  return((tot.row) / ((tot.row-1) * (tot.row-2)) * val) 
}



# 
#'@title Kurtosis (kur)
#'@description Measure of Descriptive Statistics
#'@details Kurtosis addresses the characteristics of ‘Peakedness’ & ‘Taildness’.
#'@param Numeric (column position x in the dataframe)
#'@return Gives the kurtosis value
#'@export


kur <- function(column.pos) {
  
  val <- 0
  
  temp <- my.deviate(col.pos=column.pos,Deviate = 1,MeanAbsDeviation = 0,MeanSqDev = 0)
  tot.col <- ncol(temp)
  temp.sd <- sd(temp[,tot.col])
  tot.row <- nrow(temp)
  
  for(t in 1:tot.row) {
    val <- val + ((temp[t,tot.col] / temp.sd) ^ 4)
  } 
  
  p1 <- tot.row * (tot.row+1)
  p2 <- (tot.row-1) * (tot.row-2) * (tot.row-3)
  p3 <- 3 * ((tot.row-1)^2)
  p4 <- (tot.row-2) * (tot.row-3)
  
  return(((p1/p2) * val) - (p3/p4))
  
}




optimal.lamda <- function(col.pos.x,my.alpha=0.05) {
  
  mydata <- col.pos.x # recent insertion
  col.pos <- 1
  
  options(scipen=100)
  my.row.id <- 0
  my.matrix <- data.frame(matrix(data = 0,nrow = 400 ,ncol = 2))
  
  for (j in seq(from = -2,to = 2,by = .01)) {
    
    if(j==0) {
      
    } else {
      
      my.row.id <- (my.row.id+1)
      
      my.matrix[my.row.id,1] <- as.numeric(j)
      
      
      my.matrix[my.row.id,2] <- as.numeric(norm.p((mydata[col.pos]^j)))
      
    }
  }
  
  if(round(my.matrix[which.max(my.matrix[,2]),1],digits = 4)<my.alpha) {
    return(0) } else {
      return(round(my.matrix[which.max(my.matrix[,2]),1],digits = 4))
      
    }
}



#'@title  Box-Cox Transformation (boxcox.t)
#'@description  Box-Cox Power Transformations from Non-Normal to Normality
#'@details Box-Cox transformation is usually applied in order to achieve modeling assumptions.  As stated earlier techniques such as Simple Linear Regression, Multiple Linear Regression, Logistic Regression and other Classification Techniques like Discriminant Analysis, DT, NN, etc., also requires data to be normally distributed.(E.g.) In ENERGY model building, the usage of ENERGY in any house, district or state will be right-skewed.  Using Box-Cox when the same is transformed into ‘NORMAL’ the assumption is met, at the same time after transformation that particular units will lose its original UNIT OF MEASUREMENT.
#'@param Numeric column position x in the dataframe
#'@return Returns optimal lambda value and transformed column vector
#'@export

boxcox.t <- function(your.object) {
  
  cat('\n','Optimal Lamda = ',optimal.lamda(your.object),'\n')
  return(your.object^optimal.lamda(your.object))  
  
}



#'@title Transformation of Poisson Data (poi.trans)
#'@description  Freeman Tukey Variance Stabilization transformation for POISSON DISTRIBUTION DATA
#'@details This technique widely used to Stabilize Variance during model building.if the COUNT data follows Poisson Distribution, this transformation is then applied.
#'@param Numeric column position x in the dataframe
#'@return Returns transformed column vector
#'@export 
poi.trans <- function(your.object) {
  options(digits=6)
  mydata <- your.object
  
  mydata.poi.trans <- ((mydata^0.5) +  ((mydata+1)^0.5)) / 2
  
  mydata <- data.frame(cbind(mydata,mydata.poi.trans))
  colnames(mydata) [2] <- paste(names(your.object),'.POISSONtrans',sep = "")
  
  return (mydata)
  
}



#
#'@title Transformation of Proportion Data (prop.trans)
#'@description Freeman Tukey Arcsin Transformation For Proportions
#'@details This technique is widely used to Stabilize Variance during model building.  The proportion should be arrived from a COUNT data only.  It should not be applied under scenarios like where the Proportion is arrived from Continuous Data.  Should not be applied for Profit Margins.
#'@param Numeric column position x in the dataframe
#'@return Returns transformed column vector
#'@export 

prop.trans <- function(your.object) {
  options(digits=6)
  mydata <- your.object
  
  mydata.prop.trans <- asin(mydata^0.5) 
  
  mydata <- data.frame(cbind(mydata,mydata.prop.trans))
  colnames(mydata) [2] <- paste(names(your.object),'.PROPORTIONtrans',sep = "")
  
  return (mydata)
  
}
x<-c(.3,.76,.45)
prop.trans(x)

