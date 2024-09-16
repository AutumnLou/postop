# library
library(stabm)
library(rcartocolor)
##################################################
# Colours
palette <- carto_pal(12, "Safe")



# Rate Image

# Baseline
co <- as.data.frame(co)
rownames(co) <- co[,1]
co <- co[,-1]

# Top
top <- as.data.frame(top)
rownames(top) <- top[,1]
top <- top[,-1]

# Rate image function
rate.image <- function(mat,nms,mod.nms){
  par(mar = c(3, 5, 6, 2.1))
	image(t(mat),axes=F,col=grey(c(100:0)/100))
	ats = seq(0,1,len=ncol(mat))
	axis(3,las=2,at=ats,labels=nms)
	xgrid = ats+diff(ats)[1]/2
	ats = seq(0,1,len=nrow(mat))
	axis(2,las=2,at=ats,labels=mod.nms)
	ygrid = ats+diff(ats)[1]/2
	abline(h=ygrid,lty=3,lwd=1,col='grey')
	abline(v=xgrid,lty=1,lwd=.8,col='grey')
	return(list(xgrid=xgrid,ygrid=ygrid))
}

# Baseline
mat = t(co)
rate.image(mat=mat,rownames(co),colnames(co))

# Top
mat2 = t(top)
rate.image(mat=mat2,rownames(top),colnames(top))

# Example: say you have 4 models, 10 features
# Dummy rate matrix (overall rates per model):
#rates = matrix(runif(40),nrow=4,ncol=10)
#nms = paste("X",c(1:10),sep="")
#mod.nms = paste("model",c(1:4),sep="")


#require(stabm)

###############################################
# Indices
## Baseline
co_cox <- read.csv(file = "C:\\Users\\User\\Documents\\PhD\\JournalArticle\\Post_op\\var_select_co_cox.csv")
co_cox <- co_cox[,-1]
co_las <- read.csv(file = "C:\\Users\\User\\Documents\\PhD\\JournalArticle\\Post_op\\var_select_co_las.csv")
co_las <- co_las[,-1]
co_boo <- read.csv(file = "C:\\Users\\User\\Documents\\PhD\\JournalArticle\\Post_op\\var_select_co_boo.csv")
co_boo <- co_boo[,-1]
co_rsf <- read.csv(file = "C:\\Users\\User\\Documents\\PhD\\JournalArticle\\Post_op\\var_select_co_rsf.csv")
co_rsf <- co_rsf[,-1]

# Lists of Features used
Cox = LASSO = Boosted = RSF = NULL
for(i in 1:100){
  # Create Lists
  Cox[[i]] = which(co_cox[i,]==1)
  LASSO[[i]] = which(co_las[i,]==1)
  Boosted[[i]] = which(co_boo[i,]==1)
  RSF[[i]] = which(co_rsf[i,]==1)
}

# Hamming index
ham = c(stabilityHamming(Cox, p=8),
        stabilityHamming(LASSO, p=14),
        stabilityHamming(Boosted, p=14),
        stabilityHamming(RSF, p=8))
# Jaccard index
jac = c(stabilityJaccard(Cox),
        stabilityJaccard(LASSO),
        stabilityJaccard(Boosted),
        stabilityJaccard(RSF))
# Davis index
dav = c(stabilityDavis(Cox, p=8),
        stabilityDavis(LASSO, p=14),
        stabilityDavis(Boosted, p=14),
        stabilityDavis(RSF, p=8))
# Dice index
dice = c(stabilityDice(Cox),
         stabilityDice(LASSO),
         stabilityDice(Boosted),
         stabilityDice(RSF))

# Graph
stabs = data.frame(Jaccard=jac,Davis=dav,Dice=dice)
matplot(t(stabs),t='b',col=palette[c(7,8,11,9)],lwd=3,ylab="Rate",xlab="Index",xaxt='n', ylim = 0:1)
axis(1, at=c(1:ncol(stabs)), labels=names(stabs))
#legend("bottomleft",lty=1,lwd=2,bty='n',
#       legend=c("Cox", "LASSO", "Boosted", "RSF"),col=palette[c(7,8,11,9)])

## Top
cor_cox <- read.csv(file = "C:\\Users\\User\\Documents\\PhD\\JournalArticle\\Post_op\\var_select_comb_cor_cox.csv")
cor_cox <- cor_cox[,-1]
uni_las <- read.csv(file = "C:\\Users\\User\\Documents\\PhD\\JournalArticle\\Post_op\\var_select_comb_uni_las.csv")
uni_las <- uni_las[,-1]
cor_boo <- read.csv(file = "C:\\Users\\User\\Documents\\PhD\\JournalArticle\\Post_op\\var_select_comb_cor_boo.csv")
cor_boo <- cor_boo[,-1]
cor_uni_rsf <- read.csv(file = "C:\\Users\\User\\Documents\\PhD\\JournalArticle\\Post_op\\var_select_comb_cor_uni_rsf.csv")
cor_uni_rsf <- cor_uni_rsf[,-1]

# Lists of Features used
Cox = LASSO = Boosted = RSF = NULL
for(i in 1:100){
  # Create Lists
  Cox[[i]] = which(cor_cox[i,]==1)
  LASSO[[i]] = which(uni_las[i,]==1)
  Boosted[[i]] = which(cor_boo[i,]==1)
  RSF[[i]] = which(cor_uni_rsf[i,]==1)
}

# Hamming index
#ham = c(stabilityHamming(Cox, p=26000),
#        stabilityHamming(LASSO, p=26000),
#        stabilityHamming(Boosted, p=26000),
#        stabilityHamming(RSF, p=26000))
# Jaccard index
jac = c(stabilityJaccard(Cox),
        stabilityJaccard(LASSO),
        stabilityJaccard(Boosted),
        stabilityJaccard(RSF))
# Davis index
dav = c(stabilityDavis(Cox, p=26000),
        stabilityDavis(LASSO, p=26000),
        stabilityDavis(Boosted, p=26000),
        stabilityDavis(RSF, p=26000))
# Dice index
dice = c(stabilityDice(Cox),
         stabilityDice(LASSO),
         stabilityDice(Boosted),
         stabilityDice(RSF))

# Graph
stabs = data.frame(Jaccard=jac,Davis=dav,Dice=dice)
matplot(t(stabs),t='b',col=palette[c(7,8,11,9)],lwd=3,ylab="Rate",xlab="Index",xaxt='n', ylim = 0:1)
axis(1, at=c(1:ncol(stabs)), labels=names(stabs))
legend("topright",lty=1,lwd=2,bty='n',
       legend=c("Cox", "LASSO", "Boosted", "RSF"),col=palette[c(7,8,11,9)])




################################################
# Dummy rate matrix (rates per bootstrap sample for a given model):
rates = matrix(sample(c(0,1),size=1000,replace=TRUE),nrow=100,ncol=10)


# turn rates matrix into a list of feature indices per bootstrap set:
sels.model1 = sels.model2 = sels.model3 = sels.model4 = NULL
n = nrow(rates)
for(i in 1:n){
	# these are just very dummy ways of creating toy examples
	sels.model1[[i]] = which(rates[i,]==1)
	sels.model2[[i]] = c(1:5)
	sels.model3[[i]] = unique(c(1,which(sample(1:10,size=10,replace=TRUE)<5)))
	sels.model4[[i]] = sample(1:10,1)
}

# Hamming index
P = 10 # total number of features
ham = c(stabilityHamming(sels.model1,p=P),
	stabilityHamming(sels.model2,p=P),
	stabilityHamming(sels.model3,p=P),
	stabilityHamming(sels.model4,p=P))
# Jaccard index
jac = c(stabilityJaccard(sels.model1),
	stabilityJaccard(sels.model2),
	stabilityJaccard(sels.model3),
	stabilityJaccard(sels.model4))
# Davis index
dav = c(stabilityDavis(sels.model1,p=P),
	stabilityDavis(sels.model2,p=P),
	stabilityDavis(sels.model3,p=P),
	stabilityDavis(sels.model4,p=P))
# Dice index
dice = c(stabilityDice(sels.model1),
	stabilityDice(sels.model2),
	stabilityDice(sels.model3),
	stabilityDice(sels.model4))
# Wald index
wald = c(stabilityWald(sels.model1,p=P),
	stabilityWald(sels.model2,p=P),
	stabilityWald(sels.model3,p=P),
	stabilityWald(sels.model4,p=P))

stabs = data.frame(Hamming=ham,Jaccard=jac,Davis=dav,Dice=dice,Wald=wald)
matplot(t(stabs),t='b',col=c(1:4),lwd=3,ylab="Rate",xlab="Index",xaxt='n')
axis(1, at=c(1:ncol(stabs)), labels=names(stabs))
legend("bottomleft",lty=1,lwd=2,bty='n',
	legend=paste("model",c(1:4),sep=""),col=c(1:4))
