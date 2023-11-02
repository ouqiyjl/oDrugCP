rm(list = ls())
setwd("E:\\dxm\\BLCA")

exp <- read.table("exp_all.txt",sep = "\t",header = T,stringsAsFactors = F) #导入表达谱，包括miRNA与gene
path <- read.table("risk_pathway.txt",sep = "\t",header = T,stringsAsFactors = F) #导入通路
PPI <- read.table("ppi.txt",sep = "\t",header = T,stringsAsFactors = F) #导入实验证实的protein-protein网络
mir_ge <- read.table("miRNA-gene.txt",sep = "\t",header = T,stringsAsFactors = F) #导入miRNA-靶基因网络
PPI[,3:4] <- lapply(PPI[,3:4],as.character) #将PPI网络的entrez列转为字符型
tumor <- grep("Tumor.*",colnames(exp)) #找出表达谱的癌症样本
control <- grep("Normal.*",colnames(exp)) #找出表达谱的正常样本

#列出所有的通路组合
all_path <- NULL
for(i in 1:(nrow(path)-1)){
  for (j in (i+1):nrow(path)){
    xx <- cbind(path[i,c(1:3)],path[j,c(1:3)])
    colnames(xx) <- ""
    all_path <- rbind(all_path,xx)
  }
}
all_path <- all_path[,c(1,4,2,5,3,6)]
colnames(all_path) <- c("id1","id2","name1","name2","gene1","gene2")
rownames(all_path) <- c()

#计算通路之间的crosstalk
scorefun <- function(pathway){
  p1 <- unlist(strsplit(pathway[5],"\\|"))
  p2 <- unlist(strsplit(pathway[6],"\\|"))
  
  g_id1 <- which(PPI[,3] %in% p1)
  g_id2 <- which(PPI[,4] %in% p2)
  id3 <- intersect(na.omit(g_id1),na.omit(g_id2))
  sum1 <- 0
  if(length(id3)!=0){
    for (i in 1:length(id3)) {
      x <- PPI$ENTREZID1[id3[i]]
      y <- PPI$ENTREZID2[id3[i]]
      if(x%in%rownames(exp) & y%in%rownames(exp)){
        s1 <- (-1)*log(t.test(exp[x,tumor],exp[x,control],alternative = "two.sided")$p.value)+
          (-1)*log(t.test(exp[y,tumor],exp[y,control],alternative = "two.sided")$p.value)+
          (-1)*log(cor.test(as.numeric(exp[x,]),as.numeric(exp[y,]),alternative = "two.sided")$p.value)
      }
      else{s1=0}
      sum1=sum1+s1
    }
  }
  
  g_id3 <- which(PPI[,4] %in% p1)
  g_id4 <- which(PPI[,3] %in% p2)
  id4 <- intersect(na.omit(g_id3),na.omit(g_id4))
  id5 <- setdiff(id4,id3)
  sum2 <- 0
  if(length(id5)!=0){
    for (i in 1:length(id5)) {
      m <- PPI$ENTREZID2[id5[i]]
      n <- PPI$ENTREZID1[id5[i]]
      if(m%in%rownames(exp) & n%in%rownames(exp)){
        s2 <- (-1)*log(t.test(exp[m,tumor],exp[m,control],alternative = "two.sided")$p.value)+
          (-1)*log(t.test(exp[n,tumor],exp[n,control],alternative = "two.sided")$p.value)+
          (-1)*log(cor.test(as.numeric(exp[m,]),as.numeric(exp[n,]),alternative = "two.sided")$p.value)
      }
      else{s2=0}
      sum2=sum2+s2
    }
  }
  sum3 <- sum1+sum2
  return(sum3)
}
score1 <- apply(all_path,1,function(x) scorefun(x))
int_path <- all_path[which(score1!=0&score1!=Inf),]
score2 <- score1[which(score1!=0&score1!=Inf)]
row.names(int_path) <- c()
#cross <- cbind.data.frame(int_path[,c(3:6)],score2)
#write.table(cross,"crosstalk.txt",sep = "\t",quote = F,row.names = F)

scorefun1 <- function(z){
  sum9 <- NULL
  for (k in 1:nrow(int_path)){
   p1 <- unlist(strsplit(int_path[k,5],"\\|"))
   p2 <- unlist(strsplit(int_path[k,6],"\\|"))
   p3 <- unlist(strsplit(z,"\\|"))
  
   p4 <- setdiff(p1,p3)
   p5 <- setdiff(p2,p3)
  
   g_id1 <- which(PPI[,3] %in% p4)
   g_id2 <- which(PPI[,4] %in% p5)
   id3 <- intersect(na.omit(g_id1),na.omit(g_id2))
   sum1 <- 0
   if(length(id3)!=0){
    for (i in 1:length(id3)) {
      x <- PPI$ENTREZID1[id3[i]]
      y <- PPI$ENTREZID2[id3[i]]
      if(x%in%rownames(exp) & y%in%rownames(exp)){
        s1 <- (-1)*log(t.test(exp[x,tumor],exp[x,control],alternative = "two.sided")$p.value)+
          (-1)*log(t.test(exp[y,tumor],exp[y,control],alternative = "two.sided")$p.value)+
          (-1)*log(cor.test(as.numeric(exp[x,]),as.numeric(exp[y,]),alternative = "two.sided")$p.value)
      }
      else{s1=0}
      sum1=sum1+s1
     }
   }
  
   g_id3 <- which(PPI[,4] %in% p4)
   g_id4 <- which(PPI[,3] %in% p5)
   id4 <- intersect(na.omit(g_id3),na.omit(g_id4))
   id5 <- setdiff(id4,id3)
   sum2 <- 0
   if(length(id5)!=0){
    for (i in 1:length(id5)) {
      m <- PPI$ENTREZID2[id5[i]]
      n <- PPI$ENTREZID1[id5[i]]
      if(m%in%rownames(exp) & n%in%rownames(exp)){
        s2 <- (-1)*log(t.test(exp[m,tumor],exp[m,control],alternative = "two.sided")$p.value)+
          (-1)*log(t.test(exp[n,tumor],exp[n,control],alternative = "two.sided")$p.value)+
          (-1)*log(cor.test(as.numeric(exp[m,]),as.numeric(exp[n,]),alternative = "two.sided")$p.value)
      }
      else{s2=0}
      sum2=sum2+s2
     }
   }
   sum3 <- sum1+sum2
   sum9 <- c(sum9,sum3)
 }
 des <- 1-(sum(sum9)/sum(score2))
 return(des)
}

#导入药物
drug <- read.table("drug_target.txt",sep = "\t",header = T,stringsAsFactors = F)
drug_s <- apply(drug,1,function(x) scorefun1(x[2]))
drug_ss <- cbind.data.frame(drug[which(drug_s>0.01),],drug_s[which(drug_s>0.01)])
colnames(drug_ss) <- c("drug","target","score")
drug_t <- drug[which(drug_s>0.01),]

#计算药物组合的破坏程度
d_t <- NULL
for(i in 1:(nrow(drug_t)-1)){
  for (j in (i+1):nrow(drug_t)){
    xx <- paste(union(unlist(strsplit(drug_t[i,2],"\\|")),unlist(strsplit(drug_t[j,2],"\\|"))),collapse = "|")
    yy <- data.frame(drug_t[i,1],drug_t[j,1],xx)
    colnames(yy) <- ""
    d_t <- rbind(d_t,yy)
  }
}
colnames(d_t) <- c("id1","id2","target")
rownames(d_t) <- c()

des_all <- apply(d_t,1,function(x) scorefun1(x[3]))
drug_zh <- cbind.data.frame(d_t,des_all)


#比较两药物相加和药物组合的破坏程度
drug_ll <- NULL
for(i in 1:(nrow(drug_ss)-1)){
  for (j in (i+1):nrow(drug_ss)){
    zz <- drug_ss[i,3]+drug_ss[j,3]
    yy <- data.frame(drug_ss[i,1],drug_ss[j,1],zz)
    colnames(yy) <- ""
    drug_ll <- rbind(drug_ll,yy)
  }
}
colnames(drug_ll) <- c("id1","id2","score")
rownames(drug_ll) <- c()
drug_zz <- drug_zh[which(des_all>drug$score),]

write.table(drug_ss[,c(1,3)],"drug_1.txt",sep = "\t",quote = F,row.names = F)
write.table(drug_zz[,c(1,2,4)],"drug_2.txt",sep = "\t",quote = F,row.names = F)



